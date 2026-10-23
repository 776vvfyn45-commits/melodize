/* Melodize DSP — pitch tracking, note segmentation, key detection, quantizing,
   beatbox drum-hit detection, WAV encoding.
   Pure functions, no DOM / Web Audio dependency: runs in the browser AND in Node (for tests).
   Pitch tracker: FFT-accelerated YIN with multiple candidates per frame + Viterbi
   continuity decoding (a "pYIN-lite"), which removes the classic octave-jump errors. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.MelodizeDSP = factory();
})(typeof self !== 'undefined' ? self : this, function () {
'use strict';

const NOTE_NAMES = ['C','C#','D','D#','E','F','F#','G','G#','A','A#','B'];
const freqToMidi = f => 69 + 12 * Math.log2(f / 440);
const midiToFreq = m => 440 * Math.pow(2, (m - 69) / 12);
function midiToNoteName(m) {
  const n = Math.round(m);
  return NOTE_NAMES[((n % 12) + 12) % 12] + (Math.floor(n / 12) - 1);
}

/* ---------- small helpers ---------- */
function median(arr) {
  if (!arr.length) return 0;
  const a = Array.from(arr).sort((x, y) => x - y);
  const m = a.length >> 1;
  return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
}
function percentile(arr, p) {
  if (!arr.length) return 0;
  const a = Array.from(arr).sort((x, y) => x - y);
  return a[Math.min(a.length - 1, Math.max(0, Math.floor(p * (a.length - 1))))];
}

/* ---------- FFT (iterative radix-2) ---------- */
function makeFFT(n) {
  const levels = Math.round(Math.log2(n));
  if ((1 << levels) !== n) throw new Error('FFT size must be a power of 2');
  const cos = new Float64Array(n / 2), sin = new Float64Array(n / 2);
  for (let i = 0; i < n / 2; i++) { cos[i] = Math.cos(2 * Math.PI * i / n); sin[i] = Math.sin(2 * Math.PI * i / n); }
  const rev = new Uint32Array(n);
  for (let i = 0; i < n; i++) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (levels - 1));
  function transform(re, im, inverse) {
    for (let i = 0; i < n; i++) {
      const j = rev[i];
      if (j > i) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let size = 2; size <= n; size <<= 1) {
      const half = size >> 1, step = n / size;
      for (let i = 0; i < n; i += size) {
        for (let j = i, k = 0; j < i + half; j++, k += step) {
          const c = cos[k], s = inverse ? sin[k] : -sin[k];
          const tre = re[j + half] * c - im[j + half] * s;
          const tim = re[j + half] * s + im[j + half] * c;
          re[j + half] = re[j] - tre; im[j + half] = im[j] - tim;
          re[j] += tre; im[j] += tim;
        }
      }
    }
    if (inverse) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  }
  return { n, transform };
}

/* ---------- 2nd-order high-pass biquad (RBJ), removes DC + rumble ---------- */
function highpass(data, sampleRate, fc) {
  fc = fc || 80;
  const w0 = 2 * Math.PI * fc / sampleRate, cosw = Math.cos(w0), alpha = Math.sin(w0) / (2 * 0.7071);
  const b0 = (1 + cosw) / 2, b1 = -(1 + cosw), b2 = (1 + cosw) / 2;
  const a0 = 1 + alpha, a1 = -2 * cosw, a2 = 1 - alpha;
  const out = new Float32Array(data.length);
  let x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  for (let i = 0; i < data.length; i++) {
    const x0 = data[i];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    out[i] = y0; x2 = x1; x1 = x0; y2 = y1; y1 = y0;
  }
  return out;
}

/* ---------- YIN frame analysis using FFT autocorrelation ---------- */
function makeYin(sampleRate, opts) {
  opts = opts || {};
  const W = opts.windowSize || 2048;
  const half = W >> 1;
  const minF = opts.minFreq || 70, maxF = opts.maxFreq || 1000;
  const N = 4096;
  const fft = makeFFT(N);
  const tauMin = Math.max(2, Math.floor(sampleRate / maxF));
  const tauMax = Math.min(half - 2, Math.ceil(sampleRate / minF));
  const reA = new Float64Array(N), imA = new Float64Array(N);
  const reB = new Float64Array(N), imB = new Float64Array(N);
  const cum = new Float64Array(W + 1);
  const cmnd = new Float64Array(half);

  // returns { cands:[{midi, cmnd}], rms }
  function analyze(buf, offset) {
    reA.fill(0); imA.fill(0); reB.fill(0); imB.fill(0);
    cum[0] = 0;
    let sumSq = 0;
    for (let i = 0; i < W; i++) {
      const v = buf[offset + i];
      reB[i] = v; if (i < half) reA[i] = v;
      sumSq += v * v; cum[i + 1] = sumSq;
    }
    const rms = Math.sqrt(sumSq / W);
    if (rms < 1e-6) return { cands: [], rms };
    fft.transform(reA, imA, false); fft.transform(reB, imB, false);
    // conj(A) * B  -> cross-correlation
    for (let k = 0; k < N; k++) {
      const ar = reA[k], ai = -imA[k], br = reB[k], bi = imB[k];
      reA[k] = ar * br - ai * bi; imA[k] = ar * bi + ai * br;
    }
    fft.transform(reA, imA, true);
    const e0 = cum[half];
    cmnd[0] = 1;
    let running = 0;
    for (let tau = 1; tau < half; tau++) {
      const d = e0 + (cum[tau + half] - cum[tau]) - 2 * reA[tau];
      const dd = d < 0 ? 0 : d;
      running += dd;
      cmnd[tau] = running === 0 ? 1 : dd * tau / running;
    }
    const cands = [];
    for (let tau = tauMin; tau <= tauMax; tau++) {
      if (cmnd[tau] < cmnd[tau - 1] && cmnd[tau] <= cmnd[tau + 1] && cmnd[tau] < 0.4) {
        const s0 = cmnd[tau - 1], s1 = cmnd[tau], s2 = cmnd[tau + 1];
        const denom = 2 * (2 * s1 - s2 - s0);
        const better = denom !== 0 ? tau + (s2 - s0) / denom : tau;
        const f = sampleRate / better;
        cands.push({ midi: freqToMidi(f), cmnd: s1, tau });
      }
    }
    cands.sort((a, b) => a.cmnd - b.cmnd);
    // make sure octave relatives of the strongest candidates are always available to the Viterbi decoder
    const base = cands.slice(0, 3);
    for (const c of base) {
      for (const target of [Math.round(c.tau / 2), c.tau * 2]) {
        if (target < tauMin + 1 || target > tauMax) continue;
        let bt = -1, bv = 9;
        for (let t = target - 3; t <= target + 3; t++) {
          if (t < tauMin || t > tauMax) continue;
          if (cmnd[t] < bv) { bv = cmnd[t]; bt = t; }
        }
        if (bt < 0 || bv >= 0.4) continue;
        if (cands.some(o => Math.abs(o.tau - bt) <= 3)) continue;
        const s0 = cmnd[bt - 1], s1 = cmnd[bt], s2 = cmnd[bt + 1];
        const denom = 2 * (2 * s1 - s2 - s0);
        const better = (denom !== 0 && bt > 1 && bt < half - 1) ? bt + (s2 - s0) / denom : bt;
        cands.push({ midi: freqToMidi(sampleRate / better), cmnd: s1, tau: bt });
      }
    }
    cands.sort((a, b) => a.cmnd - b.cmnd);
    // keep the shortest-lag (highest-pitch) "equally good" candidates + the 2 best by score,
    // so the decoder always has the true pitch AND its octave relatives to choose from
    if (cands.length) {
      const bestC = cands[0].cmnd;
      const keep = new Set(cands.filter(c => c.cmnd <= bestC + 0.08).sort((a, b) => a.tau - b.tau).slice(0, 4));
      cands.slice(0, 2).forEach(c => keep.add(c));
      return { cands: Array.from(keep), rms };
    }
    return { cands: [], rms };
  }
  return { analyze, W, tauMin };
}

/* ---------- Pitch track: frames -> smoothed midi per frame ---------- */
function trackPitch(data, sampleRate, opts) {
  opts = opts || {};
  const hop = opts.hop || 256;
  const yin = makeYin(sampleRate, opts);
  const W = yin.W;
  const filtered = highpass(data, sampleRate, opts.hpFreq || 80);

  const raw = [];
  for (let start = 0; start + W <= filtered.length; start += hop) {
    const r = yin.analyze(filtered, start);
    let fs = 0; const c0 = start + (W >> 1) - 256;
    for (let q = 0; q < 512; q++) { const v = filtered[Math.min(filtered.length - 1, Math.max(0, c0 + q))]; fs += v * v; }
    raw.push({ t: (start + W / 2) / sampleRate, cands: r.cands, rms: r.rms, fast: Math.sqrt(fs / 512) });
  }
  if (!raw.length) return { frames: [], hopSec: hop / sampleRate, gate: 0 };

  // adaptive gate from the recording's own noise floor / peak level
  const rmsList = raw.map(f => f.rms);
  const floor = percentile(rmsList, 0.10), peak = percentile(rmsList, 0.95);
  let gate = Math.min(Math.max(floor * 3, peak * 0.04), peak * 0.15);
  gate = Math.max(gate, 0.0015);

  raw.forEach(f => { f.voiced = f.rms >= gate && f.cands.length > 0; });

  // Viterbi inside each contiguous voiced run: obs cost = cmnd (+ tiny low-tau bias), transition = jump cost
  const midi = new Array(raw.length).fill(NaN);
  let i = 0;
  while (i < raw.length) {
    if (!raw[i].voiced) { i++; continue; }
    let j = i;
    while (j < raw.length && raw[j].voiced) j++;
    // run is [i, j)
    let prevCost = raw[i].cands.map(c => c.cmnd + 0.015 * Math.log2(c.tau / yin.tauMin));
    const back = [];
    for (let k = i + 1; k < j; k++) {
      const cur = raw[k].cands, prv = raw[k - 1].cands;
      const cost = new Array(cur.length), bp = new Array(cur.length);
      for (let a = 0; a < cur.length; a++) {
        let best = Infinity, bi = 0;
        for (let b = 0; b < prv.length; b++) {
          const jump = Math.abs(cur[a].midi - prv[b].midi);
          const trans = 0.12 * Math.min(jump, 14) + (jump > 10.5 ? 1.2 : 0);
          const c = prevCost[b] + trans;
          if (c < best) { best = c; bi = b; }
        }
        cost[a] = best + cur[a].cmnd + 0.015 * Math.log2(cur[a].tau / yin.tauMin); bp[a] = bi;
      }
      prevCost = cost; back.push(bp);
    }
    let state = 0, bestC = Infinity;
    prevCost.forEach((c, idx) => { if (c < bestC) { bestC = c; state = idx; } });
    for (let k = j - 1; k >= i; k--) {
      midi[k] = raw[k].cands[state].midi;
      if (k > i) state = back[k - i - 1][state];
    }
    i = j;
  }

  // light 3-tap median smoothing on voiced neighbours, then drop tiny voiced blips
  const sm = midi.slice();
  for (let k = 1; k < midi.length - 1; k++) {
    if (isNaN(midi[k])) continue;
    const v = [midi[k]];
    if (!isNaN(midi[k - 1])) v.push(midi[k - 1]);
    if (!isNaN(midi[k + 1])) v.push(midi[k + 1]);
    sm[k] = median(v);
  }
  // fill 1–2 frame unvoiced gaps between two similar pitches
  for (let k = 1; k < sm.length - 3; k++) {
    if (!isNaN(sm[k]) || isNaN(sm[k - 1])) continue;
    for (let g = 1; g <= 2; g++) {
      if (!isNaN(sm[k + g])) {
        if (Math.abs(sm[k + g] - sm[k - 1]) < 0.7) { for (let q = 0; q < g; q++) sm[k + q] = sm[k - 1]; }
        break;
      }
    }
  }
  let s = 0;
  while (s < sm.length) {
    if (isNaN(sm[s])) { s++; continue; }
    let e = s; while (e < sm.length && !isNaN(sm[e])) e++;
    if (e - s < 3) for (let q = s; q < e; q++) sm[q] = NaN;
    s = e;
  }

  const frames = raw.map((f, idx) => ({ t: f.t, midi: sm[idx], rms: f.rms, fast: f.fast }));
  return { frames, hopSec: hop / sampleRate, gate, peak };
}

/* ---------- tuning offset (singer consistently sharp/flat) ---------- */
function estimateTuningOffset(frames) {
  let sx = 0, sy = 0, n = 0;
  for (const f of frames) {
    if (isNaN(f.midi)) continue;
    const frac = f.midi - Math.round(f.midi);           // -0.5..0.5
    const ang = 2 * Math.PI * frac;
    sx += Math.cos(ang); sy += Math.sin(ang); n++;
  }
  if (n < 5) return 0;
  const R = Math.sqrt(sx * sx + sy * sy) / n;
  if (R < 0.35) return 0;                                // no consistent offset to trust
  const off = Math.atan2(sy, sx) / (2 * Math.PI);        // in semitones
  return Math.max(-0.45, Math.min(0.45, off));
}

/* ---------- segmentation into notes ---------- */
function segmentNotes(frames, hopSec, tuningOffset) {
  const notes = [];
  tuningOffset = tuningOffset || 0;
  let i = 0;
  while (i < frames.length) {
    if (isNaN(frames[i].midi)) { i++; continue; }
    let j = i; while (j < frames.length && !isNaN(frames[j].midi)) j++;
    // run [i, j)
    let cur = null, pending = [];
    const finalize = (n, endT) => {
      if (!n) return;
      n.end = endT;
      notes.push(n);
    };
    const newNote = (k) => ({ vals: [frames[k].midi - tuningOffset], rms: [frames[k].rms], start: frames[k].t - hopSec / 2, end: frames[k].t + hopSec / 2, peak: frames[k].fast || frames[k].rms, dipMin: Infinity, dip: 0 });
    for (let k = i; k < j; k++) {
      const m = frames[k].midi - tuningOffset, r = frames[k].rms, rf = frames[k].fast || frames[k].rms;
      if (!cur) { cur = newNote(k); continue; }
      const med = median(cur.vals.slice(-18));
      if (Math.abs(m - med) > 0.75) {
        pending.push(k);
        if (pending.length >= 3) {
          const pv = pending.map(p => frames[p].midi - tuningOffset);
          if (Math.max(...pv) - Math.min(...pv) < 0.9) {
            finalize(cur, frames[pending[0]].t - hopSec / 2);
            cur = newNote(pending[0]);
            for (let q = 1; q < pending.length; q++) {
              cur.vals.push(frames[pending[q]].midi - tuningOffset); cur.rms.push(frames[pending[q]].rms);
              cur.end = frames[pending[q]].t + hopSec / 2; cur.peak = Math.max(cur.peak, frames[pending[q]].rms);
            }
            pending = [];
          }
        }
        continue;
      }
      // back near the note's pitch: absorb any pending jitter frames
      for (const p of pending) { cur.vals.push(frames[p].midi - tuningOffset); cur.rms.push(frames[p].rms); }
      pending = [];
      // re-articulation: energy dips then rises sharply at the same pitch ("ta ta ta")
      if (rf < cur.peak * 0.4) { cur.dip++; cur.dipMin = Math.min(cur.dipMin, rf); }
      else if (cur.dip >= 2 && rf > cur.dipMin * 1.8 && rf > cur.peak * 0.55 && (cur.end - cur.start) > 0.12) {
        finalize(cur, frames[k].t - hopSec / 2);
        cur = newNote(k);
        continue;
      }
      cur.vals.push(m); cur.rms.push(r); cur.end = frames[k].t + hopSec / 2;
      if (rf > cur.peak) cur.peak = rf;
      if (rf >= cur.peak * 0.6) { cur.dip = 0; cur.dipMin = Infinity; }
    }
    if (cur) finalize(cur, cur.end);
    i = j;
  }
  let maxRms = 0; notes.forEach(n => { n.meanRms = n.rms.reduce((a, b) => a + b, 0) / n.rms.length; maxRms = Math.max(maxRms, n.meanRms); });
  return notes.map(n => ({
    rawMidi: median(n.vals),
    start: n.start, dur: n.end - n.start,
    vel: Math.max(0.45, Math.min(1, 0.45 + 0.55 * (n.meanRms / (maxRms || 1))))
  }));
}

/* ---------- key detection (Krumhansl–Schmuckler) + scale snapping ---------- */
const KS_MAJOR = [6.35,2.23,3.48,2.33,4.38,4.09,2.52,5.19,2.39,3.66,2.29,2.88];
const KS_MINOR = [6.33,2.68,3.52,5.38,2.60,3.53,2.54,4.75,3.98,2.69,3.34,3.17];
function corr(a, b) {
  const n = a.length; let ma = 0, mb = 0;
  for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  let num = 0, da = 0, db = 0;
  for (let i = 0; i < n; i++) { num += (a[i] - ma) * (b[i] - mb); da += (a[i] - ma) ** 2; db += (b[i] - mb) ** 2; }
  return da && db ? num / Math.sqrt(da * db) : 0;
}
function detectKey(notes) {
  const hist = new Array(12).fill(0);
  notes.forEach(n => { hist[((Math.round(n.rawMidi) % 12) + 12) % 12] += Math.max(0.05, n.dur); });
  let best = { root: 0, mode: 'major', score: -2 };
  for (let r = 0; r < 12; r++) {
    const rot = hist.map((_, i) => hist[(i + r) % 12]);
    const sM = corr(rot, KS_MAJOR), sm = corr(rot, KS_MINOR);
    if (sM > best.score) best = { root: r, mode: 'major', score: sM };
    if (sm > best.score) best = { root: r, mode: 'minor', score: sm };
  }
  best.name = NOTE_NAMES[best.root] + ' ' + best.mode;
  return best;
}
const SCALES = {
  major:      [0,2,4,5,7,9,11],
  minor:      [0,2,3,5,7,8,10],
  majorPent:  [0,2,4,7,9],
  minorPent:  [0,3,5,7,10],
  blues:      [0,3,5,6,7,10]
};
function snapToScale(midi, root, scaleName) {
  const sc = SCALES[scaleName];
  if (!sc) return Math.round(midi);
  let best = Math.round(midi), bestD = Infinity;
  const base = Math.floor(midi / 12) * 12;
  for (let oct = -12; oct <= 12; oct += 12) {
    for (const d of sc) {
      const c = base + oct + (((root + d) % 12 + 12) % 12);
      const dist = Math.abs(c - midi);
      if (dist < bestD) { bestD = dist; best = c; }
    }
  }
  return best;
}

/* ---------- rhythm quantize ---------- */
function quantizeNotes(notes, bpm, stepsPerBeat) {
  if (!stepsPerBeat || !bpm) return notes.map(n => Object.assign({}, n));
  const grid = 60 / bpm / stepsPerBeat;
  let out = notes.map(n => {
    const s = Math.round(n.start / grid) * grid;
    let e = Math.round((n.start + n.dur) / grid) * grid;
    if (e <= s) e = s + grid;
    return Object.assign({}, n, { start: s, dur: e - s });
  });
  out.sort((a, b) => a.start - b.start);
  // monophonic cleanup: same-start -> keep the longer; trim overlaps
  const clean = [];
  for (const n of out) {
    const prev = clean[clean.length - 1];
    if (prev && Math.abs(prev.start - n.start) < 1e-6) { if (n.dur > prev.dur) clean[clean.length - 1] = n; continue; }
    if (prev && prev.start + prev.dur > n.start) prev.dur = Math.max(grid, n.start - prev.start);
    clean.push(n);
  }
  return clean;
}

/* ---------- full voice -> notes pipeline ---------- */
function transcribe(data, sampleRate, opts) {
  opts = opts || {};
  const tracked = trackPitch(data, sampleRate, opts);
  const tuning = opts.autoTune === false ? 0 : estimateTuningOffset(tracked.frames);
  let segs = segmentNotes(tracked.frames, tracked.hopSec, tuning).filter(n => n.dur >= (opts.minNote || 0.07));
  let key = segs.length ? detectKey(segs) : null;
  let root = null;
  let scaleName = opts.scale;
  if (scaleName === 'auto') {            // follow the detected key's own major/minor scale when confident
    if (key && segs.length >= 4 && key.score > 0.55 && typeof opts.keyRoot !== 'number') { scaleName = key.mode; }
    else if (typeof opts.keyRoot === 'number') { scaleName = (key && key.mode === 'minor') ? 'minor' : 'major'; }
    else scaleName = 'off';
  }
  opts = Object.assign({}, opts, { scale: scaleName });
  if (opts.scale && opts.scale !== 'off') {
    if (typeof opts.keyRoot === 'number') root = opts.keyRoot;
    else if (key) {
      root = key.root;
      const majorish = opts.scale === 'major' || opts.scale === 'majorPent';
      const minorish = opts.scale === 'minor' || opts.scale === 'minorPent' || opts.scale === 'blues';
      if (majorish && key.mode === 'minor') root = (key.root + 3) % 12;
      if (minorish && key.mode === 'major') root = (key.root + 9) % 12;
    }
  }
  let notes = segs.map(n => ({
    midi: root !== null ? snapToScale(n.rawMidi, root, opts.scale) : Math.round(n.rawMidi),
    rawMidi: n.rawMidi, start: n.start, dur: Math.max(0.1, n.dur), vel: n.vel
  }));
  // merge consecutive identical pitches that touch (tiny gaps) only when no articulation dip split them
  if (opts.bpm && opts.grid) notes = quantizeNotes(notes, opts.bpm, opts.grid);
  return { notes, key, tuningCents: Math.round(tuning * 100), gate: tracked.gate };
}

/* ---------- beatbox: detect drum hits & classify (kick / snare / hat / open hat) ---------- */
function detectDrumHits(data, sampleRate, opts) {
  opts = opts || {};
  const x = highpass(data, sampleRate, 35);
  const block = Math.max(32, Math.round(sampleRate * 0.004));
  const nb = Math.floor(x.length / block);
  const env = new Float64Array(nb);
  for (let b = 0; b < nb; b++) {
    let s = 0; for (let i = 0; i < block; i++) { const v = x[b * block + i]; s += v * v; }
    env[b] = Math.sqrt(s / block);
  }
  // trailing max-hold removes the ripple low-frequency kick waveforms put on the envelope
  for (let b = nb - 1; b >= 0; b--) { let m = env[b]; for (let q = 1; q <= 3 && b - q >= 0; q++) if (env[b - q] > m) m = env[b - q]; env[b] = m; }
  const peak = percentile(env, 0.99), floor = percentile(env, 0.2);
  const thr = Math.max(peak * 0.06, floor * 4, 0.003);
  const refractory = Math.round(0.085 / (block / sampleRate));
  const fft = makeFFT(1024);
  const re = new Float64Array(1024), im = new Float64Array(1024);
  const hits = [];
  let last = -1e9;
  for (let b = 5; b < nb - 2; b++) {
    if (env[b] < thr) continue;
    let minPrev = Infinity; for (let q = 1; q <= 5; q++) minPrev = Math.min(minPrev, env[b - q]);
    if (!(env[b] > Math.max(minPrev, 1e-5) * 2.2)) continue;
    if (b - last < refractory) continue;
    // walk back to the true start of the rise
    let s = b; while (s > 1 && env[s - 1] < env[s]) s--;
    last = b;
    const startSample = s * block;
    if (startSample + 1024 > x.length) continue;
    re.fill(0); im.fill(0);
    for (let i = 0; i < 1024; i++) re[i] = x[startSample + i] * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / 1023));
    fft.transform(re, im, false);
    let total = 0, low = 0, high = 0, cent = 0;
    const binHz = sampleRate / 1024;
    for (let k = 1; k < 512; k++) {
      const p = re[k] * re[k] + im[k] * im[k];
      const f = k * binHz;
      total += p; cent += p * f;
      if (f < 250) low += p; else if (f > 3500) high += p;
    }
    const lowFrac = total ? low / total : 0, highFrac = total ? high / total : 0, centroid = total ? cent / total : 0;
    let voice;
    if (lowFrac > 0.68) voice = 'kick';
    else if (lowFrac > 0.10) voice = 'snare';
    else if (highFrac > 0.42 || centroid > 4200) voice = 'hat';
    else voice = 'snare';
    // decay length -> open hat
    let e = b; const pk = env[b];
    while (e < nb - 1 && env[e] > pk * 0.25 && (e - b) * block / sampleRate < 0.5) e++;
    const decay = (e - b) * block / sampleRate;
    if (voice === 'hat' && decay > 0.14) voice = 'openhat';
    hits.push({ t: startSample / sampleRate, voice, strength: Math.min(1, env[b] / (peak || 1)), lowFrac, highFrac, centroid });
  }
  return hits;
}

function hitsToPattern(hits, bpm, stepsPerBeat) {
  stepsPerBeat = stepsPerBeat || 4;
  const stepDur = 60 / bpm / stepsPerBeat;
  const stepsPerBar = stepsPerBeat * 4;
  let maxStep = 0;
  const placed = hits.map(h => { const s = Math.round(h.t / stepDur); maxStep = Math.max(maxStep, s); return { step: s, voice: h.voice, strength: h.strength }; });
  const bars = Math.max(1, Math.ceil((maxStep + 1) / stepsPerBar));
  const steps = bars * stepsPerBar;
  const voices = {};
  placed.forEach(p => { (voices[p.voice] = voices[p.voice] || new Array(steps).fill(0))[p.step] = 1; });
  return { voices, steps, bars };
}

/* ---------- WAV encoding (16-bit PCM) ---------- */
function encodeWav(channels, sampleRate) {
  const nCh = channels.length, len = channels[0].length;
  const buf = new ArrayBuffer(44 + len * nCh * 2);
  const v = new DataView(buf);
  const wr = (o, s) => { for (let i = 0; i < s.length; i++) v.setUint8(o + i, s.charCodeAt(i)); };
  wr(0, 'RIFF'); v.setUint32(4, 36 + len * nCh * 2, true); wr(8, 'WAVE'); wr(12, 'fmt ');
  v.setUint32(16, 16, true); v.setUint16(20, 1, true); v.setUint16(22, nCh, true);
  v.setUint32(24, sampleRate, true); v.setUint32(28, sampleRate * nCh * 2, true);
  v.setUint16(32, nCh * 2, true); v.setUint16(34, 16, true); wr(36, 'data'); v.setUint32(40, len * nCh * 2, true);
  let o = 44;
  for (let i = 0; i < len; i++) for (let c = 0; c < nCh; c++) {
    const s = Math.max(-1, Math.min(1, channels[c][i]));
    v.setInt16(o, s < 0 ? s * 0x8000 : s * 0x7FFF, true); o += 2;
  }
  return buf;
}

return { NOTE_NAMES, freqToMidi, midiToFreq, midiToNoteName, median, makeFFT, highpass,
  trackPitch, estimateTuningOffset, segmentNotes, detectKey, SCALES, snapToScale,
  quantizeNotes, transcribe, detectDrumHits, hitsToPattern, encodeWav };
});
