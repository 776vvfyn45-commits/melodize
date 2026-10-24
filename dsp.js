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

/* ---------- piano-roll editing operations (pure functions on [{midi,startB,durB,vel}]) ---------- */
const rollOps = {
  snap(b, step) { return step ? Math.round(b / step) * step : b; },
  hit(notes, beat, midi, edge) {                 // topmost note under the pointer; zone 'end' = right-edge resize handle
    for (let i = notes.length - 1; i >= 0; i--) {
      const n = notes[i];
      if (n.midi === midi && beat >= n.startB - 1e-6 && beat <= n.startB + n.durB + 1e-6)
        return { i, zone: (n.durB > edge * 1.6 && beat > n.startB + n.durB - edge) ? 'end' : 'body' };
    }
    return null;
  },
  add(notes, startB, midi, durB, vel) { return notes.concat([{ midi, startB: Math.max(0, startB), durB: Math.max(0.0625, durB), vel: vel || 0.85 }]); },
  move(notes, i, orig, dBeat, dMidi, step, lo, hi) {
    const out = notes.slice();
    out[i] = Object.assign({}, orig, {
      startB: Math.max(0, rollOps.snap(orig.startB + dBeat, step)),
      midi: Math.max(lo, Math.min(hi, orig.midi + dMidi))
    });
    return out;
  },
  resize(notes, i, orig, newDur, step) {
    const out = notes.slice();
    let d = step ? Math.max(step, rollOps.snap(newDur, step)) : Math.max(0.0625, newDur);
    out[i] = Object.assign({}, orig, { durB: d });
    return out;
  },
  remove(notes, i) { return notes.filter((_, k) => k !== i); },
  duplicate(notes, i) { const n = notes[i]; return notes.concat([Object.assign({}, n, { startB: n.startB + n.durB })]); },
  transpose(notes, semis, i, lo, hi) {
    return notes.map((n, k) => (i < 0 || k === i) ? Object.assign({}, n, { midi: Math.max(lo, Math.min(hi, n.midi + semis)) }) : n);
  },
  quantize(notes, step, i) {
    if (!step) return notes.slice();
    return notes.map((n, k) => {
      if (i >= 0 && k !== i) return n;
      const s = Math.round(n.startB / step) * step, e = Math.max(s + step, Math.round((n.startB + n.durB) / step) * step);
      return Object.assign({}, n, { startB: s, durB: e - s });
    });
  }
};

/* ---------- MIDI file export (format 1): one tempo track + one track per melody / drum track ---------- */
function encodeMidi(bpm, trackDefs) {
  const PPQ = 480, out = [];
  const w16 = v => [(v >> 8) & 255, v & 255], w32 = v => [(v >>> 24) & 255, (v >>> 16) & 255, (v >>> 8) & 255, v & 255];
  const vlq = v => { v = Math.max(0, Math.round(v)); const b = [v & 127]; while ((v >>= 7) > 0) b.unshift((v & 127) | 128); return b; };
  const str = s => Array.from(s).map(c => c.charCodeAt(0) & 127);
  const chunk = (id, bytes) => str(id).concat(w32(bytes.length), bytes);
  const mpq = Math.round(60000000 / bpm);
  const tempo = [0, 0xFF, 0x51, 3, (mpq >> 16) & 255, (mpq >> 8) & 255, mpq & 255, 0, 0xFF, 0x58, 4, 4, 2, 24, 8, 0, 0xFF, 0x2F, 0];
  const chunks = [chunk('MTrk', tempo)];
  let mc = 0;
  trackDefs.forEach(td => {
    const ch = td.drums ? 9 : (mc++ % 9);
    const ev = [];
    const name = str(td.name || 'Track'); ev.push({ t: 0, b: [0xFF, 0x03, name.length].concat(name), o: 0 });
    if (!td.drums) ev.push({ t: 0, b: [0xC0 | ch, td.program || 0], o: 0 });
    td.notes.forEach(n => {
      const on = Math.round(n.startB * PPQ), off = Math.max(on + 10, Math.round((n.startB + n.durB) * PPQ));
      const m = Math.max(0, Math.min(127, Math.round(n.midi))), v = Math.max(1, Math.min(127, Math.round((n.vel || 0.85) * 127)));
      ev.push({ t: on, b: [0x90 | ch, m, v], o: 1 });
      ev.push({ t: off, b: [0x80 | ch, m, 0], o: 0 });
    });
    ev.sort((a, b) => a.t - b.t || a.o - b.o);
    const bytes = []; let last = 0;
    ev.forEach(e => { bytes.push.apply(bytes, vlq(e.t - last)); bytes.push.apply(bytes, e.b); last = e.t; });
    bytes.push(0, 0xFF, 0x2F, 0);
    chunks.push(chunk('MTrk', bytes));
  });
  const head = chunk('MThd', w16(1).concat(w16(chunks.length), w16(PPQ)));
  const all = head.concat.apply(head, chunks);
  return new Uint8Array(all);
}

/* ---------- read back a 16-bit PCM WAV written by encodeWav (mono or stereo -> mono) ---------- */
function decodeWav(arrayBuffer) {
  const v = new DataView(arrayBuffer);
  const nCh = v.getUint16(22, true), sr = v.getUint32(24, true), bits = v.getUint16(34, true);
  if (bits !== 16) throw new Error('unsupported wav');
  let off = 12, dataOff = 0, dataLen = 0;
  while (off + 8 <= v.byteLength) {
    const id = String.fromCharCode(v.getUint8(off), v.getUint8(off + 1), v.getUint8(off + 2), v.getUint8(off + 3)), sz = v.getUint32(off + 4, true);
    if (id === 'data') { dataOff = off + 8; dataLen = Math.min(sz, v.byteLength - dataOff); break; }
    off += 8 + sz + (sz & 1);
  }
  const frames = Math.floor(dataLen / 2 / nCh), out = new Float32Array(frames);
  for (let i = 0; i < frames; i++) { let s = 0; for (let c = 0; c < nCh; c++) s += v.getInt16(dataOff + (i * nCh + c) * 2, true) / 32768; out[i] = s / nCh; }
  return { data: out, sr };
}

/* ---------- loudness: bring a finished mix to a consistent level without clipping ---------- */
function normalizeLoudness(chans, targetRmsDb, peakCeilDb) {
  let sum = 0, n = 0, peak = 0;
  for (const ch of chans) for (let i = 0; i < ch.length; i++) { const v = ch[i]; sum += v * v; n++; const a = Math.abs(v); if (a > peak) peak = a; }
  if (!n || peak < 1e-5) return 1;
  // measure RMS only over the audible part so silence/tails don't drag the average down
  let sum2 = 0, n2 = 0; const gate = peak * 0.03;
  for (const ch of chans) for (let i = 0; i < ch.length; i++) { const v = ch[i]; if (Math.abs(v) > gate) { sum2 += v * v; n2++; } }
  const rmsDb = 10 * Math.log10((sum2 / Math.max(1, n2)) + 1e-12);
  let gain = Math.pow(10, (targetRmsDb - rmsDb) / 20);
  const maxGain = Math.pow(10, peakCeilDb / 20) / peak;
  gain = Math.min(gain, maxGain, 8);
  if (gain < 0.25) gain = 0.25;
  for (const ch of chans) for (let i = 0; i < ch.length; i++) ch[i] *= gain;
  return gain;
}

/* ---------- note clean-up: sharper onsets, fewer ghost / split notes ---------- */
function energyEnvelope(data, sr, blockSec) {
  const block = Math.max(16, Math.round(sr * (blockSec || 0.004))), nb = Math.floor(data.length / block);
  const env = new Float64Array(nb);
  for (let b = 0; b < nb; b++) { let s = 0; for (let i = 0; i < block; i++) { const v = data[b * block + i]; s += v * v; } env[b] = Math.sqrt(s / block); }
  return { env, dt: block / sr };
}
function refineNotes(segs, data, sr) {
  if (!segs.length) return segs;
  const E = energyEnvelope(data, sr, 0.004), env = E.env, dt = E.dt;
  const at = t => Math.max(0, Math.min(env.length - 1, Math.round(t / dt)));
  const peakIn = (a, b) => { let m = 0; for (let i = at(a); i <= at(b); i++) if (env[i] > m) m = env[i]; return m; };
  const minIn = (a, b) => { let m = Infinity; for (let i = at(a); i <= at(b); i++) if (env[i] < m) m = env[i]; return m; };
  let out = segs.map(n => Object.assign({}, n));
  // 1) pull each start back to the true energy onset (pitch tracking only locks on ~30-50 ms late)
  out.forEach((n, idx) => {
    const prevEnd = idx ? out[idx - 1].start + out[idx - 1].dur : 0;
    const lo = Math.max(prevEnd, n.start - 0.07), pk = peakIn(n.start, Math.min(n.start + n.dur, n.start + 0.15));
    if (!pk) return;
    let k = at(n.start);
    while (k > at(lo) && env[k - 1] > pk * 0.22) k--;
    const t = k * dt;
    if (t < n.start) { n.dur += n.start - t; n.start = t; }
  });
  // 2) ghost notes: very short AND very quiet next to the loudest part
  const gmax = Math.max.apply(null, out.map(n => peakIn(n.start, n.start + n.dur)));
  out = out.filter(n => !(n.dur < 0.13 && peakIn(n.start, n.start + n.dur) < gmax * 0.18));
  // 3) octave blips / single-note glitches: short note far from both neighbours
  out = out.filter((n, i) => {
    const p = out[i - 1], q = out[i + 1];
    if (n.dur < 0.16 && p && q && Math.abs(n.rawMidi - p.rawMidi) >= 6 && Math.abs(n.rawMidi - q.rawMidi) >= 6) return false;
    return true;
  });
  // 4) heal splits: same pitch separated by a hole where the energy never really dipped
  const merged = [];
  for (const n of out) {
    const p = merged[merged.length - 1];
    if (p && Math.abs(Math.round(p.rawMidi) - Math.round(n.rawMidi)) === 0 && n.start - (p.start + p.dur) < 0.09) {
      const gapMin = minIn(p.start + p.dur * 0.7, n.start + 0.02);
      const lowerPeak = Math.min(peakIn(p.start, p.start + p.dur), peakIn(n.start, n.start + n.dur));
      if (gapMin > lowerPeak * 0.55) { p.dur = n.start + n.dur - p.start; p.vel = Math.max(p.vel, n.vel); continue; }
    }
    // a very short note between two same-pitch notes is a wobble, absorb it
    merged.push(n);
  }
  const final = [];
  for (let i = 0; i < merged.length; i++) {
    const n = merged[i], p = final[final.length - 1], q = merged[i + 1];
    if (p && q && n.dur < 0.1 && Math.round(p.rawMidi) === Math.round(q.rawMidi) && Math.abs(n.rawMidi - p.rawMidi) <= 2.5) {
      p.dur = q.start + q.dur - p.start; i++; continue;
    }
    final.push(n);
  }
  return final;
}

/* ---------- tempo from the way someone sang (onset times, seconds) ---------- */
function estimateTempo(onsets) {
  if (!onsets || onsets.length < 4) return null;
  const t0 = onsets[0], ts = onsets.map(t => t - t0);
  let best = null, scores = [];
  for (let bpm = 60; bpm <= 180; bpm += 0.5) {
    const grid = 60 / bpm / 2;                          // eighth-note grid
    let s = 0;
    for (const t of ts) s += Math.cos(2 * Math.PI * t / grid);
    s /= ts.length;
    const prior = Math.exp(-0.5 * Math.pow(Math.log2(bpm / 100) / 0.45, 2));
    const score = s * 0.8 + prior * 0.2;
    scores.push(s);
    if (!best || score > best.score) best = { bpm, score, fit: s };
  }
  if (!best || best.fit < 0.55) return null;          // onsets don't line up with any steady pulse
  return { bpm: Math.round(best.bpm), fit: Math.round(best.fit * 100) / 100 };
}

/* ---------- full voice -> notes pipeline ---------- */
function transcribe(data, sampleRate, opts) {
  opts = opts || {};
  const tracked = trackPitch(data, sampleRate, opts);
  const tuning = opts.autoTune === false ? 0 : estimateTuningOffset(tracked.frames);
  let segs = segmentNotes(tracked.frames, tracked.hopSec, tuning).filter(n => n.dur >= (opts.minNote || 0.07));
  if (opts.refine !== false) segs = refineNotes(segs, data, sampleRate).filter(n => n.dur >= (opts.minNote || 0.07));
  let tempoEst = null;
  if (opts.autoTempo && segs.length >= 4) {
    tempoEst = estimateTempo(segs.map(n => n.start));
    if (tempoEst) opts = Object.assign({}, opts, { bpm: tempoEst.bpm });
  }
  let leadSec = 0;
  if (opts.trimLead && segs.length) { const t0 = segs[0].start; leadSec = t0; segs.forEach(n => { n.start -= t0; }); }
  if (opts.latency) segs.forEach(n => { n.start = Math.max(0, n.start - opts.latency); });
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
  return { notes, key, tuningCents: Math.round(tuning * 100), gate: tracked.gate, tempo: tempoEst, bpm: opts.bpm || null, leadSec };
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

/* ---------- Instrument advisor: describe a melody and rank instruments for it ---------- */
// lo/hi = comfortable MIDI range, sweet = centre of its best register, sustain 0 (plucked) .. 1 (bowed/blown/held),
// agile = how comfortably it plays fast passages, moods / genres = 0..1 affinity
const INSTRUMENT_PROFILES = {
  piano:     { lo:36, hi:96,  sweet:66, sustain:.50, agile:.95, moods:{happy:.6,sad:.7,calm:.7,energetic:.5,dark:.4,epic:.5}, genres:{pop:.9,classical:.9,jazz:.8,cinematic:.7,folk:.5,edm:.4,rock:.5,hiphop:.6} },
  guitar:    { lo:40, hi:84,  sweet:60, sustain:.30, agile:.70, moods:{happy:.8,sad:.6,calm:.7,energetic:.5,dark:.2,epic:.2}, genres:{folk:1,pop:.8,rock:.5,jazz:.3,cinematic:.3} },
  eguitar:   { lo:40, hi:88,  sweet:62, sustain:.50, agile:.80, moods:{energetic:1,epic:.6,dark:.6,happy:.5,sad:.3,calm:.1}, genres:{rock:1,pop:.6,edm:.3,jazz:.3,hiphop:.2} },
  nylon:     { lo:40, hi:84,  sweet:62, sustain:.30, agile:.70, moods:{calm:1,sad:.7,happy:.5,dark:.2,energetic:.2,epic:.1}, genres:{folk:.9,jazz:.6,classical:.6,cinematic:.3,pop:.4} },
  bass:      { lo:28, hi:67,  sweet:46, sustain:.50, agile:.60, moods:{energetic:.7,dark:.7,happy:.4,calm:.3,sad:.3,epic:.3}, genres:{rock:.8,hiphop:.9,pop:.6,edm:.6,jazz:.8,folk:.3} },
  violin:    { lo:55, hi:100, sweet:76, sustain:1.0, agile:.70, moods:{sad:1,epic:.7,calm:.5,happy:.4,dark:.4,energetic:.3}, genres:{classical:1,cinematic:.9,folk:.7,pop:.3} },
  cello:     { lo:36, hi:81,  sweet:55, sustain:1.0, agile:.50, moods:{sad:1,dark:.8,epic:.8,calm:.6,happy:.2,energetic:.2}, genres:{classical:1,cinematic:1,folk:.4,pop:.3} },
  flute:     { lo:60, hi:96,  sweet:79, sustain:.80, agile:.80, moods:{calm:1,happy:.7,sad:.4,energetic:.2,dark:.1,epic:.2}, genres:{classical:.8,folk:.8,cinematic:.6,pop:.3,jazz:.3} },
  clarinet:  { lo:50, hi:91,  sweet:70, sustain:.80, agile:.70, moods:{calm:.7,happy:.6,sad:.5,dark:.3,energetic:.3,epic:.2}, genres:{jazz:1,classical:.8,folk:.4,cinematic:.3} },
  sax:       { lo:46, hi:84,  sweet:64, sustain:.80, agile:.80, moods:{sad:.8,calm:.6,energetic:.5,happy:.4,dark:.4,epic:.2}, genres:{jazz:1,pop:.7,rock:.5,hiphop:.4,cinematic:.3} },
  trumpet:   { lo:52, hi:88,  sweet:70, sustain:.60, agile:.80, moods:{epic:1,energetic:.9,happy:.8,dark:.2,sad:.3,calm:.1}, genres:{jazz:.9,cinematic:.8,classical:.6,pop:.5,rock:.3} },
  trombone:  { lo:40, hi:77,  sweet:55, sustain:.70, agile:.40, moods:{epic:.8,energetic:.6,dark:.6,happy:.4,sad:.3,calm:.1}, genres:{jazz:.8,cinematic:.7,classical:.5,rock:.2} },
  horn:      { lo:34, hi:77,  sweet:55, sustain:.90, agile:.30, moods:{epic:1,calm:.6,dark:.5,sad:.5,happy:.3,energetic:.2}, genres:{cinematic:1,classical:.9} },
  harp:      { lo:30, hi:100, sweet:66, sustain:.30, agile:.90, moods:{calm:1,happy:.6,sad:.5,epic:.3,dark:.1,energetic:.1}, genres:{classical:.8,cinematic:.8,folk:.5,pop:.2} },
  organ:     { lo:36, hi:84,  sweet:60, sustain:1.0, agile:.60, moods:{epic:.8,dark:.6,calm:.5,sad:.4,happy:.3,energetic:.3}, genres:{classical:.7,rock:.6,cinematic:.7,folk:.3,jazz:.4} },
  harmonium: { lo:36, hi:84,  sweet:60, sustain:1.0, agile:.50, moods:{calm:1,sad:.7,happy:.3,epic:.2,dark:.2,energetic:.1}, genres:{folk:1,classical:.4,cinematic:.2} },
  xylo:      { lo:67, hi:108, sweet:84, sustain:0.0, agile:1.0, moods:{happy:1,energetic:.6,calm:.2,sad:.1,dark:0,epic:.1}, genres:{pop:.5,folk:.5,classical:.3,cinematic:.2} },
  supersaw:  { lo:48, hi:96,  sweet:69, sustain:.70, agile:.80, moods:{energetic:1,happy:.9,epic:.9,calm:.2,sad:.2,dark:.3}, genres:{edm:1,pop:.8,cinematic:.3} },
  pluck:     { lo:48, hi:96,  sweet:69, sustain:.10, agile:1.0, moods:{energetic:.9,happy:.8,calm:.3,sad:.3,dark:.2,epic:.3}, genres:{edm:1,pop:.8,hiphop:.5} },
  pad:       { lo:36, hi:84,  sweet:58, sustain:1.0, agile:.20, moods:{calm:1,sad:.6,dark:.6,epic:.6,happy:.3,energetic:.1}, genres:{cinematic:1,edm:.7,pop:.5} },
  sub:       { lo:24, hi:55,  sweet:38, sustain:.60, agile:.40, moods:{dark:1,energetic:.7,epic:.4,calm:.2,sad:.2,happy:.1}, genres:{edm:1,hiphop:1,pop:.5} },
  wobble:    { lo:24, hi:60,  sweet:41, sustain:.70, agile:.50, moods:{energetic:1,dark:.9,epic:.4,calm:0,sad:0,happy:.2}, genres:{edm:1,hiphop:.4} },
  acid:      { lo:36, hi:72,  sweet:52, sustain:.30, agile:.90, moods:{energetic:1,dark:.6,happy:.3,epic:.2,calm:0,sad:.1}, genres:{edm:1} },
  bell:      { lo:55, hi:100, sweet:76, sustain:.20, agile:.70, moods:{calm:.8,happy:.7,sad:.6,epic:.3,dark:.2,energetic:.2}, genres:{pop:.6,cinematic:.7,edm:.5,jazz:.4} },
  chip:      { lo:48, hi:96,  sweet:72, sustain:.40, agile:1.0, moods:{happy:1,energetic:.9,calm:.1,sad:.1,dark:.1,epic:.3}, genres:{edm:.4,pop:.3} }
};
const clamp01 = x => Math.max(0, Math.min(1, x));

// notes: [{midi, start, dur}] in seconds
function melodyFeatures(notes) {
  const n = notes.length;
  if (!n) return null;
  const midis = notes.map(x => x.midi), sorted = midis.slice().sort((a, b) => a - b);
  const first = Math.min(...notes.map(x => x.start)), last = Math.max(...notes.map(x => x.start + x.dur));
  const span = Math.max(0.5, last - first);
  const avgDur = notes.reduce((a, x) => a + x.dur, 0) / n;
  const legato = clamp01(notes.reduce((a, x) => a + x.dur, 0) / span);
  let step = 0, leap = 0, rep = 0;
  for (let i = 1; i < n; i++) {
    const d = Math.abs(midis[i] - midis[i - 1]);
    if (d === 0) rep++; if (d <= 2) step++; if (d >= 5) leap++;
  }
  const intervals = Math.max(1, n - 1);
  return {
    count: n, median: median(midis), low: sorted[0], high: sorted[n - 1], range: sorted[n - 1] - sorted[0],
    avgDur, noteRate: n / span, legato,
    stepwise: step / intervals, leaps: leap / intervals, repeats: rep / intervals,
    sustainLevel: clamp01((avgDur - 0.15) / 0.85) * 0.7 + legato * 0.3,
    speed: clamp01((n / span) / 5)
  };
}

function describeMelody(f) {
  const reg = f.median < 48 ? 'low' : f.median < 60 ? 'low-to-mid' : f.median < 72 ? 'mid-range' : f.median < 84 ? 'high' : 'very high';
  const motion = f.stepwise > 0.65 ? 'smooth, stepwise' : f.leaps > 0.35 ? 'jumpy with big leaps' : 'balanced';
  const feel = f.noteRate < 1.2 ? 'slow' : f.noteRate < 3 ? 'medium-paced' : 'fast';
  const notesTxt = f.avgDur > 0.5 ? 'long held notes' : f.avgDur < 0.25 ? 'short quick notes' : 'medium-length notes';
  return reg + ' register (' + midiToNoteName(f.low) + '–' + midiToNoteName(f.high) + '), ' + motion + ', ' + feel + ' with ' + notesTxt;
}

function rangeFit(midis, prof, shift) {
  let inside = 0;
  for (const m of midis) {
    const x = m + 12 * shift;
    const out = x < prof.lo ? prof.lo - x : (x > prof.hi ? x - prof.hi : 0);
    inside += clamp01(1 - out / 8);
  }
  inside /= midis.length;
  const med = median(midis) + 12 * shift;
  return 0.6 * inside + 0.4 * (1 - Math.min(1, Math.abs(med - prof.sweet) / 18));
}

const MOOD_LABEL = { happy:'happy', sad:'sad', energetic:'energetic', calm:'calm', dark:'dark', epic:'epic' };
// notes: [{midi,...}], ctx: { mood, genre }  -> sorted list [{ key, pct, octave, reason }]
function recommendInstruments(notes, ctx) {
  const f = melodyFeatures(notes);
  if (!f) return [];
  ctx = ctx || {};
  const midis = notes.map(n => n.midi);
  const out = [];
  Object.keys(INSTRUMENT_PROFILES).forEach(key => {
    const p = INSTRUMENT_PROFILES[key];
    let bestShift = 0, bestReg = -1;
    for (let k = -2; k <= 2; k++) {
      const r = rangeFit(midis, p, k) - 0.1 * Math.abs(k);      // prefer no shift; never more than 2 octaves
      if (r > bestReg) { bestReg = r; bestShift = k; }
    }
    const reg = clamp01(bestReg);
    const sus = 1 - Math.abs(f.sustainLevel - p.sustain);
    const agile = f.speed <= p.agile ? 1 : clamp01(1 - (f.speed - p.agile) * 1.5);
    let score = 0.5 * reg + 0.3 * sus + 0.2 * agile, maxS = 1;
    const mw = ctx.mood ? (p.moods[ctx.mood] || 0) : 0, gw = ctx.genre ? (p.genres[ctx.genre] || 0) : 0;
    if (ctx.mood) { score += 0.25 * mw; maxS += 0.25; }
    if (ctx.genre) { score += 0.25 * gw; maxS += 0.25; }
    const bits = [];
    if (reg >= 0.8 && bestShift === 0) bits.push('your melody sits right in its natural range');
    else if (bestShift !== 0) bits.push("I'd move it " + (bestShift > 0 ? 'up' : 'down') + ' ' + Math.abs(bestShift) + ' octave' + (Math.abs(bestShift) > 1 ? 's' : '') + ' to fit its range');
    if (sus >= 0.75) bits.push(f.sustainLevel > 0.5 ? 'long, flowing notes suit its sustained tone' : 'short, punchy notes suit its plucked attack');
    if (f.speed > 0.55 && p.agile >= 0.8) bits.push('it handles fast passages cleanly');
    if (ctx.mood && mw >= 0.7) bits.push("it's a classic choice for a " + MOOD_LABEL[ctx.mood] + ' feel');
    if (ctx.genre && gw >= 0.7) bits.push('it fits ' + ctx.genre + ' really well');
    if (!bits.length) bits.push('a workable match for this tune');
    const reason = bits.slice(0, 3).join('; ');
    out.push({ key, pct: Math.round(100 * clamp01(score / maxS)), octave: bestShift, reason: reason.charAt(0).toUpperCase() + reason.slice(1) + '.' });
  });
  out.sort((a, b) => b.pct - a.pct);
  return out;
}

/* ---------- Chat brain: turn a typed / spoken sentence into structured intents (rule-based, no network) ---------- */
const INSTRUMENT_ALIASES = {
  eguitar:['electric guitar','e guitar','rock guitar','distorted guitar'], nylon:['nylon guitar','nylon','classical guitar','spanish guitar'],
  guitar:['acoustic guitar','guitar','acoustic'], bass:['bass guitar','electric bass','bass'], sub:['sub bass','subbass','sub','808 bass'],
  piano:['piano','keyboard','keys'], violin:['violin','fiddle'], cello:['cello'], flute:['flute'], clarinet:['clarinet'],
  sax:['saxophone','sax'], trumpet:['trumpet'], trombone:['trombone'], horn:['french horn','horn'], harp:['harp'],
  organ:['organ'], harmonium:['harmonium'], xylo:['xylophone','xylo','marimba'],
  supersaw:['supersaw','super saw','edm lead','synth lead','lead synth','lead'], pluck:['pluck','plucky'],
  pad:['synth pad','pad','ambient pad'], wobble:['wobble bass','wobble','dubstep'], acid:['acid','303'],
  bell:['electric piano','e piano','epiano','rhodes','bell','bells'], chip:['chiptune','8 bit','8-bit','chip','retro game']
};
const STYLE_WORDS = {
  bigroom:['big room','festival','drop'], techno:['techno'], house:['house','edm','club','dance'],
  trap:['trap'], boombap:['boom bap','hip hop','hiphop','rap'], dnb:['drum and bass','drum n bass','dnb','jungle'],
  rock:['rock'], pop:['pop'], disco:['disco','funk'], reggaeton:['reggaeton','dembow'],
  afro:['afro','latin','afrobeat'], lofi:['lofi','lo fi','lo-fi']
};
const STRONG_STYLES = ['house','techno','trap','dnb','reggaeton','disco','boombap','bigroom','lofi','afro'];
const MOOD_WORDS = {
  happy:['happy','joyful','cheerful','fun','bright','uplifting'], sad:['sad','melancholy','melancholic','emotional','heartbreak','heartbroken','crying'],
  energetic:['energetic','energy','hype','upbeat','party','pumped','workout','exciting'], calm:['calm','chill','relax','relaxing','relaxed','peaceful','soft','dreamy','sleepy'],
  dark:['dark','moody','scary','heavy','aggressive','gritty'], epic:['epic','heroic','powerful','grand','massive','dramatic']
};
const GENRE_WORDS = {
  pop:['pop'], edm:['edm','electronic','dance','house','techno'], rock:['rock','metal','punk'], hiphop:['hip hop','hiphop','rap','trap'],
  classical:['classical','orchestra','orchestral'], folk:['folk','country','acoustic'], jazz:['jazz','blues'], cinematic:['cinematic','film','movie','soundtrack','trailer']
};
function normText(t) {
  return ' ' + String(t || '').toLowerCase().replace(/[’']/g, '').replace(/[^a-z0-9\- ]+/g, ' ').replace(/\s+/g, ' ').trim() + ' ';
}
function findWord(n, list) {   // returns the matched keyword (longest first) or null
  const sorted = list.slice().sort((a, b) => b.length - a.length);
  for (const w of sorted) if (n.indexOf(' ' + w + ' ') !== -1 || n.indexOf(' ' + w + 's ') !== -1) return w;
  return null;
}
function parseChat(text) {
  const n = normText(text);
  const r = { text: String(text || '').trim(), instrument: null, style: null, mood: null, genre: null, bpm: null, tempoDelta: 0 };
  const has = re => re.test(n);
  // instrument: longest alias wins across all instruments
  let bestLen = 0;
  Object.keys(INSTRUMENT_ALIASES).forEach(key => {
    const w = findWord(n, INSTRUMENT_ALIASES[key]);
    if (w && w.length > bestLen) { bestLen = w.length; r.instrument = key; }
  });
  bestLen = 0;
  Object.keys(STYLE_WORDS).forEach(key => {
    const w = findWord(n, STYLE_WORDS[key]);
    if (w && w.length > bestLen) { bestLen = w.length; r.style = key; }
  });
  Object.keys(MOOD_WORDS).forEach(k => { if (!r.mood && findWord(n, MOOD_WORDS[k])) r.mood = k; });
  bestLen = 0;
  Object.keys(GENRE_WORDS).forEach(k => { const w = findWord(n, GENRE_WORDS[k]); if (w && w.length > bestLen) { bestLen = w.length; r.genre = k; } });
  let m = n.match(/ (\d{2,3}) ?(bpm|beats per minute) /) || n.match(/ (?:tempo|bpm) (?:of |to |at |is )?(\d{2,3}) /) || n.match(/ at (\d{2,3}) /);
  if (m) { const v = parseInt(m[1], 10); if (v >= 50 && v <= 220) r.bpm = Math.max(60, Math.min(200, v)); }
  if (has(/ (faster|speed (it )?up|quicker) /)) r.tempoDelta = 10;
  if (has(/ (slower|slow (it )?down) /)) r.tempoDelta = -10;
  r.greet = has(/^ (hi|hello|hey|hola|namaste|yo|good (morning|evening|afternoon)) /);
  r.help = has(/ (help|what can you do|what do you do|how do i|how does this work|commands|options) /);
  r.thanks = has(/ (thanks|thank you|thx|cheers|awesome|great job|nice) /) && !has(/ (make|add|play) /);
  r.sing = has(/ (sing|hum|humming|whistle|record|listen to me|let me sing|here is my tune|mic|my voice) /) && !has(/ (stop recording) /);
  r.beatbox = has(/ (beatbox|beat box|beatboxing) /);
  r.buildSong = has(/ (full song|whole song|complete song|build (me )?a song|make (me )?a song|song starter|arrange it|finish (my )?song|full track) /);
  const beatWord = has(/ (beat|beats|drum|drums|rhythm|groove|percussion|drumbeat) /);
  r.wantBeat = !r.beatbox && (beatWord || (r.style && STRONG_STYLES.indexOf(r.style) !== -1 && !r.instrument && !has(/ (song|track) /)));
  r.recommend = has(/ (which|what) (instrument|sound)s? |recommend|suggest|suits?|suitable|fits?|best (instrument|sound)|advice|advise|belongs? /) ;
  r.play = has(/ (play|start|hear it|listen) /) && !r.sing;
  r.stop = has(/ (stop|pause|silence|quiet) /);
  r.exportSong = has(/ (export|download|render|bounce|save as wav|save (the )?song) /);
  r.saveProject = has(/ save (the |my )?project /);
  r.undo = has(/ (undo|remove (the )?last|delete (the )?last|take (that )?back|scrap that) /);
  r.clearAll = has(/ (clear (everything|all)|start over|reset|delete everything|new project) /);
  r.fill = has(/ (fill|fills|drum fill) /);
  r.variation = has(/ (variation|vary|change it up|mix it up|switch it up) /);
  extraChat(n, text, r);
  return r;
}

/* ----- friendlier understanding: questions, follow-ups, small talk, unsupported topics, typo hints ----- */
function lev(a, b) {
  if (a === b) return 0;
  const m = a.length, k = b.length; if (!m) return k; if (!k) return m;
  let prev = []; for (let j = 0; j <= k; j++) prev[j] = j;
  for (let i = 1; i <= m; i++) {
    const cur = [i];
    for (let j = 1; j <= k; j++) cur[j] = Math.min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1));
    prev = cur;
  }
  return prev[k];
}
const UNSUPPORTED = [
  { id: 'lyrics', re: / (lyrics?|write (me )?a song about|words to (the|my) song|rap verse) /, say: "I can't write lyrics or sing words — I turn your voice into instruments and beats. But I can build the music under your tune, and you can add words on top." },
  { id: 'vocals', re: / (autotune|auto tune|pitch correct(ion)?|backing vocals?|harmon(y|ies)|choir) /, say: "Pitch correction, harmonies and choirs aren't built in. For a harmony feel you can layer the same tune on two instruments (like strings and a pad) — I can do that now. Your raw voice can be kept as a track with “keep my real voice”." },
  { id: 'mp3', re: / (mp3|midi file|export midi|stems|flac|aac) /, say: "Right now I export WAV files (studio quality), and you can convert them to MP3 in any free converter. MIDI and stem export aren't in yet." },
  { id: 'share', re: / (upload|share|post|publish|spotify|youtube|instagram|tiktok|soundcloud|whatsapp) /, say: "I can't post online for you, but export gives you a WAV you can upload anywhere." },
  { id: 'import', re: / (import|upload|load) (a |my |an )?(song|file|loop|sample|audio|track|mp3|wav) /, say: "Importing your own audio files isn't supported yet — right now everything starts from your voice, a beat or the built-in sounds." }
];
function extraChat(n, raw, r) {
  const has = re => re.test(n);
  r.question = /\?\s*$/.test(String(raw || '').trim()) || has(/^ (can|could|would|will|do|does|is|are|shall|should|how about|what about|whats|what is|why|how) /) || has(/ (is it possible|is there a way|are you able|are we able) /);
  r.canWe = has(/ (can|could|would|will|shall) (we|you|i|u|this|it)\b/) || has(/ (is it possible|is there a way|are you able|are we able|able to|how about|what about|what if) /) || has(/^ (any chance|maybe|perhaps) /);
  r.yes = has(/^ (yes|yeah|yea|yep|yup|sure|ok|okay|k|please|pls|go ahead|do it|sounds good|lets do it|lets go|why not|definitely|absolutely|of course|right|correct|that one|that works|perfect) /);
  r.no = has(/^ (no|nope|nah|not now|not really|dont|cancel|never mind|nevermind|skip|keep them|leave it|no thanks) /);
  r.again = has(/ (again|another|one more|different (one|style|sound|beat|instrument)|something else|other one|next one|new one|redo|re do|try (a )?(different|another)) /);
  r.surprise = has(/ (surprise me|you choose|you pick|anything|whatever|your choice|dealer.?s choice|i dont know|idk|not sure what|inspire me) /);
  r.improve = has(/ (better|cooler|nicer|professional|polished|fuller|richer|improve|more interesting|sound good|sounds good|fix it|make it pop) /) && has(/ (this|it|that|song|track|tune|melody|mix) /);
  r.whatNext = has(/ (what next|whats next|what now|what should i|what do i do|where do i start|get started|help me start|next step|what else) /);
  r.removeDrums = has(/ (remove|delete|drop|get rid of|take out|no more) (the |all )?(drums?|beats?|percussion) /);
  r.removeMelody = has(/ (remove|delete|drop|get rid of|take out) (the |my )?(melody|tune|lead) /);
  r.like = has(/ (i love (it|this|that)|love it|sounds (great|good|amazing|awesome|nice|perfect|fire)|nailed it|perfect|im happy|so good) /);
  r.dislike = has(/ (dont like|do not like|hate|not good|sounds (bad|off|weird|wrong|terrible|awful|boring)|doesnt sound|not what i|wrong|boring|meh|ugly|too (loud|quiet|fast|slow|busy|simple)) /);
  r.how = has(/ (how are you|hows it going|whats up|how do you do|how is it going) /);
  r.who = has(/ (who are you|what are you|your name|are you (an? )?(ai|bot|robot|real|human|person|chatgpt|claude)|who made you|who built you) /);
  r.bye = has(/^ (bye|goodbye|see you|later|cya|good night|gotta go) /);
  r.sorry = has(/ (sorry|my bad|oops|whoops) /);
  r.bored = has(/ (bored|nothing to do|stuck|no idea|no ideas|blank) /);
  r.unsupported = null;
  for (const u of UNSUPPORTED) if (u.re.test(n)) { r.unsupported = u; break; }
  r.vague = has(/ (something|stuff|thing|things|cool|nice|good) /) && !r.instrument && !r.style && !r.wantBeat && !r.recommend;
  // typo hints: only when nothing was recognised and the message is short
  r.didYouMean = null;
  if (!r.instrument && !r.style && !r.wantBeat && !r.mood && !r.genre) {
    const words = n.trim().split(' ').filter(w => w.length >= 4);
    let best = null;
    for (const w of words) {
      for (const key of Object.keys(INSTRUMENT_ALIASES)) for (const al of INSTRUMENT_ALIASES[key]) {
        if (al.indexOf(' ') !== -1 || al.length < 4) continue;
        const d = lev(w, al); if (d <= (al.length >= 8 ? 2 : 1) && w !== al && (!best || d < best.d)) best = { kind: 'instrument', key, word: w, to: al, d };
      }
      for (const key of Object.keys(STYLE_WORDS)) for (const al of STYLE_WORDS[key]) {
        if (al.indexOf(' ') !== -1 || al.length < 4) continue;
        const d = lev(w, al); if (d <= (al.length >= 8 ? 2 : 1) && w !== al && (!best || d < best.d)) best = { kind: 'style', key, word: w, to: al, d };
      }
    }
    r.didYouMean = best;
  }
}

return { NOTE_NAMES, freqToMidi, midiToFreq, midiToNoteName, median, makeFFT, highpass,
  trackPitch, estimateTuningOffset, segmentNotes, detectKey, SCALES, snapToScale,
  quantizeNotes, transcribe, refineNotes, estimateTempo, normalizeLoudness, rollOps, encodeMidi, decodeWav, detectDrumHits, hitsToPattern, encodeWav,
  INSTRUMENT_PROFILES, melodyFeatures, describeMelody, recommendInstruments, parseChat, INSTRUMENT_ALIASES };
});
