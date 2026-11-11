# Melodize: project notes for Claude Code

## What it is
Melodize is a free, installable web app (PWA) that turns humming, singing and beatboxing into instruments and full songs. It is for non-musicians. It has a dark studio look and a friendly voice assistant called "Nova".

Live site: https://776vvfyn45-commits.github.io/melodize/
GitHub user: 776vvfyn45-commits (hosted on GitHub Pages, repo "melodize")

## Ground rules
- Everything must stay free: no paid services and no backend to maintain. The only optional paid piece is "Real AI mode", which is off by default.
- No frameworks or build step. It is plain HTML and JavaScript in single files. Keep it that way unless asked.
- The app must work offline after the first load (service worker).
- Phone first: it must work well on Android, iPhone and desktop.
- The owner is not a developer. Explain changes in plain language, keep changes small, and tell her exactly which files to upload.

## Files
- `index.html` (about 275 KB): the whole app (UI, tracks, instruments, mixer, recorder, piano roll, playlist, chat/Nova, taste profile). This is the file that goes live.
- `dsp.js` (about 99 KB): the audio-analysis brain: pitch tracking (YIN), note segmentation, key and scale detection, tempo, quantizing, MIDI export, chords and auto-band, jam-partner phrases, scene-to-music, drum-hit detection, and the `transcribe` pipeline.
- `sw.js`: service worker. The cache name is `melodize-vNN` (currently v27). Downloaded sounds (Tone.js, instrument samples, drum kits) live in a separate cache, `melodize-assets-v1`, which is kept across updates. Only bump that name if the sound files themselves must be refreshed.
- `manifest.json` and `icon-*.png`: install-to-home-screen settings and icons.
- `Melodize.html`: an older standalone copy. Do not edit it or upload it.
- `ai-worker.js` and `AI_SETUP.md`: the optional Cloudflare Worker for "Real AI mode". Setup steps are in `AI_SETUP.md`.
- `_tests/`: Node test scripts for the DSP and logic. They are not part of the live site.

## Important rules when editing
1. **Bump the cache name in `sw.js`** (for example `melodize-v26` to `melodize-v27`) every time `index.html`, `dsp.js` or other shell files change. Otherwise people keep the old version.
2. If a new file must work offline, add it to `APP_SHELL` in `sw.js`.
3. After changes, run the tests in `_tests/` with Node where they apply, and check the app in a browser.
4. Never put API keys in the files. Real AI mode uses a Cloudflare Worker secret instead.

## What already exists
Voice recording (raw PCM) and transcription to notes; automatic speech-versus-singing detection (hands-free voice AI); many instruments, synth and sampled drum kits with fills and variations; effects and a mix chain (EQ, glue compressor, loudness, export); piano-roll note editor; song timeline and sections; MIDI export and import; real vocal track recording; auto-band (chords and bass from a hum); jam partner (answer, continue, ending); scene-to-music ("describe a feeling"); taste profile that learns what the user keeps; Instrument Advisor and an advisor chat with text, voice and singing input; an offline chat brain plus optional Real AI mode.

## Main screen: the quick song maker (added Nov 2026)
The first thing people see is a 3-step "Make a song from your voice" panel (`#quickPanel`, code section "8b. QUICK SONG MAKER" in `index.html`):
1. Hum (big mic button; free timing, so the tempo follows the singer).
2. Pick an instrument (12 favourites, "Choose for me", "More instruments", or type a name; unknown names map to the closest sound via `INST_NEAR`).
3. Pick a feel (`FEELS`: Just my tune, Calm, Happy, Sad, Party, Epic, Rock) or describe one in words (uses `DSP.sceneRecipe`).
"Make my song" builds lead + band + beat, auto-mixes, arranges an intro, and sets the master preset and reverb room for the feel. Changing any choice rebuilds right away. The finished song has Higher/Lower, Slower/Faster, Fix a wrong note (piano roll), Undo/Redo (whole-song snapshots via `projectSnapshot()`), Download .wav, and Save.
The old full studio is still there, hidden behind "Show full studio" (`#studioWrap`; the choice is remembered in localStorage).
Tracks can have `enterBar`: the part keeps its place in time but stays silent until that bar. Bass lines use it so they stay lined up with the chords during the intro.

## Tone.js guard (keep it)
At the top of the main script, `Tone.PolySynth.prototype._scheduleEvent` is wrapped so a removed synth ignores late note events. Without it, removing a synth track while it plays jams Tone's timer queue: an error repeats about 16 times a second and later note-offs get stuck. Rebuilding songs (the quick maker, Nova's band) triggers this often.

## Testing without Node
This PC has no Node or Python. To test, serve the folder with a small PowerShell HttpListener script, open http://localhost:PORT in the browser pane, and replace `navigator.mediaDevices.getUserMedia` with a synthesized "dum dum" hum (vibrato, scoops, +30 cent detune, timing jitter) to drive the real record button. The browser pane blocks service workers, so test `sw.js` by running its code against a fake `caches` object.

## Unfinished plan (the next work)
- **Stage 2: Sound quality.** Reverb, mastering, stems export, better drum kits. (Started, not finished.)
- **Stage 3: Editing power.** Quantize controls, undo, key and tempo changes, automation.
- **Stage 4: Songwriting help, studio UI and project save.**
- Verify each stage with tests and screenshots, then give the owner the upload list.

## Known good practices
- Test transcription on realistic "dum dum" style singing, not just clean tones. Accuracy was improved earlier this way.
- Keep the voice assistant Nova friendly and plain-spoken, not technical.
