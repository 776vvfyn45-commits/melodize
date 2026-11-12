// Melodize service worker
// - App shell is precached (fresh from the network, bypassing the HTTP cache)
// - HTML/JS/manifest: network-first, so updates show up as soon as you're online
// - Everything else (Tone.js, instrument samples, drum kits): cache-first, so it works offline after first use
// - The app shell cache is renamed on every update; the sound cache (ASSET_CACHE) is kept,
//   so an update never throws away instruments people already downloaded.

const CACHE_NAME = 'melodize-v28';
const ASSET_CACHE = 'melodize-assets-v1';
const APP_SHELL = [
  './',
  './index.html',
  './dsp.js',
  './manifest.json',
  './icon-192.png',
  './icon-512.png'
];

self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_NAME)
      .then(cache => Promise.all(APP_SHELL.map(url =>
        fetch(new Request(url, { cache: 'reload' })).then(res => { if (res.ok) return cache.put(url, res); }).catch(() => {})
      )))
      .then(() => self.skipWaiting())
  );
});

self.addEventListener('activate', event => {
  event.waitUntil(
    caches.keys().then(keys => Promise.all(keys.map(async k => {
      if (k === CACHE_NAME || k === ASSET_CACHE || k.indexOf('melodize-') !== 0) return;   // never touch other apps on the same site
      if (/^melodize-v\d+$/.test(k)) {
        // older versions kept sounds in the app cache: move them over before deleting it
        try {
          const old = await caches.open(k), assets = await caches.open(ASSET_CACHE);
          for (const req of await old.keys()) {
            const u = new URL(req.url);
            if (u.origin !== self.location.origin) { const r = await old.match(req); if (r) await assets.put(req, r); }
          }
        } catch (e) {}
      }
      return caches.delete(k);
    }))).then(() => self.clients.claim())
  );
});

function isAppFile(url) {
  return url.origin === self.location.origin &&
    /(\/|index\.html|dsp\.js|manifest\.json)$/.test(url.pathname);
}

self.addEventListener('fetch', event => {
  const req = event.request;
  if (req.method !== 'GET') return;
  const url = new URL(req.url);

  if (isAppFile(url)) {
    // network-first for the app itself
    event.respondWith(
      fetch(req).then(res => {
        if (res.ok) { const copy = res.clone(); caches.open(CACHE_NAME).then(c => c.put(req, copy)).catch(() => {}); }
        return res;
      }).catch(() => caches.match(req).then(r => r || caches.match('./index.html')))
    );
    return;
  }

  // cache-first for libraries and samples (kept across app updates)
  event.respondWith(
    caches.match(req).then(cached => {
      if (cached) return cached;
      return fetch(req).then(res => {
        if (res.ok || res.type === 'opaque') {
          const copy = res.clone();
          const isAsset = url.origin !== self.location.origin;
          caches.open(isAsset ? ASSET_CACHE : CACHE_NAME).then(c => c.put(req, copy)).catch(() => {});
        }
        return res;
      }).catch(() => cached);
    })
  );
});
