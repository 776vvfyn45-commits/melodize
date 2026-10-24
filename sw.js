// Melodize service worker v10
// - App shell is precached (fresh from the network, bypassing the HTTP cache)
// - HTML/JS/manifest: network-first, so updates show up as soon as you're online
// - Everything else (Tone.js, instrument samples): cache-first, so it works offline after first use

const CACHE_NAME = 'melodize-v10';
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
    caches.keys().then(keys =>
      Promise.all(keys.filter(k => k !== CACHE_NAME).map(k => caches.delete(k)))
    ).then(() => self.clients.claim())
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

  // cache-first for libraries and samples
  event.respondWith(
    caches.match(req).then(cached => {
      if (cached) return cached;
      return fetch(req).then(res => {
        if (res.ok || res.type === 'opaque') {
          const copy = res.clone();
          caches.open(CACHE_NAME).then(c => c.put(req, copy)).catch(() => {});
        }
        return res;
      }).catch(() => cached);
    })
  );
});
