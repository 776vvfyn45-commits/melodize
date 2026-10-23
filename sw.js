// Melodize service worker — app-shell precache + runtime cache-first for everything else
// (including the instrument samples and Tone.js, so the app keeps working offline after first load).

const CACHE_NAME = 'melodize-v1';
const APP_SHELL = [
  './',
  './index.html',
  './manifest.json',
  './icon-192.png',
  './icon-512.png'
];

self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_NAME)
      .then(cache => cache.addAll(APP_SHELL))
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

self.addEventListener('fetch', event => {
  if (event.request.method !== 'GET') return;

  event.respondWith(
    caches.match(event.request).then(cached => {
      if (cached) return cached;

      return fetch(event.request).then(response => {
        // Cache a copy of same-origin AND cross-origin (opaque) responses for offline reuse —
        // this is what lets instrument samples and Tone.js keep working without a network.
        const copy = response.clone();
        caches.open(CACHE_NAME).then(cache => {
          cache.put(event.request, copy).catch(() => {});
        });
        return response;
      }).catch(() => cached);
    })
  );
});
