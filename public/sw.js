// ChemCalc Service Worker v0.1 - Network First
var CACHE = "chemcalc-v0.1";
var STATIC = ["/", "/index.html"];
self.addEventListener("install", function(e){
  e.waitUntil(caches.open(CACHE).then(function(c){return c.addAll(STATIC);}));
  self.skipWaiting();
});
self.addEventListener("activate", function(e){
  e.waitUntil(caches.keys().then(function(ks){
    return Promise.all(ks.filter(function(k){return k!==CACHE;}).map(function(k){return caches.delete(k);}));
  }));
  self.clients.claim();
});
self.addEventListener("fetch", function(e){
  if(e.request.method!=="GET") return;
  e.respondWith(
    fetch(e.request).then(function(resp){
      if(resp&&resp.status===200){
        var clone=resp.clone();
        caches.open(CACHE).then(function(c){c.put(e.request,clone);});
      }
      return resp;
    }).catch(function(){
      return caches.match(e.request).then(function(cached){
        return cached||new Response("\u79bb\u7ebf\u4e2d\uff0c\u8bf7\u8fde\u63a5\u7f51\u7edc",{status:503});
      });
    })
  );
});
