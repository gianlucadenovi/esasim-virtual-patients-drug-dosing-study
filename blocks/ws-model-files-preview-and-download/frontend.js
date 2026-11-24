// ── EsaSim block ──────────────────────────────────────────────────────
// Model files (preview and download) (ws-model-files-preview-and-download) · version 1 · view, JavaScript: function render(ctx)
// Puts back together the files a block streams on a telemetry channel, shows the 3D model among them and offers every file to download
// Workspace: Gianluca De Novi's workspace · Author: Gianluca De Novi
// Saved: 2026-11-24
//
// Inputs   ctx.input("id") gives { files, stream }; ctx.rows("id"), ctx.json("id"), ctx.text("id")
//   stream      Stream                    Live telemetry · The files, in pieces, on the channel
// Settings ctx.params.key
//   channel     text      Channel = "model_files"
//   show        text      Show = "patient.glb"
// EsaSim rewrites these lines from the block's definition each time it is saved: change the definition, not them.
// ── end of EsaSim block header ────────────────────────────────────────

// Model files: puts back together the files a block streams on a telemetry channel (gzip, base64, in pieces), shows
// the 3D model among them and offers every file to download, for when a worker's result files do not reach EsaSim.
function render(ctx) {
  var THREE = ctx.THREE, p = ctx.params || {};
  var channel = p.channel || 'model_files', show = (p.show || 'patient.glb').toLowerCase();
  var files = {}, done = {}, order = [], shown = null, lastNew = Date.now(), asked = 0, told = false;

  ctx.el.innerHTML = '';
  var root = document.createElement('div');
  root.style.cssText = 'position:relative;width:100%;height:100%;min-height:420px;font:13px system-ui,sans-serif;color:#dfe6f3;background:#0b0f17;overflow:hidden';
  ctx.el.appendChild(root);
  var canvasBox = document.createElement('div');
  canvasBox.style.cssText = 'position:absolute;inset:0';
  root.appendChild(canvasBox);
  var panel = document.createElement('div');
  panel.style.cssText = 'position:absolute;top:10px;left:10px;max-width:340px;background:rgba(14,19,30,.86);border:1px solid #2a3550;border-radius:8px;padding:10px 12px';
  root.appendChild(panel);
  var title = document.createElement('div');
  title.style.cssText = 'font-weight:600;margin-bottom:6px';
  title.textContent = 'Model files';
  panel.appendChild(title);
  var status = document.createElement('div');
  status.style.cssText = 'color:#9aa3b2;margin-bottom:8px';
  status.textContent = 'Waiting for the files (run the flow)…';
  panel.appendChild(status);
  var list = document.createElement('div');
  panel.appendChild(list);
  var all = document.createElement('button');
  all.textContent = 'Download all';
  all.style.cssText = btnCss() + ';margin-top:8px;display:none';
  panel.appendChild(all);
  var hint = document.createElement('div');
  hint.style.cssText = 'color:#9aa3b2;margin-top:8px;font-size:12px;display:none';
  hint.textContent = 'Put them in the project\'s Assets, in the folder models/.';
  panel.appendChild(hint);

  function btnCss() {
    return 'background:#2d5bd7;color:#fff;border:0;border-radius:5px;padding:4px 10px;cursor:pointer;font:12px system-ui,sans-serif';
  }
  function kb(n) { return n > 1048576 ? (n / 1048576).toFixed(1) + ' MB' : Math.max(1, Math.round(n / 1024)) + ' KB'; }

  // ── the 3D preview ──
  var renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
  canvasBox.appendChild(renderer.domElement);
  var scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0b0f17);
  var camera = new THREE.PerspectiveCamera(35, 1, 0.01, 100);
  camera.up.set(0, 0, 1);
  scene.add(new THREE.HemisphereLight(0xdfe8ff, 0x202430, 1.1));
  var sun = new THREE.DirectionalLight(0xffffff, 1.6);
  sun.position.set(1.5, 3, 4);
  scene.add(sun);
  var back = new THREE.DirectionalLight(0x5c8dff, 0.8);
  back.position.set(-2, -3, 2);
  scene.add(back);
  var orbit = { az: 0.35, el: 0.15, dist: 2, target: new THREE.Vector3(0, 0, 1) };
  function place() {
    var c = Math.cos(orbit.el);
    camera.position.set(orbit.target.x + orbit.dist * c * Math.sin(orbit.az), orbit.target.y + orbit.dist * c * Math.cos(orbit.az),
                        orbit.target.z + orbit.dist * Math.sin(orbit.el));
    camera.lookAt(orbit.target);
    renderer.render(scene, camera);
  }
  function resize() {
    var w = root.clientWidth || 600, h = root.clientHeight || 420;
    renderer.setSize(w, h);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    place();
  }
  var drag = null;
  renderer.domElement.addEventListener('pointerdown', function (e) { drag = [e.clientX, e.clientY]; renderer.domElement.setPointerCapture(e.pointerId); });
  renderer.domElement.addEventListener('pointerup', function () { drag = null; });
  renderer.domElement.addEventListener('pointermove', function (e) {
    if (!drag) return;
    orbit.az -= (e.clientX - drag[0]) * 0.008;
    orbit.el = Math.max(-1.4, Math.min(1.4, orbit.el + (e.clientY - drag[1]) * 0.008));
    drag = [e.clientX, e.clientY];
    place();
  });
  renderer.domElement.addEventListener('wheel', function (e) {
    e.preventDefault();
    orbit.dist *= Math.exp(e.deltaY * 0.001);
    place();
  }, { passive: false });
  ctx.on('resize', resize);
  resize();

  function parseGlb(buf) {
    var dv = new DataView(buf), jl = dv.getUint32(12, true);
    var js = JSON.parse(new TextDecoder().decode(new Uint8Array(buf, 20, jl)));
    var bin = 20 + jl + 8;
    var types = { 5126: Float32Array, 5125: Uint32Array, 5123: Uint16Array, 5121: Uint8Array, 5122: Int16Array, 5120: Int8Array };
    var sizes = { SCALAR: 1, VEC2: 2, VEC3: 3, VEC4: 4 };
    function acc(i) {
      var a = js.accessors[i], v = js.bufferViews[a.bufferView], T = types[a.componentType];
      var off = bin + (v.byteOffset || 0) + (a.byteOffset || 0);
      var arr = new T(buf.slice(off, off + a.count * sizes[a.type] * T.BYTES_PER_ELEMENT));
      if (a.normalized || T === Int16Array || T === Int8Array) {
        var f = new Float32Array(arr.length), s = T === Int16Array ? 32767 : T === Int8Array ? 127 : 1;
        for (var k = 0; k < arr.length; k++) f[k] = a.normalized ? Math.max(arr[k] / s, -1) : arr[k];
        return f;
      }
      return arr;
    }
    var mats = (js.materials || []).map(function (m) {
      var pb = m.pbrMetallicRoughness || {}, c = pb.baseColorFactor || [0.8, 0.8, 0.8, 1];
      var op = c[3] == null ? 1 : c[3];
      return new THREE.MeshStandardMaterial({
        color: new THREE.Color(c[0], c[1], c[2]), roughness: pb.roughnessFactor == null ? 0.5 : pb.roughnessFactor,
        metalness: pb.metallicFactor || 0, transparent: op < 1, opacity: op, depthWrite: op >= 1,
        side: m.doubleSided ? THREE.DoubleSide : THREE.FrontSide,
        emissive: m.emissiveFactor ? new THREE.Color(m.emissiveFactor[0], m.emissiveFactor[1], m.emissiveFactor[2]) : new THREE.Color(0, 0, 0)
      });
    });
    function node(i) {
      var n = js.nodes[i], g = new THREE.Group();
      if (n.translation) g.position.fromArray(n.translation);
      if (n.rotation) g.quaternion.fromArray(n.rotation);
      if (n.scale) g.scale.fromArray(n.scale);
      if (n.mesh != null) {
        js.meshes[n.mesh].primitives.forEach(function (pr) {
          var geo = new THREE.BufferGeometry();
          geo.setAttribute('position', new THREE.BufferAttribute(acc(pr.attributes.POSITION), 3));
          if (pr.attributes.NORMAL != null) geo.setAttribute('normal', new THREE.BufferAttribute(acc(pr.attributes.NORMAL), 3));
          else geo.computeVertexNormals();
          if (pr.indices != null) geo.setIndex(new THREE.BufferAttribute(acc(pr.indices), 1));
          var mesh = new THREE.Mesh(geo, pr.material != null ? mats[pr.material] : new THREE.MeshStandardMaterial());
          if (mesh.material.transparent) mesh.renderOrder = 1;
          g.add(mesh);
        });
      }
      (n.children || []).forEach(function (c) { g.add(node(c)); });
      return g;
    }
    var out = new THREE.Group();
    (js.scenes[js.scene || 0].nodes || []).forEach(function (i) { out.add(node(i)); });
    return out;
  }

  function preview(name, blob) {
    blob.arrayBuffer().then(function (buf) {
      if (shown) scene.remove(shown);
      shown = parseGlb(buf);
      scene.add(shown);
      var box = new THREE.Box3().setFromObject(shown), size = box.getSize(new THREE.Vector3());
      box.getCenter(orbit.target);
      orbit.dist = Math.max(size.x, size.y, size.z) * 2.2 + 0.05;
      place();
    }).catch(function (e) { status.textContent = 'Could not show ' + name + ': ' + e.message; });
  }

  function b64(s) {
    var bin = atob(s), u = new Uint8Array(bin.length);
    for (var i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
    return u;
  }

  function finish(name) {
    var f = files[name];
    var blob = new Blob(f.chunks.map(b64));
    var ready = f.gz ? new Response(blob.stream().pipeThrough(new DecompressionStream('gzip'))).blob() : Promise.resolve(blob);
    return ready.then(function (b) {
      done[name] = b;
      var row = document.createElement('div');
      row.style.cssText = 'display:flex;align-items:center;justify-content:space-between;gap:10px;margin:3px 0';
      var lab = document.createElement('span');
      lab.textContent = name + ' · ' + kb(b.size);
      var btn = document.createElement('button');
      btn.textContent = 'Download';
      btn.style.cssText = btnCss();
      btn.onclick = function () { save(name, b); };
      row.appendChild(lab);
      row.appendChild(btn);
      list.appendChild(row);
      if (name.toLowerCase() === show || (!shown && /\.glb$/i.test(name) && name.toLowerCase().indexOf(show.replace(/\.glb$/, '')) >= 0)) preview(name, b);
    });
  }

  function save(name, b) {
    if (typeof ctx.save === 'function') return ctx.save(name, b);
    var url = URL.createObjectURL(b), a = document.createElement('a');
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 20000);
  }
  all.onclick = function () {
    var names = order.filter(function (n) { return done[n]; }), k = 0;
    (function next() { if (k < names.length) { save(names[k], done[names[k++]]); setTimeout(next, 700); } })();
  };

  function update() {
    var x = ctx.input(ctx.ports()[0]) || {}, fresh = false;
    (x.stream || []).forEach(function (it) {
      var d = it.data;
      if (it.channel !== channel || !d || !d.file || d.b64 == null) return;
      var f = files[d.file];
      if (!f) { f = files[d.file] = { parts: d.parts, gz: !!d.gz, size: d.size, chunks: new Array(d.parts), have: 0, busy: false }; order.push(d.file); }
      if (f.chunks[d.part] == null) { f.chunks[d.part] = d.b64; f.have++; fresh = true; }
    });
    if (fresh) lastNew = Date.now();
    var got = 0, want = 0;
    order.forEach(function (n) {
      var f = files[n];
      got += f.have; want += f.parts;
      if (f.have === f.parts && !f.busy) { f.busy = true; finish(n); }
    });
    if (order.length) {
      var complete = order.filter(function (n) { return files[n].have === files[n].parts; }).length;
      status.textContent = complete === order.length ? order.length + ' files ready. Drag to turn the model, scroll to zoom.'
        : 'Receiving: ' + got + ' of ' + want + ' pieces (' + complete + ' of ' + order.length + ' files)…';
      all.style.display = hint.style.display = complete === order.length ? '' : 'none';
      if (complete === order.length && !told && ctx.send) { told = true; ctx.send(channel, { complete: true }).catch(function () {}); }
    }
  }

  // pieces still missing a while after the last one came: ask the block (while it runs) to send them again
  setInterval(function () {
    if (!order.length || Date.now() - lastNew < 8000 || asked > 5) return;
    var missing = {};
    order.forEach(function (n) {
      var f = files[n], m = [];
      for (var i = 0; i < f.parts; i++) if (f.chunks[i] == null) m.push(i);
      if (m.length) missing[n] = m;
    });
    if (!Object.keys(missing).length) return;
    asked++;
    lastNew = Date.now();
    if (ctx.send) ctx.send(channel, { missing: missing }).catch(function () {});
  }, 3000);

  ctx.on('update', update);
  update();
}
