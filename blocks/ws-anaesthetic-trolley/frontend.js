// ── EsaSim block ──────────────────────────────────────────────────────
// Anaesthetic trolley (ws-anaesthetic-trolley) · version 1 · view, JavaScript: function render(ctx)
// The anaesthetist's controls as a trolley: pick a drug, set how much, choose a push or an infusion, administer; the airway and the oxygen; the job, the clock and the demo
// Workspace: Gianluca De Novi's workspace · Author: Gianluca De Novi
// Saved: 2026-11-25
//
// Inputs   ctx.input("id") gives { files, stream }; ctx.rows("id"), ctx.json("id"), ctx.text("id")
//   stream      Stream                    Live telemetry · The run's telemetry: the patient channel carries the doses, the infusions running and the advice
// Settings ctx.params.key
//   channel     text      Commands go on = "control"
//   data_channel text      Patient channel = "patient"
//   to          text      Only to the block = ""
// EsaSim rewrites these lines from the block's definition each time it is saved: change the definition, not them.
// ── end of EsaSim block header ────────────────────────────────────────

// The anaesthetist's trolley: pick a drug, set how much, choose how it goes in, and administer it. Everything the
// room can be told is here, in groups, instead of a row of controls under the 3D patient: the drugs, the airway and
// the oxygen, and the case itself (the job, the patient, the clock, the demo). It reads the patient channel for the
// doses this patient's weight comes to and for what is already running, and sends one command per action.
function render(ctx) {
  var p = ctx.params || {};
  var CH = p.channel || 'control';
  var DATA = p.data_channel || 'patient';
  var TO = p.to || '';
  var GROUPS = [
    ['Asleep', ['Propofol', 'Midazolam']],
    ['Pain and the knife', ['Remifentanil', 'Fentanyl', 'Morphine']],
    ['Not moving', ['Rocuronium']],
    ['Undo it', ['Sugammadex', 'Neostigmine']],
    ['The circulation', ['Phenylephrine', 'Ephedrine', 'Atropine']],
    ['Afterwards', ['Ondansetron']]
  ];
  var FAVOURITES = ['Propofol', 'Remifentanil', 'Rocuronium', 'Sugammadex', 'Phenylephrine'];
  var JOBS = ['Operation: a four-hour laparotomy', 'Intensive care: sedation and the wake-up test',
              'Therapeutic coma for status epilepticus'];
  var CLOCKS = ['Real time', '1 min/s', '5 min/s', '15 min/s', '1 h/s'];
  var AIRWAYS = ['No support', 'Bag-mask', 'Intubated'];

  var drug = 'Propofol', route = 'bolus', x = 1, V = null, busy = false;

  ctx.el.innerHTML = '';
  var root = el('div', 'display:flex;flex-direction:column;gap:8px;width:100%;height:100%;padding:10px;box-sizing:border-box;' +
    'background:#070a11;color:#d8dee9;font:13px system-ui,sans-serif;overflow:auto');
  ctx.el.appendChild(root);

  var hint = el('div', 'padding:7px 10px;border-radius:6px;background:#1b1706;border:1px solid #4a3c0c;color:#ffe08a;font-size:12.5px', root);
  hint.textContent = 'Run the flow to open the room.';

  // the drugs
  var card = box(root, 'Give a drug');
  var favs = el('div', 'display:flex;flex-wrap:wrap;gap:5px;margin-bottom:8px', card);
  var favB = {};
  FAVOURITES.forEach(function (n) { favB[n] = chip(favs, n, function () { pick(n); }); });

  var pickRow = el('div', 'display:flex;gap:8px;align-items:center;margin-bottom:9px', card);
  el('span', 'color:#9aa3b2;min-width:52px', pickRow).textContent = 'or any';
  var sel = el('select', 'flex:1;background:#121a2b;color:#e5e9f0;border:1px solid #2b3754;' +
    'border-radius:6px;padding:5px 7px;font:13px system-ui,sans-serif', pickRow);
  GROUPS.forEach(function (g) {
    var og = document.createElement('optgroup'); og.label = g[0];
    g[1].forEach(function (n) { var o = document.createElement('option'); o.value = o.textContent = n; og.appendChild(o); });
    sel.appendChild(og);
  });
  sel.onchange = function () { pick(sel.value); };

  var routeRow = el('div', 'display:flex;gap:6px;margin-bottom:9px', card);
  var rB = { bolus: tab(routeRow, 'One push (IV)', function () { setRoute('bolus'); }),
             infusion: tab(routeRow, 'Infusion', function () { setRoute('infusion'); }) };

  var amount = el('div', 'display:flex;align-items:center;gap:10px;margin-bottom:4px', card);
  btn(amount, '−', function () { step(-1); }, 'width:36px;font-size:16px');
  var big = el('div', 'flex:1;text-align:center;font:700 22px ui-monospace,Menlo,monospace;color:#7dffb0', amount);
  big.textContent = '—';
  btn(amount, '+', function () { step(1); }, 'width:36px;font-size:16px');
  var sub = el('div', 'text-align:center;color:#9aa3b2;font-size:11.5px;margin-bottom:9px', card);

  var go = btn(card, 'Administer', administer,
    'width:100%;padding:9px;background:#16643a;border-color:#1f8a50;color:#eaffee;font-weight:700;font-size:14px');
  var running = el('div', 'margin-top:8px;font-size:12px;color:#e879f9;min-height:16px', card);
  btn(card, 'Stop every infusion', function () { send({ control: 'stop infusions', value: true }); },
    'width:100%;margin-top:6px');

  // the airway
  var air = box(root, 'The airway and the oxygen');
  var preox = btn(air, 'Preoxygenate · fill the lungs with oxygen', function () {
    send({ control: 'preoxygenate', value: true });
  }, 'width:100%;margin-bottom:7px');
  var airRow = el('div', 'display:flex;gap:6px', air);
  var airB = {};
  AIRWAYS.forEach(function (a) { airB[a] = tab(airRow, a, function () { send({ control: 'airway', value: a }); }); });

  // the case
  var job = box(root, 'The case');
  var jobSel = line(job, 'Job', JOBS, function (v) { send({ control: 'scenario', value: v }); });
  line(job, 'Clock', CLOCKS, function (v) { send({ control: 'clock', value: v }); });
  var bottom = el('div', 'display:flex;gap:6px;margin-top:8px', job);
  btn(bottom, 'Demo: watch it played', function () { send({ control: 'demo', value: true }); }, 'flex:1');
  btn(bottom, 'Start again', function () { send({ control: 'start again', value: true }); }, 'flex:1');
  el('div', 'color:#6b7280;font-size:11px;margin-top:7px', job).textContent =
    'The patient is chosen on the block before the run.';

  function el(tag, css, parent) { var e = document.createElement(tag); if (css) e.style.cssText = css; if (parent) parent.appendChild(e); return e; }
  function box(parent, title) {
    var b = el('div', 'background:#0c111c;border:1px solid #1d2638;border-radius:8px;padding:10px', parent);
    el('div', 'font-weight:700;font-size:12px;letter-spacing:.03em;color:#9fb4ff;margin-bottom:8px', b).textContent = title.toUpperCase();
    return b;
  }
  function btn(parent, text, fn, css) {
    var b = el('button', 'background:#1c2742;color:#e5e9f0;border:1px solid #2f3d5e;border-radius:6px;padding:6px 10px;' +
      'cursor:pointer;font:13px system-ui,sans-serif;' + (css || ''), parent);
    b.textContent = text; b.onclick = fn; return b;
  }
  function chip(parent, text, fn) {
    var b = el('button', 'background:#121a2b;color:#c9d1dc;border:1px solid #2b3754;border-radius:999px;padding:4px 11px;' +
      'cursor:pointer;font:12px system-ui,sans-serif', parent);
    b.textContent = text; b.onclick = fn; return b;
  }
  function tab(parent, text, fn) {
    var b = el('button', 'flex:1;background:#121a2b;color:#c9d1dc;border:1px solid #2b3754;border-radius:6px;padding:6px 4px;' +
      'cursor:pointer;font:12px system-ui,sans-serif', parent);
    b.textContent = text; b.onclick = fn; return b;
  }
  function line(parent, label, options, fn) {
    var r = el('div', 'display:flex;gap:8px;align-items:center;margin-bottom:6px', parent);
    el('span', 'color:#9aa3b2;min-width:46px;font-size:12px', r).textContent = label;
    var s = el('select', 'flex:1;background:#121a2b;color:#e5e9f0;border:1px solid #2b3754;border-radius:6px;padding:5px 7px;' +
      'font:12.5px system-ui,sans-serif', r);
    options.forEach(function (o) { var e = document.createElement('option'); e.value = e.textContent = o; s.appendChild(e); });
    s.onchange = function () { fn(s.value); };
    return s;
  }
  function mark(b, on) {
    b.style.background = on ? '#24406e' : '#121a2b';
    b.style.borderColor = on ? '#4c7bd6' : '#2b3754';
    b.style.color = on ? '#eaf1ff' : '#c9d1dc';
  }

  function pick(n) { drug = n; if (sel.value !== n) sel.value = n; draw(); }
  function setRoute(r) { route = r; x = r === 'infusion' ? (infusionNow(drug) || 1) : 1; draw(); }
  function step(d) {
    var hi = route === 'infusion' ? 3 : 4;
    x = Math.max(0, Math.min(hi, Math.round((x + d * 0.25) * 100) / 100));
    draw();
  }
  function doses() { return (V && V.doses) || {}; }
  function infusionNow(n) { return (V && V.inf && V.inf[n]) || 0; }
  function send(msg) {
    if (busy) return;
    busy = true;
    var done = function () { busy = false; };
    try { ctx.send(CH, msg, TO ? { node: TO } : null).then(done, done); } catch (e) { done(); }
  }
  function administer() { send({ control: 'administer', value: { drug: drug, x: x, route: route } }); }

  function draw() {
    var d = doses()[drug] || {};
    FAVOURITES.forEach(function (n) { mark(favB[n], n === drug); });
    mark(rB.bolus, route === 'bolus'); mark(rB.infusion, route === 'infusion');
    if (route === 'infusion' && !d.rate) {
      big.textContent = 'no infusion';
      sub.textContent = drug + ' is given as a push, not as an infusion.';
      go.disabled = true; go.style.opacity = 0.45;
      return;
    }
    go.disabled = false; go.style.opacity = 1;
    var has = route === 'infusion' ? d.rate : d.bolus;
    big.textContent = has ? scaled(has, x) : (x + ' ×');
    sub.textContent = (has ? 'the usual ' + has + ' × ' + x.toFixed(2) : '× ' + x.toFixed(2)) +
      (V && V.weight ? ' · ' + V.weight + ' kg' : '');
    go.textContent = route === 'infusion' ? (x > 0 ? 'Run this infusion' : 'Stop this infusion') : 'Administer';
  }
  function scaled(text, k) {
    var m = /^([\d.]+)\s*(.*)$/.exec(text);
    if (!m) return text;
    var n = parseFloat(m[1]) * k;
    return (n >= 100 ? n.toFixed(0) : n >= 10 ? n.toFixed(1) : n.toFixed(2)) + ' ' + m[2];
  }

  function update() {
    var last = null;
    ctx.ports().forEach(function (port) {
      (ctx.input(port).stream || []).forEach(function (it) {
        if (it.channel === DATA && it.data && it.data.phase != null) last = it.data;
      });
    });
    if (!last) return;
    var first = !V;
    V = Object.assign({}, V, last);          // (the doses come every few seconds; keep the last we saw)
    hint.textContent = V.hint || '';
    if (first && V.scenario && JOBS.indexOf(V.scenario) >= 0) jobSel.value = V.scenario;
    if (first && V.drug) pick(V.drug);
    AIRWAYS.forEach(function (a) { mark(airB[a], a === V.airway); });
    mark(preox, false);
    var inf = V.inf || {};
    var names = Object.keys(inf);
    running.textContent = names.length
      ? 'Running: ' + names.map(function (n) { return n + ' ×' + inf[n]; }).join(' · ')
      : 'No infusion running.';
    draw();
  }

  ctx.on('update', update);
  update();
  draw();
}
