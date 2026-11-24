// ── EsaSim block ──────────────────────────────────────────────────────
// Bedside monitor (ws-bedside-monitor) · version 2 · view, JavaScript: function render(ctx)
// A bedside patient monitor: ECG, pulse oximetry and breathing waveforms, vitals (HR, SpO2, resp, BP, temperature, QTc, pain, sedation), alarms and the drugs on board
// Workspace: Gianluca De Novi's workspace · Author: Gianluca De Novi
// Saved: 2026-11-24
//
// Inputs   ctx.input("id") gives { files, stream }; ctx.rows("id"), ctx.json("id"), ctx.text("id")
//   stream      Stream                    Live telemetry · The vitals
// Settings ctx.params.key
//   channel     text      Vitals channel = "vitals"
// EsaSim rewrites these lines from the block's definition each time it is saved: change the definition, not them.
// ── end of EsaSim block header ────────────────────────────────────────

// Bedside monitor: the patient's vitals as a bedside monitor shows them (ECG, pulse oximetry and breathing waveforms,
// the numbers, the alarms, an optional beep). Reads the telemetry channel "vitals". The patient and the drugs are
// chosen and given from the 3D view's controls.
function render(ctx) {
  var p = ctx.params || {};
  var CH = p.channel || 'vitals';
  var V = null, sound = false, audio = null;

  ctx.el.innerHTML = '';
  var root = el('div', 'display:flex;flex-direction:column;width:100%;height:100%;min-height:380px;background:#05070b;color:#d8dee9;font:13px system-ui,sans-serif;overflow:hidden');
  ctx.el.appendChild(root);

  var head = el('div', 'display:flex;gap:10px;align-items:center;padding:6px 12px;background:#0d1220;border-bottom:1px solid #1f2a40', root);
  var status = el('div', 'flex:1;color:#c9d1dc;font-size:12px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis', head);
  status.textContent = 'Run the flow to open the bedside.';
  var bSound = btn(head, '🔇 Beep', function () {
    sound = !sound;
    bSound.textContent = (sound ? '🔊' : '🔇') + ' Beep';
    if (sound && !audio) { try { audio = new (window.AudioContext || window.webkitAudioContext)(); } catch (e) { audio = null; } }
  });

  var mon = el('div', 'flex:1;display:flex;min-height:0;position:relative', root);
  var waves = el('div', 'flex:1;position:relative;min-width:0', mon);
  var canvas = document.createElement('canvas');
  canvas.style.cssText = 'position:absolute;inset:0;width:100%;height:100%';
  waves.appendChild(canvas);
  var alarm = el('div', 'position:absolute;top:6px;left:8px;right:8px;padding:4px 10px;border-radius:4px;font-weight:700;letter-spacing:.04em;display:none;z-index:2', waves);
  var side = el('div', 'width:min(46%,320px);display:grid;grid-template-columns:1fr 1fr;grid-auto-rows:1fr;gap:1px;background:#141b2b;border-left:1px solid #141b2b', mon);
  var N = {
    hr: box('HR', '#39e75f', 'bpm', 2),
    spo2: box('SpO₂', '#22d3ee', '%', 1), rr: box('RESP', '#facc15', '/min', 1),
    bp: box('NIBP', '#f87171', 'mmHg', 2),
    temp: box('TEMP', '#e5e7eb', '°C', 1), qtc: box('QTc', '#fb923c', 'ms', 1),
    pain: box('PAIN', '#f9a8d4', '/10', 1), sed: box('SEDATION', '#a5b4fc', 'Ramsay', 1)
  };
  var foot = el('div', 'padding:5px 12px;background:#080b13;border-top:1px solid #141b2b;font-size:12px;color:#9aa3b2;display:flex;flex-direction:column;gap:2px', root);
  var fClin = el('div', '', foot), fDrugs = el('div', 'color:#e879f9', foot);

  function el(tag, css, parent) { var e = document.createElement(tag); if (css) e.style.cssText = css; if (parent) parent.appendChild(e); return e; }
  function btn(parent, text, fn) {
    var b = el('button', 'background:#1c2742;color:#e5e9f0;border:1px solid #2f3d5e;border-radius:5px;padding:2px 9px;cursor:pointer;font:12px system-ui,sans-serif', parent);
    b.textContent = text; b.onclick = fn; return b;
  }
  function box(label, color, unit, span) {
    var b = el('div', 'background:#05070b;padding:5px 10px;display:flex;flex-direction:column;justify-content:center;min-height:0;grid-column:span ' + span, side);
    var top = el('div', 'display:flex;justify-content:space-between;font-size:11px;font-weight:600;color:' + color, b);
    el('span', '', top).textContent = label;
    var u = el('span', 'opacity:.75', top); u.textContent = unit;
    var v = el('div', 'font:700 28px ui-monospace,Menlo,monospace;color:' + color + ';line-height:1.1', b); v.textContent = '--';
    var s = el('div', 'font-size:11px;color:' + color + ';opacity:.85;white-space:nowrap;overflow:hidden;text-overflow:ellipsis', b);
    return { b: b, v: v, s: s, u: u };
  }

  // ── the latest vitals ──
  function update() {
    var last = null;
    ctx.ports().forEach(function (port) {
      (ctx.input(port).stream || []).forEach(function (it) { if (it.channel === CH && it.data && it.data.hr != null) last = it.data; });
    });
    if (!last) return;
    V = last;
    status.textContent = [V.who, V.day].filter(Boolean).join(' · ');
    N.hr.v.textContent = V.hr; N.hr.s.textContent = V.rhythm || '';
    N.spo2.v.textContent = V.spo2;
    N.rr.v.textContent = V.rr; N.rr.s.textContent = V.airway && V.airway !== 'clear' ? V.airway : '';
    N.bp.v.textContent = V.sys + '/' + V.dia; N.bp.s.textContent = '(' + V.map + ')';
    N.temp.v.textContent = (+V.temp).toFixed(1);
    N.qtc.v.textContent = V.qtc;
    N.pain.v.textContent = V.pain != null ? Math.round(+V.pain) : '--';
    N.sed.v.textContent = V.sedation != null ? V.sedation : '--'; N.sed.s.textContent = V.sedation_txt || '';
    N.qtc.v.style.color = +V.qtc >= 500 ? '#f87171' : '#fb923c';
    fClin.textContent = 'Nausea: ' + (V.nausea || '–') + ' · Airways: ' + (V.airway || '–') + ' · Other effects: ' + (V.effects || 'none');
    fDrugs.textContent = (V.drugs && V.drugs.length) ? 'On board: ' + V.drugs.join(' · ') : 'No drug given yet';
    var red = (V.alarms || []).filter(function (a) { return /LOW|ARRHYTHMIA|QTc/.test(a) || a === 'HR HIGH'; });
    if ((V.alarms || []).length) {
      alarm.style.display = '';
      alarm.textContent = '⚠ ' + V.alarms.join(' · ');
      alarm.style.background = red.length ? '#b91c1c' : '#a16207';
      alarm.style.color = '#fff';
    } else alarm.style.display = 'none';
  }

  // ── waveforms: drawn in real time at the patient's current rates ──
  var g = canvas.getContext('2d'), W = 0, H = 0, dpr = 1, x = 0, last = performance.now();
  var beatT = 0, beatLen = 0.83, beatKind = 'n', breathPh = 0;
  var prev = { ecg: null, pl: null, rs: null };
  function resize() {
    dpr = Math.min(2, window.devicePixelRatio || 1);
    W = Math.max(50, waves.clientWidth); H = Math.max(50, waves.clientHeight);
    canvas.width = W * dpr; canvas.height = H * dpr; g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.fillStyle = '#05070b'; g.fillRect(0, 0, W, H); x = 0; prev = { ecg: null, pl: null, rs: null };
    labels();
  }
  function labels() {
    g.font = '600 11px system-ui,sans-serif';
    g.fillStyle = '#39e75f'; g.fillText('II', 8, H * 0.08 + 22);
    g.fillStyle = '#22d3ee'; g.fillText('Pleth', 8, H * 0.52 + 12);
    g.fillStyle = '#facc15'; g.fillText('Resp', 8, H * 0.76 + 12);
  }
  function gauss(t, m, s) { return Math.exp(-0.5 * ((t - m) / s) * ((t - m) / s)); }
  function ecg(t, kind, f, q) {
    if (kind === 'v') return 1.25 * gauss(t, 0.06 * f + 0.1, 0.035) - 0.35 * gauss(t, 0.2 * f + 0.1, 0.03) - 0.45 * gauss(t, 0.42 * f, 0.06);
    var v = 0.14 * gauss(t, 0.09 * f, 0.022);
    if (kind === 'd') return v;
    return v - 0.1 * gauss(t, 0.19 * f, 0.008) + 1.0 * gauss(t, 0.21 * f, 0.009) - 0.28 * gauss(t, 0.235 * f, 0.009) + 0.3 * gauss(t, 0.42 * f * q, 0.04);
  }
  function pleth(t, kind, f) {
    if (kind === 'd') return 0;
    var a = kind === 'v' ? 0.5 : 1;
    return a * (gauss(t, 0.3 * f, 0.07) + 0.35 * gauss(t, 0.52 * f, 0.06));
  }
  function nextBeat() {
    var hr = V ? Math.max(25, +V.hr) : 72;
    beatLen = 60 / hr * (1 + (Math.random() - 0.5) * 0.04);
    var r = Math.random(), ect = V ? +V.ectopy || 0 : 0, drop = V ? +V.dropped || 0 : 0;
    beatKind = r < ect ? 'v' : (r < ect + drop ? 'd' : 'n');
    if (beatKind === 'v') beatLen *= 1.35;
    if (beatKind !== 'd' && sound && audio) beep();
  }
  function beep() {
    var sp = V ? +V.spo2 : 98, o = audio.createOscillator(), gn = audio.createGain();
    o.frequency.value = 440 + (sp - 80) * 22; o.type = 'sine';
    gn.gain.setValueAtTime(0.0001, audio.currentTime); gn.gain.exponentialRampToValueAtTime(0.12, audio.currentTime + 0.01);
    gn.gain.exponentialRampToValueAtTime(0.0001, audio.currentTime + 0.12);
    o.connect(gn); gn.connect(audio.destination); o.start(); o.stop(audio.currentTime + 0.14);
  }
  function frame(now) {
    var dt = Math.min(0.1, (now - last) / 1000); last = now;
    if (!W) { requestAnimationFrame(frame); return; }
    var pxs = W / 6;                                    // (6 s across)
    var steps = Math.max(1, Math.round(dt * pxs));
    var hdt = dt / steps;
    var rr = V ? Math.max(3, +V.rr) : 14, deep = V && V.sedation ? Math.max(0.35, 1 - 0.12 * (V.sedation - 1)) : 1;
    var q = V && V.qtc ? Math.max(0.9, Math.min(1.5, +V.qtc / 420)) : 1;    // (a long QT: the T wave later)
    for (var k = 0; k < steps; k++) {
      beatT += hdt;
      if (beatT >= beatLen) { beatT -= beatLen; nextBeat(); }
      breathPh += hdt * rr / 60 * 2 * Math.PI;
      var f = Math.min(1, Math.sqrt(beatLen / 0.9));
      var ye = H * 0.30 - ecg(beatT, beatKind, f, q) * H * 0.17;
      var yp = H * 0.68 - pleth(beatT, beatKind, f) * H * 0.11 * (V ? Math.max(0.4, (+V.spo2 - 70) / 28) : 1);
      var yr = H * 0.90 - Math.sin(breathPh) * H * 0.05 * deep;
      g.fillStyle = '#05070b'; g.fillRect(x + 1, 0, 14, H);
      line('#39e75f', prev.ecg, ye); line('#22d3ee', prev.pl, yp); line('#facc15', prev.rs, yr);
      prev = { ecg: ye, pl: yp, rs: yr };
      x += 1;
      if (x >= W) { x = 0; prev = { ecg: null, pl: null, rs: null }; labels(); }
    }
    requestAnimationFrame(frame);
  }
  function line(c, y0, y1) {
    if (y0 == null) return;
    g.strokeStyle = c; g.lineWidth = 1.6; g.beginPath(); g.moveTo(x, y0); g.lineTo(x + 1, y1); g.stroke();
  }

  ctx.on('resize', resize);
  ctx.on('update', update);
  setTimeout(resize, 30);
  update();
  requestAnimationFrame(frame);
}
