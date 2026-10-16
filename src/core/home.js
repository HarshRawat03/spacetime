/* Home page: the Levenspiel plot as a live instrument, then the full equipment list. */
(function (G) {
  'use strict';
  const ST = G.ST;
  const N = ST.N, K = ST.K, D = ST.D, UI = ST.UI;
  const $ = (s, r = document) => r.querySelector(s);

  const KIN = {
    n: { label: 'nth order, −r = kCₐⁿ', p: 'Order n', min: -0.5, max: 3, step: 0.05, v: 1 },
    auto: { label: 'Autocatalytic, A + R → 2R', p: 'Seed C_R0/C_A0', min: 0.005, max: 0.5, step: 0.005, v: 0.05 },
    rev: { label: 'Reversible first order, A ⇌ R', p: 'Equilibrium X_Ae', min: 0.6, max: 0.995, step: 0.005, v: 0.9 },
    sat: { label: 'Saturation, kC/(1+KC)', p: 'K·C_A0', min: 0, max: 20, step: 0.1, v: 4 },
    inh: { label: 'Substrate inhibited, kC/(1+KC)²', p: 'K·C_A0', min: 0.2, max: 12, step: 0.1, v: 4 }
  };
  const ARR = [
    ['pfr', 'Plug flow'], ['mfr', 'Mixed flow'], ['m2', 'Two mixed, equal'], ['m4', 'Four mixed, equal'],
    ['mp', 'Mixed then plug, best split'], ['rec', 'Plug with recycle, best R']
  ];
  const st = { kin: 'n', p: 1, X: 0.8, arr: 'pfr' };
  try { const s = JSON.parse(localStorage.getItem('st:home') || 'null'); if (s && KIN[s.kin]) Object.assign(st, s); } catch (e) { /* no storage */ }

  function model() {
    const base = { CA0: 1, eps: 0, b: 1, MB: 1, r: 1, MR: 0, s: 0, MS: 0, a: 1, ob: 0, k: 1, orr: 1, os: 0 };
    let v;
    if (st.kin === 'n') v = Object.assign(base, { law: 'power', a: st.p });
    else if (st.kin === 'auto') v = Object.assign(base, { law: 'auto', MR: st.p });
    else if (st.kin === 'rev') v = Object.assign(base, { law: 'rev', k1: 1, k2: (1 - st.p) / st.p });
    else if (st.kin === 'sat') v = Object.assign(base, { law: 'sat', sk1: 1 + st.p, sk2: st.p });
    else v = Object.assign(base, { law: 'custom', expr: 'k*CA/(1+K*CA)^2', cpar: 'k = ' + (1 + st.p) ** 2 + ', K = ' + st.p });
    return K.make(v);
  }
  function tanksTau(m, Xf, n) {
    const f = t => { let X = 0; for (let i = 0; i < n; i++) X = K.XMFR(m, t, X); return X - Xf; };
    let hi = K.tauMFR(m, Xf) || 10;
    if (!isFinite(hi)) return { tau: Infinity };
    for (let i = 0; i < 60 && f(hi) < 0; i++) hi *= 2;
    const t = N.brent(f, 1e-9, hi, { tol: 1e-12 });
    const xs = [0];
    for (let i = 0; i < n; i++) xs.push(K.XMFR(m, t, xs[i]));
    return { tau: t * n, each: t, xs };
  }
  function compute() {
    const m = model();
    const top = Math.min(m.Xe, m.Xmax);
    const Xf = Math.min(st.X, top * 0.995);
    const r = { m, Xf, top };
    r.pfr = K.tauPFR(m, Xf);
    r.mfr = K.tauMFR(m, Xf);
    r.m2 = tanksTau(m, Xf, 2);
    r.m4 = tanksTau(m, Xf, 4);
    const mp = N.fminScan(x1 => K.tauMFR(m, x1) + K.tauPFR(m, Xf, x1), 1e-6, Xf * 0.999, 80);
    r.mp = { tau: mp.f, X1: mp.x };
    const recTau = lr => {
      const R = Math.pow(10, lr), X1 = R / (R + 1) * Xf;
      return (R + 1) * K.tauPFR(m, Xf, X1);
    };
    const rc = N.fminScan(recTau, -3, 2.5, 80);
    r.rec = { tau: rc.f, R: Math.pow(10, rc.x) };
    if (r.rec.R < 0.0015 && r.pfr <= r.rec.tau) r.rec = { tau: r.pfr, R: 0 };
    return r;
  }

  function plotSpec(r) {
    const c = K.curve(r.m, Math.min(0.995, r.top * 0.999), 240);
    const y = c.y.map(v => v);
    const Xf = r.Xf;
    const finite = y.filter(isFinite);
    const ymaxData = Math.max(...finite.filter((v, i) => c.x[i] <= Xf + 1e-9));
    const ycap = Math.min(Math.max(...finite), ymaxData * 1.6 + 0.2);
    const spec = {
      h: 300,
      x: { label: 'X_A', min: 0, max: 1 },
      y: { label: 'C_{A0}/(−r_A)', min: 0, max: ycap },
      series: [{ name: 'C_{A0}/(−r_A)', x: c.x, y, color: 1, w: 2.4, hold: false }],
      marks: [],
      legend: false,
      drag: { x: Xf, min: 0.02, max: Math.min(0.99, r.top * 0.99), onDrag: x => { st.X = x; save(); draw(); } }
    };
    const inv = x => 1 / r.m.rate(x);
    const fill = (a, b, color) => {
      const xs = N.linspace(a, b, 80);
      spec.series.push({ x: xs, y: xs.map(inv), type: 'area', fill: 0.22, color, hover: false, legend: false, noRange: true });
    };
    const rect = (a, b, h, color) => spec.marks.push({ t: 'rect', x0: a, x1: b, y0: 0, y1: h, color, fill: 0.2, noRange: true });
    if (st.arr === 'pfr') fill(0, Xf, 3);
    else if (st.arr === 'mfr') rect(0, Xf, inv(Xf), 2);
    else if (st.arr === 'm2' || st.arr === 'm4') { const xs = (st.arr === 'm2' ? r.m2 : r.m4).xs || []; for (let i = 1; i < xs.length; i++) rect(xs[i - 1], xs[i], inv(xs[i]), i % 2 ? 2 : 4); }
    else if (st.arr === 'mp') { rect(0, r.mp.X1, inv(r.mp.X1), 2); fill(r.mp.X1, Xf, 3); }
    else if (st.arr === 'rec') {
      const R = r.rec.R, X1 = R / (R + 1) * Xf;
      rect(0, Xf, r.rec.tau / Xf, 6);
      if (R > 0) spec.marks.push({ t: 'v', x: X1, label: 'X_{A1}', color: 6, noRange: true });
    }
    spec.marks.push({ t: 'v', x: Xf, label: 'X_{Af} ' + N.fmt(Xf, 3), color: 8, dash: false, noRange: true, dy: -2 });
    return spec;
  }
  function diagram(r) {
    const m = r.m, Xf = r.Xf;
    const prof = (a, b) => N.linspace(0, 1, 16).map(s => {
      const t = s * K.tauPFR(m, b, a);
      return isFinite(t) ? K.XPFR(m, t, a).X : a + s * (b - a);
    });
    let units = [];
    if (st.arr === 'pfr') units = [{ type: 'pfr', tag: 'R-101', X: Xf, prof: prof(0, Xf) }];
    else if (st.arr === 'mfr') units = [{ type: 'mfr', tag: 'R-101', X: Xf }];
    else if (st.arr === 'm2' || st.arr === 'm4') { const xs = (st.arr === 'm2' ? r.m2 : r.m4).xs || []; units = xs.slice(1).map((x, i) => ({ type: 'mfr', tag: 'R-10' + (i + 1), X: x })); }
    else if (st.arr === 'mp') units = [{ type: 'mfr', tag: 'R-101', X: r.mp.X1 }, { type: 'pfr', tag: 'R-102', X: Xf, prof: prof(r.mp.X1, Xf) }];
    else if (st.arr === 'rec') {
      const R = r.rec.R, X1 = R / (R + 1) * Xf;
      let s = D.pipe([[6, 50], [70, 50]], { label: 'feed' });
      s += D.node(74, 50);
      s += D.pipe([[77, 50], [104, 50]]);
      s += D.tube(106, 39, 190, 22, { prof: prof(X1, Xf), tag: 'R-101' });
      s += D.pipe([[298, 50], [356, 50]]);
      s += D.node(360, 50);
      s += D.pipe([[363, 50], [440, 50]], { label: 'product' , at: 1, dx: -50});
      s += D.pipe([[360, 54], [360, 104], [74, 104], [74, 56]], { cls: 'rec' });
      s += D.text(217, 118, 'recycle R = ' + N.fmt(R, 3), { cls: 'pid-r' });
      s += D.text(330, 42, 'X ' + N.fmt(Xf, 3), { cls: 'pid-r' });
      return D.svg(450, 128, s, 'plug flow reactor with recycle');
    }
    return D.train(units, { h: 118, aria: 'reactor arrangement' }).svg;
  }

  let host = null;
  const save = () => { try { localStorage.setItem('st:home', JSON.stringify(st)); } catch (e) { /* ignore */ } };
  function draw() {
    if (!host) return;
    let r;
    try { r = compute(); } catch (e) { $('.inst-out', host).innerHTML = '<div class="note bad">' + UI.esc(e.message) + '</div>'; return; }
    const vals = [['pfr', r.pfr], ['mfr', r.mfr], ['m2', r.m2.tau], ['m4', r.m4.tau], ['mp', r.mp.tau], ['rec', r.rec.tau]];
    const best = vals.reduce((a, b) => (b[1] < a[1] ? b : a));
    const names = { pfr: 'Plug', mfr: 'Mixed', m2: '2 mixed', m4: '4 mixed', mp: 'Mixed→plug', rec: 'Recycle' };
    const kp = vals.map(([k, v]) => '<div class="kpi' + (k === st.arr ? ' main' : '') + '"><div class="kl">' + names[k] + '</div><div class="kv">' + N.fmt(v, 3, { html: true }) + '</div><div class="ks">' + (k === best[0] ? 'smallest' : (v / best[1]).toFixed(2) + '× smallest') + '</div></div>').join('');
    const extra = st.arr === 'mp' ? 'Switch to plug flow at X_A = ' + N.fmt(r.mp.X1, 3) + '.' : st.arr === 'rec' ? 'Best recycle ratio R = ' + N.fmt(r.rec.R, 3) + '.' : st.arr === 'm2' || st.arr === 'm4' ? 'Each tank: kτ = ' + N.fmt((st.arr === 'm2' ? r.m2 : r.m4).each, 3) + '.' : '';
    $('.inst-out', host).innerHTML = '<figure class="diag">' + diagram(r) + '</figure><div class="plot-host" id="hp"></div><div class="kpis" aria-label="Dimensionless space-time kτ for each arrangement">' + kp + '</div><p class="cap">Numbers are k·τ for C<sub>A0</sub> = 1: the shaded area under the curve. Drag the handle to change the target conversion. ' + UI.lab(extra) + '</p>';
    ST.Plot.render($('#hp', host), plotSpec(r));
  }

  ST.Home = {
    html() {
      const nTools = ST.TOOLS.length, nEx = ST.TOOLS.reduce((s, t) => s + (t.presets || []).filter(p => p.book).length, 0);
      const chs = new Set(ST.TOOLS.map(t => t.ch)), app = chs.has(31);
      const nCh = chs.size - (app ? 1 : 0);
      const kinOpts = Object.keys(KIN).map(k => '<option value="' + k + '"' + (st.kin === k ? ' selected' : '') + '>' + KIN[k].label + '</option>').join('');
      const arrOpts = ARR.map(a => '<option value="' + a[0] + '"' + (st.arr === a[0] ? ' selected' : '') + '>' + a[1] + '</option>').join('');
      const kd = KIN[st.kin];
      return '<section class="hero"><div><span class="tag tag-l">CRE · 3RD EDITION</span><h1>Every chapter of Levenspiel, as a <em>working</em> instrument.</h1>' +
        '<p>Calculators for the whole of <i>Chemical Reaction Engineering</i>: rate data analysis, ideal and non-ideal reactors, temperature effects, solid catalysts, fluid–fluid and fluid–particle systems, and fermentation.</p>' +
        '<p>Load a worked example from the book and the tool checks its answer against the printed one. Change any input and every plot, reactor sketch and table follows.</p>' +
        '<div class="facts"><span><b>' + nTools + '</b>tools</span><span><b>' + nCh + '</b>chapters' + (app ? ' + appendix' : '') + '</span><span><b>' + nEx + '</b>book answers checked</span></div>' +
        '<div class="go"><a class="btn" href="#part-H">Browse the chapters</a><a class="btn" href="#mfr">Start with a reactor</a></div></div>' +
        '<div class="inst" id="inst"><div class="inst-h"><b>Plug or mixed? The design equation, drawn.</b><span class="tag">H-000</span></div><div class="inst-b"><div class="inst-ctl">' +
        '<label for="h-kin">Kinetics<select id="h-kin">' + kinOpts + '</select></label>' +
        '<label for="h-p"><span class="lh"><span id="h-pl">' + UI.lab(kd.p) + '</span><output id="h-po">' + N.fmt(st.p, 3) + '</output></span><input id="h-p" class="rng" type="range" min="' + kd.min + '" max="' + kd.max + '" step="' + kd.step + '" value="' + st.p + '"></label>' +
        '<label for="h-arr">Arrangement<select id="h-arr">' + arrOpts + '</select></label>' +
        '</div><div class="inst-out"></div></div></div></section>' + ST.App.plantHTML();
    },
    mount(main) {
      host = $('#inst', main);
      const kin = $('#h-kin', host), p = $('#h-p', host), arr = $('#h-arr', host);
      kin.addEventListener('change', () => {
        st.kin = kin.value;
        const kd = KIN[st.kin];
        st.p = kd.v;
        p.min = kd.min; p.max = kd.max; p.step = kd.step; p.value = kd.v;
        $('#h-pl', host).innerHTML = UI.lab(kd.p);
        $('#h-po', host).textContent = N.fmt(st.p, 3);
        if (st.kin === 'auto' && st.arr === 'pfr') { st.arr = 'mp'; arr.value = 'mp'; }
        save(); draw();
      });
      p.addEventListener('input', () => { st.p = +p.value; $('#h-po', host).textContent = N.fmt(st.p, 3); save(); draw(); });
      arr.addEventListener('change', () => { st.arr = arr.value; save(); draw(); });
      draw();
    }
  };
})(window);
