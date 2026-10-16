/* Single-reaction kinetics block shared by the ideal-reactor tools.
   Stoichiometry: A + bB -> rR + sS, conversion X of key reactant A.          */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N, X = ST.X;
  const K = ST.K = {};

  K.LAWS = [
    ['power', 'Power law', '-r_A = k\\,C_A^{a}C_B^{\\beta}'],
    ['rev', 'Reversible', '-r_A = k_1C_A^{a}C_B^{\\beta} - k_2C_R^{\\rho}C_S^{\\sigma}'],
    ['sat', 'Saturation', '-r_A = \\dfrac{k_1C_A}{1+k_2C_A}'],
    ['auto', 'Autocatalytic', '-r_A = k\\,C_AC_R'],
    ['lh', 'Langmuir–Hinshelwood', '-r_A = \\dfrac{k\\,C_AC_B}{(1+K_AC_A+K_BC_B)^2}'],
    ['inhib', 'Product inhibited', '-r_A = \\dfrac{k\\,C_A^{a}}{1+K_RC_R}'],
    ['custom', 'Your own expression', '-r_A = f(C_A, C_B, C_R, C_S, X_A)']
  ];

  /* parse "k = 0.5, K = 2" */
  K.parseParams = function (s) {
    const out = {};
    String(s || '').split(/[,;\n]/).map(t => t.trim()).filter(Boolean).forEach(kv => {
      const m = /^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([-+0-9.eE]+)$/.exec(kv);
      if (!m) throw new Error('Cannot read parameter "' + kv + '" (use name = value)');
      out[m[1]] = parseFloat(m[2]);
    });
    return out;
  };

  const usesB = v => (v.law === 'power' && +v.ob !== 0) || (v.law === 'rev' && (+v.ob !== 0)) || v.law === 'lh' || v.law === 'custom';
  const usesR = v => v.law === 'auto' || v.law === 'inhib' || v.law === 'rev' || v.law === 'custom';
  const usesS = v => (v.law === 'rev' && +v.os !== 0) || v.law === 'custom';
  K.usesB = usesB; K.usesR = usesR;

  /* input specs; o.prefix lets two blocks live in one tool; o.gas shows the expansion factor */
  K.inputs = function (o = {}) {
    const C = o.cu || 'mol/L', t = o.tu || 'min';
    const d = Object.assign({ law: 'power', k: 0.5, a: 1, ob: 0, k1: 1, k2: 0.2, orr: 1, os: 0, sk1: 1, sk2: 0.5,
      KA: 0.5, KB: 0.5, KR: 1, CA0: 1, eps: 0, b: 1, MB: 1, r: 1, MR: 0, s: 0, MS: 0,
      expr: 'k*CA^2/(1+K*CA)', cpar: 'k = 2, K = 0.5' }, o.d || {});
    const L = [
      { group: o.title || 'Rate law' },
      { k: 'law', t: 'select', label: 'Form', opts: K.LAWS.filter(l => !o.only || o.only.indexOf(l[0]) >= 0).map(l => [l[0], l[1]]), v: d.law },
      { k: 'k', label: '$k$', unit: 'rate constant', v: d.k, min: 1e-4, max: 1e3, scale: 'log', if: v => ['power', 'auto', 'lh', 'inhib'].indexOf(v.law) >= 0 },
      { k: 'a', label: 'Order in A, $a$', v: d.a, min: -1, max: 4, step: 0.05, if: v => ['power', 'rev', 'inhib'].indexOf(v.law) >= 0 },
      { k: 'ob', label: 'Order in B, $\\beta$', v: d.ob, min: 0, max: 3, step: 0.05, if: v => ['power', 'rev'].indexOf(v.law) >= 0 },
      { k: 'k1', label: '$k_1$', unit: 'forward', v: d.k1, min: 1e-4, max: 1e3, scale: 'log', if: v => v.law === 'rev' },
      { k: 'k2', label: '$k_2$', unit: 'reverse', v: d.k2, min: 0, max: 1e3, scale: 'log', if: v => v.law === 'rev' },
      { k: 'orr', label: 'Order in R, $\\rho$', v: d.orr, min: 0, max: 3, step: 0.05, if: v => v.law === 'rev' },
      { k: 'os', label: 'Order in S, $\\sigma$', v: d.os, min: 0, max: 3, step: 0.05, if: v => v.law === 'rev' },
      { k: 'sk1', label: '$k_1$', v: d.sk1, min: 1e-4, max: 1e3, scale: 'log', if: v => v.law === 'sat' },
      { k: 'sk2', label: '$k_2$', unit: C + '⁻¹', v: d.sk2, min: 0, max: 1e3, scale: 'log', if: v => v.law === 'sat' },
      { k: 'KA', label: '$K_A$', unit: C + '⁻¹', v: d.KA, min: 0, max: 100, if: v => v.law === 'lh' },
      { k: 'KB', label: '$K_B$', unit: C + '⁻¹', v: d.KB, min: 0, max: 100, if: v => v.law === 'lh' },
      { k: 'KR', label: '$K_R$', unit: C + '⁻¹', v: d.KR, min: 0, max: 100, if: v => v.law === 'inhib' },
      { k: 'expr', t: 'text', label: '$-r_A =$', v: d.expr, mono: true, hint: 'Use CA, CB, CR, CS, XA, CA0 and your parameters', if: v => v.law === 'custom' },
      { k: 'cpar', t: 'text', label: 'Parameters', v: d.cpar, mono: true, if: v => v.law === 'custom' },
      { group: 'Feed and stoichiometry' },
      { k: 'CA0', label: '$C_{A0}$', unit: C, v: d.CA0, min: 1e-4, max: 1e4, scale: 'log' }
    ];
    if (o.gas !== false) L.push({ k: 'eps', label: 'Expansion factor $\\varepsilon_A$', v: d.eps, min: -0.95, max: 4, step: 0.01, hint: '0 for liquids; (moles at X=1 − moles at X=0)/moles at X=0' });
    L.push(
      { k: 'b', label: 'Moles B per mole A, $b$', v: d.b, min: 0, max: 5, step: 0.1, if: usesB },
      { k: 'MB', label: '$M = C_{B0}/C_{A0}$', v: d.MB, min: 0, max: 50, step: 0.05, if: usesB },
      { k: 'r', label: 'Moles R formed per mole A, $r$', v: d.r, min: 0, max: 5, step: 0.1, if: usesR },
      { k: 'MR', label: '$C_{R0}/C_{A0}$', v: d.MR, min: 0, max: 50, step: 0.01, if: usesR },
      { k: 's', label: 'Moles S formed per mole A, $s$', v: d.s, min: 0, max: 5, step: 0.1, if: usesS },
      { k: 'MS', label: '$C_{S0}/C_{A0}$', v: d.MS, min: 0, max: 50, step: 0.01, if: usesS }
    );
    if (o.prefix) L.forEach(f => { if (f.k) f.k = o.prefix + f.k; });
    return L;
  };

  const texNum = v => {
    if (!isFinite(v)) return '\\infty';
    const s = N.fmt(v, 4, { plain: true });
    if (s.indexOf('e') > 0) { const [m, e] = s.split('e'); return m + '\\times10^{' + e + '}'; }
    return s.replace('\u2212', '-');
  };
  K.texNum = texNum;
  const pw = (sym, n) => (n === 0 ? '' : n === 1 ? sym : sym + '^{' + texNum(n) + '}');

  /* build the model from values v (with optional prefix) */
  K.make = function (v, prefix = '') {
    const g = k => v[prefix + k];
    const law = g('law'), CA0 = +g('CA0'), eps = +(g('eps') || 0);
    const b = usesB({ law, ob: g('ob') }) ? +g('b') : 0, MB = +(g('MB') || 0);
    const r = usesR({ law }) ? +g('r') : 0, MR = +(g('MR') || 0);
    const s = usesS({ law, os: g('os') }) ? +g('s') : 0, MS = +(g('MS') || 0);
    if (!(CA0 > 0)) throw new Error('C_A0 must be positive');
    const conc = X => {
      const den = 1 + eps * X;
      return { A: CA0 * (1 - X) / den, B: CA0 * (MB - b * X) / den, R: CA0 * (MR + r * X) / den, S: CA0 * (MS + s * X) / den };
    };
    let rate, tex;
    const P = (c, n) => (n === 0 ? 1 : Math.pow(Math.max(c, 0), n));
    switch (law) {
      case 'power': {
        const k = +g('k'), a = +g('a'), ob = +g('ob');
        rate = x => { const c = conc(x); return k * P(c.A, a) * P(c.B, ob); };
        tex = '-r_A = ' + texNum(k) + '\\,' + pw('C_A', a) + pw('C_B', ob);
        break;
      }
      case 'rev': {
        const k1 = +g('k1'), k2 = +g('k2'), a = +g('a'), ob = +g('ob'), orr = +g('orr'), os = +g('os');
        rate = x => { const c = conc(x); return k1 * P(c.A, a) * P(c.B, ob) - k2 * P(c.R, orr) * P(c.S, os); };
        tex = '-r_A = ' + texNum(k1) + '\\,' + pw('C_A', a) + pw('C_B', ob) + ' - ' + texNum(k2) + '\\,' + (pw('C_R', orr) + pw('C_S', os) || '1');
        break;
      }
      case 'sat': {
        const k1 = +g('sk1'), k2 = +g('sk2');
        rate = x => { const c = conc(x); return k1 * Math.max(c.A, 0) / (1 + k2 * Math.max(c.A, 0)); };
        tex = '-r_A = \\dfrac{' + texNum(k1) + '\\,C_A}{1+' + texNum(k2) + '\\,C_A}';
        break;
      }
      case 'auto': {
        const k = +g('k');
        rate = x => { const c = conc(x); return k * Math.max(c.A, 0) * Math.max(c.R, 0); };
        tex = '-r_A = ' + texNum(k) + '\\,C_AC_R';
        break;
      }
      case 'lh': {
        const k = +g('k'), KA = +g('KA'), KB = +g('KB');
        rate = x => { const c = conc(x), A = Math.max(c.A, 0), B = Math.max(c.B, 0); return k * A * B / Math.pow(1 + KA * A + KB * B, 2); };
        tex = '-r_A = \\dfrac{' + texNum(k) + '\\,C_AC_B}{(1+' + texNum(KA) + 'C_A+' + texNum(KB) + 'C_B)^2}';
        break;
      }
      case 'inhib': {
        const k = +g('k'), a = +g('a'), KR = +g('KR');
        rate = x => { const c = conc(x); return k * P(c.A, a) / (1 + KR * Math.max(c.R, 0)); };
        tex = '-r_A = \\dfrac{' + texNum(k) + '\\,' + pw('C_A', a) + '}{1+' + texNum(KR) + 'C_R}';
        break;
      }
      case 'custom': {
        const pars = K.parseParams(g('cpar'));
        const names = ['CA', 'CB', 'CR', 'CS', 'XA', 'CA0'].concat(Object.keys(pars));
        const fn = X.compile(g('expr'), names);
        const pv = Object.values(pars);
        rate = x => { const c = conc(x); return fn([Math.max(c.A, 0), Math.max(c.B, 0), Math.max(c.R, 0), Math.max(c.S, 0), x, CA0].concat(pv)); };
        tex = '-r_A = ' + String(g('expr')).replace(/\*/g, '\\,').replace(/CA0/g, 'C_{A0}').replace(/C([ABRS])/g, 'C_$1').replace(/XA/g, 'X_A');
        break;
      }
      default: throw new Error('Unknown rate form');
    }
    let Xmax = 1;
    if (b > 0) Xmax = Math.min(Xmax, MB / b);
    if (eps < 0) Xmax = Math.min(Xmax, -1 / eps * 0.999999);
    const m = { law, CA0, eps, b, MB, r, MR, s, MS, conc, rate, tex, Xmax,
      useB: usesB({ law, ob: g('ob') }), useR: usesR({ law }), useS: usesS({ law, os: g('os') }) };
    /* equilibrium (or stall) conversion: first X where rate drops to zero */
    m.Xe = Xmax;
    const r0 = rate(0);
    if (law === 'rev' || law === 'custom') {
      const roots = N.roots(x => rate(x), 0, Xmax * (1 - 1e-12), 400);
      const pos = roots.filter(x => x > 1e-12);
      if (pos.length) m.Xe = pos[0];
    }
    if (!(r0 > 0) && law !== 'auto') m.stall = true;
    return m;
  };

  /* ---------- performance of ideal reactors for one reaction ---------- */
  K.tauPFR = function (m, X1, X0 = 0) {
    if (X1 <= X0) return 0;
    if (X1 >= m.Xe - 1e-12) return Infinity;
    const v = m.CA0 * N.quad(x => 1 / m.rate(x), X0, X1, { rtol: 1e-10 });
    return v > 0 ? v : Infinity;
  };
  K.tBatch = function (m, X1, X0 = 0) {
    if (X1 <= X0) return 0;
    if (X1 >= m.Xe - 1e-12) return Infinity;
    return m.CA0 * N.quad(x => 1 / ((1 + m.eps * x) * m.rate(x)), X0, X1, { rtol: 1e-10 });
  };
  K.tauMFR = function (m, X1, X0 = 0) {
    if (X1 <= X0) return 0;
    if (X1 >= m.Xe - 1e-12) return Infinity;
    const r = m.rate(X1);
    return r > 0 ? m.CA0 * (X1 - X0) / r : Infinity;
  };
  /* conversion leaving a PFR of space time tau with entering conversion X0 */
  K.XPFR = function (m, tau, X0 = 0, o = {}) {
    if (tau <= 0) return { X: X0, prof: { t: [0], x: [X0] } };
    const cap = Math.min(m.Xe, m.Xmax) * (1 - 1e-12);
    const sol = N.ode((t, y) => [Math.max(0, m.rate(Math.min(y[0], cap))) / m.CA0], 0, [X0], tau,
      { rtol: 1e-10, atol: 1e-13, minPts: o.minPts || 12, stop: (t, y) => y[0] >= cap });
    const X = Math.min(sol.last()[0], cap);
    const n = o.n || 100, ts = N.linspace(0, tau, n);
    const tEnd = sol.t[sol.t.length - 1];
    return { X, prof: { t: ts, x: ts.map(t => (t <= tEnd ? Math.min(sol.at(t)[0], cap) : X)) } };
  };
  K.XBatch = function (m, tEnd, o = {}) {
    const cap = Math.min(m.Xe, m.Xmax) * (1 - 1e-12);
    const sol = N.ode((t, y) => [(1 + m.eps * y[0]) * Math.max(0, m.rate(Math.min(y[0], cap))) / m.CA0], 0, [0], tEnd,
      { rtol: 1e-10, atol: 1e-13, minPts: 60, stop: (t, y) => y[0] >= cap });
    const X = Math.min(sol.last()[0], cap), n = o.n || 100, ts = N.linspace(0, tEnd, n), tl = sol.t[sol.t.length - 1];
    return { X, prof: { t: ts, x: ts.map(t => (t <= tl ? Math.min(sol.at(t)[0], cap) : X)) } };
  };
  /* all steady states of a single MFR with entering conversion X0 */
  K.XMFRall = function (m, tau, X0 = 0) {
    if (tau <= 0) return [X0];
    const top = Math.min(m.Xe, m.Xmax);
    const f = x => m.CA0 * (x - X0) - tau * m.rate(x);
    const roots = N.roots(f, X0, top, 240, { tol: 1e-14 }).filter(x => x >= X0 - 1e-12);
    if (!roots.length) return [top];
    return roots;
  };
  K.XMFR = function (m, tau, X0 = 0) {
    const r = K.XMFRall(m, tau, X0);
    return r[r.length - 1];
  };

  /* Levenspiel curve 1/(-rA) vs X sampled for plotting */
  K.curve = function (m, Xhi, n = 200) {
    const top = Math.min(Xhi, m.Xe * 0.999);
    const xs = N.linspace(0, top, n);
    return { x: xs, y: xs.map(x => { const r = m.rate(x); return r > 0 ? 1 / r : NaN; }) };
  };

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
