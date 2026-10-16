/* Helpers shared by the tool definitions: units, flow basis, standard plots and tables. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N, K = ST.K;
  const KIT = ST.KIT = {};

  /* units of an nth-order rate constant, e.g. (L/mol)^0.5/min */
  KIT.kUnit = function (n, cu = 'mol/L', tu = 'min') {
    const [num, den] = cu.split('/');
    const e = 1 - n;
    if (Math.abs(e) < 1e-9) return '1/' + tu;
    const pw = Math.abs(e), ps = Math.abs(pw - 1) < 1e-9 ? '' : '^{' + N.fmt(pw, 3) + '}';
    if (e > 0) return '(' + cu + ')' + ps + '/' + tu;
    if (!den) return '(1/' + cu + ')' + ps + '/' + tu;   /* a pressure standing in for concentration */
    return '(' + den + '/' + num + ')' + ps + '/' + tu;
  };
  KIT.lawOrder = v => (v.law === 'power' || v.law === 'rev' ? (+v.a || 0) + (+v.ob || 0) : v.law === 'auto' || v.law === 'lh' ? 2 : 1);
  KIT.kUnitFor = v => KIT.kUnit(KIT.lawOrder(v), v._cu, v._tu);

  /* attach dynamic units to the K block fields */
  KIT.withUnits = function (fields) {
    fields.forEach(f => {
      if (f.k === 'k' || f.k === 'k1') f.unit = v => KIT.kUnitFor(v);
      if (f.k === 'k2' && f.label && f.label.indexOf('k_2') >= 0) f.unit = v => KIT.kUnit((+v.orr || 0) + (+v.os || 0), v._cu, v._tu);
      if (f.k === 'sk1') f.unit = v => '1/' + (v._tu || 'min');
      if (f.k === 'sk2' || f.k === 'KA' || f.k === 'KB' || f.k === 'KR') f.unit = v => '(' + (v._cu || 'mol/L') + ')⁻¹';
      if (f.k === 'CA0') f.unit = '[C]';
    });
    return fields;
  };

  /* throughput inputs */
  KIT.flowInputs = function (d = {}) {
    return [
      { group: 'Throughput (optional)' },
      { k: 'basis', t: 'select', label: 'Size the reactor from', opts: [['none', 'space-time only'], ['v0', 'volumetric feed rate v₀'], ['FA0', 'molar feed rate F_A0']], v: d.basis || 'none' },
      { k: 'v0', label: 'Feed rate $v_0$', unit: v => (v._vu || 'L') + '/' + (v._tu || 'min'), v: d.v0 || 1, min: 1e-6, max: 1e6, scale: 'log', if: v => v.basis === 'v0' },
      { k: 'FA0', label: 'Feed rate $F_{A0}$', unit: v => 'mol/' + (v._tu || 'min'), v: d.FA0 || 1, min: 1e-6, max: 1e6, scale: 'log', if: v => v.basis === 'FA0' }
    ];
  };
  /* volume units implied by the concentration unit */
  KIT.volUnit = cu => (/L$/.test(cu) ? 'L' : /m³$/.test(cu) ? 'm³' : /cm³$/.test(cu) ? 'cm³' : 'L');
  KIT.v0of = (v, CA0) => (v.basis === 'v0' ? +v.v0 : v.basis === 'FA0' ? +v.FA0 / CA0 : NaN);

  /* Levenspiel plot of CA0/(-rA) against X with optional shading */
  KIT.levPlot = function (m, Xf, o = {}) {
    const top = Math.min(m.Xe, m.Xmax);
    const hi = Math.min(o.xmax || Math.max(Xf * 1.15, Math.min(0.99, top * 0.98)), top * 0.995);
    const xs = N.linspace(0, hi, 220);
    const ys = xs.map(x => { const r = m.rate(x); return r > 0 ? m.CA0 / r : NaN; });
    const atX = xs.filter(x => x <= Xf).map(x => m.CA0 / m.rate(x)).filter(isFinite);
    const ymax = o.ymax || Math.max(...atX, m.CA0 / m.rate(Xf)) * 1.35;
    const spec = {
      title: o.title || 'Levenspiel plot',
      x: { label: 'X_A', min: 0, max: hi },
      y: { label: 'C_{A0}/(−r_A)', unit: o.tu || '', min: 0, max: isFinite(ymax) ? ymax : undefined },
      series: [{ name: 'C_{A0}/(−r_A)', x: xs, y: ys, color: 1 }],
      marks: [], legend: false
    };
    const inv = x => m.CA0 / m.rate(x);
    (o.shade || []).forEach(s => {
      if (s.type === 'pfr') {
        const a = s.a || 0, b = s.b !== undefined ? s.b : Xf, sx = N.linspace(a, b, 90);
        spec.series.push({ x: sx, y: sx.map(inv), type: 'area', fill: 0.22, color: s.color || 3, hover: false, legend: false, noRange: true });
      } else if (s.type === 'mfr') {
        const a = s.a || 0, b = s.b !== undefined ? s.b : Xf;
        spec.marks.push({ t: 'rect', x0: a, x1: b, y0: 0, y1: s.h !== undefined ? s.h : inv(b), color: s.color || 2, fill: 0.2, noRange: true, label: s.label });
      }
    });
    if (o.markX !== false) spec.marks.push({ t: 'v', x: Xf, label: 'X_{Af}', color: 8, noRange: true });
    if (m.Xe < m.Xmax - 1e-9 && m.Xe < 1) spec.marks.push({ t: 'v', x: m.Xe, label: 'X_{Ae}', color: 2, noRange: true });
    return spec;
  };

  KIT.concTable = function (m, rows, title) {
    const usesB = m.useB, usesR = m.useR, usesS = m.useS && (m.s > 0 || m.MS > 0);
    const cols = ['', 'X_A', 'C_A'].concat(usesB ? ['C_B'] : [], usesR ? ['C_R'] : [], usesS ? ['C_S'] : [], ['−r_A']);
    return {
      title: title || 'Compositions', cols,
      rows: rows.map(([lab, X]) => {
        const c = m.conc(X);
        return [lab, X, c.A].concat(usesB ? [c.B] : [], usesR ? [c.R] : [], usesS ? [c.S] : [], [m.rate(X)]);
      })
    };
  };

  /* safe K.make with friendly message */
  KIT.model = function (v, prefix) {
    const m = K.make(v, prefix);
    const r0 = m.rate(0);
    if (!isFinite(r0)) throw new Error('The rate is not finite at the feed composition. Check orders and concentrations.');
    return m;
  };

  /* standard rate-law + feed inputs with units wired */
  KIT.kin = (o = {}) => KIT.withUnits(K.inputs(o));

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
