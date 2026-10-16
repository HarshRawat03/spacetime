/* Chapters 28–30 — microbial fermentation: the generalized Monod equation with product poisoning
   (Eq. 28.8), batch, plug and mixed fermentors, washout and the optimum, rate constants from runs;
   and the best arrangement of fermentors for a feed (Examples 29.1–29.2, 30.1). */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, D = ST.D;
  const fmt = (x, s = 4) => N.fmt(x, s);
  const BIO = ST.BIO = {};

  /* p: {k, CM, CRs (≤0 or ∞: none), n, yCA, yRA, CA0, CC0, CR0}; state by the extent Δ = A consumed */
  BIO.model = function (p) {
    const k = +p.k, CM = +p.CM, n = +p.n, yCA = +p.yCA, yRA = +p.yRA, CA0 = +p.CA0, CC0 = +p.CC0 || 0, CR0 = +p.CR0 || 0;
    const CRs = +p.CRs > 0 && isFinite(+p.CRs) ? +p.CRs : Infinity;
    const Dmax = Math.min(CA0, isFinite(CRs) ? Math.max((CRs - CR0) / yRA, 0) : Infinity);
    const CA = d => CA0 - d, CC = d => CC0 + yCA * d, CR = d => CR0 + yRA * d;
    const rC = d => {
      const a = CA0 - d, c = CC0 + yCA * d, rr = CR0 + yRA * d;
      if (!(a > 0) || !(c > 0)) return 0;
      const pois = isFinite(CRs) ? Math.pow(Math.max(1 - rr / CRs, 0), n) : 1;
      return k * pois * c * (CM > 0 ? a / (a + CM) : 1);
    };
    /* mixed flow: Δ = τ r_C(Δ)/y; the largest root, or 0 (washout) */
    const mixed = (tau, d0 = 0) => {
      /* feed already reacted to extent d0 (cells present); exit extent d ≥ d0 */
      const F = d => (d - d0) * yCA - tau * rC(d);
      const hi = Dmax * (1 - 1e-12);
      if (!(hi > d0)) return d0;
      const xs = [];
      const nS = 240;
      for (let i = 0; i <= nS; i++) { const u = i / nS; xs.push(d0 + (hi - d0) * (u < 0.5 ? Math.pow(2 * u, 3) / 2 : 1 - Math.pow(2 * (1 - u), 3) / 2)); }
      let best = null;
      for (let i = 1; i < xs.length; i++) {
        const a = F(xs[i - 1]), b = F(xs[i]);
        if (i === 1 && Math.abs(a) < 1e-300 && rC(xs[0]) === 0) continue;
        if (a * b < 0) { const r = N.brent(F, xs[i - 1], xs[i], { fa: a, fb: b, tol: 1e-14 * Math.max(hi, 1) }); if (isFinite(r)) best = r; }
      }
      if (best === null) return CC0 > 0 || d0 > 0 ? hi : d0;
      return best;
    };
    /* plug or batch: dΔ/dt = r_C/y from Δ0 over time t */
    const plug = (t, d0 = 0) => {
      if (!(t > 0)) return d0;
      if (rC(d0) <= 0) return d0;
      const sol = N.ode((tt, y) => [rC(Math.min(Math.max(y[0], 0), Dmax)) / yCA], 0, [d0], t, { rtol: 1e-9, atol: 1e-12 * Math.max(Dmax, 1), minPts: 30 });
      return Math.min(sol.last()[0], Dmax);
    };
    const tPlug = (d0, d1) => (d1 <= d0 ? 0 : N.quad(d => yCA / rC(d), d0, d1, { rtol: 1e-9 }));
    /* the composition of maximum rate */
    const dMax = () => {
      const f = d => -rC(d);
      const sc = N.fminScan(f, 0, Dmax * (1 - 1e-9), 200);
      return sc.x;
    };
    return { k, CM, n, yCA, yRA, CA0, CC0, CR0, CRs, Dmax, CA, CC, CR, rC, mixed, plug, tPlug, dMax };
  };
  BIO.inputs = function (d = {}, cond) {
    const c = cond || (() => true);
    return [
      { group: 'Kinetics: r_C = k(1 − C_R/C_R*)ⁿ C_AC_C/(C_A + C_M)', if: c },
      { k: 'k', label: 'k', unit: '1/[t]', v: d.k !== undefined ? d.k : 2, min: 1e-12, max: 1e9, scale: 'log', if: c },
      { k: 'CM', label: 'Monod constant C_M (0: substrate never limiting)', v: d.CM !== undefined ? d.CM : 1, min: 0, max: 1e12, if: c },
      { k: 'pois', t: 'check', label: 'Product poisoning', v: !!d.pois, if: c },
      { k: 'CRs', label: 'C_R* (all activity stops)', v: d.CRs || 12, min: 1e-12, max: 1e12, if: v => c(v) && v.pois },
      { k: 'n', label: 'Order of poisoning n', v: d.n || 1, min: 0.05, max: 5, step: 0.05, if: v => c(v) && v.pois },
      { k: 'yCA', label: 'φ(C/A), cells formed per A used', v: d.yCA || 1, min: 1e-9, max: 1e9, if: c },
      { k: 'yRA', label: 'φ(R/A), product formed per A used', v: d.yRA !== undefined ? d.yRA : 1, min: 0, max: 1e9, if: c },
      { group: 'Feed', if: c },
      { k: 'CA0', label: 'C_A0', v: d.CA0 || 3, min: 1e-12, max: 1e12, if: c },
      { k: 'CC0', label: 'C_C0 (cells in the feed)', v: d.CC0 || 0, min: 0, max: 1e12, if: c },
      { k: 'CR0', label: 'C_R0', v: d.CR0 || 0, min: 0, max: 1e12, if: c }
    ];
  };
  /* one batch run sampled on a time grid (a single integration, dense output); growth starts after the lag */
  BIO.batchPath = function (m, T, n = 160, lag = 0) {
    const ts = N.linspace(0, T, n);
    if (!(m.CC0 > 0) || !(T > lag)) return { t: ts, d: ts.map(() => 0) };
    const sol = N.ode((tt, y) => [m.rC(Math.min(Math.max(y[0], 0), m.Dmax)) / m.yCA], 0, [0], T - lag, { rtol: 1e-8, atol: 1e-12 * Math.max(m.Dmax, 1), minPts: 60 });
    return { t: ts, d: ts.map(t => (t <= lag ? 0 : Math.min(Math.max(sol.at(t - lag)[0], 0), m.Dmax))) };
  };
  /* with food in great excess C_A hardly moves: show what does */
  BIO.axis = m => {
    const excess = m.Dmax < 0.2 * m.CA0;
    if (!excess) return { excess, lab: 'C_A', f: m.CA, max: m.CA0 };
    return m.yRA > 0 ? { excess, lab: 'C_R', f: m.CR, max: m.CR(m.Dmax) } : { excess, lab: 'C_A0 − C_A', f: d => d, max: m.Dmax };
  };
  BIO.params = v => ({ k: +v.k, CM: +v.CM, CRs: v.pois ? +v.CRs : Infinity, n: v.pois ? +v.n : 1, yCA: +v.yCA, yRA: +v.yRA, CA0: +v.CA0, CC0: +v.CC0, CR0: +v.CR0 });

  /* layout sketch: a first stage on part of the feed (m mixed tanks in parallel, or a vessel run
     as a recycle loop), the rest bypassed, then stages in series (tanks or plug flow).
     o = {feed, first:{kind:'tanks'|'loop', m, tag, fill, flow, loop}, bypass (label or null), series:[{kind:'tank'|'tube', tag, fill}], out} */
  BIO.sketch = function (o) {
    const W = 50, H = 46, LW = 96, LH = 18, PW = 96, pitch = 104;
    const loop = o.first.kind === 'loop' || o.first.kind === 'tube', m = loop ? 1 : Math.max(1, o.first.m || 1);
    const top = loop ? 58 : 34;
    const y0 = top + (m - 1) * pitch / 2 + (loop ? 0 : H / 2);
    const ys = [];
    for (let i = 0; i < m; i++) ys.push(y0 + (i - (m - 1) / 2) * pitch);
    const xs = 52;
    let s = D.pipe([[4, y0], [xs, y0]], { label: o.feed, dy: -8, arrow: false });
    let xa = xs + 30, xj, yLow;
    if (loop) {
      s += D.pipe([[xs, y0], [xa - 2, y0]], { arrow: true });
      s += D.tube(xa, y0 - LH / 2, LW, LH, { tag: o.first.tag, prof: o.first.prof });
      const xr = xa + LW + 14;
      s += D.pipe([[xa + LW, y0], [xr + 18, y0]], { arrow: false });
      if (o.first.kind === 'loop') s += D.pipe([[xr, y0], [xr, y0 - 32], [xa - 14, y0 - 32], [xa - 14, y0 - 3]], { label: o.first.loop || 'R = ∞', at: 1, dx: -(xr - xa + 14) / 2, dy: -6, anchor: 'middle', flow: false });
      if (o.first.flow) s += D.text(xa + LW / 2, y0 - 14, o.first.flow, { cls: 'pid-r', fs: 10 });
      xj = xr + 18; yLow = y0 + LH / 2 + 30;
    } else {
      const xm = xs + 14, xo = xa + W + 16;
      if (m > 1) { s += D.pipe([[xs, y0], [xm, y0]], { arrow: false }); s += D.pipe([[xm, ys[0]], [xm, ys[m - 1]]], { arrow: false, flow: false }); s += D.pipe([[xo, ys[0]], [xo, ys[m - 1]]], { arrow: false, flow: false }); }
      ys.forEach((y, i) => {
        s += D.pipe([[m > 1 ? xm : xs, y], [xa - 2, y]], { arrow: true });
        s += D.tank(xa, y - H / 2, W, H, { fill: o.first.fill, tag: (o.first.tag || 'F') + '-' + (i + 1) });
        s += D.pipe([[xa + W, y], [xo + (m > 1 ? 0 : 14), y]], { arrow: false });
      });
      if (o.first.flow) s += D.text((xs + xa) / 2 + 4, ys[0] - 8, o.first.flow, { cls: 'pid-r', fs: 10 });
      xj = xo + 14; yLow = ys[m - 1] + H / 2 + 40;
    }
    if (o.bypass) {
      s += D.pipe([[xs, y0], [xs, yLow], [xj, yLow], [xj, y0 + 2]], { label: o.bypass, at: 1, dx: (xj - xs) / 2, dy: -5, anchor: 'middle' });
      s += D.node(xs, y0) + D.node(xj, y0);
    }
    let x = xj;
    (o.series || []).forEach(u => {
      if (u.kind === 'tube') { s += D.pipe([[x, y0], [x + 20, y0]], { arrow: true }); s += D.tube(x + 22, y0 - LH / 2, PW, LH, { tag: u.tag, prof: u.prof }); x += 22 + PW; }
      else { s += D.pipe([[x, y0], [x + 18, y0]], { arrow: true }); s += D.tank(x + 20, y0 - H / 2, W, H, { fill: u.fill, tag: u.tag }); x += 20 + W; }
    });
    s += D.pipe([[x, y0], [x + 78, y0]], { label: o.out, dy: -8 });
    const h = Math.max(o.bypass ? yLow + 16 : 0, (loop ? y0 + LH / 2 : ys[m - 1] + H / 2) + 34);
    return D.svg(x + 86, h, s, o.aria || 'fermentor arrangement');
  };

  /* ===================================================================== */
  /* B-2901  microbial fermentors: Monod and product poisoning              */
  /* ===================================================================== */
  ST.tool({
    id: 'monod', tag: 'B-2901', ch: 29, sec: '29.1–29.2, 30.2', title: 'Microbial fermentation: Monod kinetics, product poisoning, washout and the optimum',
    lede: 'Cells growing on a substrate with the generalized Monod equation of Chapter 28 — food limitation, product poisoning of any order, or both at once — in batch, plug flow and mixed flow. It finds washout, the composition of maximum rate, the best mixed-flow residence time and the most cells or product a fermentor can make; and from mixed-flow or batch runs it extracts k, C_M, C_R* and n with the book’s straight lines, and tells substrate limitation from product limitation.',
    keys: 'microbial fermentation monod equation substrate limiting product limiting poisoning chemostat mixed flow fermentor batch plug washout optimum residence time maximum production rate cells yield fractional yield lag time E coli glucose mannitol lactose yeast alcohol fruit fly example 29.3 30.1 problem 29.1 29.2 29.3 29.6 29.18 30.1 30.2 30.8 30.9 30.10 30.11 30.12',
    units: { t: 'h' },
    inputs: [
      { group: 'Task' },
      { k: 'mode', t: 'seg', label: 'Do', opts: [['reactor', 'fermentor performance'], ['fit', 'rate constants from runs']], v: 'reactor' }
    ].concat(BIO.inputs({ k: 4, CM: 0.4, yCA: 0.1, yRA: 0, CA0: 6 }, v => v.mode === 'reactor')).concat([
      { group: 'Fermentor', if: v => v.mode === 'reactor' },
      { k: 'Vm', label: 'Volume V', unit: 'm³', v: 1, min: 1e-12, max: 1e9, if: v => v.mode === 'reactor' },
      { k: 'vf', label: 'Feed rate v (mixed flow)', unit: 'm³/[t]', v: 3, min: 1e-12, max: 1e12, if: v => v.mode === 'reactor' },
      { k: 'tb', label: 'Batch time (or plug τ)', unit: '[t]', v: 1, min: 0, max: 1e9, if: v => v.mode === 'reactor' },
      { k: 'lag', label: 'Lag time before growth (batch)', unit: '[t]', v: 0, min: 0, max: 1e9, if: v => v.mode === 'reactor' },
      { group: 'Runs', if: v => v.mode === 'fit' },
      { k: 'ftype', t: 'select', label: 'Runs', opts: [['mfs', 'mixed flow: C_A (and C_C) at several τ — substrate limiting'], ['mfp', 'mixed flow: C_R at several τ — product limiting'], ['mfx', 'mixed flow: which limits? (fit both)'], ['bat', 'batch: intervals with Δt, mean C_A, C_C at start and end']], v: 'mfs', if: v => v.mode === 'fit' },
      { k: 'fCA0', label: 'C_A0 of the feed', v: 160, min: 1e-12, max: 1e12, if: v => v.mode === 'fit' && v.ftype !== 'bat' && v.ftype !== 'mfp' },
      { k: 'fV', label: 'Volume (if flows are given)', v: 1, min: 1e-12, max: 1e12, if: v => v.mode === 'fit' && v.ftype !== 'bat' },
      { k: 'fdata', t: 'table', label: 'Runs: flow v (or τ), C_A, C_C', cols: ['v or τ', 'C_A', 'C_C'], v: [[0.2, 4, 15.6], [0.4, 10, 15], [0.8, 40, 12], [1, 100, 6]], if: v => v.mode === 'fit' && (v.ftype === 'mfs' || v.ftype === 'mfx') },
      { k: 'pdata', t: 'table', label: 'Runs: flow v (or τ), C_R', cols: ['v or τ', 'C_R'], v: [[1 / 1.5, 24], [1 / 3, 30]], if: v => v.mode === 'fit' && v.ftype === 'mfp' },
      { k: 'isFlow', t: 'check', label: 'First column is a flow rate (τ = V/v)', v: true, if: v => v.mode === 'fit' && v.ftype !== 'bat' },
      { k: 'fyRA', label: 'φ(R/A) (to turn C_A into C_R)', v: 0.5, min: 0, max: 1e9, if: v => v.mode === 'fit' && v.ftype === 'mfx' },
      { k: 'CRsKnown', label: 'C_R* if known (0 = fit it, n = 1)', v: 36, min: 0, max: 1e12, if: v => v.mode === 'fit' && v.ftype === 'mfp' },
      { k: 'bdata', t: 'table', label: 'Intervals: Δt, mean C_A, C_C at start, C_C at end', cols: ['Δt', 'C̄_A', 'C_C start', 'C_C end'], v: [[0.54, 137, 15.5, 23], [0.36, 114, 23, 30], [0.33, 90, 30, 38.8], [0.35, 43, 38.8, 48.5], [0.37, 29, 48.5, 58.3], [0.38, 9, 58.3, 61.3], [0.37, 2, 61.3, 62.5]], if: v => v.mode === 'fit' && v.ftype === 'bat' }
    ]),
    presets: [
      { name: 'Example 29.3 — E. coli on glucose', v: { mode: 'reactor', k: 4, CM: 0.4, pois: false, yCA: 0.1, yRA: 0, CA0: 6, CC0: 0, CR0: 0, Vm: 1, vf: 3, tb: 1, lag: 0, _tu: 'h' }, book: { vopt: { v: 3, tol: 0.01 }, FAx: { v: 14.4, tol: 0.01 }, FCmax: { v: 1.44, tol: 0.01 } }, note: 'φ(C/A) = 0.1, C_M = 0.4 and C_A0 = 6 kg/m³, V = 1 m³. The rate is printed as 40C_AC_C/(C_A + 0.4), but the solution divides kτ_opt = 1.33 by k = 4, so k = 4 hr⁻¹ is used here; with 40 every flow and rate is ten times larger. Book: N = 4, τ_opt = 1/3 hr, v_opt = 3 m³/hr, 14.4 kg/hr of glucose used and 1.44 kg/hr of cells made.' },
      { name: 'Example 30.1 — fruit fly cocktail', v: { mode: 'reactor', k: Math.sqrt(3), CM: 0, pois: true, CRs: 120, n: 1, yCA: 1, yRA: 1, CA0: 1e9, CC0: 0, CR0: 0, Vm: 30, vf: 25.98, tb: 1, lag: 0, _tu: 'h' }, book: { vopt: { v: 25.98, tol: 0.005 }, FRmax: { v: 1558, tol: 0.005 }, CRopt: { v: 60, tol: 0.005 } }, note: 'Product-limited, k = √3 hr⁻¹, n = 1, C_R* = 120 kg alcohol/m³, a 30 m³ mixed fermentor with food in excess. Book: 6% alcohol, v_opt = 25.98 m³/hr and 1558 kg/hr of alcohol.' },
      { name: 'Problem 30.9 — winter, n = ½', v: { mode: 'reactor', k: Math.sqrt(3), CM: 0, pois: true, CRs: 120, n: 0.5, yCA: 1, yRA: 1, CA0: 1e9, CC0: 0, CR0: 0, Vm: 30, vf: 25.98, tb: 1, lag: 0, _tu: 'h' }, expect: { CRopt: { v: 120 / 1.5, tol: 1e-4 } } },
      { name: 'Problem 30.10 — summer, n = 2', v: { mode: 'reactor', k: Math.sqrt(3), CM: 0, pois: true, CRs: 120, n: 2, yCA: 1, yRA: 1, CA0: 1e9, CC0: 0, CR0: 0, Vm: 30, vf: 25.98, tb: 1, lag: 0, _tu: 'h' }, expect: { CRopt: { v: 40, tol: 1e-4 }, FRmax: { v: Math.sqrt(3) * 30 * 120 * 4 / 27, tol: 1e-3 } } },
      { name: 'Problem 29.2 — mannitol, V = 5 m³', v: { mode: 'reactor', k: 1.2, CM: 2, pois: false, yCA: 0.1, yRA: 0, CA0: 6, CC0: 0, CR0: 0, Vm: 5, vf: 1, tb: 5, lag: 0, _tu: 'h' }, check: { CC: { v: 0.56, tol: 1e-4 } } },
      { name: 'Problem 29.3 — mannitol, V = 1 m³ (washout)', v: { mode: 'reactor', k: 1.2, CM: 2, pois: false, yCA: 0.1, yRA: 0, CA0: 6, CC0: 0, CR0: 0, Vm: 1, vf: 1, tb: 1, lag: 0, _tu: 'h' }, check: { CC: { v: 0, tol: 1e-9 } } },
      { name: 'Problems 30.1–30.2 — poison-limited, one tank', v: { mode: 'reactor', k: 2, CM: 50, pois: true, CRs: 12, n: 1, yCA: 0.01, yRA: 0.1, CA0: 1e6, CC0: 0, CR0: 0, Vm: 1, vf: 1, tb: 1, lag: 0, _tu: 'h' }, check: { CR: { v: 6, tol: 1e-3 } }, note: 'v = 1 gives C_R = 6. At v = 4 (Problem 30.2) the whole tank washes out; B-2902 finds the bypass that saves it (C_R = 1.5).' },
      { name: 'Batch with a lag', v: { mode: 'reactor', k: 1.2, CM: 2, pois: false, yCA: 0.1, yRA: 0.4, CA0: 6, CC0: 0.05, CR0: 0, Vm: 1, vf: 0.5, tb: 6, lag: 1.5, _tu: 'h' } },
      { name: 'Problem 29.1 — E. coli on lactose', v: { mode: 'fit', ftype: 'mfs', fCA0: 160, fV: 1, isFlow: true, fdata: [[0.2, 4, 15.6], [0.4, 10, 15], [0.8, 40, 12], [1, 100, 6]], _tu: 'h' }, check: { k: { v: 1.2, tol: 0.01 }, CM: { v: 20, tol: 0.01 } }, note: 'V = 1 L, C_A0 = 160 mg/L. The points fall exactly on Eq. 29.8.' },
      { name: 'Problem 29.6 — yeast, two flow rates', v: { mode: 'fit', ftype: 'mfs', fCA0: 500, fV: 1, isFlow: true, fdata: [[0.5, 100, 200], [1, 300, 100]], _tu: 'h' }, check: { k: { v: 2, tol: 0.01 }, CM: { v: 300, tol: 0.01 }, yC: { v: 0.5, tol: 0.01 } }, note: 'Both runs make 100 g/hr of cells in 1 m³ from C_A0 = 500 mol/m³, so C_C = 200 and 100 g/m³.' },
      { name: 'Problem 29.18 — Monod’s batch runs (method b)', v: { mode: 'fit', ftype: 'bat', bdata: [[0.54, 137, 15.5, 23], [0.36, 114, 23, 30], [0.33, 90, 30, 38.8], [0.35, 43, 38.8, 48.5], [0.37, 29, 48.5, 58.3], [0.38, 9, 58.3, 61.3], [0.37, 2, 61.3, 62.5]], _tu: 'h' } },
      { name: 'Problem 30.8 — poison limit at 36', v: { mode: 'fit', ftype: 'mfp', fV: 1, isFlow: false, CRsKnown: 36, pdata: [[1.5, 24], [3, 30]], _tu: 'h' }, check: { k: { v: 2, tol: 0.01 } } },
      { name: 'Problem 30.11 — books to glucose', v: { mode: 'fit', ftype: 'mfp', fV: 50, isFlow: true, CRsKnown: 90, pdata: [[16, 54], [4, 75]], _tu: 'h' }, check: { n: { v: 1.58, tol: 0.01 }, k: { v: 1.37, tol: 0.01 } }, note: 'V = 50 L. C_R → 90 as v → 0 gives C_R* directly; two runs then fix k and n (Eq. 30.15).' },
      { name: 'Problem 30.12 — Microbe or Ferment?', v: { mode: 'fit', ftype: 'mfx', fCA0: 150, fV: 46.4, isFlow: true, fyRA: 0.5, fdata: [[4.64, 5, 0], [20, 125, 0], [22, 150, 0]], _tu: 'h' }, book: { k: { v: 0.5, tol: 0.005 }, CM: { v: 20, tol: 0.005 }, kp: { v: 0.5, tol: 0.005 }, CRsP: { v: 90.6, tol: 0.005 } }, check: { verdict: 'substrate limiting' }, note: 'The constants each side claims are printed in the problem (Microbe: k = 0.50, C_R* = 90.6, n = 1; Ferment: C_M = 20, k = 0.50). Both fit the two operating runs; only Monod’s also predicts washout at v = 22.' }
    ],
    eqs: [
      ['r_C = k_{obs}\\frac{C_AC_C}{C_A+C_M}', 'Eq. 28.8', 'The generalized Monod equation.'],
      ['k_{obs} = k\\left(1-\\frac{C_R}{C_R^*}\\right)^n', 'Eq. 28.7'],
      ['C_C - C_{C0} = \\varphi(C/A)(C_{A0}-C_A)', 'Eq. 28.4'],
      ['C_R - C_{R0} = \\varphi(R/A)(C_{A0}-C_A)', 'Eq. 28.4'],
      ['C_{A,max\\ rate} = \\sqrt{C_M^2 + C_M\\left(C_{A0} + \\varphi(A/C)C_{C0}\\right)} - C_M', 'Eq. 29.2'],
      ['k\\tau_m = \\frac{C_M + C_A}{C_A}', 'Eq. 29.7', 'Monod kinetics, no cells in the feed.'],
      ['N = \\sqrt{1+\\frac{C_{A0}}{C_M}},\\qquad k\\tau_{opt} = \\frac{N}{N-1}', 'Eqs. 29.9–29.10'],
      ['k\\tau_{washout} = \\frac{N^2}{N^2-1}', 'Eq. 29.11'],
      ['F_{C,max} = kV\\varphi(C/A)\\,C_{A0}\\,\\frac{N-1}{N+1}', 'Fig. 29.5'],
      ['k\\tau_m = \\frac{(C_{A0}-C_A)(C_A+C_M)}{\\varphi(A/C)C_{C0}C_A + C_A(C_{A0}-C_A)}', 'Eq. 29.12', 'Cells in the feed.'],
      ['\\frac{1}{C_A} = \\frac{k}{C_M}\\tau_m - \\frac{1}{C_M}', 'Eq. 29.8'],
      ['\\frac{C_C}{r_C} = \\frac{1}{k}+\\frac{C_M}{k}\\frac{1}{C_A}', 'Eq. 29.5'],
      ['k\\tau_m = \\frac{1}{(1-C_R/C_R^*)^n}', 'Eqs. 30.10, 30.14', 'Product limiting, C_C0 = C_R0 = 0.'],
      ['C_{R,opt} = \\frac{C_R^*}{n+1},\\qquad F_{R,max} = \\frac{kVC_R^*n^n}{(n+1)^{n+1}}', 'Fig. 30.6'],
      ['C_R = C_R^* - \\frac{C_R^*}{k}\\frac{1}{\\tau_m}', 'Eq. 30.11'],
      ['\\log\\tau_m = -\\log k + n\\log\\frac{C_R^*}{C_R^*-C_R}', 'Eq. 30.15']
    ],
    method: ['Constant fractional yields reduce every composition to one variable, the amount of A consumed. Mixed flow is then the root of (C_A0 − C_A)φ(C/A) = τ r_C (searched over the whole range, so washout and multiple states are found); batch and plug flow integrate dC_A/dt = −r_C/φ(C/A). The production rate of a mixed fermentor per unit volume equals r_C at its exit composition, so the optimum is simply the composition of maximum rate, found by golden-section search — Eqs. 29.10 and 30.4 are its closed forms.', 'Fits use the book’s linearizations by least squares. To decide between substrate and product limitation, both are fitted to the operating runs and then asked to predict every run, including any washout.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'h';
      if (v.mode === 'fit') return fitRun(v, kpi, warn, info, plots, tables, tu);
      const m = BIO.model(BIO.params(v));
      const V = +v.Vm, vf = +v.vf, tau = V / vf;
      const dM = m.mixed(tau);
      kpi.push({ k: 'CA', l: 'Mixed flow: C_A leaving', v: m.CA(dM), main: true, s: dM <= 1e-12 * m.CA0 && m.CC0 === 0 ? 'washout: the cells cannot keep up' : 'τ = ' + fmt(tau, 4) + ' ' + tu + ', kτ = ' + fmt(m.k * tau, 4) });
      kpi.push({ k: 'CC', l: 'Cells leaving C_C', v: m.CC(dM) }, { k: 'CR', l: 'Product leaving C_R', v: m.CR(dM) });
      kpi.push({ k: 'FC', l: 'Cells made, v(C_C − C_C0)', v: vf * (m.CC(dM) - m.CC0), u: 'per ' + tu });
      /* washout */
      let tauW = NaN;
      if (m.CC0 === 0) {
        const r0 = m.rC(1e-9 * m.Dmax) / (m.yCA * 1e-9 * m.Dmax);
        tauW = 1 / r0;
        kpi.push({ k: 'tauW', l: 'Washout below τ', v: tauW, u: tu, s: 'v above ' + fmt(V / tauW, 4) + ' m³/' + tu });
      }
      /* optimum */
      const dO = m.dMax(), rO = m.rC(dO), tauO = (dO) * m.yCA / rO;
      if (m.CC0 === 0) {
        kpi.push({ k: 'tauopt', l: 'Best mixed-flow τ (maximum production)', v: tauO, u: tu, s: 'C_A = ' + fmt(m.CA(dO), 4) + ', C_C = ' + fmt(m.CC(dO), 4) + (isFinite(m.CRs) ? ', C_R = ' + fmt(m.CR(dO), 4) : '') });
        kpi.push({ k: 'vopt', l: 'Best feed rate v_opt = V/τ_opt', v: V / tauO, u: 'm³/' + tu });
        kpi.push({ k: 'FAx', l: 'Food used at the optimum, v(C_A0 − C_A)', v: V / tauO * dO, u: 'per ' + tu });
        kpi.push({ k: 'FCmax', l: 'Most cells per unit time', v: V * rO, u: 'per ' + tu });
        if (m.yRA > 0) { kpi.push({ k: 'FRmax', l: 'Most product per unit time', v: V * rO * m.yRA / m.yCA, u: 'per ' + tu }); kpi.push({ k: 'CRopt', l: 'C_R at the optimum', v: m.CR(dO) }); }
        if (!isFinite(m.CRs) && m.CM > 0) {
          const Nn = Math.sqrt(1 + m.CA0 / m.CM);
          info.push('Monod kinetics: N = √(1 + C_A0/C_M) = ' + fmt(Nn, 4) + ', kτ_opt = N/(N − 1) = ' + fmt(Nn / (Nn - 1), 4) + ' and kτ_washout = N²/(N² − 1) = ' + fmt(Nn * Nn / (Nn * Nn - 1), 4) + ' (Eqs. 29.9–29.11).');
        }
      }
      kpi.push({ k: 'CAmr', l: BIO.axis(m).lab + ' at the maximum rate (Eq. 29.2 / 30.4)', v: BIO.axis(m).f(dO), s: 'r_C,max = ' + fmt(rO, 4) + ' per ' + tu });
      /* batch */
      const tb = +v.tb, lag = +v.lag;
      if (m.CC0 > 0) {
        const dB = m.plug(Math.max(tb - lag, 0));
        kpi.push({ k: 'CAb', l: 'Batch (or plug) after ' + fmt(tb, 4) + ' ' + tu, v: m.CA(dB), s: 'C_C = ' + fmt(m.CC(dB), 4) + (lag > 0 ? '; first ' + fmt(lag, 3) + ' ' + tu + ' lag' : '') });
      } else info.push('Batch and plug flow need cells in the feed (C_C0 > 0): with none, nothing grows. Mixed flow seeds itself.');
      /* charts: mixed flow summary */
      const tmax = Math.max(4 * tauO, 3 * (isFinite(tauW) ? tauW : tauO), 1.5 * tau);
      const tg = N.linspace(tmax / 400, tmax, 220);
      const dg = tg.map(t => m.mixed(t));
      const xs1 = BIO.axis(m);
      const ser = [{ name: xs1.excess ? 'C_A0 − C_A' : 'C_A', x: tg.map(t => m.k * t), y: dg.map(xs1.excess ? (d => d) : m.CA), color: 1, w: 2.2 }, { name: 'C_C', x: tg.map(t => m.k * t), y: dg.map(m.CC), color: 3, w: 2.2 }];
      if (m.yRA > 0) ser.push({ name: 'C_R', x: tg.map(t => m.k * t), y: dg.map(m.CR), color: 4, w: 2.2 });
      ser.push({ name: 'production per volume, r_C (scaled)', x: tg.map(t => m.k * t), y: dg.map(d => m.rC(d) / Math.max(rO, 1e-300) * Math.max(m.CC(dO), m.CR(dO), 1e-300)), color: 'var(--muted)', dash: '4 4' });
      const marks = [{ t: 'v', x: m.k * tau, color: 8, label: 'this fermentor', dy: 30 }];
      if (m.CC0 === 0) { marks.push({ t: 'v', x: m.k * tauO, color: 2, label: 'optimum', dy: 15 }); if (isFinite(tauW)) marks.push({ t: 'v', x: m.k * tauW, color: 'var(--muted)', label: 'washout' }); }
      plots.push({ title: 'Mixed fermentor against kτ (Figs. 29.5, 30.4, 30.6)', x: { label: 'kτ_m', min: 0, max: m.k * tmax }, y: { label: 'concentration', min: 0 }, series: ser, marks, hmax: 460, cap: 'Below washout no cells survive; the maximum production rate comes at the composition of maximum rate.' });
      /* 1/r against C_A */
      const dd = N.linspace(m.Dmax * 1e-4, m.Dmax * 0.999, 240);
      plots.push({ title: '1/(−r_A) against ' + xs1.lab + ': the U-shaped curve (Fig. 29.6)', x: { label: xs1.lab, min: 0 }, y: { label: '1/(−r_A)', log: true }, series: [{ name: '1/(−r_A)', x: dd.map(xs1.f), y: dd.map(d => m.yCA / m.rC(d)), color: 1, w: 2.2 }], marks: [{ t: 'v', x: xs1.f(dO), color: 2, label: 'max rate' }], cap: 'To reach any C_A between the minimum and C_A0, mix part of the feed at the minimum; to go below it, mixed flow to the minimum and plug flow beyond (Ch. 29.3).' });
      /* batch profile */
      if (m.CC0 > 0) {
        const tEnd = Math.max(tb, lag + 2 * m.tPlug(0, m.Dmax * 0.99)) || tb;
        const bp = BIO.batchPath(m, tEnd, 200, lag), ts = bp.t, ds = bp.d;
        const bs = [{ name: xs1.excess ? 'C_A0 − C_A' : 'C_A', x: ts, y: ds.map(xs1.excess ? (d => d) : m.CA), color: 1, w: 2.2 }, { name: 'C_C', x: ts, y: ds.map(m.CC), color: 3, w: 2.2 }];
        if (m.yRA > 0) bs.push({ name: 'C_R', x: ts, y: ds.map(m.CR), color: 4, w: 2.2 });
        plots.push({ title: 'Batch or plug flow fermentor (Figs. 29.1, 30.1)', x: { label: 't', unit: tu, min: 0 }, y: { label: 'concentration', min: 0 }, series: bs, marks: lag > 0 ? [{ t: 'band', x0: 0, x1: lag, color: 'var(--grid)' }, { t: 'text', x: lag / 2, y: (xs1.excess ? m.Dmax : m.CA0) * 0.5, s: 'lag', anchor: 'middle' }] : [] });
      }
      if (isFinite(m.CRs)) {
        const cr = N.linspace(0, 1, 100);
        plots.push({ title: 'Poisoning by product: k_obs/k against C_R/C_R* (Fig. 28.4)', x: { label: 'C_R/C_R*', min: 0, max: 1 }, y: { label: 'k_obs/k', min: 0, max: 1 }, series: [0.5, 1, 2].map((nn, i) => ({ name: 'n = ' + nn, x: cr, y: cr.map(q => Math.pow(1 - q, nn)), color: i + 2 })).concat([{ name: 'this n = ' + fmt(m.n, 3), x: cr, y: cr.map(q => Math.pow(1 - q, m.n)), color: 1, w: 2.4 }]) });
      }
      return { kpi, warn, info, plots, tables };
    }
  });

  function fitRun(v, kpi, warn, info, plots, tables, tu) {
    if (v.ftype === 'bat') {
      const rows = (v.bdata || []).filter(r => r && r.slice(0, 4).every(z => z !== '' && isFinite(+z))).map(r => {
        const dt = +r[0], CA = +r[1], c1 = +r[2], c2 = +r[3];
        const Cm = Math.abs(c2 - c1) > 1e-12 ? (c2 - c1) / Math.log(c2 / c1) : c1;
        return { dt, CA, Cm, r: (c2 - c1) / dt };
      }).filter(q => q.r > 0 && q.CA > 0);
      if (rows.length < 2) throw new Error('Enter at least two intervals with growth.');
      const L = N.linreg(rows.map(q => 1 / q.CA), rows.map(q => q.Cm / q.r));
      const k = 1 / L.a, CM = L.b * k;
      kpi.push({ k: 'k', l: 'k (method b, Eq. 29.5)', v: k, u: '1/' + tu, main: true, s: 'R² = ' + fmt(L.r2, 3) });
      kpi.push({ k: 'CM', l: 'C_M', v: CM });
      if (!(L.a > 0)) warn.push('The intercept is not positive: these intervals do not determine k.');
      tables.push({ title: 'Intervals', cols: ['Δt', 'C̄_A', 'C̄_C (log mean)', 'r_C = ΔC_C/Δt', 'C_C/r_C', '1/C_A'], rows: rows.map(q => [q.dt, q.CA, q.Cm, q.r, q.Cm / q.r, 1 / q.CA]) });
      const xs = [0, Math.max(...rows.map(q => 1 / q.CA)) * 1.05];
      plots.push({ title: 'Method (b): C_C/r_C against 1/C_A (Fig. 29.3)', x: { label: '1/C_A', min: 0 }, y: { label: 'C_C/r_C', unit: tu, min: 0 }, series: [{ name: 'intervals', x: rows.map(q => 1 / q.CA), y: rows.map(q => q.Cm / q.r), type: 'pts', color: 1 }, { name: 'slope C_M/k, intercept 1/k', x: xs, y: xs.map(x => L.a + L.b * x), color: 2 }] });
      info.push('Rates from differences are noisy, especially near the end of a run where little changes (Ch. 29.1 comments); mixed-flow runs give these constants far more directly.');
      return { kpi, warn, info, plots, tables };
    }
    const V = +v.fV;
    if (v.ftype === 'mfp') {
      const rows = (v.pdata || []).filter(r => r && isFinite(+r[0]) && isFinite(+r[1]) && r[0] !== '' && r[1] !== '').map(r => ({ tau: v.isFlow ? V / +r[0] : +r[0], CR: +r[1] }));
      if (rows.length < 2) throw new Error('Enter at least two runs.');
      const CRk = +v.CRsKnown;
      if (CRk > 0) {
        const pts = rows.filter(q => q.CR < CRk).map(q => ({ x: Math.log10(CRk / (CRk - q.CR)), y: Math.log10(q.tau) }));
        const L = N.linreg(pts.map(p => p.x), pts.map(p => p.y));
        const n = L.b, k = Math.pow(10, -L.a);
        kpi.push({ k: 'k', l: 'k', v: k, u: '1/' + tu, main: true, s: 'from Eq. 30.15 with C_R* = ' + fmt(CRk, 4) });
        kpi.push({ k: 'n', l: 'Order of poisoning n', v: n, s: pts.length === 2 ? 'two runs: an exact fit' : 'R² = ' + fmt(L.r2, 4) });
        kpi.push({ k: 'CRopt', l: 'Best C_R in mixed flow C_R*/(n + 1)', v: CRk / (n + 1) });
        kpi.push({ k: 'vopt', l: 'Best feed rate kV(n/(n + 1))ⁿ', v: k * V * Math.pow(n / (n + 1), n), u: 'per ' + tu });
        kpi.push({ k: 'FRmax', l: 'Most product per unit time', v: k * V * CRk * Math.pow(n, n) / Math.pow(n + 1, n + 1), u: 'per ' + tu });
        const xs = [0, Math.max(...pts.map(p => p.x)) * 1.15];
        plots.push({ title: 'Eq. 30.15: log τ against log C_R*/(C_R* − C_R) (Fig. 30.7)', x: { label: 'log₁₀ C_R*/(C_R* − C_R)', min: 0 }, y: { label: 'log₁₀ τ' }, series: [{ name: 'runs', x: pts.map(p => p.x), y: pts.map(p => p.y), type: 'pts', color: 1 }, { name: 'slope n', x: xs, y: xs.map(x => L.a + L.b * x), color: 2 }] });
      } else {
        const L = N.linreg(rows.map(q => 1 / q.tau), rows.map(q => q.CR));
        const CRs = L.a, k = -CRs / L.b;
        kpi.push({ k: 'k', l: 'k (n = 1, Eq. 30.11)', v: k, u: '1/' + tu, main: true }, { k: 'CRs', l: 'C_R*', v: CRs });
        const xs = [0, Math.max(...rows.map(q => 1 / q.tau)) * 1.2];
        plots.push({ title: 'Eq. 30.11: C_R against 1/τ (Fig. 30.3)', x: { label: '1/τ', unit: '1/' + tu, min: 0 }, y: { label: 'C_R', min: 0 }, series: [{ name: 'runs', x: rows.map(q => 1 / q.tau), y: rows.map(q => q.CR), type: 'pts', color: 1 }, { name: 'intercept C_R*, washout at 1/τ = k', x: xs, y: xs.map(x => L.a + L.b * x), color: 2 }] });
      }
      return { kpi, warn, info, plots, tables };
    }
    /* substrate limiting (and the discrimination) */
    const CA0 = +v.fCA0;
    const rows = (v.fdata || []).filter(r => r && isFinite(+r[0]) && isFinite(+r[1]) && r[0] !== '' && r[1] !== '').map(r => ({ tau: v.isFlow ? V / +r[0] : +r[0], CA: +r[1], CC: r[2] !== '' && isFinite(+r[2]) ? +r[2] : NaN }));
    const op = rows.filter(q => q.CA < CA0 * (1 - 1e-9)), wo = rows.filter(q => q.CA >= CA0 * (1 - 1e-9));
    if (op.length < 2) throw new Error('Enter at least two runs that are not washed out.');
    const L = N.linreg(op.map(q => q.tau), op.map(q => 1 / q.CA));
    const CM = -1 / L.a, k = L.b * CM;
    kpi.push({ k: 'k', l: 'k (Eq. 29.8)', v: k, u: '1/' + tu, main: v.ftype === 'mfs' });
    kpi.push({ k: 'CM', l: 'C_M', v: CM });
    const ys = op.filter(q => isFinite(q.CC) && q.CC > 0).map(q => q.CC / (CA0 - q.CA));
    if (ys.length) kpi.push({ k: 'yC', l: 'φ(C/A) from the cell concentrations', v: N.mean(ys), s: ys.length > 1 ? 'spread ' + fmt((Math.max(...ys) - Math.min(...ys)) / N.mean(ys) * 100, 3) + '%' : '' });
    if (CM > 0 && k > 0) {
      const Nn = Math.sqrt(1 + CA0 / CM), tauO = Nn / ((Nn - 1) * k);
      kpi.push({ k: 'tauW', l: 'Washout below τ = N²/((N² − 1)k)', v: Nn * Nn / ((Nn * Nn - 1) * k), u: tu });
      kpi.push({ k: 'tauopt', l: 'Best τ = N/((N − 1)k)', v: tauO, u: tu, s: 'v_opt = ' + fmt(V / tauO, 4) });
      if (ys.length) kpi.push({ k: 'FCmax', l: 'Most cells per unit time', v: V / tauO * N.mean(ys) * CA0 * Nn / (Nn + 1), u: 'per ' + tu });
    }
    const xs = [0, Math.max(...rows.map(q => q.tau)) * 1.1];
    plots.push({ title: 'Eq. 29.8: 1/C_A against τ (Fig. 29.4)', x: { label: 'τ_m', unit: tu, min: 0 }, y: { label: '1/C_A' }, series: [{ name: 'runs', x: op.map(q => q.tau), y: op.map(q => 1 / q.CA), type: 'pts', color: 1 }, { name: 'slope k/C_M, intercept −1/C_M', x: xs, y: xs.map(x => L.a + L.b * x), color: 2 }].concat(wo.length ? [{ name: 'washed out', x: wo.map(q => q.tau), y: wo.map(q => 1 / q.CA), type: 'pts', color: 5, marker: 's' }] : []), marks: [{ t: 'h', y: 1 / CA0, color: 'var(--muted)', label: '1/C_A0' }] });
    if (v.ftype === 'mfx') {
      const yRA = +v.fyRA;
      const prA = op.map(q => ({ tau: q.tau, CR: yRA * (CA0 - q.CA) }));
      const Lp = N.linreg(prA.map(q => 1 / q.tau), prA.map(q => q.CR));
      const CRs = Lp.a, kp = -CRs / Lp.b;
      kpi.push({ k: 'kp', l: 'Product-limited fit (n = 1): k', v: kp, u: '1/' + tu }, { k: 'CRsP', l: 'Product-limited fit: C_R*', v: CRs });
      /* predictions for every run */
      const predS = q => { const a = k * q.tau - 1; return a > 0 && CM / a < CA0 ? CM / a : CA0; };
      const predP = q => { const kt = kp * q.tau; if (kt <= 1) return CA0; return CA0 - CRs * (1 - 1 / kt) / yRA; };
      const errS = rows.map(q => Math.abs(predS(q) - q.CA) / CA0), errP = rows.map(q => Math.abs(predP(q) - q.CA) / CA0);
      const eS = Math.max(...errS), eP = Math.max(...errP);
      const verdict = eS < 0.02 && eP >= 0.02 ? 'substrate limiting' : eP < 0.02 && eS >= 0.02 ? 'product limiting' : eS < 0.02 ? 'cannot tell' : 'neither fits';
      kpi.unshift({ k: 'verdict', l: 'Which limits growth', v: verdict, main: true, s: 'largest misfit: Monod ' + fmt(eS * 100, 3) + '%, poisoning ' + fmt(eP * 100, 3) + '% of C_A0' });
      tables.push({ title: 'Each run against both rate equations', cols: ['τ', 'C_A measured', 'Monod predicts', 'poisoning predicts'], rows: rows.map(q => [q.tau, q.CA, predS(q), predP(q)]) });
      info.push('The two forms give look-alike curves in mixed flow, so two operating runs fit either one (Ch. 30, Discussion). They differ at washout: for Monod kinetics it depends on C_A0, for product poisoning it comes at kτ = 1 whatever the feed.');
    }
    return { kpi, warn, info, plots, tables };
  }

  /* ===================================================================== */
  /* B-2902  the best arrangement of fermentors                             */
  /* ===================================================================== */
  /* m mixed tanks in parallel on a fraction β of the feed, bypass mixed in, then the rest in series */
  function mfrNet(m, Ntot, V, vf, beta, M) {
    const t1 = m * V / (beta * vf);
    let d = M.mixed(t1, 0);
    d = beta * d;
    const tauS = V / vf;
    for (let i = m; i < Ntot; i++) d = M.mixed(tauS, d);
    return d;
  }
  function bestMfr(Ntot, V, vf, M) {
    let best = null;
    for (let m = 1; m <= Ntot; m++) {
      const f = b => -mfrNet(m, Ntot, V, vf, b, M);
      const bs = N.linspace(0.005, 1, 80);
      let bi = 0, bf = Infinity;
      bs.forEach((b, i) => { const q = f(b); if (q < bf) { bf = q; bi = i; } });
      let b = bs[bi];
      if (bi > 0 && bi < bs.length - 1) { const r = N.fmin(f, bs[bi - 1], bs[bi + 1]); if (r.f < bf) { b = r.x; bf = r.f; } }
      if (!best || -bf > best.d + 1e-12 * M.Dmax) best = { m, beta: b, d: -bf };
    }
    return best;
  }
  /* PFR with a recycle ratio R on a fraction β of the feed (no side tap) */
  function pfrRec(V, vf, R, beta, M) {
    const tauR = V / ((R + 1) * beta * vf);
    if (R <= 1e-9) return beta * M.plug(tauR, 0);
    const G = d => M.plug(tauR, R * d / (R + 1)) - d;
    const hi = M.Dmax * (1 - 1e-10);
    const xs = N.linspace(1e-9 * hi, hi, 120);
    let root = 0;
    for (let i = 1; i < xs.length; i++) { const a = G(xs[i - 1]), b = G(xs[i]); if (a * b < 0) root = N.brent(G, xs[i - 1], xs[i], { fa: a, fb: b, tol: 1e-12 * hi }); }
    return beta * root;
  }

  ST.tool({
    id: 'fermnet', tag: 'B-2902', ch: 29, sec: '29.3', title: 'The best arrangement of fermentors: bypass, recycle, staging',
    lede: 'For autocatalytic growth the rate curve is U-shaped, so the best layout is never obvious. Given a feed and the vessels you have — one to four mixed tanks, or a plug-flow fermentor with or without a side tap — this tool finds the arrangement of bypass, parallel and series flow, or recycle, that drives C_A lowest (and C_C or C_R highest); or, for a target, the smallest volume by Chapter 29’s rule: mixed flow to the composition of maximum rate, plug flow beyond, part-bypass if the target lies above it.',
    keys: 'optimum arrangement fermentors bypass recycle side tap mixed flow reactors in parallel series plug flow fermentor with recycle smallest reactor size lowest concentration autocatalytic U shaped curve example 29.1 29.2 problem 29.7 29.8 29.9 29.10 29.11 29.12 29.13 29.14 29.15 29.16 29.17 29.19 30.1 30.2 30.3 30.4 30.5 30.6 30.7',
    units: { t: 'h' }, slow: true,
    inputs: BIO.inputs({ k: 2, CM: 1, yCA: 1, yRA: 0, CA0: 3 }).concat([
      { group: 'Equipment and feed' },
      { k: 'eq', t: 'select', label: 'Problem', opts: [['mfr', 'lowest C_A with N equal mixed tanks'], ['pfr', 'lowest C_A with one plug-flow vessel'], ['tgt', 'smallest volume for a target']], v: 'mfr' },
      { k: 'vf', label: 'Feed rate v', v: 3, min: 1e-12, max: 1e12 },
      { k: 'Nt', label: 'Number of mixed tanks', v: 1, min: 1, max: 4, step: 1, if: v => v.eq === 'mfr' },
      { k: 'Vt', label: v => (v.eq === 'pfr' ? 'Plug-flow volume V_p' : 'Volume of each tank V_m'), v: 1, min: 1e-12, max: 1e12, if: v => v.eq !== 'tgt' },
      { k: 'tap', t: 'check', label: 'Side tap allowed (part of the vessel can run as a recycle loop)', v: true, if: v => v.eq === 'pfr' || v.eq === 'tgt' },
      { k: 'tq', t: 'select', label: 'Target', opts: [['CA', 'C_A leaving'], ['CC', 'C_C leaving'], ['CR', 'C_R leaving']], v: 'CC', if: v => v.eq === 'tgt' },
      { k: 'tval', label: 'Target value', v: 1.25, min: 0, max: 1e12, if: v => v.eq === 'tgt' }
    ]),
    presets: [
      { name: 'Example 29.1a — one tank, v = 3', v: { k: 2, CM: 1, pois: false, yCA: 1, yRA: 0, CA0: 3, CC0: 0, CR0: 0, eq: 'mfr', vf: 3, Nt: 1, Vt: 1 }, book: { through: { v: 1, tol: 0.03 } }, note: 'k = 2, C_A0 = 3, C_M = 1, φ(C/A) = 1, V_m = 1. Book: τ = 1/3 would wash out; send v = 1 through the tank (τ_opt = 1) and bypass the rest.' },
      { name: 'Example 29.1b — one tank, v = 1', v: { k: 2, CM: 1, pois: false, yCA: 1, yRA: 0, CA0: 3, CC0: 0, CR0: 0, eq: 'mfr', vf: 1, Nt: 1, Vt: 1 }, book: { through: { v: 1, tol: 0.01 } }, note: 'τ = 1 is the optimum: pass everything.' },
      { name: 'Example 29.1c — one tank, v = 1/3', v: { k: 2, CM: 1, pois: false, yCA: 1, yRA: 0, CA0: 3, CC0: 0, CR0: 0, eq: 'mfr', vf: 1 / 3, Nt: 1, Vt: 1 }, book: { through: { v: 1 / 3, tol: 0.01 } } },
      { name: 'Example 29.1d — two tanks, v = 3', v: { k: 2, CM: 1, pois: false, yCA: 1, yRA: 0, CA0: 3, CC0: 0, CR0: 0, eq: 'mfr', vf: 3, Nt: 2, Vt: 1 }, book: { through: { v: 2, tol: 0.03 } }, note: 'Book: both tanks in parallel at τ = 1, each taking v = 1, with 1 bypassed.' },
      { name: 'Example 29.1e — two tanks, v = 1.5', v: { k: 2, CM: 1, pois: false, yCA: 1, yRA: 0, CA0: 3, CC0: 0, CR0: 0, eq: 'mfr', vf: 1.5, Nt: 2, Vt: 1 }, book: { through: { v: 1, tol: 0.03 } }, note: 'Book: the first tank at its optimum (v = 1) with 0.5 bypassed, then everything through the second. The table answers Problem 29.19, comparing it with the all-series and all-parallel schemes.' },
      { name: 'Example 29.1f — two tanks, v = 0.5', v: { k: 2, CM: 1, pois: false, yCA: 1, yRA: 0, CA0: 3, CC0: 0, CR0: 0, eq: 'mfr', vf: 0.5, Nt: 2, Vt: 1 }, book: { through: { v: 0.5, tol: 0.01 } }, note: 'Book: both in series.' },
      { name: 'Example 29.2a — plug flow, V_p = 3, v = 2', v: { k: 2, CM: 1, pois: false, yCA: 1, yRA: 0, CA0: 3, CC0: 0, CR0: 0, eq: 'pfr', vf: 2, Vt: 3, tap: true }, book: { Vmix: { v: 2, tol: 0.03 } }, note: 'Book: run the first 2 volume units as a recycle loop at the optimum (approaching mixed flow), then 1 unit in plug flow.' },
      { name: 'Example 29.2b — plug flow, V_p = 2, v = 3', v: { k: 2, CM: 1, pois: false, yCA: 1, yRA: 0, CA0: 3, CC0: 0, CR0: 0, eq: 'pfr', vf: 3, Vt: 2, tap: true }, book: { through: { v: 2, tol: 0.03 } }, note: 'Book: all of it as a recycle loop at τ = 1 with v = 2 through it, the rest bypassed.' },
      { name: 'Problem 29.4 — mannitol, V = 5: can bypass help?', v: { k: 1.2, CM: 2, pois: false, yCA: 0.1, yRA: 0, CA0: 6, CC0: 0, CR0: 0, eq: 'mfr', vf: 1, Nt: 1, Vt: 5 }, check: { CC: { v: 0.56, tol: 1e-3 } }, note: 'τ = 5 is already beyond τ_opt, so nothing beats sending it all through.' },
      { name: 'Problem 29.5 — mannitol, V = 1: bypass saves it', v: { k: 1.2, CM: 2, pois: false, yCA: 0.1, yRA: 0, CA0: 6, CC0: 0, CR0: 0, eq: 'mfr', vf: 1, Nt: 1, Vt: 1 }, check: { CC: { v: 0.24, tol: 1e-3 } }, note: 'The whole feed washes out; 0.6 through the tank at τ_opt leaves C_A = 2 there, and the mixed stream carries C_C = 0.24.' },
      { name: 'Problems 29.7–29.8 — lowest C_A, one tank', v: { k: 2, CM: 1, pois: false, yCA: 0.5, yRA: 1, CA0: 3, CC0: 0, CR0: 0, eq: 'mfr', vf: 1 / 3, Nt: 1, Vt: 1 }, note: 'φ(C/A) = 0.5; v = 1/3 for 29.7 and v = 3 for 29.8.' },
      { name: 'Problems 29.9–29.10 — two tanks', v: { k: 2, CM: 1, pois: false, yCA: 0.5, yRA: 1, CA0: 3, CC0: 0, CR0: 0, eq: 'mfr', vf: 2, Nt: 2, Vt: 1 } },
      { name: 'Problems 29.11–29.12 — three tanks', v: { k: 2, CM: 1, pois: false, yCA: 0.5, yRA: 1, CA0: 3, CC0: 0, CR0: 0, eq: 'mfr', vf: 6, Nt: 3, Vt: 1 } },
      { name: 'Problem 29.13 — smallest plug vessel for C_C = 0.5', v: { k: 2, CM: 1, pois: false, yCA: 0.5, yRA: 1, CA0: 3, CC0: 0, CR0: 0, eq: 'tgt', vf: 3, tap: true, tq: 'CC', tval: 0.5 }, check: { V: { v: 1.5, tol: 1e-3 } } },
      { name: 'Problem 29.14 — C_C = 1.25, no side tap', v: { k: 2, CM: 1, pois: false, yCA: 0.5, yRA: 1, CA0: 3, CC0: 0, CR0: 0, eq: 'tgt', vf: 3, tap: false, tq: 'CC', tval: 1.25 }, check: { Vrec: { v: 3.9034, tol: 1e-3 } }, note: 'Recycle over the whole vessel, R ≈ 2.1; with a side tap it would be 3.79.' },
      { name: 'Problem 29.15 — C_C = 1.44, side tap', v: { k: 2, CM: 1, pois: false, yCA: 0.5, yRA: 1, CA0: 3, CC0: 0, CR0: 0, eq: 'tgt', vf: 3, tap: true, tq: 'CC', tval: 1.44 }, check: { V: { v: 4.7894, tol: 1e-3 } } },
      { name: 'Problem 29.16 — lowest C_A, V_p = 4, v = 6', v: { k: 2, CM: 1, pois: false, yCA: 0.5, yRA: 1, CA0: 3, CC0: 0, CR0: 0, eq: 'pfr', vf: 6, Vt: 4, tap: true } },
      { name: 'Problems 30.1–30.4 — poison-limited, tanks (v = 4)', v: { k: 2, CM: 50, pois: true, CRs: 12, n: 1, yCA: 0.01, yRA: 0.1, CA0: 1e6, CC0: 0, CR0: 0, eq: 'mfr', vf: 4, Nt: 1, Vt: 1 }, check: { CR: { v: 1.5, tol: 1e-3 } }, note: 'Problem 30.2 (v = 4): one unit through the tank at its optimum, three bypassed. Set v = 1 for 30.1, two tanks with v = 1 or 3 for 30.3–30.4.' },
      { name: 'Problems 30.5–30.7 — plug vessel for C_R = 6, 4, 9', v: { k: 2, CM: 50, pois: true, CRs: 12, n: 1, yCA: 0.01, yRA: 0.1, CA0: 1e6, CC0: 0, CR0: 0, eq: 'tgt', vf: 3, tap: true, tq: 'CR', tval: 6 }, check: { V: { v: 3, tol: 1e-3 } }, note: 'C_R = 6 is the maximum-rate point: a recycle loop of V = 3 at τ = 1. Targets 4 and 9 need 2 and 4.65.' }
    ],
    eqs: [
      ['\\tau_m = \\frac{C_{A0}-C_A}{-r_A},\\qquad \\tau_p = \\int_{C_A}^{C_{A0}}\\frac{dC_A}{-r_A}', 'Ch. 5'],
      ['\\frac{V}{v} = (R+1)\\int_{C_{A1}}^{C_{Af}}\\frac{-dC_A}{-r_A},\\qquad C_{A1} = \\frac{C_{A0}+RC_{Af}}{R+1}', 'Ch. 6', 'Plug flow with recycle.'],
      ['C_{A,opt} \\to \\text{mixed flow},\\quad C_A < C_{A,opt} \\to \\text{then plug flow}', 'Ch. 29.3', 'The rule for U-shaped 1/r curves.']
    ],
    method: ['Mixed tanks: every arrangement of the form “m tanks in parallel on a fraction β of the feed, the remainder bypassed, then the other tanks in series” is evaluated, β optimized by a scan and golden-section search for each m. This family contains all of the book’s answers for one and two tanks. A plug-flow vessel with a side tap is split into a recycle loop (treated as mixed flow) on part of the feed, then plug flow for everything; without a side tap the whole vessel runs with one recycle ratio R, whose steady state is the fixed point of the plug-flow integral, and R and the bypass are optimized together.', 'Smallest volume for a target: with a side tap the book’s rule is exact — mixed flow to the composition of maximum rate (with the right fraction bypassed if the target lies above it) and plug flow beyond. Without one, V(R) = (R + 1)v∫dC_A/(−r_A) is minimized over R, with bypass when that helps.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'h';
      const M = BIO.model(BIO.params(v));
      const vf = +v.vf;
      const out = d => ({ CA: M.CA(d), CC: M.CC(d), CR: M.CR(d) });
      const dO = M.dMax(), rO = M.rC(dO), tauO = dO * M.yCA / rO;
      kpi.push({ k: 'CAopt', l: BIO.axis(M).lab + ' at the maximum rate', v: BIO.axis(M).f(dO), s: 'τ_opt in one mixed tank = ' + fmt(tauO, 4) + ' ' + tu });
      if (v.eq === 'mfr') {
        const Nt = Math.max(1, Math.min(4, Math.round(+v.Nt))), V = +v.Vt;
        const b = bestMfr(Nt, V, vf, M), o = out(b.d);
        leave(kpi, o, M);
        kpi.push({ k: 'through', l: 'Feed sent through the first stage', v: b.beta * vf, s: (b.m > 1 ? b.m + ' tanks in parallel' : 'one tank') + (b.beta < 0.999 ? ', ' + fmt((1 - b.beta) * vf, 4) + ' bypassed' : ', nothing bypassed') + (Nt > b.m ? '; then ' + (Nt - b.m) + ' in series' : '') });
        const allSeries = mfrNet(1, Nt, V, vf, 1, M), allPar = mfrNet(Nt, Nt, V, vf, 1, M);
        tables.push({ title: 'Compared with the obvious layouts', cols: ['layout', 'C_A leaving', 'C_C leaving'], rows: [['best found', o.CA, o.CC], ['all in series, no bypass', M.CA(allSeries), M.CC(allSeries)], ['all in parallel, no bypass', M.CA(allPar), M.CC(allPar)]] });
        const ax = BIO.axis(M), d1 = M.mixed(b.m * V / (b.beta * vf), 0);
        const ser = [];
        let dd = b.beta * d1;
        for (let i = b.m; i < Nt; i++) { dd = M.mixed(V / vf, dd); ser.push({ kind: 'tank', tag: 'F-' + (i + 1), fill: dd / M.Dmax }); }
        const diagram = BIO.sketch({ feed: 'v = ' + fmt(vf, 3), first: { kind: 'tanks', m: b.m, tag: 'F', fill: d1 / M.Dmax, flow: (b.m > 1 ? fmt(b.beta * vf / b.m, 3) + ' each' : fmt(b.beta * vf, 3)) }, bypass: b.beta < 0.999 ? 'bypass ' + fmt((1 - b.beta) * vf, 3) : null, series: ser, out: ax.lab + ' = ' + fmt(ax.f(b.d), 3) });
        return { kpi, warn, info, plots: [curve(M, dO, undefined, b.d)], tables, diagram };
      }
      if (v.eq === 'pfr') {
        const V = +v.Vt;
        let best;
        if (v.tap) {
          const f = z => {
            const fr = Math.min(Math.max(z[0], 0), 1), beta = Math.min(Math.max(z[1], 1e-4), 1);
            const V1 = fr * V;
            let d = V1 > 0 ? M.mixed(V1 / (beta * vf), 0) : 0;
            d = beta * d;
            d = M.plug((V - V1) / vf, d);
            return -d;
          };
          let bf = Infinity, bz = null;
          N.linspace(0.02, 1, 25).forEach(fr => N.linspace(0.04, 1, 25).forEach(be => { const q = f([fr, be]); if (q < bf) { bf = q; bz = [fr, be]; } }));
          const r = N.nelderMead(f, bz, { tol: 1e-12, maxIt: 400, step: [0.04, 0.04] });
          const z = r.f < bf ? r.x : bz;
          best = { Vmix: Math.min(Math.max(z[0], 0), 1) * V, beta: Math.min(Math.max(z[1], 1e-4), 1), d: -Math.min(r.f, bf) };
          leave(kpi, out(best.d), M);
          kpi.push({ k: 'Vmix', l: 'Volume run as a recycle loop (≈ mixed flow)', v: best.Vmix, s: fmt(V - best.Vmix, 4) + ' then in plug flow' });
          kpi.push({ k: 'through', l: 'Feed sent through the loop', v: best.beta * vf, s: best.beta < 0.999 ? fmt((1 - best.beta) * vf, 4) + ' bypassed' : 'nothing bypassed' });
        } else {
          const f = z => { const R = Math.pow(10, z[0]), beta = Math.min(Math.max(z[1], 1e-4), 1); return -pfrRec(V, vf, R, beta, M); };
          let bf = Infinity, bz = null;
          N.linspace(-2, 3, 16).forEach(lr => N.linspace(0.05, 1, 16).forEach(be => { const q = f([lr, be]); if (q < bf) { bf = q; bz = [lr, be]; } }));
          const r = N.nelderMead(f, bz, { tol: 1e-10, maxIt: 200, step: [0.3, 0.06] });
          const z = r.f < bf ? r.x : bz;
          best = { R: Math.pow(10, z[0]), beta: Math.min(Math.max(z[1], 1e-4), 1), d: -Math.min(r.f, bf) };
          leave(kpi, out(best.d), M);
          kpi.push({ k: 'R', l: 'Best recycle ratio', v: best.R }, { k: 'through', l: 'Feed sent through', v: best.beta * vf, s: best.beta < 0.999 ? fmt((1 - best.beta) * vf, 4) + ' bypassed' : '' });
        }
        const ax = BIO.axis(M);
        const sk = best.Vmix !== undefined
          ? { first: { kind: 'loop', tag: 'loop ' + fmt(best.Vmix, 3), flow: fmt(best.beta * vf, 3) + ' in' }, series: V - best.Vmix > 1e-6 * V ? [{ kind: 'tube', tag: 'plug ' + fmt(V - best.Vmix, 3) }] : [] }
          : { first: { kind: 'loop', tag: 'V = ' + fmt(V, 3), loop: 'R = ' + fmt(best.R, 3), flow: fmt(best.beta * vf, 3) + ' in' }, series: [] };
        const diagram = BIO.sketch(Object.assign(sk, { feed: 'v = ' + fmt(vf, 3), bypass: best.beta < 0.999 ? 'bypass ' + fmt((1 - best.beta) * vf, 3) : null, out: ax.lab + ' = ' + fmt(ax.f(best.d), 3) }));
        return { kpi, warn, info, plots: [curve(M, dO, undefined, best.d)], tables, diagram };
      }
      /* smallest volume for a target */
      const tq = v.tq, tv = +v.tval;
      const dT = tq === 'CA' ? M.CA0 - tv : tq === 'CC' ? (tv - M.CC0) / M.yCA : (tv - M.CR0) / M.yRA;
      if (!(dT > 0 && dT < M.Dmax)) throw new Error('The target lies outside what the feed can reach (0 < C_A0 − C_A < ' + fmt(M.Dmax, 4) + ').');
      const Vfree = dT <= dO ? vf * (dT / dO) * tauO : vf * tauO + vf * M.tPlug(dO, dT);
      if (v.tap) {
        kpi.unshift({ k: 'V', l: 'Smallest volume (side tap, the book’s rule)', v: Vfree, main: true, s: dT <= dO ? 'mixed flow at the maximum rate on ' + fmt(dT / dO * 100, 3) + '% of the feed, the rest bypassed' : 'recycle loop of ' + fmt(vf * tauO, 4) + ' at the maximum rate, then ' + fmt(vf * M.tPlug(dO, dT), 4) + ' in plug flow' });
      }
      /* no side tap: recycle over the whole vessel, optimize R (and bypass) */
      const Vrec = (R, beta) => {
        const dOut = dT / beta;
        if (!(dOut < M.Dmax)) return Infinity;
        const d1 = R * dOut / (R + 1);
        if (M.rC(d1) <= 0) return Infinity;
        return (R + 1) * beta * vf * M.tPlug(d1, dOut);
      };
      let bv = Infinity, bR = 0, bB = 1;
      N.logspace(1e-3, 1e4, 120).forEach(R => N.linspace(Math.max(dT / M.Dmax * 1.0001, 0.02), 1, 30).forEach(be => { const q = Vrec(R, be); if (q < bv) { bv = q; bR = R; bB = be; } }));
      {
        const g = z => Vrec(Math.pow(10, z[0]), Math.min(Math.max(z[1], 1e-4), 1));
        const r = N.nelderMead(g, [Math.log10(Math.max(bR, 1e-3)), bB], { tol: 1e-12, maxIt: 300, step: [0.05, 0.02] });
        if (r.f < bv) { bv = r.f; bR = Math.pow(10, r.x[0]); bB = Math.min(Math.max(r.x[1], 1e-4), 1); }
      }
      kpi.push({ k: 'Vrec', l: 'Smallest volume without a side tap (one recycle ratio)', v: bv, main: !v.tap, s: 'R = ' + fmt(bR, 3) + (bB < 0.999 ? ', ' + fmt((1 - bB) * 100, 3) + '% bypassed' : '') });
      kpi.push({ k: 'Vmixed', l: 'One mixed tank, no bypass', v: vf * dT * M.yCA / M.rC(dT) });
      if (M.CC0 > 0) kpi.push({ k: 'Vplug', l: 'Plain plug flow', v: vf * M.tPlug(0, dT) });
      if (!v.tap) info.push('With a side tap the same target needs ' + fmt(Vfree, 4) + '.');
      const ax = BIO.axis(M);
      let sk;
      if (v.tap) sk = dT <= dO ? { first: { kind: 'loop', tag: 'V = ' + fmt(Vfree, 3), flow: fmt(dT / dO * vf, 3) + ' in' }, bypass: 'bypass ' + fmt((1 - dT / dO) * vf, 3) } : { first: { kind: 'loop', tag: 'loop ' + fmt(vf * tauO, 3) }, series: [{ kind: 'tube', tag: 'plug ' + fmt(vf * M.tPlug(dO, dT), 3) }] };
      else sk = { first: { kind: 'loop', tag: 'V = ' + fmt(bv, 3), loop: 'R = ' + fmt(bR, 3), flow: fmt(bB * vf, 3) + ' in' }, bypass: bB < 0.999 ? 'bypass ' + fmt((1 - bB) * vf, 3) : null };
      const diagram = BIO.sketch(Object.assign(sk, { feed: 'v = ' + fmt(vf, 3), out: ax.lab + ' = ' + fmt(ax.f(dT), 3) }));
      return { kpi, warn, info, plots: [curve(M, dO, dT)], tables, diagram };
    }
  });
  function leave(kpi, o, M) {
    const prod = isFinite(M.CRs) && M.yRA > 0;
    const k = [{ k: 'CA', l: prod ? 'C_A leaving' : 'Lowest C_A leaving', v: o.CA, main: !prod }, { k: 'CC', l: 'C_C leaving', v: o.CC }];
    if (M.yRA > 0) k.push({ k: 'CR', l: prod ? 'Highest C_R leaving' : 'C_R leaving', v: o.CR, main: prod });
    if (prod) k.unshift(k.pop());
    kpi.unshift(...k);
  }
  function curve(M, dO, dT, dOut) {
    const ax = BIO.axis(M);
    const dd = N.linspace(M.Dmax * 1e-4, M.Dmax * 0.999, 240);
    const marks = [{ t: 'v', x: ax.f(dO), color: 2, label: 'maximum rate' }];
    if (dT !== undefined) marks.push({ t: 'v', x: ax.f(dT), color: 8, label: 'target' });
    if (dOut !== undefined) marks.push({ t: 'v', x: ax.f(dOut), color: 8, label: 'leaving' });
    return { title: '1/(−r_A) against ' + ax.lab + ' (Fig. 29.6)', x: { label: ax.lab, min: 0, max: ax.max }, y: { label: '1/(−r_A)', log: true }, series: [{ name: '1/(−r_A)', x: dd.map(ax.f), y: dd.map(d => M.yCA / M.rC(d)), color: 1, w: 2.2 }], marks, cap: 'Between the minimum and C_A0, mixed flow with bypass; below the minimum, mixed flow to it and plug flow beyond.' };
  }
})(typeof window !== 'undefined' ? window : globalThis);
