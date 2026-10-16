/* Chapter 3 — Interpretation of batch reactor data; Chapter 4 — Introduction to reactor design */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, KIT = ST.KIT, D = ST.D;
  const RL = 0.082057;
  const PU = { atm: 1, mmHg: 1 / 760, kPa: 1 / 101.325, Pa: 1 / 101325, bar: 1 / 1.01325 };

  /* ---- convert the measured column into C_A(t) ---- */
  function batchData(v) {
    const rows = (v.data || []).filter(r => isFinite(r[0]) && isFinite(r[1])).map(r => [+r[0], +r[1]]).sort((a, b) => a[0] - b[0]);
    if (rows.length < 3) throw new Error('Enter at least three (t, value) pairs.');
    let CA0 = +v.CA0, unit = v._cu, t = rows.map(r => r[0]), C;
    const y = rows.map(r => r[1]);
    switch (v.dtype) {
      case 'CA': C = y.slice(); if (v.CA0from === 'data') CA0 = rows[0][0] === 0 ? y[0] : CA0; break;
      case 'XA': C = y.map(x => CA0 * (1 - x)); break;
      case 'CR': C = y.map(c => CA0 - (+v.aA / +v.rR) * (c - (+v.CR0 || 0))); break;
      case 'lam': { const l0 = +v.lam0, li = +v.lamInf; C = y.map(l => CA0 * (1 - (l - l0) / (li - l0))); break; }
      case 'pA': case 'pi': {
        let p = y.slice(), p0 = +v.pA0;
        if (v.dtype === 'pi') { const pi0 = +v.pi0; p = y.map(P => p0 - (+v.aA / +v.dn) * (P - pi0)); }
        if (v.toC) { const f = PU[v.pu] / (RL * (+v.T)); C = p.map(q => q * f); CA0 = p0 * f; unit = 'mol/L'; }
        else { C = p; CA0 = p0; unit = v.pu; }
        break;
      }
      default: throw new Error('Unknown data type');
    }
    if (v.dtype !== 'CA' && v.dtype !== 'pA' && v.dtype !== 'pi') { /* keep given CA0 */ }
    const pts = t.map((ti, i) => ({ t: ti, C: C[i] }));
    return { pts, CA0, unit };
  }

  /* ---- candidate rate forms for the integral method (constant volume) ---- */
  const lg = x => Math.log(x);
  const MODELS = {
    n0: { name: 'Zero order', tex: '-r_A = k', p: ['k'], tf: ['log'],
      C: (t, p, c) => Math.max(c.CA0 - p[0] * t, 0),
      lin: (pt, p, c) => ({ x: pt.t, y: c.CA0 - pt.C }), linLab: 'C_{A0} − C_A', slope: p => p[0], eq: 'Eq. 3.31' },
    n1: { name: 'First order', tex: '-r_A = kC_A', p: ['k'], tf: ['log'],
      C: (t, p, c) => c.CA0 * Math.exp(-p[0] * t),
      lin: (pt, p, c) => ({ x: pt.t, y: lg(c.CA0 / pt.C) }), linLab: 'ln(C_{A0}/C_A)', slope: p => p[0], eq: 'Eq. 3.12' },
    n2: { name: 'Second order (2A, or A+B with M = 1)', tex: '-r_A = kC_A^2', p: ['k'], tf: ['log'],
      C: (t, p, c) => c.CA0 / (1 + c.CA0 * p[0] * t),
      lin: (pt, p, c) => ({ x: pt.t, y: 1 / pt.C - 1 / c.CA0 }), linLab: '1/C_A − 1/C_{A0}', slope: p => p[0], eq: 'Eq. 3.16' },
    n2AB: { name: 'Second order A + bB, M ≠ b', tex: '-r_A = kC_AC_B', p: ['k'], tf: ['log'], needB: true,
      C: (t, p, c) => {
        const M = c.M, b = c.b;
        if (Math.abs(M - b) < 1e-9) { const z = b * p[0] * c.CA0 * t; return c.CA0 / (1 + z); }
        const E = Math.exp(c.CA0 * (M - b) * p[0] * t), X = M * (E - 1) / (M * E - b);
        return c.CA0 * (1 - X);
      },
      lin: (pt, p, c) => { const X = 1 - pt.C / c.CA0; return { x: pt.t, y: lg((c.M - c.b * X) / (c.M * (1 - X))) }; },
      linLab: 'ln[(M − bX_A)/M(1 − X_A)]', slope: (p, c) => c.CA0 * (c.M - c.b) * p[0], eq: 'Eq. 3.14, 3.18' },
    n3: { name: 'Third order (3A, or equimolar)', tex: '-r_A = kC_A^3', p: ['k'], tf: ['log'],
      C: (t, p, c) => c.CA0 / Math.sqrt(1 + 2 * p[0] * c.CA0 * c.CA0 * t),
      lin: (pt, p, c) => ({ x: pt.t, y: 1 / (pt.C * pt.C) - 1 / (c.CA0 * c.CA0) }), linLab: '1/C_A² − 1/C_{A0}²', slope: p => 2 * p[0], eq: 'Eq. 3.27 form' },
    nth: { name: 'nth order', tex: '-r_A = kC_A^n', p: ['k', 'n'], tf: ['log', 'none'],
      C: (t, p, c) => {
        const k = p[0], n = p[1];
        if (Math.abs(n - 1) < 1e-6) return c.CA0 * Math.exp(-k * t);
        const base = Math.pow(c.CA0, 1 - n) - (1 - n) * k * t;
        if (base <= 0) return 0;
        return Math.pow(base, 1 / (1 - n));
      },
      lin: (pt, p, c) => ({ x: pt.t, y: Math.abs(p[1] - 1) < 1e-6 ? lg(c.CA0 / pt.C) : (Math.pow(pt.C, 1 - p[1]) - Math.pow(c.CA0, 1 - p[1])) / (p[1] - 1) }),
      linLab: '(C_A^{1−n} − C_{A0}^{1−n})/(n − 1)', slope: p => p[0], eq: 'Eq. 3.29' },
    rev1: { name: 'Reversible first order A ⇌ R', tex: '-r_A = k_1C_A - k_2C_R', p: ['k1', 'XAe'], tf: ['log', 'logit'],
      C: (t, p, c) => { const Xe = p[1], X = Xe * (1 - Math.exp(-(c.MR + 1) / (c.MR + Xe) * p[0] * t)); return c.CA0 * (1 - X); },
      lin: (pt, p, c) => { const X = 1 - pt.C / c.CA0; return { x: pt.t, y: -lg(1 - X / p[1]) }; },
      linLab: '−ln(1 − X_A/X_{Ae})', slope: (p, c) => (c.MR + 1) / (c.MR + p[1]) * p[0], eq: 'Eq. 3.54' },
    rev2: { name: 'Reversible second order, C_{A0} = C_{B0}', tex: '-r_A = k_1C_AC_B - k_2C_RC_S', p: ['k1', 'XAe'], tf: ['log', 'logit'],
      C: (t, p, c) => {
        const Xe = p[1], Gx = Math.exp(2 * p[0] * (1 / Xe - 1) * c.CA0 * t);
        const X = Xe * (1 - Gx) / (2 * Xe - 1 - Gx);
        return c.CA0 * (1 - X);
      },
      lin: (pt, p, c) => { const X = 1 - pt.C / c.CA0, Xe = p[1]; return { x: pt.t, y: lg((Xe - (2 * Xe - 1) * X) / (Xe - X)) }; },
      linLab: 'ln[(X_{Ae} − (2X_{Ae} − 1)X_A)/(X_{Ae} − X_A)]', slope: (p, c) => 2 * p[0] * (1 / p[1] - 1) * c.CA0, eq: 'Eq. 3.56' },
    auto: { name: 'Autocatalytic A + R → 2R', tex: '-r_A = kC_AC_R', p: ['k', 'CR0'], tf: ['log', 'log'],
      C: (t, p, c) => { const C0 = c.CA0 + p[1]; return C0 / (1 + (p[1] / c.CA0) * Math.exp(C0 * p[0] * t)); },
      lin: (pt, p, c) => { const CR = c.CA0 + p[1] - pt.C; return { x: pt.t, y: lg((CR / p[1]) / (pt.C / c.CA0)) }; },
      linLab: 'ln[(C_R/C_{R0})/(C_A/C_{A0})]', slope: (p, c) => (c.CA0 + p[1]) * p[0], eq: 'Eq. 3.42' },
    sat: { name: 'Shifting order k₁C_A/(1 + k₂C_A)', tex: '-r_A = \\frac{k_1C_A}{1+k_2C_A}', p: ['k1', 'k2'], tf: ['log', 'log'],
      C: (t, p, c) => { const k1 = p[0], k2 = p[1]; if (k2 < 1e-12) return c.CA0 * Math.exp(-k1 * t); return N.lambertW(k2 * c.CA0 * Math.exp(k2 * c.CA0 - k1 * t)) / k2; },
      lin: (pt, p, c) => ({ x: pt.t / (c.CA0 - pt.C), y: lg(c.CA0 / pt.C) / (c.CA0 - pt.C) }), linLab: 'ln(C_{A0}/C_A)/(C_{A0} − C_A)', linX: 't/(C_{A0} − C_A)', slope: p => p[0], icept: p => -p[1], eq: 'Eq. 3.58c' }
  };
  const TF = { log: [Math.exp, Math.log], logit: [u => 1 / (1 + Math.exp(-u)), x => Math.log(x / (1 - x))], none: [u => u, x => x] };

  function guess(id, pts, c) {
    const last = pts[pts.length - 1], tl = last.t || 1, Cl = Math.max(last.C, 1e-9 * c.CA0);
    const k1 = Math.max(lg(c.CA0 / Cl) / tl, 1e-9);
    switch (id) {
      case 'n0': return [Math.max((c.CA0 - Cl) / tl, 1e-12)];
      case 'n1': return [k1];
      case 'n2': return [Math.max((1 / Cl - 1 / c.CA0) / tl, 1e-12)];
      case 'n2AB': return [k1 / Math.max(c.CA0 * c.M, 1e-9)];
      case 'n3': return [Math.max((1 / (Cl * Cl) - 1 / (c.CA0 * c.CA0)) / (2 * tl), 1e-12)];
      case 'nth': return [k1 / Math.pow(c.CA0, 0.0), 1.2];
      case 'rev1': { const Xe = Math.min(0.995, Math.max(0.05, (1 - Cl / c.CA0) * 1.05)); return [k1, Xe]; }
      case 'rev2': { const Xe = Math.min(0.995, Math.max(0.55, (1 - Cl / c.CA0) * 1.05)); return [k1 / c.CA0, Xe]; }
      case 'auto': return [k1 / c.CA0 * 2, Math.max(c.CR0 || c.CA0 * 0.02, c.CA0 * 1e-4)];
      case 'sat': return [k1 * 1.5, 0.5 / c.CA0];
    }
    return [k1];
  }
  function fitModel(id, pts, c, fixCR0) {
    const M = MODELS[id];
    const tf = M.tf.map(x => TF[x]);
    let p0 = guess(id, pts, c);
    if (id === 'auto' && fixCR0) p0[1] = c.CR0;
    const u0 = p0.map((pv, i) => tf[i][1](pv));
    const fixed = id === 'auto' && fixCR0 ? [false, true] : [];
    const scale = Math.max(...pts.map(p => Math.abs(p.C)), 1e-30);
    const resid = u => { const p = u.map((x, i) => tf[i][0](x)); return pts.map(pt => (M.C(pt.t, p, c) - pt.C) / scale); };
    let best = null;
    const starts = [u0];
    if (id === 'nth') [0.5, 2, 3].forEach(n => { const g = u0.slice(); g[1] = n; starts.push(g); });
    if (id === 'sat') [0.05, 5].forEach(f => { const g = u0.slice(); g[1] = Math.log(f / c.CA0); starts.push(g); });
    starts.forEach(s => {
      try {
        const r = N.lm(resid, s, { fixed });
        if (r.ok && isFinite(r.ssr) && (!best || r.ssr < best.ssr)) best = r;
      } catch (e) { /* try next start */ }
    });
    if (!best) return null;
    const p = best.p.map((x, i) => tf[i][0](x));
    const se = best.se.map((s, i) => (M.tf[i] === 'log' ? p[i] * s : M.tf[i] === 'logit' ? p[i] * (1 - p[i]) * s : s));
    const fitC = pts.map(pt => M.C(pt.t, p, c));
    const ssr = N.sum(pts.map((pt, i) => (pt.C - fitC[i]) ** 2));
    const n = pts.length, np = M.p.length - (fixed[1] ? 1 : 0);
    const aic = n * Math.log(Math.max(ssr, 1e-300) / n) + 2 * np + (n - np - 1 > 0 ? 2 * np * (np + 1) / (n - np - 1) : 0);
    return { id, p, se, ssr, r2: N.r2(pts.map(q => q.C), fitC), aic, fitC, np };
  }

  const dataInputs = (d = {}) => [
    { group: 'Data' },
    { k: 'dtype', t: 'select', label: 'Measured quantity', opts: [['CA', 'concentration of A'], ['XA', 'conversion of A'], ['CR', 'concentration of product R'], ['pA', 'partial pressure of A'], ['pi', 'total pressure (constant volume)'], ['lam', 'a physical property linear in conversion']], v: d.dtype || 'CA' },
    { k: 'data', t: 'table', label: 'Run', cols: ['t', 'measured'], v: d.data || [[0, 10], [20, 8], [40, 6], [60, 5], [120, 3], [180, 2], [300, 1]] },
    { k: 'CA0', label: '$C_{A0}$', unit: '[C]', v: d.CA0 || 10, min: 1e-9, max: 1e9, scale: 'log', if: v => !(v.dtype === 'pA' || v.dtype === 'pi') },
    { k: 'aA', label: 'Moles of A in the equation, a', v: 1, min: 0.1, max: 10, if: v => v.dtype === 'CR' || v.dtype === 'pi' },
    { k: 'rR', label: 'Moles of R in the equation, r', v: 1, min: 0.1, max: 10, if: v => v.dtype === 'CR' },
    { k: 'CR0', label: '$C_{R0}$', unit: '[C]', v: 0, min: 0, max: 1e9, if: v => v.dtype === 'CR' },
    { k: 'lam0', label: 'Property at t = 0, λ₀', v: 0, min: -1e9, max: 1e9, slider: false, if: v => v.dtype === 'lam' },
    { k: 'lamInf', label: 'Property at completion, λ∞', v: 1, min: -1e9, max: 1e9, slider: false, if: v => v.dtype === 'lam' },
    { k: 'dn', label: 'Change in moles in the equation, Δn', v: -1, min: -10, max: 10, step: 0.5, if: v => v.dtype === 'pi', hint: 'Δn = r + s + … − a − b − …; 2A → B gives −1.' },
    { k: 'pi0', label: 'Initial total pressure π₀', v: 1, min: 1e-9, max: 1e9, slider: false, if: v => v.dtype === 'pi' },
    { k: 'pA0', label: 'Initial partial pressure of A, p_{A0}', v: 1, min: 1e-9, max: 1e9, slider: false, if: v => v.dtype === 'pA' || v.dtype === 'pi' },
    { k: 'pu', t: 'select', label: 'Pressure unit', opts: Object.keys(PU).map(k => [k, k]), v: 'atm', if: v => v.dtype === 'pA' || v.dtype === 'pi' },
    { k: 'toC', t: 'check', label: 'Convert pressures to concentrations ($p = CRT$)', v: true, if: v => v.dtype === 'pA' || v.dtype === 'pi' },
    { k: 'T', label: 'Temperature', unit: 'K', v: 373.15, min: 1, max: 3000, if: v => (v.dtype === 'pA' || v.dtype === 'pi') && v.toC }
  ];

  /* ---------------- integral method ---------------- */
  ST.tool({
    id: 'integral', tag: 'H-301', ch: 3, sec: '3.1', title: 'Batch data: the integral method',
    lede: 'Test ten rate forms against one batch run at once. Each is fitted by nonlinear least squares and checked on its own linearized plot, then ranked by a corrected Akaike criterion.',
    keys: 'integral method of analysis batch kinetics first order second order zero order nth order reversible autocatalytic michaelis shifting order total pressure conversion fit rate equation',
    units: { C: 'mol/L', t: 's' }, slow: true, sweep: false,
    inputs: dataInputs().concat([
      { group: 'Rate forms' },
      { k: 'MB', label: 'For A + bB: $M = C_{B0}/C_{A0}$', v: 2, min: 0, max: 100, step: 0.1 },
      { k: 'b', label: 'b (moles B per mole A)', v: 1, min: 0.1, max: 5, step: 0.1 },
      { k: 'MR', label: 'For A ⇌ R: $C_{R0}/C_{A0}$', v: 0, min: 0, max: 10, step: 0.01 },
      { k: 'CR0a', label: 'Autocatalytic seed $C_{R0}$ (0 = fit it)', unit: '[C]', v: 0, min: 0, max: 1e6 },
      { k: 'show', t: 'select', label: 'Linear test to plot', opts: [['best', 'best ranked']].concat(Object.keys(MODELS).map(k => [k, MODELS[k].name])), v: 'best' }
    ]),
    presets: [
      { name: 'Example 3.1', v: { dtype: 'CA', data: [[0, 10], [20, 8], [40, 6], [60, 5], [120, 3], [180, 2], [300, 1]], CA0: 10, _cu: 'mol/L', _tu: 's', show: 'nth' }, book: { n: { v: 1.4, tol: 0.06 }, k: { v: 0.005, tol: 0.15 } }, note: 'First- and second-order straight lines both fail; the book finds n = 1.4 and k = 0.005 by the fractional-life method.' },
      { name: 'Problem 3.11', v: { dtype: 'CA', data: [[0, 1000], [100, 500], [200, 333], [300, 250], [400, 200]], CA0: 1000, _cu: 'mol/m³', _tu: 'min', show: 'n2' }, note: 'A → R in water. Which rate form fits, and what conversion follows after 5 h from 500 mol/m³?' },
      { name: 'Problem 3.15', v: { dtype: 'CA', data: [[0, 1.0], [1, 0.84], [2, 0.68], [3, 0.53], [4, 0.38], [5, 0.27], [6, 0.16], [7, 0.09], [8, 0.04], [9, 0.018], [10, 0.006], [11, 0.0025]], CA0: 1, _cu: 'mmol/L', _tu: 'h', show: 'sat' }, note: 'Sucrose hydrolysis by sucrase, $C_{E0}$ = 0.01 mmol/L. Michaelis–Menten is the shifting-order form with $k_3C_{E0}$ = k₁/k₂ and $C_M$ = 1/k₂.' },
      { name: 'Problem 3.20', v: { dtype: 'CR', data: [[0, 0], [41, 1.18], [48, 1.38], [55, 1.63], [75, 2.24], [96, 2.75], [127, 3.31], [146, 3.76], [162, 3.81], [180, 4.11], [194, 4.31], [212, 4.45], [267, 4.86], [318, 5.15], [368, 5.32], [379, 5.35], [410, 5.42]], CA0: 5.5, aA: 1, rR: 2, CR0: 0, _cu: 'mol/L', _tu: 'min', show: 'rev2' }, note: 'H₂SO₄ + (C₂H₅)₂SO₄ ⇌ 2C₂H₅SO₄H at 22.9 °C, both reactants at 5.5 mol/L; the product levels off near 5.80 mol/L.' },
      { name: 'Problem 3.21', v: { dtype: 'pi', data: [[1, 1.14], [2, 1.04], [3, 0.982], [4, 0.940], [5, 0.905], [6, 0.870], [7, 0.850], [8, 0.832], [9, 0.815], [10, 0.800], [15, 0.754], [20, 0.728]], aA: 2, dn: -1, pi0: 373.15 / 298.15, pA0: 373.15 / 298.15, pu: 'atm', toC: true, T: 373.15, _tu: 'min' }, note: 'Pure A, filled at 1 atm and 25 °C, heated to 100 °C; 2A → B. π₀ at 100 °C is taken as 1.252 atm.' },
      { name: 'Problem 3.25', v: { dtype: 'pA', data: [[0, 760], [2, 600], [4, 475], [6, 390], [8, 320], [10, 275], [12, 240], [14, 215]], pA0: 760, pu: 'mmHg', toC: true, T: 273.15, _tu: 'min', show: 'rev1', MR: 0 }, note: 'Pure gaseous A at 0 °C, A → 2.5R; at very long times p_A falls only to 150 mm Hg, so try the reversible forms.' }
    ],
    eqs: [
      ['-\\ln\\frac{C_A}{C_{A0}} = kt', 'Eq. 3.11–3.12'],
      ['\\frac{1}{C_A}-\\frac{1}{C_{A0}} = kt', 'Eq. 3.16'],
      ['\\ln\\frac{C_BC_{A0}}{C_{B0}C_A} = C_{A0}(M-1)kt', 'Eq. 3.14'],
      ['C_A^{1-n}-C_{A0}^{1-n} = (n-1)kt', 'Eq. 3.29'],
      ['-\\ln\\left(1-\\frac{X_A}{X_{Ae}}\\right) = \\frac{M+1}{M+X_{Ae}}k_1t', 'Eq. 3.54'],
      ['\\ln\\frac{X_{Ae}-(2X_{Ae}-1)X_A}{X_{Ae}-X_A} = 2k_1\\left(\\frac{1}{X_{Ae}}-1\\right)C_{A0}t', 'Eq. 3.56'],
      ['\\ln\\frac{C_{A0}(C_0-C_A)}{C_A(C_0-C_{A0})} = C_0kt,\\quad C_0 = C_{A0}+C_{R0}', 'Eq. 3.42'],
      ['\\ln\\frac{C_{A0}}{C_A}+k_2(C_{A0}-C_A) = k_1t', 'Eq. 3.58a'],
      ['p_A = p_{A0}-\\frac{a}{\\Delta n}(\\pi-\\pi_0)', 'Eq. 3.5', 'Total-pressure data from a constant-volume run.']
    ],
    method: ['Each candidate is fitted to the concentrations themselves (not to its linearized form) by Levenberg–Marquardt, with rate constants kept positive and equilibrium conversions between 0 and 1. The book’s straight-line test is then drawn with the fitted constants, so curvature you can see is curvature the data really have.', 'Ranking uses the small-sample Akaike criterion AICc = n ln(SSR/n) + 2p + 2p(p+1)/(n−p−1), which rewards fit and penalizes extra constants. Differences under about 2 are not meaningful.', 'The shifting-order form integrates exactly through the Lambert W function: $C_A = W(k_2C_{A0}e^{k_2C_{A0}-k_1t})/k_2$.'],
    run(v) {
      const d = batchData(v);
      const pts = d.pts.filter(p => p.C > 0);
      const c = { CA0: d.CA0, M: +v.MB, b: +v.b, MR: +v.MR, CR0: +v.CR0a > 0 ? +v.CR0a : d.CA0 * 0.02 };
      const warn = [];
      if (pts.length < 3) throw new Error('Need at least three positive concentrations.');
      const fits = [];
      Object.keys(MODELS).forEach(id => {
        if (id === 'rev2' && v.dtype !== 'CR' && !(pts[pts.length - 1].C > 0.05 * d.CA0)) { /* still try */ }
        const f = fitModel(id, pts, c, +v.CR0a > 0);
        if (f) fits.push(f);
      });
      fits.sort((a, b) => a.aic - b.aic);
      const best = fits[0];
      const pick = v.show === 'best' ? best : fits.find(f => f.id === v.show) || best;
      const M = MODELS[pick.id];
      const pname = (id, i) => MODELS[id].p[i].replace('XAe', 'X_{Ae}').replace('CR0', 'C_{R0}').replace('k1', 'k_1').replace('k2', 'k_2');
      const kpi = [{ k: 'model', l: 'Best by AICc', v: MODELS[best.id].name, main: true }];
      pick.p.forEach((p, i) => kpi.push({ k: MODELS[pick.id].p[i], l: pname(pick.id, i) + (pick === best ? '' : ' (' + M.name + ')'), v: p, s: isFinite(pick.se[i]) ? '± ' + N.fmt(1.96 * pick.se[i], 2) + ' (95%)' : '' }));
      kpi.push({ k: 'r2', l: 'R² on C_A', v: pick.r2 });
      if (pick.id === 'sat') kpi.push({ k: 'CM', l: 'Michaelis C_M = 1/k₂', v: 1 / pick.p[1], u: d.unit });
      const tmax = Math.max(...pts.map(p => p.t));
      const ts = N.linspace(0, tmax * 1.05, 160);
      const series = [{ name: 'data', type: 'pts', x: pts.map(p => p.t), y: pts.map(p => p.C), color: 8 }];
      fits.slice(0, 4).forEach((f, i) => series.push({ name: MODELS[f.id].name, x: ts, y: ts.map(t => MODELS[f.id].C(t, f.p, c)), color: i + 1, w: i === 0 ? 2.6 : 1.6, dash: i ? '6 4' : undefined }));
      if (fits.slice(0, 4).indexOf(pick) < 0) series.push({ name: M.name, x: ts, y: ts.map(t => M.C(t, pick.p, c)), color: 5 });
      const lin = pts.map(pt => M.lin(pt, pick.p, c)).filter(q => isFinite(q.x) && isFinite(q.y));
      const sl = M.slope(pick.p, c), ic = M.icept ? M.icept(pick.p, c) : 0;
      const lx = N.linspace(0, Math.max(...lin.map(q => q.x)) * 1.05, 20);
      const lr = N.linreg(lin.map(q => q.x), lin.map(q => q.y), { origin: !M.icept });
      const linPlot = { title: 'Straight-line test: ' + M.name + ' (' + M.eq + ')', x: { label: M.linX || 't', unit: M.linX ? '' : v._tu, min: 0 }, y: { label: M.linLab }, series: [{ name: 'data', type: 'pts', x: lin.map(q => q.x), y: lin.map(q => q.y), color: 2 }, { name: 'fitted constants', x: lx, y: lx.map(x => ic + sl * x), color: 1 }], cap: 'Least-squares line on these axes: slope ' + N.fmt(lr.b, 4) + (M.icept ? ', intercept ' + N.fmt(lr.a, 4) : ' (through the origin)') + ', R² = ' + N.fmt(lr.r2, 4) + '.' };
      const res = { title: 'Residuals, ' + M.name, x: { label: 't', unit: v._tu, min: 0 }, y: { label: 'C_A − fit', unit: d.unit }, series: [{ name: 'residual', type: 'bar', x: pts.map(p => p.t), y: pts.map((p, i) => p.C - pick.fitC[i]), color: 2 }], marks: [{ t: 'h', y: 0, color: 8, dash: false }] };
      const rank = { title: 'All rate forms, best first', cols: ['rate form', 'constants', 'R²', 'SSR', 'AICc', 'ΔAICc'], rows: fits.map(f => [MODELS[f.id].name, f.p.map((p, i) => MODELS[f.id].p[i] + ' = ' + N.fmt(p, 4)).join(', '), f.r2, f.ssr, f.aic, f.aic - best.aic]) };
      /* k from individual points for one-constant forms */
      const tabs = [rank];
      if (M.p.length === 1 || pick.id === 'nth' || pick.id === 'rev1' || pick.id === 'rev2' || pick.id === 'auto') {
        const rows = pts.filter(p => p.t > 0).map(p => {
          const f = k => M.C(p.t, [k].concat(pick.p.slice(1)), c) - p.C;
          const kk = N.root(f, pick.p[0] * 1e-3, pick.p[0] * 1e3, { log: true, n: 120 });
          return [p.t, p.C, kk];
        });
        tabs.push({ title: 'k from each point, other constants held (Fig. 3.25)', cols: ['t', 'C_A', MODELS[pick.id].p[0]], rows });
      }
      if (d.unit !== v._cu && (v.dtype === 'pA' || v.dtype === 'pi')) warn.push('Concentrations are in ' + d.unit + '.');
      if (best.id === 'n0' && pts.some(p => p.C <= 1e-9)) warn.push('Zero order ends abruptly when A runs out; check the last points.');
      return { kpi, warn, plots: [{ title: 'Data and the best rate forms', wide: true, x: { label: 't', unit: v._tu, min: 0 }, y: { label: 'C_A', unit: d.unit, min: 0 }, series }, linPlot, res], tables: tabs };
    }
  });

  /* ---------------- differential method ---------------- */
  ST.tool({
    id: 'differential', tag: 'H-302', ch: 3, sec: '3.1', title: 'Batch data: the differential method',
    lede: 'Estimate rates from a concentration–time curve by several smoothing choices, then fit nth-order or saturation kinetics; or start from tabulated rates, including runs varying several concentrations at once.',
    keys: 'differential method rate from slope tangent derivative log log plot order initial rates multiple regression lineweaver burk graphical integration batch time from rate data',
    units: { C: 'mol/L', t: 's' }, slow: true, sweep: false,
    inputs: [
      { group: 'Data' },
      { k: 'kind', t: 'select', label: 'What you have', opts: [['ct', 'C_A against t'], ['rate', 'rates at several C_A'], ['multi', 'rates at several C_A, C_B, C_C']], v: 'ct' },
      { k: 'data', t: 'table', label: 'C_A against t', cols: ['t', 'C_A'], v: [[0, 10], [20, 8], [40, 6], [60, 5], [120, 3], [180, 2], [300, 1]], if: v => v.kind === 'ct' },
      { k: 'rdata', t: 'table', label: 'Rates', cols: ['C_A', '−r_A'], v: [[1, 0.06], [2, 0.1], [4, 0.25], [6, 1.0], [7, 2.0], [9, 1.0], [12, 0.5]], if: v => v.kind === 'rate' },
      { k: 'mdata', t: 'table', label: 'Rate runs', cols: ['C_A', 'C_B', 'C_C', 'rate'], v: [[1, 3, 0.02, 9], [3, 1, 0.02, 5], [4, 4, 0.04, 32], [2, 2, 0.01, 6], [2, 4, 0.03, 20], [1, 2, 0.05, 12]], if: v => v.kind === 'multi' },
      { k: 'how', t: 'select', label: 'Slopes from', opts: [['local', 'local quadratic through neighbours'], ['logpoly', 'polynomial in ln C_A (smooth decay)'], ['poly', 'global polynomial'], ['central', 'finite differences']], v: 'logpoly', if: v => v.kind === 'ct' },
      { k: 'deg', label: 'Polynomial degree', v: 3, min: 1, max: 6, step: 1, if: v => v.kind === 'ct' && (v.how === 'poly' || v.how === 'logpoly') },
      { k: 'half', label: 'Neighbours each side', v: 2, min: 1, max: 5, step: 1, if: v => v.kind === 'ct' && v.how === 'local' },
      { group: 'Fit and use', if: v => v.kind !== 'multi' },
      { k: 'form', t: 'select', label: 'Rate form', opts: [['n', 'nth order: log(−r) against log C'], ['sat', 'k₁C/(1+k₂C): both linear tests']], v: 'n', if: v => v.kind !== 'multi' },
      { k: 'C1', label: 'Batch time from C_A =', unit: '[C]', v: 10, min: 0, max: 1e9, slider: false, if: v => v.kind === 'rate' },
      { k: 'C2', label: 'down to C_A =', unit: '[C]', v: 2, min: 0, max: 1e9, slider: false, if: v => v.kind === 'rate' }
    ],
    presets: [
      { name: 'Example 3.2', v: { kind: 'ct', data: [[0, 10], [20, 8], [40, 6], [60, 5], [120, 3], [180, 2], [300, 1]], how: 'logpoly', deg: 3, form: 'n', _cu: 'mol/L', _tu: 's' }, book: { n: { v: 1.43, tol: 0.08 }, k: { v: 0.005, tol: 0.25 } }, note: 'The book draws tangents by eye and gets n = 1.43, k = 0.005. Try each smoothing choice to see how much the slopes move.' },
      { name: 'Problem 3.24', v: { kind: 'rate', rdata: [[1, 0.06], [2, 0.1], [4, 0.25], [6, 1.0], [7, 2.0], [9, 1.0], [12, 0.5]], C1: 10, C2: 2, _cu: 'mol/L', _tu: 'h' }, note: 'Rates measured with a homogeneous catalyst; batch time to go from 10 to 2 mol/L by integrating 1/(−r_A).' },
      { name: 'Problem 3.28', v: { kind: 'multi', mdata: [[1, 3, 0.02, 9], [3, 1, 0.02, 5], [4, 4, 0.04, 32], [2, 2, 0.01, 6], [2, 4, 0.03, 20], [1, 2, 0.05, 12]] }, note: 'Rate of A + B → AB catalysed by C: regress log r on log [A], log [B] and log [C].' },
      { name: 'Problem 3.16', v: { kind: 'ct', data: [[0, 1.0], [1, 0.84], [2, 0.68], [3, 0.53], [4, 0.38], [5, 0.27], [6, 0.16], [7, 0.09], [8, 0.04], [9, 0.018], [10, 0.006]], how: 'local', half: 2, form: 'sat', _cu: 'mmol/L', _tu: 'h' }, note: 'Sucrose data again, now by the differential method.' }
    ],
    eqs: [
      ['-r_A = -\\frac{dC_A}{dt},\\qquad \\log(-r_A) = \\log k + n\\log C_A', 'Fig. 3.18'],
      ['\\frac{1}{-r_A} = \\frac{1}{k_1C_A}+\\frac{k_2}{k_1}', 'Eq. 3.61'],
      ['-r_A = \\frac{k_1}{k_2} - \\frac{1}{k_2}\\frac{-r_A}{C_A}', 'Eq. 3.62'],
      ['t = \\int_{C_{A}}^{C_{A0}}\\frac{dC_A}{-r_A}', 'Eq. 5.4', 'Batch time by integrating tabulated rates.']
    ],
    method: ['Slopes: the local option fits a quadratic through each point and its neighbours (a Savitzky–Golay estimate); the ln-polynomial option fits ln C_A as a polynomial in t, which follows decaying curves without the wiggles the book warns about in Fig. E3.2d; the global polynomial and finite differences are there to show how easily slopes go wrong.', 'The multi-component option solves log r = log k + a log C_A + b log C_B + c log C_C by least squares with 95% intervals on each order.'],
    run(v) {
      const kpi = [], plots = [], tables = [], warn = [];
      if (v.kind === 'multi') {
        const rows = (v.mdata || []).filter(r => r.every(x => isFinite(x) && x > 0));
        if (rows.length < 4) throw new Error('Need at least four complete runs.');
        const A = rows.map(r => [1, Math.log(r[0]), Math.log(r[1]), Math.log(r[2])]);
        const y = rows.map(r => Math.log(r[3]));
        const s = N.lstsq(A, y), p = s.x;
        const fit = A.map(a => a.reduce((t, x, j) => t + x * p[j], 0));
        const res = y.map((q, i) => q - fit[i]), dof = rows.length - 4, s2 = dof > 0 ? N.sum(res.map(r => r * r)) / dof : NaN;
        const tq = dof > 0 ? N.tinv(0.975, dof) : NaN;
        const se = j => Math.sqrt(s.covU[j][j] * s2);
        ['k', 'a (order in A)', 'b (order in B)', 'c (order in C)'].forEach((l, j) => kpi.push({ k: 'p' + j, l, v: j ? p[j] : Math.exp(p[0]), main: j === 0, s: isFinite(tq) ? (j ? '± ' + N.fmt(tq * se(j), 2) : '×÷ ' + N.fmt(Math.exp(tq * se(0)), 3)) : 'too few runs for intervals' }));
        kpi.push({ k: 'r2', l: 'R² (log)', v: N.r2(y, fit) });
        tables.push({ title: 'Runs and fit', cols: ['C_A', 'C_B', 'C_C', 'rate', 'fitted'], rows: rows.map((r, i) => r.concat([Math.exp(fit[i])])) });
        plots.push({ title: 'Parity: fitted against measured rate', x: { label: 'measured', log: true }, y: { label: 'fitted', log: true }, series: [{ name: 'runs', type: 'pts', x: rows.map(r => r[3]), y: fit.map(Math.exp), color: 2 }, { name: 'y = x', x: [Math.min(...rows.map(r => r[3])) * 0.8, Math.max(...rows.map(r => r[3])) * 1.2], y: [Math.min(...rows.map(r => r[3])) * 0.8, Math.max(...rows.map(r => r[3])) * 1.2], color: 8, dash: '4 4' }] });
        warn.push('Rounded orders suggest a mechanism: with a ≈ 1, b ≈ 1 and c ≈ 1 the rate is k[A][B][C] — consistent with C combining with one reactant in a fast pre-equilibrium.');
        return { kpi, plots, tables, warn };
      }
      let pairs, ct = null;
      if (v.kind === 'ct') {
        const rows = (v.data || []).filter(r => isFinite(r[0]) && isFinite(r[1]) && r[1] > 0).sort((a, b) => a[0] - b[0]);
        if (rows.length < 4) throw new Error('Need at least four points.');
        const t = rows.map(r => r[0]), C = rows.map(r => r[1]);
        let dC, smooth;
        if (v.how === 'logpoly') {
          const deg = Math.min(+v.deg, t.length - 1), f = N.polyfit(t, C.map(Math.log), deg), dc = N.polyder(f.c);
          smooth = x => Math.exp(N.polyval(f.c, x));
          dC = t.map(x => smooth(x) * N.polyval(dc, x));
        } else if (v.how === 'poly') {
          const deg = Math.min(+v.deg, t.length - 1), f = N.polyfit(t, C, deg);
          smooth = x => N.polyval(f.c, x);
          dC = N.diffData(t, C, 'poly', { deg });
        } else {
          dC = N.diffData(t, C, v.how === 'local' ? 'local' : 'central', { half: +v.half });
          smooth = N.pchip(t, C);
        }
        pairs = t.map((x, i) => ({ t: x, C: C[i], r: -dC[i] }));
        const tt = N.linspace(t[0], t[t.length - 1], 160);
        const tan = [];
        pairs.forEach(p => { const L = (t[t.length - 1] - t[0]) * 0.06; tan.push({ t: 'seg', x0: p.t - L, x1: p.t + L, y0: p.C + p.r * L, y1: p.C - p.r * L, color: 2, w: 1.6 }); });
        ct = { title: 'Smooth curve and tangents', x: { label: 't', unit: v._tu, min: 0 }, y: { label: 'C_A', unit: v._cu, min: 0 }, series: [{ name: 'data', type: 'pts', x: t, y: C, color: 8 }, { name: 'smooth curve', x: tt, y: tt.map(smooth), color: 1 }], marks: tan };
        tables.push({ title: 'Slopes', cols: ['t', 'C_A', 'dC_A/dt', 'log₁₀(−dC_A/dt)', 'log₁₀ C_A'], rows: pairs.map(p => [p.t, p.C, -p.r, p.r > 0 ? Math.log10(p.r) : NaN, Math.log10(p.C)]) });
      } else {
        pairs = (v.rdata || []).filter(r => isFinite(r[0]) && isFinite(r[1]) && r[0] > 0 && r[1] > 0).map(r => ({ C: r[0], r: r[1] })).sort((a, b) => a.C - b.C);
        if (pairs.length < 2) throw new Error('Need at least two rate points.');
      }
      const good = pairs.filter(p => p.r > 0 && p.C > 0);
      if (good.length < 2) throw new Error('No positive rates; check the data.');
      if (v.form === 'sat') {
        const l1 = N.linreg(good.map(p => 1 / p.C), good.map(p => 1 / p.r));
        const k1a = 1 / l1.b, k2a = l1.a * k1a;
        const l2 = N.linreg(good.map(p => p.r / p.C), good.map(p => p.r));
        const k2b = -1 / l2.b, k1b = l2.a * k2b;
        kpi.push({ k: 'k1', l: 'k₁ (Eq. 3.61)', v: k1a, u: '1/' + v._tu, main: true }, { k: 'k2', l: 'k₂ (Eq. 3.61)', v: k2a, u: '(' + v._cu + ')⁻¹' }, { k: 'k1b', l: 'k₁ (Eq. 3.62)', v: k1b, u: '1/' + v._tu }, { k: 'k2b', l: 'k₂ (Eq. 3.62)', v: k2b, u: '(' + v._cu + ')⁻¹' });
        const xi = N.linspace(0, Math.max(...good.map(p => 1 / p.C)) * 1.05, 20), xj = N.linspace(0, Math.max(...good.map(p => p.r / p.C)) * 1.05, 20);
        plots.push({ title: '1/(−r_A) against 1/C_A (Eq. 3.61)', x: { label: '1/C_A', min: 0 }, y: { label: '1/(−r_A)' }, series: [{ name: 'data', type: 'pts', x: good.map(p => 1 / p.C), y: good.map(p => 1 / p.r), color: 2 }, { name: 'line', x: xi, y: xi.map(x => l1.a + l1.b * x), color: 1 }] });
        plots.push({ title: '−r_A against (−r_A)/C_A (Eq. 3.62)', x: { label: '(−r_A)/C_A', min: 0 }, y: { label: '−r_A' }, series: [{ name: 'data', type: 'pts', x: good.map(p => p.r / p.C), y: good.map(p => p.r), color: 2 }, { name: 'line', x: xj, y: xj.map(x => l2.a + l2.b * x), color: 1 }] });
      } else {
        const lr = N.linreg(good.map(p => Math.log10(p.C)), good.map(p => Math.log10(p.r)));
        const tq = lr.dof > 0 ? N.tinv(0.975, lr.dof) : NaN;
        kpi.push({ k: 'n', l: 'Order n', v: lr.b, main: true, s: isFinite(tq) ? '± ' + N.fmt(tq * lr.sb, 2) + ' (95%)' : '' }, { k: 'k', l: 'k', v: Math.pow(10, lr.a), u: KIT.kUnit(lr.b, v._cu, v._tu) }, { k: 'r2', l: 'R² (log–log)', v: lr.r2 });
        const xs = good.map(p => p.C), xr = N.logspace(Math.min(...xs) * 0.8, Math.max(...xs) * 1.25, 30);
        plots.push({ title: 'log(−r_A) against log C_A (Fig. 3.18)', x: { label: 'C_A', unit: v._cu, log: true }, y: { label: '−r_A', unit: v._cu + '/' + v._tu, log: true }, series: [{ name: 'slopes', type: 'pts', x: xs, y: good.map(p => p.r), color: 2 }, { name: 'fit', x: xr, y: xr.map(x => Math.pow(10, lr.a) * Math.pow(x, lr.b)), color: 1 }] });
      }
      if (ct) plots.unshift(ct);
      if (v.kind === 'rate') {
        const xs = pairs.map(p => p.C), inv = N.pchip(xs, pairs.map(p => 1 / p.r));
        const a = Math.min(+v.C1, +v.C2), b = Math.max(+v.C1, +v.C2);
        if (a < xs[0] || b > xs[xs.length - 1]) warn.push('The integration range goes beyond the data; values outside are extrapolated.');
        const tB = N.quad(inv, a, b);
        kpi.push({ k: 'tB', l: 'Batch time ' + N.fmt(b, 3) + ' → ' + N.fmt(a, 3), v: tB, u: v._tu, main: true });
        const cx = N.linspace(xs[0], xs[xs.length - 1], 200), sx = N.linspace(a, b, 80);
        plots.push({ title: 'Area = batch time', x: { label: 'C_A', unit: v._cu, min: 0 }, y: { label: '1/(−r_A)', min: 0 }, series: [{ name: '1/(−r_A)', x: cx, y: cx.map(inv), color: 1 }, { name: 'data', type: 'pts', x: xs, y: pairs.map(p => 1 / p.r), color: 2 }, { x: sx, y: sx.map(inv), type: 'area', fill: 0.22, color: 3, hover: false, legend: false }] });
      }
      return { kpi, plots, tables, warn };
    }
  });

  /* ---------------- fractional life ---------------- */
  ST.tool({
    id: 'fraclife', tag: 'H-303', ch: 3, sec: '3.1', title: 'Half-life and fractional-life method',
    lede: 'Overall order from how the time for a fixed fractional drop depends on the starting concentration, from separate runs or from points picked along one smooth batch curve.',
    keys: 'half life fractional life order of reaction method of excess log t vs log C nitric oxide hydrogen',
    units: { C: 'mol/L', t: 's' }, sweep: false,
    inputs: [
      { group: 'Data' },
      { k: 'mode', t: 'seg', label: 'Source', opts: [['runs', 'separate runs'], ['curve', 'one batch curve']], v: 'curve' },
      { k: 'F', label: 'Fraction remaining F = C_A/C_{A0}', v: 0.8, min: 0.05, max: 0.95, step: 0.05 },
      { k: 'runs', t: 'table', label: 'Runs', cols: ['C_{A0}', 't_F'], v: [[200, 265], [240, 186], [280, 115], [320, 104], [360, 67]], if: v => v.mode === 'runs' },
      { k: 'data', t: 'table', label: 'Batch curve', cols: ['t', 'C_A'], v: [[0, 10], [20, 8], [40, 6], [60, 5], [120, 3], [180, 2], [300, 1]], if: v => v.mode === 'curve' },
      { k: 'starts', t: 'text', label: 'Starting concentrations to use', v: '10, 5, 2', mono: true, if: v => v.mode === 'curve' },
      { k: 'smooth', t: 'select', label: 'Smooth curve', opts: [['logpoly', 'polynomial in ln C_A'], ['pchip', 'monotone curve through the points']], v: 'logpoly', if: v => v.mode === 'curve' }
    ],
    presets: [
      { name: 'Example 3.1', v: { mode: 'curve', F: 0.8, data: [[0, 10], [20, 8], [40, 6], [60, 5], [120, 3], [180, 2], [300, 1]], starts: '10, 5, 2', smooth: 'logpoly', _cu: 'mol/L', _tu: 's' }, book: { n: { v: 1.4, tol: 0.08 }, k: { v: 0.005, tol: 0.2 } }, note: 'The book reads t_F from a smooth curve drawn by eye (18.5 s from 10 mol/L) and finds n = 1.4, k = 0.005.' },
      { name: 'Problem 3.26', v: { mode: 'curve', F: 0.5, data: [[0, 10], [20, 8], [40, 6], [60, 5], [120, 3], [180, 2], [300, 1]], starts: '10, 6, 2', smooth: 'logpoly' }, note: 'Same data with the half-life method.' },
      { name: 'Problem 3.8', v: { mode: 'runs', F: 0.5, runs: [[200, 265], [240, 186], [280, 115], [320, 104], [360, 67]], _cu: 'mmHg', _tu: 's' }, note: 'Equimolar H₂ + NO: half-life against initial total pressure. The pressure stands in for concentration; the order comes from the slope.' }
    ],
    eqs: [
      ['t_F = \\frac{F^{1-n}-1}{k(n-1)}C_{A0}^{1-n}', 'Eq. 3.33b'],
      ['t_{1/2} = \\frac{0.5^{1-n}-1}{k(n-1)}C_{A0}^{1-n}', 'Eq. 3.33a'],
      ['\\log t_F = \\text{const} + (1-n)\\log C_{A0}', 'Fig. 3.5']
    ],
    run(v) {
      const F = +v.F;
      let runs = [];
      const plots = [], warn = [];
      if (v.mode === 'runs') runs = (v.runs || []).filter(r => r[0] > 0 && r[1] > 0).map(r => ({ C0: r[0], tF: r[1] }));
      else {
        const rows = (v.data || []).filter(r => isFinite(r[0]) && r[1] > 0).sort((a, b) => a[0] - b[0]);
        if (rows.length < 3) throw new Error('Need at least three points on the curve.');
        const t = rows.map(r => r[0]), C = rows.map(r => r[1]);
        let Cof;
        if (v.smooth === 'logpoly') { const f = N.polyfit(t, C.map(Math.log), Math.min(3, t.length - 1)); Cof = x => Math.exp(N.polyval(f.c, x)); }
        else Cof = N.pchip(t, C);
        const tOf = c => N.root(x => Cof(x) - c, t[0], t[t.length - 1] * 3, { n: 300 });
        const starts = String(v.starts).split(/[,\s]+/).map(Number).filter(x => x > 0);
        const marks = [];
        starts.forEach(c0 => {
          const ta = tOf(c0), tb = tOf(F * c0);
          if (isFinite(ta) && isFinite(tb)) { runs.push({ C0: c0, tF: tb - ta, ta, tb }); marks.push({ t: 'seg', x0: ta, x1: tb, y0: c0, y1: c0, color: 2, w: 2 }, { t: 'seg', x0: tb, x1: tb, y0: c0, y1: F * c0, color: 2, w: 2 }); }
          else warn.push('Could not locate ' + N.fmt(c0, 3) + ' or ' + N.fmt(F * c0, 3) + ' on the curve.');
        });
        const tt = N.linspace(t[0], t[t.length - 1], 160);
        plots.push({ title: 'Times read from the smooth curve', x: { label: 't', unit: v._tu, min: 0 }, y: { label: 'C_A', unit: v._cu, min: 0 }, series: [{ name: 'data', type: 'pts', x: t, y: C, color: 8 }, { name: 'smooth curve', x: tt, y: tt.map(Cof), color: 1 }], marks });
      }
      if (runs.length < 2) throw new Error('Need at least two runs.');
      const lr = N.linreg(runs.map(r => Math.log10(r.C0)), runs.map(r => Math.log10(r.tF)));
      const n = 1 - lr.b;
      const kOf = r => (Math.abs(n - 1) < 1e-9 ? -Math.log(F) / r.tF : (Math.pow(F, 1 - n) - 1) * Math.pow(r.C0, 1 - n) / ((n - 1) * r.tF));
      const ks = runs.map(kOf), k = Math.exp(N.mean(ks.map(Math.log)));
      const kpi = [{ k: 'n', l: 'Order n', v: n, main: true, s: lr.dof > 0 ? '± ' + N.fmt(N.tinv(0.975, lr.dof) * lr.sb, 2) + ' (95%)' : 'from two runs' }, { k: 'k', l: 'k', v: k, u: KIT.kUnit(n, v._cu, v._tu) }, { k: 'r2', l: 'R²', v: lr.r2 }];
      const xr = N.logspace(Math.min(...runs.map(r => r.C0)) * 0.8, Math.max(...runs.map(r => r.C0)) * 1.25, 20);
      plots.push({ title: 'log t_F against log C_{A0}', x: { label: 'C_{A0}', unit: v._cu, log: true }, y: { label: 't_F', unit: v._tu, log: true }, series: [{ name: 'runs', type: 'pts', x: runs.map(r => r.C0), y: runs.map(r => r.tF), color: 2 }, { name: 'slope 1 − n', x: xr, y: xr.map(x => Math.pow(10, lr.a) * Math.pow(x, lr.b)), color: 1 }] });
      return { kpi, warn, plots, tables: [{ title: 'Runs', cols: ['C_{A0}', 't_F', 'k from this run'], rows: runs.map((r, i) => [r.C0, r.tF, ks[i]]) }] };
    }
  });

  /* ---------------- multiple reactions in a batch ---------------- */
  ST.tool({
    id: 'multibatch', tag: 'H-304', ch: 3, sec: '3.1', title: 'Parallel, series and catalysed reactions from batch data',
    lede: 'Rate constants for A → R and A → S in parallel, for A → R → S in series (with the maximum of the intermediate), and for a homogeneous catalysed reaction from runs at several catalyst levels.',
    keys: 'parallel reactions series reactions consecutive intermediate maximum k1 k2 catalysed homogeneous catalyst k observed',
    units: { C: 'mol/L', t: 'min' }, sweep: false,
    inputs: [
      { group: 'Scheme' },
      { k: 'mode', t: 'select', label: 'Case', opts: [['par', 'A → R, A → S (parallel, first order)'], ['ser', 'A → R → S (series, first order)'], ['serk', 'A → R → S: curves from k₁, k₂'], ['cat', 'catalysed: k_obs at several catalyst levels']], v: 'par' },
      { k: 'pdata', t: 'table', label: 'Run', cols: ['t', 'C_A', 'C_R', 'C_S'], v: [[0, 1, 0, 0], [1, 0.74, 0.17, 0.09], [2, 0.55, 0.30, 0.15], [4, 0.30, 0.46, 0.24], [6, 0.165, 0.555, 0.28], [10, 0.05, 0.63, 0.32]], if: v => v.mode === 'par' },
      { k: 'sdata', t: 'table', label: 'Run', cols: ['t', 'C_A', 'C_R'], v: [[0, 1, 0], [0.5, 0.74, 0.24], [1, 0.55, 0.38], [2, 0.30, 0.47], [3, 0.165, 0.45], [5, 0.05, 0.33], [8, 0.008, 0.17]], if: v => v.mode === 'ser' },
      { k: 'k1', label: 'k₁', unit: '1/[t]', v: 0.6, min: 1e-4, max: 1e3, scale: 'log', if: v => v.mode === 'serk' },
      { k: 'k2', label: 'k₂', unit: '1/[t]', v: 0.3, min: 1e-4, max: 1e3, scale: 'log', if: v => v.mode === 'serk' },
      { k: 'CA0', label: '$C_{A0}$', unit: '[C]', v: 1, min: 1e-6, max: 1e6, scale: 'log', if: v => v.mode === 'serk' },
      { k: 'cdata', t: 'table', label: 'Runs', cols: ['C_C', 'k_{obs}'], v: [[0, 0.11], [0.01, 0.16], [0.02, 0.215], [0.04, 0.31], [0.08, 0.52]], if: v => v.mode === 'cat' }
    ],
    presets: [
      { name: 'Problem 3.14', v: { mode: 'serk', k1: 1, k2: 1, CA0: 1 }, note: 'k₁ = k₂: the intermediate peaks at C_{A0}/e when t = 1/k.' },
      { name: 'Sample: parallel', v: { mode: 'par' }, note: 'Synthetic data with k₁ = 0.2, k₂ = 0.1 and a little noise.' },
      { name: 'Sample: series', v: { mode: 'ser' }, note: 'Synthetic data with k₁ = 0.6, k₂ = 0.3.' },
      { name: 'Sample: catalysed', v: { mode: 'cat' }, note: 'Synthetic runs with k₁ = 0.11 and k₂ = 5.1.' }
    ],
    eqs: [
      ['-\\ln\\frac{C_A}{C_{A0}} = (k_1+k_2)t', 'Eq. 3.37'],
      ['\\frac{C_R-C_{R0}}{C_S-C_{S0}} = \\frac{k_1}{k_2}', 'Eq. 3.38'],
      ['\\frac{C_R}{C_{A0}} = \\frac{k_1}{k_2-k_1}\\left(e^{-k_1t}-e^{-k_2t}\\right)', 'Eq. 3.49'],
      ['t_{max} = \\frac{\\ln(k_2/k_1)}{k_2-k_1},\\qquad \\frac{C_{R,max}}{C_{A0}} = \\left(\\frac{k_1}{k_2}\\right)^{k_2/(k_2-k_1)}', 'Eqs. 3.51–3.52'],
      ['k_{obs} = k_1 + k_2C_C', 'Eq. 3.40']
    ],
    method: ['Parallel: k₁ + k₂ from the slope of −ln(C_A/C_{A0}) against t, k₁/k₂ from the slope of C_R against C_S; a joint least-squares fit to all three curves refines both.', 'Series: a joint fit to C_A and C_R. With only C_R measured, (k₁, k₂) and (k₂, k₁) give curves of the same shape and different heights, which is why C_A must also be followed.'],
    run(v) {
      const kpi = [], plots = [], tables = [], warn = [];
      if (v.mode === 'par') {
        const rows = (v.pdata || []).filter(r => r.slice(0, 4).every(isFinite)).sort((a, b) => a[0] - b[0]);
        if (rows.length < 3) throw new Error('Need at least three times.');
        const CA0 = rows[0][1], CR0 = rows[0][2], CS0 = rows[0][3];
        const l1 = N.linreg(rows.map(r => r[0]), rows.map(r => Math.log(CA0 / r[1])), { origin: true });
        const l2 = N.linreg(rows.map(r => r[3] - CS0), rows.map(r => r[2] - CR0), { origin: true });
        const ksum = l1.b, ratio = l2.b;
        let k2 = ksum / (1 + ratio), k1 = ksum - k2;
        const model = (p, t) => { const e = Math.exp(-(p[0] + p[1]) * t), f = 1 - e; return [CA0 * e, CR0 + CA0 * p[0] / (p[0] + p[1]) * f, CS0 + CA0 * p[1] / (p[0] + p[1]) * f]; };
        const fit = N.lm(p => { const out = []; rows.forEach(r => { const m = model(p, r[0]); out.push(m[0] - r[1], m[1] - r[2], m[2] - r[3]); }); return out; }, [k1, k2], { pos: [true, true] });
        k1 = fit.p[0]; k2 = fit.p[1];
        kpi.push({ k: 'k1', l: 'k₁ (A → R)', v: k1, u: '1/' + v._tu, main: true }, { k: 'k2', l: 'k₂ (A → S)', v: k2, u: '1/' + v._tu }, { k: 'sum', l: 'k₁ + k₂ from slope', v: ksum, u: '1/' + v._tu }, { k: 'ratio', l: 'k₁/k₂ from C_R vs C_S', v: ratio });
        const ts = N.linspace(0, rows[rows.length - 1][0] * 1.05, 120);
        plots.push({ title: 'Concentrations', x: { label: 't', unit: v._tu, min: 0 }, y: { label: 'C', unit: v._cu, min: 0 }, series: [{ name: 'C_A', type: 'pts', x: rows.map(r => r[0]), y: rows.map(r => r[1]), color: 1 }, { name: 'C_R', type: 'pts', x: rows.map(r => r[0]), y: rows.map(r => r[2]), color: 2 }, { name: 'C_S', type: 'pts', x: rows.map(r => r[0]), y: rows.map(r => r[3]), color: 3 }, { name: 'fit A', x: ts, y: ts.map(t => model(fit.p, t)[0]), color: 1, legend: false }, { name: 'fit R', x: ts, y: ts.map(t => model(fit.p, t)[1]), color: 2, legend: false }, { name: 'fit S', x: ts, y: ts.map(t => model(fit.p, t)[2]), color: 3, legend: false }] });
        plots.push({ title: 'C_R against C_S (Fig. 3.6)', x: { label: 'C_S − C_{S0}', min: 0 }, y: { label: 'C_R − C_{R0}', min: 0 }, series: [{ name: 'data', type: 'pts', x: rows.map(r => r[3] - CS0), y: rows.map(r => r[2] - CR0), color: 2 }, { name: 'slope k₁/k₂', x: [0, Math.max(...rows.map(r => r[3] - CS0))], y: [0, ratio * Math.max(...rows.map(r => r[3] - CS0))], color: 1 }] });
      } else if (v.mode === 'ser' || v.mode === 'serk') {
        let k1, k2, CA0, rows = null;
        if (v.mode === 'ser') {
          rows = (v.sdata || []).filter(r => r.slice(0, 3).every(isFinite)).sort((a, b) => a[0] - b[0]);
          if (rows.length < 3) throw new Error('Need at least three times.');
          CA0 = rows[0][1];
          const l1 = N.linreg(rows.map(r => r[0]), rows.map(r => Math.log(CA0 / Math.max(r[1], 1e-12))), { origin: true });
          const CRm = (p, t) => (Math.abs(p[1] - p[0]) < 1e-9 ? CA0 * p[0] * t * Math.exp(-p[0] * t) : CA0 * p[0] / (p[1] - p[0]) * (Math.exp(-p[0] * t) - Math.exp(-p[1] * t)));
          let best = null;
          [0.3, 1, 3].forEach(f => {
            const r = N.lm(p => { const out = []; rows.forEach(rw => { out.push(CA0 * Math.exp(-p[0] * rw[0]) - rw[1], CRm(p, rw[0]) - rw[2]); }); return out; }, [l1.b, l1.b * f * 1.01], { pos: [true, true] });
            if (!best || r.ssr < best.ssr) best = r;
          });
          k1 = best.p[0]; k2 = best.p[1];
          kpi.push({ k: 'k1', l: 'k₁', v: k1, u: '1/' + v._tu, main: true }, { k: 'k2', l: 'k₂', v: k2, u: '1/' + v._tu });
        } else { k1 = +v.k1; k2 = +v.k2; CA0 = +v.CA0; }
        const same = Math.abs(k1 - k2) < 1e-9 * Math.max(k1, k2);
        const tmax = same ? 1 / k1 : Math.log(k2 / k1) / (k2 - k1);
        const CRmax = same ? CA0 / Math.E : CA0 * Math.pow(k1 / k2, k2 / (k2 - k1));
        kpi.push({ k: 'tmax', l: 't at C_{R,max}', v: tmax, u: v._tu, main: v.mode === 'serk' }, { k: 'CRmax', l: 'C_{R,max}', v: CRmax, u: v._cu });
        const CR = t => (same ? CA0 * k1 * t * Math.exp(-k1 * t) : CA0 * k1 / (k2 - k1) * (Math.exp(-k1 * t) - Math.exp(-k2 * t)));
        const tEnd = Math.max(tmax * 4, rows ? rows[rows.length - 1][0] * 1.05 : 0);
        const ts = N.linspace(0, tEnd, 200);
        const series = [{ name: 'C_A', x: ts, y: ts.map(t => CA0 * Math.exp(-k1 * t)), color: 1 }, { name: 'C_R', x: ts, y: ts.map(CR), color: 2 }, { name: 'C_S', x: ts, y: ts.map(t => CA0 - CA0 * Math.exp(-k1 * t) - CR(t)), color: 3 }];
        if (rows) series.push({ name: 'A data', type: 'pts', x: rows.map(r => r[0]), y: rows.map(r => r[1]), color: 1, legend: false }, { name: 'R data', type: 'pts', x: rows.map(r => r[0]), y: rows.map(r => r[2]), color: 2, legend: false });
        plots.push({ title: 'A → R → S in a batch (Fig. 3.11)', x: { label: 't', unit: v._tu, min: 0 }, y: { label: 'C', unit: v._cu, min: 0 }, series, marks: [{ t: 'pt', x: tmax, y: CRmax, color: 2, label: 'max' }] });
        if (v.mode === 'ser') warn.push('The swapped pair k₁ = ' + N.fmt(k2, 3) + ', k₂ = ' + N.fmt(k1, 3) + ' gives an R curve of the same shape but a different height; the A data decide between them.');
      } else {
        const rows = (v.cdata || []).filter(r => isFinite(r[0]) && isFinite(r[1]));
        const lr = N.linreg(rows.map(r => r[0]), rows.map(r => r[1]));
        kpi.push({ k: 'k1', l: 'k₁ (uncatalysed)', v: lr.a, u: '1/' + v._tu, main: true }, { k: 'k2', l: 'k₂ (catalysed)', v: lr.b, u: '1/(' + v._cu + '·' + v._tu + ')' }, { k: 'r2', l: 'R²', v: lr.r2 });
        const xr = [0, Math.max(...rows.map(r => r[0])) * 1.05];
        plots.push({ title: 'k_obs against catalyst concentration (Fig. 3.8)', x: { label: 'C_C', unit: v._cu, min: 0 }, y: { label: 'k_{obs}', unit: '1/' + v._tu, min: 0 }, series: [{ name: 'runs', type: 'pts', x: rows.map(r => r[0]), y: rows.map(r => r[1]), color: 2 }, { name: 'k₁ + k₂C_C', x: xr, y: xr.map(x => lr.a + lr.b * x), color: 1 }] });
      }
      return { kpi, plots, tables, warn };
    }
  });

  /* ---------------- varying-volume batch ---------------- */
  ST.tool({
    id: 'varvol', tag: 'H-305', ch: 3, sec: '3.2', title: 'Varying-volume batch reactor',
    lede: 'Rate constants from volume–time data in a constant-pressure batch reactor, with the expansion factor built from stoichiometry and inerts.',
    keys: 'varying volume constant pressure batch expansion factor epsilon capillary bead volume change first order zero order',
    units: { C: 'mol/L', t: 'min' }, sweep: false,
    inputs: [
      { group: 'Expansion factor' },
      { k: 'epsMode', t: 'seg', label: 'ε_A from', opts: [['stoich', 'stoichiometry'], ['given', 'given value']], v: 'stoich' },
      { k: 'a', label: 'a in aA → products', v: 2, min: 0.1, max: 10, if: v => v.epsMode === 'stoich' },
      { k: 'pm', label: 'Moles of products from a moles of A', v: 1, min: 0, max: 20, if: v => v.epsMode === 'stoich' },
      { k: 'yA0', label: 'Mole fraction of A in the feed', v: 0.8, min: 0.01, max: 1, step: 0.01, if: v => v.epsMode === 'stoich' },
      { k: 'eps', label: '$\\varepsilon_A$', v: -0.4, min: -0.95, max: 5, step: 0.01, if: v => v.epsMode === 'given' },
      { group: 'Data' },
      { k: 'data', t: 'table', label: 'Volume against time', cols: ['t', 'V/V_0'], v: [[0, 1], [3, 0.8]] },
      { k: 'CA0', label: '$C_{A0}$ (for zero and second order)', unit: '[C]', v: 1, min: 1e-9, max: 1e9, scale: 'log' }
    ],
    presets: [
      { name: 'Problem 3.29', v: { epsMode: 'stoich', a: 2, pm: 1, yA0: 0.8, data: [[0, 1], [3, 0.8]] }, note: '2A → R at constant pressure, 80% A; the volume falls 20% in 3 min.' },
      { name: 'Problem 3.30', v: { epsMode: 'stoich', a: 1, pm: 1.6, yA0: 1, data: [[0, 1], [4, 1.5]], CA0: 1.2 / (0.082057 * 298.15) }, note: 'A → 1.6R from pure A at 1.2 atm and 25 °C; the volume rises 50% in 4 min.' }
    ],
    eqs: [
      ['V = V_0(1+\\varepsilon_AX_A),\\qquad \\varepsilon_A = \\frac{V_{X_A=1}-V_{X_A=0}}{V_{X_A=0}}', 'Eqs. 3.63–3.64'],
      ['\\frac{C_{A0}}{k\\varepsilon_A}\\ln\\frac{V}{V_0} = t', 'Eq. 3.70', 'Zero order.'],
      ['-\\ln\\left(1-\\frac{\\Delta V}{\\varepsilon_AV_0}\\right) = kt', 'Eq. 3.72', 'First order.'],
      ['\\frac{(1+\\varepsilon_A)\\Delta V}{\\varepsilon_AV_0-\\Delta V}+\\varepsilon_A\\ln\\left(1-\\frac{\\Delta V}{\\varepsilon_AV_0}\\right) = kC_{A0}t', 'Eq. 3.73', 'Second order, 2A → products or equimolar A + B.']
    ],
    run(v) {
      const eps = v.epsMode === 'stoich' ? +v.yA0 * ((+v.pm - +v.a) / +v.a) : +v.eps;
      if (Math.abs(eps) < 1e-9) throw new Error('ε_A is zero: the volume does not change, use the constant-volume tools.');
      const rows = (v.data || []).filter(r => isFinite(r[0]) && isFinite(r[1]) && r[0] > 0);
      if (!rows.length) throw new Error('Enter at least one time after zero.');
      const X = rows.map(r => (r[1] - 1) / eps);
      if (X.some(x => x <= 0 || x >= 1)) throw new Error('Volume changes imply conversions outside 0–1 for this ε_A.');
      const t = rows.map(r => r[0]);
      const k0 = N.linreg(t, rows.map(r => +v.CA0 / eps * Math.log(r[1])), { origin: true }).b;
      const k1 = N.linreg(t, X.map(x => -Math.log(1 - x)), { origin: true }).b;
      const k2 = N.linreg(t, X.map(x => ((1 + eps) * x / (1 - x) + eps * Math.log(1 - x)) / +v.CA0), { origin: true }).b;
      const kpi = [{ k: 'eps', l: 'ε_A', v: eps }, { k: 'k1', l: 'k, first order', v: k1, u: '1/' + v._tu, main: true }, { k: 'k0', l: 'k, zero order', v: k0, u: v._cu + '/' + v._tu }, { k: 'k2', l: 'k, second order', v: k2, u: KIT.kUnit(2, v._cu, v._tu) }];
      const tab = { title: 'Points', cols: ['t', 'V/V_0', 'X_A', 'k₁ from point'], rows: rows.map((r, i) => [r[0], r[1], X[i], -Math.log(1 - X[i]) / r[0]]) };
      const ts = N.linspace(0, Math.max(...t) * 1.4, 60);
      return { kpi, tables: [tab], plots: [{ title: 'Volume predicted by each order', x: { label: 't', unit: v._tu, min: 0 }, y: { label: 'V/V_0' }, series: [{ name: 'data', type: 'pts', x: [0].concat(t), y: [1].concat(rows.map(r => r[1])), color: 8 }, { name: 'first order', x: ts, y: ts.map(tt => 1 + eps * (1 - Math.exp(-k1 * tt))), color: 1 }, { name: 'zero order', x: ts, y: ts.map(tt => Math.exp(k0 * eps * tt / +v.CA0)), color: 2, dash: '5 4' }] }] };
    }
  });

  /* ---------------- Ch 4: stoichiometric relations ---------------- */
  ST.tool({
    id: 'stoich', tag: 'H-401', ch: 4, sec: '4.1', title: 'Conversion and concentration relations',
    lede: 'For aA + bB → rR + sS with inerts: expansion factors, both conversions and every concentration from any one measured quantity, for liquids, isothermal gases, or gases whose temperature and pressure also change.',
    keys: 'stoichiometric table expansion factor epsilon conversion concentration relation key reactant limiting reactant special case popcorn',
    units: { C: 'mol/L', t: 'min' }, sweep: false,
    inputs: [
      { group: 'Reaction aA + bB → rR + sS' },
      { k: 'a', label: 'a', v: 1, min: 0.1, max: 10 }, { k: 'b', label: 'b', v: 3, min: 0, max: 10 },
      { k: 'r', label: 'r', v: 6, min: 0, max: 20 }, { k: 's', label: 's', v: 0, min: 0, max: 20 },
      { group: 'Feed' },
      { k: 'phase', t: 'select', label: 'System', opts: [['liq', 'liquid, constant density'], ['gas', 'gas, constant T and π'], ['gasTP', 'gas, T and π change']], v: 'gas' },
      { k: 'CA0', label: '$C_{A0}$', unit: '[C]', v: 100, min: 0, max: 1e9 },
      { k: 'CB0', label: '$C_{B0}$', unit: '[C]', v: 200, min: 0, max: 1e9 },
      { k: 'CR0', label: '$C_{R0}$', unit: '[C]', v: 0, min: 0, max: 1e9 },
      { k: 'CS0', label: '$C_{S0}$', unit: '[C]', v: 0, min: 0, max: 1e9 },
      { k: 'CI0', label: 'Inerts $C_{I0}$', unit: '[C]', v: 100, min: 0, max: 1e9, if: v => v.phase !== 'liq' },
      { k: 'epsOver', t: 'check', label: 'Give ε_A directly instead', v: false },
      { k: 'epsIn', label: '$\\varepsilon_A$', v: 30, min: -0.99, max: 100, if: v => v.epsOver },
      { k: 'T0', label: 'Feed temperature T₀', unit: 'K', v: 400, min: 1, max: 3000, if: v => v.phase === 'gasTP' },
      { k: 'T', label: 'Temperature where measured, T', unit: 'K', v: 300, min: 1, max: 3000, if: v => v.phase === 'gasTP' },
      { k: 'p0', label: 'Feed pressure π₀', v: 4, min: 1e-3, max: 1e4, if: v => v.phase === 'gasTP' },
      { k: 'p', label: 'Pressure where measured, π', v: 3, min: 1e-3, max: 1e4, if: v => v.phase === 'gasTP' },
      { group: 'Known at the point of interest' },
      { k: 'known', t: 'select', label: 'Measured', opts: [['XA', 'conversion of A'], ['CA', 'C_A'], ['XB', 'conversion of B'], ['CB', 'C_B'], ['CR', 'C_R'], ['vr', 'volumetric flow ratio v/v₀']], v: 'CA' },
      { k: 'kv', label: 'Its value', v: 40, min: 0, max: 1e9, slider: false }
    ],
    presets: [
      { name: 'Example 4.1', v: { a: 1, b: 3, r: 6, s: 0, phase: 'gas', CA0: 100, CB0: 200, CR0: 0, CS0: 0, CI0: 100, epsOver: false, known: 'CA', kv: 40 }, book: { XA: 0.5, XB: 0.75, CB: 40, eA: 0.5 } },
      { name: 'Problem 4.1', v: { a: 1, b: 1, r: 1, s: 1, phase: 'gas', CA0: 100, CB0: 200, CR0: 0, CS0: 0, CI0: 0, epsOver: false, known: 'XA', kv: 0.8 } },
      { name: 'Problem 4.2', v: { a: 1, b: 2, r: 1, s: 1, phase: 'liq', CA0: 100, CB0: 100, CR0: 0, CS0: 0, CI0: 0, epsOver: false, known: 'CA', kv: 20 } },
      { name: 'Problem 4.3', v: { a: 1, b: 1, r: 1, s: 0, phase: 'gas', CA0: 200, CB0: 100, CR0: 0, CS0: 0, CI0: 0, epsOver: false, known: 'CA', kv: 50 } },
      { name: 'Problem 4.4', v: { a: 1, b: 2, r: 1, s: 0, phase: 'gas', CA0: 100, CB0: 100, CR0: 0, CS0: 0, CI0: 0, epsOver: false, known: 'CB', kv: 20 } },
      { name: 'Problem 4.5', v: { a: 1, b: 1, r: 2, s: 0, phase: 'gasTP', CA0: 100, CB0: 200, CR0: 0, CS0: 0, CI0: 0, epsOver: false, T0: 400, T: 300, p0: 4, p: 3, known: 'CA', kv: 20 } },
      { name: 'Problem 4.6', v: { a: 1, b: 1, r: 5, s: 0, phase: 'gasTP', CA0: 100, CB0: 200, CR0: 0, CS0: 0, CI0: 0, epsOver: false, T0: 1000, T: 400, p0: 5, p: 4, known: 'CA', kv: 20 } },
      { name: 'Problem 4.7 popcorn', v: { a: 1, b: 0, r: 1, s: 0, phase: 'gas', CA0: 1, CB0: 0, CR0: 0, CS0: 0, CI0: 0, epsOver: true, epsIn: 30, known: 'vr', kv: 28 }, note: 'Popping takes corn from 1 volume to 31, so ε = 30; 1 L/min in, 28 L/min out.' }
    ],
    eqs: [
      ['\\varepsilon_A = \\frac{V_{X_A=1}-V_{X_A=0}}{V_{X_A=0}},\\qquad \\frac{a\\varepsilon_A}{C_{A0}} = \\frac{b\\varepsilon_B}{C_{B0}},\\qquad \\frac{C_{A0}X_A}{a} = \\frac{C_{B0}X_B}{b}', 'Ch. 4'],
      ['X_A = \\frac{C_{A0}-C_A}{C_{A0}+\\varepsilon_AC_A},\\qquad \\frac{C_A}{C_{A0}} = \\frac{1-X_A}{1+\\varepsilon_AX_A}', 'Special case 2'],
      ['\\frac{C_A}{C_{A0}} = \\frac{1-X_A}{1+\\varepsilon_AX_A}\\left(\\frac{T_0\\pi}{T\\pi_0}\\right)', 'Special case 3'],
      ['\\frac{C_B}{C_{A0}} = \\frac{C_{B0}/C_{A0}-(b/a)X_A}{1+\\varepsilon_AX_A}\\left(\\frac{T_0\\pi}{T\\pi_0}\\right)', 'Special case 3']
    ],
    run(v) {
      const a = +v.a, b = +v.b, r = +v.r, s = +v.s;
      const CA0 = +v.CA0, CB0 = +v.CB0, CR0 = +v.CR0, CS0 = +v.CS0, CI0 = v.phase === 'liq' ? 0 : +v.CI0;
      const tot = CA0 + CB0 + CR0 + CS0 + CI0;
      let eA = v.phase === 'liq' ? 0 : (CA0 / tot) * ((r + s - a - b) / a);
      if (v.epsOver) eA = +v.epsIn;
      const F = v.phase === 'gasTP' ? (+v.T0 * +v.p) / (+v.T * +v.p0) : 1;
      const conc = X => { const den = (1 + eA * X) / F; return { A: CA0 * (1 - X) / den, B: (CB0 - (b / a) * CA0 * X) / den, R: (CR0 + (r / a) * CA0 * X) / den, S: (CS0 + (s / a) * CA0 * X) / den, I: CI0 / den }; };
      const Xmax = b > 0 && CB0 > 0 ? Math.min(1, (CB0 * a) / (b * CA0)) : (b > 0 ? 0 : 1);
      const kv = +v.kv;
      let XA;
      const solve = f => N.root(x => f(x) - kv, 0, 1 - 1e-12, { n: 600 });
      switch (v.known) {
        case 'XA': XA = kv; break;
        case 'CA': XA = solve(x => conc(x).A); break;
        case 'XB': XA = kv * CB0 * a / (b * CA0); break;
        case 'CB': XA = solve(x => conc(x).B); break;
        case 'CR': XA = solve(x => conc(x).R); break;
        case 'vr': XA = (kv / F - 1) / eA; break;
      }
      if (!isFinite(XA)) throw new Error('No conversion of A between 0 and 1 gives that value.');
      const warn = [];
      if (XA > Xmax + 1e-9) warn.push('Not physically possible: X_A = ' + N.fmt(XA, 4) + ' would consume more B than the feed holds. B runs out at X_A = ' + N.fmt(Xmax, 4) + ', so the measured value must refer to something else (check the stoichiometry).');
      const c = conc(XA), XB = b > 0 && CB0 > 0 ? b * CA0 * XA / (a * CB0) : NaN;
      const eB = b > 0 && CB0 > 0 ? eA * CB0 * a / (b * CA0) : NaN;
      const kpi = [{ k: 'XA', l: 'X_A', v: XA, main: true }, { k: 'XB', l: 'X_B', v: XB }, { k: 'eA', l: 'ε_A', v: eA }, { k: 'eB', l: 'ε_B', v: eB }, { k: 'CA', l: 'C_A', v: c.A, u: v._cu }, { k: 'CB', l: 'C_B', v: c.B, u: v._cu }, { k: 'CR', l: 'C_R', v: c.R, u: v._cu }, { k: 'vr', l: 'v/v₀', v: (1 + eA * XA) / F }];
      const xs = N.linspace(0, Math.max(0.05, Math.min(Xmax, 1)), 80);
      const cs = xs.map(conc);
      const series = [{ name: 'C_A', x: xs, y: cs.map(q => q.A), color: 1 }];
      if (b > 0) series.push({ name: 'C_B', x: xs, y: cs.map(q => q.B), color: 2 });
      if (r > 0) series.push({ name: 'C_R', x: xs, y: cs.map(q => q.R), color: 3 });
      if (s > 0) series.push({ name: 'C_S', x: xs, y: cs.map(q => q.S), color: 4 });
      if (CI0 > 0) series.push({ name: 'inerts', x: xs, y: cs.map(q => q.I), color: 8, dash: '4 4' });
      const tab = { title: 'Stoichiometric table at the point of interest', cols: ['', 'feed', 'here', 'change'], rows: [['A', CA0, c.A, c.A - CA0], ['B', CB0, c.B, c.B - CB0], ['R', CR0, c.R, c.R - CR0], ['S', CS0, c.S, c.S - CS0], ['inerts', CI0, c.I, c.I - CI0]] };
      return { kpi, warn, plots: [{ title: 'Concentrations against conversion of A', x: { label: 'X_A', min: 0 }, y: { label: 'C', unit: v._cu }, series, marks: [{ t: 'v', x: XA, label: 'here', color: 8 }] }], tables: [tab] };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
