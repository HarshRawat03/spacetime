/* Chapter 21 — deactivating catalysts: rate and deactivation equations from batch-solids
   experiments (with and without pore diffusion), and reactor performance over a run under
   constant flow, falling flow or rising temperature, with deactivation fronts and run-length choice */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, PEL = ST.PEL;
  const CAT = ST.CAT = ST.CAT || {};
  const fmt = (x, s = 4) => N.fmt(x, s);
  const Rg = 8.314462618;

  /* activity under concentration-independent deactivation of order d: da/dt = −kd a^d */
  CAT.act = (kd, d, t) => {
    if (Math.abs(d - 1) < 1e-12) return Math.exp(-kd * t);
    const base = 1 + (d - 1) * kd * t;
    if (d < 1) return base > 0 ? Math.pow(base, 1 / (1 - d)) : 0;
    return Math.pow(base, 1 / (1 - d));
  };
  /* conversion of an nth-order reaction (ε_A = 0) for a group g = k′a τ′ C_A0^{n−1} */
  CAT.Xplug = (n, g) => {
    if (!(g > 0)) return 0;
    if (Math.abs(n - 1) < 1e-12) return 1 - Math.exp(-g);
    if (n < 1) { const v = 1 - (1 - n) * g; return v <= 0 ? 1 : 1 - Math.pow(v, 1 / (1 - n)); }
    return 1 - Math.pow(1 + (n - 1) * g, -1 / (n - 1));
  };
  CAT.Xmix = (n, g) => {
    if (!(g > 0)) return 0;
    if (Math.abs(n - 1) < 1e-12) return g / (1 + g);
    return N.brent(X => g * Math.pow(1 - X, n) - X, 0, 1, { tol: 1e-14 });
  };
  const TAUU = [['kg·s/m³', 'kg·s/m³', 'm³/(kg·s)'], ['kg·h/m³', 'kg·h/m³', 'm³/(kg·h)'], ['g·min/L', 'g·min/L', 'L/(g·min)'], ['g·h/L', 'g·h/L', 'L/(g·h)'], ['kg·s/L', 'kg·s/L', 'L/(kg·s)']];
  const tauUnit = v => (TAUU.find(q => q[0] === v.tauU) || TAUU[0])[1];
  const kUnit = v => (TAUU.find(q => q[0] === v.tauU) || TAUU[0])[2];

  /* ===================================================================== */
  /* C-2101  rate and deactivation equations from experiment                */
  /* ===================================================================== */
  ST.tool({
    id: 'deactdata', tag: 'C-2101', ch: 21, sec: 'Ch. 21.2', title: 'Deactivation kinetics from batch-solids experiments, with pore diffusion',
    lede: 'Find the rate of reaction and the rate of deactivation from experiments on a batch of catalyst: batch fluid, steady mixed or plug flow, or flow lowered to hold the gas composition. The book’s straight-line tests (Eqs. 21.22, 21.26, 21.27, 21.30) are drawn and fitted; any reaction order and deactivation order can be fitted directly; and with the pellet’s size, density and diffusivity the same data are reinterpreted for strong pore resistance, with the Thiele modulus checked for consistency, as in Example 21.1.',
    keys: 'catalyst deactivation activity rate equation batch solids batch fluid mixed flow plug flow basket reactor recycle constant composition order of deactivation kd poisoning fouling strong pore diffusion thiele modulus example 21.1 problem 21.1 21.2 21.3 21.4 21.5 21.6 21.12 21.13 21.14',
    units: { t: 'h' },
    inputs: [
      { group: 'Experiment' },
      { k: 'mode', t: 'select', label: 'Experimental device', opts: [['mc', 'batch of solids, steady mixed flow of fluid'], ['pc', 'batch of solids, steady plug flow of fluid'], ['bb', 'batch of solids, batch of fluid'], ['var', 'basket reactor, flow lowered to hold the composition']], v: 'mc' },
      { k: 'tauU', t: 'select', label: 'Units of τ′ = W C_A0/F_A0', opts: TAUU.map(q => [q[0], q[1]]), v: 'kg·s/m³', if: v => v.mode !== 'bb' },
      { k: 'tau', label: 'Weight-time τ′', unit: tauUnit, v: 4000, min: 1e-12, max: 1e15, if: v => v.mode === 'mc' || v.mode === 'pc' },
      { k: 'runT', t: 'table', label: 'Run: t and X_A', cols: ['t', 'X_A'], units: ['[t]', ''], v: [[0, 0.75], [2, 0.64], [4, 0.52], [6, 0.39]], if: v => v.mode === 'mc' || v.mode === 'pc' },
      { k: 'WV', label: 'Catalyst per volume of fluid W/V', unit: 'kg/m³ (or g/L)', v: 1, min: 1e-12, max: 1e12, if: v => v.mode === 'bb' },
      { k: 'bbT', t: 'table', label: 'Run: t and C_A', cols: ['t', 'C_A'], units: ['[t]', ''], v: [[0, 1], [0.25, 0.802], [0.5, 0.675], [1, 0.532], [2, 0.422]], if: v => v.mode === 'bb' },
      { k: 'Xe', label: 'Equilibrium conversion X_Ae (1 = irreversible)', v: 1, min: 0.01, max: 1, if: v => v.mode === 'bb' },
      { k: 'varT', t: 'table', label: 'Runs: run no., C_A0, X_A (held), t, τ′', cols: ['run', 'C_A0', 'X_A', 't', 'τ′'], v: [[1, 1, 0.5, 0, 1], [1, 1, 0.5, 1, 2.71828], [1, 1, 0.5, 2, 7.38906], [1, 1, 0.5, 3, 20.0855], [2, 2, 0.667, 0, 2], [2, 2, 0.667, 1, 5.43656], [2, 2, 0.667, 2, 14.7781], [2, 2, 0.667, 3, 40.1711]], if: v => v.mode === 'var' },
      { group: 'Rate forms to fit', if: v => v.mode !== 'var' },
      { k: 'nfit', t: 'seg', label: 'Reaction order', opts: [['1', 'first'], ['2', 'second'], ['free', 'fit it']], v: '1', if: v => v.mode !== 'var' },
      { k: 'dfit', t: 'seg', label: 'Deactivation order', opts: [['1', 'first'], ['2', 'second'], ['free', 'fit it']], v: '1', if: v => v.mode !== 'var' },
      { group: 'Pellets (to test for pore resistance)', if: v => v.mode === 'mc' || v.mode === 'pc' },
      { k: 'pel', t: 'check', label: 'Check and correct for pore diffusion (first order)', v: true, if: v => v.mode === 'mc' || v.mode === 'pc' },
      { k: 'dp', label: 'Pellet diameter (sphere)', unit: 'm', v: 2.4e-3, min: 1e-9, max: 1, if: v => (v.mode === 'mc' || v.mode === 'pc') && v.pel },
      { k: 'rhos', label: 'Pellet density ρ_s', unit: 'kg/m³ cat', v: 1500, min: 1, max: 3e4, if: v => (v.mode === 'mc' || v.mode === 'pc') && v.pel },
      { k: 'De', label: 'Effective diffusivity D_e', unit: 'm³/(m cat·s)', v: 5e-10, min: 1e-20, max: 1, scale: 'log', if: v => (v.mode === 'mc' || v.mode === 'pc') && v.pel }
    ],
    presets: [
      { name: 'Example 21.1 — mixed flow with pore resistance', v: { mode: 'mc', tauU: 'kg·s/m³', tau: 4000, runT: [[0, 0.75], [2, 0.64], [4, 0.52], [6, 0.39]], nfit: '1', dfit: '1', pel: true, dp: 2.4e-3, rhos: 1500, De: 5e-10, _tu: 'h' }, book: { kFree: { v: 7.5e-4, tol: 0.02 }, kdFree: { v: 0.26, tol: 0.02 }, kStrong: { v: 0.27, tol: 0.03 }, kdStrong: { v: 0.52, tol: 0.02 }, MTstrong: { v: 360, tol: 0.03 } }, note: '2.4 mm pellets, D_e = 5 × 10⁻¹⁰ m³/m cat·s, ρ_s = 1500 kg/m³, τ′ = 4000 kg·s/m³. Assuming no pore resistance gives k′ = 7.5 × 10⁻⁴ and k_d = 0.26 hr⁻¹ — but then M_T = 18.9, a contradiction. Assuming strong resistance gives k′ = 0.27 m³/kg·s, k_d = 0.52 hr⁻¹ and M_T = 360, which is consistent.' },
      { name: 'Problem 21.1 — composition held, two runs', v: { mode: 'var', varT: [[1, 1, 0.5, 0, 1], [1, 1, 0.5, 1, Math.E], [1, 1, 0.5, 2, Math.E ** 2], [1, 1, 0.5, 3, Math.E ** 3], [2, 2, 0.667, 0, 2], [2, 2, 0.667, 1, 2 * Math.E], [2, 2, 0.667, 2, 2 * Math.E ** 2], [2, 2, 0.667, 3, 2 * Math.E ** 3]], tauU: 'g·min/L', _tu: 'h' }, note: 'C_A0 = 1 and 2 mol/L, conversions 50% and 66.7% held by lowering the flow; τ′ grows as eᵗ in both runs.' },
      { name: 'Problem 21.2 — composition held, two runs', v: { mode: 'var', varT: [[1, 2, 0.5, 0, 1], [1, 2, 0.5, 1, Math.E], [1, 2, 0.5, 2, Math.E ** 2], [1, 2, 0.5, 3, Math.E ** 3], [2, 20, 0.8, 0, 1], [2, 20, 0.8, 0.5, Math.E], [2, 20, 0.8, 1, Math.E ** 2], [2, 20, 0.8, 1.5, Math.E ** 3]], tauU: 'g·min/L', _tu: 'h' }, note: 'C_A0 = 2 and 20 mol/L at 50% and 80%: the second run deactivates twice as fast.' },
      { name: 'Problem 21.3 — automobile converter aging', v: { mode: 'pc', tauU: 'kg·s/m³', tau: 1, runT: [[5, 0.57], [10, 0.53], [15, 0.52], [20, 0.5], [25, 0.48], [30, 0.45], [35, 0.43], [40, 0.41]], nfit: '1', dfit: 'free', pel: false, _tu: 'h' }, note: 'τ′ is not reported, so it is set to 1 and k′ stands for k′τ′.' },
      { name: 'Problem 21.4 — recycle reactor, τ′ = 2 kg·s/L', v: { mode: 'mc', tauU: 'kg·s/L', tau: 2, runT: [[1, 0.889], [2, 0.865], [4, 0.804]], nfit: '1', dfit: '1', pel: false, _tu: 'h' } },
      { name: 'Problem 21.5 — batch, reversible (X_Ae = 0.5)', v: { mode: 'bb', WV: 1, Xe: 0.5, bbT: [[0, 1], [0.25, 0.901], [0.5, 0.83], [1, 0.766], [2, 0.711]], nfit: '1', dfit: 'free', _tu: 'h' }, note: 'C_A falls to 0.684 at long times, above the equilibrium value 0.5: the catalyst dies first.' },
      { name: 'Problem 21.6 — batch, irreversible', v: { mode: 'bb', WV: 1, Xe: 1, bbT: [[0, 1], [0.25, 0.802], [0.5, 0.675], [1, 0.532], [2, 0.422]], nfit: '1', dfit: 'free', _tu: 'h' }, note: 'C_A levels off at 0.363 because the catalyst dies.' },
      { name: 'Problem 21.12 — butane dehydrogenation, plug flow', v: { mode: 'pc', tauU: 'kg·h/m³', tau: 1.1, runT: [[0, 0.89], [1, 0.78], [2, 0.63], [3, 0.47], [4, 0.34], [5, 0.26]], nfit: '1', dfit: '1', pel: false, _tu: 'h' } },
      { name: 'Problem 21.13 — catalase, run E (small pellets)', v: { mode: 'pc', tauU: 'kg·s/m³', tau: 4100, runT: [[0, 0.795], [1.25, 0.635], [2, 0.51], [3, 0.397], [4.25, 0.255], [5, 0.22], [6, 0.15], [7, 0.104]], nfit: '1', dfit: '1', pel: true, dp: 72e-6, rhos: 630, De: 5e-10, _tu: 'h' } },
      { name: 'Problem 21.14 — catalase, run B (large pellets)', v: { mode: 'pc', tauU: 'kg·s/m³', tau: 4560, runT: [[0.25, 0.57], [1, 0.475], [2, 0.39], [3, 0.3], [4, 0.23], [5, 0.186], [6, 0.14], [7, 0.115]], nfit: '1', dfit: '1', pel: true, dp: 1.45e-3, rhos: 630, De: 5e-10, _tu: 'h' } }
    ],
    eqs: [
      ['-r_A^\\prime = k^\\prime C_A^n\\mathbf{a},\\qquad -\\frac{d\\mathbf{a}}{dt} = k_d C_A^m\\mathbf{a}^d', 'Eqs. 21.7–21.8'],
      ['\\ln\\ln\\frac{C_A}{C_{A\\infty}} = \\ln\\frac{Wk^\\prime}{Vk_d} - k_dt', 'Eq. 21.22', 'Batch solids, batch fluid.'],
      ['\\ln\\left(\\frac{C_{A0}}{C_A}-1\\right) = \\ln(k^\\prime\\tau^\\prime) - k_dt', 'Eq. 21.26', 'Batch solids, steady mixed flow.'],
      ['\\ln\\tau^\\prime = k_dt + \\ln\\frac{C_{A0}-C_A}{k^\\prime C_A}', 'Eq. 21.27', 'Mixed flow lowered to hold C_A.'],
      ['\\ln\\ln\\frac{C_{A0}}{C_A} = \\ln(k^\\prime\\tau^\\prime) - k_dt', 'Eq. 21.30', 'Batch solids, steady plug flow.'],
      ['\\mathscr{E} = \\frac{1}{M_{Td}},\\qquad M_{Td} = L\\sqrt{\\frac{k^{\\prime\\prime\\prime}\\mathbf{a}}{\\mathscr{D}_e}} = M_T\\mathbf{a}^{1/2}', 'Eq. 21.38', 'Strong pore resistance: observed decay rate halves.'],
      ['\\ln\\left(\\frac{C_{A0}}{C_A}-1\\right) = \\ln\\frac{k^\\prime\\tau^\\prime}{M_T} - \\frac{k_d}{2}t', 'Ex. 21.1 (ii)']
    ],
    method: ['The straight-line tests of the book are fitted by least squares. For any reaction order n and deactivation order d the model is also fitted directly to the measured conversions or concentrations: the activity follows its integrated decay law, and the conversion at each time comes from the exact steady mixed or plug-flow equation (or, for a batch of fluid, from integrating dC/dt = −(W/V)k′Cⁿa). Parameters are found by Levenberg–Marquardt with standard errors.', 'Pore diffusion (first order, spheres, L = d_p/6): if the diffusion-free interpretation gives M_T > 0.4 it is inconsistent; the strong-resistance interpretation takes the intercept as k′τ′/M_T and the slope as k_d/2, and solves k′ = (k′τ′/M_T)² L²ρ_s/(τ′²D_e); its M_T is then checked against 4.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'h', kU = kUnit(v);
      if (v.mode === 'var') {
        const rows = (v.varT || []).filter(r => r && r.slice(0, 5).every(z => z !== '' && isFinite(+z)));
        const runs = {};
        rows.forEach(r => { (runs[+r[0]] = runs[+r[0]] || []).push({ CA0: +r[1], X: +r[2], t: +r[3], tau: +r[4] }); });
        const res = Object.keys(runs).map(id => {
          const pts = runs[id].sort((a, b) => a.t - b.t);
          if (pts.length < 2) return null;
          const L = N.linreg(pts.map(p => p.t), pts.map(p => Math.log(p.tau)));
          const CA0 = pts[0].CA0, C = CA0 * (1 - pts[0].X), tau0 = Math.exp(L.a);
          return { id, CA0, C, kd: L.b, r0: (CA0 - C) / tau0, tau0, pts, r2: L.r2 };
        }).filter(Boolean);
        if (!res.length) throw new Error('Each run needs at least two readings.');
        tables.push({ title: 'Each run: ln τ′ against t (Eq. 21.27)', cols: ['run', 'C_A0', 'C_A held', 'slope = k_d at this C_A', 'τ′ at t = 0', 'fresh rate (C_A0 − C_A)/τ′₀', 'R²'], rows: res.map(q => [q.id, q.CA0, q.C, q.kd, q.tau0, q.r0, q.r2]) });
        if (res.length >= 2) {
          const Ln = N.linreg(res.map(q => Math.log(q.C)), res.map(q => Math.log(q.r0)));
          const Lm = N.linreg(res.map(q => Math.log(q.C)), res.map(q => Math.log(Math.max(q.kd, 1e-300))));
          const n = Ln.b, k = Math.exp(Ln.a), m = Lm.b, kd0 = Math.exp(Lm.a);
          kpi.push({ k: 'n', l: 'Order of reaction n', v: n, main: true, s: '−r′ = ' + fmt(k, 4) + ' C_A^{' + fmt(n, 3) + '} a' });
          kpi.push({ k: 'k', l: 'k′', v: k, u: kU + (Math.abs(n - 1) > 1e-6 ? ' (with C in the units entered)' : '') });
          kpi.push({ k: 'm', l: 'Concentration dependence of deactivation m', v: m, s: Math.abs(m) < 0.05 ? 'deactivation independent of concentration' : '−da/dt = ' + fmt(kd0, 4) + ' C_A^{' + fmt(m, 3) + '} a' });
          kpi.push({ k: 'kd', l: 'k_d', v: kd0, u: '1/' + tu });
          if (res.length > 2 && Ln.r2 < 0.98) warn.push('The runs do not fall on one power law for the rate: try another rate form.');
        } else info.push('One run gives k_d at that composition and one point of the fresh rate; add a run at another composition to find the orders.');
        plots.push({ title: 'ln τ′ against time (Fig. 21.4)', x: { label: 't', unit: tu, min: 0 }, y: { label: 'ln τ′' }, series: res.flatMap((q, i) => [{ name: 'run ' + q.id + ' (C_A = ' + fmt(q.C, 3) + ')', x: q.pts.map(p => p.t), y: q.pts.map(p => Math.log(p.tau)), type: 'pts', color: i + 1 }, { x: [0, Math.max(...q.pts.map(p => p.t))], y: [Math.log(q.tau0), Math.log(q.tau0) + q.kd * Math.max(...q.pts.map(p => p.t))], color: i + 1, legend: false }]), cap: 'Parallel lines: deactivation does not depend on concentration. Different slopes: it does.' });
        return { kpi, warn, info, plots, tables };
      }
      const nOpt = v.nfit, dOpt = v.dfit;
      if (v.mode === 'mc' || v.mode === 'pc') {
        const pts = (v.runT || []).filter(r => r && isFinite(+r[0]) && isFinite(+r[1]) && r[0] !== '' && r[1] !== '').map(r => ({ t: +r[0], X: +r[1] })).sort((a, b) => a.t - b.t);
        if (pts.length < 2) throw new Error('Enter at least two readings.');
        const tau = +v.tau;
        const yT = v.mode === 'mc' ? pts.map(p => Math.log(p.X / (1 - p.X))) : pts.map(p => Math.log(Math.log(1 / (1 - p.X))));
        const L = N.linreg(pts.map(p => p.t), yT);
        const kFree = Math.exp(L.a) / tau, kdFree = -L.b;
        kpi.push({ k: 'kFree', l: 'k′ (no pore resistance, first order / first order)', v: kFree, u: kU, main: true, s: 'intercept ln(k′τ′) = ' + fmt(L.a, 4) });
        kpi.push({ k: 'kdFree', l: 'k_d (no pore resistance)', v: kdFree, u: '1/' + tu, s: 'R² = ' + fmt(L.r2, 4) });
        /* general fit */
        const ts = pts.map(p => p.t), Xs = pts.map(p => p.X);
        const fixedN = nOpt !== 'free', fixedD = dOpt !== 'free';
        const pred = p => ts.map(t => { const a = CAT.act(p[1], p[3], t), g = p[0] * a * tau; return v.mode === 'mc' ? CAT.Xmix(p[2], g) : CAT.Xplug(p[2], g); });
        const p0 = [kFree, Math.max(kdFree, 1e-6), fixedN ? +nOpt : 1, fixedD ? +dOpt : 1];
        const fit = N.lm(p => pred(p).map((x, i) => x - Xs[i]), p0, { pos: [true, true, false, false], fixed: [false, false, fixedN, fixedD], maxIt: 200 });
        if (!fixedN || !fixedD || nOpt !== '1' || dOpt !== '1') {
          kpi.push({ k: 'kFit', l: 'Direct fit: k′', v: fit.p[0], s: 'n = ' + fmt(fit.p[2], 3) + ', d = ' + fmt(fit.p[3], 3) + ', k_d = ' + fmt(fit.p[1], 4) + ' /' + tu + (fixedN || fixedD ? '' : '') + ' (k′ in (conc)^{1−n}·' + kU + ')' });
        }
        const tg = N.linspace(Math.min(0, ts[0]), ts[ts.length - 1] * 1.1, 120);
        plots.push({ title: v.mode === 'mc' ? 'Test of Eq. 21.26: ln(C_A0/C_A − 1) against t (Fig. 21.3)' : 'Test of Eq. 21.30: ln ln(C_A0/C_A) against t (Fig. 21.5)', x: { label: 't', unit: tu }, y: { label: v.mode === 'mc' ? 'ln(C_A0/C_A − 1)' : 'ln ln(C_A0/C_A)' }, series: [{ name: 'runs', x: ts, y: yT, type: 'pts', color: 1 }, { name: 'slope −k_d = ' + fmt(L.b, 3), x: tg, y: tg.map(t => L.a + L.b * t), color: 2 }] });
        plots.push({ title: 'Conversion against time', x: { label: 't', unit: tu }, y: { label: 'X_A', min: 0, max: 1 }, series: [{ name: 'runs', x: ts, y: Xs, type: 'pts', color: 1 }, { name: 'fit: n = ' + fmt(fit.p[2], 3) + ', d = ' + fmt(fit.p[3], 3), x: tg, y: (() => { const keep = ts.slice(); return tg.map(t => { const a = CAT.act(fit.p[1], fit.p[3], t), g = fit.p[0] * a * tau; return v.mode === 'mc' ? CAT.Xmix(fit.p[2], g) : CAT.Xplug(fit.p[2], g); }); })(), color: 2 }] });
        if (v.pel) {
          const Lp = +v.dp / 6, rs = +v.rhos, De = +v.De;
          const tauS = v.tauU === 'kg·s/m³' ? tau : NaN;
          if (!isFinite(tauS)) warn.push('The pore-diffusion check needs τ′ in kg·s/m³ (with D_e in m³/m cat·s).');
          else {
            const MTfree = Lp * Math.sqrt(kFree * rs / De);
            kpi.push({ k: 'MTfree', l: 'M_T if there were no pore resistance', v: MTfree, s: MTfree < 0.4 ? 'consistent: no pore resistance' : 'inconsistent: pores matter', state: MTfree < 0.4 ? 'ok' : 'bad' });
            const kStrong = Math.exp(2 * L.a) * Lp * Lp * rs / (De * tauS * tauS), MTs = Lp * Math.sqrt(kStrong * rs / De);
            kpi.push({ k: 'kStrong', l: 'k′ assuming strong pore resistance', v: kStrong, u: kU, s: 'intercept = ln(k′τ′/M_T)' });
            kpi.push({ k: 'kdStrong', l: 'k_d assuming strong pore resistance', v: 2 * kdFree, u: '1/' + tu, s: 'slope = −k_d/2' });
            kpi.push({ k: 'MTstrong', l: 'M_T on that assumption', v: MTs, s: MTs > 4 ? 'consistent: strong resistance' : 'inconsistent', state: MTs > 4 ? 'ok' : 'bad' });
            if (MTfree >= 0.4 && MTs > 4) info.push('Only the strong-resistance interpretation is self-consistent (Example 21.1): the true catalyst is ' + fmt(kStrong / kFree, 3) + ' times as active as it appears, and it deactivates twice as fast as the observed decay suggests.');
            else if (MTfree < 0.4) info.push('The pellets are free of pore resistance; the first interpretation stands.');
            else info.push('Neither limit is consistent: the pellets are in the intermediate regime. Use C-1801 with M_T between 0.4 and 4.');
          }
        }
        return { kpi, warn, info, plots, tables };
      }
      /* batch of solids and fluid */
      const pts = (v.bbT || []).filter(r => r && isFinite(+r[0]) && isFinite(+r[1]) && r[0] !== '' && r[1] !== '').map(r => ({ t: +r[0], C: +r[1] })).sort((a, b) => a.t - b.t);
      if (pts.length < 3) throw new Error('Enter at least three readings.');
      const WV = +v.WV, C0 = pts[0].C, Xe = +v.Xe, Ce = C0 * (1 - Xe);
      const fixedN = nOpt !== 'free', fixedD = dOpt !== 'free';
      const tmax = pts[pts.length - 1].t;
      const model = (p, tg) => {
        const sol = N.ode((t, y) => { const a = CAT.act(p[1], p[3], t), c = Math.max(y[0] - Ce, 0); return [-WV * p[0] * Math.pow(c, p[2]) * a]; }, 0, [C0], Math.max(...tg), { rtol: 1e-9, atol: 1e-14 });
        return tg.map(t => sol.at(t)[0]);
      };
      const ts = pts.map(p => p.t), Cs = pts.map(p => p.C);
      const k0 = Math.log((C0 - Ce) / Math.max(Cs[1] - Ce, 1e-9 * C0)) / Math.max(ts[1], 1e-9) / WV;
      const fit = N.lm(p => model(p, ts).map((c, i) => (c - Cs[i]) / C0), [k0, 1, fixedN ? +nOpt : 1, fixedD ? +dOpt : 1], { pos: [true, true, false, false], fixed: [false, false, fixedN, fixedD], maxIt: 200 });
      const p = fit.p;
      const Cinf = (() => { const sol = N.ode((t, y) => { const a = CAT.act(p[1], p[3], t), c = Math.max(y[0] - Ce, 0); return [-WV * p[0] * Math.pow(c, p[2]) * a]; }, 0, [C0], tmax * 200 + 50 / Math.max(p[1], 1e-9), { rtol: 1e-9, atol: 1e-14 }); return sol.last()[0]; })();
      kpi.push({ k: 'k', l: 'k′ (fit)', v: p[0], main: true, s: 'n = ' + fmt(p[2], 3) + (fixedN ? ' (fixed)' : '') });
      kpi.push({ k: 'kd', l: 'k_d (fit)', v: p[1], u: '1/' + tu, s: 'd = ' + fmt(p[3], 3) + (fixedD ? ' (fixed)' : '') });
      kpi.push({ k: 'Cinf', l: 'C_A at very long times', v: Cinf, s: Xe < 1 ? 'equilibrium would be ' + fmt(Ce, 4) : 'the catalyst dies before A is used up' });
      if (Math.abs(p[2] - 1) < 1e-9 && Math.abs(p[3] - 1) < 1e-9 && Xe === 1) {
        const yy = Cs.map(c => Math.log(Math.log(c / Cinf)));
        plots.push({ title: 'Test of Eq. 21.22 (Fig. 21.2)', x: { label: 't', unit: tu }, y: { label: 'ln ln(C_A/C_A∞)' }, series: [{ name: 'readings', x: ts, y: yy.map(q => (isFinite(q) ? q : NaN)), type: 'pts', color: 1 }, { name: 'slope −k_d', x: [0, tmax], y: [Math.log(WV * p[0] / p[1]), Math.log(WV * p[0] / p[1]) - p[1] * tmax], color: 2 }] });
      }
      const tg = N.linspace(0, tmax * 1.5, 150);
      plots.push({ title: 'Concentration in the batch', x: { label: 't', unit: tu, min: 0 }, y: { label: 'C_A', min: 0 }, series: [{ name: 'readings', x: ts, y: Cs, type: 'pts', color: 1 }, { name: 'fit', x: tg, y: model(p, tg), color: 2 }, { name: 'C_A∞', x: [0, tmax * 1.5], y: [Cinf, Cinf], color: 'var(--muted)', dash: '4 4' }] });
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* C-2102  reactor performance with a deactivating catalyst               */
  /* ===================================================================== */
  ST.tool({
    id: 'deactdesign', tag: 'C-2102', ch: 21, sec: 'Ch. 21.3', title: 'Reactors with deactivating catalyst: runs, policies and regeneration',
    lede: 'What a batch of deactivating catalyst delivers over a run: conversion and activity against time, the mean conversion of Eq. 21.40, for plug or mixed flow, any reaction and deactivation order, with or without pore resistance. Deactivation that depends on reactant or product moves through the bed as a front, followed cell by cell. Three operating policies are compared — constant flow, flow cut to hold conversion, temperature raised to hold conversion up to the limit — and the best run length before regeneration is found.',
    keys: 'deactivation design mean conversion run length regeneration policy constant flow decreasing flow rising temperature isothermal optimal temperature policy activity front poisoning parallel series side by side deactivation strong pore diffusion example 21.2 problem 21.7 21.8 21.9 21.10 21.11 21.15 cycle time',
    units: { t: 'day' },
    inputs: [
      { group: 'Reaction (fresh catalyst, constant flow)' },
      { k: 'pat', t: 'seg', label: 'Fluid flow', opts: [['plug', 'plug'], ['mixed', 'mixed']], v: 'plug' },
      { k: 'n', label: 'Reaction order n', v: 2, min: 0, max: 4, step: 0.5 },
      { k: 'kt', label: 'k′ τ′ C_A0ⁿ⁻¹ (dimensionless group, fresh)', v: 100, min: 1e-9, max: 1e9, scale: 'log', hint: 'Example 21.2: 0.2 × 10 × 50 = 100.' },
      { k: 'CA0', label: 'C_A0 (for concentration-dependent deactivation)', unit: 'mol/m³', v: 50, min: 1e-12, max: 1e7 },
      { group: 'Deactivation' },
      { k: 'dtype', t: 'select', label: 'Deactivation depends on', opts: [['indep', 'nothing (independent)'], ['sum', 'C_A + C_R (constant for A ⇌ R isomerization)'], ['par', 'C_A (parallel: A → P↓)'], ['ser', 'C_R (series: R → P↓)']], v: 'indep' },
      { k: 'd', label: 'Deactivation order d', v: 0, min: 0, max: 4, step: 0.5 },
      { k: 'kd', label: 'k_d', unit: v => (v.dtype === 'indep' ? '1/[t]' : '(m³/mol)^m/[t]'), v: 8.3125e-3, min: 0, max: 1e9 },
      { k: 'm', label: 'Concentration order m', v: 1, min: 0, max: 3, step: 0.5, if: v => v.dtype !== 'indep' },
      { k: 'trun', label: 'Run length', unit: '[t]', v: 120, min: 1e-6, max: 1e7 },
      { group: 'Pore diffusion' },
      { k: 'pore', t: 'seg', label: 'Pellets', opts: [['none', 'no resistance'], ['strong', 'strong resistance (ℰ = 1/M_T a^½)'], ['exact', 'first order, exact ℰ']], v: 'none' },
      { k: 'MT0', label: 'M_T of the fresh catalyst', v: 10, min: 1e-4, max: 1e4, scale: 'log', if: v => v.pore === 'exact' },
      { k: 'shape', t: 'seg', label: 'Shape', opts: [['slab', 'flat plate'], ['cyl', 'cylinder'], ['sph', 'sphere']], v: 'sph', if: v => v.pore === 'exact' },
      { group: 'Operating policy' },
      { k: 'pol', t: 'select', label: 'Policy (Fig. 21.8)', opts: [['K', 'K: constant flow, isothermal'], ['L', 'L: flow lowered to hold X_A'], ['M', 'M: temperature raised to hold X_A']], v: 'K' },
      { k: 'phase', t: 'seg', label: 'Phase (for policy M)', opts: [['liq', 'liquid'], ['gas', 'ideal gas']], v: 'gas', if: v => v.pol === 'M' },
      { k: 'T0', label: 'Starting temperature', unit: 'K', v: 600, min: 100, max: 3000, if: v => v.pol === 'M' },
      { k: 'Tmax', label: 'Highest allowable temperature T*', unit: 'K', v: 800, min: 100, max: 3000, if: v => v.pol === 'M' },
      { k: 'ER', label: 'E/R of the reaction', unit: 'K', v: 7200, min: 0, max: 1e5, if: v => v.pol === 'M' },
      { k: 'EdR', label: 'E_d/R of deactivation', unit: 'K', v: 0, min: 0, max: 1e5, if: v => v.pol === 'M' },
      { k: 'EDR', label: 'E_diff/R (strong resistance only)', unit: 'K', v: 0, min: 0, max: 1e5, if: v => v.pol === 'M' && v.pore !== 'none' },
      { group: 'Regeneration' },
      { k: 'treg', label: 'Time to regenerate or replace', unit: '[t]', v: 0, min: 0, max: 1e7 }
    ],
    presets: [
      { name: 'Example 21.2a — zero-order deactivation', v: { pat: 'plug', n: 2, kt: 100, CA0: 50, dtype: 'indep', d: 0, kd: 8.3125e-3, trun: 120, pore: 'none', pol: 'K', treg: 0, _tu: 'day' }, book: { Xbar: { v: 0.96, tol: 0.01 }, X0: { v: 0.99, tol: 0.003 }, Xend: { v: 0.2, tol: 0.02 } }, note: 'Second-order isomerization in a packed bed: X_A = 100a/(1 + 100a). Over 120 days the mean conversions are 0.96, 0.73, 0.40 and 0.30 for deactivation orders 0, 1, 2 and 3.' },
      { name: 'Example 21.2b — poisoning by A and R', v: { pat: 'plug', n: 2, kt: 100, CA0: 50, dtype: 'sum', d: 1, kd: 1e-3, m: 1, trun: 120, pore: 'none', pol: 'K', treg: 0, _tu: 'day' }, book: { Xbar: { v: 0.73, tol: 0.01 } } },
      { name: 'Example 21.2c — second-order deactivation', v: { pat: 'plug', n: 2, kt: 100, CA0: 50, dtype: 'indep', d: 2, kd: 3.325, trun: 120, pore: 'none', pol: 'K', treg: 0, _tu: 'day' }, book: { Xbar: { v: 0.4, tol: 0.02 } } },
      { name: 'Example 21.2d — third-order deactivation', v: { pat: 'plug', n: 2, kt: 100, CA0: 50, dtype: 'indep', d: 3, kd: 666.5, trun: 120, pore: 'none', pol: 'K', treg: 0, _tu: 'day' }, book: { Xbar: { v: 0.3, tol: 0.02 } } },
      { name: 'Problem 21.15 — 12-day run', v: { pat: 'plug', n: 2, kt: 100, CA0: 0.05, dtype: 'sum', d: 1, kd: 10, m: 1, trun: 12, pore: 'none', pol: 'K', treg: 0, _tu: 'day' }, note: 'k′ = 200 L²/mol·hr·g, τ′ = 10 g·hr/L, C_A0 = 0.05 mol/L: k′τ′C_A0 = 100; −da/dt = 10(C_A + C_R)a = 0.5a per day.' },
      { name: 'Problem 21.7 — temperature raised from 600 to 800 K', v: { pat: 'plug', n: 2, kt: 10, CA0: 1, dtype: 'indep', d: 1, kd: 0.1, trun: 28, pore: 'none', pol: 'M', phase: 'gas', T0: 600, Tmax: 800, ER: 7200, EdR: 0, treg: 0, _tu: 'day' }, note: 'Second-order gas reaction with k = k₀e^{−7200/T}: holding k′a/T² constant (Eq. 21.48) until 800 K leaves an activity of (800/600)² e^{−3} = 0.089 at reactivation, whatever the deactivation rate.' },
      { name: 'Problem 21.8 — strong pores, 88% falling to 64%', v: { pat: 'plug', n: 1, kt: 2.12, CA0: 1, dtype: 'indep', d: 1, kd: 0.005840, trun: 383, pore: 'strong', pol: 'K', treg: 0, _tu: 'day' }, note: 'First order in the strong pore regime: ln ln(1/(1 − X)) falls by k_d t/2. From 88% to 64% in 250 days gives k_d = 0.00584 per day; 50% comes after about 383 days.' },
      { name: 'Problem 21.9 — the same catalyst as fine particles', v: { pat: 'plug', n: 1, kt: 2.12, CA0: 1, dtype: 'indep', d: 1, kd: 0.005840, trun: 125, pore: 'none', pol: 'K', treg: 0, _tu: 'day' }, note: 'Without pore resistance the decay is felt at full speed: 88% to 64% in about 125 days.' },
      { name: 'Problem 21.10 — fine particles in mixed flow', v: { pat: 'mixed', n: 1, kt: 7.333, CA0: 1, dtype: 'indep', d: 1, kd: 0.005840, trun: 243, pore: 'none', pol: 'K', treg: 0, _tu: 'day' }, note: 'Mixed flow: ln(X/(1 − X)) falls by k_d t: about 243 days from 88% to 64%.' },
      { name: 'Problem 21.11 — third-order deactivation, run to 10% activity', v: { pat: 'plug', n: 1, kt: 3, CA0: 1, dtype: 'indep', d: 3, kd: 3, trun: 16.5, pore: 'strong', pol: 'K', treg: 0, _tu: 'day' }, note: 'a = (1 + 6t)^{−1/2} reaches 0.1 after 16.5 days. For a first-order rate k′τ′ = k′W/v = 0.030 × 10 000/100 = 3 with W in grams and v in L/min; under strong pore resistance the rate falls only as a^½.' },
      { name: 'Poisoning front: parallel deactivation in a packed bed', v: { pat: 'plug', n: 1, kt: 4, CA0: 10, dtype: 'par', d: 1, kd: 0.02, m: 1, trun: 60, pore: 'none', pol: 'K', treg: 5, _tu: 'day' }, note: 'A poisons the catalyst where A is richest: activity dies first at the inlet and the active zone moves down the bed.' }
    ],
    eqs: [
      ['\\bar X_A = \\frac{1}{t_{run}}\\int_0^{t_{run}}X_A(t)\\,dt', 'Eq. 21.40'],
      ['\\mathbf{a} = e^{-k_dt}\\ (d = 1),\\qquad \\mathbf{a} = [1+(d-1)k_dt]^{1/(1-d)}\\ (d\\neq 1)', 'Eq. 21.39'],
      ['\\frac{W}{F_{A0}} = \\int_0^{X_A}\\frac{dX_A}{k^\\prime\\mathbf{a}C_A^n}', 'Eq. 21.43', 'Activity uniform through the bed.'],
      ['\\frac{F_{A0}}{\\mathbf{a}} = \\text{const}\\ \\text{(policy L)},\\qquad \\frac{k^\\prime\\mathbf{a}}{T^n} = \\text{const}\\ \\text{(policy M, ideal gas)}', 'Eqs. 21.45, 21.48', 'Liquids: k′a = const (Eq. 21.47).'],
      ['k^\\prime\\mathbf{a}\\mathscr{D}_e = \\text{const}\\ \\text{(liquid)},\\qquad \\frac{k^\\prime\\mathbf{a}\\mathscr{D}_e}{T^{n+1}} = \\text{const}\\ \\text{(gas)}', 'Eqs. 21.49–21.50', 'Policy M with strong pore resistance.'],
      ['\\ln\\frac{C_{A0}}{C_A} = \\frac{k^\\prime\\tau^\\prime}{M_T}\\exp\\left(-\\frac{k_dt}{2}\\right)', 'Eq. 21.43 (Table 21.1)', 'First order, strong pore resistance, plug flow.']
    ],
    method: ['With concentration-independent decay the activity is the same everywhere and the conversion at each instant comes from the exact plug or mixed equation for order n with the group k′aτ′C_A0ⁿ⁻¹. When decay depends on C_A or C_R the bed is divided into 60 cells, each with its own activity; at each time step the steady concentration profile is solved cell by cell and each cell’s activity is advanced with its local concentration (pseudo-steady state), which shows the poisoning front moving through the bed.', 'Pore resistance multiplies the rate by ℰ: 1/(M_T a^½) in the strong limit, or the exact first-order effectiveness for the chosen shape. Policy L lowers the feed so that the group stays constant; policy M raises the temperature to keep k′a/Tⁿ (gas) or k′a (liquid) constant, multiplied by D_e and one more power of T with strong pore resistance, while the decay speeds up with E_d, until T* ends the run. The best run length maximizes the mean conversion per cycle, t_run + t_reg.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'day';
      const n = +v.n, g0 = +v.kt, d = +v.d, kd = +v.kd, trun = +v.trun, CA0 = +v.CA0, m = +v.m;
      const plug = v.pat === 'plug';
      const etaF = a => (v.pore === 'strong' ? 1 / Math.sqrt(Math.max(a, 1e-300)) : v.pore === 'exact' ? PEL.eta1(v.shape, +v.MT0 * Math.sqrt(Math.max(a, 1e-300))) / PEL.eta1(v.shape, +v.MT0) : 1);
      const Xof = a => { const g = g0 * a * etaF(a); return plug ? CAT.Xplug(n, g) : CAT.Xmix(n, g); };
      const nt = 241, ts = N.linspace(0, trun, nt);
      let Xs, As, front = null;
      if (v.pol === 'K' && (v.dtype === 'par' || v.dtype === 'ser') && plug) {
        /* deactivation front: cells along the bed */
        const Nc = 60, a = new Float64Array(Nc).fill(1);
        Xs = []; As = []; front = [];
        const dt = trun / (nt - 1);
        const profile = () => {
          let X = 0; const Cloc = [];
          for (let i = 0; i < Nc; i++) {
            const gi = g0 * a[i] * etaF(a[i]) / Nc;
            const Xin = X;
            let Xo;
            if (Math.abs(n - 1) < 1e-12) Xo = 1 - (1 - Xin) * Math.exp(-gi);
            else if (n < 1) { const q = Math.pow(1 - Xin, 1 - n) - (1 - n) * gi; Xo = q <= 0 ? 1 : 1 - Math.pow(q, 1 / (1 - n)); }
            else Xo = 1 - Math.pow(Math.pow(1 - Xin, 1 - n) + (n - 1) * gi, 1 / (1 - n));
            const Xm = 0.5 * (Xin + Xo);
            Cloc.push(v.dtype === 'par' ? CA0 * (1 - Xm) : CA0 * Xm);
            X = Xo;
          }
          return { X, Cloc };
        };
        for (let j = 0; j < nt; j++) {
          const P = profile();
          Xs.push(P.X); As.push(N.mean(Array.from(a)));
          if (j % 60 === 0 || j === nt - 1) front.push({ t: ts[j], a: Array.from(a) });
          if (j === nt - 1) break;
          /* RK2 on each cell with its local concentration (held over the step) */
          for (let i = 0; i < Nc; i++) {
            const r = c => kd * Math.pow(Math.max(P.Cloc[i], 0), m) * Math.pow(Math.max(c, 0), d);
            const k1 = r(a[i]), k2 = r(a[i] - dt * k1 / 2);
            a[i] = Math.max(a[i] - dt * k2, 0);
          }
        }
      } else {
        const kdEff = v.dtype === 'sum' ? kd * Math.pow(CA0, m) : kd;
        if ((v.dtype === 'par' || v.dtype === 'ser') && !plug) {
          /* mixed flow: uniform activity, decay with the exit concentration */
          As = [1]; Xs = [Xof(1)];
          const dt = trun / (nt - 1);
          for (let j = 1; j < nt; j++) {
            const aPrev = As[j - 1], X = Xs[j - 1], C = v.dtype === 'par' ? CA0 * (1 - X) : CA0 * X;
            const aN = Math.max(aPrev - dt * kd * Math.pow(C, m) * Math.pow(aPrev, d), 0);
            As.push(aN); Xs.push(Xof(aN));
          }
        } else if (v.pol === 'M') {
          /* policy M: hold the performance group by raising T */
          const T0 = +v.T0, Tm = +v.Tmax, ER = +v.ER, EdR = +v.EdR, EDR = +v.EDR || 0;
          const gas = v.phase === 'gas', strong = v.pore !== 'none';
          const pw = gas ? (strong ? n + 1 : n) : 0;
          const lhs = T => Math.exp(-ER / T + ER / T0) * (strong ? Math.exp(-EDR / T + EDR / T0) : 1) / Math.pow(T / T0, pw);
          const Tof = a => (a >= 1 ? T0 : N.brent(T => lhs(T) * a - 1, T0, 50 * T0, { tol: 1e-10 }));
          const aEnd = 1 / lhs(Tm);
          const sol = N.ode((t, y) => { const T = Tof(Math.min(y[0], 1)); return [-kd * Math.exp(-EdR / T + EdR / T0) * Math.pow(Math.max(y[0], 0), d)]; }, 0, [1], trun, { rtol: 1e-9, atol: 1e-14, event: (t, y) => y[0] - aEnd, minPts: 200 });
          const tEnd = sol.t[sol.t.length - 1];
          const reached = sol.status === 'event';
          kpi.push({ k: 'aEnd', l: 'Activity when T reaches T*', v: aEnd, main: true, s: (gas ? 'ideal gas: k′a/T^' + fmt(pw, 3) : 'liquid: k′a') + (strong ? '·D_e' : '') + ' held constant' });
          kpi.push({ k: 'tEnd', l: reached ? 'Run ends at T* after' : 'T still below T* at the end of the run, t =', v: tEnd, u: tu });
          kpi.push({ k: 'X', l: 'Conversion held throughout', v: Xof(1) });
          const tt = sol.t, aa = sol.y.map(y => y[0]);
          plots.push({ title: 'Policy M: activity and temperature during the run', x: { label: 't', unit: tu, min: 0 }, y: { label: 'a', min: 0, max: 1 }, series: [{ name: 'activity', x: tt, y: aa, color: 1, w: 2.2 }, { name: '(T − T₀)/(T* − T₀)', x: tt, y: aa.map(q => (Tof(Math.min(q, 1)) - T0) / (Tm - T0)), color: 2, dash: '5 4' }], cap: 'Raising T to hold conversion is the usual optimum when decay is less temperature-sensitive than reaction (Fig. 21.9).' });
          return { kpi, warn, info, plots, tables };
        } else {
          As = ts.map(t => CAT.act(kdEff, d, t));
          Xs = As.map(Xof);
        }
      }
      const Xbar = N.trapz(ts, Xs) / trun;
      if (v.pol === 'L') {
        /* flow lowered: feed ∝ a·ℰ(a) for first-order-like scaling, conversion held */
        const thr = N.trapz(ts, As.map(a => a * etaF(a))) / trun;
        kpi.push({ k: 'thr', l: 'Mean throughput / initial throughput (policy L)', v: thr, main: true, s: 'conversion held at ' + fmt(Xs[0], 4) });
        kpi.push({ k: 'aEnd', l: 'Activity at the end of the run', v: As[As.length - 1] });
        plots.push({ title: 'Policy L: feed rate needed to hold conversion', x: { label: 't', unit: tu, min: 0 }, y: { label: 'F_A0/F_A0,initial', min: 0, max: 1 }, series: [{ name: 'feed rate', x: ts, y: As.map(a => a * etaF(a)), color: 1, w: 2.2 }] });
        return { kpi, warn, info, plots, tables };
      }
      kpi.push({ k: 'Xbar', l: 'Mean conversion over the run (Eq. 21.40)', v: Xbar, main: true, s: 'from ' + fmt(Xs[0], 4) + ' to ' + fmt(Xs[Xs.length - 1], 4) });
      kpi.push({ k: 'X0', l: 'Conversion at the start', v: Xs[0] }, { k: 'Xend', l: 'Conversion at the end', v: Xs[Xs.length - 1] });
      kpi.push({ k: 'aEnd', l: 'Activity at the end' + (front ? ' (bed average)' : ''), v: As[As.length - 1] });
      /* best run length with regeneration */
      const treg = +v.treg;
      if (treg > 0 || true) {
        const cum = N.cumtrapz(ts, Xs);
        const yld = ts.map((t, i) => (t > 0 ? cum[i] / (t + treg) : (treg > 0 ? 0 : Xs[0])));
        let bi = 1; yld.forEach((q, i) => { if (i > 0 && q > yld[bi]) bi = i; });
        if (treg > 0) {
          kpi.push({ k: 'tbest', l: 'Run length giving the most product per cycle', v: ts[bi], u: tu, s: bi === nt - 1 ? 'at or beyond the run length entered: run longer' : 'mean conversion per cycle ' + fmt(yld[bi], 4) });
          plots.push({ title: 'Product per unit of cycle time against run length', x: { label: 'run length', unit: tu, min: 0 }, y: { label: '∫X dt/(t_run + t_reg)', min: 0 }, series: [{ name: 't_reg = ' + fmt(treg, 3) + ' ' + tu, x: ts, y: yld, color: 1, w: 2.2 }], marks: [{ t: 'pt', x: ts[bi], y: yld[bi], color: 8, label: 'best' }] });
        }
      }
      plots.push({ title: 'Conversion and activity during the run (Fig. E21.2)', x: { label: 't', unit: tu, min: 0, max: trun }, y: { label: 'X_A, a', min: 0, max: 1 }, series: [{ name: 'X_A', x: ts, y: Xs, color: 1, w: 2.2 }, { name: front ? 'activity, bed average' : 'activity a', x: ts, y: As, color: 2, dash: '5 4' }, { name: 'mean X̄_A', x: [0, trun], y: [Xbar, Xbar], color: 'var(--muted)', dash: '2 3' }] });
      if (front) {
        const zz = N.linspace(0.5 / 60, 1 - 0.5 / 60, 60);
        plots.push({ title: 'Activity along the bed: the deactivation front', x: { label: 'fraction of the bed from the inlet', min: 0, max: 1 }, y: { label: 'a', min: 0, max: 1 }, series: front.map((f, i) => ({ name: 't = ' + fmt(f.t, 3) + ' ' + tu, x: zz, y: f.a, color: (i % 7) + 1 })), cap: v.dtype === 'par' ? 'Parallel poisoning starts at the inlet, where A is richest, and the dead zone grows downstream.' : 'Series poisoning attacks the outlet first, where the product is richest.' });
      }
      if (v.pore === 'strong') info.push('With strong pore resistance the rate falls only as a^½: the catalyst seems to deactivate at half its true rate (Table 21.1).');
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
