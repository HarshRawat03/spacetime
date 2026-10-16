/* Chapter 8 — Potpourri of multiple reactions */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, KIT = ST.KIT, D = ST.D;

  /* ---------- series A -> R -> S, orders n1 (A) and n2 (R) ---------- */
  function seriesSys(k1, k2, n1, n2) {
    const P = (c, n) => (n === 0 ? (c > 0 ? Math.min(1, c * 1e6) : 0) : Math.pow(Math.max(c, 0), n));
    return { r1: ca => k1 * P(ca, n1), r2: cr => k2 * P(cr, n2) };
  }
  function seriesPFR(sys, CA0, tEnd, n = 200) {
    const sol = N.ode((t, y) => { const a = sys.r1(y[0]), b = sys.r2(y[1]); return [-a, a - b, b]; }, 0, [CA0, 0, 0], tEnd, { rtol: 1e-9, atol: 1e-12, minPts: 40 });
    const ts = N.linspace(0, tEnd, n);
    return { ts, y: ts.map(sol.at), sol };
  }
  function seriesMFR(sys, CA0, tau) {
    const CA = N.root(c => CA0 - c - tau * sys.r1(c), 0, CA0, { n: 200 });
    const gen = CA0 - CA;
    const CR = N.root(c => gen - c - tau * sys.r2(c), 0, gen, { n: 200 });
    return [CA, isFinite(CR) ? CR : 0, gen - (isFinite(CR) ? CR : 0)];
  }
  /* product distribution, first-order steps, kappa = k2/k1 */
  const crPFR = (x, ka) => (Math.abs(ka - 1) < 1e-9 ? x * Math.log(1 / x) : (Math.pow(x, ka) - x) / (1 - ka)); /* x = CA/CA0 */
  const crMFR = (x, ka) => x * (1 - x) / (x + ka * (1 - x));

  ST.tool({
    id: 'series', tag: 'H-801', ch: 8, sec: '8.1–8.4', title: 'Reactions in series A → R → S',
    lede: 'Concentration histories, the maximum of the intermediate and when it occurs, and product-distribution charts for plug, mixed and staged flow; first-order steps exactly, zero-order and other orders numerically. Also works backwards: k₂/k₁ from one measured composition.',
    keys: 'series reactions consecutive intermediate maximum optimum residence time kappa k2/k1 product distribution zero order first order figure 8.3 8.5 8.6 8.7 8.8',
    units: { C: 'mol/L', t: 'min' },
    inputs: [
      { group: 'Kinetics' },
      { k: 'k1', label: 'k₁ (A → R)', v: 0.1, min: 1e-5, max: 1e4, scale: 'log' },
      { k: 'n1', label: 'Order of step 1', v: 1, min: 0, max: 3, step: 0.5 },
      { k: 'k2', label: 'k₂ (R → S)', v: 0.05, min: 0, max: 1e4, scale: 'log' },
      { k: 'n2', label: 'Order of step 2', v: 1, min: 0, max: 3, step: 0.5 },
      { k: 'CA0', label: '$C_{A0}$', unit: '[C]', v: 1, min: 1e-6, max: 1e6, scale: 'log' },
      { k: 'Nt', label: 'Equal tanks for the staged case', v: 3, min: 2, max: 30, step: 1 },
      { group: 'Work backwards (first-order steps)' },
      { k: 'inv', t: 'check', label: 'Find k₂/k₁ from one measured composition', v: false },
      { k: 'invR', t: 'seg', label: 'Measured in', opts: [['pfr', 'batch or plug flow'], ['mfr', 'mixed flow']], v: 'mfr', if: v => v.inv },
      { k: 'obsA', label: 'Measured C_A/C_{A0}', v: 0.5, min: 0.001, max: 0.999, step: 0.001, if: v => v.inv },
      { k: 'obsR', label: 'Measured C_R/C_{A0}', v: 0.333, min: 0.0001, max: 0.999, step: 0.001, if: v => v.inv },
      { k: 'v0', label: 'Feed rate (for volumes)', unit: v => KIT.volUnit(v._cu) + '/' + v._tu, v: 1000 / 60, min: 1e-9, max: 1e9, scale: 'log' }
    ],
    presets: [
      { name: 'Problem 8.2', v: { k1: 0.1, n1: 1, k2: 0.1, n2: 1, CA0: 1, v0: 1000 / 60, inv: false, _cu: 'mol/L', _tu: 'min' }, note: 'k₁ = k₂ = 0.1 min⁻¹, 1000 L/h: the plug reactor peaks at 1/e of C_{A0}, the mixed one at 1/4.' },
      { name: 'Problem 8.4', v: { k1: 0.1, n1: 1, k2: 0.05, n2: 1, CA0: 100, inv: true, invR: 'mfr', obsA: 0.5, obsR: 1 / 3 }, note: 'Two mixed flow runs from 100: (50, 33⅓, 16⅔) and (25, 30, 45). Enter either run; both give the same k₂/k₁ if A → R → S is the right scheme.' },
      { name: 'First then zero order', v: { k1: 0.2, n1: 1, k2: 0.05, n2: 0, CA0: 1, inv: false }, note: '§8.2: R is removed at a constant rate while any is present.' },
      { name: 'Zero then first order', v: { k1: 0.1, n1: 0, k2: 0.2, n2: 1, CA0: 1, inv: false }, note: '§8.3: A disappears at a constant rate until it runs out.' }
    ],
    eqs: [
      ['\\frac{C_R}{C_{A0}} = \\frac{k_1}{k_2-k_1}\\left(e^{-k_1\\tau}-e^{-k_2\\tau}\\right)', 'Eq. 8.7', 'Plug flow or batch.'],
      ['\\frac{C_{R,max}}{C_{A0}} = \\left(\\frac{k_1}{k_2}\\right)^{k_2/(k_2-k_1)},\\qquad \\tau_{p,opt} = \\frac{\\ln(k_2/k_1)}{k_2-k_1}', 'Eqs. 8.8–8.9'],
      ['\\frac{C_R}{C_{A0}} = \\frac{k_1\\tau}{(1+k_1\\tau)(1+k_2\\tau)}', 'Eq. 8.13', 'Mixed flow.'],
      ['\\frac{C_{R,max}}{C_{A0}} = \\frac{1}{\\left[(k_2/k_1)^{1/2}+1\\right]^2},\\qquad \\tau_{m,opt} = \\frac{1}{\\sqrt{k_1k_2}}', 'Eqs. 8.15–8.16'],
      ['\\frac{C_R}{C_{A0}} = \\frac{1}{1-\\kappa}\\left[\\left(\\frac{C_A}{C_{A0}}\\right)^{\\kappa}-\\frac{C_A}{C_{A0}}\\right],\\quad \\kappa = \\frac{k_2}{k_1}', 'Eq. 8.37 form', 'Time-free product distribution, plug flow.']
    ],
    run(v) {
      const k1 = +v.k1, k2 = +v.k2, n1 = +v.n1, n2 = +v.n2, CA0 = +v.CA0;
      const sys = seriesSys(k1, k2, n1, n2);
      const kpi = [], warn = [], plots = [], tables = [];
      const first = n1 === 1 && n2 === 1;
      /* plug flow */
      let tP, cP, tM, cM;
      if (first) {
        if (Math.abs(k2 - k1) < 1e-12 * k1) { tP = 1 / k1; cP = CA0 / Math.E; } else if (k2 === 0) { tP = Infinity; cP = CA0; } else { tP = Math.log(k2 / k1) / (k2 - k1); cP = CA0 * Math.pow(k1 / k2, k2 / (k2 - k1)); }
        tM = k2 > 0 ? 1 / Math.sqrt(k1 * k2) : Infinity; cM = k2 > 0 ? CA0 / Math.pow(Math.sqrt(k2 / k1) + 1, 2) : CA0;
      } else {
        const horizon = 20 / Math.max(Math.min(k1 * Math.pow(CA0, n1 - 1), k2 > 0 ? k2 * Math.pow(CA0, Math.max(n2 - 1, 0)) : k1), 1e-9);
        const run = seriesPFR(sys, CA0, horizon, 800);
        let bi = 0; run.y.forEach((y, i) => { if (y[1] > run.y[bi][1]) bi = i; });
        const o = N.fmin(t => -run.sol.at(t)[1], run.ts[Math.max(0, bi - 1)], run.ts[Math.min(run.ts.length - 1, bi + 1)]);
        tP = o.x; cP = -o.f;
        const om = N.fminScan(t => -seriesMFR(sys, CA0, t)[1], horizon * 1e-4, horizon, 120, { log: true });
        tM = om.x; cM = -om.f;
      }
      kpi.push({ k: 'CRmaxP', l: 'C_{R,max}, plug or batch', v: cP, u: v._cu, main: true }, { k: 'tP', l: 'τ for it, plug', v: tP, u: v._tu }, { k: 'CRmaxM', l: 'C_{R,max}, mixed', v: cM, u: v._cu }, { k: 'tM', l: 'τ for it, mixed', v: tM, u: v._tu });
      const v0 = +v.v0;
      if (isFinite(v0)) kpi.push({ k: 'VP', l: 'Plug volume at the optimum', v: tP * v0, u: KIT.volUnit(v._cu) }, { k: 'VM', l: 'Mixed volume at the optimum', v: tM * v0, u: KIT.volUnit(v._cu) });
      /* staged: N equal tanks, optimize each tau */
      const Nt = Math.round(+v.Nt);
      const stage = te => { let c = [CA0, 0, 0]; for (let i = 0; i < Nt; i++) { const CA = N.root(x => c[0] - x - te * sys.r1(x), 0, c[0], { n: 120 }); const gen = c[0] - CA + c[1]; const CR = N.root(x => gen - x - te * sys.r2(x), 0, gen, { n: 120 }); c = [CA, CR, CA0 - CA - CR]; } return c; };
      const os = N.fminScan(t => -stage(t)[1], (tM / Nt) * 0.02, tM * 2, 60, { log: true });
      kpi.push({ k: 'CRmaxN', l: 'C_{R,max}, ' + Nt + ' equal tanks', v: -os.f, u: v._cu, s: 'total τ = ' + N.fmt(os.x * Nt, 3) });
      /* histories */
      const hor = Math.max(3 * tP, 3 * tM, 1e-6);
      const run = seriesPFR(sys, CA0, isFinite(hor) ? hor : 50 / k1);
      const ts = run.ts;
      plots.push({ title: 'Batch or plug flow (Fig. 8.3a)', x: { label: 'τ', unit: v._tu, min: 0 }, y: { label: 'C', unit: v._cu, min: 0 }, series: [{ name: 'C_A', x: ts, y: run.y.map(y => y[0]), color: 1 }, { name: 'C_R', x: ts, y: run.y.map(y => y[1]), color: 2 }, { name: 'C_S', x: ts, y: run.y.map(y => y[2]), color: 3 }], marks: [{ t: 'pt', x: tP, y: cP, color: 2, label: 'max' }] });
      const tm = N.logspace(Math.max(tM, 1e-6) * 0.01, Math.max(tM, 1e-6) * 100, 160), cm = tm.map(t => seriesMFR(sys, CA0, t));
      plots.push({ title: 'Mixed flow (Fig. 8.5a)', x: { label: 'τ', unit: v._tu, log: true }, y: { label: 'C', unit: v._cu, min: 0 }, series: [{ name: 'C_A', x: tm, y: cm.map(c => c[0]), color: 1 }, { name: 'C_R', x: tm, y: cm.map(c => c[1]), color: 2 }, { name: 'C_S', x: tm, y: cm.map(c => c[2]), color: 3 }], marks: [{ t: 'pt', x: tM, y: cM, color: 2, label: 'max' }] });
      /* product distribution charts */
      if (first) {
        const ka = k2 / k1, xs = N.linspace(0.001, 1, 200);
        const fam = [0.1, 0.25, 0.5, 1, 2, 4];
        const series = [];
        fam.forEach((f, i) => series.push({ name: 'κ = ' + f, x: xs.map(x => 1 - x), y: xs.map(x => crPFR(x, f)), color: 8, w: 1, legend: false, hover: false }));
        series.push({ name: 'plug, κ = ' + N.fmt(ka, 3), x: xs.map(x => 1 - x), y: xs.map(x => crPFR(x, ka)), color: 2, w: 2.6 }, { name: 'mixed, κ = ' + N.fmt(ka, 3), x: xs.map(x => 1 - x), y: xs.map(x => crMFR(x, ka)), color: 6, w: 2.6, dash: '6 4' });
        const marks = fam.map(f => { const x = 0.5; return { t: 'text', x: 1 - 0.12, y: crPFR(0.12, f), s: String(f), color: 8, noRange: true }; });
        plots.push({ title: 'Product distribution, time-free (Figs. 8.7–8.8)', x: { label: 'X_A = 1 − C_A/C_{A0}', min: 0, max: 1 }, y: { label: 'C_R/C_{A0}', min: 0, max: 1 }, series, marks });
        const phiP = xs.map(x => crPFR(x, ka) / (1 - x)), phiM = xs.map(x => crMFR(x, ka) / (1 - x));
        plots.push({ title: 'Fractional yield of R (Fig. 8.6)', x: { label: 'X_A', min: 0, max: 1 }, y: { label: 'C_R/(C_{A0} − C_A)', min: 0, max: 1 }, series: [{ name: 'plug', x: xs.map(x => 1 - x), y: phiP, color: 2 }, { name: 'mixed', x: xs.map(x => 1 - x), y: phiM, color: 6, dash: '6 4' }] });
      }
      if (v.inv) {
        const x = +v.obsA, r = +v.obsR;
        const f = v.invR === 'pfr' ? ka => crPFR(x, ka) - r : ka => crMFR(x, ka) - r;
        const ka = N.root(f, 1e-6, 1e4, { n: 400, log: true });
        kpi.unshift({ k: 'kappa', l: 'k₂/k₁ from the measurement', v: ka, main: true, s: isFinite(ka) ? '' : 'no κ matches: check the scheme' });
        if (!isFinite(ka)) warn.push('No value of k₂/k₁ reproduces that composition for first-order A → R → S in ' + (v.invR === 'pfr' ? 'plug flow' : 'mixed flow') + '. A different scheme (parallel steps, other orders) is needed.');
      }
      const rows = [];
      N.linspace(0, isFinite(hor) ? hor : 10, 9).forEach(t => { const y = run.sol.at(t), m = seriesMFR(sys, CA0, Math.max(t, 1e-12)); rows.push([t, y[0], y[1], y[2], m[0], m[1], m[2]]); });
      tables.push({ title: 'Plug against mixed flow at the same τ', cols: ['τ', 'C_A plug', 'C_R plug', 'C_S plug', 'C_A mixed', 'C_R mixed', 'C_S mixed'], rows });
      return { kpi, warn, plots, tables };
    }
  });

  /* ---------- series-parallel A + B -> R, R + B -> S ---------- */
  ST.tool({
    id: 'seriespar', tag: 'H-802', ch: 8, sec: '8.6', title: 'Series–parallel reactions A + B → R, R + B → S',
    lede: 'Product distribution charts of Figs. 8.13 and 8.14 for any k₂/k₁, the effect of adding B slowly or all at once, a third step (benzene chlorination), and k₂/k₁ worked back from any one observation.',
    keys: 'series parallel reactions competitive consecutive A+B R+B figure 8.13 8.14 k2/k1 mixing beaker chlorination benzene ethylene oxide glycol',
    units: { C: 'mol/L', t: 'min' },
    inputs: [
      { group: 'Kinetics' },
      { k: 'ka', label: 'k₂/k₁', v: 0.5, min: 0.001, max: 100, scale: 'log' },
      { k: 'three', t: 'check', label: 'Third step U + B → … (k₃/k₂)', v: false },
      { k: 'kb', label: 'k₃/k₂', v: 1 / 240, min: 1e-5, max: 100, scale: 'log', if: v => v.three },
      { group: 'Operating point' },
      { k: 'how', t: 'seg', label: 'Contacting', opts: [['pfr', 'batch or plug in A'], ['mfr', 'mixed flow']], v: 'pfr' },
      { k: 'dB', label: 'B consumed per A fed, ΔB/A₀', v: 1, min: 0, max: 3, step: 0.01 },
      { group: 'Work backwards' },
      { k: 'inv', t: 'select', label: 'Observation', opts: [['none', 'none'], ['AdB', 'A left and B consumed'], ['AR', 'A left and R present'], ['dBR', 'B consumed and R present'], ['dBS', 'B consumed and S present']], v: 'none' },
      { k: 'o1', label: 'first value (per mole A fed)', v: 0.67, min: 0, max: 3, step: 0.001, if: v => v.inv !== 'none' },
      { k: 'o2', label: 'second value (per mole A fed)', v: 0.5, min: 0, max: 3, step: 0.001, if: v => v.inv !== 'none' }
    ],
    presets: [
      { name: 'Example 8.2(a)', v: { ka: 4, how: 'pfr', inv: 'AdB', o1: 0.67, o2: 0.5 }, book: { kinv: { v: 4, tol: 0.08 } }, note: 'Half a mole of B poured in slowly onto a mole of A; when B is used up 0.67 mol A remain.' },
      { name: 'Example 8.2(b)', v: { ka: 0.4, how: 'pfr', inv: 'dBR', o1: 1.25, o2: 0.5 }, book: { kinv: { v: 0.4, tol: 0.08 } }, note: 'One mole A and 1.25 mol B mixed fast; at the end 0.5 mol R.' },
      { name: 'Example 8.2(c)', v: { ka: 1.45, how: 'pfr', inv: 'dBS', o1: 0.9, o2: 0.3 }, book: { kinv: { v: 1.45, tol: 0.08 } }, note: 'One mole A and 1.25 mol B; when 0.9 mol B has reacted, 0.3 mol S is present.' },
      { name: 'Fig. 8.15 chlorination', v: { ka: 1 / 8, three: true, kb: 1 / 240, how: 'pfr', dB: 1, inv: 'none' }, note: 'Benzene + Cl₂ to mono-, di- and trichlorobenzene with k₂/k₁ = 1/8 and k₃/k₂ = 1/240 (MacMullin).' }
    ],
    eqs: [
      ['\\frac{dC_R}{dC_A} = -1 + \\frac{k_2}{k_1}\\frac{C_R}{C_A}', 'Eq. 8.34', 'B cancels: the distribution does not depend on how B is added, as long as A and R are not mixed back.'],
      ['\\frac{C_R}{C_{A0}} = \\frac{1}{1-k_2/k_1}\\left[\\left(\\frac{C_A}{C_{A0}}\\right)^{k_2/k_1}-\\frac{C_A}{C_{A0}}\\right]', 'Eq. 8.37', 'Batch or plug flow in A.'],
      ['\\frac{C_R}{C_{A0}} = \\frac{C_A(C_{A0}-C_A)}{C_{A0}\\left[C_A+(k_2/k_1)(C_{A0}-C_A)\\right]}', 'Eq. 8.41', 'Mixed flow.'],
      ['\\Delta C_B = (C_{A0}-C_A) + C_S = 2(C_{A0}-C_A) - C_R', 'Eq. 8.38–8.39']
    ],
    run(v) {
      const ka = +v.ka, how = v.how, three = !!v.three, kb = +v.kb * ka;
      const kpi = [], plots = [], warn = [];
      /* trajectory parameterized by x = C_A/C_A0 */
      const crOf = x => (how === 'pfr' ? crPFR(x, ka) : crMFR(x, ka));
      let traj;
      if (!three) traj = x => { const r = crOf(x), s = 1 - x - r; return { a: x, r, s, u: 0, dB: (1 - x) + s }; };
      else {
        /* three steps: A -> R (1), R -> S (ka), S -> U (kb), B taking part in each; integrate in x = C_A/C_A0 */
        if (how === 'pfr') {
          const sol = N.ode((x, y) => { const r = y[0], s2 = y[1]; return [-(1 - ka * r / x), -(ka * r - kb * s2) / x]; }, 1, [0, 0], 1e-6, { rtol: 1e-9, atol: 1e-12 });
          traj = x => { const y = sol.at(x), r = y[0], s2 = y[1], u = Math.max(0, 1 - x - r - s2); return { a: x, r, s: s2, u, dB: (1 - x) + s2 + 2 * u }; };
        } else {
          traj = x => { const th = (1 - x) / x, R = th * x / (1 + ka * th), S = ka * th * R / (1 + kb * th), U = Math.max(0, 1 - x - R - S); return { a: x, r: R, s: S, u: U, dB: (1 - x) + S + 2 * U }; };
        }
      }
      /* find x at the given dB */
      const dB = +v.dB;
      const xAt = N.root(x => traj(x).dB - dB, 1e-6, 1, { n: 400 });
      if (isFinite(xAt)) {
        const p = traj(xAt);
        kpi.push({ k: 'A', l: 'A left per A fed', v: p.a, main: true }, { k: 'R', l: 'R per A fed', v: p.r }, { k: 'S', l: 'S per A fed', v: p.s });
        if (three) kpi.push({ k: 'U', l: 'U per A fed', v: p.u });
      } else warn.push('That much B cannot be consumed: the reaction stops when A, R' + (three ? ' and S' : '') + ' are used up.');
      /* inverse */
      if (v.inv !== 'none') {
        const o1 = +v.o1, o2 = +v.o2;
        let A, R;
        if (v.inv === 'AdB') { A = o1; R = 2 * (1 - A) - o2; }
        else if (v.inv === 'AR') { A = o1; R = o2; }
        else if (v.inv === 'dBR') { R = o2; A = 1 - (o1 + R) / 2; }
        else { const S = o2; A = 1 - (o1 - S); R = 1 - A - S; }
        const f = k => (how === 'pfr' ? crPFR(A, k) : crMFR(A, k)) - R;
        const kinv = N.root(f, 1e-5, 1e4, { n: 500, log: true });
        kpi.unshift({ k: 'kinv', l: 'k₂/k₁ from the observation', v: kinv, main: true, s: 'A = ' + N.fmt(A, 3) + ', R = ' + N.fmt(R, 3) + ' per A fed' });
        if (!isFinite(kinv)) warn.push('No k₂/k₁ fits that observation for this contacting.');
      }
      /* charts */
      const xsF = N.linspace(1, 0.0005, 260);
      const fam = [0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10];
      const series = fam.map(f => ({ x: xsF, y: xsF.map(x => (how === 'pfr' ? (f === 0 ? 1 - x : crPFR(x, f)) : (f === 0 ? 1 - x : crMFR(x, f)))), color: 8, w: 1, legend: false, hover: false }));
      const tx = xsF.map(traj);
      series.push({ name: 'R, k₂/k₁ = ' + N.fmt(ka, 3), x: xsF, y: tx.map(p => p.r), color: 2, w: 2.6 }, { name: 'S', x: xsF, y: tx.map(p => p.s), color: 3, w: 2 });
      if (three) series.push({ name: 'U', x: xsF, y: tx.map(p => p.u), color: 4, w: 2 });
      const marks = fam.filter(f => f > 0).map(f => ({ t: 'text', x: 0.08, y: (how === 'pfr' ? crPFR(0.08, f) : crMFR(0.08, f)), s: String(f), color: 8, noRange: true }));
      if (isFinite(xAt)) marks.push({ t: 'v', x: xAt, label: 'ΔB/A₀ = ' + N.fmt(dB, 3), color: 6 });
      plots.push({ title: (how === 'pfr' ? 'Fig. 8.13, batch or plug flow' : 'Fig. 8.14, mixed flow') + ': light lines are k₂/k₁ = ' + fam.filter(f => f > 0).join(', '), wide: true, h: 380, x: { label: 'C_A/C_{A0}', min: 0, max: 1 }, y: { label: 'C/C_{A0}', min: 0, max: 1 }, series, marks });
      const bs = N.linspace(0, 2, 200).map(b => { const x = N.root(xx => traj(xx).dB - b, 1e-6, 1, { n: 200 }); return isFinite(x) ? traj(x) : null; });
      const bx = N.linspace(0, 2, 200);
      plots.push({ title: 'Composition against B consumed', x: { label: 'ΔB/A₀', min: 0, max: 2 }, y: { label: 'per mole A fed', min: 0, max: 1 }, series: [{ name: 'A', x: bx, y: bs.map(p => (p ? p.a : NaN)), color: 1 }, { name: 'R', x: bx, y: bs.map(p => (p ? p.r : NaN)), color: 2 }, { name: 'S', x: bx, y: bs.map(p => (p ? p.s : NaN)), color: 3 }].concat(three ? [{ name: 'U', x: bx, y: bs.map(p => (p ? p.u : NaN)), color: 4 }] : []) });
      return { kpi, warn, plots };
    }
  });

  /* ---------- Denbigh scheme ---------- */
  ST.tool({
    id: 'denbigh', tag: 'H-803', ch: 8, sec: '8.7', title: 'The Denbigh reactions',
    lede: 'A →₁ R →₃ S with A →₂ T and R →₄ U, all first order: compositions in plug and mixed flow, the maximum of R and its space-time, and every special case by setting rate constants to zero.',
    keys: 'denbigh reaction scheme first order network special cases maximum intermediate plug mixed sodium sulfide',
    units: { C: 'mol/L', t: 'min' },
    inputs: [
      { group: 'Rate constants' },
      { k: 'k1', label: 'k₁ (A → R)', unit: '1/[t]', v: 1, min: 0, max: 1e4, scale: 'log' },
      { k: 'k2', label: 'k₂ (A → T)', unit: '1/[t]', v: 0.3, min: 0, max: 1e4, scale: 'log' },
      { k: 'k3', label: 'k₃ (R → S)', unit: '1/[t]', v: 0.5, min: 0, max: 1e4, scale: 'log' },
      { k: 'k4', label: 'k₄ (R → U)', unit: '1/[t]', v: 0.1, min: 0, max: 1e4, scale: 'log' },
      { k: 'CA0', label: '$C_{A0}$', unit: '[C]', v: 1, min: 1e-6, max: 1e6, scale: 'log' }
    ],
    presets: [
      { name: 'Series only (k₂ = k₄ = 0)', v: { k1: 1, k2: 0, k3: 0.5, k4: 0 } },
      { name: 'Parallel only (k₃ = k₄ = 0)', v: { k1: 1, k2: 0.5, k3: 0, k4: 0 } }
    ],
    eqs: [
      ['k_{12} = k_1+k_2,\\qquad k_{34} = k_3+k_4', '§8.7'],
      ['\\frac{C_{R,max}}{C_{A0}} = \\frac{k_1}{k_{12}}\\left(\\frac{k_{12}}{k_{34}}\\right)^{k_{34}/(k_{34}-k_{12})},\\quad \\tau_{p,opt} = \\frac{\\ln(k_{34}/k_{12})}{k_{34}-k_{12}}', 'Eqs. 8.51–8.52', 'Plug flow.'],
      ['\\frac{C_{R,max}}{C_{A0}} = \\frac{k_1}{\\left(k_{12}^{1/2}+k_{34}^{1/2}\\right)^2},\\quad \\tau_{m,opt} = \\frac{1}{\\sqrt{k_{12}k_{34}}}', 'Eqs. 8.56–8.57', 'Mixed flow.']
    ],
    run(v) {
      const k1 = +v.k1, k2 = +v.k2, k3 = +v.k3, k4 = +v.k4, C0 = +v.CA0;
      const k12 = k1 + k2, k34 = k3 + k4;
      if (!(k12 > 0)) throw new Error('A must react: k₁ + k₂ > 0.');
      const pf = t => { const A = C0 * Math.exp(-k12 * t); const R = Math.abs(k34 - k12) < 1e-12 ? C0 * k1 * t * Math.exp(-k12 * t) : C0 * k1 / (k34 - k12) * (Math.exp(-k12 * t) - Math.exp(-k34 * t)); const T = k2 / k12 * (C0 - A); const RS = C0 * k1 / k12 * (1 - Math.exp(-k12 * t)) - R; return [A, R, k34 > 0 ? RS * k3 / k34 : 0, T, k34 > 0 ? RS * k4 / k34 : 0]; };
      const mf = t => { const A = C0 / (1 + k12 * t), R = k1 * t * A / (1 + k34 * t); return [A, R, k3 * t * R, k2 * t * A, k4 * t * R]; };
      let tp, cp, tm, cm;
      if (k34 === 0) { tp = tm = Infinity; cp = cm = C0 * k1 / k12; }
      else {
        tp = Math.abs(k34 - k12) < 1e-12 ? 1 / k12 : Math.log(k34 / k12) / (k34 - k12);
        cp = Math.abs(k34 - k12) < 1e-12 ? C0 * k1 / k12 / Math.E : C0 * k1 / k12 * Math.pow(k12 / k34, k34 / (k34 - k12));
        tm = 1 / Math.sqrt(k12 * k34); cm = C0 * k1 / Math.pow(Math.sqrt(k12) + Math.sqrt(k34), 2);
      }
      const kpi = [{ k: 'cp', l: 'C_{R,max}, plug', v: cp, u: v._cu, main: true }, { k: 'tp', l: 'τ for it, plug', v: tp, u: v._tu }, { k: 'cm', l: 'C_{R,max}, mixed', v: cm, u: v._cu }, { k: 'tm', l: 'τ for it, mixed', v: tm, u: v._tu }];
      const hor = isFinite(tp) ? 4 * Math.max(tp, tm) : 6 / k12;
      const ts = N.linspace(0, hor, 200), tl = N.logspace(hor * 1e-3, hor * 10, 200);
      const names = ['A', 'R', 'S', 'T', 'U'];
      const mk = (f, xs) => names.map((n, i) => ({ name: 'C_' + n, x: xs, y: xs.map(t => f(t)[i]), color: i + 1 }));
      return {
        kpi,
        plots: [{ title: 'Plug flow (Fig. 8.17)', x: { label: 'τ', unit: v._tu, min: 0 }, y: { label: 'C', unit: v._cu, min: 0 }, series: mk(pf, ts), marks: [{ t: 'pt', x: tp, y: cp, color: 2, label: 'max' }] }, { title: 'Mixed flow (Fig. 8.18)', x: { label: 'τ', unit: v._tu, log: true }, y: { label: 'C', unit: v._cu, min: 0 }, series: mk(mf, tl), marks: [{ t: 'pt', x: tm, y: cm, color: 2, label: 'max' }] }],
        tables: [{ title: 'At the plug-flow optimum and the mixed-flow optimum', cols: ['', 'C_A', 'C_R', 'C_S', 'C_T', 'C_U'], rows: [['plug, τ = ' + N.fmt(tp, 3)].concat(isFinite(tp) ? pf(tp) : [NaN, NaN, NaN, NaN, NaN]), ['mixed, τ = ' + N.fmt(tm, 3)].concat(isFinite(tm) ? mf(tm) : [NaN, NaN, NaN, NaN, NaN])] }]
      };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
