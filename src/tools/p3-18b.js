/* Chapter 18 (II) — rate equations from catalytic reactors; reactor size with pore diffusion;
   product distribution when pores interfere with multiple reactions */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, PEL = ST.PEL, KIT = ST.KIT;
  const CAT = ST.CAT = ST.CAT || {};
  const fmt = (x, s = 4) => N.fmt(x, s);

  /* units of a rate constant per kg of catalyst for order n: vol^n mol^(1−n) / (kg·t) */
  CAT.kpU = function (n, cu, tu, per = 'kg') {
    const V = KIT.volUnit(cu || 'mol/L'), e = 1 - n, num = [], den = [];
    const pw = (b, p) => (Math.abs(p - 1) < 1e-9 ? b : b + '^{' + fmt(p, 3) + '}');
    if (n > 1e-9) num.push(pw(V, n));
    if (e > 1e-9) num.push(pw('mol', e));
    if (e < -1e-9) den.push(pw('mol', -e));
    den.push(per, tu || 'h');
    return (num.length ? num.join('·') : '1') + '/(' + den.join('·') + ')';
  };
  const Cof = (CA0, eps) => X => CA0 * (1 - X) / (1 + eps * X);
  const Xof = (CA0, eps) => C => (1 - C / CA0) / (1 + eps * C / CA0);

  /* candidate rate forms, rates per unit mass of catalyst */
  const MOD = {
    p0: { name: 'zero order, k′', p: ['k′'], pos: [true], f: (p, C) => (C > 0 ? p[0] : 0), order: () => 0 },
    p1: { name: 'first order, k′C', p: ['k′'], pos: [true], f: (p, C) => p[0] * Math.max(C, 0), order: () => 1 },
    p2: { name: 'second order, k′C²', p: ['k′'], pos: [true], f: (p, C) => p[0] * Math.max(C, 0) * Math.max(C, 0), order: () => 2 },
    pn: { name: 'nth order, k′Cⁿ', p: ['k′', 'n'], pos: [true, false], f: (p, C) => (C > 0 ? p[0] * Math.pow(C, p[1]) : 0), order: p => p[1] },
    lh: { name: 'saturation, k′C/(1 + KC)', p: ['k′', 'K'], pos: [true, true], f: (p, C) => (C > 0 ? p[0] * C / (1 + p[1] * C) : 0), order: () => 1 },
    lh2: { name: 'dual site, k′C/(1 + KC)²', p: ['k′', 'K'], pos: [true, true], f: (p, C) => (C > 0 ? p[0] * C / Math.pow(1 + p[1] * C, 2) : 0), order: () => 1 },
    p1r: { name: 'reversible, k′(C − C_e)', p: ['k′', 'C_e'], pos: [true, true], f: (p, C) => p[0] * (C - p[1]), order: () => 1 }
  };
  CAT.MOD = MOD;
  const aicc = (n, ssr, k) => n * Math.log(Math.max(ssr, 1e-300) / n) + 2 * k + (n - k - 1 > 0 ? 2 * k * (k + 1) / (n - k - 1) : 0);

  /* fit candidates to direct rate data (C, r) */
  function fitRates(C, r, ids) {
    const sc = Math.max(...r.map(Math.abs)) || 1;
    return ids.map(id => {
      const M = MOD[id];
      const k0 = N.mean(r.map((ri, i) => ri / Math.max(M.f(M.p.length > 1 ? [1, id === 'pn' ? 1 : id === 'p1r' ? 0 : 1e-9] : [1], C[i]), 1e-300)));
      let p0 = [Math.abs(k0) || 1];
      if (id === 'pn') p0 = [Math.abs(N.mean(r.map((ri, i) => ri / C[i]))) || 1, 1];
      if (id === 'lh' || id === 'lh2') p0 = [Math.abs(N.mean(r.map((ri, i) => ri / C[i]))) * 1.5 || 1, 1 / N.mean(C)];
      if (id === 'p1r') p0 = [Math.abs(N.mean(r.map((ri, i) => ri / C[i]))) || 1, Math.min(...C) * 0.2];
      const fit = N.lm(p => C.map((c, i) => (M.f(p, c) - r[i]) / sc), p0, { pos: M.pos, maxIt: 200 });
      const pred = C.map(c => M.f(fit.p, c));
      return { id, name: M.name, p: fit.p, se: fit.se, ssr: fit.ssr * sc * sc, aic: aicc(C.length, fit.ssr, M.p.length), r2: N.r2(r, pred), pred };
    });
  }
  /* fit candidates to integral data via a solver pred(p) → predicted y */
  function fitInt(ids, y, predFor, p0For, sc) {
    return ids.map(id => {
      const M = MOD[id];
      const fit = N.lm(p => { const yp = predFor(id, p); return yp.map((q, i) => (isFinite(q) ? (q - y[i]) / sc : 10)); }, p0For(id), { pos: M.pos, maxIt: 120 });
      const pred = predFor(id, fit.p);
      return { id, name: M.name, p: fit.p, se: fit.se, ssr: fit.ssr * sc * sc, aic: aicc(y.length, fit.ssr, M.p.length), r2: N.r2(y, pred), pred };
    });
  }
  const parStr = f => MOD[f.id].p.map((nm, i) => nm + ' = ' + fmt(f.p[i], 4) + (isFinite(f.se[i]) ? ' ± ' + fmt(f.se[i], 2) : '')).join(', ');

  /* ===================================================================== */
  /* C-1803  rate equation from catalytic reactor data                      */
  /* ===================================================================== */
  ST.tool({
    id: 'catdata', tag: 'C-1803', ch: 18, sec: 'Ch. 18.6', title: 'Rate equation from differential, integral, mixed, recycle and batch catalytic reactors',
    lede: 'Turn catalytic reactor runs into a rate equation, with the analysis each experimental reactor calls for: rates straight from a differential reactor or a basket mixed reactor; the integral and the differential method for an integral (plug) reactor; single runs with recycle; and the closed-loop batch G/S reactor read by total pressure. Every candidate rate form is fitted with its uncertainty and ranked.',
    keys: 'experimental catalytic reactor differential reactor integral reactor plug flow mixed flow basket carberry berty recycle reactor batch closed loop rate equation fitting example 18.2 18.3 problem 18.4 18.5 18.6 18.7 18.8 18.15 18.17 18.18 18.19 weight time first order',
    units: { C: 'mol/L', t: 'h' },
    inputs: [
      { group: 'Reactor and feed' },
      { k: 'rt', t: 'select', label: 'Experimental reactor', opts: [['diff', 'differential reactor: C_in and C_out per run'], ['int', 'integral (plug flow) reactor: W and C_out'], ['mfr', 'mixed flow (basket, Berty, high recycle): flow and C_out'], ['rec', 'single runs with any recycle ratio'], ['batch', 'batch of gas over a batch of catalyst']], v: 'diff' },
      { k: 'CA0', label: 'C_A0 (feed basis: pure feed, unconverted)', unit: '[C]', v: 0.1, min: 1e-12, max: 1e6, if: v => v.rt !== 'batch' && v.rt !== 'rec' },
      { k: 'eps', label: 'Expansion factor ε_A', v: 3, min: -0.99, max: 20, if: v => v.rt !== 'batch' },
      { k: 'W', label: 'Catalyst W', unit: 'kg', v: 0.01, min: 1e-12, max: 1e9, if: v => v.rt === 'diff' || v.rt === 'mfr' || v.rt === 'batch' },
      { k: 'fb', t: 'seg', label: 'Feed given as', opts: [['v0', 'volumetric v₀'], ['FA0', 'molar F_A0']], v: 'v0', if: v => v.rt === 'diff' || v.rt === 'int' || v.rt === 'mfr' },
      { k: 'v0', label: 'Feed rate v₀', unit: v => KIT.volUnit(v._cu) + '/[t]', v: 20, min: 1e-12, max: 1e12, if: v => (v.rt === 'diff' || v.rt === 'int') && v.fb === 'v0' },
      { k: 'FA0', label: 'Feed rate F_A0', unit: 'mol/[t]', v: 2, min: 1e-12, max: 1e12, if: v => (v.rt === 'diff' || v.rt === 'int') && v.fb === 'FA0' },
      { k: 'diffT', t: 'table', label: 'Runs: C_A,in and C_A,out', cols: ['C_A,in', 'C_A,out'], units: ['[C]', '[C]'], v: [[0.1, 0.084], [0.08, 0.07], [0.06, 0.055], [0.04, 0.038]], if: v => v.rt === 'diff' },
      { k: 'intT', t: 'table', label: 'Runs: W and C_A,out', cols: ['W, kg', 'C_A,out'], units: ['', '[C]'], v: [[0.02, 0.074], [0.04, 0.06], [0.08, 0.044], [0.16, 0.029]], if: v => v.rt === 'int' },
      { k: 'ratio', t: 'check', label: 'Second column is C_A,out/C_A0 (or p_A,out/p_A,in)', v: false, if: v => v.rt === 'mfr' },
      { k: 'mfrT', t: 'table', label: v => 'Runs: ' + (v.fb === 'FA0' ? 'F_A0' : 'v₀') + ' and the leaving C_A', cols: ['flow', 'C_A,out'], v: [[0.14, 8], [0.42, 6], [1.67, 4], [2.5, 2], [1.25, 1]], if: v => v.rt === 'mfr' },
      { k: 'recT', t: 'table', label: 'Runs: recycle ratio R (0 = none, ≥1000 = very large), W, v₀, C_A0, C_A,out', cols: ['R', 'W', 'v₀', 'C_A0', 'C_A,out'], v: [[1000, 3, 1, 2, 0.5]], if: v => v.rt === 'rec' },
      { k: 'nrec', label: 'Order n to evaluate (single runs)', v: 1, min: 0, max: 4, step: 0.5, if: v => v.rt === 'rec' },
      { k: 'Vg', label: 'Gas volume of the loop V', unit: v => KIT.volUnit(v._cu), v: 0.1, min: 1e-12, max: 1e9, if: v => v.rt === 'batch' },
      { k: 'yk', t: 'select', label: 'Measured', opts: [['C', 'C_A directly'], ['pA', 'partial pressure of A'], ['pi', 'total pressure π']], v: 'pi', if: v => v.rt === 'batch' },
      { k: 'TK', label: 'Temperature', unit: 'K', v: 609, min: 1, max: 5000, if: v => v.rt === 'batch' && v.yk !== 'C' },
      { k: 'pu', t: 'select', label: 'Pressure unit', opts: [['atm', 'atm'], ['mmHg', 'mm Hg'], ['kPa', 'kPa'], ['bar', 'bar']], v: 'atm', if: v => v.rt === 'batch' && v.yk !== 'C' },
      { k: 'nu', label: 'Moles of gas formed per mole of A reacted (2A → R: ½)', v: 0.5, min: 0, max: 20, if: v => v.rt === 'batch' && v.yk === 'pi' },
      { k: 'pI', label: 'Partial pressure of inerts in the loop', v: 0, min: 0, max: 1e6, if: v => v.rt === 'batch' && v.yk === 'pi' },
      { k: 'batT', t: 'table', label: 'Readings: t and the measured quantity', cols: ['t', 'reading'], units: ['[t]', ''], v: [[0, 1], [4, 0.75], [8, 0.67], [16, 0.6], [36, 0.55]], if: v => v.rt === 'batch' },
      { group: 'Rate forms to try' },
      { k: 'cands', t: 'text', label: 'Candidates (p0 p1 p2 pn lh lh2 p1r)', v: 'p1 p2 pn lh', mono: true }
    ],
    presets: [
      { name: 'Example 18.2 — differential reactor', v: { rt: 'diff', CA0: 0.1, eps: 3, W: 0.01, fb: 'v0', v0: 20, diffT: [[0.1, 0.084], [0.08, 0.07], [0.06, 0.055], [0.04, 0.038]], cands: 'p1 p2 pn lh', _cu: 'mol/L', _tu: 'h' }, book: { k1: { v: 96, tol: 0.03 } }, note: 'A → 4R at 3.2 atm and 117 °C, W = 0.01 kg, 20 L/hr of feed, C_A0 = 0.1 mol/L, ε_A = 3. The book draws −r′_A against C_A by eye: k′ = 96 L/hr·kg cat.' },
      { name: 'Example 18.3a — integral reactor, integral method', v: { rt: 'int', CA0: 0.1, eps: 3, fb: 'v0', v0: 20, intT: [[0.02, 0.074], [0.04, 0.06], [0.08, 0.044], [0.16, 0.029]], cands: 'p1 p2 pn lh', _cu: 'mol/L', _tu: 'h' }, book: { kInt: { v: 95, tol: 0.03 }, kDiff: { v: 93, tol: 0.03 } }, note: 'Same reaction, plug flow with 0.02 to 0.16 kg of catalyst. The integral test of Eq. (i) gives k′ = 95; differentiating the X_A against W/F_A0 curve (Example 18.3b) gives 93 L/hr·kg.' },
      { name: 'Problem 18.9 data — basket reactor', v: { rt: 'mfr', CA0: 10, eps: 0, W: 0.1, fb: 'FA0', ratio: false, mfrT: [[0.14, 8], [0.42, 6], [1.67, 4], [2.5, 2], [1.25, 1]], cands: 'p1 p2 pn lh lh2', _cu: 'mol/m³', _tu: 'min' }, note: 'Carberry basket reactor, 100 g of catalyst, C_A0 = 10 mol/m³, F_A0 from 0.14 to 2.5 mol/min. The rate rises as C_A falls: strong adsorption of A.' },
      { name: 'Problem 18.18 data — high recycle', v: { rt: 'mfr', CA0: 10, eps: 0, W: 1, fb: 'v0', ratio: false, mfrT: [[0.005, 1], [0.02, 2], [0.065, 3], [0.133, 6], [0.54, 9]], cands: 'p1 p2 pn lh lh2', _cu: 'mol/m³', _tu: 'h' }, note: 'W = 1 kg, very large recycle, C_A0 = 10 mol/m³; the book’s v₀ of 5 to 540 L/hr is entered in m³/hr to match the unit bar.' },
      { name: 'Problem 18.15 — basket reactor, A → 3R at 8 atm', v: { rt: 'mfr', CA0: 0.1002, eps: 2, W: 0.001, fb: 'v0', ratio: true, mfrT: [[100, 0.8], [22, 0.5], [4, 0.2], [1, 0.1], [0.6, 0.05]], cands: 'p1 p2 pn lh', _cu: 'mol/L', _tu: 'h' }, note: 'Pure A at 8 atm and 700 °C (C_A0 = 0.100 mol/L), 1 g of 3 mm catalyst, feed 0.6 to 100 L/hr; the exit partial-pressure ratio is given.' },
      { name: 'Problem 18.4 — very large recycle, n = ½', v: { rt: 'rec', eps: 0, nrec: 0.5, recT: [[1000, 3, 1, 2, 0.5]], _cu: 'mol/m³', _tu: 'h' }, note: 'A (2 mol/m³) fed at 1 m³/hr through 3 kg of catalyst with a recycle loop; 0.5 mol/m³ leaves. A → R, half order.' },
      { name: 'Problem 18.5 — very large recycle, A → 3R, half inerts', v: { rt: 'rec', eps: 1, nrec: 1, recT: [[1000, 3, 1, 2, 0.5]], _cu: 'mol/m³', _tu: 'h' }, note: 'ε_A = (3 − 1)(0.5) = 1; first order.' },
      { name: 'Problem 18.6 — no recycle, A → 3R, 25% A', v: { rt: 'rec', eps: 0.5, nrec: 2, recT: [[0, 3, 1, 2, 0.5]], _cu: 'mol/m³', _tu: 'h' }, note: 'ε_A = (3 − 1)(0.25) = 0.5; second order; plug flow.' },
      { name: 'Problem 18.17a — second order in a recycle reactor', v: { rt: 'rec', eps: 0, nrec: 2, recT: [[1000, 3, 1, 2, 0.5]], _cu: 'mol/L', _tu: 'h' }, note: '3 g of catalyst (enter W in the unit you want k in: here grams), 1 L/hr of feed at 2 mol/L leaving at 0.5 mol/L: k′ = 2 L²/mol·g·hr.' },
      { name: 'Problem 18.19 — closed-loop batch G/S reactor', v: { rt: 'batch', W: 2, Vg: 0.1, yk: 'pi', TK: 609, pu: 'atm', nu: 0.5, pI: 0, batT: [[0, 1], [4, 0.75], [8, 0.67], [16, 0.6], [36, 0.55]], cands: 'p1 p2 pn lh p1r', _cu: 'mol/L', _tu: 'min' }, note: 'Pure A at 609 K, 2A → R, 100 cm³ of gas and 2 g of nonporous catalyst (enter W in g to get rates per gram). π falls from 1 to 0.55 atm: p_A = 2π − π₀.' }
    ],
    eqs: [
      ['\\frac{W}{F_{A0}} = \\int_{X_{A,in}}^{X_{A,out}}\\frac{dX_A}{-r_A^\\prime} \\approx \\frac{X_{A,out}-X_{A,in}}{(-r_A^\\prime)_{avg}}', 'Eq. 18.52', 'Differential reactor.'],
      ['\\frac{W}{F_{A0}} = \\int_0^{X_A}\\frac{dX_A}{-r_A^\\prime},\\qquad -r_A^\\prime = \\frac{dX_A}{d(W/F_{A0})}', 'Eqs. 18.53–18.54', 'Integral reactor: integral and differential analysis.'],
      ['k^\\prime\\frac{C_{A0}W}{F_{A0}} = (1+\\varepsilon_A)\\ln\\frac{1}{1-X_A} - \\varepsilon_A X_A', 'Eq. 18.44'],
      ['-r_{A,out}^\\prime = \\frac{F_{A0}X_{A,out}}{W}', 'Eq. 18.55', 'Mixed flow: each run gives a rate.'],
      ['\\frac{W}{F_{A0}} = (R+1)\\int_{\\frac{R}{R+1}X_{Af}}^{X_{Af}}\\frac{dX_A}{-r_A^\\prime}', 'Eq. 6.21', 'Recycle reactor.'],
      ['\\frac{t}{C_{A0}} = \\frac{V}{W}\\int\\frac{dX_A}{-r_A^\\prime}', 'Eq. 18.56', 'Batch of gas, batch of catalyst.'],
      ['\\frac{C_A}{C_{A0}} = \\frac{1-X_A}{1+\\varepsilon_A X_A}', 'Eq. 3.73']
    ],
    method: ['Rates are formed exactly as in the book. Candidate rate forms are fitted by Levenberg–Marquardt: directly to the rates for differential and mixed reactors, and to the measured conversions (plug flow) or concentrations (batch) by integrating each candidate, so no rate is ever estimated by hand differentiation. Standard errors come from the Jacobian; forms are ranked by AICc.', 'The integral reactor is also analysed both of the book’s ways: the first-order linearization of Eq. 18.44 through the origin, and the differential method with a cubic through the origin for X_A against W/F_A0, whose slopes are the rates.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const cu = v._cu || 'mol/L', tu = v._tu || 'h', V = KIT.volUnit(cu);
      const ids = String(v.cands || '').split(/[\s,;]+/).filter(s => MOD[s]);
      if (!ids.length) throw new Error('List at least one candidate rate form.');
      const eps = +v.eps;
      const rateU = 'mol/(kg·' + tu + ')';
      const rankTable = (fits, title) => {
        fits.sort((a, b) => a.aic - b.aic);
        tables.push({ title: title || 'Rate forms, best first', cols: ['rate form', 'constants', 'R²', 'ΔAICc'], rows: fits.map(f => [f.name, parStr(f), f.r2, f.aic - fits[0].aic]) });
        kpi.push({ k: 'best', l: 'Best rate form (AICc)', v: fits[0].name, s: parStr(fits[0]) });
        return fits;
      };
      const ratePlot = (C, r, fits, title, extra) => {
        const cmax = Math.max(...C) * 1.15, cg = N.linspace(0, cmax, 160);
        const ser = [{ name: 'runs', x: C, y: r, type: 'pts', color: 1 }].concat(fits.slice(0, 4).map((f, i) => ({ name: f.name, x: cg, y: cg.map(c => MOD[f.id].f(f.p, c)), color: i + 2, dash: i ? '5 4' : undefined })));
        if (extra) ser.push(extra);
        plots.push({ title: title || 'Rate against concentration', x: { label: 'C_A', unit: cu, min: 0 }, y: { label: '−r′_A', unit: rateU, min: 0 }, series: ser });
      };
      if (v.rt === 'diff' || v.rt === 'mfr') {
        const CA0 = +v.CA0, W = +v.W;
        let C = [], r = [], rows = [];
        if (v.rt === 'diff') {
          const FA0 = v.fb === 'FA0' ? +v.FA0 : +v.v0 * CA0;
          (v.diffT || []).filter(q => q && isFinite(+q[0]) && isFinite(+q[1]) && q[0] !== '' && q[1] !== '').forEach(q => {
            const ci = +q[0], co = +q[1], xi = Xof(CA0, eps)(ci), xo = Xof(CA0, eps)(co), rr = (xo - xi) / (W / FA0);
            C.push((ci + co) / 2); r.push(rr);
            rows.push([ci, co, (ci + co) / 2, xi, xo, xo - xi, rr, (ci - co) / ((ci + co) / 2)]);
          });
          if (C.length < 1) throw new Error('Enter at least one run.');
          tables.push({ title: 'Runs (Table E18.2)', cols: ['C_A,in', 'C_A,out', 'C_A,avg', 'X_A,in', 'X_A,out', 'ΔX_A', '−r′_A = ΔX_A/(W/F_A0)', 'spread / mean'], rows });
          kpi.push({ k: 'FA0', l: 'F_A0', v: FA0, u: 'mol/' + tu });
          const spread = Math.max(...rows.map(q => q[7]));
          if (spread > 0.25) warn.push('Concentration changes by ' + fmt(spread * 100, 3) + '% in one run: too much for a differential reactor; treat it as an integral reactor.');
        } else {
          (v.mfrT || []).filter(q => q && isFinite(+q[0]) && isFinite(+q[1]) && q[0] !== '' && q[1] !== '').forEach(q => {
            const co = v.ratio ? +q[1] * CA0 : +q[1], xo = Xof(CA0, eps)(co);
            const FA0 = v.fb === 'FA0' ? +q[0] : +q[0] * CA0, rr = FA0 * xo / W;
            C.push(co); r.push(rr); rows.push([+q[0], co, xo, W / FA0, rr]);
          });
          if (C.length < 1) throw new Error('Enter at least one run.');
          tables.push({ title: 'Runs: each gives the rate at the leaving composition (Eq. 18.55)', cols: [v.fb === 'FA0' ? 'F_A0' : 'v₀', 'C_A,out', 'X_A,out', 'W/F_A0', '−r′_A'], rows });
        }
        const lin = N.linreg(C, r, { origin: true });
        kpi.push({ k: 'k1', l: 'First order through the origin: k′', v: lin.b, u: CAT.kpU(1, cu, tu), main: true, s: '± ' + fmt(lin.sb, 2) + ' (R² = ' + fmt(lin.r2, 4) + ')' });
        if (C.length >= 2) { const fits = rankTable(fitRates(C, r, ids.filter(id => MOD[id].p.length < C.length || C.length > 2))); ratePlot(C, r, fits, v.rt === 'diff' ? 'Rate against mean concentration (Fig. E18.2)' : 'Rate against the leaving concentration'); }
        else ratePlot(C, r, [], 'Rate against concentration');
        return { kpi, warn, info, plots, tables };
      }
      if (v.rt === 'int') {
        const CA0 = +v.CA0, FA0 = v.fb === 'FA0' ? +v.FA0 : +v.v0 * CA0;
        const pts = (v.intT || []).filter(q => q && isFinite(+q[0]) && isFinite(+q[1]) && q[0] !== '' && q[1] !== '').map(q => ({ W: +q[0], C: +q[1] })).sort((a, b) => a.W - b.W);
        if (pts.length < 2) throw new Error('Enter at least two runs.');
        const tau = pts.map(q => q.W / FA0), X = pts.map(q => Xof(CA0, eps)(q.C));
        /* integral test for first order (Eq. 18.44) */
        const yv = X.map(x => (1 + eps) * Math.log(1 / (1 - x)) - eps * x), xv = tau.map(t => CA0 * t);
        const L1 = N.linreg(xv, yv, { origin: true });
        kpi.push({ k: 'kInt', l: 'Integral method, first order: k′', v: L1.b, u: CAT.kpU(1, cu, tu), main: true, s: 'R² = ' + fmt(L1.r2, 4) });
        /* differential method: cubic through the origin */
        const A = tau.map(t => [t, t * t, t * t * t]), cf = N.lstsq(A, X).x;
        const dX = t => cf[0] + 2 * cf[1] * t + 3 * cf[2] * t * t;
        const tD = [0].concat(tau), CD = [CA0].concat(pts.map(q => q.C)), rD = tD.map(dX);
        const L2 = N.linreg(CD, rD, { origin: true });
        kpi.push({ k: 'kDiff', l: 'Differential method, first order: k′', v: L2.b, u: CAT.kpU(1, cu, tu), s: 'slopes of the fitted X_A(W/F_A0) curve' });
        tables.push({ title: 'Integral and differential analysis (Tables E18.3a, b)', cols: ['W', 'W/F_A0', 'C_A,out', 'X_A', '(1+ε)ln 1/(1−X) − εX', 'C_A0 W/F_A0', 'slope dX/d(W/F_A0)'], rows: [[0, 0, CA0, 0, 0, 0, dX(0)]].concat(pts.map((q, i) => [q.W, tau[i], q.C, X[i], yv[i], xv[i], dX(tau[i])])) });
        /* model fits on X by integration */
        const tmax = Math.max(...tau);
        const Cx = Cof(CA0, eps);
        const predFor = (id, p) => {
          const M = MOD[id];
          const sol = N.ode((t, y) => [Math.max(M.f(p, Cx(Math.min(y[0], 0.999999999))), 0)], 0, [0], tmax, { rtol: 1e-8, atol: 1e-12 });
          return tau.map(t => sol.at(t)[0]);
        };
        const p0For = id => (id === 'pn' ? [L1.b, 1] : id === 'lh' || id === 'lh2' ? [L1.b * 1.5, 1 / CA0] : id === 'p1r' ? [L1.b, CA0 * 0.05] : id === 'p2' ? [L1.b / CA0] : id === 'p0' ? [L1.b * CA0] : [L1.b]);
        const fits = rankTable(fitInt(ids, X, predFor, p0For, 1), 'Rate forms fitted to the conversions, best first');
        const tg = N.linspace(0, tmax * 1.1, 120);
        plots.push({ title: 'Conversion against weight-time (Fig. E18.3b)', x: { label: 'W/F_A0', unit: 'kg·' + tu + '/mol', min: 0 }, y: { label: 'X_A', min: 0 }, series: [{ name: 'runs', x: [0].concat(tau), y: [0].concat(X), type: 'pts', color: 1 }, { name: 'cubic through the origin', x: tg, y: tg.map(t => cf[0] * t + cf[1] * t * t + cf[2] * t * t * t), color: 'var(--muted)', dash: '3 3' }].concat(fits.slice(0, 3).map((f, i) => { const yy = (() => { const M = MOD[f.id]; const sol = N.ode((t, y) => [Math.max(M.f(f.p, Cx(Math.min(y[0], 0.999999999))), 0)], 0, [0], tmax * 1.1, { rtol: 1e-8, atol: 1e-12 }); return tg.map(t => sol.at(t)[0]); })(); return { name: f.name, x: tg, y: yy, color: i + 2 }; })) });
        const xl = N.linspace(0, Math.max(...xv) * 1.1, 50);
        plots.push({ title: 'Integral test for first order (Fig. E18.3a)', x: { label: 'C_A0 W/F_A0', unit: 'kg·' + tu + '/' + V, min: 0 }, y: { label: '(1+ε_A) ln 1/(1−X_A) − ε_A X_A', min: 0 }, series: [{ name: 'runs', x: xv, y: yv, type: 'pts', color: 1 }, { name: 'slope k′ = ' + fmt(L1.b, 3), x: xl, y: xl.map(x => L1.b * x), color: 2 }] });
        ratePlot(CD, rD, fits, 'Rates from the slopes (Fig. E18.3c)', { name: 'k′ = ' + fmt(L2.b, 3) + ' through the origin', x: [0, CA0 * 1.1], y: [0, L2.b * CA0 * 1.1], color: 6, dash: '2 3' });
        return { kpi, warn, info, plots, tables };
      }
      if (v.rt === 'rec') {
        const n = +v.nrec;
        const rows = [], ks = [];
        (v.recT || []).filter(q => q && q.slice(0, 5).every(z => z !== '' && isFinite(+z))).forEach(q => {
          const R = +q[0], W = +q[1], v0 = +q[2], CA0 = +q[3], C = +q[4];
          const Cx = Cof(CA0, eps), Xf = Xof(CA0, eps)(C), FA0 = v0 * CA0;
          let k, pat;
          if (R >= 1000) { k = FA0 * Xf / (W * Math.pow(C, n)); pat = 'mixed (R → ∞)'; }
          else {
            const Xi = R * Xf / (R + 1);
            const I = N.quad(x => 1 / Math.pow(Cx(x), n), Xi, Xf, { rtol: 1e-11 });
            k = FA0 * (R + 1) * I / W; pat = R === 0 ? 'plug' : 'recycle R = ' + fmt(R, 3);
          }
          ks.push(k); rows.push([pat, W, v0, CA0, C, Xf, k]);
        });
        if (!ks.length) throw new Error('Enter at least one run.');
        tables.push({ title: 'Rate constant from each run', cols: ['flow pattern', 'W', 'v₀', 'C_A0', 'C_A,out', 'X_A', 'k′'], rows });
        const km = N.mean(ks);
        kpi.push({ k: 'k', l: 'Rate constant k′ for order ' + fmt(n, 3), v: km, u: CAT.kpU(n, cu, tu), main: true, s: '−r′_A = k′ C_A' + (n !== 1 ? '^{' + fmt(n, 3) + '}' : '') + (ks.length > 1 ? '; spread ' + fmt((Math.max(...ks) - Math.min(...ks)) / km * 100, 3) + '%' : '') });
        if (ks.length > 1 && (Math.max(...ks) - Math.min(...ks)) / km > 0.1) warn.push('The runs do not agree on k′ for this order: try another order.');
        const r0 = rows[0], CA0 = r0[3];
        const cg = N.linspace(0, CA0, 100);
        plots.push({ title: 'The rate equation found', x: { label: 'C_A', unit: cu, min: 0 }, y: { label: '−r′_A', unit: rateU, min: 0 }, series: [{ name: '−r′ = k′Cⁿ', x: cg, y: cg.map(c => km * Math.pow(c, n)), color: 2 }], marks: rows.map(q => ({ t: 'pt', x: q[4], y: km * Math.pow(q[4], n), color: 8, label: 'exit' })) });
        return { kpi, warn, info, plots, tables };
      }
      /* batch of gas over a batch of catalyst */
      const W = +v.W, Vg = +v.Vg;
      const raw = (v.batT || []).filter(q => q && isFinite(+q[0]) && isFinite(+q[1]) && q[0] !== '' && q[1] !== '').map(q => [+q[0], +q[1]]).sort((a, b) => a[0] - b[0]);
      if (raw.length < 3) throw new Error('Enter at least three readings.');
      const pf = { atm: 1, mmHg: 1 / 760, kPa: 1 / 101.325, bar: 1 / 1.01325 }[v.pu] || 1;
      const Rl = 0.082057366; /* L·atm/(mol·K) */
      const volF = { L: 1, 'm³': 1000, 'cm³': 0.001 }[V] || 1;
      const toC = pA => pA * pf / (Rl * +v.TK) / volF * (cu === 'mmol/L' ? 1000 : 1);
      let t = raw.map(q => q[0]), C;
      if (v.yk === 'C') C = raw.map(q => q[1]);
      else if (v.yk === 'pA') C = raw.map(q => toC(q[1]));
      else {
        const nu = +v.nu, pI = +v.pI, pA0 = raw[0][1] - pI;
        if (Math.abs(nu - 1) < 1e-9) throw new Error('With ν = 1 the total pressure does not change; measure p_A instead.');
        const pA = raw.map(q => (nu * pA0 - (q[1] - pI)) / (nu - 1));
        C = pA.map(toC);
        tables.push({ title: 'Partial pressure of A from the total pressure', cols: ['t', 'π', 'p_A', 'C_A'], rows: raw.map((q, i) => [q[0], q[1], pA[i], C[i]]) });
      }
      const C0 = C[0], tmax = Math.max(...t);
      const predFor = (id, p) => {
        const M = MOD[id];
        const sol = N.ode((tt, y) => [-(W / Vg) * M.f(p, Math.max(y[0], 0))], 0, [C0], tmax, { rtol: 1e-8, atol: 1e-14 * Math.max(C0, 1e-300) });
        return t.map(tt => sol.at(tt)[0]);
      };
      const kGuess = Vg / W * Math.log(C0 / Math.max(C[C.length - 1], C0 * 1e-6)) / tmax;
      const p0For = id => (id === 'pn' ? [kGuess, 1] : id === 'lh' || id === 'lh2' ? [kGuess * 1.5, 1 / C0] : id === 'p1r' ? [kGuess * 1.3, Math.min(...C) * 0.8] : id === 'p2' ? [kGuess / C0] : id === 'p0' ? [kGuess * C0] : [kGuess]);
      const fits = rankTable(fitInt(ids, C, predFor, p0For, C0), 'Rate forms fitted to C_A(t), best first');
      const best = fits[0];
      kpi.unshift({ k: 'kbest', l: 'Best fit: ' + best.name, v: best.p[0], u: CAT.kpU(MOD[best.id].order(best.p), cu, tu, 'kg'), main: true, s: parStr(best) });
      const tg = N.linspace(0, tmax * 1.05, 150);
      plots.push({ title: 'Concentration of A in the loop', x: { label: 't', unit: tu, min: 0 }, y: { label: 'C_A', unit: cu, min: 0 }, series: [{ name: 'readings', x: t, y: C, type: 'pts', color: 1 }].concat(fits.slice(0, 3).map((f, i) => { const M = MOD[f.id]; const sol = N.ode((tt, y) => [-(W / Vg) * M.f(f.p, Math.max(y[0], 0))], 0, [C0], tmax * 1.05, { rtol: 1e-8, atol: 1e-14 * C0 }); return { name: f.name, x: tg, y: tg.map(tt => sol.at(tt)[0]), color: i + 2, dash: i ? '5 4' : undefined }; })) });
      const dC = N.diffData(t, C, 'local', { half: 1 });
      ratePlot(C, dC.map(d => -d * Vg / W), fits, 'Rates from local slopes −(V/W) dC_A/dt');
      info.push('Rates are per unit of the catalyst amount entered as W (kg or g).');
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* C-1804  catalytic reactor size: plug, mixed, recycle, batch, f(z)      */
  /* ===================================================================== */
  const lawOrder = v => (v.law === 'nth' ? +v.n : v.law === 'ab' ? 2 : 1);
  function rateModel(v) {
    const CA0 = +v.CA0, eps = +v.eps, Cx = Cof(CA0, eps);
    const rhoS = +v.rhoS, rhoB = +v.rhoB;
    /* constant per kg from the chosen basis */
    const kmass = v.kb === 'pellet' ? +v.k / rhoS : v.kb === 'bed' ? +v.k / rhoB : +v.k;
    let free, rC;
    if (v.law === 'nth') rC = C => (C > 0 ? kmass * Math.pow(C, +v.n) : 0);
    else if (v.law === 'lh') rC = C => (C > 0 ? kmass * C / (1 + +v.K * C) : 0);
    else if (v.law === 'rev') rC = C => kmass * (C - (CA0 - C) / +v.Keq);
    else if (v.law === 'ab') {
      /* A + bB → products, −r′ = k′ C_A C_B; C_B follows from stoichiometry */
      const M = +v.CB0 / CA0, b = +v.bB;
      free = X => { const CA = Cx(X), CB = CA0 * (M - b * X) / (1 + eps * X); return CB > 0 && CA > 0 ? kmass * CA * CB : 0; };
    } else if (v.law === 'data') {
      const rows = (v.rateT || []).filter(q => q && +q[0] > 0 && +q[1] > 0).map(q => [+q[0], +q[1]]).sort((a, b) => a[0] - b[0]);
      if (rows.length < 2) throw new Error('Enter at least two (C_A, −r′_A) points.');
      const lx = rows.map(q => Math.log(q[0])), ly = rows.map(q => Math.log(q[1]));
      const pc = N.pchip(lx, ly);
      const s0 = (ly[1] - ly[0]) / (lx[1] - lx[0]), s1 = (ly[ly.length - 1] - ly[ly.length - 2]) / (lx[lx.length - 1] - lx[lx.length - 2]);
      rC = C => {
        if (!(C > 0)) return 0;
        const l = Math.log(C);
        if (l < lx[0]) return Math.exp(ly[0] + s0 * (l - lx[0]));
        if (l > lx[lx.length - 1]) return Math.exp(ly[ly.length - 1] + s1 * (l - lx[lx.length - 1]));
        return Math.exp(pc(l));
      };
      rC.dataRange = [rows[0][0], rows[rows.length - 1][0]];
    }
    if (!free) free = X => rC(Cx(X));
    return { CA0, eps, Cx, free, rC, kmass };
  }
  /* effectiveness factor along the reactor: exact for nth order, generalized modulus otherwise */
  function etaModel(v, RM) {
    if (!v.pore) return () => 1;
    const sh = v.shape, L = CAT.Lof(sh, +v.dp / 1000), De = +v.De, rhoS = +v.rhoS;
    if (v.law === 'nth') {
      const n = +v.n, tr = CAT.trace({ shape: sh, kind: 'nth', n });
      const k3 = RM.kmass * rhoS;
      return (X, Lc) => { const C = RM.Cx(X); if (!(C > 0)) return 1; const MT = (Lc || L) * Math.sqrt(k3 * Math.pow(C, n - 1) / De); return PEL.etaAt(tr, MT)[0].eta; };
    }
    /* generalized modulus with the local rate form against concentration (Eq. 18.26) on the exact first-order curve */
    const rOfC = v.law === 'ab' ? (C, X) => RM.free(X) * C / Math.max(RM.Cx(X), 1e-300) : (C) => RM.rC(C);
    const Ce = v.law === 'rev' ? (CA0X => CA0X) : null;
    return (X, Lc) => {
      const Cs = RM.Cx(X);
      if (!(Cs > 0)) return 1;
      const rs = v.law === 'ab' ? RM.free(X) : RM.rC(Cs);
      if (!(rs > 0)) return 1;
      let lo = 0;
      if (v.law === 'rev') { const Keq = +v.Keq, Ctot = RM.CA0 * (1 - X) / (1 + RM.eps * X) + RM.CA0 * X / (1 + RM.eps * X); lo = Ctot / (1 + Keq); }
      const integ = N.quad(c => Math.max(v.law === 'ab' ? rs * (c / Cs) : RM.rC(c), 0), lo, Cs, { rtol: 1e-7 });
      const Mg = rs * rhoS * (Lc || L) / Math.sqrt(2 * De * rhoS * integ);
      return PEL.eta1(sh, Mg);
    };
  }
  const sizeReactor = (v, RM, etaF, pattern, X, R) => {
    const FA0 = v.fb === 'FA0' ? +v.FA0 : +v.v0 * RM.CA0;
    const rate = x => RM.free(x) * etaF(x);
    if (pattern === 'mixed') { const r = rate(X); return { W: FA0 * X / r, rOut: r }; }
    if (pattern === 'rec') { const Xi = R * X / (R + 1); return { W: FA0 * (R + 1) * N.quad(x => 1 / rate(x), Xi, X, { rtol: 1e-8 }) }; }
    return { W: FA0 * N.quad(x => 1 / rate(x), 0, X, { rtol: 1e-9 }) };
  };

  ST.tool({
    id: 'catdesign', tag: 'C-1804', ch: 18, sec: 'Ch. 18.5', title: 'Catalytic reactor size: plug, mixed, recycle, batch, with pore diffusion',
    lede: 'How much catalyst — and how much bed — for a given duty, in plug flow, mixed flow, with any recycle ratio, in a batch G/S reactor, or in a vertical reactor whose solid fraction varies with height. The rate can be a law or a table of measured rates, and pore diffusion can be switched on with the effectiveness factor recomputed at every point of the reactor: exactly for nth-order kinetics and through the generalized modulus for the rest. Compare flow patterns and pellet sizes side by side.',
    keys: 'catalytic reactor design packed bed weight of catalyst plug flow mixed flow recycle batch fluidized bed solid fraction f(z) effectiveness factor along the reactor pellet size comparison bed volume bulk density example 18.4 18.5 18.6 problem 18.3 18.9 18.10 18.11 18.12 18.13 18.17 18.31 18.32 18.33 18.34 18.35 18.36',
    units: { C: 'mol/L', t: 'h' },
    inputs: [
      { group: 'Feed and duty' },
      { k: 'fb', t: 'seg', label: 'Feed given as', opts: [['FA0', 'molar F_A0'], ['v0', 'volumetric v₀']], v: 'FA0' },
      { k: 'FA0', label: 'F_A0', unit: 'mol/[t]', v: 2000, min: 1e-12, max: 1e15, scale: 'log', if: v => v.fb === 'FA0' },
      { k: 'v0', label: 'v₀', unit: v => KIT.volUnit(v._cu) + '/[t]', v: 20000, min: 1e-12, max: 1e15, scale: 'log', if: v => v.fb === 'v0' },
      { k: 'CA0', label: 'C_A0', unit: '[C]', v: 0.1, min: 1e-12, max: 1e7 },
      { k: 'eps', label: 'ε_A', v: 3, min: -0.99, max: 20 },
      { k: 'solve', t: 'seg', label: 'Find', opts: [['W', 'catalyst for a conversion'], ['X', 'conversion for a catalyst load']], v: 'W' },
      { k: 'X', label: 'Conversion X_A', v: 0.35, min: 0.0001, max: 0.9999, if: v => v.solve === 'W' },
      { k: 'Wgiven', label: 'Catalyst W', unit: 'kg', v: 100, min: 1e-12, max: 1e12, scale: 'log', if: v => v.solve === 'X' },
      { group: 'Rate, free of pore diffusion' },
      { k: 'law', t: 'select', label: 'Rate form', opts: [['nth', 'k Cⁿ'], ['lh', 'k C/(1 + K C)'], ['ab', 'k C_A C_B  (A + bB)'], ['rev', 'k (C_A − C_R/K), A ⇌ R'], ['data', 'measured rates (C_A, −r′_A)']], v: 'nth' },
      { k: 'kb', t: 'seg', label: 'Constant per', opts: [['mass', 'kg catalyst'], ['pellet', 'm³ of pellets'], ['bed', 'm³ of bed']], v: 'mass', if: v => v.law !== 'data' },
      { k: 'k', label: 'Rate constant', unit: v => (v.kb === 'mass' ? CAT.kpU(lawOrder(v), v._cu, v._tu) : CAT.kpU(lawOrder(v), v._cu, v._tu, 'm³ ' + (v.kb === 'pellet' ? 'pellet' : 'bed'))), v: 96, min: 1e-15, max: 1e15, scale: 'log', if: v => v.law !== 'data' },
      { k: 'n', label: 'Order n', v: 1, min: 0, max: 4, step: 0.5, if: v => v.law === 'nth' },
      { k: 'K', label: 'K', unit: v => KIT.volUnit(v._cu) + '/mol', v: 0.02, min: 0, max: 1e9, if: v => v.law === 'lh' },
      { k: 'CB0', label: 'C_B0', unit: '[C]', v: 10, min: 0, max: 1e7, if: v => v.law === 'ab' },
      { k: 'bB', label: 'b (mol B per mol A)', v: 1, min: 0, max: 10, if: v => v.law === 'ab' },
      { k: 'Keq', label: 'Equilibrium constant K = C_R/C_A', v: 4, min: 1e-6, max: 1e9, if: v => v.law === 'rev' },
      { k: 'rateT', t: 'table', label: 'Measured rates', cols: ['C_A', '−r′_A'], units: ['[C]', 'mol/(kg·[t])'], v: [[0.039, 3.4], [0.0575, 5.4], [0.075, 7.6], [0.092, 9.1]], if: v => v.law === 'data' },
      { group: 'Reactor' },
      { k: 'pat', t: 'select', label: 'Flow pattern', opts: [['plug', 'plug flow (packed bed)'], ['mixed', 'mixed flow (very large recycle, fluidized gas)'], ['rec', 'recycle with ratio R'], ['batch', 'batch of gas over the catalyst'], ['fz', 'vertical reactor, solids fraction f(z) (first order)']], v: 'plug' },
      { k: 'R', label: 'Recycle ratio R', v: 1, min: 0, max: 1e4, if: v => v.pat === 'rec' },
      { k: 'Vg', label: 'Gas volume V (batch)', unit: v => KIT.volUnit(v._cu), v: 1, min: 1e-12, max: 1e12, if: v => v.pat === 'batch' },
      { k: 'u0', label: 'Superficial gas velocity u₀', unit: 'm/[t]', v: 1, min: 1e-9, max: 1e9, if: v => v.pat === 'fz' },
      { k: 'fz', t: 'select', label: 'Solids profile', opts: [['const', 'uniform f over height H'], ['exp', 'f = f₀ e^{−a z} (Fig. 18.10 left)'], ['two', 'dense bottom zone + lean zone (Fig. 18.10 right)']], v: 'const', if: v => v.pat === 'fz' },
      { k: 'Hz', label: 'Height H', unit: 'm', v: 2, min: 1e-6, max: 1e4, if: v => v.pat === 'fz' },
      { k: 'f0', label: v => (v.fz === 'two' ? 'f in the dense zone' : 'f (at the bottom)'), v: 0.4, min: 1e-6, max: 0.7, if: v => v.pat === 'fz' },
      { k: 'az', label: 'Decay constant a', unit: '1/m', v: 1, min: 0, max: 1e3, if: v => v.pat === 'fz' && v.fz === 'exp' },
      { k: 'Hd', label: 'Height of the dense zone', unit: 'm', v: 0.5, min: 0, max: 1e4, if: v => v.pat === 'fz' && v.fz === 'two' },
      { k: 'fl', label: 'f in the lean zone', v: 0.02, min: 0, max: 0.7, if: v => v.pat === 'fz' && v.fz === 'two' },
      { group: 'Pellets and pore diffusion' },
      { k: 'pore', t: 'check', label: 'Account for pore diffusion', v: false },
      { k: 'shape', t: 'seg', label: 'Shape', opts: [['slab', 'flat plate'], ['cyl', 'cylinder'], ['sph', 'sphere']], v: 'sph', if: v => v.pore },
      { k: 'dp', label: v => 'Pellet ' + CAT.sizeWord(v.shape), unit: 'mm', v: 6, min: 1e-6, max: 1e5, scale: 'log', if: v => v.pore },
      { k: 'De', label: 'Effective diffusivity D_e', unit: 'm³/(m cat·[t])', v: 1e-6, min: 1e-20, max: 1e3, scale: 'log', if: v => v.pore },
      { k: 'dp2', label: 'Compare with another size (0 = none)', unit: 'mm', v: 0, min: 0, max: 1e5, if: v => v.pore },
      { k: 'rhoS', label: 'Pellet density ρ_s', unit: 'kg/m³ pellet', v: 2000, min: 1, max: 30000, if: v => v.pore || v.kb === 'pellet' || v.pat === 'fz' },
      { group: 'Bed' },
      { k: 'rhoB', label: 'Bulk density of the bed ρ_b', unit: 'kg/m³ bed', v: 1000, min: 1, max: 30000 },
      { k: 'dil', label: 'Inert solids per volume of catalyst (dilution)', v: 0, min: 0, max: 100 }
    ],
    presets: [
      { name: 'Example 18.4 — plug flow from the rate equation', v: { fb: 'FA0', FA0: 2000, CA0: 0.1, eps: 3, solve: 'W', X: 0.35, law: 'nth', kb: 'mass', k: 96, n: 1, pat: 'plug', pore: false, _cu: 'mol/L', _tu: 'h' }, book: { W: { v: 140, tol: 0.01 } }, note: 'A → 4R, 2000 mol/hr of pure A at 3.2 atm and 117 °C (C_A0 = 0.1 mol/L, ε_A = 3), −r′ = 96 C_A L/hr·kg, 35% conversion: W = 140 kg.' },
      { name: 'Example 18.5 — plug flow from rate data', v: { fb: 'FA0', FA0: 2000, CA0: 0.1, eps: 3, solve: 'W', X: 0.35, law: 'data', rateT: [[0.039, 3.4], [0.0575, 5.4], [0.075, 7.6], [0.092, 9.1]], pat: 'plug', pore: false, _cu: 'mol/L', _tu: 'h' }, book: { W: { v: 147, tol: 0.03 } }, note: 'Four measured rates and no rate equation: graphical integration of 1/(−r′) against X_A gives 0.0735 hr·kg/mol and W = 147 kg.' },
      { name: 'Example 18.6 — mixed flow (very large recycle)', v: { fb: 'FA0', FA0: 2000, CA0: 0.1, eps: 3, solve: 'W', X: 0.35, law: 'nth', kb: 'mass', k: 96, n: 1, pat: 'mixed', pore: false, _cu: 'mol/L', _tu: 'h' }, book: { W: { v: 230, tol: 0.01 } }, note: 'C_A = 0.0317 mol/L at the exit; W = 2000 × 0.35/(96 × 0.0317) = 230 kg.' },
      { name: 'Problem 18.11 — saturation kinetics', v: { fb: 'v0', v0: 1000, CA0: 100, eps: 0, solve: 'W', X: 0.8, law: 'lh', kb: 'mass', k: 50, K: 0.02, pat: 'plug', pore: false, _cu: 'mol/m³', _tu: 'h' }, note: '1000 m³/hr of pure A at 100 mol/m³, A → R, −r′ = 50C_A/(1 + 0.02C_A) mol/kg·hr; 80% conversion in a packed bed.' },
      { name: 'Problem 18.12 — second order', v: { fb: 'v0', v0: 1000, CA0: 100, eps: 0, solve: 'W', X: 0.8, law: 'nth', kb: 'mass', k: 8, n: 2, pat: 'plug', pore: false, _cu: 'mol/m³', _tu: 'h' }, note: 'Same duty with −r′ = 8C_A² mol/kg·hr.' },
      { name: 'Problem 18.13 — A + B, conversion from W', v: { fb: 'v0', v0: 10, CA0: 0.1, eps: 0, solve: 'X', Wgiven: 4, law: 'ab', kb: 'mass', k: 0.6, CB0: 10, bB: 1, pat: 'plug', pore: false, _cu: 'mol/m³', _tu: 'h' }, note: '10 m³/hr over 4 kg of catalyst, −r′ = 0.6 C_A C_B mol/kg·hr, C_A0 = 0.1 and C_B0 = 10 mol/m³.' },
      { name: 'Problem 18.17b,c — second order, diluted bed', v: { fb: 'v0', v0: 1000, CA0: 1, eps: 0, solve: 'W', X: 0.8, law: 'nth', kb: 'mass', k: 2000, n: 2, pat: 'plug', pore: false, rhoB: 1000, dil: 4, _cu: 'mol/L', _tu: 'h' }, note: 'k′ = 2 L²/mol·g·hr = 2000 L²/mol·kg·hr from part (a); 1000 L/hr at 1 mol/L, 80% conversion. Diluting 1 : 4 with inert solids leaves W unchanged and makes the bed five times larger.' },
      { name: 'Problem 18.31 — doubling the pellet size, strong pores', v: { fb: 'v0', v0: 1, CA0: 1, eps: 0, solve: 'X', Wgiven: 42.56, law: 'nth', kb: 'pellet', k: 5000, n: 1, pat: 'plug', pore: true, shape: 'sph', dp: 9, De: 1e-6, dp2: 18, rhoS: 2000, _cu: 'mol/m³', _tu: 's' }, note: 'Constants chosen so that 9 mm pellets are deep in the strong-diffusion regime and give 63.2% conversion; with 18 mm pellets the rate per pellet halves and the conversion falls to 1 − e^{−1/2} = 39.3%.' },
      { name: 'Problem 18.32 — second order, 1.2 cm pellets', v: { fb: 'v0', v0: 1, CA0: 6.67, eps: 0, solve: 'W', X: 0.8, law: 'nth', kb: 'mass', k: 0.01, n: 2, pat: 'plug', pore: true, shape: 'sph', dp: 12, De: 2e-6, dp2: 0, rhoS: 2000, _cu: 'mol/m³', _tu: 's' }, note: '1 m³/s of 1/3 A, 1/3 B, 1/3 inert at 336 °C and 1 atm (C_A0 = C_B0 = 6.67 mol/m³), A + B → R + S with k′ = 0.01 m⁶/mol·kg·s, ρ_s = 2000 kg/m³, D_e = 2 × 10⁻⁶ m³/m cat·s. Because C_A = C_B throughout, the rate is k′C_A².' },
      { name: 'Problem 18.33 — fluidized 1 mm or packed 1.5 cm', v: { fb: 'v0', v0: 1, CA0: 20, eps: 0, solve: 'W', X: 0.8, law: 'nth', kb: 'mass', k: 1.667e-3, n: 1, pat: 'mixed', pore: true, shape: 'sph', dp: 1, De: 1e-6, dp2: 15, rhoS: 2000, _cu: 'mol/m³', _tu: 's' }, note: 'The 1.2 mm experiment (80% in mixed flow, 10 g, 4 cm³/s) gives k′ ≈ 1.67 × 10⁻³ m³/kg·s once its small pore effect is removed (C-1802). The comparison table puts a mixed fluidized bed of 1 mm particles against a packed bed of 15 mm pellets.' },
      { name: 'Problem 18.34 — liquid, slow reaction', v: { fb: 'FA0', FA0: 1e4, CA0: 1000, eps: 0, solve: 'W', X: 0.9, law: 'nth', kb: 'bed', k: 8e-4 * 3600, n: 1, pat: 'plug', pore: true, shape: 'sph', dp: 6, De: 4e-8 * 3600, dp2: 0, rhoS: 4000, rhoB: 2000, _cu: 'mol/m³', _tu: 'h' }, note: 'k‴ = 8 × 10⁻⁴ m³/m³ bed·s, 6 mm pellets, D_e = 4 × 10⁻⁸ m³/m cat·s, bed voidage 0.5 and bulk density 2000 kg/m³ (so ρ_s = 4000 kg/m³); 10⁴ mol A/hr at 10³ mol/m³ to 90%. Per-second constants are entered per hour.' },
      { name: 'Problem 18.35 — liquid, fast reaction', v: { fb: 'FA0', FA0: 1e4, CA0: 1000, eps: 0, solve: 'W', X: 0.9, law: 'nth', kb: 'bed', k: 2 * 3600, n: 1, pat: 'plug', pore: true, shape: 'sph', dp: 6, De: 4e-8 * 3600, dp2: 0, rhoS: 4000, rhoB: 2000, _cu: 'mol/m³', _tu: 'h' }, note: 'As Problem 18.34 with k‴ = 2 m³/m³ bed·s: now strongly diffusion-limited.' },
      { name: 'Problem 18.36 — friable catalyst in a liquid fluidized bed', v: { fb: 'v0', v0: 1, CA0: 1, eps: 0, solve: 'X', Wgiven: 0.216, law: 'nth', kb: 'pellet', k: 2.996 * 2000, n: 1, pat: 'plug', pore: false, rhoS: 2000, _cu: 'mol/m³', _tu: 's' }, note: 'The bed gave 95% with 5 mm spheres (here W = 1, k chosen for that). After the spheres wear to 3 mm the same number of particles weighs 0.216 of the original: with no pore resistance X falls to 47.6%. Turn on pore diffusion in the strong regime to see the 66% of part (b).' },
      { name: 'Suspended solids, uniform f (Eq. 18.51)', v: { fb: 'v0', v0: 1, CA0: 1, eps: 0, solve: 'X', law: 'nth', kb: 'pellet', k: 2, n: 1, pat: 'fz', u0: 1, fz: 'const', Hz: 2, f0: 0.4, pore: false, rhoS: 2000, _cu: 'mol/m³', _tu: 's' }, note: 'First order with constant solids fraction: ln(C_A0/C_A) = k‴ f H/u₀ when ε_A = 0.' }
    ],
    eqs: [
      ['\\frac{W}{F_{A0}} = \\int_0^{X_{A,out}}\\frac{dX_A}{-r_A^\\prime},\\qquad \\tau^\\prime = \\frac{WC_{A0}}{F_{A0}}', 'Eqs. 18.41–18.42', 'Plug flow.'],
      ['k^\\prime\\tau^\\prime = (1+\\varepsilon_A)\\ln\\frac{1}{1-X_A} - \\varepsilon_A X_A', 'Eq. 18.44', 'Plug flow, first order.'],
      ['\\frac{W}{F_{A0}} = \\frac{X_{A,out}-X_{A,in}}{-r_{A,out}^\\prime},\\qquad k^\\prime\\tau^\\prime = \\frac{X_A(1+\\varepsilon_A X_A)}{1-X_A}', 'Eqs. 18.45–18.46', 'Mixed flow.'],
      ['\\frac{t}{C_{A0}} = \\frac{V}{W}\\int\\frac{dX_A}{-r_A^\\prime}', 'Eq. 18.47', 'Batch of gas, batch of catalyst.'],
      ['C_{A0}\\int_0^{X_A}\\frac{dX_A}{-r_A^{\\prime\\prime\\prime}} = \\frac{1}{u_0}\\int_0^H f\\,dz', 'Eq. 18.49', 'Solids fraction f varying with height.'],
      ['(1+\\varepsilon_A)\\ln\\frac{1}{1-X_A}-\\varepsilon_A X_A = \\frac{k^{\\prime\\prime\\prime}}{u_0}\\int_0^H f\\,dz', 'Eq. 18.50'],
      ['-r_A^{\\prime\\prime\\prime} = k^{\\prime\\prime\\prime}C_A^n\\,\\mathscr{E}(M_T),\\qquad k^{\\prime\\prime\\prime} = \\rho_s k^\\prime', 'Eq. 18.19', 'Effectiveness re-evaluated at every point of the reactor, with M_T from Eq. 18.28 at the local C_A.']
    ],
    method: ['Plug-flow integrals are evaluated by adaptive Gauss–Kronrod quadrature of 1/(−r′) against X_A; for a given catalyst load the conversion is found by Brent’s method on the same integral (or on the mixed-flow balance, reporting every root when the rate passes through a maximum). Measured rates are interpolated with a monotone cubic in log C – log r and extended as power laws beyond the data.', 'With pore diffusion the effectiveness factor is recomputed at the local concentration: from the exact nth-order pellet solution of C-1801, or, for other rate forms, from the generalized modulus of Eq. 18.26 evaluated by quadrature of the local rate law and placed on the exact first-order curve of the chosen shape.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const cu = v._cu || 'mol/L', tu = v._tu || 'h', V = KIT.volUnit(cu);
      const RM = rateModel(v), etaF = etaModel(v, RM);
      const CA0 = RM.CA0, eps = RM.eps;
      const FA0 = v.fb === 'FA0' ? +v.FA0 : +v.v0 * CA0, v0 = FA0 / CA0;
      const rhoB = +v.rhoB, dil = +v.dil || 0;
      let Xeq = 1;
      if (v.law === 'rev') { const Keq = +v.Keq; Xeq = Keq / (1 + Keq); kpi.push({ k: 'Xeq', l: 'Equilibrium conversion', v: Xeq }); }
      if (v.law === 'ab' && +v.bB > 0) Xeq = Math.min(1, (+v.CB0 / CA0) / +v.bB);
      const rate = x => RM.free(x) * etaF(x);
      const bedV = W => W / rhoB * (1 + dil);
      const volUnitBed = 'm³';
      /* ----- suspended solids f(z) ----- */
      if (v.pat === 'fz') {
        if (!(v.law === 'nth' && +v.n === 1)) throw new Error('Eqs. 18.49–18.51 are written for first-order reactions: choose order 1.');
        const k3 = v.kb === 'pellet' ? +v.k : v.kb === 'bed' ? +v.k / rhoB * +v.rhoS : +v.k * +v.rhoS;
        const H = +v.Hz, f0 = +v.f0;
        const fz = v.fz === 'exp' ? z => f0 * Math.exp(-(+v.az) * z) : v.fz === 'two' ? z => (z <= +v.Hd ? f0 : +v.fl) : () => f0;
        const I = v.fz === 'exp' ? (+v.az > 0 ? f0 * (1 - Math.exp(-(+v.az) * H)) / (+v.az) : f0 * H) : v.fz === 'two' ? f0 * Math.min(+v.Hd, H) + (+v.fl) * Math.max(H - +v.Hd, 0) : f0 * H;
        let eta = 1;
        if (v.pore) { const L = CAT.Lof(v.shape, +v.dp / 1000); eta = PEL.eta1(v.shape, L * Math.sqrt(k3 / +v.De)); kpi.push({ k: 'eta', l: 'Effectiveness factor', v: eta }); }
        const rhs = k3 * eta * I / +v.u0;
        const X = N.brent(x => (1 + eps) * Math.log(1 / (1 - x)) - eps * x - rhs, 0, 1 - 1e-15, { tol: 1e-14 });
        kpi.push({ k: 'X', l: 'Conversion X_A', v: X, main: true }, { k: 'If', l: '∫₀ᴴ f dz', v: I, u: 'm' }, { k: 'group', l: 'k‴ ℰ ∫f dz / u₀', v: rhs });
        const zs = N.linspace(0, H, 200);
        const Ic = N.cumtrapz(zs, zs.map(fz));
        const Xz = Ic.map(ii => { const g = k3 * eta * ii / +v.u0; return N.brent(x => (1 + eps) * Math.log(1 / (1 - x)) - eps * x - g, 0, 1 - 1e-15, { tol: 1e-13 }); });
        plots.push({ title: 'Solids fraction and conversion up the reactor', x: { label: 'height z', unit: 'm', min: 0, max: H }, y: { label: 'f, X_A', min: 0, max: 1 }, series: [{ name: 'solids fraction f', x: zs, y: zs.map(fz), color: 3, type: 'area', line: true }, { name: 'X_A', x: zs, y: Xz, color: 1, w: 2.2 }] });
        return { kpi, warn, info, plots, tables };
      }
      /* ----- solve for W or for X ----- */
      const pat = v.pat, R = +v.R;
      const Wfor = (X, etaFn, pp = pat) => {
        const rr = x => RM.free(x) * etaFn(x);
        if (pp === 'mixed') return FA0 * X / rr(X);
        if (pp === 'rec') { const Xi = R * X / (R + 1); return FA0 * (R + 1) * N.quad(x => 1 / rr(x), Xi, X, { rtol: 1e-8 }); }
        if (pp === 'batch') return NaN;
        return FA0 * N.quad(x => 1 / rr(x), 0, X, { rtol: 1e-9 });
      };
      const tBatch = (X, W) => {
        /* constant-volume batch of gas: C_A = C_A0 (1 − X) */
        const rr = x => { const C = CA0 * (1 - x); return (v.law === 'ab' ? RM.free(x) : RM.rC ? RM.rC(C) : RM.free(x)) * etaF(x); };
        return CA0 * (+v.Vg / W) * N.quad(x => 1 / rr(x), 0, X, { rtol: 1e-9 });
      };
      let X, W, roots = null;
      if (pat === 'batch') {
        W = v.solve === 'X' ? +v.Wgiven : +v.Wgiven || 1;
        X = v.solve === 'W' ? +v.X : +v.X;
        const t = tBatch(X, W);
        kpi.push({ k: 't', l: 'Time for X_A = ' + fmt(X, 3) + ' with W = ' + fmt(W, 4) + ' kg', v: t, u: tu, main: true });
        const xs = N.linspace(0, Math.min(X * 1.3, Xeq * 0.995), 120);
        plots.push({ title: 'Conversion against time in the batch', x: { label: 't', unit: tu, min: 0 }, y: { label: 'X_A', min: 0 }, series: [{ name: 'X_A(t)', x: xs.map(x => tBatch(x, W)), y: xs, color: 1, w: 2.2 }] });
        return { kpi, warn, info, plots, tables };
      }
      if (v.solve === 'W') {
        X = +v.X;
        if (!(X < Xeq)) throw new Error('The target conversion is beyond equilibrium (' + fmt(Xeq, 4) + ').');
        W = Wfor(X, etaF);
      } else {
        W = +v.Wgiven;
        if (pat === 'mixed') {
          const g = x => FA0 * x - W * rate(x);
          roots = N.roots(g, 1e-9, Math.min(Xeq, 1) * (1 - 1e-9), 600);
          if (!roots.length) throw new Error('No steady state found for this load.');
          X = roots[roots.length - 1];
          if (roots.length > 1) warn.push('With this rate form the mixed reactor has ' + roots.length + ' steady states: X_A = ' + roots.map(r => fmt(r, 4)).join(', ') + '. The highest is shown.');
        } else {
          const g = x => Wfor(x, etaF) - W;
          X = N.brent(g, 1e-12, Math.min(Xeq, 1) * (1 - 1e-12), { tol: 1e-13 });
        }
      }
      const tauP = W * CA0 / FA0;
      if (v.solve === 'W') kpi.push({ k: 'W', l: 'Catalyst needed W', v: W, u: 'kg', main: true, s: pat === 'mixed' ? 'mixed flow' : pat === 'rec' ? 'recycle R = ' + fmt(R, 3) : 'plug flow' });
      else kpi.push({ k: 'X', l: 'Conversion X_A', v: X, main: true, s: 'W = ' + fmt(W, 4) + ' kg' });
      kpi.push({ k: 'Vbed', l: 'Bed volume' + (dil ? ' (with inerts)' : ''), v: bedV(W), u: volUnitBed, s: 'ρ_b = ' + fmt(rhoB, 4) + ' kg/m³' + (dil ? ', ' + fmt(1 + dil, 3) + '× diluted' : '') });
      kpi.push({ k: 'tauP', l: 'Weight-time τ′ = W C_A0/F_A0', v: tauP, u: 'kg·' + tu + '/' + V });
      kpi.push({ k: 'CAout', l: 'C_A leaving', v: RM.Cx(X), u: cu });
      if (v.pore) {
        const e0 = etaF(0), e1 = etaF(X);
        kpi.push({ k: 'eta', l: 'Effectiveness factor, inlet → outlet', v: e0, s: '→ ' + fmt(e1, 4) + (v.law === 'nth' && +v.n === 1 ? ' (constant for first order)' : '') });
        const Wfree = Wfor(Math.min(X, Xeq * 0.999999), () => 1);
        kpi.push({ k: 'Wfree', l: 'Catalyst if the pellets were free of pore resistance', v: Wfree, u: 'kg' });
        if (+v.dp2 > 0) {
          const L2 = CAT.Lof(v.shape, +v.dp2 / 1000);
          const eta2 = x => etaF(x, L2);
          if (v.solve === 'W') kpi.push({ k: 'W2', l: 'Catalyst with ' + fmt(+v.dp2, 3) + ' mm pellets', v: Wfor(X, eta2), u: 'kg' });
          else {
            let X2;
            if (pat === 'mixed') { const rs = N.roots(x => FA0 * x - W * RM.free(x) * eta2(x), 1e-9, Xeq * (1 - 1e-9), 600); X2 = rs.length ? rs[rs.length - 1] : NaN; }
            else X2 = N.brent(x => Wfor(x, eta2) - W, 1e-12, Xeq * (1 - 1e-12), { tol: 1e-13 });
            kpi.push({ k: 'X2', l: 'Conversion with ' + fmt(+v.dp2, 3) + ' mm pellets', v: X2 });
          }
        }
      }
      /* comparison of patterns (and sizes) */
      const pats = [['plug', 'plug flow'], ['mixed', 'mixed flow']].concat(pat === 'rec' ? [['rec', 'recycle R = ' + fmt(R, 3)]] : []);
      const crow = [];
      const Xc = Math.min(X, Xeq * 0.999999);
      pats.forEach(([p, nm]) => {
        const Wf = Wfor(Xc, etaF, p); const W1 = v.pore ? Wf : NaN;
        const Wn = v.pore ? Wfor(Xc, () => 1, p) : Wf;
        let W2 = NaN;
        if (v.pore && +v.dp2 > 0) { const L2 = CAT.Lof(v.shape, +v.dp2 / 1000); W2 = Wfor(Xc, x => etaF(x, L2), p); }
        crow.push([nm, Wn, v.pore ? W1 : '—', v.pore && +v.dp2 > 0 ? W2 : '—', bedV(v.pore ? W1 : Wn)]);
      });
      tables.push({ title: 'Catalyst for X_A = ' + fmt(Xc, 4) + ' by flow pattern' + (v.pore ? ' and pellet' : ''), cols: ['flow pattern', 'W, no pore resistance', v.pore ? 'W, ' + fmt(+v.dp, 3) + ' mm pellets' : 'W with pellets', v.pore && +v.dp2 > 0 ? 'W, ' + fmt(+v.dp2, 3) + ' mm pellets' : 'W, second size', 'bed volume, m³'], rows: crow });
      /* Levenspiel plot */
      const xmax = Math.min(Xeq * 0.995, Math.max(X * 1.15, 0.05));
      const xs = N.linspace(0, xmax, 220);
      const inv = xs.map(x => 1 / rate(x)), invFree = xs.map(x => 1 / RM.free(x));
      const ymax = Math.max(...inv.filter(isFinite).slice(0, Math.max(2, Math.round(220 * Math.min(1, X / xmax)) + 1))) * 1.25;
      const ser = [{ name: v.pore ? '1/(−r′_A) with pore diffusion' : '1/(−r′_A)', x: xs, y: inv, color: 1, w: 2.2 }];
      if (v.pore) ser.push({ name: 'free of pore resistance', x: xs, y: invFree, color: 2, dash: '5 4' });
      const marks = [{ t: 'v', x: X, color: 8, label: 'X_A = ' + fmt(X, 3) }];
      if (pat === 'mixed') marks.push({ t: 'rect', x0: 0, x1: X, y0: 0, y1: 1 / rate(X), color: 3, fill: 0.14 });
      else if (pat === 'rec') { const Xi = R * X / (R + 1); const xa = xs.filter(x => x >= Xi && x <= X); ser.push({ name: 'recycle: area × (R + 1)', x: xa, y: xa.map(x => 1 / rate(x)), type: 'area', color: 3, fill: 0.16, legend: false }); }
      else { const xa = xs.filter(x => x <= X).concat([X]); ser.push({ name: 'area = W/F_A0', x: xa, y: xa.map(x => 1 / rate(x)), type: 'area', color: 3, fill: 0.16, legend: false }); }
      if (v.law === 'data' && RM.rC.dataRange) {
        const [c0, c1] = RM.rC.dataRange, x0 = Xof(CA0, eps)(c1), x1 = Xof(CA0, eps)(c0);
        marks.push({ t: 'band', x0, x1, color: 'var(--grid)' });
        info.push('The shaded band is covered by the measured rates; outside it the rate is extrapolated as a power law from the nearest two points.');
        const rows = (v.rateT || []).filter(q => q && +q[0] > 0 && +q[1] > 0);
        ser.push({ name: 'measured', x: rows.map(q => Xof(CA0, eps)(+q[0])), y: rows.map(q => 1 / +q[1]), type: 'pts', color: 4 });
      }
      plots.push({ title: 'Performance chart: 1/(−r′_A) against X_A (Fig. E18.5)', x: { label: 'X_A', min: 0, max: xmax }, y: { label: '1/(−r′_A)', unit: 'kg·' + tu + '/mol', min: 0, max: isFinite(ymax) ? ymax : undefined }, series: ser, marks, cap: pat === 'mixed' ? 'Mixed flow: the rectangle, rate at the exit composition.' : 'Plug flow: the shaded area is W/F_A0.' });
      /* profile along a plug bed */
      if (pat === 'plug' || pat === 'rec') {
        const Xp = N.linspace(0, X, 120), Wp = [0];
        for (let i = 1; i < Xp.length; i++) Wp.push(Wp[i - 1] + FA0 * N.quad(x => 1 / rate(x), Xp[i - 1], Xp[i], { rtol: 1e-8 }) * (pat === 'rec' ? 1 : 1));
        const sp = [{ name: 'X_A', x: Wp, y: Xp, color: 1, w: 2.2 }];
        if (v.pore) sp.push({ name: 'ℰ', x: Wp, y: Xp.map(x => etaF(x)), color: 5, dash: '5 4' });
        if (pat === 'plug') plots.push({ title: 'Along the packed bed', x: { label: 'catalyst passed W', unit: 'kg', min: 0 }, y: { label: 'X_A' + (v.pore ? ', ℰ' : ''), min: 0, max: 1 }, series: sp });
      }
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* C-1805  product distribution with pore diffusion                       */
  /* ===================================================================== */
  /* series A → R → S with linear pellet rates −r_A = kA C_A, r_R = a C_A − kR C_R (per pellet volume) */
  const serPlugMax = (kA, a, kR) => {
    if (Math.abs(kR - kA) < 1e-12 * kA) return { y: a / kA * Math.exp(-1), t: 1 / kA };
    const t = Math.log(kR / kA) / (kR - kA);
    return { y: a / (kR - kA) * (Math.exp(-kA * t) - Math.exp(-kR * t)), t };
  };
  const serMixMax = (kA, a, kR) => { const t = 1 / Math.sqrt(kA * kR); return { y: a * t / ((1 + kA * t) * (1 + kR * t)), t }; };
  const serPlugAt = (kA, a, kR, X) => { const t = -Math.log(1 - X) / kA; return Math.abs(kR - kA) < 1e-12 * kA ? a * t * Math.exp(-kA * t) : a / (kR - kA) * (Math.exp(-kA * t) - Math.exp(-kR * t)); };
  const serMixAt = (kA, a, kR, X) => { const t = X / ((1 - X) * kA); return a * t * (1 - X) / (1 + kR * t); };

  ST.tool({
    id: 'catselect', tag: 'C-1805', ch: 18, sec: 'Ch. 18.7', title: 'Product distribution when pores interfere: parallel and series reactions',
    lede: 'How strong pore diffusion changes what comes out of a catalytic reactor. For side-by-side reactions of any orders the pellet equation is solved with both paths inside it, giving the exact observed fractional yield at every concentration — not just Eq. 18.62’s estimate. For A → R → S the pellet is solved exactly for any Thiele modulus, so the yield of R moves continuously from the diffusion-free value to Wheeler’s strong-resistance limit; Table 18.2 is reproduced, and the falsified order, activation energy and k ratio of branched pore structures are computed.',
    keys: 'multiple reactions selectivity fractional yield pore diffusion parallel side by side series consecutive A R S wheeler table 18.2 maximum intermediate plug mixed strong resistance pore branching alpha falsified kinetics problem 18.28 18.29 18.30',
    inputs: [
      { group: 'Reaction scheme' },
      { k: 'sch', t: 'seg', label: 'Scheme', opts: [['ser', 'A → R → S'], ['par', 'A → R, A → S'], ['alpha', 'falsified kinetics (Eqs. 18.70–18.72)']], v: 'ser' },
      { k: 'k21', label: 'k₂/k₁ (true, at the same basis)', v: 1, min: 1e-4, max: 1e4, scale: 'log', if: v => v.sch === 'ser' },
      { k: 'DR', label: 'D_R/D_A (diffusivity ratio in the pores)', v: 1, min: 0.01, max: 100, scale: 'log', if: v => v.sch === 'ser' },
      { k: 'k1', label: 'k₁ for A → R', v: 2, min: 1e-12, max: 1e12, scale: 'log', if: v => v.sch === 'par' },
      { k: 'a1', label: 'Order a₁', v: 1, min: 0, max: 4, step: 0.5, if: v => v.sch === 'par' },
      { k: 'k2', label: 'k₂ for A → S', v: 3, min: 1e-12, max: 1e12, scale: 'log', if: v => v.sch === 'par' },
      { k: 'a2', label: 'Order a₂', v: 1, min: 0, max: 4, step: 0.5, if: v => v.sch === 'par' },
      { k: 'CA0', label: 'Feed C_A0', unit: 'mol/m³', v: 10, min: 1e-12, max: 1e7, if: v => v.sch === 'par' },
      { k: 'XA', label: 'Conversion of A', v: 0.9, min: 0.001, max: 0.9999, if: v => v.sch === 'par' },
      { group: 'Pellets', if: v => v.sch !== 'alpha' },
      { k: 'shape', t: 'seg', label: 'Shape', opts: [['slab', 'flat plate'], ['cyl', 'cylinder'], ['sph', 'sphere']], v: 'slab', if: v => v.sch !== 'alpha' },
      { k: 'MT1', label: v => (v.sch === 'ser' ? 'Thiele modulus of the first step, M_T1 = L√(k₁/D_A)' : 'Thiele modulus at the feed, L√((r₁ + r₂)/(D_A C_A0))'), v: 10, min: 1e-4, max: 1e4, scale: 'log', if: v => v.sch !== 'alpha' },
      { group: 'Inference from an experiment', if: v => v.sch === 'ser' },
      { k: 'inf', t: 'check', label: 'Find k₂/k₁ from an observed C_R,max/C_A0', v: false, if: v => v.sch === 'ser' },
      { k: 'obs', label: 'Observed C_R,max/C_A0', v: 0.23, min: 0.001, max: 0.999, if: v => v.sch === 'ser' && v.inf },
      { k: 'obsPat', t: 'seg', label: 'Observed in', opts: [['plug', 'plug flow'], ['mixed', 'mixed flow']], v: 'mixed', if: v => v.sch === 'ser' && v.inf },
      { k: 'obsReg', t: 'seg', label: 'Regime of that experiment', opts: [['free', 'diffusion-free'], ['strong', 'strong pore resistance'], ['given', 'the M_T above']], v: 'strong', if: v => v.sch === 'ser' && v.inf },
      { group: 'Pore structure', if: v => v.sch === 'alpha' },
      { k: 'alpha', label: 'Degree of branching α (0 nonporous, 1 one pore size, 2 two sizes …)', v: 1, min: 0, max: 6, step: 0.5, if: v => v.sch === 'alpha' },
      { k: 'nTrue', label: 'True order n', v: 2, min: 0, max: 4, step: 0.5, if: v => v.sch === 'alpha' },
      { k: 'Etrue', label: 'True activation energy E', unit: 'kJ/mol', v: 120, min: 0, max: 1000, if: v => v.sch === 'alpha' },
      { k: 'Ediff', label: 'E_diff', unit: 'kJ/mol', v: 10, min: 0, max: 200, if: v => v.sch === 'alpha' },
      { k: 'kr', label: '(k₂/k₁) true, side-by-side reactions', v: 9, min: 1e-6, max: 1e6, scale: 'log', if: v => v.sch === 'alpha' }
    ],
    presets: [
      { name: 'Table 18.2 — k₂/k₁ = 1', v: { sch: 'ser', k21: 1, DR: 1, shape: 'slab', MT1: 1000, inf: false }, book: { pNone: { v: 0.368, tol: 0.003 }, pStrong: { v: 0.184, tol: 0.005 }, mNone: { v: 0.25, tol: 0.003 }, mStrong: { v: 0.125, tol: 0.005 } }, note: 'Strong pore resistance halves the best yield of R in both flow patterns.' },
      { name: 'Table 18.2 — k₂/k₁ = 1/4', v: { sch: 'ser', k21: 0.25, DR: 1, shape: 'slab', MT1: 1000, inf: false }, book: { pNone: { v: 0.63, tol: 0.003 }, pStrong: { v: 0.333, tol: 0.005 }, mNone: { v: 0.444, tol: 0.003 }, mStrong: { v: 0.229, tol: 0.005 } } },
      { name: 'Table 18.2 — k₂/k₁ = 1/64', v: { sch: 'ser', k21: 1 / 64, DR: 1, shape: 'slab', MT1: 1000, inf: false }, book: { pNone: { v: 0.936, tol: 0.003 }, mNone: { v: 0.79, tol: 0.003 }, mStrong: { v: 0.486, tol: 0.005 } }, note: 'The book prints 0.650 for plug flow with strong resistance; Eq. 18.68 itself gives 0.660, which is what is computed here.' },
      { name: 'Table 18.2 — k₂/k₁ = 16', v: { sch: 'ser', k21: 16, DR: 1, shape: 'slab', MT1: 1000, inf: false }, book: { pNone: { v: 0.051, tol: 0.03 }, pStrong: { v: 0.031, tol: 0.03 }, mNone: { v: 0.04, tol: 0.01 }, mStrong: { v: 0.022, tol: 0.03 } } },
      { name: 'Intermediate pore resistance, sphere', v: { sch: 'ser', k21: 0.25, DR: 1, shape: 'sph', MT1: 1, inf: false }, note: 'Between the two columns of Table 18.2: the exact pellet solution gives the yield for any M_T.' },
      { name: 'Problem 18.28 — size does not matter', v: { sch: 'ser', k21: 1, DR: 1, shape: 'sph', MT1: 0.05, inf: true, obs: 0.5, obsPat: 'plug', obsReg: 'free' }, note: '4 mm and 8 mm pellets both give C_R,max/C_A0 = 0.5 in plug flow with W/F_A0 doubling with size — so diffusion-free, plug flow already, and k₂/k₁ = 0.5: nothing better is possible.' },
      { name: 'Problem 18.29 — strong resistance in a basket reactor', v: { sch: 'ser', k21: 1, DR: 1, shape: 'sph', MT1: 0.05, inf: true, obs: 0.23, obsPat: 'mixed', obsReg: 'strong' }, note: '6 mm and 12 mm pellets give the same 23% with W/F_A0 doubled: strong pore resistance. Infer k₂/k₁, then the best attainable — fine particles in plug flow.' },
      { name: 'Problem 18.30 — equal orders side by side', v: { sch: 'par', k1: 2, a1: 1, k2: 3, a2: 1, CA0: 10, XA: 0.9, shape: 'sph', MT1: 10 }, note: 'R and S both first order: the fractional yield is 0.4 whatever the regime or flow pattern, so C_R = 3.6 mol/m³ at 90% decomposition.' },
      { name: 'Different orders side by side', v: { sch: 'par', k1: 1, a1: 1, k2: 0.1, a2: 2, CA0: 10, XA: 0.9, shape: 'slab', MT1: 30 }, note: 'Here pore diffusion lowers the concentration inside the pellet and favours the lower-order path; the exact result differs from Eq. 18.62’s estimate.' },
      { name: 'Branched pores (Eq. 18.72)', v: { sch: 'alpha', alpha: 2, nTrue: 2, Etrue: 120, Ediff: 10, kr: 9 } }
    ],
    eqs: [
      ['\\varphi_{true} = \\frac{1}{1+(k_2/k_1)C_A^{a_2-a_1}}', 'Eq. 18.60', 'Equal orders: no change (Eq. 18.63).'],
      ['\\varphi_{obs} \\cong \\frac{1}{1+(k_2/k_1)C_{Ag}^{(a_2-a_1)/2}}', 'Eq. 18.62', 'Strong pore resistance, estimate.'],
      ['\\frac{dC_{Rg}}{dC_{Ag}} = -\\frac{1}{1+\\gamma} + \\gamma\\frac{C_{Rg}}{C_{Ag}},\\qquad \\gamma = (k_2/k_1)^{1/2}', 'Eq. 18.65', 'Strong pore resistance, A → R → S.'],
      ['\\frac{C_{R,max}}{C_{A0}} = \\frac{\\gamma^{\\gamma/(1-\\gamma)}}{1+\\gamma}', 'Eq. 18.68', 'Plug flow, strong resistance.'],
      ['\\frac{C_{R,max}}{C_{A0}} = \\frac{1}{(1+\\gamma)(\\gamma^{1/2}+1)^2}', 'Eq. 18.69', 'Mixed flow, strong resistance.'],
      ['-r_A = k_1\\mathscr{E}_1 C_{As},\\qquad r_R = aC_{As} - k_2\\mathscr{E}_2C_{Rs}', 'Exact pellet', 'Any modulus; the coupling a is below.'],
      ['a = k_1\\mathscr{E}_1 - \\frac{k_1k_2(\\mathscr{E}_1-\\mathscr{E}_2)}{k_2-k_1\\mathscr{D}_R/\\mathscr{D}_A}', 'Exact pellet', 'Reduces to Eq. 18.65 when both moduli are large.'],
      ['\\mathbf{E}_{obs} = \\left(1-\\frac{1}{2^\\alpha}\\right)\\mathbf{E}_{diff} + \\frac{\\mathbf{E}}{2^\\alpha}', 'Eq. 18.72'],
      ['n_{obs} = 1+\\frac{n-1}{2^\\alpha},\\qquad \\left(\\frac{k_2}{k_1}\\right)_{obs} = \\left(\\frac{k_2}{k_1}\\right)^{1/2^\\alpha}', 'Eqs. 18.70–18.72']
    ],
    method: ['Series reactions: because both steps are first order, the pellet equations for A and R are linear and solved exactly: C_A follows the first-order profile and C_R is a combination of the profiles for k₁ and k₂. The pellet then acts on the bulk stream through two effectiveness factors and one coupling coefficient, and the plug and mixed reactor equations are solved in closed form, maximized analytically over space time.', 'Parallel reactions: the pellet equation with both rates is shot from the centre with a second quadrature carrying the formation of R, so the observed fractional yield is exact at any modulus; in the strong limit it tends to the thin-zone value ∫w₁/√(2F) dc / √(2F(1)). Reactor yields integrate the local fractional yield over the concentration path (plug) or take it at the exit (mixed).'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      if (v.sch === 'alpha') {
        const a = +v.alpha, p = Math.pow(2, a), n = +v.nTrue, E = +v.Etrue, Ed = +v.Ediff, kr = +v.kr;
        kpi.push({ k: 'Eobs', l: 'Observed activation energy', v: (1 - 1 / p) * Ed + E / p, u: 'kJ/mol', main: true });
        kpi.push({ k: 'nobs', l: 'Observed order', v: 1 + (n - 1) / p });
        kpi.push({ k: 'krobs', l: 'Observed (k₂/k₁), side by side', v: Math.pow(kr, 1 / p) });
        const as = N.linspace(0, 5, 101);
        plots.push({ title: 'What you would measure as the pore structure branches', x: { label: 'α', min: 0, max: 5 }, y: { label: 'E_obs/E, n_obs, (k₂/k₁)_obs / (k₂/k₁)', min: 0 }, series: [{ name: 'E_obs / E', x: as, y: as.map(q => ((1 - Math.pow(2, -q)) * Ed + E * Math.pow(2, -q)) / E), color: 1 }, { name: 'n_obs', x: as, y: as.map(q => 1 + (n - 1) * Math.pow(2, -q)), color: 2 }, { name: '(k₂/k₁)_obs / (k₂/k₁)', x: as, y: as.map(q => Math.pow(kr, Math.pow(2, -q)) / kr), color: 3 }], marks: [{ t: 'v', x: a, color: 8, label: 'α = ' + fmt(a, 2) }], cap: 'As α grows the measured E falls toward E_diff and the order toward one: only experiments free of pore resistance reveal the true kinetics (Ch. 18.7).' });
        return { kpi, warn, info, plots, tables };
      }
      const sh = v.shape;
      if (v.sch === 'ser') {
        let k21 = +v.k21;
        const DR = +v.DR, M1 = +v.MT1;
        /* pellet coefficients for a given modulus; L = 1, D_A = 1, k1 = M1² */
        const coef = (M, kr) => { const k1 = M * M; return PEL.series(sh, 1, k1, k1 * kr, 1, DR); };
        const best = (M, kr) => { const c = coef(M, kr); return { plug: serPlugMax(c.kA, c.a, c.kR), mixed: serMixMax(c.kA, c.a, c.kR), c }; };
        if (v.inf) {
          const target = +v.obs, pat = v.obsPat, reg = v.obsReg;
          const Mexp = reg === 'free' ? 1e-4 : reg === 'strong' ? 1e4 : M1;
          const g = lk => best(Mexp, Math.exp(lk))[pat].y - target;
          const lk = N.brent(g, Math.log(1e-6), Math.log(1e6), { tol: 1e-12 });
          if (!isFinite(lk)) throw new Error('No k₂/k₁ gives that yield in this regime.');
          k21 = Math.exp(lk);
          kpi.push({ k: 'k21', l: 'k₂/k₁ inferred from the experiment', v: k21, main: true, s: 'from C_R,max/C_A0 = ' + fmt(target, 3) + ' in ' + pat + ' flow, ' + (reg === 'free' ? 'diffusion-free' : reg === 'strong' ? 'strong pore resistance' : 'M_T1 = ' + fmt(M1, 3)) });
          const bf = best(1e-4, k21);
          kpi.push({ k: 'bestPossible', l: 'Best attainable: fine particles, plug flow', v: bf.plug.y, s: 'mixed flow, fine particles: ' + fmt(bf.mixed.y, 3) });
        }
        const free = best(1e-4, k21), strong = best(1e4, k21), here = best(M1, k21);
        kpi.push({ k: 'pHere', l: 'C_R,max/C_A0, plug flow, these pellets', v: here.plug.y, main: !v.inf, s: 'M_T1 = ' + fmt(M1, 3) + ', M_T2 = ' + fmt(here.c.M2, 3) });
        kpi.push({ k: 'mHere', l: 'C_R,max/C_A0, mixed flow, these pellets', v: here.mixed.y });
        kpi.push({ k: 'pNone', l: 'Plug flow, no pore resistance (Eq. 8.8)', v: free.plug.y });
        kpi.push({ k: 'pStrong', l: 'Plug flow, strong resistance (Eq. 18.68)', v: strong.plug.y, s: 'decrease ' + fmt((1 - strong.plug.y / free.plug.y) * 100, 3) + '%' });
        kpi.push({ k: 'mNone', l: 'Mixed flow, no pore resistance (Eq. 8.15)', v: free.mixed.y });
        kpi.push({ k: 'mStrong', l: 'Mixed flow, strong resistance (Eq. 18.69)', v: strong.mixed.y, s: 'decrease ' + fmt((1 - strong.mixed.y / free.mixed.y) * 100, 3) + '%' });
        if (Math.abs(DR - 1) > 1e-9) info.push('With D_R ≠ D_A the strong-resistance limit is not exactly Eq. 18.68: the exact pellet solution is used throughout.');
        /* charts */
        const Xs = N.linspace(0.001, 0.999, 200);
        const curves = [['no pore resistance', free.c, 2, '5 4'], ['these pellets', here.c, 1, undefined], ['strong resistance', strong.c, 3, '2 3']];
        plots.push({ title: 'C_R/C_A0 against conversion of A', x: { label: 'X_A', min: 0, max: 1 }, y: { label: 'C_R/C_A0', min: 0 }, series: curves.map(([nm, c, col, d]) => ({ name: 'plug, ' + nm, x: Xs, y: Xs.map(x => serPlugAt(c.kA, c.a, c.kR, x)), color: col, dash: d, w: 2 })).concat(curves.map(([nm, c, col]) => ({ name: 'mixed, ' + nm, x: Xs, y: Xs.map(x => serMixAt(c.kA, c.a, c.kR, x)), color: col, dash: '1 3', w: 1.4, legend: false }))), cap: 'Solid and dashed: plug flow; dotted: mixed flow, in the same colours.' });
        const Ms = N.logspace(0.01, 100, 90);
        plots.push({ title: 'Best yield of R against the Thiele modulus', x: { label: 'M_T1 = L√(k₁/D_A)', log: true, min: 0.01, max: 100 }, y: { label: 'C_R,max/C_A0', min: 0 }, series: [{ name: 'plug flow', x: Ms, y: Ms.map(M => best(M, k21).plug.y), color: 1, w: 2.2 }, { name: 'mixed flow', x: Ms, y: Ms.map(M => best(M, k21).mixed.y), color: 2, w: 2.2 }], marks: [{ t: 'v', x: M1, color: 8, label: 'these pellets' }, { t: 'h', y: free.plug.y / 2, color: 'var(--muted)', label: 'half the diffusion-free plug value' }], cap: 'The yield falls by about half between the diffusion-free and the strong-resistance regimes (Table 18.2).' });
        /* Table 18.2 */
        const tb = [1 / 64, 1 / 16, 1 / 4, 1, 4, 16].map(kr => { const f = best(1e-4, kr), s = best(1e4, kr); return [kr, f.plug.y, s.plug.y, (1 - s.plug.y / f.plug.y) * 100, f.mixed.y, s.mixed.y, (1 - s.mixed.y / f.mixed.y) * 100]; });
        tables.push({ title: 'Table 18.2 recomputed (first-order reactions in series)', cols: ['k₂/k₁', 'plug, none', 'plug, strong', '% decrease', 'mixed, none', 'mixed, strong', '% decrease'], rows: tb });
        return { kpi, warn, info, plots, tables };
      }
      /* parallel */
      const k1 = +v.k1, a1 = +v.a1, k2 = +v.k2, a2 = +v.a2, CA0 = +v.CA0, XA = +v.XA, MT0 = +v.MT1;
      /* fix D and L from the modulus at the feed: L √((r1+r2)/(D C)) = MT0 with L = 1 */
      const r0 = k1 * Math.pow(CA0, a1) + k2 * Math.pow(CA0, a2), D = r0 / (CA0 * MT0 * MT0);
      const Cs = N.logspace(CA0 * (1 - XA) * 0.98, CA0, 24);
      const obs = Cs.map(C => PEL.parallel(sh, C, [k1, k2], [a1, a2], D, 1));
      const phiT = C => 1 / (1 + (k2 / k1) * Math.pow(C, a2 - a1));
      const phi62 = C => 1 / (1 + (k2 / k1) * Math.pow(C, (a2 - a1) / 2));
      const phiO = N.pchip(Cs.map(Math.log), obs.map(q => q.phiObs));
      const phiOf = C => phiO(Math.log(C));
      const Cout = CA0 * (1 - XA);
      const plugY = f => N.quad(c => f(c), Cout, CA0, { rtol: 1e-9 });
      const CRp = { free: plugY(phiT), here: plugY(phiOf), est: plugY(phi62) };
      const CRm = { free: phiT(Cout) * (CA0 - Cout), here: phiOf(Cout) * (CA0 - Cout), est: phi62(Cout) * (CA0 - Cout) };
      kpi.push({ k: 'CRp', l: 'C_R leaving, plug flow, these pellets', v: CRp.here, u: 'mol/m³', main: true, s: 'diffusion-free ' + fmt(CRp.free, 4) + '; Eq. 18.62 estimate ' + fmt(CRp.est, 4) });
      kpi.push({ k: 'CRm', l: 'C_R leaving, mixed flow, these pellets', v: CRm.here, u: 'mol/m³', s: 'diffusion-free ' + fmt(CRm.free, 4) + '; Eq. 18.62 estimate ' + fmt(CRm.est, 4) });
      kpi.push({ k: 'phiIn', l: 'Fractional yield at the feed: true → observed', v: obs[obs.length - 1].phiObs, s: 'true ' + fmt(phiT(CA0), 4) + ', ℰ = ' + fmt(obs[obs.length - 1].eta, 4) });
      if (a1 === a2) info.push('Equal orders: the fractional yield is the same inside the pellet as outside, so pore diffusion and flow pattern do not change the product distribution (Eq. 18.63).');
      else info.push(a1 < a2 ? 'The desired path is of lower order: the low concentration inside pellets favours it, so here pore resistance helps the yield.' : 'The desired path is of higher order: the low concentration inside pellets hurts it; use small pellets and plug flow.');
      plots.push({ title: 'Fractional yield of R against concentration of A', x: { label: 'C_A at the pellet surface', unit: 'mol/m³', log: true }, y: { label: 'φ(R/A)', min: 0, max: 1 }, series: [{ name: 'true (no pore resistance)', x: Cs, y: Cs.map(phiT), color: 2, dash: '5 4' }, { name: 'observed, exact pellet', x: Cs, y: obs.map(q => q.phiObs), color: 1, w: 2.2 }, { name: 'Eq. 18.62 estimate', x: Cs, y: Cs.map(phi62), color: 3, dash: '2 3' }] });
      tables.push({ title: 'Along the concentration path', cols: ['C_A', 'φ true', 'φ observed (exact)', 'φ Eq. 18.62', 'ℰ', 'M_T (surface rate)'], rows: obs.filter((q, i) => i % 3 === 0 || i === obs.length - 1).map((q, i) => { const C = Cs[obs.indexOf(q)]; return [C, phiT(C), q.phiObs, phi62(C), q.eta, q.MT]; }) });
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
