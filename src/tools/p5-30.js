/* Chapter 28 — the overall picture: fractional yields, what ends a batch, the Monod and poisoning factors.
   Chapter 30 — product-limiting fermentation: batch and plug flow (Eq. 30.5), mixed flow (Eqs. 30.9–30.14),
   the optimum and the best way to run a plug-flow fermentor (Fig. 30.2). Uses ST.BIO from p5-29. */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, BIO = ST.BIO;
  const fmt = (x, s = 4) => N.fmt(x, s);

  /* ===================================================================== */
  /* B-2801  yields, growth and what stops it                               */
  /* ===================================================================== */
  ST.tool({
    id: 'microbe', tag: 'B-2801', ch: 28, title: 'Microbial fermentation: fractional yields, batch growth and what stops it',
    lede: 'Give the fractional yields, or a stoichiometric equation such as 10A + Cell → 18R + 2 Cells, and get every φ, the most cells and product a batch can make, and whether it ends because the food runs out or because its own product poisons it. Batch curves for feeds from a quarter to four times C_A0 reproduce Fig. 28.3; the semi-log growth curve shows the lag and the slope k of Fig. 28.1.',
    keys: 'microbial fermentation introduction overall picture fractional yield phi stoichiometry cells substrate food product waste lag phase exponential growth stationary phase batch fermentor maximum cell concentration famine pollution generalized monod equation product poisoning order k observed doubling time figure 28.1 28.3 28.4 wine yeast alcohol 10A + cell 18R 2 cells',
    units: { t: 'h' },
    inputs: [
      { group: 'Stoichiometry' },
      { k: 'ysrc', t: 'seg', label: 'Yields from', opts: [['phi', 'fractional yields'], ['eq', 'an equation']], v: 'phi' },
      { k: 'yCA', label: 'φ(C/A), cells formed per A used', v: 0.1, min: 1e-9, max: 1e9, if: v => v.ysrc === 'phi' },
      { k: 'yRA', label: 'φ(R/A), product formed per A used', v: 0.4, min: 0, max: 1e9, if: v => v.ysrc === 'phi' },
      { k: 'sa', label: 'aA + C → rR + cC: a, A used', v: 10, min: 1e-9, max: 1e9, if: v => v.ysrc === 'eq' },
      { k: 'sr', label: 'r, product formed', v: 18, min: 0, max: 1e9, if: v => v.ysrc === 'eq' },
      { k: 'sc', label: 'c, cells leaving per cell entering', v: 2, min: 1, max: 1e9, if: v => v.ysrc === 'eq' }
    ].concat(BIO.inputs({ k: 1.2, CM: 2, CA0: 6, CC0: 0.05 }).filter(i => i.k !== 'yCA' && i.k !== 'yRA')).concat([
      { group: 'Batch' },
      { k: 'lag', label: 'Lag before growth starts', unit: '[t]', v: 1, min: 0, max: 1e9 },
      { k: 'fam', t: 'check', label: 'Compare feeds from ¼ to 4 times C_A0', v: true }
    ]),
    presets: [
      { name: 'Fig. 28.3 — the food runs out', v: { ysrc: 'phi', yCA: 0.1, yRA: 0.4, k: 1.2, CM: 2, pois: false, CA0: 6, CC0: 0.05, CR0: 0, lag: 1, fam: true, _tu: 'h' }, expect: { CAmr: { v: Math.sqrt(4 + 2 * (6 + 0.05 / 0.1)) - 2, tol: 1e-6 } }, check: { stop: 'food', CCmax: { v: 0.65, tol: 1e-9 } }, note: 'The final cell concentration rises in proportion to C_A0 (Fig. 28.3, left). The C_A of fastest growth is Eq. 29.2.' },
      { name: 'Fig. 28.3 — the product poisons the cells', v: { ysrc: 'phi', yCA: 0.1, yRA: 0.5, k: 1, CM: 0.5, pois: true, CRs: 2, n: 1, CA0: 10, CC0: 0.05, CR0: 0, lag: 0.5, fam: true, _tu: 'h' }, check: { stop: 'poison', CCmax: { v: 0.45, tol: 1e-9 } }, note: 'Once C_A0 exceeds C_R*/φ(R/A) = 4, every batch stops at the same cell concentration (Fig. 28.3, right).' },
      { name: 'Problem 30.8 — 10A + Cell → 18R + 2 Cells', v: { ysrc: 'eq', sa: 10, sr: 18, sc: 2, k: 2, CM: 0, pois: true, CRs: 36, n: 1, CA0: 250, CC0: 0.1, CR0: 0, lag: 0, fam: false, _tu: 'h' }, expect: { yCA: { v: 0.1, tol: 1e-12 }, yRA: { v: 1.8, tol: 1e-12 }, yRC: { v: 18, tol: 1e-12 } }, check: { stop: 'poison' }, note: 'One cell eats 10 A and becomes two, making 18 R. With C_R* = 36, only 20 of the 250 mol/m³ of A can ever be used.' },
      { name: 'Wine — the yeast quits near 12% alcohol', v: { ysrc: 'phi', yCA: 0.05, yRA: 0.51, k: 0.25, CM: 2, pois: true, CRs: 95, n: 1, CA0: 220, CC0: 0.5, CR0: 0, lag: 4, fam: true, _tu: 'h' }, note: 'Illustrative numbers, kg/m³: 220 of sugar, ethanol at 0.51 kg per kg sugar, C_R* = 95 kg/m³ (about 12% by volume). The richest musts stop sweet; below about 186 kg/m³ of sugar the wine ferments dry.' }
    ],
    eqs: [
      ['\\varphi(R/A) = \\varphi(R/C)\\,\\varphi(C/A)', 'Eq. 28.2'],
      ['\\varphi(A/C) = 1/\\varphi(C/A)', 'Eq. 28.2'],
      ['r_C = (-r_A)\\,\\varphi(C/A),\\qquad r_R = r_C\\,\\varphi(R/C)', 'Eq. 28.3'],
      ['C_C = C_{C0}+\\varphi(C/A)(C_{A0}-C_A)', 'Eq. 28.4'],
      ['C_R = C_{R0}+\\varphi(R/A)(C_{A0}-C_A)', 'Eq. 28.4'],
      ['r_C = k\\frac{C_AC_C}{C_A+C_M}', 'Eq. 28.6', 'Monod.'],
      ['k_{obs} = k\\left(1-\\frac{C_R}{C_R^*}\\right)^n', 'Eq. 28.7'],
      ['C_{C,max} = C_{C0} + \\varphi(C/A)\\min\\left(C_{A0},\\ \\frac{C_R^*-C_{R0}}{\\varphi(R/A)}\\right)', 'Fig. 28.3', 'Food or poison, whichever gives out first.']
    ],
    method: ['Constant yields tie every concentration to the amount of A used (Eq. 28.4), so a batch is one differential equation, integrated after the lag with an embedded Runge–Kutta method. A batch can go on only until A is gone or C_R reaches C_R*, whichever comes first; that fixes the final composition without integrating anything.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'h';
      const pv = Object.assign({}, v);
      if (v.ysrc === 'eq') {
        pv.yCA = (+v.sc - 1) / +v.sa; pv.yRA = +v.sr / +v.sa;
        if (!(pv.yCA > 0)) throw new Error('Cells must multiply: c has to exceed 1.');
      }
      const m = BIO.model(BIO.params(pv));
      kpi.push({ k: 'yCA', l: 'φ(C/A), cells per A', v: m.yCA, s: 'φ(A/C) = ' + fmt(1 / m.yCA, 4) });
      kpi.push({ k: 'yRA', l: 'φ(R/A), product per A', v: m.yRA });
      kpi.push({ k: 'yRC', l: 'φ(R/C), product per cell', v: m.yRA / m.yCA });
      const dPois = isFinite(m.CRs) ? (m.CRs - m.CR0) / Math.max(m.yRA, 1e-300) : Infinity;
      if (isFinite(m.CRs) && m.CR0 >= m.CRs) throw new Error('C_R0 is already at or above C_R*: nothing can grow.');
      const pois = dPois < m.CA0;
      kpi.unshift({ k: 'stop', l: 'What ends a batch', v: pois ? 'poison' : 'food', main: true, s: pois ? 'the product reaches C_R* with C_A = ' + fmt(m.CA0 - dPois, 4) + ' uneaten' : isFinite(m.CRs) ? 'the food runs out; C_R ends at ' + fmt(m.CR(m.CA0), 4) + ', below C_R* = ' + fmt(m.CRs, 4) : 'the food runs out (no poisoning in this rate form)' });
      kpi.push({ k: 'CCmax', l: 'Most cells a batch can reach, C_C,max', v: m.CC(m.Dmax) });
      if (m.yRA > 0) kpi.push({ k: 'CRmax', l: 'Most product, C_R,max', v: m.CR(m.Dmax) });
      const kobs = d => { const a = m.CA(d); return m.k * (isFinite(m.CRs) ? Math.pow(Math.max(1 - m.CR(d) / m.CRs, 0), m.n) : 1) * (m.CM > 0 ? a / (a + m.CM) : 1); };
      kpi.push({ k: 'kobs0', l: 'Specific growth rate in the fresh feed', v: kobs(0), u: '1/' + tu, s: 'cells double every ' + fmt(Math.LN2 / kobs(0), 4) + ' ' + tu });
      const dO = m.dMax();
      kpi.push({ k: 'CAmr', l: 'C_A where cells grow fastest (Eq. 29.2 without poison)', v: m.CA(dO), s: 'r_C,max = ' + fmt(m.rC(dO), 4) + ' per ' + tu });
      if (!(m.CC0 > 0)) {
        warn.push('A batch with no cells never starts: set C_C0 above zero to see it grow.');
        return { kpi, warn, info, plots, tables };
      }
      const lag = +v.lag;
      const t99 = m.tPlug(0, 0.99 * m.Dmax);
      kpi.push({ k: 't99', l: 'Batch time to 99% of the way', v: lag + t99, u: tu, s: (lag > 0 ? fmt(lag, 3) + ' lag + ' : '') + fmt(t99, 4) + ' growing' });
      const fac = v.fam ? [0.25, 0.5, 1, 2, 4] : [1];
      const models = fac.map(f => BIO.model(Object.assign(BIO.params(pv), { CA0: m.CA0 * f })));
      const T = lag + 1.35 * Math.max(...models.map(mm => mm.tPlug(0, 0.99 * mm.Dmax)));
      const paths = models.map(mm => BIO.batchPath(mm, T, 180, lag));
      plots.push({ title: 'Batch growth for a range of feeds (Fig. 28.3)', x: { label: 't', unit: tu, min: 0, max: T }, y: { label: 'C_C', min: 0 },
        series: paths.map((p, i) => ({ name: 'C_A0 = ' + fmt(models[i].CA0, 3), x: p.t, y: p.d.map(models[i].CC), color: fac[i] === 1 ? 1 : i + 2, w: fac[i] === 1 ? 2.6 : 1.6 })),
        marks: lag > 0 ? [{ t: 'band', x0: 0, x1: lag, color: 'var(--grid)' }] : [],
        cap: pois ? 'Poison-limited: past C_A0 = ' + fmt(m.CR0 + dPois, 3) + ' every batch levels off at the same C_C.' : 'Food-limited: the plateau rises with C_A0.' });
      const ca = N.linspace(0, 4 * m.CA0, 200);
      const ccMax = c => m.CC0 + m.yCA * Math.min(c, dPois);
      const mk = [];
      if (isFinite(dPois) && dPois < 4 * m.CA0) mk.push({ t: 'v', x: dPois, color: 'var(--muted)', label: 'C_R* reached' });
      plots.push({ title: 'Largest cell concentration a batch can reach (Fig. 28.3)', x: { label: 'C_A0', min: 0 }, y: { label: 'C_C,max', min: 0 },
        series: [{ name: 'C_C,max', x: ca, y: ca.map(ccMax), color: 1, w: 2.4 }, { name: 'if food were the only limit', x: ca, y: ca.map(c => m.CC0 + m.yCA * c), color: 'var(--muted)', dash: '4 4' }, { name: 'these feeds', x: models.map(mm => mm.CA0), y: models.map(mm => ccMax(mm.CA0)), type: 'pts', color: 2 }],
        marks: mk });
      const base = paths[fac.indexOf(1)];
      const lnC = base.d.map(d => Math.log(m.CC(d) / m.CC0));
      const Y = Math.max(...lnC), k0 = kobs(0);
      plots.push({ title: 'ln(C_C/C_C0) against time: lag, then the steepest slope k_obs (Fig. 28.1)', x: { label: 't', unit: tu, min: 0, max: T }, y: { label: 'ln(C_C/C_C0)', min: 0 },
        series: [{ name: 'ln(C_C/C_C0)', x: base.t, y: lnC, color: 1, w: 2.4 }],
        marks: [{ t: 'seg', x0: lag, y0: 0, x1: lag + Y / k0, y1: Y, color: 2, dash: '5 4' }, { t: 'text', x: lag + 0.55 * Y / k0, y: 0.5 * Y, s: 'slope ' + fmt(k0, 3), color: 2 }].concat(lag > 0 ? [{ t: 'band', x0: 0, x1: lag, color: 'var(--grid)' }, { t: 'text', x: lag / 2, y: 0.85 * Y, s: 'lag', anchor: 'middle' }] : []) });
      if (m.CM > 0) {
        const cg = N.linspace(0, Math.max(m.CA0, 10 * m.CM), 200);
        plots.push({ title: 'Growth constant against food (Fig. 28.1)', x: { label: 'C_A', min: 0 }, y: { label: 'k_obs', unit: '1/' + tu, min: 0, max: m.k * 1.08 },
          series: [{ name: 'k C_A/(C_A + C_M)', x: cg, y: cg.map(c => m.k * c / (c + m.CM)), color: 1, w: 2.4 }],
          marks: [{ t: 'h', y: m.k, color: 'var(--muted)', label: 'k' }, { t: 'h', y: m.k / 2, color: 'var(--muted)', label: 'k/2' }, { t: 'v', x: m.CM, color: 2, label: 'C_M' }] });
      }
      if (isFinite(m.CRs)) {
        const q = N.linspace(0, 1, 101);
        plots.push({ title: 'Poisoning by product (Fig. 28.4)', x: { label: 'C_R/C_R*', min: 0, max: 1 }, y: { label: 'k_obs/k', min: 0, max: 1.05 },
          series: [0.5, 1, 2].map((nn, i) => ({ name: 'n = ' + nn, x: q, y: q.map(z => Math.pow(1 - z, nn)), color: i + 2, w: 1.5 })).concat([{ name: 'n = ' + fmt(m.n, 3), x: q, y: q.map(z => Math.pow(1 - z, m.n)), color: 1, w: 2.6 }]) });
      }
      tables.push({ title: 'Final state of each batch', cols: ['C_A0', 'ends by', 'C_A left', 'C_C', 'C_R', 'time to 99%'], rows: models.map(mm => [mm.CA0, mm.Dmax < mm.CA0 * (1 - 1e-12) ? 'poison' : 'food', mm.CA(mm.Dmax), mm.CC(mm.Dmax), mm.CR(mm.Dmax), lag + mm.tPlug(0, 0.99 * mm.Dmax)]) });
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* B-3001  product-limiting fermentation                                  */
  /* ===================================================================== */
  function pmodel(v) {
    const yRC = +v.yRC;
    if (v.food) return BIO.model({ k: +v.k, CM: +v.CM, CRs: +v.CRs, n: +v.n, yCA: +v.yRA / yRC, yRA: +v.yRA, CA0: +v.CA0, CC0: +v.CC0, CR0: +v.CR0 });
    return BIO.model({ k: +v.k, CM: 0, CRs: +v.CRs, n: +v.n, yCA: 1 / yRC, yRA: 1, CA0: 1e300, CC0: +v.CC0, CR0: +v.CR0 });
  }
  /* closed forms for plenty of food (Eqs. 30.4, 30.5, 30.13) */
  const PL = {
    a: v => +v.CR0 - +v.yRC * +v.CC0,
    CRmr: v => Math.max(+v.CR0, (+v.CRs + +v.n * PL.a(v)) / (+v.n + 1)),
    plug1: (v, kt) => { const C = +v.CRs, a = PL.a(v), R0 = +v.CR0; if (!(R0 - a > 0)) return R0; const Q = Math.exp(kt * (C - a) / C) * (R0 - a) / (C - R0); return (a + Q * C) / (1 + Q); },
    tplug1: (v, x) => { const C = +v.CRs, a = PL.a(v), R0 = +v.CR0; return C / (C - a) * Math.log((x - a) * (C - R0) / ((R0 - a) * (C - x))) / +v.k; },
    ktMix: (v, x) => (x - +v.CR0) / ((x - PL.a(v)) * Math.pow(1 - x / +v.CRs, +v.n))
  };
  ST.tool({
    id: 'poison', tag: 'B-3001', ch: 30, sec: '30.1–30.2', title: 'Product-limiting fermentation: batch, plug and mixed fermentors',
    lede: 'Cells slowed by their own product, r_R = k(1 − C_R/C_R*)ⁿ(C_R − C_R0 + φ(R/C)C_C0), with food in excess or limiting too. Batch and plug flow follow Eq. 30.5, mixed flow Eqs. 30.9–30.14 with washout at kτ = 1 for any feed; the optimum sits at C_R*/(n + 1). For a target C_R the tool lays out the smallest plug-flow fermentor the way Fig. 30.2 does: bypass the excess fluid below the maximum-rate point, a recycle loop then plug flow above it.',
    keys: 'product limiting microbial fermentation product poisoning alcohol order of poisoning n batch plug flow fermentor mixed flow fermentor washout optimum residence time maximum production rate fruit fly cocktail recycle bypass side tap smallest volume equation 30.5 30.9 30.10 30.13 30.14 figure 30.1 30.2 30.4 30.5 30.6 example 30.1 problem 30.1 30.2 30.5 30.6 30.7 30.9 30.10',
    units: { t: 'h' },
    inputs: [
      { group: 'Kinetics' },
      { k: 'k', label: 'k', unit: '1/[t]', v: 2, min: 1e-12, max: 1e9, scale: 'log' },
      { k: 'n', label: 'Order of poisoning n', v: 1, min: 0.05, max: 5, step: 0.05 },
      { k: 'CRs', label: 'C_R* (all activity stops)', v: 12, min: 1e-12, max: 1e12 },
      { k: 'yRC', label: 'φ(R/C), product per cell formed', v: 10, min: 1e-9, max: 1e12 },
      { k: 'food', t: 'check', label: 'Food can run short too (Eq. 30.1)', v: false },
      { k: 'CA0', label: 'C_A0', v: 100, min: 1e-12, max: 1e12, if: v => v.food },
      { k: 'CM', label: 'C_M', v: 20, min: 0, max: 1e12, if: v => v.food },
      { k: 'yRA', label: 'φ(R/A)', v: 0.5, min: 1e-12, max: 1e9, if: v => v.food },
      { group: 'Feed' },
      { k: 'CC0', label: 'C_C0 (cells in the feed)', v: 0, min: 0, max: 1e12 },
      { k: 'CR0', label: 'C_R0', v: 0, min: 0, max: 1e12 },
      { group: 'Task' },
      { k: 'task', t: 'seg', label: 'Find', opts: [['mixed', 'mixed fermentor'], ['plug', 'batch or plug flow'], ['best', 'smallest vessel for a target']], v: 'mixed' },
      { k: 'V', label: 'Volume V', v: 1, min: 1e-12, max: 1e12, if: v => v.task === 'mixed' },
      { k: 'vf', label: 'Feed rate v', v: 1, min: 1e-12, max: 1e12, if: v => v.task !== 'plug' },
      { k: 'tp', label: 'Batch time or plug-flow τ_p', unit: '[t]', v: 1, min: 0, max: 1e9, if: v => v.task === 'plug' },
      { k: 'tgt', label: v => (v.task === 'best' ? 'Target C_R' : 'Target C_R (0: none)'), v: 6, min: 0, max: 1e12, if: v => v.task !== 'mixed' }
    ],
    presets: [
      { name: 'Example 30.1 — fruit fly cocktail', v: { k: Math.sqrt(3), n: 1, CRs: 120, yRC: 1, food: false, CC0: 0, CR0: 0, task: 'mixed', V: 30, vf: 25.98, _tu: 'h' }, book: { vopt: { v: 25.98, tol: 0.005 }, FRmax: { v: 1558, tol: 0.005 }, CRopt: { v: 60, tol: 0.005 } }, note: 'k = √3 hr⁻¹, n = 1, C_R* = 0.12 kg alcohol/kg = 120 kg/m³, V_m = 30 m³. Book: C_R,opt = 60 kg/m³ (6% alcohol), v_opt = 25.98 m³/hr, F_R = 1558 kg/hr.' },
      { name: 'Problem 30.9 — winter, n = ½', v: { k: Math.sqrt(3), n: 0.5, CRs: 120, yRC: 1, food: false, CC0: 0, CR0: 0, task: 'mixed', V: 30, vf: 30, _tu: 'h' }, expect: { CRopt: { v: 80, tol: 1e-6 }, vopt: { v: Math.sqrt(3) * 30 * Math.sqrt(1 / 3), tol: 1e-6 }, FRmax: { v: Math.sqrt(3) * 30 * 120 * Math.sqrt(0.5) / Math.pow(1.5, 1.5), tol: 1e-6 } } },
      { name: 'Problem 30.10 — summer, n = 2', v: { k: Math.sqrt(3), n: 2, CRs: 120, yRC: 1, food: false, CC0: 0, CR0: 0, task: 'mixed', V: 30, vf: 23.09, _tu: 'h' }, expect: { CRopt: { v: 40, tol: 1e-6 }, vopt: { v: Math.sqrt(3) * 30 * 4 / 9, tol: 1e-6 }, FRmax: { v: Math.sqrt(3) * 30 * 120 * 4 / 27, tol: 1e-6 } } },
      { name: 'Problems 30.1–30.2 — one tank, v = 1', v: { k: 2, n: 1, CRs: 12, yRC: 10, food: false, CC0: 0, CR0: 0, task: 'mixed', V: 1, vf: 1, _tu: 'h' }, check: { CR: { v: 6, tol: 1e-6 } }, note: 'At v = 4 (Problem 30.2) kτ = 0.5 < 1 and the tank washes out; B-2902 finds the bypass that rescues it.' },
      { name: 'Eq. 30.9 — seed cells and product in the feed', v: { k: 1, n: 1, CRs: 10, yRC: 1, food: false, CC0: 1, CR0: 2, task: 'mixed', V: 1, vf: 2, _tu: 'h' }, expect: { CR: { v: (-9 + Math.sqrt(201)) / 2, tol: 1e-6 } }, note: 'kτ = 0.5 is below 1, yet the cells in the feed keep the tank going.' },
      { name: 'Fig. 30.1 — plug flow with seed cells (Eq. 30.5)', v: { k: 1, n: 1, CRs: 10, yRC: 1, food: false, CC0: 0.5, CR0: 0, task: 'plug', tp: 2, tgt: 8, _tu: 'h' }, expect: { CR: { v: (function () { const Q = Math.exp(2 * 1.05) * 0.05; return (-0.5 + 10 * Q) / (1 + Q); })(), tol: 1e-6 }, tT: { v: 10 / 10.5 * Math.log(85), tol: 1e-6 } } },
      { name: 'Problem 30.5 — C_R = 6 at v = 3', v: { k: 2, n: 1, CRs: 12, yRC: 10, food: false, CC0: 0, CR0: 0, task: 'best', vf: 3, tgt: 6, _tu: 'h' }, check: { V: { v: 3, tol: 1e-6 } }, note: 'C_R*/2 is the maximum-rate point: one recycle loop of V = 3, nothing else.' },
      { name: 'Problem 30.6 — C_R = 4: bypass the excess', v: { k: 2, n: 1, CRs: 12, yRC: 10, food: false, CC0: 0, CR0: 0, task: 'best', vf: 3, tgt: 4, _tu: 'h' }, check: { V: { v: 2, tol: 1e-6 } } },
      { name: 'Problem 30.7 — C_R = 9: loop, then plug flow', v: { k: 2, n: 1, CRs: 12, yRC: 10, food: false, CC0: 0, CR0: 0, task: 'best', vf: 3, tgt: 9, _tu: 'h' }, check: { V: { v: 3 + 1.5 * Math.log(3), tol: 1e-6 } } },
      { name: 'Food runs short as well', v: { k: 1, n: 1, CRs: 40, yRC: 5, food: true, CA0: 100, CM: 20, yRA: 0.5, CC0: 0, CR0: 0, task: 'mixed', V: 10, vf: 4, _tu: 'h' }, note: 'With C_A0 = 100 the product can reach at most φ(R/A)C_A0 = 50 > C_R* = 40, but the Monod factor slows things long before.' }
    ],
    eqs: [
      ['r_C = k\\left(1-\\frac{C_R}{C_R^*}\\right)^n C_C', 'Eqs. 30.2, 30.12', 'Food in excess.'],
      ['r_R = \\varphi(R/C)\\,r_C = k\\left(1-\\frac{C_R}{C_R^*}\\right)^n\\left(C_R - C_{R0} + \\varphi(R/C)C_{C0}\\right)', 'Eq. 30.3'],
      ['C_{R,max\\ rate} = \\tfrac{1}{2}\\left(C_{R0} + C_R^* - \\varphi(R/C)C_{C0}\\right)\\quad (n=1)', 'Eq. 30.4'],
      ['k\\tau_p = \\frac{C_R^*}{C_R^*-C_{R0}+\\varphi(R/C)C_{C0}}\\ln\\frac{C_C(C_R^*-C_{R0})}{C_{C0}(C_R^*-C_R)}', 'Eq. 30.5'],
      ['k\\tau_m = \\frac{C_R - C_{R0}}{\\left(C_R - C_{R0} + \\varphi(R/C)C_{C0}\\right)\\left(1-C_R/C_R^*\\right)^n}', 'Eqs. 30.9, 30.13'],
      ['k\\tau_m = \\frac{1}{(1-C_R/C_R^*)^n}\\quad(k\\tau_m>1)', 'Eqs. 30.10, 30.14', 'C_C0 = C_R0 = 0.'],
      ['C_{R,opt} = \\frac{C_R^*}{n+1},\\qquad v_{opt} = kV\\left(\\frac{n}{n+1}\\right)^n', 'Fig. 30.6'],
      ['F_{R,max} = \\frac{kVC_R^*\\,n^n}{(n+1)^{n+1}}\\quad\\left(= kVC_R^*/4\\ \\text{for}\\ n=1\\right)', 'Fig. 30.6']
    ],
    method: ['The same engine as B-2901, written in terms of C_R: with food in excess the rate needs only C_R, and the closed forms of the chapter (Eqs. 30.4, 30.5 and 30.13) are shown beside the numerical answers as a check. The production rate of a mixed fermentor per unit volume is r_R at its exit, so the optimum always sits at the composition of maximum rate; the smallest vessel for a target follows Fig. 30.2.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'h';
      if (+v.CR0 >= +v.CRs) throw new Error('C_R0 must be below C_R*.');
      const M = pmodel(v);
      const k = +v.k, C = +v.CRs, nn = +v.n;
      const rR = d => M.rC(d) * M.yRA / M.yCA;
      const dO = M.dMax(), tauO = dO * M.yCA / M.rC(dO);
      const plain = !v.food;
      const CRof = d => M.CR(d);
      let shade = null, dMark, diagram;
      if (v.task === 'mixed') {
        const V = +v.V, vf = +v.vf, tau = V / vf;
        const d = M.mixed(tau, 0);
        const wash = d <= 1e-12 * M.Dmax && !(M.CC0 > 0);
        kpi.push({ k: 'CR', l: 'C_R leaving', v: M.CR(d), main: true, s: wash ? 'washout: kτ = ' + fmt(k * tau, 4) + ' is too small' : 'kτ = ' + fmt(k * tau, 4) });
        kpi.push({ k: 'CC', l: 'C_C leaving', v: M.CC(d) });
        if (v.food) kpi.push({ k: 'CA', l: 'C_A leaving', v: M.CA(d) });
        kpi.push({ k: 'FR', l: 'Product made, v(C_R − C_R0)', v: vf * (M.CR(d) - M.CR0), u: 'per ' + tu });
        if (plain && !wash) kpi.push({ k: 'ktEq', l: 'Eq. 30.13 at this C_R gives kτ', v: PL.ktMix(v, M.CR(d)) });
        if (!(M.CC0 > 0)) {
          const kw = M.rC(1e-9 * M.Dmax) / (M.yCA * 1e-9 * M.Dmax);
          kpi.push({ k: 'tauW', l: 'Washout below τ', v: 1 / kw, u: tu, s: plain && M.CR0 === 0 ? 'kτ = 1 whatever the feed' : 'kτ = ' + fmt(k / kw, 4) });
        }
        kpi.push({ k: 'CRopt', l: 'C_R at the optimum', v: M.CR(dO), s: plain ? 'C_R*/(n + 1) = ' + fmt(C / (nn + 1), 4) + (M.CC0 > 0 || M.CR0 > 0 ? ' for a clean feed' : '') : '' });
        kpi.push({ k: 'tauopt', l: 'Best τ', v: tauO, u: tu });
        kpi.push({ k: 'vopt', l: 'Best feed rate v_opt = V/τ_opt', v: V / tauO, u: 'per ' + tu });
        kpi.push({ k: 'FRmax', l: 'Most product per unit time', v: V * rR(dO), u: 'per ' + tu });
        shade = { kind: 'rect', d0: 0, d1: d };
        dMark = d;
        plots.push(mixedChart(M, v, tau, tu));
        info.push('Production per unit volume equals r_R at the exit, so the best tank sits at the maximum-rate composition: the tangent from the origin on the C_R–kτ curve (Fig. 30.4).');
        diagram = BIO.sketch({ feed: 'v = ' + fmt(vf, 3), first: { kind: 'tanks', m: 1, tag: 'F', fill: d / M.Dmax }, bypass: null, series: [], out: 'C_R = ' + fmt(M.CR(d), 3) });
      } else if (v.task === 'plug') {
        const tp = +v.tp;
        if (!(M.CC0 > 0)) warn.push('No cells in the feed: a batch or a plug-flow vessel never starts. Seed it (C_C0 > 0), or use mixed flow or recycle.');
        const d = M.plug(tp, 0);
        kpi.push({ k: 'CR', l: 'C_R after τ_p', v: M.CR(d), main: true, s: 'kτ_p = ' + fmt(k * tp, 4) });
        kpi.push({ k: 'CC', l: 'C_C', v: M.CC(d) });
        if (v.food) kpi.push({ k: 'CA', l: 'C_A', v: M.CA(d) });
        if (plain && nn === 1 && M.CC0 > 0) kpi.push({ k: 'CReq', l: 'Eq. 30.5 gives', v: PL.plug1(v, k * tp) });
        const tg = +v.tgt;
        if (tg > M.CR0 && M.CC0 > 0) {
          const dT = (tg - M.CR0) / M.yRA;
          if (dT < M.Dmax) {
            kpi.push({ k: 'tT', l: 'Time to reach C_R = ' + fmt(tg, 4), v: M.tPlug(0, dT), u: tu, s: plain && nn === 1 ? 'Eq. 30.5: ' + fmt(PL.tplug1(v, tg), 5) : '' });
          } else warn.push('C_R = ' + fmt(tg, 4) + ' is beyond reach: the product stops at ' + fmt(M.CR(M.Dmax), 4) + '.');
        }
        shade = { kind: 'area', d0: 0, d1: d };
        dMark = d;
        plots.push(plugChart(M, v, tp, tu));
        diagram = BIO.sketch({ feed: 'C_C0 = ' + fmt(M.CC0, 3), first: { kind: 'tube', tag: 'kτ_p ' + fmt(k * tp, 3), prof: N.linspace(0, tp, 12).map(t => M.plug(t, 0) / M.Dmax) }, bypass: null, series: [], out: 'C_R = ' + fmt(M.CR(d), 3) });
      } else {
        const vf = +v.vf, tg = +v.tgt;
        const dT = (tg - M.CR0) / M.yRA;
        if (!(dT > 0 && dT < M.Dmax)) throw new Error('The target must lie between C_R0 and ' + fmt(M.CR(M.Dmax), 4) + '.');
        let V, how, sk;
        if (dO <= 1e-9 * M.Dmax) {
          V = vf * M.tPlug(0, dT); how = 'plug flow all the way: the rate only falls from the feed on';
          sk = { first: { kind: 'tube', tag: 'V = ' + fmt(V, 3) } };
        } else if (dT <= dO) {
          V = vf * (dT / dO) * tauO; how = 'one loop at R = ∞ kept at the maximum rate, C_R = ' + fmt(M.CR(dO), 4) + ', on ' + fmt(dT / dO * vf, 4) + ' of the feed; bypass the rest';
          sk = { first: { kind: 'loop', tag: 'V = ' + fmt(V, 3), flow: fmt(dT / dO * vf, 3) + ' in' }, bypass: 'bypass ' + fmt((1 - dT / dO) * vf, 3) };
        } else {
          const V1 = vf * tauO, V2 = vf * M.tPlug(dO, dT);
          V = V1 + V2; how = 'a loop at R = ∞ to C_R = ' + fmt(M.CR(dO), 4) + ' (V = ' + fmt(V1, 4) + '), then plug flow (V = ' + fmt(V2, 4) + ')';
          sk = { first: { kind: 'loop', tag: 'loop ' + fmt(V1, 3) }, series: [{ kind: 'tube', tag: 'plug ' + fmt(V2, 3) }] };
        }
        kpi.push({ k: 'V', l: 'Smallest vessel', v: V, main: true, s: how });
        const Vm = vf * dT * M.yCA / M.rC(dT);
        const Vp = M.CC0 > 0 ? vf * M.tPlug(0, dT) : Infinity;
        kpi.push({ k: 'Vm', l: 'One mixed tank instead', v: Vm, s: fmt(Vm / V, 3) + '× larger' });
        tables.push({ title: 'Ways to reach C_R = ' + fmt(tg, 4) + ' at v = ' + fmt(vf, 4), cols: ['arrangement', 'volume', 'relative'], rows: [['best (Fig. 30.2)', V, 1], ['one mixed tank, no bypass', Vm, Vm / V], ['plain plug flow', isFinite(Vp) ? Vp : 'never starts', isFinite(Vp) ? Vp / V : '—']] });
        shade = { kind: 'best', dO, dT };
        dMark = dT;
        diagram = BIO.sketch(Object.assign(sk, { feed: 'v = ' + fmt(vf, 3), out: 'C_R = ' + fmt(tg, 3) }));
      }
      /* 1/r_R against C_R with the area that is τ */
      const dd = N.linspace(M.Dmax * 1e-3, M.Dmax * 0.995, 260);
      const inv = dd.map(d => 1 / rR(d));
      const ymin = Math.min(...inv), yTop = 7 * ymin;
      const ser = [{ name: '1/r_R', x: dd.map(CRof), y: inv, color: 1, w: 2.4 }];
      const marks = [{ t: 'v', x: M.CR(dO), color: 2, label: plain && nn === 1 && !(M.CC0 > 0) && M.CR0 === 0 ? 'C_R*/2' : 'max rate' }];
      if (shade) {
        if (shade.kind === 'rect') marks.push({ t: 'rect', x0: M.CR(0), x1: M.CR(shade.d1), y0: 0, y1: Math.min(1 / rR(shade.d1), yTop), color: 3, fill: 0.2, label: 'τ_m' });
        if (shade.kind === 'area') { const z = N.linspace(0, shade.d1, 120); ser.push({ name: 'area = τ_p', x: z.map(CRof), y: z.map(d => Math.min(1 / rR(d), 1e6 * yTop)), type: 'area', color: 3, fill: 0.22, legend: false }); }
        if (shade.kind === 'best') {
          marks.push({ t: 'rect', x0: M.CR(0), x1: M.CR(shade.dO), y0: 0, y1: Math.min(1 / rR(dO), yTop), color: 3, fill: 0.2, label: 'loop' });
          if (shade.dT > dO) { const z = N.linspace(dO, shade.dT, 120); ser.push({ name: 'plug section', x: z.map(CRof), y: z.map(d => 1 / rR(d)), type: 'area', color: 4, fill: 0.22, legend: false }); }
        }
      }
      if (dMark !== undefined && dMark > 0) marks.push({ t: 'v', x: M.CR(dMark), color: 8, label: v.task === 'best' ? 'target' : 'leaving', dy: 15 });
      plots.unshift({ title: '1/r_R against C_R: areas are space-times (Figs. 30.1, 30.2)', x: { label: 'C_R', min: M.CR0, max: M.CR(M.Dmax) }, y: { label: '1/r_R', unit: tu, min: 0, max: yTop }, series: ser, marks, cap: 'Below the maximum-rate point keep the fluid there with a loop and bypass the excess; above it, reach the point in one step and continue in plug flow.' });
      const q = N.linspace(0, 1, 101);
      plots.push({ title: 'Slowdown with the order of poisoning (Fig. 30.5)', x: { label: 'C_R/C_R*', min: 0, max: 1 }, y: { label: 'k_obs/k', min: 0, max: 1.05 },
        series: [{ name: 'n = 0', x: [0, 0.999, 1], y: [1, 1, 0], color: 'var(--muted)', w: 1.4 }].concat([0.5, 1, 2].map((m2, i) => ({ name: 'n = ' + m2, x: q, y: q.map(z => Math.pow(1 - z, m2)), color: i + 2, w: 1.5 }))).concat([{ name: 'this n = ' + fmt(nn, 3), x: q, y: q.map(z => Math.pow(1 - z, nn)), color: 1, w: 2.6 }]) });
      return { kpi, warn, info, plots, tables, diagram };
    }
  });
  function mixedChart(M, v, tau, tu) {
    const k = +v.k;
    const kw = M.CC0 > 0 ? 0 : 1 / (M.yCA * 1e-9 * M.Dmax / M.rC(1e-9 * M.Dmax));
    const dO = M.dMax(), tauO = dO * M.yCA / M.rC(dO);
    const top = Math.max(4 * tauO, 1.4 * tau) * k;
    const kt = N.linspace(top / 500, top, 260);
    const ds = kt.map(z => M.mixed(z / k, 0));
    const CRo = M.CR(dO), slope = (CRo - M.CR0) / (k * tauO);
    const marks = [{ t: 'h', y: M.CR(M.Dmax), color: 'var(--muted)', label: 'C_R*' }, { t: 'seg', x0: 0, y0: M.CR0, x1: Math.min(top, (M.CR(M.Dmax) - M.CR0) / slope), y1: M.CR0 + slope * Math.min(top, (M.CR(M.Dmax) - M.CR0) / slope), color: 2, dash: '5 4' }, { t: 'pt', x: k * tauO, y: CRo, color: 2, label: 'optimum' }, { t: 'v', x: k * tau, color: 8, label: 'this tank', dy: 15 }];
    if (kw > 0) marks.push({ t: 'v', x: k / kw, color: 'var(--muted)', label: 'washout' });
    const ser = [{ name: 'C_R', x: kt, y: ds.map(M.CR), color: 1, w: 2.4 }, { name: 'C_C', x: kt, y: ds.map(M.CC), color: 3, w: 1.8 }];
    return { title: 'Mixed fermentor: C_R against kτ_m (Figs. 30.4, 30.6)', x: { label: 'kτ_m', min: 0, max: top }, y: { label: 'concentration', min: 0, max: 1.12 * M.CR(M.Dmax) }, series: ser, marks, cap: 'The dashed tangent from the origin touches the curve at the highest production rate v(C_R − C_R0).' };
  }
  function plugChart(M, v, tp, tu) {
    const k = +v.k, phi = M.yRA / M.yCA;
    const seeds = M.CC0 > 0 ? [M.CC0 / 4, M.CC0, M.CC0 * 4] : [0.01, 0.05, 0.2].map(f => f * M.CRs / phi);
    const T = Math.max(3 * tp, 1e-9);
    const ser = seeds.map((c0, i) => {
      const mm = BIO.model({ k: M.k, CM: M.CM, CRs: M.CRs, n: M.n, yCA: M.yCA, yRA: M.yRA, CA0: M.CA0, CC0: c0, CR0: M.CR0 });
      const p = BIO.batchPath(mm, T, 180, 0);
      return { name: 'φ(R/C)C_C0 = ' + fmt(phi * c0, 3), x: p.t.map(t => k * t), y: p.d.map(mm.CR), color: Math.abs(c0 - M.CC0) < 1e-12 * Math.max(c0, 1) ? 1 : i + 2, w: Math.abs(c0 - M.CC0) < 1e-12 * Math.max(c0, 1) ? 2.6 : 1.6 };
    });
    return { title: 'Batch or plug flow: C_R against kτ_p for several seedings (Fig. 30.1)', x: { label: 'kτ_p', min: 0, max: k * T }, y: { label: 'C_R', min: 0, max: M.CR(M.Dmax) * 1.05 }, series: ser, marks: [{ t: 'h', y: M.CR(M.Dmax), color: 'var(--muted)', label: 'C_R*' }, { t: 'v', x: k * tp, color: 8, label: 'this τ_p' }] };
  }
})(typeof window !== 'undefined' ? window : globalThis);
