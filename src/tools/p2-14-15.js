/* Chapters 14 and 15 — the tanks-in-series model; the convection model for laminar flow */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, RT = ST.RT, KIT = ST.KIT, K = ST.K, D = ST.D, RTDI = ST.RTDI;

  /* ===================================================================== */
  /* F-1401  tanks in series                                                */
  /* ===================================================================== */
  const Emax = Nk => (Nk <= 1 ? 1 : Math.exp(Math.log(Nk) + (Nk - 1) * Math.log(Nk - 1) - (Nk - 1) - N.lgamma(Nk)));
  ST.tool({
    id: 'tanks', tag: 'F-1401', ch: 14, sec: 'Ch. 14', title: 'Tanks in series: RTD, fitting N, recirculation, and conversion for any kinetics',
    lede: 'The MacMullin–Weber curves for any N, including non-integer N; N from a variance, from a peak’s width or from data; what adding tanks does to the peak and the spread; tracer in closed recirculating loops with or without throughflow; the fraction of material staying longer than any time; and conversion as a microfluid tank by tank or as a macrofluid, against plug flow and the dispersion model of equal variance.',
    keys: 'tanks in series model macmullin weber gamma erlang N number of tanks variance recirculation closed loop throughflow peak width dollar bills problem 14.2 14.3 14.6 example 14.3 conversion micro macro',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'Task' },
      { k: 'what', t: 'select', label: 'Show', opts: [['curves', 'E and F curves for several N'], ['fit', 'N from measurements'], ['add', 'adding or removing tanks'], ['recirc', 'recirculating systems'], ['frac', 'fraction staying longer than t'], ['react', 'conversion in N tanks']], v: 'curves' },
      { k: 'nlist', t: 'text', label: 'N values to draw', v: '1, 2, 3, 5, 10, 30', mono: true, if: v => v.what === 'curves' },
      { k: 'fitHow', t: 'select', label: 'From', opts: [['var', 't̄ and σ²'], ['peak', 'a peak’s position and width (Fig. 14.3)'], ['data', 'pulse data (table)']], v: 'peak', if: v => v.what === 'fit' },
      { k: 'tb', label: 'Mean time t̄', unit: '[t]', v: 1, min: 1e-12, max: 1e12, if: v => (v.what === 'fit' && v.fitHow === 'var') || v.what === 'frac' || v.what === 'react' },
      { k: 's2', label: 'Variance σ²', unit: '[t]²', v: 0.25, min: 1e-15, max: 1e15, if: v => v.what === 'fit' && v.fitHow === 'var' },
      { k: 'thMax', label: 'Peak position θ_max', v: 1 / 6, min: 1e-6, max: 10, if: v => v.what === 'fit' && v.fitHow === 'peak', hint: 'In any consistent time unit; only the ratio matters.' },
      { k: 'dTh', label: 'Width at the inflection points Δθ', v: 0.1, min: 1e-9, max: 10, if: v => v.what === 'fit' && v.fitHow === 'peak' },
      { k: 'data', t: 'table', label: 'Pulse response (t, C)', cols: ['t', 'C'], v: [[1, 0], [3, 0], [5, 10], [7, 10], [9, 10], [11, 10], [13, 0], [15, 0]], if: v => v.what === 'fit' && v.fitHow === 'data' },
      { k: 'N1', label: 'Tanks now', v: 10, min: 1, max: 1000, step: 1, if: v => v.what === 'add' },
      { k: 'N2', label: 'Tanks after the change', v: 20, min: 1, max: 1000, step: 1, if: v => v.what === 'add' },
      { k: 'Cm1', label: 'Peak concentration now', v: 100, min: 0, max: 1e12, if: v => v.what === 'add' },
      { k: 'sp1', label: 'Spread now', v: 1, min: 0, max: 1e12, if: v => v.what === 'add' },
      { k: 'Nr', label: 'Tanks per circuit', v: 5, min: 1, max: 500, step: 1, if: v => v.what === 'recirc' },
      { k: 'tc', label: 'Mean circulation time', unit: '[t]', v: 1, min: 1e-9, max: 1e9, if: v => v.what === 'recirc' },
      { k: 'tsys', label: 'Mean residence time of the whole system (0 = closed)', unit: '[t]', v: 0, min: 0, max: 1e12, if: v => v.what === 'recirc' },
      { k: 'Nf', label: 'N', v: 2, min: 0.2, max: 1000, scale: 'log', if: v => v.what === 'frac' || v.what === 'react' },
      { k: 'tq', label: 'Time t', unit: '[t]', v: 2, min: 0, max: 1e12, if: v => v.what === 'frac' },
      { k: 'pop', label: 'Total population (optional)', v: 0, min: 0, max: 1e18, if: v => v.what === 'frac' }
    ].concat(KIT.kin({ gas: false, d: { law: 'power', k: 1, a: 2, CA0: 1 } }).map(f => Object.assign({}, f, { if: f.group ? (vv => vv.what === 'react') : f.if ? (vv => vv.what === 'react' && f.if(vv)) : (vv => vv.what === 'react') }))),
    presets: [
      { name: 'RTD curves (Fig. 14.2)', v: { what: 'curves', nlist: '1, 2, 3, 5, 10, 30' } },
      { name: 'Example 14.3 (N from a peak)', v: { what: 'fit', fitHow: 'peak', thMax: 1 / 6, dTh: 0.1 }, book: { N: { v: 12, tol: 0.02 } }, note: 'The first peak sits at θ = 1/6 with width 1/10 at its inflection points: Δθ/θ_max = 2/√(N−1) gives N = 12.' },
      { name: 'Problem 14.8 (flat response)', v: { what: 'fit', fitHow: 'data', data: [[1, 0], [3, 0], [5, 10], [7, 10], [9, 10], [11, 10], [13, 0], [15, 0]] }, note: 'A flat-topped response: how many tanks?' },
      { name: 'Problem 14.2 (ten more tanks)', v: { what: 'add', N1: 10, N2: 20, Cm1: 100, sp1: 1 }, note: 'Ten well-behaved tanks give a peak of 100 mmol and a spread of 1 min. Add ten more: peak and spread?' },
      { name: 'Recirculation (Fig. 14.5)', v: { what: 'recirc', Nr: 5, tc: 1, tsys: 0 }, note: 'Five tanks per circuit in a closed loop: each pass broadens the peaks until the tracer is uniform.' },
      { name: 'Recirculation with throughflow (Fig. 14.6)', v: { what: 'recirc', Nr: 5, tc: 1, tsys: 4 }, note: 'Slow throughflow: the circulation pattern rides on the exponential decay of one big stirred tank.' },
      { name: 'Problem 14.3 (dollar bills)', v: { what: 'frac', Nf: 1, tb: 1.25, tq: 21, pop: 1.25e9 }, book: { left: { v: 63, tol: 0.03 } }, note: '1.25 billion bills in circulation, 1 billion replaced a year, withdrawn at random: a mixed tank with t̄ = 1.25 years. Old bills left after 21 years?' },
      { name: 'Problem 14.6 (particles, two tanks)', v: { what: 'frac', Nf: 2, tb: 2, tq: 2, pop: 0 }, book: { frac: { v: 0.406, tol: 0.01 } }, note: 'Particles need 2 min to convert fully; two 1 m³ tanks in series at 1 m³/min: the fraction staying at least 2 min.' },
      { name: 'Second-order reaction in N tanks', v: { what: 'react', Nf: 4, tb: 10, law: 'power', k: 1, a: 2, ob: 0, CA0: 1 }, note: 'kC_A0τ = 10 in four tanks: microfluid tank by tank, macrofluid through the Erlang RTD, and the equivalents.' }
    ],
    eqs: [
      ['E_\\theta = N\\frac{(N\\theta)^{N-1}}{(N-1)!}e^{-N\\theta},\\qquad \\sigma^2_\\theta = \\frac{1}{N}', 'Eq. 14.3'],
      ['\\theta_{max} = \\frac{N-1}{N},\\qquad \\frac{\\Delta\\theta}{\\theta_{max}} = \\frac{2}{\\sqrt{N-1}}', 'Fig. 14.3'],
      ['F = 1 - e^{-N\\theta}\\left[1 + N\\theta + \\frac{(N\\theta)^2}{2!} + \\cdots + \\frac{(N\\theta)^{N-1}}{(N-1)!}\\right]', 'Eq. 14.8'],
      ['C_{\\theta,pulse} = Ne^{-N\\theta}\\sum_{m=1}^{\\infty}\\frac{(N\\theta)^{mN-1}}{(mN-1)!}', 'Eq. 14.5', 'Closed recirculation.'],
      ['\\frac{C_A}{C_{A0}} = \\frac{1}{(1+k\\tau_i)^N},\\qquad \\frac{C_{N\\,tanks}}{C_{plug}} \\approx 1 + \\frac{(k\\tau)^2}{2N}', 'Eqs. 14.9–14.11'],
      ['\\frac{\\bar C_A}{C_{A0}} = \\int_0^\\infty \\left(\\frac{C_A}{C_{A0}}\\right)_{batch}\\frac{t^{N-1}}{\\bar t^{N}}\\frac{N^N}{(N-1)!}e^{-tN/\\bar t}\\,dt', 'Eq. 14.15', 'Macrofluid.']
    ],
    method: ['Non-integer N uses the gamma distribution, which reduces to the Erlang form for whole N; F comes from the regularized incomplete gamma function.', 'Microfluid conversion solves the mixed-flow balance tank by tank (whole N); the macrofluid result integrates the batch history over the gamma RTD.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'min';
      if (v.what === 'curves') {
        const ns = String(v.nlist).split(/[,;\s]+/).map(parseFloat).filter(x => x > 0).slice(0, 8);
        const th = N.linspace(0, 3, 400);
        plots.push({ title: 'E_θ for N tanks in series (Fig. 14.2)', x: { label: 'θ = t/t̄', min: 0, max: 3 }, y: { label: 'E_θ', min: 0 }, series: ns.map((n, i) => ({ name: 'N = ' + N.fmt(n, 3), x: th, y: th.map(q => RT.Etanks(q, 1, n)), color: i + 1 })), hmax: 480 });
        plots.push({ title: 'F curves (Fig. 14.7)', x: { label: 'θ', min: 0, max: 3 }, y: { label: 'F', min: 0, max: 1 }, series: ns.map((n, i) => ({ name: 'N = ' + N.fmt(n, 3), x: th, y: th.map(q => RT.Ftanks(q, 1, n)), color: i + 1 })) });
        tables.push({ title: 'Properties (Fig. 14.3)', cols: ['N', 'σ²_θ', 'θ_max', 'E_θ,max', 'D/uL with the same variance'], rows: ns.map(n => [n, 1 / n, n > 1 ? (n - 1) / n : 0, Emax(n), RTDI.dFromVar(1 / n, 'closed')]) });
        kpi.push({ k: 'n', l: 'Curves', v: ns.length });
        return { kpi, warn, info, plots, tables };
      }
      if (v.what === 'fit') {
        let Nq, tb = NaN;
        if (v.fitHow === 'var') { tb = +v.tb; Nq = tb * tb / +v.s2; }
        else if (v.fitHow === 'peak') { const r = 2 * +v.thMax / +v.dTh; Nq = 1 + r * r; }
        else { const R = RTDI.pulse(RTDI.rows(v.data, 'readings'), { rule: 'pl' }); tb = R.mean; Nq = R.mean * R.mean / R.var; kpi.push({ k: 'tbar', l: 't̄ of the data', v: R.mean }, { k: 'var', l: 'σ² of the data', v: R.var }); const th = N.linspace(0, 3, 300); plots.push({ title: 'Data against N tanks of the same variance', x: { label: 'θ', min: 0 }, y: { label: 'E_θ', min: 0 }, series: [{ name: 'data', x: R.t.map(t => t / R.mean), y: R.E.map(e => e * R.mean), type: 'pts', color: 1 }, { name: 'N = ' + N.fmt(Nq, 3), x: th, y: th.map(q => RT.Etanks(q, 1, Nq)), color: 2 }] }); }
        kpi.push({ k: 'N', l: 'Number of tanks N', v: Nq, main: true, s: 'nearest whole number ' + Math.max(1, Math.round(Nq)) });
        kpi.push({ k: 'd', l: 'Equivalent D/uL (closed vessel)', v: RTDI.dFromVar(1 / Nq, 'closed') });
        if (Nq < 1) warn.push('N < 1: the curve is wider than one mixed tank — look for bypassing or stagnant regions (Ch. 12) instead.');
        return { kpi, warn, info, plots, tables };
      }
      if (v.what === 'add') {
        const N1 = +v.N1, N2 = +v.N2;
        /* same tank size: t̄ ∝ N; peak height of C ∝ E_max(N)/t̄; spread ∝ √N */
        const peak2 = +v.Cm1 * (Emax(N2) / N2) / (Emax(N1) / N1), sp2 = +v.sp1 * Math.sqrt(N2 / N1);
        kpi.push({ k: 'Cm2', l: 'Peak concentration after', v: peak2, main: true, s: 'large-N estimate ×√(N₁/N₂) = ' + N.fmt(+v.Cm1 * Math.sqrt(N1 / N2), 4) }, { k: 'sp2', l: 'Spread after', v: sp2, main: true });
        kpi.push({ k: 'rel', l: 'Relative spread σ/t̄ after / before', v: Math.sqrt(N1 / N2), s: 'the relative spread falls as 1/√N' });
        const t = N.linspace(0, 2.5 * Math.max(N1, N2), 400);
        plots.push({ title: 'Output for the same pulse, same tank size', x: { label: 't/τ_i' }, y: { label: 'C (relative)', min: 0 }, series: [{ name: N1 + ' tanks', x: t, y: t.map(q => RT.Etanks(q, N1, N1) * +v.Cm1 / (Emax(N1) / N1)), color: 1 }, { name: N2 + ' tanks', x: t, y: t.map(q => RT.Etanks(q, N2, N2) * +v.Cm1 / (Emax(N1) / N1)), color: 2 }] });
        return { kpi, warn, info, plots, tables };
      }
      if (v.what === 'recirc') {
        const Nr = Math.round(+v.Nr), tc = +v.tc, ts = +v.tsys;
        const tmax = tc * 6, t = N.linspace(1e-9, tmax, 600);
        const Ci = q => { let s = 0; for (let mm = 1; mm < 400; mm++) { const term = RT.Etanks(q, mm * tc, mm * Nr); s += term; if (mm * tc > q * 1.5 + 5 * tc) break; } return s * tc; };
        const y = t.map(q => Ci(q) * (ts > 0 ? Math.exp(-q / ts) : 1));
        plots.push({ title: ts > 0 ? 'Recirculation with throughflow (Fig. 14.6)' : 'Closed recirculating system (Fig. 14.5)', x: { label: 't/t̄_circuit', min: 0 }, y: { label: 'C/C_∞', min: 0 }, series: [{ name: 'recorder signal', x: t.map(q => q / tc), y, color: 1 }].concat(ts > 0 ? [{ name: 'well-mixed system, e^{−t/t̄}', x: t.map(q => q / tc), y: t.map(q => Math.exp(-q / ts)), color: 8, dash: '5 4' }] : [{ name: 'uniform', x: [0, 6], y: [1, 1], color: 8, dash: '5 4' }]), cap: 'C_∞ is the concentration once the tracer is spread evenly through the system.' });
        const pk = []; for (let i = 1; i < y.length - 1; i++) if (y[i] > y[i - 1] && y[i] >= y[i + 1]) pk.push([t[i] / tc, y[i]]);
        kpi.push({ k: 'peaks', l: 'Visible peaks', v: pk.length, s: pk.slice(0, 4).map(p => 'θ ' + N.fmt(p[0], 3)).join(', ') });
        tables.push({ title: 'Peaks', cols: ['t/t̄_circuit', 'C/C_∞'], rows: pk });
        return { kpi, warn, info, plots, tables };
      }
      if (v.what === 'frac') {
        const Nq = +v.Nf, tb = +v.tb, tq = +v.tq;
        const F = RT.Ftanks(tq, tb, Nq);
        kpi.push({ k: 'frac', l: 'Fraction staying longer than t', v: 1 - F, main: true }, { k: 'fracS', l: 'Fraction leaving before t', v: F });
        if (+v.pop > 0) kpi.push({ k: 'left', l: 'Of the population, still inside after t', v: +v.pop * (1 - F), main: true });
        const t = N.linspace(0, Math.max(tq * 1.5, tb * 4), 300);
        plots.push({ title: 'Fraction still inside, 1 − F', x: { label: 't', min: 0 }, y: { label: '1 − F', min: 0, max: 1 }, series: [{ name: 'N = ' + N.fmt(Nq, 3), x: t, y: t.map(q => 1 - RT.Ftanks(q, tb, Nq)), color: 1 }], marks: [{ t: 'v', x: tq, color: 2, label: 't' }] });
        return { kpi, warn, info, plots, tables };
      }
      /* react */
      const m = KIT.model(v), tb = +v.tb, Nq = +v.Nf, nI = Math.max(1, Math.round(Nq));
      let X = 0; const stages = [];
      for (let i = 0; i < nI; i++) { const Xo = K.XMFR(m, tb / nI, X); stages.push([i + 1, X, Xo]); X = Xo; }
      const Xb = RT.batchX(m, tb * (1 + 30 / Math.sqrt(Nq)));
      const macro = N.quad(t => Xb(t) * RT.Etanks(t, tb, Nq), 0, tb * (1 + 30 / Math.sqrt(Nq)), { rtol: 1e-9 });
      const Xp = K.XPFR(m, tb).X, Xm = K.XMFR(m, tb);
      const dEq = RTDI.dFromVar(1 / Nq, 'closed');
      const Xd = isFinite(dEq) ? RT.dispersionX(m, tb, dEq).X : NaN;
      kpi.push({ k: 'Xmicro', l: nI + ' tanks, microfluid', v: X, main: true }, { k: 'Xmacro', l: 'Macrofluid through the same RTD', v: macro }, { k: 'Xp', l: 'Plug flow', v: Xp }, { k: 'Xm', l: 'One mixed tank', v: Xm }, { k: 'Xd', l: 'Dispersion model, D/uL = ' + N.fmt(dEq, 3), v: Xd });
      if (Math.abs(nI - Nq) > 1e-9) info.push('Tank-by-tank results use N = ' + nI + '; the macrofluid integral uses N = ' + N.fmt(Nq, 4) + ' exactly.');
      const Ns = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 50];
      const xN = Ns.map(n => { let x = 0; for (let i = 0; i < n; i++) x = K.XMFR(m, tb / n, x); return x; });
      plots.push({ title: 'Conversion against the number of tanks (same total τ)', x: { label: 'N', log: true, min: 1, max: 50 }, y: { label: 'X_A', min: 0, max: 1 }, series: [{ name: 'N tanks (microfluid)', x: Ns, y: xN, color: 1, marker: 'o' }, { name: 'plug flow', x: [1, 50], y: [Xp, Xp], color: 8, dash: '5 4' }] });
      tables.push({ title: 'Tank by tank', cols: ['tank', 'X_A in', 'X_A out'], rows: stages });
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* F-1501  laminar flow: the convection model                             */
  /* ===================================================================== */
  /* RTD from a velocity profile u(ξ) over ξ ∈ [0,1] with flux weight w(ξ) (pipe: 2ξ, slit/film: 1) */
  function profileRTD(kind, n) {
    /* returns {u(ξ) normalized to mean 1, w(ξ), umax} */
    if (kind === 'pipe') return { u: x => 2 * (1 - x * x), w: x => 2 * x, label: 'Newtonian, circular pipe' };
    if (kind === 'power') { const a = (3 * n + 1) / (n + 1), p = (n + 1) / n; return { u: x => a * (1 - Math.pow(x, p)), w: x => 2 * x, label: 'power-law fluid, n = ' + N.fmt(n, 3) }; }
    if (kind === 'slit') return { u: x => 1.5 * (1 - x * x), w: () => 1, label: 'between parallel plates' };
    return { u: x => 1.5 * (1 - x * x), w: () => 1, label: 'falling film' };
  }
  /* E_θ by binning flux-weighted residence times θ = 1/u */
  function Etheta(prof, thG) {
    const M = 40000, ths = [], ws = [];
    for (let i = 0; i < M; i++) { const x = (i + 0.5) / M, u = prof.u(x); if (u <= 0) continue; ths.push(1 / u); ws.push(prof.w(x) * u / M); }
    const tot = N.sum(ws);
    const E = new Array(thG.length).fill(0), dth = thG[1] - thG[0];
    ths.forEach((t, i) => { const j = Math.floor((t - thG[0]) / dth); if (j >= 0 && j < thG.length) E[j] += ws[i] / tot / dth; });
    return { E, ths, ws: ws.map(w => w / tot), min: Math.min(...ths) };
  }
  ST.tool({
    id: 'laminar', tag: 'F-1501', ch: 15, sec: 'Ch. 15', title: 'Laminar flow reactors: the convection model, its RTD and conversion',
    lede: 'Which model applies — pure convection, dispersion or something between — from Re, Sc and L/d. The E, E* and E** curves that different tracer injection and detection give, and their F curves. Conversion for any kinetics in Newtonian pipe flow, power-law fluids, slits and falling films, against plug flow, with the zero-, first- and second-order closed forms, and the yield of an intermediate in laminar flow.',
    keys: 'laminar flow convection model parabolic velocity profile E star E double star flux planar measurement non newtonian power law falling film zero first second order bosworth cleland wilhelm denbigh breakthrough time regime map intermediate series reaction',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'Task' },
      { k: 'what', t: 'select', label: 'Show', opts: [['regime', 'which model applies'], ['curves', 'E, E*, E** and F curves'], ['react', 'conversion in laminar flow'], ['series', 'intermediate yield, A → R → S']], v: 'react' },
      { k: 'Re', label: 'Reynolds number', v: 100, min: 1e-4, max: 2100, scale: 'log', if: v => v.what === 'regime' },
      { k: 'Sc', label: 'Schmidt number', v: 1000, min: 0.1, max: 1e7, scale: 'log', if: v => v.what === 'regime' },
      { k: 'LD', label: 'L/d_t', v: 100, min: 1e-2, max: 1e7, scale: 'log', if: v => v.what === 'regime' },
      { k: 'prof', t: 'select', label: 'Channel and fluid', opts: [['pipe', 'Newtonian, circular pipe'], ['power', 'power-law fluid in a pipe'], ['slit', 'Newtonian, between parallel plates'], ['film', 'falling film']], v: 'pipe', if: v => v.what !== 'regime' },
      { k: 'n', label: 'Flow index n (1 = Newtonian)', v: 0.5, min: 0.05, max: 3, if: v => v.what !== 'regime' && v.prof === 'power' },
      { k: 'tau', label: 'Mean residence time t̄', unit: '[t]', v: 1, min: 1e-9, max: 1e9, scale: 'log', if: v => v.what === 'react' || v.what === 'series' },
      { k: 'kAR', label: 'k₁ (A → R)', unit: '1/[t]', v: 1, min: 1e-9, max: 1e9, scale: 'log', if: v => v.what === 'series' },
      { k: 'kRS', label: 'k₂ (R → S)', unit: '1/[t]', v: 0.5, min: 1e-9, max: 1e9, scale: 'log', if: v => v.what === 'series' }
    ].concat(KIT.kin({ gas: false, d: { law: 'power', k: 1.6094, a: 1, CA0: 1 } }).map(f => Object.assign({}, f, { if: f.group ? (vv => vv.what === 'react') : f.if ? (vv => vv.what === 'react' && f.if(vv)) : (vv => vv.what === 'react') }))),
    presets: [
      { name: 'Problem 15.1 (plug would give 80%)', v: { what: 'react', prof: 'pipe', tau: 1, law: 'power', a: 1, ob: 0, k: 1.6094, CA0: 1 }, note: 'A viscous liquid, first order, in laminar flow; plug flow of the same t̄ would convert 80% (kt̄ = ln 5).' },
      { name: 'Zero order (Bosworth)', v: { what: 'react', prof: 'pipe', tau: 1, law: 'power', a: 0, ob: 0, k: 0.8, CA0: 1 }, note: 'Zero order with kt̄/C_A0 = 0.8: compare with (1 − kt̄/2C_A0)².' },
      { name: 'Second order (Denbigh)', v: { what: 'react', prof: 'pipe', tau: 1, law: 'power', a: 2, ob: 0, k: 4, CA0: 1 }, note: 'kC_A0t̄ = 4.' },
      { name: 'Shear-thinning polymer melt', v: { what: 'react', prof: 'power', n: 0.4, tau: 1, law: 'power', a: 1, ob: 0, k: 1.6094, CA0: 1 }, note: 'A flatter velocity profile narrows the RTD: closer to plug flow than a Newtonian liquid.' },
      { name: 'E, E*, E** (Fig. 15.7)', v: { what: 'curves', prof: 'pipe' } },
      { name: 'Regime map (Fig. 15.2)', v: { what: 'regime', Re: 100, Sc: 1000, LD: 100 } },
      { name: 'A → R → S in laminar flow (Fig. 15.11)', v: { what: 'series', prof: 'pipe', tau: 1, kAR: 1, kRS: 1 } }
    ],
    eqs: [
      ['E = \\frac{\\bar t^2}{2t^3},\\quad E^* = \\frac{\\bar t}{2t^2},\\quad E^{**} = \\frac{1}{2t}\\qquad (t \\ge \\bar t/2)', 'Eqs. 15.3–15.5'],
      ['F = 1-\\frac{1}{4\\theta^2},\\qquad F^* = 1-\\frac{1}{2\\theta}\\qquad (\\theta \\ge \\tfrac12)', 'Eq. 15.6'],
      ['\\frac{C_A}{C_{A0}} = \\left(1-\\frac{k\\bar t}{2C_{A0}}\\right)^2', 'Eq. 15.9', 'Zero order (Bosworth).'],
      ['\\frac{C_A}{C_{A0}} = y^2\\,\\mathrm{ei}(y) + (1-y)e^{-y},\\quad y = \\frac{k\\bar t}{2}', 'Eq. 15.10', 'First order (Cleland and Wilhelm).'],
      ['\\frac{C_A}{C_{A0}} = 1 - kC_{A0}\\bar t\\left[1-\\frac{kC_{A0}\\bar t}{2}\\ln\\left(1+\\frac{2}{kC_{A0}\\bar t}\\right)\\right]', 'Eq. 15.11', 'Second order (Denbigh).'],
      ['\\sigma \\propto L \\ \\text{(convection)}\\quad\\text{vs}\\quad \\sigma \\propto \\sqrt{L}\\ \\text{(dispersion)}', 'Eqs. 15.1–15.2']
    ],
    method: ['In pure convection each streamline is its own batch reactor (a macrofluid), so the outlet is the flux-weighted average of batch results over the velocity profile — integrated directly over the radius rather than over E, which avoids the singular tail of E.', 'Non-Newtonian and other channels use their velocity profiles: power-law u/ū = (3n+1)/(n+1)·[1 − (r/R)^{(n+1)/n}], parallel plates and falling films 1.5[1 − ξ²]. E curves are built by binning flux-weighted residence times.', 'The regime test uses the radial-diffusion time ratio 𝒟t̄/R² = 4(L/d)/(Re·Sc): above about 0.8 Taylor dispersion holds, below about 0.003 molecular diffusion is negligible and the convection model holds.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      if (v.what === 'regime') {
        const Re = +v.Re, Sc = +v.Sc, LD = +v.LD, Bo = Re * Sc, tD = 4 * LD / Bo;
        const reg = tD > 0.8 ? (Bo < 1 ? 'pure diffusion' : 'dispersion (Taylor–Aris)') : tD < 0.003 ? 'pure convection' : 'between convection and dispersion';
        kpi.push({ k: 'Bo', l: 'Bodenstein number Re·Sc', v: Bo }, { k: 'tD', l: '𝒟t̄/R² = 4(L/d)/(Re·Sc)', v: tD, main: true }, { k: 'reg', l: 'Model to use', v: reg, main: true });
        if (Re > 2100) warn.push('Re above 2100: turbulent, so use the dispersion correlations of Chapter 13.');
        if (reg === 'between convection and dispersion') info.push('In the no-man’s-land: compute with both bounding models and average, as the book suggests.');
        if (reg === 'dispersion (Taylor–Aris)') kpi.push({ k: 'dL', l: 'D/uL from Taylor–Aris', v: (1 / Bo + Bo / 192) / LD });
        /* map */
        const LDs = N.logspace(0.01, 1e7, 120);
        plots.push({ title: 'Which model applies (after Fig. 15.2)', x: { label: 'L/d_t', log: true }, y: { label: 'Re·Sc', log: true, min: 1e-2, max: 1e8 }, series: [{ name: 'pure convection below this line (𝒟t̄/R² = 0.003)', x: LDs, y: LDs.map(x => 4 * x / 0.003), color: 1 }, { name: 'dispersion above this line (𝒟t̄/R² = 0.8)', x: LDs, y: LDs.map(x => 4 * x / 0.8), color: 2 }, { name: 'Re·Sc = 1 (diffusion dominates below)', x: [0.01, 1e7], y: [1, 1], color: 8, dash: '5 4' }], marks: [{ t: 'pt', x: LD, y: Bo, color: 3, label: 'this flow' }], cap: 'Convection regime lies above-left of the first line (short tubes, high Re·Sc); dispersion regime below-right of the second.' });
        return { kpi, warn, info, plots, tables };
      }
      const prof = profileRTD(v.prof, +v.n);
      if (v.what === 'curves') {
        const th = N.linspace(0.3, 3, 300);
        const Eb = Etheta(prof, N.linspace(0, 4, 801));
        const thB = N.linspace(0, 4, 801).map(x => x + 0.0025);
        const ser = [{ name: 'E_θ (flux in, flux out): the RTD', x: thB, y: Eb.E, color: 1 }];
        if (v.prof === 'pipe') ser.push({ name: 'E*_θ (one planar end)', x: th, y: th.map(q => (q >= 0.5 ? 1 / (2 * q * q) : 0)), color: 2 }, { name: 'E**_θ (planar in and out)', x: th, y: th.map(q => (q >= 0.5 ? 1 / (2 * q) : 0)), color: 3 });
        plots.push({ title: 'Pulse responses (Fig. 15.7) — ' + prof.label, x: { label: 'θ', min: 0, max: 3 }, y: { label: 'E_θ', min: 0, max: 4.5 }, series: ser });
        if (v.prof === 'pipe') plots.push({ title: 'Step responses (Fig. 15.9)', x: { label: 'θ', min: 0, max: 3 }, y: { label: 'F', min: 0, max: 1 }, series: [{ name: 'F (flux)', x: th, y: th.map(q => (q >= 0.5 ? 1 - 1 / (4 * q * q) : 0)), color: 1 }, { name: 'F* (planar measurement)', x: th, y: th.map(q => (q >= 0.5 ? 1 - 1 / (2 * q) : 0)), color: 2 }] });
        kpi.push({ k: 'thb', l: 'Breakthrough θ (first fluid out)', v: Eb.min, main: true, s: 'velocity ratio u_max/ū = ' + N.fmt(1 / Eb.min, 4) });
        info.push('E is the RTD for reactor work. E* and E** come from planar injection or detection and have infinite variance: convert them before use (E = E*/θ = E**/θ²).');
        return { kpi, warn, info, plots, tables };
      }
      /* conversion by integrating over the cross-section */
      const tb = +v.tau;
      const flowAvg = f => { const g = x => prof.w(x) * prof.u(x) * f(tb / prof.u(x)); return N.quad(g, 0, 1 - 1e-12, { rtol: 1e-9 }) / N.quad(x => prof.w(x) * prof.u(x), 0, 1, { rtol: 1e-12 }); };
      if (v.what === 'series') {
        const k1 = +v.kAR, k2 = +v.kRS;
        const R = t => (Math.abs(k1 - k2) < 1e-12 * k1 ? k1 * t * Math.exp(-k1 * t) : k1 / (k2 - k1) * (Math.exp(-k1 * t) - Math.exp(-k2 * t)));
        const CR = flowAvg(R), CA = flowAvg(t => Math.exp(-k1 * t));
        const RP = R(tb), RM = k1 * tb / ((1 + k1 * tb) * (1 + k2 * tb));
        kpi.push({ k: 'CR', l: 'C_R/C_A0, laminar flow', v: CR, main: true, s: 'X_A = ' + N.fmt(1 - CA, 4) }, { k: 'CRp', l: 'C_R/C_A0, plug flow, same t̄', v: RP }, { k: 'CRm', l: 'C_R/C_A0, mixed flow, same t̄', v: RM });
        /* best over the mean residence time */
        const lamAt = t2 => flowAvg2(t2), bestL = N.fminScan(q => -lamAt(Math.exp(q)), Math.log(0.02 / k1), Math.log(50 / k1), 60);
        function flowAvg2(t2) { const g = x => prof.w(x) * prof.u(x) * R(t2 / prof.u(x)); return N.quad(g, 0, 1 - 1e-12, { rtol: 1e-9 }) / N.quad(x => prof.w(x) * prof.u(x), 0, 1, { rtol: 1e-12 }); }
        const RPmax = Math.abs(k1 - k2) < 1e-12 * k1 ? Math.exp(-1) : Math.pow(k1 / k2, k2 / (k2 - k1)), RMmax = 1 / Math.pow(Math.sqrt(k2 / k1) + 1, 2);
        kpi.push({ k: 'CRmaxL', l: 'Most R, laminar flow', v: -bestL.f, main: true, s: 'at t̄ = ' + N.fmt(Math.exp(bestL.x), 4) }, { k: 'CRmaxP', l: 'Most R, plug flow (Eq. 8.8)', v: RPmax }, { k: 'CRmaxM', l: 'Most R, mixed flow (Eq. 8.15)', v: RMmax });
        kpi.push({ k: 'pos', l: 'Laminar maximum, from plug (0) to mixed (1)', v: (RPmax + bestL.f) / (RPmax - RMmax), s: 'the book reads about 0.2 from Fig. 15.11' });
        const ts = N.logspace(0.01, 50, 60).map(q => q / k1);
        const pts = ts.map(t => { const tbb = t; const g = x => prof.w(x) * prof.u(x) * (Math.abs(k1 - k2) < 1e-12 * k1 ? k1 * (tbb / prof.u(x)) * Math.exp(-k1 * tbb / prof.u(x)) : k1 / (k2 - k1) * (Math.exp(-k1 * tbb / prof.u(x)) - Math.exp(-k2 * tbb / prof.u(x)))); const ga = x => prof.w(x) * prof.u(x) * Math.exp(-k1 * tbb / prof.u(x)); const norm = N.quad(x => prof.w(x) * prof.u(x), 0, 1); return [1 - N.quad(ga, 0, 1 - 1e-12) / norm, N.quad(g, 0, 1 - 1e-12) / norm]; });
        const pp = ts.map(t => [1 - Math.exp(-k1 * t), R(t)]), mm = ts.map(t => [k1 * t / (1 + k1 * t), k1 * t / ((1 + k1 * t) * (1 + k2 * t))]);
        plots.push({ title: 'Intermediate against conversion (Fig. 15.11)', x: { label: 'X_A', min: 0, max: 1 }, y: { label: 'C_R/C_A0', min: 0 }, series: [{ name: 'plug flow', x: pp.map(p => p[0]), y: pp.map(p => p[1]), color: 1 }, { name: 'laminar flow', x: pts.map(p => p[0]), y: pts.map(p => p[1]), color: 2 }, { name: 'mixed flow', x: mm.map(p => p[0]), y: mm.map(p => p[1]), color: 3 }], marks: [{ t: 'pt', x: 1 - CA, y: CR, color: 2, label: 'this t̄' }] });
        return { kpi, warn, info, plots, tables };
      }
      const m = KIT.model(v);
      const umin = prof.u(0.999999), tl = tb / Math.max(umin, 1e-9);
      const Xb = RT.batchX(m, Math.min(tl, tb * 1e4));
      const C = flowAvg(t => 1 - Xb(t));
      const Cp = 1 - K.XPFR(m, tb).X, Cm = 1 - K.XMFR(m, tb);
      kpi.push({ k: 'C', l: 'C_A/C_A0 leaving, laminar flow', v: C, main: true, s: 'X_A = ' + N.fmt(1 - C, 4) }, { k: 'Cp', l: 'Plug flow, same t̄', v: Cp }, { k: 'Cm', l: 'Mixed flow, same t̄', v: Cm });
      if (v.prof === 'pipe' && v.law === 'power') {
        const a = +v.a, k = +v.k, C0 = +v.CA0;
        let cf = NaN;
        if (a === 0) { const z = k * tb / (2 * C0); cf = z <= 1 ? (1 - z) * (1 - z) : 0; }
        else if (a === 1) { const y = k * tb / 2; cf = y * y * N.E1(y) + (1 - y) * Math.exp(-y); }
        else if (a === 2) { const z = k * C0 * tb; cf = 1 - z * (1 - z / 2 * Math.log(1 + 2 / z)); }
        if (isFinite(cf)) kpi.push({ k: 'Cform', l: 'Closed form (Eqs. 15.9–15.11)', v: cf });
      }
      /* Fig. 15.10: laminar vs plug — one batch history serves every mean residence time */
      const top = Math.min(m.Xe, m.Xmax);
      const tpLo = K.tauPFR(m, Math.min(0.05, top * 0.5)), tpHi = K.tauPFR(m, Math.min(0.995, top * 0.999));
      const tbs = N.logspace(tpLo * 0.5, tpHi * 30, 40);
      const XbAll = RT.batchX(m, tbs[tbs.length - 1] * 60);
      const norm = N.quad(x => prof.w(x) * prof.u(x), 0, 1, { rtol: 1e-12 });
      const Cof = t2 => N.quad(x => prof.w(x) * prof.u(x) * (1 - XbAll(t2 / prof.u(x))), 0, 1 - 1e-12, { rtol: 1e-8 }) / norm;
      const Cs = tbs.map(Cof);
      const fr = N.logspace(0.01, 0.95, 24);
      const pts = fr.map(cf => { if (!(1 - cf < top)) return [cf, NaN]; const tp = K.tauPFR(m, 1 - cf); const lc = Cs.map(c => Math.log(Math.max(c, 1e-300))); if (Math.log(cf) > lc[0] || Math.log(cf) < lc[lc.length - 1]) return [cf, NaN]; const t2 = Math.exp(N.interp(lc.slice().reverse(), tbs.map(Math.log).reverse(), Math.log(cf))); return [cf, t2 / tp]; });
      plots.push({ title: 'Laminar flow against plug flow (Fig. 15.10)', x: { label: 'C_A/C_A0 leaving', log: true, min: 0.01, max: 1 }, y: { label: 'V_laminar/V_plug', min: 1 }, series: [{ name: prof.label, x: pts.map(p => p[0]), y: pts.map(p => p[1]), color: 1, marker: 'o' }], cap: 'Even at high conversion convection costs little compared with dispersion or backmixing.' });
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
