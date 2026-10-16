/* Chapter 9 — multiple steady states in mixed flow; temperature and product distribution */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, TH = ST.TH, D = ST.D, XTI = ST.XTI;
  const C2K = TH.C2K, K2C = TH.K2C, R = TH.R;

  /* ===================================================================== */
  /* H-904  multiple steady states, ignition, extinction, stability          */
  /* ===================================================================== */
  const EX95 = { form: 'A=R', rev: true, phase: 'liq', CA0: 4, MR: 0, ksrc: 'arr', A1: Math.exp(17.34), E1: 48900, Ksrc: 'arr2', A2: Math.exp(42.04), E2: 124200, dH: -75300, dCp: 0, Cp: 1046, Tlo: -20, Thi: 110 };
  const inp904 = XTI.inputs(EX95).map(f => {
    if (f.k === 'Tlo') return Object.assign({}, f, { label: 'Chart from', hint: undefined });
    if (f.k === 'Thi') return Object.assign({}, f, { label: 'Chart to' });
    if (f.group === 'Energy and temperature limits') return { group: 'Energy and chart range' };
    if (f.k === 'phase') return Object.assign({}, f, { opts: [['liq', 'liquid, constant density']], v: 'liq' });
    return f;
  });

  ST.tool({
    id: 'mfrstates', tag: 'H-904', ch: 9, sec: 'Ch. 9.1', title: 'Exothermic mixed flow: multiple steady states, ignition, extinction, stability',
    lede: 'van Heerden’s construction: the S-shaped material balance against the energy-balance line, with wall cooling if you want it. Every steady state with its slope and dynamic stability, the ignition and extinction feed temperatures, hysteresis against feed temperature and space time, and start-up trajectories in the phase plane.',
    keys: 'multiple steady states van heerden ignition extinction hysteresis stability slope condition limit cycle oscillation cstr mixed flow exothermic phase plane start up heat generation removal semenov s-curve runaway',
    units: { t: 'min' }, slow: true,
    inputs: inp904.concat([
      { group: 'Mixed reactor' },
      { k: 'tau', label: 'Space time τ', unit: '[t]', v: 8.105, min: 1e-6, max: 1e7, scale: 'log' },
      { k: 'T1', label: 'Feed temperature', unit: '°C', v: 4.46, min: -150, max: 1500 },
      { k: 'Uv', label: 'Wall cooling UA/V', unit: 'W/(L·K)', v: 0, min: 0, max: 1e4, hint: 'Zero for an adiabatic reactor.' },
      { k: 'Tc', label: 'Coolant temperature', unit: '°C', v: 20, min: -150, max: 1500, if: v => +v.Uv > 0 },
      { group: 'Start-up and sweeps' },
      { k: 'X0', label: 'Tank starts at conversion', v: 0, min: 0, max: 0.999 },
      { k: 'T0', label: 'Tank starts at', unit: '°C', v: 25, min: -150, max: 1500 },
      { k: 'tEnd', label: 'Simulate for', unit: 'τ', v: 20, min: 1, max: 1000, scale: 'log' },
      { k: 'portrait', t: 'check', label: 'Phase portrait', v: true },
      { k: 'sTa', label: 'Space-time sweep from', unit: '[t]', v: 0.05, min: 1e-6, max: 1e7, scale: 'log' },
      { k: 'sTb', label: 'Space-time sweep to', unit: '[t]', v: 500, min: 1e-6, max: 1e8, scale: 'log' }
    ]),
    presets: [
      { name: 'Example 9.5 design point', v: Object.assign({}, EX95, { tau: 8.105, T1: 4.457, Uv: 0, X0: 0, T0: 25, tEnd: 20, sTa: 0.05, sTb: 500 }), book: { Xup: { v: 0.8, tol: 0.01 }, Tup: { v: 62, tol: 0.02 } },
        note: 'The 2000 L mixed reactor of Example 9.5 fed at its adiabatic inlet temperature: one steady state, the design point at X = 0.8 and about 62 °C. Its hysteresis band lies just below, for feeds between about −7 and −4.7 °C.' },
      { name: 'Same reactor inside the hysteresis band', v: Object.assign({}, EX95, { tau: 8.105, T1: -6, Uv: 0, X0: 0, T0: -6, tEnd: 40, sTa: 0.05, sTb: 500 }), note: 'Fed at −6 °C the tank has three steady states. Started cold it stays cold; set the start-up temperature to 40 °C and it lights off to the hot branch.' },
      { name: 'Irreversible, classic S-curve', v: { form: 'A=R', rev: false, oA: 1, phase: 'liq', CA0: 2, MR: 0, ksrc: 'ref', k1r: 0.002, Tk: 25, E1: 80000, dH: -100000, dCp: 0, Cp: 1000, Tlo: -40, Thi: 160, tau: 10, T1: 25, Uv: 0, X0: 0, T0: 25, tEnd: 30, sTa: 0.01, sTb: 1000 }, note: 'Adiabatic temperature rise 100 K, E = 80 kJ/mol: three steady states at a 25 °C feed. The feed must be raised past about 32 °C to ignite, and once lit the tank stays lit down to about 3 °C.' },
      { name: 'Cooled tank that oscillates', v: { form: 'A=R', rev: false, oA: 1, phase: 'liq', CA0: 2, MR: 0, ksrc: 'ref', k1r: 0.002, Tk: 25, E1: 100000, dH: -100000, dCp: 0, Cp: 1000, Tlo: -10, Thi: 200, tau: 3, T1: 0, Uv: 20, Tc: 80, X0: 0, T0: 80, tEnd: 60, sTa: 0.01, sTb: 1000 }, note: 'Strong wall cooling at 80 °C with a cold feed: the single steady state is an unstable focus, so the tank cycles between ignition and quench indefinitely.' }
    ],
    eqs: [
      ['X_A = \\frac{\\tau\\,(-r_A)}{C_{A0}}', 'Eq. 5.11', 'Material balance: the S-shaped curve of Fig. 9.13.'],
      ['X_A = \\frac{C\'_p(T-T_1) + \\dfrac{UA}{F_{A0}}(T-T_c)}{-\\Delta H_r}', 'Eq. 25', 'Energy balance with wall cooling.'],
      ['\\frac{dX_A}{dt} = -\\frac{X_A}{\\tau} + \\frac{-r_A}{C_{A0}},\\quad C\'_p\\frac{dT}{dt} = \\frac{C\'_p(T_1-T)}{\\tau} + (-\\Delta H_r)\\frac{-r_A}{C_{A0}} - \\frac{UA}{V C_{A0}}(T-T_c)', '', 'Transient balances; their Jacobian decides stability.'],
      ['\\text{saddle if } \\det J < 0 \\iff \\left.\\frac{dX}{dT}\\right|_{\\text{MB}} > \\left.\\frac{dX}{dT}\\right|_{\\text{EB}};\\quad \\text{unstable focus if } \\operatorname{tr} J > 0', '', 'Slope and dynamic conditions.']
    ],
    method: ['The material-balance curve is solved at each temperature by bracketed root-finding (unique for positive-order kinetics); steady states are sign changes of X_MB − X_EB on a fine grid, refined by Brent’s method.', 'Feed-temperature hysteresis is traced exactly as a parametric curve: at each reactor temperature the material balance fixes X and the energy balance then gives the feed temperature, so all three branches appear without root-finding. Ignition and extinction are its turning points.', 'Stability comes from the eigenvalues of the 2×2 Jacobian of the transient balances (finite differences). Start-up and the phase portrait integrate the same balances with an adaptive Runge–Kutta method.'],
    run(v) {
      const m = XTI.model(v), kpi = [], warn = [], info = [];
      if (m.gas) throw new Error('This analysis assumes constant density.');
      const tau = +v.tau, T1 = C2K(+v.T1), Uv = +v.Uv, Tc = C2K(+v.Tc), CA0 = m.CA0f(300), tsec = m.tsec;
      const uw = Uv * tsec / CA0, UAF = uw * tau;
      const top = T => (m.rev ? m.Xe(T) : m.Xmax) * (1 - 1e-12);
      const XMBt = (T, ta) => { const xe = top(T), g = X => ta * m.rate(X, T) / CA0 - X; const g0 = g(0); if (!(g0 > 0)) return 0; return N.brent(g, 0, xe, { fa: g0, fb: g(xe), tol: 1e-14 }); };
      const XMB = T => XMBt(T, tau);
      const XEB = T => (m.Cp * (T - T1) + UAF * (T - Tc)) / (-m.dHA(T));
      const f = (t, y) => { const r = m.rate(y[0], y[1]) / CA0; return [-y[0] / tau + r, ((T1 - y[1]) / tau * m.Cp + (-m.dHA(y[1])) * r - uw * (y[1] - Tc)) / m.Cp]; };
      const Ta = m.Tlo, Tb = m.Thi;
      if (!(Tb > Ta + 1)) throw new Error('Widen the chart range.');
      /* steady states */
      const Ts = N.linspace(Ta, Tb, 700), Fv = Ts.map(T => XMB(T) - XEB(T));
      const ss = [];
      for (let i = 1; i < Ts.length; i++) if (Fv[i - 1] * Fv[i] <= 0 && Fv[i - 1] !== Fv[i]) {
        const T = N.brent(T => XMB(T) - XEB(T), Ts[i - 1], Ts[i], { tol: 1e-10 });
        if (isFinite(T) && !ss.some(q => Math.abs(q.T - T) < 1e-6)) ss.push({ T, X: XMB(T) });
      }
      ss.forEach(q => {
        const hx = 1e-6, ht = 1e-4;
        const a = f(0, [q.X + hx, q.T]), b = f(0, [q.X - hx, q.T]), c = f(0, [q.X, q.T + ht]), d = f(0, [q.X, q.T - ht]);
        const J = [[(a[0] - b[0]) / (2 * hx), (c[0] - d[0]) / (2 * ht)], [(a[1] - b[1]) / (2 * hx), (c[1] - d[1]) / (2 * ht)]];
        const tr = J[0][0] + J[1][1], det = J[0][0] * J[1][1] - J[0][1] * J[1][0], disc = tr * tr - 4 * det;
        q.tr = tr; q.det = det;
        q.ev = disc >= 0 ? [(tr + Math.sqrt(disc)) / 2, (tr - Math.sqrt(disc)) / 2] : [{ re: tr / 2, im: Math.sqrt(-disc) / 2 }];
        q.type = det < 0 ? 'saddle (unstable)' : tr < 0 ? (disc >= 0 ? 'stable node' : 'stable focus') : (disc >= 0 ? 'unstable node' : 'unstable focus');
        q.stable = det > 0 && tr < 0;
        q.sMB = (XMB(q.T + 0.01) - XMB(q.T - 0.01)) / 0.02; q.sEB = (XEB(q.T + 0.01) - XEB(q.T - 0.01)) / 0.02;
      });
      if (!ss.length) warn.push('No steady state inside the chart range; widen it.');
      const st = ss.filter(q => q.stable);
      kpi.push({ k: 'nSS', l: 'Steady states', v: ss.length, main: true, s: ss.map(q => q.stable ? 'stable' : 'unstable').join(', ') || 'none in range' });
      if (ss.length) {
        const up = ss[ss.length - 1], lo = ss[0];
        kpi.push({ k: 'Xup', l: 'Upper steady state: X_A', v: up.X, main: true, s: up.type }, { k: 'Tup', l: 'Upper steady state: T', v: K2C(up.T), u: '°C' });
        if (ss.length > 1) kpi.push({ k: 'Xlo', l: 'Lower steady state: X_A', v: lo.X, s: lo.type + ', ' + N.fmt(K2C(lo.T), 4) + ' °C' });
        if (!st.length) warn.push('No steady state is stable: expect sustained oscillations (a limit cycle).');
        if (ss.length > 1) info.push('With ' + ss.length + ' steady states, where the reactor ends up depends on how it is started. The middle state is the ignition point: the start-up trajectory decides which side it falls to.');
      }
      /* start-up */
      const tEnd = +v.tEnd * tau;
      const sim = N.ode(f, 0, [+v.X0, C2K(+v.T0)], tEnd, { rtol: 1e-8, atol: 1e-10, minPts: 300 });
      const yE = sim.last();
      kpi.push({ k: 'Xdyn', l: 'After ' + N.fmt(+v.tEnd, 3) + 'τ from the start-up state', v: yE[0], s: N.fmt(K2C(yE[1]), 4) + ' °C' });
      /* feed-temperature hysteresis (parametric) */
      const Tg = N.linspace(Ta, Tb, 600);
      const curve = Tg.map(T => { const X = XMB(T); return { T, X, T1: T - ((-m.dHA(T)) * X - UAF * (T - Tc)) / m.Cp }; });
      const turns = [];
      for (let i = 1; i < curve.length - 1; i++) {
        const a = curve[i].T1 - curve[i - 1].T1, b = curve[i + 1].T1 - curve[i].T1;
        if (a * b < 0) turns.push({ T1: curve[i].T1, T: curve[i].T, X: curve[i].X, kind: a > 0 ? 'max' : 'min' });
      }
      const ign = turns.find(q => q.kind === 'max'), ext = turns.find(q => q.kind === 'min' && (!ign || q.T > ign.T));
      if (ign) kpi.push({ k: 'T1ign', l: 'Ignition feed temperature', v: K2C(ign.T1), u: '°C', s: 'lower branch ends at ' + N.fmt(K2C(ign.T), 4) + ' °C' });
      if (ext) kpi.push({ k: 'T1ext', l: 'Extinction feed temperature', v: K2C(ext.T1), u: '°C', s: 'upper branch ends at ' + N.fmt(K2C(ext.T), 4) + ' °C' });
      /* plots */
      const C = K2C;
      const tg = N.linspace(Ta, Tb, 400);
      const ser = [{ name: 'material balance', x: tg.map(C), y: tg.map(XMB), color: 1, w: 2.2 }, { name: 'energy balance', x: tg.map(C), y: tg.map(XEB), color: 6, w: 2 }];
      if (m.rev) ser.push({ name: 'equilibrium', x: tg.map(C), y: tg.map(T => m.Xe(T)), color: 2, dash: '5 4' });
      if (v.portrait) {
        const starts = [];
        N.linspace(Ta, Tb, 7).forEach(T => { starts.push([0, T]); starts.push([Math.min(0.98, top(T) * 0.98), T]); });
        starts.forEach((s0, i) => {
          const so = N.ode(f, 0, s0, tEnd, { rtol: 1e-6, atol: 1e-9, minPts: 120 });
          ser.push({ x: so.y.map(y => C(y[1])), y: so.y.map(y => y[0]), color: 'var(--muted)', w: 0.8, hover: false, legend: false, noRange: true });
        });
      }
      ser.push({ name: 'start-up path', x: sim.y.map(y => C(y[1])), y: sim.y.map(y => y[0]), color: 3, w: 2 });
      const mk = ss.map((q, i) => ({ t: 'pt', x: C(q.T), y: q.X, color: q.stable ? 4 : 2, open: !q.stable, r: 5.5, label: q.stable ? 'stable' : 'unstable' }));
      mk.push({ t: 'pt', x: +v.T0, y: +v.X0, color: 3, r: 3.5 });
      const plots = [{ title: 'Material and energy balances (Figs. 9.13–9.15)', x: { label: 'T', unit: '°C', min: C(Ta), max: C(Tb) }, y: { label: 'X_A', min: 0, max: 1 }, series: ser, marks: mk, aspect: 0.5, hmax: 560, wide: true, cap: 'Filled points are stable, open points unstable. Thin grey lines are trajectories from the edges of the chart.' }];
      const dH = T => -m.dHA(T);
      plots.push({ title: 'Heat generated and removed, per mol A fed', x: { label: 'T', unit: '°C', min: C(Ta), max: C(Tb) }, y: { label: 'heat', unit: 'kJ/mol A' }, series: [{ name: 'heat generated', x: tg.map(C), y: tg.map(T => dH(T) * XMB(T) / 1000), color: 1 }, { name: 'heat removed', x: tg.map(C), y: tg.map(T => (m.Cp * (T - T1) + UAF * (T - Tc)) / 1000), color: 6 }], cap: 'Generated: (−ΔH_r)·X from the material balance. Removed: C′_p(T − T₁) plus wall cooling. Crossings are steady states.', marks: ss.map(q => ({ t: 'pt', x: C(q.T), y: dH(q.T) * q.X / 1000, color: q.stable ? 4 : 2, open: !q.stable })) });
      /* hysteresis in feed temperature: split stable / unstable */
      const stab = curve.map((p, i) => i > 0 && i < curve.length - 1 ? (curve[i + 1].T1 - curve[i - 1].T1) > 0 : true);
      const sA = curve.map((p, i) => (stab[i] ? C(p.T) : NaN)), sB = curve.map((p, i) => (!stab[i] || (i > 0 && !stab[i - 1]) || (i < curve.length - 1 && !stab[i + 1]) ? C(p.T) : NaN));
      const x1 = curve.map(p => C(p.T1));
      plots.push({ title: 'Reactor temperature against feed temperature', x: { label: 'feed T₁', unit: '°C' }, y: { label: 'reactor T', unit: '°C' }, series: [{ name: 'stable branches', x: x1, y: sA, color: 1, w: 2 }, { name: 'middle (unstable) branch', x: x1, y: sB, color: 2, dash: '5 4' }], marks: [{ t: 'v', x: +v.T1, color: 8, label: 'this feed' }].concat(ign ? [{ t: 'pt', x: C(ign.T1), y: C(ign.T), color: 3, label: 'ignition' }] : [], ext ? [{ t: 'pt', x: C(ext.T1), y: C(ext.T), color: 5, label: 'extinction' }] : []), cap: 'With a fixed space time. On the middle branch dT₁/dT < 0: those states are saddles.' });
      /* space-time sweep */
      const tA = +v.sTa, tB = +v.sTb;
      if (tB > tA) {
        if (Uv === 0) {
          const pts = Tg.map(T => { const X = XEB(T); if (!(X > 0 && X < top(T))) return null; const r = m.rate(X, T); if (!(r > 0)) return null; return { tau: X * CA0 / r, X, T }; });
          const xs = pts.map(p => (p && p.tau >= tA * 0.5 && p.tau <= tB * 2 ? p.tau : NaN)), ys = pts.map(p => (p ? p.X : NaN));
          plots.push({ title: 'Conversion against space time (feed at ' + N.fmt(+v.T1, 4) + ' °C)', x: { label: 'τ', unit: v._tu || 'min', log: true, min: tA, max: tB }, y: { label: 'X_A', min: 0, max: 1 }, series: [{ name: 'steady states', x: xs, y: ys, color: 1, w: 2 }], marks: [{ t: 'v', x: tau, color: 8, label: 'this τ' }], cap: 'Folds in the curve are ignition and extinction space times.' });
        } else {
          const taus = N.logspace(tA, tB, 70), px = [], py = [];
          taus.forEach(ta => {
            const uaf = uw * ta;
            const g = T => XMBt(T, ta) - (m.Cp * (T - T1) + uaf * (T - Tc)) / (-m.dHA(T));
            const tt = N.linspace(Ta, Tb, 260); let g0 = g(tt[0]);
            for (let i = 1; i < tt.length; i++) { const g1 = g(tt[i]); if (g0 * g1 <= 0 && g0 !== g1) { const T = N.brent(g, tt[i - 1], tt[i], { tol: 1e-8 }); px.push(ta); py.push(XMBt(T, ta)); } g0 = g1; }
          });
          plots.push({ title: 'Steady-state conversion against space time', x: { label: 'τ', unit: v._tu || 'min', log: true, min: tA, max: tB }, y: { label: 'X_A', min: 0, max: 1 }, series: [{ name: 'steady states', x: px, y: py, type: 'pts', color: 1, r: 2.6 }], marks: [{ t: 'v', x: tau, color: 8, label: 'this τ' }] });
        }
      }
      const tt = sim.t.map(t => t / tau);
      plots.push({ title: 'Start-up transient', x: { label: 't/τ', min: 0 }, y: { label: 'T', unit: '°C' }, series: [{ name: 'T', x: tt, y: sim.y.map(y => C(y[1])), color: 6 }], marks: ss.map(q => ({ t: 'h', y: C(q.T), color: q.stable ? 4 : 2 })) });
      plots.push({ title: 'Start-up transient: conversion', x: { label: 't/τ', min: 0 }, y: { label: 'X_A', min: 0, max: 1 }, series: [{ name: 'X_A', x: tt, y: sim.y.map(y => y[0]), color: 1 }] });
      const tables = [{ title: 'Steady states', cols: ['#', 'X_A', 'T, °C', 'dX/dT, MB', 'dX/dT, EB', 'trace J', 'det J', 'eigenvalues, 1/' + (v._tu || 'min'), 'type'], rows: ss.map((q, i) => [i + 1, q.X, C(q.T), q.sMB, q.sEB, q.tr, q.det, q.ev.map(e => typeof e === 'number' ? N.fmt(e, 3) : N.fmt(e.re, 3) + ' ± ' + N.fmt(e.im, 3) + 'i').join(', '), q.type]) }];
      const diag = (() => {
        let s = D.pipe([[4, 58], [70, 58]], { label: 'feed ' + N.fmt(+v.T1, 3) + ' °C', dy: -8 });
        s += D.tank(74, 26, 70, 72, { fill: yE[0], tag: 'R-906', label: N.fmt(K2C(yE[1]), 3) + ' °C, X ' + N.fmt(yE[0], 3), jacket: Uv > 0 });
        s += D.pipe([[146, 58], [230, 58]], { label: 'product', dy: -8, at: 0, dx: 14 });
        if (Uv > 0) s += D.text(109, 120, 'coolant ' + N.fmt(+v.Tc, 3) + ' °C', { fs: 10 }).replace('y="120"', 'y="10"');
        return D.svg(240, 150, s, 'mixed reactor with feed and product');
      })();
      return { kpi, warn, info, plots, tables, diagram: diag };
    }
  });

  /* ===================================================================== */
  /* H-905  temperature and product distribution                            */
  /* ===================================================================== */
  const P911 = 'A -> R ; k0 = 1e9, ER = 6000\nA -> T ; k0 = 1e7, ER = 4000\nR -> S ; k0 = 1e8, ER = 9000\nR -> U ; k0 = 1e12, ER = 12000';
  const SER = 'A -> R ; k0 = 3e6, ER = 6000\nR -> S ; k0 = 3e13, ER = 12000';

  function kf(net, r, T) { return r.k0 !== undefined ? r.k0 * Math.exp(-r.ER / T) : r.k * (r.E ? Math.exp(-r.E / R * (1 / T - 1 / net.Tref)) : 1); }
  function kb(net, r, T) { return r.kr0 !== undefined ? r.kr0 * Math.exp(-r.ERr / T) : r.kr * (r.Er ? Math.exp(-r.Er / R * (1 / T - 1 / net.Tref)) : 1); }
  const isLinear = net => net.rx.every(r => !r.fn && r.reac.length === 1 && r.reac[0][1] === 1 && (!r.rev || (r.prod.length === 1 && r.prod[0][1] === 1)));
  function Kmat(net, T) {
    const n = net.ns, K = N.range(n).map(() => new Array(n).fill(0));
    net.rx.forEach(r => {
      const i = r.reac[0][0], k = kf(net, r, T);
      K[i][i] -= k; r.prod.forEach(([j, nu]) => { K[j][i] += nu * k; });
      if (r.rev) { const j = r.prod[0][0], b = kb(net, r, T); K[j][j] -= b; K[i][j] += b; }
    });
    return K;
  }
  /* exp(A)·x by scaling and squaring with a Taylor core; plain loops (small dense matrices) */
  function expmv(A, x) {
    const n = A.length;
    let nrm = 0;
    for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += Math.abs(A[i][j]); if (s > nrm) nrm = s; }
    const sq = Math.max(0, Math.ceil(Math.log2(Math.max(nrm, 1e-300) / 0.25)));
    const sc = Math.pow(2, -sq);
    const As = new Float64Array(n * n), E = new Float64Array(n * n), T = new Float64Array(n * n), W = new Float64Array(n * n);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { As[i * n + j] = A[i][j] * sc; E[i * n + j] = i === j ? 1 : 0; T[i * n + j] = E[i * n + j]; }
    const mul = (P, Q, out) => { for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { let s = 0; for (let k = 0; k < n; k++) s += P[i * n + k] * Q[k * n + j]; out[i * n + j] = s; } };
    for (let k = 1; k <= 12; k++) { mul(T, As, W); for (let q = 0; q < n * n; q++) { T[q] = W[q] / k; E[q] += T[q]; } }
    for (let i = 0; i < sq; i++) { mul(E, E, W); E.set(W); }
    const y = new Array(n).fill(0);
    for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += E[i * n + j] * x[j]; y[i] = s; }
    return y;
  }
  /* one stage: 'M' mixed or 'P' plug at temperature T, space time tau */
  function stage(net, lin, type, C, T, tau) {
    if (tau <= 0) return C.slice();
    if (lin) {
      const K = Kmat(net, T);
      if (type === 'M') return N.solve(K.map((r, i) => r.map((q, j) => (i === j ? 1 : 0) - tau * q)), C).map(x => Math.max(x, 0));
      return expmv(K.map(r => r.map(q => q * tau)), C).map(x => Math.max(x, 0));
    }
    if (type === 'M') return ST.NR.mfr(net, C, tau, { T, guess: C }).C;
    return ST.NR.batch(net, C, tau, { T, n: 4 }).end;
  }

  ST.NETT = { stage, isLinear, kf, kb, Kmat, expmv };

  ST.tool({
    id: 'tempdist', tag: 'H-905', ch: 9, sec: 'Ch. 9.2', title: 'Temperature and product distribution: best temperatures for multiple reactions',
    lede: 'Rule 6 made quantitative. For parallel paths: the optimum temperature (Eq. 32) and the fractional-yield map over temperature and concentration. For any network with Arrhenius constants: the best temperature for every stage of a reactor train, isothermal or stagewise, rising or falling, with space times fixed or free.',
    keys: 'temperature product distribution activation energy optimum temperature parallel series falling rising temperature progression problem 9.11 fractional yield staged reactors optimize temperature profile multiple reactions arrhenius',
    units: { t: 'min' }, slow: true,
    inputs: [
      { group: 'Problem' },
      { k: 'mode', t: 'seg', label: 'Analyse', opts: [['par', 'parallel paths from A'], ['net', 'any network, staged reactors']], v: 'net' },
      { group: 'Parallel paths: rate of A → product i is k_{0i} e^{−E_i/RT} C_A^{n_i}', if: v => v.mode === 'par' },
      { k: 'use3', t: 'check', label: 'Third path (→ T)', v: true, if: v => v.mode === 'par' },
      { k: 'pk1', label: 'k₀ for → R (wanted)', v: 1e8, min: 1e-10, max: 1e30, scale: 'log', if: v => v.mode === 'par' },
      { k: 'pE1', label: 'E for → R', unit: 'kJ/mol', v: 60, min: 0, max: 400, if: v => v.mode === 'par' },
      { k: 'pn1', label: 'order for → R', v: 1, min: 0, max: 4, step: 0.05, if: v => v.mode === 'par' },
      { k: 'pk2', label: 'k₀ for → S', v: 1e5, min: 1e-10, max: 1e30, scale: 'log', if: v => v.mode === 'par' },
      { k: 'pE2', label: 'E for → S', unit: 'kJ/mol', v: 40, min: 0, max: 400, if: v => v.mode === 'par' },
      { k: 'pn2', label: 'order for → S', v: 1, min: 0, max: 4, step: 0.05, if: v => v.mode === 'par' },
      { k: 'pk3', label: 'k₀ for → T', v: 1e13, min: 1e-10, max: 1e30, scale: 'log', if: v => v.mode === 'par' && v.use3 },
      { k: 'pE3', label: 'E for → T', unit: 'kJ/mol', v: 100, min: 0, max: 400, if: v => v.mode === 'par' && v.use3 },
      { k: 'pn3', label: 'order for → T', v: 1, min: 0, max: 4, step: 0.05, if: v => v.mode === 'par' && v.use3 },
      { k: 'pC', label: 'C_A at which to compare', unit: 'mol/L', v: 1, min: 1e-8, max: 1e4, scale: 'log', if: v => v.mode === 'par' },
      { k: 'pCa', label: 'Map: lowest C_A', unit: 'mol/L', v: 0.01, min: 1e-8, max: 1e4, scale: 'log', if: v => v.mode === 'par' },
      { k: 'pCb', label: 'Map: highest C_A', unit: 'mol/L', v: 10, min: 1e-8, max: 1e5, scale: 'log', if: v => v.mode === 'par' },
      { group: 'Network and reactor train', if: v => v.mode === 'net' },
      { k: 'net', t: 'area', label: 'Reactions (k0 and ER = E/R in K, or k with E and Tref)', v: P911, rows: 6, if: v => v.mode === 'net' },
      { k: 'feed', t: 'text', label: 'Feed concentrations', v: 'A = 1', mono: true, if: v => v.mode === 'net' },
      { k: 'key', t: 'text', label: 'Key reactant', v: 'A', mono: true, if: v => v.mode === 'net' },
      { k: 'want', t: 'text', label: 'Wanted product', v: 'S', mono: true, if: v => v.mode === 'net' },
      { k: 'obj', t: 'seg', label: 'Maximize', opts: [['phi', 'fractional yield Φ'], ['prod', 'production C/C_key,0']], v: 'phi', if: v => v.mode === 'net' },
      { k: 'train', t: 'text', label: 'Reactor train (M = mixed, P = plug)', v: 'M M', mono: true, if: v => v.mode === 'net' },
      { k: 'tauMode', t: 'seg', label: 'Space times', opts: [['free', 'free (each up to a limit)'], ['total', 'fixed total, best split'], ['equal', 'fixed total, equal split']], v: 'free', if: v => v.mode === 'net' },
      { k: 'tauTot', label: 'Total space time', unit: '[t]', v: 10, min: 1e-8, max: 1e10, scale: 'log', if: v => v.mode === 'net' && v.tauMode !== 'free' },
      { k: 'tauMax', label: 'Largest space time per stage', unit: '[t]', v: 1e12, min: 1e-6, max: 1e15, scale: 'log', if: v => v.mode === 'net' && v.tauMode === 'free', hint: 'Space times range over 18 decades below this.' },
      { group: 'Temperature range' },
      { k: 'Tmin', label: 'Lowest allowed', unit: '°C', v: 10, min: -200, max: 2000 },
      { k: 'Tmax', label: 'Highest allowed', unit: '°C', v: 90, min: -200, max: 2500 }
    ],
    presets: [
      { name: 'Problem 9.11', v: { mode: 'net', net: P911, feed: 'A = 1', key: 'A', want: 'S', obj: 'phi', train: 'M M', tauMode: 'free', tauMax: 1e12, Tmin: 10, Tmax: 90 }, note: 'A → R → S (wanted) with side reactions A → T and R → U, first order, k = k₀e^{−E/RT} as given. Two mixed reactors anywhere between 10 and 90 °C: what temperatures, and what fractional yield of S?' },
      { name: 'Series A → R → S, falling progression', v: { mode: 'net', net: SER, feed: 'A = 1', key: 'A', want: 'R', obj: 'prod', train: 'P P P P P P', tauMode: 'equal', tauTot: 10, Tmin: 10, Tmax: 120 }, note: 'E₁ < E₂ for R’s formation and destruction, but at a fixed total space time a hot start pays: six plug-flow segments find a falling temperature profile.' },
      { name: 'Eq. 32: wanted path between two others', v: { mode: 'par', use3: true, pk1: 1e8, pE1: 60, pn1: 1, pk2: 1e5, pE2: 40, pn2: 1, pk3: 1e13, pE3: 100, pn3: 1, pC: 1, pCa: 0.01, pCb: 10, Tmin: 0, Tmax: 200 }, note: 'E₂ < E₁ < E₃: a low temperature favours R over T, a high one favours R over S, so an intermediate temperature is best.' },
      { name: 'Different orders and energies', v: { mode: 'par', use3: false, pk1: 1e9, pE1: 70, pn1: 2, pk2: 1e6, pE2: 45, pn2: 1, pC: 1, pCa: 0.001, pCb: 10, Tmin: 0, Tmax: 200 }, note: 'Wanted path second order with the larger E: high T and high C_A both help, so both sit at their limits.' }
    ],
    eqs: [
      ['\\frac{k_1}{k_2} = \\frac{k_{10}}{k_{20}}\\,e^{-(E_1-E_2)/RT}', 'Eq. 30', 'High T favours the step with the larger E.'],
      ['\\frac{1}{T_{opt}} = \\frac{R}{E_3-E_2}\\ln\\left[\\frac{E_3-E_1}{E_1-E_2}\\,\\frac{k_{30}}{k_{20}}\\right]', 'Eq. 32', 'A → R wanted, E₂ < E₁ < E₃, same orders.'],
      ['\\varphi\\left(\\frac{R}{A}\\right) = \\frac{k_1 C_A^{n_1}}{\\sum_i k_i C_A^{n_i}}', 'Eq. 7.7']
    ],
    method: ['Parallel paths: the fractional yield is mapped over T and C_A and maximized by bounded one-dimensional searches. Eq. 32 is evaluated when it applies, and compared with the numerical optimum.', 'Networks: every reactor stage has its own temperature (and space time unless fixed). For first-order networks the stages are solved exactly — a linear solve for mixed flow and a matrix exponential (scaling and squaring) for plug flow — otherwise by Newton continuation and adaptive Runge–Kutta integration. Temperatures are mapped through a logistic transform onto the allowed range, space times through logarithms, and the objective is maximized by Nelder–Mead from several starts, then compared with the best single temperature for the whole train.'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const Tmin = C2K(+v.Tmin), Tmax = C2K(+v.Tmax);
      if (!(Tmax > Tmin)) throw new Error('The highest temperature must be above the lowest.');
      if (v.mode === 'par') {
        const P = [[+v.pk1, +v.pE1 * 1000, +v.pn1], [+v.pk2, +v.pE2 * 1000, +v.pn2]].concat(v.use3 ? [[+v.pk3, +v.pE3 * 1000, +v.pn3]] : []);
        const nm = ['R', 'S', 'T'];
        const k = (i, T) => P[i][0] * Math.exp(-P[i][1] / (R * T));
        const phi = (T, C) => { const r = P.map((p, i) => k(i, T) * Math.pow(C, p[2])); return r[0] / N.sum(r); };
        const C0 = +v.pC;
        const best = N.fminScan(T => -phi(T, C0), Tmin, Tmax, 200);
        kpi.push({ k: 'phimax', l: 'Best φ(R/A) at C_A = ' + N.fmt(C0, 3), v: -best.f, main: true }, { k: 'Topt', l: 'at temperature', v: K2C(best.x), u: '°C', main: true, s: Math.abs(best.x - Tmax) < 0.05 ? 'upper limit' : Math.abs(best.x - Tmin) < 0.05 ? 'lower limit' : 'interior optimum' });
        const sameOrder = P.every(p => Math.abs(p[2] - P[0][2]) < 1e-12);
        if (P.length === 3 && sameOrder) {
          const [E1, E2, E3] = [P[0][1], P[1][1], P[2][1]];
          const lo = Math.min(E2, E3), hi = Math.max(E2, E3), iLo = E2 < E3 ? 1 : 2, iHi = E2 < E3 ? 2 : 1;
          if (E1 > lo && E1 < hi) {
            const Tq = (hi - lo) / (R * Math.log((hi - E1) / (E1 - lo) * P[iHi][0] / P[iLo][0]));
            kpi.push({ k: 'Teq32', l: 'Eq. 32 gives', v: K2C(Tq), u: '°C', s: Tq < Tmin || Tq > Tmax ? 'outside the allowed range' : 'inside the allowed range' });
          }
        }
        /* rules */
        P.slice(1).forEach((p, j) => info.push('R against ' + nm[j + 1] + ': E_R ' + (P[0][1] > p[1] ? '>' : P[0][1] < p[1] ? '<' : '=') + ' E_' + nm[j + 1] + ' → ' + (P[0][1] > p[1] ? 'high T helps.' : P[0][1] < p[1] ? 'low T helps.' : 'temperature does not matter.') + (Math.abs(P[0][2] - p[2]) > 1e-12 ? ' Orders ' + N.fmt(P[0][2], 3) + ' vs ' + N.fmt(p[2], 3) + ' → ' + (P[0][2] > p[2] ? 'high C_A helps.' : 'low C_A helps.') : '')));
        /* map */
        const Ts = N.linspace(Tmin, Tmax, 90), Cs = N.logspace(+v.pCa, +v.pCb, 80);
        const Z = Cs.map(c => Ts.map(T => phi(T, c)));
        let bm = { f: -1 };
        Cs.forEach((c, j) => Ts.forEach((T, i) => { if (Z[j][i] > bm.f) bm = { f: Z[j][i], T, c }; }));
        const lv = N.linspace(0.1, 0.9, 9);
        const cont = ST.Plot.contour(Ts.map(K2C), Cs.map(Math.log10), Z, lv);
        const ser = [], mks = [];
        cont.forEach(c => {
          if (!c.x.length) return;
          ser.push({ x: c.x, y: c.y.map(q => (isFinite(q) ? Math.pow(10, q) : NaN)), color: 'var(--muted)', w: 0.9, hover: false, legend: false, noRange: true });
          let bi = -1, bx = -Infinity; for (let i = 0; i < c.x.length; i++) if (isFinite(c.x[i]) && c.x[i] > bx) { bx = c.x[i]; bi = i; }
          if (bi >= 0) mks.push({ t: 'text', x: c.x[bi], y: Math.pow(10, c.y[bi]), s: N.fmt(c.level, 2), anchor: 'end', color: 'var(--muted)', noRange: true });
        });
        mks.push({ t: 'pt', x: K2C(bm.T), y: bm.c, color: 2, label: 'best φ = ' + N.fmt(bm.f, 3) });
        kpi.push({ k: 'phiMap', l: 'Best φ over the whole map', v: bm.f, s: N.fmt(K2C(bm.T), 4) + ' °C, C_A = ' + N.fmt(bm.c, 3) });
        const plots = [{ title: 'Fractional yield of R over temperature and concentration', x: { label: 'T', unit: '°C', min: K2C(Tmin), max: K2C(Tmax) }, y: { label: 'C_A', unit: 'mol/L', log: true, min: +v.pCa, max: +v.pCb }, series: ser, marks: mks, aspect: 0.5, wide: true, cap: 'Contours of φ(R/A). With equal orders they run vertically: concentration then does not matter.' }];
        const tg = N.linspace(Tmin, Tmax, 260);
        plots.push({ title: 'Instantaneous fractional yields against temperature', x: { label: 'T', unit: '°C' }, y: { label: 'φ', min: 0, max: 1 }, series: P.map((p, i) => ({ name: 'φ(' + nm[i] + '/A)', x: tg.map(K2C), y: tg.map(T => { const r = P.map((q, j) => k(j, T) * Math.pow(C0, q[2])); return r[i] / N.sum(r); }), color: i + 1 })), marks: [{ t: 'v', x: K2C(best.x), color: 8, label: 'best for R' }] });
        plots.push({ title: 'Arrhenius plot of the paths', x: { label: '1000/T', unit: '1/K' }, y: { label: 'ln k' }, series: P.map((p, i) => ({ name: 'k(→' + nm[i] + ')', x: tg.map(T => 1000 / T), y: tg.map(T => Math.log(k(i, T))), color: i + 1 })), cap: 'Steeper line = larger activation energy = gains more from heating.' });
        return { kpi, warn, info, plots, tables: [{ title: 'Yields at the compared concentration', cols: ['T, °C'].concat(P.map((p, i) => 'φ(' + nm[i] + '/A)')), rows: N.linspace(Tmin, Tmax, 9).map(T => [K2C(T)].concat(P.map((p, i) => { const r = P.map((q, j) => k(j, T) * Math.pow(C0, q[2])); return r[i] / N.sum(r); }))) }] };
      }
      /* ---------- network mode ---------- */
      const net = ST.parseNet(v.net);
      const lin = isLinear(net);
      const C0 = new Array(net.ns).fill(0);
      String(v.feed).split(/[,;\n]+/).forEach(t => {
        const mm = /^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*([-+0-9.eE]+)\s*$/.exec(t);
        if (!t.trim()) return;
        if (!mm) throw new Error('Feed: write it as “A = 1, B = 0.5”.');
        const i = net.idx(mm[1]);
        if (i < 0) throw new Error('Feed: species ' + mm[1] + ' is not in the network.');
        C0[i] = parseFloat(mm[2]);
      });
      const iK = net.idx(String(v.key).trim()), iW = net.idx(String(v.want).trim());
      if (iK < 0) throw new Error('Key reactant not found in the network.');
      if (iW < 0) throw new Error('Wanted product not found in the network.');
      if (!(C0[iK] > 0)) throw new Error('Feed some of the key reactant.');
      const types = String(v.train).toUpperCase().replace(/[^MP]/g, '').split('');
      if (!types.length) throw new Error('Describe the train with letters M and P, e.g. “M M” or “P P P”.');
      if (types.length > 8) throw new Error('Up to eight stages.');
      const ns = types.length, mode = v.tauMode, tauTot = +v.tauTot, tauMax = +v.tauMax;
      const lg = x => 1 / (1 + Math.exp(-x)), Tof = z => Tmin + (Tmax - Tmin) * lg(z);
      const unpack = (x, shared) => {
        const Ts = shared ? new Array(ns).fill(Tof(x[0])) : x.slice(0, ns).map(Tof);
        const rest = x.slice(shared ? 1 : ns);
        let taus;
        if (mode === 'free') taus = rest.map(q => tauMax * Math.exp(-Math.log(1e18) * lg(-q)));
        else if (mode === 'equal') taus = new Array(ns).fill(tauTot / ns);
        else { const e = rest.concat([0]).map(Math.exp), s = N.sum(e); taus = e.map(q => tauTot * q / s); }
        return { Ts, taus };
      };
      const run = (Ts, taus, keep) => {
        let C = C0.slice();
        const outs = [C.slice()];
        for (let i = 0; i < ns; i++) { C = stage(net, lin, types[i], C, Ts[i], taus[i]); if (keep) outs.push(C.slice()); }
        return { C, outs };
      };
      const score = C => { const made = C[iW] - C0[iW], used = C0[iK] - C[iK]; return v.obj === 'phi' ? (used > 1e-12 * C0[iK] ? made / used : 0) : made / C0[iK]; };
      const nTau = mode === 'free' ? ns : mode === 'total' ? ns - 1 : 0;
      const solve = shared => {
        const nT = shared ? 1 : ns;
        const f = x => { const u = unpack(x, shared); const s = score(run(u.Ts, u.taus).C); return isFinite(s) ? -s : 1e9; };
        const starts = [];
        [-2.5, 0, 2.5].forEach(zt => (nTau ? [0, mode === 'free' ? 3 : -1, mode === 'free' ? -3 : 1] : [0]).forEach(zq => starts.push(new Array(nT).fill(zt).concat(new Array(nTau).fill(zq)))));
        if (!shared && ns > 1) { starts.push(N.linspace(2.5, -2.5, nT).concat(new Array(nTau).fill(0))); starts.push(N.linspace(-2.5, 2.5, nT).concat(new Array(nTau).fill(0))); }
        let best = null;
        starts.forEach(x0 => {
          if (x0.length === 0) return;
          const r = N.nelderMead(f, x0, { tol: 1e-10, maxIt: 260 * x0.length, step: x0.map(() => 1.2) });
          if (!best || r.f < best.f) best = r;
        });
        if (!best) best = { x: [], f: f([]) };
        best = N.nelderMead(f, best.x, { tol: 1e-13, maxIt: 400 * Math.max(1, best.x.length), step: best.x.map(() => 0.3) });
        return Object.assign(unpack(best.x, shared), { val: -best.f });
      };
      const iso = solve(true), opt = ns > 1 ? solve(false) : iso;
      const better = opt.val >= iso.val ? opt : iso;
      const lab = v.obj === 'phi' ? 'Φ(' + v.want + '/' + v.key + ')' : 'C_' + v.want + '/C_' + v.key + '0';
      kpi.push({ k: 'best', l: 'Best ' + lab + ', stagewise temperatures', v: better.val, main: true, s: better.Ts.map(T => N.fmt(K2C(T), 3)).join(' → ') + ' °C' });
      kpi.push({ k: 'iso', l: 'Best ' + lab + ', one temperature throughout', v: iso.val, s: 'at ' + N.fmt(K2C(iso.Ts[0]), 4) + ' °C' });
      better.Ts.forEach((T, i) => kpi.push({ k: 'T' + (i + 1), l: 'Stage ' + (i + 1) + ' (' + (types[i] === 'M' ? 'mixed' : 'plug') + ') temperature', v: K2C(T), u: '°C', s: 'τ = ' + N.fmt(better.taus[i], 4) + ' ' + (v._tu || 'min') }));
      const edge = T => Math.abs(T - Tmax) < 0.02 * (Tmax - Tmin) ? 'high' : Math.abs(T - Tmin) < 0.02 * (Tmax - Tmin) ? 'low' : 'mid';
      const shape = better.Ts.map(edge);
      if (ns > 1) {
        const d = better.Ts[ns - 1] - better.Ts[0];
        info.push('Optimal progression: ' + (Math.abs(d) < 0.02 * (Tmax - Tmin) ? 'flat' : d < 0 ? 'falling' : 'rising') + ' (' + shape.join(' → ') + '). Gain over the best single temperature: ' + N.fmt(100 * (better.val - iso.val) / Math.max(Math.abs(iso.val), 1e-12), 3) + '%.');
      }
      if (mode === 'free' && better.taus.some(t => t > 0.5 * tauMax)) info.push('A stage wants as long a space time as allowed: the yield keeps improving toward complete reaction there.');
      if (ST.ADV && ST.ADV.analyze) { try { ST.ADV.analyze(net, String(v.want).trim(), String(v.key).trim()).lines.forEach(s => info.push(s)); } catch (e) { /* advisory only */ } }
      /* profiles through the train */
      const res = run(better.Ts, better.taus, true);
      const tables = [{ title: 'Stages', cols: ['stage', 'type', 'T, °C', 'τ'].concat(net.species.map(s => 'C_' + s + ' out')), rows: better.Ts.map((T, i) => [i + 1, types[i] === 'M' ? 'mixed' : 'plug', K2C(T), better.taus[i]].concat(res.outs[i + 1])) }];
      const cum = [0]; better.taus.forEach(t => cum.push(cum[cum.length - 1] + t));
      const plots = [];
      plots.push({ title: 'Stage temperatures', x: { label: mode === 'free' ? 'stage' : 'cumulative τ', unit: mode === 'free' ? '' : (v._tu || 'min') }, y: { label: 'T', unit: '°C', min: K2C(Tmin), max: K2C(Tmax) }, series: [{ name: 'stagewise optimum', x: mode === 'free' ? N.range(ns + 1) : cum, y: better.Ts.concat([better.Ts[ns - 1]]).map(K2C), type: 'step', color: 1, w: 2 }, { name: 'best single temperature', x: mode === 'free' ? [0, ns] : [0, cum[ns]], y: [K2C(iso.Ts[0]), K2C(iso.Ts[0])], color: 8, dash: '5 4' }] });
      plots.push({ title: 'Concentrations leaving each stage', x: { label: 'stage' }, y: { label: 'C', unit: 'mol/L', min: 0 }, series: net.species.map((s, j) => ({ name: 'C_' + s, x: N.range(ns + 1), y: res.outs.map(c => c[j]), color: j + 1, marker: 'o' })) });
      /* isothermal sweep */
      const tg = N.linspace(Tmin, Tmax, 41);
      const sweep = tg.map(T => { const f = x => { const u = unpack([0].concat(x), true); u.Ts = new Array(ns).fill(T); const s = score(run(u.Ts, u.taus).C); return isFinite(s) ? -s : 1e9; }; if (!nTau) return -f([]); let b = null; [-2, 0, 2].forEach(z => { const r = N.nelderMead(f, new Array(nTau).fill(z), { tol: 1e-8, maxIt: 120 * nTau, step: new Array(nTau).fill(1.5) }); if (!b || r.f < b.f) b = r; }); return -b.f; });
      plots.push({ title: 'Best ' + lab + ' with the whole train at one temperature', x: { label: 'T', unit: '°C' }, y: { label: lab, min: 0 }, series: [{ name: lab, x: tg.map(K2C), y: sweep, color: 1 }], marks: [{ t: 'h', y: better.val, color: 2, label: 'stagewise optimum' }] });
      if (!lin) info.push('Nonlinear kinetics: stages are integrated numerically, so the search is slower and may settle on a local optimum.');
      return { kpi, warn, info, plots, tables, diagram: (() => { let s = D.pipe([[4, 58], [30, 58]], { label: 'feed', dy: -8 }); let x = 34; types.forEach((t, i) => { if (t === 'M') { s += D.tank(x, 28, 54, 62, { fill: Math.min(1, (C0[iK] - res.outs[i + 1][iK]) / C0[iK]), tag: 'R-' + (i + 1), label: N.fmt(K2C(better.Ts[i]), 3) + ' °C', jacket: true }); x += 54; } else { s += D.tube(x, 47, 90, 22, { prof: [(C0[iK] - res.outs[i][iK]) / C0[iK], (C0[iK] - res.outs[i + 1][iK]) / C0[iK]], tag: 'R-' + (i + 1), label: N.fmt(K2C(better.Ts[i]), 3) + ' °C' }); x += 90; } s += D.pipe([[x + 2, 58], [x + 34, 58]], {}); x += 38; }); return D.svg(x + 10, 130, s, 'reactor train with stage temperatures'); })() };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
