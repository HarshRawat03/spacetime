/* Chapter 5 — Ideal reactors for a single reaction */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, K = ST.K, KIT = ST.KIT, D = ST.D;

  const capX = m => Math.min(m.Xe, m.Xmax) * (1 - 1e-9);
  const eqWarn = (m, X, out) => {
    if (m.Xe < m.Xmax && X >= m.Xe * 0.999) out.push('The rate falls to zero at $X_A$ = ' + N.fmt(m.Xe, 4) + ' (equilibrium). Pick a conversion below it.');
    else if (X >= m.Xmax * 0.999 && m.Xmax < 1) out.push('B runs out at $X_A$ = ' + N.fmt(m.Xmax, 4) + ' (limiting reactant).');
  };
  const sizeKPIs = (v, m, tau, kpi) => {
    const v0 = KIT.v0of(v, m.CA0);
    if (isFinite(v0) && isFinite(tau)) {
      kpi.push({ k: 'V', l: 'Reactor volume V', v: tau * v0, u: KIT.volUnit(v._cu) + ' (with v₀ in ' + KIT.volUnit(v._cu) + '/' + v._tu + ')', main: true });
      kpi.push({ k: 'v0', l: 'Feed rate v₀', v: v0, u: KIT.volUnit(v._cu) + '/' + v._tu });
    }
  };

  /* ---------------- batch ---------------- */
  ST.tool({
    id: 'batch', tag: 'H-501', ch: 5, sec: '5.1', title: 'Ideal batch reactor',
    lede: 'Reaction time to reach a conversion, or conversion after a given time, for any single rate law at constant or changing volume.',
    keys: 'batch reactor time constant volume constant pressure variable volume productivity half-life',
    units: { C: 'mol/L', t: 'min' },
    inputs: [
      { group: 'Question' },
      { k: 'mode', t: 'seg', label: 'Solve for', opts: [['t', 'time for $X_A$'], ['X', '$X_A$ after time $t$']], v: 't' },
      { k: 'Xf', label: 'Target conversion $X_A$', v: 0.9, min: 0.001, max: 0.999, step: 0.001, lo: 0, hi: 1, if: v => v.mode === 't' },
      { k: 'tEnd', label: 'Reaction time $t$', unit: '[t]', v: 10, min: 1e-3, max: 1e5, scale: 'log', lo: 0, if: v => v.mode === 'X' },
      { k: 'td', label: 'Turnaround time per batch', unit: '[t]', v: 0, min: 0, max: 1e3, lo: 0, hint: 'Filling, emptying and cleaning. Used for the average production rate.' }
    ].concat(KIT.kin({ d: { k: 0.2 } })),
    presets: [
      { name: 'Problem 3.18', v: { mode: 't', Xf: 0.9975, law: 'sat', sk1: 0.1, sk2: 0.5, CA0: 10, _cu: 'mol/L', _tu: 'min' }, note: 'Enzyme kinetics $-r_A$ = 200 $C_AC_{E0}$/(2 + $C_A$) with $C_{E0}$ = 0.001, so $k_1$ = 0.1 min⁻¹ and $k_2$ = 0.5 L/mol; from 10 down to 0.025 mol/L.' },
      { name: 'Problem 3.19', v: { mode: 'X', tEnd: 1, law: 'power', k: 3, a: 0.5, ob: 0, CA0: 1, _cu: 'mol/L', _tu: 'h' }, note: 'Half-order reaction runs to completion in finite time.' },
      { name: 'Problem 3.2', v: { mode: 't', Xf: 0.75, law: 'power', k: Math.log(2) / 5, a: 1, ob: 0, CA0: 1, _tu: 'min' }, note: 'First order with 50% in 5 min: 75% needs 10 min, five more minutes.' },
      { name: 'Problem 3.3', v: { mode: 't', Xf: 0.75, law: 'power', k: 0.2, a: 2, ob: 0, CA0: 1, _tu: 'min' }, note: 'Second order with 50% in 5 min ($kC_{A0}$ = 0.2 min⁻¹): 75% needs 15 min.' }
    ],
    eqs: [
      ['t = N_{A0}\\int_0^{X_A}\\frac{dX_A}{(-r_A)V}', 'Eq. 5.3', 'General: any temperature history, any volume change.'],
      ['t = C_{A0}\\int_0^{X_A}\\frac{dX_A}{-r_A} = -\\int_{C_{A0}}^{C_A}\\frac{dC_A}{-r_A}\\qquad(\\varepsilon_A=0)', 'Eq. 5.4'],
      ['t = C_{A0}\\int_0^{X_A}\\frac{dX_A}{(-r_A)(1+\\varepsilon_AX_A)}', 'Eq. 5.5', 'Volume changing linearly with conversion, $V = V_0(1+\\varepsilon_AX_A)$.']
    ],
    method: ['The integral is evaluated by adaptive Gauss–Kronrod quadrature; the time course comes from an adaptive Runge–Kutta (Dormand–Prince) integration of $dX_A/dt = (1+\\varepsilon_AX_A)(-r_A)/C_{A0}$.', 'Concentrations follow $C_A = C_{A0}(1-X_A)/(1+\\varepsilon_AX_A)$ and the matching stoichiometric relations for B, R and S.'],
    run(v) {
      const m = KIT.model(v), warn = [], info = [];
      const cap = capX(m);
      let Xf, t;
      if (v.mode === 't') {
        Xf = Math.min(+v.Xf, cap);
        eqWarn(m, +v.Xf, warn);
        t = K.tBatch(m, Xf);
        if (!isFinite(t)) throw new Error('This conversion is never reached with these kinetics.');
      } else {
        t = +v.tEnd;
        Xf = K.XBatch(m, t).X;
      }
      if (m.rate(0) <= 0) warn.push('The rate is zero at the start (for an autocatalytic reaction, seed some R).');
      const prof = K.XBatch(m, t, { n: 140 });
      const ts = prof.prof.t, xs = prof.prof.x;
      const c = xs.map(m.conc);
      const series = [{ name: 'C_A', x: ts, y: c.map(q => q.A), color: 1 }];
      if (m.b > 0) series.push({ name: 'C_B', x: ts, y: c.map(q => q.B), color: 2 });
      if (m.r > 0 || m.MR > 0) series.push({ name: 'C_R', x: ts, y: c.map(q => q.R), color: 3 });
      if (m.s > 0) series.push({ name: 'C_S', x: ts, y: c.map(q => q.S), color: 4 });
      const half = Xf > 0.5 || v.mode === 'X' ? K.tBatch(m, Math.min(0.5, cap)) : NaN;
      const cyc = t + (+v.td || 0);
      const kpi = [
        { k: 't', l: 'Reaction time t', v: t, u: v._tu, main: v.mode === 't' },
        { k: 'X', l: 'Conversion X_A', v: Xf, main: v.mode === 'X' },
        { k: 'CA', l: 'C_A at end', v: m.conc(Xf).A, u: v._cu },
        { k: 'thalf', l: 'Half-life t_{1/2}', v: 0.5 < cap ? half : NaN, u: v._tu },
        { k: 'prod', l: 'Production per volume', v: m.CA0 * Xf / cyc, u: 'mol A/(' + KIT.volUnit(v._cu) + '·' + v._tu + ')', s: 'averaged over reaction + turnaround' }
      ];
      if (m.eps) kpi.push({ k: 'VV0', l: 'V/V_0 at end', v: 1 + m.eps * Xf });
      /* Levenspiel-type plot for the batch: CA0/((1+eps X)(-rA)) */
      const hi = Math.min(Math.max(Xf * 1.1, 0.5), cap);
      const lx = N.linspace(0, hi, 200), ly = lx.map(x => m.CA0 / ((1 + m.eps * x) * m.rate(x)));
      const sx = N.linspace(0, Xf, 90);
      const lev = {
        title: 'Area under the curve = reaction time',
        x: { label: 'X_A', min: 0, max: hi }, y: { label: 'C_{A0}/[(1+ε_AX_A)(−r_A)]', unit: v._tu, min: 0, max: Math.max(...ly.filter((q, i) => lx[i] <= Xf && isFinite(q))) * 1.4 },
        series: [{ name: 'integrand', x: lx, y: ly, color: 1 }, { x: sx, y: sx.map(x => m.CA0 / ((1 + m.eps * x) * m.rate(x))), type: 'area', fill: 0.22, color: 3, hover: false, legend: false, noRange: true }],
        marks: [{ t: 'v', x: Xf, label: 'X_A', color: 8, noRange: true }], legend: false
      };
      const tab = { title: 'Time course', cols: ['t', 'X_A', 'C_A'].concat(m.b > 0 ? ['C_B'] : [], m.r > 0 || m.MR > 0 ? ['C_R'] : [], m.eps ? ['V/V_0'] : []), rows: [] };
      N.linspace(0, t, 11).forEach(tt => {
        const X = N.interpClamp(ts, xs, tt), q = m.conc(X);
        tab.rows.push([tt, X, q.A].concat(m.b > 0 ? [q.B] : [], m.r > 0 || m.MR > 0 ? [q.R] : [], m.eps ? [1 + m.eps * X] : []));
      });
      const diag = D.svg(260, 120, D.tank(100, 22, 64, 74, { fill: Xf, tag: 'R-501', jacket: false, level: m.eps ? Math.min(0.95, 0.72 * (1 + m.eps * Xf)) : 0.78 }) + D.text(200, 60, 'X ' + N.fmt(Xf, 3), { cls: 'pid-r', anchor: 'start' }) + D.text(200, 76, 't ' + N.fmt(t, 3) + ' ' + v._tu, { cls: 'pid-r', anchor: 'start' }), 'batch reactor');
      return {
        kpi, warn, info, diagram: diag,
        plots: [{ title: 'Concentrations', x: { label: 't', unit: v._tu, min: 0 }, y: { label: 'C', unit: v._cu, min: 0 }, series }, lev],
        tables: [tab]
      };
    }
  });

  /* ---------------- mixed flow ---------------- */
  ST.tool({
    id: 'mfr', tag: 'H-502', ch: 5, sec: '5.2', title: 'Mixed flow reactor',
    lede: 'Size a steady-state mixed flow reactor for a target conversion, or find every steady state for a given space-time. Any rate law, liquids or expanding gases.',
    keys: 'cstr mfr backmix stirred tank steady state space time space velocity multiple steady states',
    units: { C: 'mol/L', t: 'min' },
    inputs: [
      { group: 'Question' },
      { k: 'mode', t: 'seg', label: 'Solve for', opts: [['tau', 'space-time for $X_A$'], ['X', '$X_A$ for space-time']], v: 'tau' },
      { k: 'Xf', label: 'Outlet conversion $X_{Af}$', v: 0.9, min: 0.001, max: 0.999, step: 0.001, lo: 0, hi: 1, if: v => v.mode === 'tau' },
      { k: 'tau', label: 'Space-time τ', unit: '[t]', v: 10, min: 1e-3, max: 1e5, scale: 'log', lo: 0, if: v => v.mode === 'X' },
      { k: 'Xi', label: 'Conversion entering $X_{Ai}$', v: 0, min: 0, max: 0.99, step: 0.01, lo: 0, hi: 1, hint: 'Non-zero when the feed was partly converted upstream (Eq. 5.12).' }
    ].concat(KIT.kin({ d: { k: 0.2 } }), KIT.flowInputs()),
    presets: [
      { name: 'Example 5.3', v: { mode: 'tau', Xf: 0.75, law: 'rev', k1: 25, a: 2, ob: 1, k2: 3, orr: 1, os: 0, CA0: 0.8, eps: 0, b: 0.5, MB: 1.75, r: 0.5, MR: 0, basis: 'none', _cu: 'mol/L', _tu: 'min' },
        book: { rate: 0.2 },
        note: 'Rewritten with B (the limiting reactant) as the key component: ½A + B → ½R, $-r_B$ = 25$C_AC_B^2$ − 3$C_R$. The printed solution then divides by 0.2 (= $C_B$) where the balance needs $C_{B0}-C_B$ = 0.6, and reports v = 6 L/min. The balance gives τ = 3 min, so a 6 L reactor takes 2 L/min in total, 1 L/min of each feed.' },
      { name: 'Problem 6.7 tank', v: { mode: 'tau', Xf: 0.99, law: 'sat', sk1: 1 / 0.2, sk2: 1 / 0.2, CA0: 1, basis: 'v0', v0: 10, _tu: 'min' }, note: '$-r_A$ = $C_A$/(0.2 + $C_A$) in one tank to 99% conversion, 10 L/min.' },
      { name: 'Autocatalytic, three states', v: { mode: 'X', tau: 1.2, law: 'custom', expr: 'k*CA/(1+K*CA)^2', cpar: 'k = 10, K = 2.2', CA0: 2, basis: 'none', Xi: 0 }, note: 'Substrate-inhibited kinetics give three steady states for this space-time.' }
    ],
    eqs: [
      ['\\frac{V}{F_{A0}} = \\frac{\\tau}{C_{A0}} = \\frac{X_{Af}-X_{Ai}}{(-r_A)_f}', 'Eq. 5.11–5.12'],
      ['\\tau = \\frac{V}{v_0} = \\frac{C_{A0}-C_A}{-r_A}\\qquad(\\varepsilon_A=0)', 'Eq. 5.13'],
      ['k\\tau = \\frac{X_A(1+\\varepsilon_AX_A)}{1-X_A}', 'Eq. 5.14b', 'First order, any expansion factor.'],
      ['C_A = \\frac{-1+\\sqrt{1+4k\\tau C_{A0}}}{2k\\tau}', 'Eq. 5.15', 'Second order, constant density.']
    ],
    method: ['For a given space-time every root of $C_{A0}(X_A-X_{Ai}) = \\tau(-r_A)$ is located by scanning [0, $X_{max}$] for sign changes and refining each with Brent’s method. Nonlinear kinetics (substrate inhibition, autocatalysis) can give several roots: each is a possible steady state, and the middle one of three is unstable.', 'The graphical solution plots the rate curve against the material-balance line of slope $C_{A0}/\\tau$; intersections are the operating points.'],
    run(v) {
      const m = KIT.model(v), warn = [], info = [];
      const cap = capX(m), Xi = +v.Xi || 0;
      let Xf, tau, roots = null;
      if (v.mode === 'tau') {
        Xf = Math.min(+v.Xf, cap);
        eqWarn(m, +v.Xf, warn);
        if (Xf <= Xi) throw new Error('Outlet conversion must exceed the entering conversion.');
        tau = K.tauMFR(m, Xf, Xi);
      } else {
        tau = +v.tau;
        roots = K.XMFRall(m, tau, Xi);
        Xf = roots[roots.length - 1];
        if (roots.length > 1) warn.push(roots.length + ' steady states satisfy the balance: $X_A$ = ' + roots.map(x => N.fmt(x, 4)).join(', ') + '. Which one you get depends on start-up; with three, the middle one is unstable.');
      }
      const kpi = [
        { k: 'tau', l: 'Space-time τ', v: tau, u: v._tu, main: v.mode === 'tau' && v.basis === 'none' },
        { k: 'X', l: 'Outlet X_A', v: Xf, main: v.mode === 'X' },
        { k: 'CA', l: 'Outlet C_A', v: m.conc(Xf).A, u: v._cu },
        { k: 'rate', l: '−r_A in reactor', v: m.rate(Xf), u: v._cu + '/' + v._tu },
        { k: 's', l: 'Space velocity s', v: 1 / tau, u: '1/' + v._tu }
      ];
      sizeKPIs(v, m, tau, kpi);
      if (v.law === 'power' && +v.ob === 0) {
        const n = +v.a, k = +v.k;
        if (Math.abs(n - 1) < 1e-12 && Xi === 0) kpi.push({ k: 'chk', l: 'Closed form, Eq. 5.14b', v: Xf * (1 + m.eps * Xf) / (1 - Xf) / k, u: v._tu, s: 'τ, agrees with the balance' });
        if (Math.abs(n - 2) < 1e-12 && !m.eps) kpi.push({ k: 'chk', l: 'Eq. 5.15 check', v: (-1 + Math.sqrt(1 + 4 * k * tau * m.CA0)) / (2 * k * tau), u: v._cu, s: 'C_A from closed form' });
      }
      const lev = KIT.levPlot(m, Xf, { tu: v._tu, shade: [{ type: 'mfr', a: Xi, b: Xf }], title: 'Rectangle = τ' });
      const xs = N.linspace(0, Math.min(cap, 1), 240);
      const gr = {
        title: 'Graphical solution',
        x: { label: 'X_A', min: 0, max: Math.min(1, cap) }, y: { label: '−r_A', unit: v._cu + '/' + v._tu, min: 0 },
        series: [{ name: 'rate curve', x: xs, y: xs.map(m.rate), color: 1 }, { name: 'balance line', x: [Xi, Math.min(1, cap)], y: [0, m.CA0 * (Math.min(1, cap) - Xi) / tau], color: 2, dash: '6 4', noRange: true }],
        marks: (roots || [Xf]).map(x => ({ t: 'pt', x, y: m.rate(x), color: 2, label: N.fmt(x, 3) }))
      };
      const tab = KIT.concTable(m, [['feed', 0], ...(Xi > 0 ? [['entering', Xi]] : []), ['outlet', Xf]], 'Feed and outlet');
      const diag = D.train([{ type: 'mfr', tag: 'R-502', X: Xf }], { h: 118 }).svg;
      return { kpi, warn, info, diagram: diag, plots: [lev, gr], tables: [tab] };
    }
  });

  /* ---------------- plug flow ---------------- */
  ST.tool({
    id: 'pfr', tag: 'H-503', ch: 5, sec: '5.3', title: 'Plug flow reactor',
    lede: 'Space-time, volume and holding time for a target conversion, or conversion for a given space-time, with profiles along the reactor. Any rate law, any constant expansion factor.',
    keys: 'pfr tubular plug flow space time holding time residence time expansion gas',
    units: { C: 'mol/L', t: 'min' },
    inputs: [
      { group: 'Question' },
      { k: 'mode', t: 'seg', label: 'Solve for', opts: [['tau', 'space-time for $X_A$'], ['X', '$X_A$ for space-time']], v: 'tau' },
      { k: 'Xf', label: 'Outlet conversion $X_{Af}$', v: 0.9, min: 0.001, max: 0.999, step: 0.001, lo: 0, hi: 1, if: v => v.mode === 'tau' },
      { k: 'tau', label: 'Space-time τ', unit: '[t]', v: 10, min: 1e-3, max: 1e5, scale: 'log', lo: 0, if: v => v.mode === 'X' },
      { k: 'Xi', label: 'Conversion entering $X_{Ai}$', v: 0, min: 0, max: 0.99, step: 0.01, lo: 0, hi: 1 }
    ].concat(KIT.kin({ d: { k: 0.2 } }), KIT.flowInputs()),
    presets: [
      { name: 'Example 5.4', v: { mode: 'tau', Xf: 0.8, law: 'power', k: 0.01, a: 0.5, ob: 0, CA0: 0.0625, eps: 1, basis: 'none', Xi: 0, _cu: 'mol/L', _tu: 's' }, book: { tau: 33.2 }, note: 'A → 3R with 50% inerts, so $\\varepsilon_A$ = 1.' },
      { name: 'Example 5.5', v: { mode: 'tau', Xf: 0.8, law: 'power', k: 10, a: 1, ob: 0, CA0: 60, eps: 0.75, basis: 'FA0', FA0: 40, Xi: 0, _cu: 'mol/m³', _tu: 'h' }, book: { V: 0.148 }, note: 'Phosphine decomposition 4PH₃ → P₄ + 6H₂ at 649 °C and 460 kPa; 40 mol/h of feed. Volume in m³ (148 L).' }
    ],
    eqs: [
      ['\\frac{V}{F_{A0}} = \\frac{\\tau}{C_{A0}} = \\int_{X_{Ai}}^{X_{Af}}\\frac{dX_A}{-r_A}', 'Eq. 5.17–5.18'],
      ['k\\tau = C_{A0}X_A', 'Eq. 5.20', 'Zero order, any $\\varepsilon_A$.'],
      ['k\\tau = -(1+\\varepsilon_A)\\ln(1-X_A) - \\varepsilon_AX_A', 'Eq. 5.21', 'First order.'],
      ['C_{A0}k\\tau = 2\\varepsilon_A(1+\\varepsilon_A)\\ln(1-X_A) + \\varepsilon_A^2X_A + (\\varepsilon_A+1)^2\\frac{X_A}{1-X_A}', 'Eq. 5.23', 'Second order, 2A → products or equimolar A + B.'],
      ['\\bar t = C_{A0}\\int_0^{X_A}\\frac{dX_A}{(-r_A)(1+\\varepsilon_AX_A)}', 'Eq. 5.24', 'Holding time. Equals τ only when the density is constant.']
    ],
    method: ['τ is the integral of the Levenspiel curve, evaluated by adaptive quadrature. The profile along the reactor comes from integrating $dX_A/d\\tau = (-r_A)/C_{A0}$, so position $V/V_{total}$ maps onto τ.', 'For a gas whose moles change with reaction the holding time $\\bar t$ differs from τ; the tool reports both, as in the popcorn-popper discussion of §5.3.'],
    run(v) {
      const m = KIT.model(v), warn = [], info = [];
      const cap = capX(m), Xi = +v.Xi || 0;
      let Xf, tau;
      if (v.mode === 'tau') {
        Xf = Math.min(+v.Xf, cap);
        eqWarn(m, +v.Xf, warn);
        if (Xf <= Xi) throw new Error('Outlet conversion must exceed the entering conversion.');
        tau = K.tauPFR(m, Xf, Xi);
      } else { tau = +v.tau; Xf = K.XPFR(m, tau, Xi).X; }
      if (m.rate(Xi) <= 0) warn.push('The rate is zero at the inlet, so a plug flow reactor cannot start this reaction. Seed with product or use recycle.');
      const tbar = m.CA0 * N.quad(x => 1 / ((m.rate(x)) * (1 + m.eps * x)), Xi, Xf);
      const kpi = [
        { k: 'tau', l: 'Space-time τ', v: tau, u: v._tu, main: v.mode === 'tau' && v.basis === 'none' },
        { k: 'X', l: 'Outlet X_A', v: Xf, main: v.mode === 'X' },
        { k: 'CA', l: 'Outlet C_A', v: m.conc(Xf).A, u: v._cu },
        { k: 'tbar', l: 'Holding time t̄', v: tbar, u: v._tu, s: m.eps ? 'differs from τ: density changes' : 'equal to τ' }
      ];
      sizeKPIs(v, m, tau, kpi);
      if (v.law === 'power' && +v.ob === 0 && Xi === 0) {
        const n = +v.a, k = +v.k, e = m.eps;
        let tc = NaN, eq = '';
        if (Math.abs(n) < 1e-12) { tc = m.CA0 * Xf / k; eq = '5.20'; }
        if (Math.abs(n - 1) < 1e-12) { tc = (-(1 + e) * Math.log(1 - Xf) - e * Xf) / k; eq = '5.21'; }
        if (Math.abs(n - 2) < 1e-12) { tc = (2 * e * (1 + e) * Math.log(1 - Xf) + e * e * Xf + (e + 1) ** 2 * Xf / (1 - Xf)) / (k * m.CA0); eq = '5.23'; }
        if (eq) kpi.push({ k: 'chk', l: 'Closed form, Eq. ' + eq, v: tc, u: v._tu, s: 'agrees with the integral' });
      }
      const pr = K.XPFR(m, tau, Xi, { n: 120 });
      const pos = pr.prof.t.map(t => t / tau);
      const cs = pr.prof.x.map(m.conc);
      const series = [{ name: 'X_A', x: pos, y: pr.prof.x, color: 1 }];
      const cser = [{ name: 'C_A', x: pos, y: cs.map(c => c.A), color: 1 }];
      if (m.b > 0) cser.push({ name: 'C_B', x: pos, y: cs.map(c => c.B), color: 2 });
      if (m.r > 0 || m.MR > 0) cser.push({ name: 'C_R', x: pos, y: cs.map(c => c.R), color: 3 });
      const lev = KIT.levPlot(m, Xf, { tu: v._tu, shade: [{ type: 'pfr', a: Xi, b: Xf }], title: 'Area = τ' });
      const tab = { title: 'Profile along the reactor', cols: ['V/V_{total}', 'τ', 'X_A', 'C_A'], rows: N.linspace(0, 1, 11).map(s => { const X = N.interpClamp(pos, pr.prof.x, s); return [s, s * tau, X, m.conc(X).A]; }) };
      const diag = D.train([{ type: 'pfr', tag: 'R-503', X: Xf, prof: N.linspace(0, 1, 20).map(s => N.interpClamp(pos, pr.prof.x, s)) }], { h: 104, pl: 260 }).svg;
      return {
        kpi, warn, info, diagram: diag,
        plots: [lev, { title: 'Profiles', x: { label: 'V/V_{total}', min: 0, max: 1 }, y: { label: 'C', unit: v._cu, min: 0 }, series: cser }],
        tables: [tab]
      };
    }
  });

  /* ---------------- rate from flow reactor data ---------------- */
  ST.tool({
    id: 'flowdata', tag: 'H-504', ch: 5, sec: '5.2–5.3', title: 'Rate equation from flow reactor runs',
    lede: 'Turn steady-state runs from a mixed flow reactor (rates read directly) or from a plug flow reactor (integral data) into a fitted rate equation, with expansion accounted for.',
    keys: 'kinetics from mixed flow reactor data plug flow integral analysis order fit log log expansion factor',
    units: { C: 'mol/L', t: 'min' },
    inputs: [
      { group: 'Data' },
      { k: 'kind', t: 'seg', label: 'Reactor', opts: [['mfr', 'Mixed flow'], ['pfr', 'Plug flow']], v: 'mfr' },
      { k: 'inp', t: 'select', label: 'Each run gives', opts: [['tauC', 'τ and outlet C_A'], ['tauX', 'τ and outlet X_A'], ['vC', 'feed rate v₀ and outlet C_A'], ['tau3', 'τ, its own feed C_A0, outlet C_A']], v: 'vC' },
      { k: 'V', label: 'Reactor volume V', unit: v => KIT.volUnit(v._cu), v: 0.1, min: 1e-6, max: 1e6, scale: 'log', if: v => v.inp === 'vC' },
      { k: 'CA0', label: 'Feed $C_{A0}$', unit: '[C]', v: 100, min: 1e-6, max: 1e6, scale: 'log', if: v => v.inp !== 'tau3' },
      { k: 'eps', label: 'Expansion factor $\\varepsilon_A$', v: -0.5, min: -0.95, max: 4, step: 0.01 },
      { k: 'runs', t: 'table', label: 'Runs', cols: ['τ or v₀', 'C_A or X_A'], v: [[10, 85.7], [3, 66.7], [1.2, 50], [0.5, 33.4]], if: v => v.inp !== 'tau3' },
      { k: 'runs3', t: 'table', label: 'Runs', cols: ['τ', 'C_{A0}', 'C_A'], v: [[60, 50, 20], [35, 100, 40], [11, 100, 60], [20, 200, 80], [11, 200, 100]], if: v => v.inp === 'tau3' },
      { group: 'Model' },
      { k: 'model', t: 'select', label: 'Fit', opts: [['n', 'nth order, k C_A^n'], ['sat', 'saturation, k₁C_A/(1+k₂C_A)'], ['n1', 'first order'], ['n2', 'second order']], v: 'n' }
    ],
    presets: [
      { name: 'Example 5.2', v: { kind: 'mfr', inp: 'vC', V: 0.1, CA0: 100, eps: -0.5, runs: [[10, 85.7], [3, 66.7], [1.2, 50], [0.5, 33.4]], model: 'n', _cu: 'mmol/L', _tu: 'h' }, book: { n: { v: 2, tol: 0.05 } }, note: 'Dimerization 2A → R of pure gaseous A ($\\varepsilon_A$ = −0.5). The book reads the slope as 2 and then fixes k = 0.36 L/(mmol·h) from the lowest point; with n held at 2 a least-squares k over all four runs is 0.34.' },
      { name: 'Example 5.2, ignoring expansion', v: { kind: 'mfr', inp: 'vC', V: 0.1, CA0: 100, eps: 0, runs: [[10, 85.7], [3, 66.7], [1.2, 50], [0.5, 33.4]], model: 'n', _cu: 'mmol/L', _tu: 'h' }, note: 'Setting $\\varepsilon_A$ = 0 gives the wrong order, about 1.6, as the book warns.' },
      { name: 'Problem 6.8', v: { kind: 'mfr', inp: 'tau3', eps: 0, runs3: [[60, 50, 20], [35, 100, 40], [11, 100, 60], [20, 200, 80], [11, 200, 100]], model: 'n', _cu: 'mmol/L', _tu: 's' }, note: 'Five mixed flow runs with different feeds. Compare an nth-order fit with the saturation form.' }
    ],
    eqs: [
      ['-r_A = \\frac{C_{A0}X_A}{\\tau},\\qquad X_A = \\frac{1-C_A/C_{A0}}{1+\\varepsilon_AC_A/C_{A0}}', 'Eqs. 5.11, 4.5', 'Each mixed flow run gives one rate at the outlet composition.'],
      ['\\log(-r_A) = \\log k + n\\log C_A', 'Fig. E5.2'],
      ['\\tau = C_{A0}\\int_0^{X_A}\\frac{dX_A}{-r_A(X_A)}', 'Eq. 5.17', 'Plug flow runs are fitted by integrating the trial rate law for every run (nonlinear least squares).']
    ],
    run(v) {
      const three = v.inp === 'tau3';
      const rows = (three ? v.runs3 : v.runs || []).filter(r => r.slice(0, three ? 3 : 2).every(isFinite) && r[0] > 0);
      if (rows.length < 2) throw new Error('Enter at least two runs.');
      const e = +v.eps;
      const pts = rows.map(r => {
        const tau = v.inp === 'vC' ? (+v.V) / r[0] : r[0];
        const C0 = three ? r[1] : +v.CA0;
        let X, CA;
        if (v.inp === 'tauX') { X = r[1]; CA = C0 * (1 - X) / (1 + e * X); }
        else { CA = three ? r[2] : r[1]; X = (1 - CA / C0) / (1 + e * CA / C0); }
        return { tau, X, CA, C0 };
      });
      const CA0 = three ? pts[0].C0 : +v.CA0;
      const kpi = [], warn = [], tables = [], plots = [];
      if (v.kind === 'mfr') {
        pts.forEach(p => { p.r = p.C0 * p.X / p.tau; });
        const good = pts.filter(p => p.r > 0 && p.CA > 0);
        let fit, fitR;
        if (v.model === 'sat') {
          const f = N.lm(q => good.map(p => Math.log(q[0] * p.CA / (1 + q[1] * p.CA)) - Math.log(p.r)), [good[0].r / good[0].CA, 0.01], { pos: [true, true] });
          fit = { k1: f.p[0], k2: f.p[1], se: f.se };
          fitR = c => fit.k1 * c / (1 + fit.k2 * c);
          kpi.push({ k: 'k1', l: 'k_1', v: fit.k1, u: '1/' + v._tu, main: true }, { k: 'k2', l: 'k_2', v: fit.k2, u: '(' + v._cu + ')⁻¹' });
        } else {
          const lr = N.linreg(good.map(p => Math.log10(p.CA)), good.map(p => Math.log10(p.r)));
          let n = lr.b, k = Math.pow(10, lr.a);
          if (v.model === 'n1' || v.model === 'n2') {
            n = v.model === 'n1' ? 1 : 2;
            k = Math.pow(10, N.mean(good.map(p => Math.log10(p.r) - n * Math.log10(p.CA))));
          }
          fitR = c => k * Math.pow(c, n);
          kpi.push({ k: 'n', l: 'Order n', v: n, main: true, s: v.model === 'n' ? '± ' + N.fmt(lr.sb * N.tinv(0.975, Math.max(1, lr.dof)), 2) + ' (95%)' : 'fixed' }, { k: 'k', l: 'k', v: k, u: KIT.kUnit(n, v._cu, v._tu) }, { k: 'r2', l: 'R² (log–log)', v: lr.r2 });
        }
        tables.push({ title: 'Rates from the runs', cols: ['τ', 'C_{A0}', 'C_A', 'X_A', '−r_A', 'fitted −r_A'], rows: pts.map(p => [p.tau, p.C0, p.CA, p.X, p.r, fitR(p.CA)]) });
        const cs = N.logspace(Math.min(...good.map(p => p.CA)) * 0.8, Math.max(...good.map(p => p.CA)) * 1.2, 60);
        plots.push({ title: 'Rate against concentration', x: { label: 'C_A', unit: v._cu, log: true }, y: { label: '−r_A', unit: v._cu + '/' + v._tu, log: true }, series: [{ name: 'runs', type: 'pts', x: good.map(p => p.CA), y: good.map(p => p.r), color: 2 }, { name: 'fit', x: cs, y: cs.map(fitR), color: 1 }] });
      } else {
        const tauModel = (rate, X, C0 = CA0) => C0 * N.quad(x => 1 / rate(C0 * (1 - x) / (1 + e * x)), 0, X);
        let p0, rateOf, names;
        if (v.model === 'sat') { p0 = [1, 0.1]; rateOf = q => c => q[0] * c / (1 + q[1] * c); names = ['k_1', 'k_2']; }
        else if (v.model === 'n1') { p0 = [0.1]; rateOf = q => c => q[0] * c; names = ['k']; }
        else if (v.model === 'n2') { p0 = [0.1]; rateOf = q => c => q[0] * c * c; names = ['k']; }
        else { p0 = [0.1, 1]; rateOf = q => c => q[0] * Math.pow(Math.max(c, 1e-300), q[1]); names = ['k', 'n']; }
        const k1g = -Math.log(Math.max(1e-6, 1 - pts[pts.length - 1].X)) / pts[pts.length - 1].tau;
        p0[0] = v.model === 'n2' ? k1g / CA0 : k1g * Math.pow(CA0, v.model === 'n' ? 0 : 0);
        const res = N.lm(q => pts.map(p => (tauModel(rateOf(q), p.X, p.C0) - p.tau) / p.tau), p0, { pos: p0.map((_, i) => !(names[i] === 'n')) });
        const q = res.p;
        names.forEach((nm, i) => kpi.push({ k: nm.replace(/[^a-z0-9]/gi, ''), l: nm, v: q[i], main: i === 0, s: isFinite(res.se[i]) ? '± ' + N.fmt(1.96 * res.se[i], 2) : '' }));
        const fr = rateOf(q);
        tables.push({ title: 'Runs and fit', cols: ['τ', 'C_{A0}', 'X_A', 'τ from fit'], rows: pts.map(p => [p.tau, p.C0, p.X, tauModel(fr, p.X, p.C0)]) });
        const xs = N.linspace(0, Math.max(...pts.map(p => p.X)) * 1.05, 60);
        plots.push({ title: 'Integral test', x: { label: 'X_A', min: 0 }, y: { label: 'τ', unit: v._tu, min: 0 }, series: [{ name: 'runs', type: 'pts', x: pts.map(p => p.X), y: pts.map(p => p.tau), color: 2 }, { name: 'fitted model', x: xs, y: xs.map(x => tauModel(fr, x)), color: 1 }] });
      }
      return { kpi, warn, plots, tables };
    }
  });

  /* ---------------- rates of several species from one MFR run ---------------- */
  ST.tool({
    id: 'mfrrates', tag: 'H-505', ch: 5, sec: '5.2', title: 'Rates of every component in a mixed flow run',
    lede: 'One steady-state run in a mixed flow reactor gives the rate of formation or disappearance of every measured component, even when the stoichiometry is unknown.',
    keys: 'mixed flow rate of reaction several components complex stoichiometry example 5.1',
    units: { C: 'mol/L', t: 'min' },
    inputs: [
      { group: 'Run' },
      { k: 'V', label: 'Reactor volume V', unit: v => KIT.volUnit(v._cu), v: 1, min: 1e-6, max: 1e6, scale: 'log' },
      { k: 'v0', label: 'Flow rate v', unit: v => KIT.volUnit(v._cu) + '/' + v._tu, v: 1, min: 1e-6, max: 1e6, scale: 'log', hint: 'Liquid: the same flow in and out.' },
      { k: 'sp', t: 'table', label: 'Components (feed, outlet)', cols: ['C_{in}', 'C_{out}'], v: [[0.1, 0.02], [0.01, 0.03], [0, 0.04]] },
      { k: 'names', t: 'text', label: 'Names, in table order', v: 'A, B, C' }
    ],
    presets: [{ name: 'Example 5.1', v: { V: 1, v0: 1, sp: [[0.1, 0.02], [0.01, 0.03], [0, 0.04]], names: 'A, B, C', _cu: 'mol/L', _tu: 'min' }, book: { r0: 0.08, r1: -0.02, r2: -0.04 }, note: 'Rates of disappearance: A 0.08, B −0.02, C −0.04 mol/L·min (B and C are being formed).' }],
    eqs: [['-r_i = \\frac{C_{i0}-C_i}{\\tau},\\qquad \\tau = V/v', 'Eq. 5.13']],
    run(v) {
      const tau = +v.V / +v.v0, names = String(v.names || '').split(',').map(s => s.trim());
      const rows = (v.sp || []).filter(r => isFinite(r[0]) && isFinite(r[1]));
      const kpi = [{ k: 'tau', l: 'Space-time τ', v: tau, u: v._tu }];
      rows.forEach((r, i) => kpi.push({ k: 'r' + i, l: '−r_{' + (names[i] || 'species ' + (i + 1)) + '}', v: (r[0] - r[1]) / tau, u: v._cu + '/' + v._tu, s: r[0] - r[1] >= 0 ? 'disappearing' : 'being formed' }));
      return { kpi, tables: [{ title: 'Component rates', cols: ['component', 'C_{in}', 'C_{out}', '−r_i', 'r_i'], rows: rows.map((r, i) => [names[i] || String(i + 1), r[0], r[1], (r[0] - r[1]) / tau, (r[1] - r[0]) / tau]) }] };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
