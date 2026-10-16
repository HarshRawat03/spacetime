/* Chapter 24 — fluid–fluid reactors: design. Towers (countercurrent, cocurrent), tanks (mixed/mixed,
   plug gas through mixed liquid), and batches of liquid under mixed or plug flow of gas; straight
   absorption, every reaction regime, exhaustion of B inside the contactor, pure or diluted gas. */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, FF = ST.FF, D = ST.D;
  const fmt = (x, s = 4) => N.fmt(x, s);

  ST.tool({
    id: 'gldesign', tag: 'N-2401', ch: 24, sec: 'Ch. 24', title: 'Gas–liquid reactors: towers, tanks and batches of liquid',
    lede: 'Size a gas–liquid contactor for absorption with or without reaction. Countercurrent and cocurrent towers, agitated tanks with both phases mixed, bubble tanks with gas in plug flow, and batches of liquid under mixed or plug flow of gas — each with the material balance of the book and the rate of Eq. 23.5 re-evaluated at every point, so the regime may change along the contactor and B may run out inside it. Height or volume for a given duty, or the outlet for a given size, the batch time with its efficiency of A use, and the profiles along the way.',
    keys: 'gas liquid reactor design absorber tower countercurrent cocurrent packed column agitated tank bubble tank semibatch batch liquid height of tower straight mass transfer absorption with reaction instantaneous fast slow reaction zones B exhausted material balance operating line example 24.1 24.2 24.3 24.4 24.5 24.6 problem 24.2 24.3 24.4 24.5 24.6 24.7 24.8 24.9 24.10 24.11 24.12 24.13 24.14 24.15 24.16 chlorine caustic CO2 carbonate buffer',
    units: { t: 'h' }, slow: true,
    inputs: [
      { group: 'Contactor' },
      { k: 'pat', t: 'select', label: 'Contacting pattern (Fig. 24.2)', opts: [['cc', 'tower: plug G / plug L, countercurrent'], ['co', 'tower: plug G / plug L, cocurrent'], ['mm', 'agitated tank: mixed G / mixed L'], ['pm', 'bubble tank: plug G / mixed L'], ['mb', 'mixed G / batch of L'], ['pb', 'plug G / batch of L']], v: 'cc' },
      { k: 'find', t: 'seg', label: 'Find', opts: [['V', 'size for the duty'], ['out', 'outlet for a given size']], v: 'V', if: v => v.pat !== 'mb' && v.pat !== 'pb' },
      { k: 'Vr', label: 'Contactor volume V_r', unit: 'm³', v: 1, min: 1e-9, max: 1e9, scale: 'log', if: v => v.find === 'out' || v.pat === 'mb' || v.pat === 'pb' },
      { k: 'Acs', label: 'Cross-section (to report a height)', unit: 'm²', v: 1, min: 1e-9, max: 1e6, scale: 'log', if: v => v.pat === 'cc' || v.pat === 'co' },
      { group: 'Gas' },
      { k: 'pure', t: 'check', label: 'Gas is pure A (p_A = π throughout)', v: false },
      { k: 'pi', label: 'Total pressure π', unit: 'Pa', v: 1e5, min: 1, max: 1e9, scale: 'log' },
      { k: 'dil', t: 'check', label: 'Dilute-gas balance, p_A ≪ π (the book’s Eqs. 24.10, 24.12, 24.16, 24.18)', v: false, if: v => !v.pure },
      { k: 'Fg', label: 'Inert gas flow F_g', unit: 'mol/[t]', v: 1e5, min: 1e-12, max: 1e12, scale: 'log', if: v => !v.pure },
      { k: 'pin', label: 'p_A entering', unit: 'Pa', v: 100, min: 0, max: 1e9, if: v => !v.pure },
      { k: 'pout', label: 'p_A leaving (target)', unit: 'Pa', v: 20, min: 0, max: 1e9, if: v => !v.pure && v.find !== 'out' && v.pat !== 'mb' && v.pat !== 'pb' },
      { k: 'FA0', label: 'Feed of pure A (for the fraction absorbed)', unit: 'mol/[t]', v: 0, min: 0, max: 1e12, if: v => v.pure },
      { group: 'Liquid' },
      { k: 'Fl', label: 'Inert liquid flow F_l', unit: 'mol/[t]', v: 7e5, min: 1e-12, max: 1e14, scale: 'log', if: v => v.pat !== 'mb' && v.pat !== 'pb' },
      { k: 'CT', label: 'Molar density of the liquid C_T ≈ C_U', unit: 'mol/m³', v: 56000, min: 1, max: 1e6 },
      { k: 'CBin', label: v => (v.pat === 'mb' || v.pat === 'pb' ? 'C_B at the start, C_B0' : 'C_B entering'), unit: 'mol/m³', v: 0, min: 0, max: 1e7 },
      { k: 'CBf', label: 'C_B at the end of the batch', unit: 'mol/m³', v: 55.6, min: 0, max: 1e7, if: v => v.pat === 'mb' || v.pat === 'pb' },
      { k: 'CAin', label: 'Dissolved A entering with the liquid', unit: 'mol/m³', v: 0, min: 0, max: 1e7, if: v => v.pat !== 'mb' && v.pat !== 'pb' },
      { k: 'hold', t: 'check', label: 'Liquid recirculated: composition held at C_B entering', v: false, if: v => v.pat !== 'mb' && v.pat !== 'pb' }
    ].concat(FF.inputs({ k: 0, b: 1, kAga: 0.32, kAla: 0.1, a: 100, fl: 0.1, H: 12.5, DA: 3.6e-6, rB: 1 }).map(f => (f.k === 'theory' ? Object.assign({}, f, { opts: f.opts.slice(0, 2) }) : f))),
    presets: [
      { name: 'Example 24.1 — tower for straight absorption', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 1e5, pin: 100, pout: 20, Fl: 7e5, CT: 56000, CBin: 0, CAin: 0, hold: false, inst: false, k: 0, b: 1, gf: true, kAga: 0.32, kAla: 0.1, a: 100, fl: 0.1, H: 12.5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' }, book: { h: { v: 512.5, tol: 0.01 }, CAout: { v: 6.4, tol: 0.01 }, shL: { v: 0.976, tol: 0.003 } }, note: 'Impurity in air reduced from 100 to 20 Pa by pure water, F_g/A = 10⁵ and F_l/A = 7 × 10⁵ mol/hr·m², k_Ag a = 0.32 mol/hr·m³·Pa, k_Al a = 0.1 hr⁻¹, H_A = 12.5 Pa·m³/mol. Book: 97.6% of the resistance in the liquid film and a 512.5 m tower.' },
      { name: 'Example 24.2 — high concentration of B (800)', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 1e5, pin: 100, pout: 20, Fl: 7e5, CT: 56000, CBin: 800, CAin: 0, hold: false, inst: true, b: 1, gf: true, kAga: 0.32, kAla: 0.1, a: 100, fl: 0.1, H: 12.5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' }, book: { h: { v: 5.03, tol: 0.01 }, CBout: { v: 793.6, tol: 0.001 } }, note: 'Instantaneous reaction with plenty of B: the gas film alone controls everywhere (Case B), 5.03 m instead of 512.5 m. Example 24.5 reaches the same height with the general rate expression, which is what is done here.' },
      { name: 'Example 24.3 — low concentration of B (32)', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 1e5, pin: 100, pout: 20, Fl: 7e5, CT: 56000, CBin: 32, CAin: 0, hold: false, inst: true, b: 1, gf: true, kAga: 0.32, kAla: 0.1, a: 100, fl: 0.1, H: 12.5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' }, book: { h: { v: 24.4, tol: 0.01 }, CBout: { v: 25.6, tol: 0.001 } }, note: 'The reaction plane lies inside the liquid film throughout (Case A): 24.4 m.' },
      { name: 'Example 24.4 — intermediate B (128): two zones', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 1e5, pin: 100, pout: 20, Fl: 7e5, CT: 56000, CBin: 128, CAin: 0, hold: false, inst: true, b: 1, gf: true, kAga: 0.32, kAla: 0.1, a: 100, fl: 0.1, H: 12.5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' }, book: { h: { v: 6.91, tol: 0.01 }, pswitch: { v: 39.5, tol: 0.01 } }, note: 'Case B near the top, Case A near the bottom; the plane reaches the interface where p_A = 39.5 Pa. Height 6.91 m.' },
      { name: 'Example 24.6 — batch of liquid, gas in mixed flow', v: { pat: 'mb', Vr: 1.62, pure: false, dil: true, pi: 1e5, Fg: 9000, pin: 1000, CT: 55555.6, CBin: 555.6, CBf: 55.6, inst: false, k: 2.6e5, b: 1, gf: true, kAga: 0.72, kAla: 144, a: 100, fl: 0.9, H: 1e3, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' }, book: { t: { v: 9.13, tol: 0.02 }, r0: { v: 50.44, tol: 0.02 }, rf: { v: 48.13, tol: 0.03 } }, note: 'B from 555.6 to 55.6 mol/m³ by bubbling 9000 mol/hr of gas at 1000 Pa of A. The example states V_l = 1.62 m³ and f_l = 0.9, i.e. V_r = 1.8 m³, but its arithmetic uses 1.62 m³ for V_r; that is followed here, with the book’s dilute-gas balance, to reproduce 9.13 hr (with V_r = 1.8 m³ it is about 9.8 hr; with exact gas ratios 1–2% less).' },
      { name: 'Problem 24.2 — very soluble gas, tank (H = 0)', v: { pat: 'mm', find: 'V', pure: false, pi: 1e5, Fg: 90000, pin: 1000, pout: 100, Fl: 9000, CT: 55556, CBin: 5556, CAin: 0, hold: false, inst: false, k: 0, b: 1, gf: true, kAga: 0.72, kAla: 144, a: 200, fl: 0.9, H: 0, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.3 — straight absorption, tower (H = 18)', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 90000, pin: 1000, pout: 100, Fl: 900000, CT: 55556, CBin: 55.56, CAin: 0, hold: false, inst: false, k: 0, b: 1, gf: true, kAga: 0.36, kAla: 72, a: 100, fl: 0.08, H: 18, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.4 — straight absorption, tower (H = 1.8)', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 90000, pin: 1000, pout: 100, Fl: 900000, CT: 55556, CBin: 55.56, CAin: 0, hold: false, inst: false, k: 0, b: 1, gf: true, kAga: 0.36, kAla: 72, a: 100, fl: 0.08, H: 1.8, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.5 — instantaneous, tower (H = 10⁵)', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 90000, pin: 1000, pout: 100, Fl: 900000, CT: 55556, CBin: 55.56, CAin: 0, hold: false, inst: true, b: 1, gf: true, kAga: 0.36, kAla: 72, a: 100, fl: 0.08, H: 1e5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.6 — fast reaction, tank (k = 2.6 × 10⁷)', v: { pat: 'mm', find: 'V', pure: false, pi: 1e5, Fg: 90000, pin: 1000, pout: 100, Fl: 9000, CT: 55556, CBin: 5556, CAin: 0, hold: false, inst: false, k: 2.6e7, b: 1, gf: true, kAga: 0.72, kAla: 144, a: 200, fl: 0.9, H: 1e5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.7 — tank (k = 2.6 × 10⁵)', v: { pat: 'mm', find: 'V', pure: false, pi: 1e5, Fg: 90000, pin: 1000, pout: 100, Fl: 9000, CT: 55556, CBin: 5556, CAin: 0, hold: false, inst: false, k: 2.6e5, b: 1, gf: true, kAga: 0.72, kAla: 144, a: 200, fl: 0.9, H: 1e5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.8 — tower (H = 10³, k = 2.6 × 10³)', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 90000, pin: 1000, pout: 100, Fl: 900000, CT: 55556, CBin: 55.56, CAin: 0, hold: false, inst: false, k: 2.6e3, b: 1, gf: true, kAga: 0.36, kAla: 72, a: 100, fl: 0.08, H: 1e3, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.9 — tower (H = 10⁵, k = 2.6 × 10⁷)', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 90000, pin: 1000, pout: 100, Fl: 900000, CT: 55556, CBin: 55.56, CAin: 0, hold: false, inst: false, k: 2.6e7, b: 1, gf: true, kAga: 0.36, kAla: 72, a: 100, fl: 0.08, H: 1e5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.10 — tower (H = 10³, k = 2.6 × 10⁵)', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 1e5, Fg: 90000, pin: 1000, pout: 100, Fl: 900000, CT: 55556, CBin: 55.56, CAin: 0, hold: false, inst: false, k: 2.6e5, b: 1, gf: true, kAga: 0.36, kAla: 72, a: 100, fl: 0.08, H: 1e3, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.11 — pure CO₂ into a carbonate buffer', v: { pat: 'pm', find: 'out', Vr: 0.6041, pure: true, pi: 101325, FA0: 101325 * 0.0363 / (8.314462618 * 293.15), Fl: 1, CT: 55556, CBin: 300, CAin: 0, hold: true, inst: false, k: 0.433, b: 1, gf: false, kAla: 0.025, a: 120, fl: 0.08, H: 3500, DA: 1.4e-9, rB: 1, theory: 'vkh', _tu: 's' }, note: 'Pure CO₂ at 1 atm and 20 °C (0.0363 m³/s) through 0.6041 m³ of packing irrigated with recirculated buffer at C_B ≈ 300 mol/m³: no gas-film resistance and uniform liquid, so the rate is the same everywhere.' },
      { name: 'Problem 24.12 — chlorine into caustic, 99% removal', v: { pat: 'cc', find: 'V', Acs: 1, pure: false, pi: 101325, Fg: 100 * 0.9764, pin: 0.0236 * 101325, pout: 0.000236 * 101325, Fl: 250, CT: 55556, CBin: 2736, CAin: 0, hold: false, inst: true, b: 2, gf: true, kAga: 133 / 101325, kAla: 45 / 3600, a: 55, fl: 0.1, H: 125e-6 * 101325, DA: 1.5e-9, rB: 1, theory: 'vkh', _tu: 's' }, note: 'G = 100 mol/s·m² with 2.36% Cl₂, L = 250 mol/s·m² of 10% NaOH (2736 mol/m³), H_A = 125 × 10⁻⁶ atm·m³/mol. k_g a = 133 is read as mol/s·m³·atm (per hour would give a tower kilometres high); k_l a = 45 hr⁻¹. The gas film controls; the book suggests doubling the height for deviations from plug flow.' },
      { name: 'Problem 24.13 — Example 24.6 with k = 2.6 × 10⁵, H = 10⁵', v: { pat: 'mb', Vr: 1.62, pure: false, dil: true, pi: 1e5, Fg: 9000, pin: 1000, CT: 55555.6, CBin: 555.6, CBf: 55.6, inst: false, k: 2.6e5, b: 1, gf: true, kAga: 0.72, kAla: 144, a: 100, fl: 0.9, H: 1e5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.14 — k = 2.6 × 10⁹, H = 10⁵', v: { pat: 'mb', Vr: 1.62, pure: false, dil: true, pi: 1e5, Fg: 9000, pin: 1000, CT: 55555.6, CBin: 555.6, CBf: 55.6, inst: false, k: 2.6e9, b: 1, gf: true, kAga: 0.72, kAla: 144, a: 100, fl: 0.9, H: 1e5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.15 — k = 2.6 × 10³, H = 10⁵', v: { pat: 'mb', Vr: 1.62, pure: false, dil: true, pi: 1e5, Fg: 9000, pin: 1000, CT: 55555.6, CBin: 555.6, CBf: 55.6, inst: false, k: 2.6e3, b: 1, gf: true, kAga: 0.72, kAla: 144, a: 100, fl: 0.9, H: 1e5, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 24.16 — k = 2.6 × 10¹¹, H = 10³', v: { pat: 'mb', Vr: 1.62, pure: false, dil: true, pi: 1e5, Fg: 9000, pin: 1000, CT: 55555.6, CBin: 555.6, CBf: 55.6, inst: false, k: 2.6e11, b: 1, gf: true, kAga: 0.72, kAla: 144, a: 100, fl: 0.9, H: 1e3, DA: 3.6e-6, rB: 1, theory: 'vkh', _tu: 'h' } }
    ],
    eqs: [
      ['F_g\\,d\\mathbf{Y}_A = -\\frac{F_l}{b}\\,d\\mathbf{X}_B = (-r_A^{\\prime\\prime\\prime})\\,dV_r,\\qquad \\mathbf{Y}_A = \\frac{p_A}{\\pi - p_A}', 'Eqs. 24.7, 24.9'],
      ['V_r = F_g\\int_{\\mathbf{Y}_{A1}}^{\\mathbf{Y}_{A2}}\\frac{d\\mathbf{Y}_A}{-r_A^{\\prime\\prime\\prime}},\\qquad F_g(\\mathbf{Y}_{A2}-\\mathbf{Y}_{A1}) = \\frac{F_l}{b}(\\mathbf{X}_{B1}-\\mathbf{X}_{B2})', 'Eq. 24.9', 'Countercurrent tower; cocurrent: F_l → −F_l.'],
      ['V_r = hA_{cs} = \\frac{F_g}{\\pi K_{Ag}a}\\int_{p_{A1}}^{p_{A2}}\\frac{dp_A}{p_A-p_A^*},\\qquad \\frac{1}{K_{Ag}} = \\frac{1}{k_{Ag}}+\\frac{H_A}{k_{Al}}', 'Eq. 24.6', 'Straight mass transfer.'],
      ['F_g(\\mathbf{Y}_{A,in}-\\mathbf{Y}_{A,out}) = \\frac{F_l}{b}(\\mathbf{X}_{B,in}-\\mathbf{X}_{B,out}) = (-r_A^{\\prime\\prime\\prime})_{exit}V_r', 'Eq. 24.11', 'Mixed G / mixed L.'],
      ['V_r = F_g\\int_{\\mathbf{Y}_{A,out}}^{\\mathbf{Y}_{A,in}}\\frac{d\\mathbf{Y}_A}{(-r_A^{\\prime\\prime\\prime})_{C_{B,out}}}', 'Eq. 24.15', 'Plug G / mixed L.'],
      ['t = \\frac{f_l}{b}\\int_{C_{Bf}}^{C_{B0}}\\frac{dC_B}{-r_A^{\\prime\\prime\\prime}},\\qquad t_{min} = \\frac{V_l(C_{B0}-C_{Bf})}{bF_g\\,p_{A,in}/(\\pi-p_{A,in})}', 'Eqs. 24.19–24.20', 'Batch of liquid; the ratio is the fraction of A used.']
    ],
    method: ['Towers: the material balance fixes C_B (or the dissolved C_A for straight absorption) at every gas composition, so the volume is a single quadrature of F_g dY/(−r‴) over the gas composition, done adaptively with the rate of Eq. 23.5 evaluated self-consistently at each point (N-2301). Where the balance says B is exhausted, the remainder of the contactor is treated as straight absorption into the B-free liquid with its dissolved A followed by the balance. Switches between regimes are located exactly (for instantaneous reactions, where k_Ag p_A = k_Bl C_B/b).', 'Tanks: exit compositions from the overall balance; for plug gas through mixed liquid the gas path is integrated at the exit liquid composition. Given a size, the outlet is found by Brent’s method on the duty. Batches: at each C_B the gas leaving is found from the gas balance (mixed gas) or by integrating along the bubble path (plug gas), and the time is the quadrature of V_l dC_B/(b × rate of uptake). Ratios Y = p/(π − p) and X = C/C_U are used throughout, so rich gases are handled; for dilute systems they reduce to the book’s Eqs. 24.10, 24.12, 24.16 and 24.18.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'h';
      const o = FF.params(v);
      const pi = +v.pi, CU = +v.CT, b = o.b;
      const react = v.inst || +v.k > 0;
      const dil = !!v.dil;
      const Y = dil ? (p => p / pi) : (p => p / (pi - p)), pOfY = dil ? (y => pi * y) : (y => pi * y / (1 + y)), dYdp = dil ? (() => 1 / pi) : (p => pi / ((pi - p) * (pi - p)));
      const rate = (p, CB, CA) => FF.local(o, p, react ? Math.max(CB, 0) : 0, CA || 0);
      const ru = 'mol/(m³·' + tu + ')';
      const pat = v.pat;
      /* ---------- batch of liquid ---------- */
      if (pat === 'mb' || pat === 'pb') {
        if (!react) throw new Error('A batch of liquid is followed through the consumption of B: give a reaction.');
        const Vr = +v.Vr, Vl = o.fl * Vr, CB0 = +v.CBin, CBf = +v.CBf;
        if (!(CBf < CB0)) throw new Error('C_B at the end must be below C_B0.');
        let uptake, gasOut;
        if (v.pure) {
          uptake = CB => rate(pi, CB).r * Vr; gasOut = () => pi;
        } else {
          const Fg = +v.Fg, pin = +v.pin, Yin = Y(pin);
          if (pat === 'mb') {
            gasOut = CB => { const g = p => Fg * (Yin - Y(p)) - rate(p, CB).r * Vr; return g(0) <= 0 ? 0 : N.brent(g, 0, pin, { tol: 1e-12 * pin }); };
          } else {
            gasOut = CB => {
              const Vof = pout => Fg * N.quad(p => dYdp(p) / rate(p, CB).r, pout, pin, { rtol: 1e-7 });
              if (Vof(pin * 1e-9) <= Vr) return pin * 1e-9;
              return N.brent(pout => Vof(pout) - Vr, pin * 1e-9, pin * (1 - 1e-12), { tol: 1e-12 * pin });
            };
          }
          uptake = CB => Fg * (Yin - Y(gasOut(CB)));
        }
        const t = N.quad(CB => Vl / (b * uptake(CB)), CBf, CB0, { rtol: 1e-7 });
        kpi.push({ k: 't', l: 'Time to bring C_B from ' + fmt(CB0, 4) + ' to ' + fmt(CBf, 4), v: t, u: tu, main: true });
        const r0 = uptake(CB0) / Vr, rf = uptake(CBf) / Vr;
        kpi.push({ k: 'r0', l: '−r‴_A at the start', v: r0, u: ru }, { k: 'rf', l: '−r‴_A at the end', v: rf, u: ru });
        if (!v.pure) {
          const tmin = Vl * (CB0 - CBf) / (b * +v.Fg * Y(+v.pin));
          kpi.push({ k: 'tmin', l: 'Least time, if every A entering reacted (Eq. 24.20)', v: tmin, u: tu });
          kpi.push({ k: 'eff', l: 'Fraction of the entering A that reacts (Eq. 24.21)', v: tmin / t, s: fmt((1 - tmin / t) * 100, 3) + '% leaves unreacted' });
          kpi.push({ k: 'pout0', l: 'p_A leaving at the start → end', v: gasOut(CB0), u: 'Pa', s: '→ ' + fmt(gasOut(CBf), 4) + ' Pa' });
        } else if (+v.FA0 > 0) kpi.push({ k: 'frac0', l: 'Fraction of the A fed absorbed, start → end', v: Math.min(1, uptake(CB0) / +v.FA0), s: '→ ' + fmt(Math.min(1, uptake(CBf) / +v.FA0), 4) });
        const R0 = rate(v.pure ? pi : gasOut(CB0), CB0), Rf = rate(v.pure ? pi : gasOut(CBf), CBf);
        tables.push({ title: 'Start and end of the batch', cols: ['', 'C_B', 'M_H', 'E_i', 'E', 'regime', '−r‴_A'], rows: [['start', CB0, R0.MH, R0.Ei, R0.E, R0.case, R0.r], ['end', CBf, Rf.MH, Rf.Ei, Rf.E, Rf.case, Rf.r]] });
        const Cg = N.linspace(CBf, CB0, 80), inv = Cg.map(c => Vr / uptake(c));
        plots.push({ title: 'Processing time (Fig. 24.6)', x: { label: 'C_B', unit: 'mol/m³', min: 0 }, y: { label: '1/(−r‴_A)', unit: 'm³·' + tu + '/mol', min: 0 }, series: [{ name: '1/(−r‴_A)', x: Cg, y: inv, color: 1, w: 2.2, type: 'area', line: true, fill: 0.16 }], cap: 'Area × f_l/b is the batch time.' });
        /* C_B against time */
        const ts = [0], cs = [CB0];
        const Cq = N.linspace(CB0, CBf, 60);
        for (let i = 1; i < Cq.length; i++) { ts.push(ts[i - 1] + N.quad(c => Vl / (b * uptake(c)), Cq[i], Cq[i - 1], { rtol: 1e-6 })); cs.push(Cq[i]); }
        plots.push({ title: 'B in the batch against time', x: { label: 't', unit: tu, min: 0 }, y: { label: 'C_B', unit: 'mol/m³', min: 0 }, series: [{ name: 'C_B', x: ts, y: cs, color: 1, w: 2.2 }] });
        return { kpi, warn, info, plots, tables, diagram: D.svg(300, 150, D.tank(110, 30, 80, 90, { fill: 0.5, tag: 'R-2401', label: pat === 'mb' ? 'mixed gas' : 'gas in plug flow' }) + D.pipe([[150, 140], [150, 122]], { label: 'gas in', dx: 8, anchor: 'start' }) + D.pipe([[150, 28], [150, 8]], { label: 'gas out', dx: 8, anchor: 'start' }), 'batch of liquid') };
      }
      /* ---------- steady flows ---------- */
      const Fl = +v.Fl, CBin = +v.CBin, CAin = +v.CAin || 0;
      const XBin = CBin / CU, XAin = CAin / CU;
      if (v.pure) {
        /* p_A = π everywhere; follow the liquid */
        const FA0 = +v.FA0;
        const rAt = CB => rate(pi, CB, 0).r;
        let Vr, CBout;
        if (v.hold || !react) {
          if (!react) throw new Error('With pure gas and no reaction, give the liquid flow and use a diluted-gas description.');
          const r = rAt(CBin);
          if (v.find === 'out') { Vr = +v.Vr; CBout = CBin; }
          else throw new Error('With the liquid composition held constant, give the contactor volume.');
          const abs = r * Vr;
          kpi.push({ k: 'abs', l: 'A absorbed', v: abs, u: 'mol/' + tu, main: true, s: '−r‴_A = ' + fmt(r, 4) + ' ' + ru });
          if (FA0 > 0) kpi.push({ k: 'frac', l: 'Fraction of the entering A absorbed', v: Math.min(1, abs / FA0), s: abs > FA0 ? 'the gas would be used up before the outlet' : '' });
          const R = rate(pi, CBin, 0);
          kpi.push({ k: 'MH', l: 'M_H', v: R.MH }, { k: 'E', l: 'E', v: R.E, s: 'E_i = ' + fmt(R.Ei, 4) + ', regime ' + R.case });
          kpi.push({ k: 'shB', l: 'Resistance in the main body of liquid', v: R.share[2] });
          return { kpi, warn, info, plots, tables };
        }
        const Vfor = CBo => {
          if (pat === 'mm' || pat === 'pm') return Fl / CU / b * (CBin - CBo) / rAt(CBo);
          return Fl / (b * CU) * N.quad(c => 1 / rAt(c), CBo, CBin, { rtol: 1e-8 });
        };
        if (v.find === 'out') { Vr = +v.Vr; CBout = N.brent(c => Vfor(c) - Vr, CBin * 1e-9, CBin * (1 - 1e-12), { tol: 1e-12 * CBin }); }
        else throw new Error('With pure gas, give the contactor volume and find the outlet.');
        kpi.push({ k: 'CBout', l: 'C_B leaving', v: CBout, u: 'mol/m³', main: true });
        const abs = Fl / CU * (CBin - CBout) / b;
        kpi.push({ k: 'abs', l: 'A absorbed', v: abs, u: 'mol/' + tu });
        if (FA0 > 0) kpi.push({ k: 'frac', l: 'Fraction of the entering A absorbed', v: Math.min(1, abs / FA0) });
        return { kpi, warn, info, plots, tables };
      }
      const Fg = +v.Fg, pin = +v.pin, Yin = Y(pin);
      if (!(pin > 0 && pin < pi)) throw new Error('p_A entering must lie between 0 and π.');
      /* liquid state from the balance at a gas composition p (towers): returns {CB, CA} */
      const cc = pat === 'cc';
      const stateAt = (p, pTop) => {
        if (v.hold) return { CB: CBin, CA: CAin };
        /* countercurrent: liquid enters at the top (gas leaves at pTop); cocurrent: both enter at the bottom */
        const dY = cc ? Y(p) - Y(pTop) : Yin - Y(p);
        if (!react) return { CB: 0, CA: CU * (XAin + Fg * dY / Fl) };
        const XB = XBin - b * Fg * dY / Fl;
        if (XB > 0) return { CB: CU * XB, CA: 0 };
        /* B exhausted: beyond this point straight absorption with dissolved A */
        const dYx = Fl * XBin / (b * Fg);
        return { CB: 0, CA: CU * (Fg * (dY - dYx) / Fl), exhausted: true };
      };
      const towerV = pout => {
        let pts = [pout, pin];
        if (!v.hold && react) {
          const dYx = Fl * XBin / (b * Fg);
          const px = cc ? pOfY(Y(pout) + dYx) : pOfY(Yin - dYx);
          if (px > pout && px < pin) pts = [pout, px, pin];
        }
        let Vt = 0;
        for (let i = 1; i < pts.length; i++) Vt += Fg * N.quad(p => { const s = stateAt(p, pout), r = rate(p, s.CB, s.CA).r; return r > 0 ? dYdp(p) / r : 1e300; }, pts[i - 1], pts[i], { rtol: 1e-8 });
        return Vt;
      };
      const tankState = pout => {
        const dY = Yin - Y(pout);
        if (v.hold) return { CB: CBin, CA: CAin };
        if (!react) return { CB: 0, CA: CU * (XAin + Fg * dY / Fl) };
        const XB = XBin - b * Fg * dY / Fl;
        return XB > 0 ? { CB: CU * XB, CA: 0 } : { CB: 0, CA: CU * (Fg * (dY - Fl * XBin / (b * Fg)) / Fl), exhausted: true };
      };
      const tankV = pout => {
        const s = tankState(pout);
        if (pat === 'mm') { const r = rate(pout, s.CB, s.CA).r; return r > 0 ? Fg * (Yin - Y(pout)) / r : Infinity; }
        return Fg * N.quad(p => { const r = rate(p, s.CB, s.CA).r; return r > 0 ? dYdp(p) / r : 1e300; }, pout, pin, { rtol: 1e-8 });
      };
      const Vof = (pat === 'cc' || pat === 'co') ? towerV : tankV;
      let pout, Vr;
      if (v.find === 'out') {
        Vr = +v.Vr;
        const lo = pin * 1e-9;
        const glo = Vof(lo);
        if (!(glo > Vr)) pout = lo;
        else pout = N.brent(p => { const q = Vof(p); return (isFinite(q) && q < 1e290 ? q : 1e290) - Vr; }, lo, pin * (1 - 1e-12), { tol: 1e-12 * pin });
      } else {
        pout = +v.pout;
        if (!(pout < pin && pout >= 0)) throw new Error('The target p_A must lie below the entering value.');
        Vr = Vof(pout);
      }
      if (!(isFinite(Vr) && Vr < 1e250)) throw new Error('This duty cannot be met: the driving force vanishes before p_A falls to the target (equilibrium or B exhausted).');
      const tower = pat === 'cc' || pat === 'co';
      if (v.find === 'out') kpi.push({ k: 'pout', l: 'p_A leaving', v: pout, u: 'Pa', main: true, s: 'removal ' + fmt((1 - Y(pout) / Yin) * 100, 4) + '%' });
      else kpi.push({ k: 'V', l: 'Contactor volume V_r', v: Vr, u: 'm³', main: !tower });
      if (tower) kpi.push({ k: 'h', l: 'Height of tower', v: Vr / +v.Acs, u: 'm', main: v.find !== 'out', s: 'cross-section ' + fmt(+v.Acs, 3) + ' m²' });
      /* outlet liquid */
      const sOut = tower ? (cc ? stateAt(pin, pout) : stateAt(pout, pout)) : tankState(pout);
      if (react) kpi.push({ k: 'CBout', l: 'C_B in the liquid leaving', v: sOut.CB, u: 'mol/m³', s: sOut.exhausted ? 'B is used up inside the contactor' : '' });
      if (!react || sOut.exhausted) kpi.push({ k: 'CAout', l: 'Dissolved A in the liquid leaving', v: sOut.CA, u: 'mol/m³' });
      /* profiles */
      if (tower) {
        const ps = N.logspace(Math.max(pout, pin * 1e-6), pin, 120);
        const cum = [0];
        for (let i = 1; i < ps.length; i++) cum.push(cum[i - 1] + Fg * N.quad(p => { const s = stateAt(p, pout), r = rate(p, s.CB, s.CA).r; return r > 0 ? dYdp(p) / r : 0; }, ps[i - 1], ps[i], { rtol: 1e-6 }));
        const hs = cum.map(c => c / +v.Acs), Htot = hs[hs.length - 1];
        /* height measured from the gas inlet (bottom) for both towers */
        const z = hs.map(hh => Htot - hh);
        const st = ps.map(p => stateAt(p, pout)), Rs = ps.map((p, i) => rate(p, st[i].CB, st[i].CA));
        const ser = [{ name: 'p_A / p_A,in', x: z, y: ps.map(p => p / pin), color: 1, w: 2.2 }];
        if (react) ser.push({ name: 'C_B / max C_B', x: z, y: st.map(s => s.CB / Math.max(...st.map(q => q.CB), 1e-300)), color: 3, w: 2.2 });
        if (!react || st.some(s => s.CA > 0)) ser.push({ name: 'dissolved C_A / max', x: z, y: st.map(s => s.CA / Math.max(...st.map(q => q.CA), 1e-300)), color: 4, dash: '5 4' });
        const marks = [];
        let lastCase = Rs[0].case;
        for (let i = 1; i < Rs.length; i++) if (Rs[i].case !== lastCase) { marks.push({ t: 'v', x: z[i], color: 'var(--muted)', label: lastCase + ' | ' + Rs[i].case }); lastCase = Rs[i].case; }
        plots.push({ title: 'Along the tower, from the gas inlet', x: { label: 'height above the gas inlet', unit: 'm', min: 0, max: Htot }, y: { label: 'relative concentration', min: 0, max: 1.05 }, series: ser, marks, cap: 'Regime changes along the tower are marked with their case letters (Fig. 23.3).' });
        /* switch point for instantaneous reactions */
        if (v.inst) {
          const sw = [];
          for (let i = 1; i < Rs.length; i++) if ((Rs[i].case === 'A') !== (Rs[i - 1].case === 'A')) {
            const f = p => { const s = stateAt(p, pout), R = rate(p, s.CB, s.CA); return o.kAga * p - R.kBla * s.CB / b; };
            sw.push(N.brent(f, ps[i - 1], ps[i], { tol: 1e-10 * pin }));
          }
          if (sw.length) kpi.push({ k: 'pswitch', l: 'Reaction plane reaches the interface at p_A', v: sw[0], u: 'Pa', s: 'Case B above (lower p_A), Case A below' });
        }
        const inv = ps.map((p, i) => 1 / Rs[i].r);
        plots.push({ title: 'Performance chart: 1/(−r‴_A) against p_A (Fig. 24.3)', x: { label: 'p_A', unit: 'Pa', min: 0 }, y: { label: '1/(−r‴_A)', unit: 'm³·' + tu + '/mol', min: 0 }, series: [{ name: '1/(−r‴_A)', x: ps, y: inv, color: 1, w: 2.2, type: 'area', line: true, fill: 0.15 }], cap: 'For dilute gas the area times F_g/π is the volume.' });
        const sh = ['gas film', 'liquid film', 'main body'];
        plots.push({ title: 'Where the resistance lies along the tower', x: { label: 'height above the gas inlet', unit: 'm', min: 0, max: Htot }, y: { label: 'share of resistance', min: 0, max: 1 }, series: sh.map((nm, j) => ({ name: nm, x: z, y: Rs.map(R => R.share[j]), color: j + 1 })) });
        const Rm = Rs[Math.floor(Rs.length / 2)];
        kpi.push({ k: 'shL', l: 'Liquid-film share of resistance (mid-tower)', v: Rm.share[1], s: 'gas film ' + fmt(Rm.share[0], 3) });
        const liqOut = react && !sOut.exhausted ? 'C_B ' + fmt(sOut.CB, 4) : 'C_A ' + fmt(sOut.CA, 4);
        const liqIn = react ? 'C_B ' + fmt(CBin, 4) : 'C_A ' + fmt(CAin, 4);
        const gx0 = 20, cx = 170, cw = 64, gx1 = 400;
        let body = D.column(cx, 18, cw, 200, { tag: 'T-2401', prof: N.linspace(0, 1, 12).map(q => (cc ? 1 - q : q)) });
        body += D.pipe([[gx0, 200], [cx - 2, 200]], { label: 'gas in · ' + fmt(pin, 4) + ' Pa', dy: -7, dx: 2 });
        body += D.pipe([[cx + cw + 2, 34], [gx1, 34]], { label: 'gas out · ' + fmt(pout, 4) + ' Pa', dy: -7, dx: 8 });
        if (cc) {
          body += D.pipe([[gx0, 52], [cx - 2, 52]], { label: 'liquid in · ' + liqIn, dy: 16, dx: 2 });
          body += D.pipe([[cx + cw + 2, 186], [gx1, 186]], { label: 'liquid out · ' + liqOut, dy: 16, dx: 8 });
        } else {
          body += D.pipe([[gx0, 182], [cx - 2, 182]], { label: 'liquid in · ' + liqIn, dy: 16, dx: 2 });
          body += D.pipe([[cx + cw + 2, 52], [gx1, 52]], { label: 'liquid out · ' + liqOut, dy: 16, dx: 8 });
        }
        body += D.text(cx + cw / 2, 252, 'h = ' + fmt(Vr / +v.Acs, 4) + ' m');
        const diag = D.svg(420, 262, body, cc ? 'countercurrent absorption tower' : 'cocurrent absorption tower', 1.15);
        return { kpi, warn, info, plots, tables, diagram: diag };
      }
      const R = rate(pout, sOut.CB, sOut.CA);
      kpi.push({ k: 'r', l: '−r‴_A at the exit composition', v: R.r, u: ru, s: 'regime ' + R.case + ', M_H = ' + fmt(R.MH, 4) + ', E = ' + fmt(R.E, 4) });
      kpi.push({ k: 'shL', l: 'Liquid-film share of resistance', v: R.share[1], s: 'gas film ' + fmt(R.share[0], 3) + ', main body ' + fmt(R.share[2], 3) });
      if (pat === 'pm') {
        const ps = N.linspace(pout, pin, 100);
        plots.push({ title: '1/(−r‴_A) along the rising bubbles (Eq. 24.16)', x: { label: 'p_A', unit: 'Pa' }, y: { label: '1/(−r‴_A)', min: 0 }, series: [{ name: 'at C_B,out', x: ps, y: ps.map(p => 1 / rate(p, sOut.CB, sOut.CA).r), color: 1, type: 'area', line: true, fill: 0.15 }] });
      }
      return { kpi, warn, info, plots, tables, diagram: D.svg(300, 150, D.tank(110, 30, 80, 90, { fill: 0.5, tag: 'R-2401', label: pat === 'mm' ? 'both mixed' : 'gas in plug flow' }), 'gas–liquid tank') };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
