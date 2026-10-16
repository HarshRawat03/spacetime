/* Chapter 1 — Overview; Chapter 2 — Kinetics of homogeneous reactions */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, KIT = ST.KIT, D = ST.D;
  const RG = 8.314;

  /* ---------------- Ch 1: rate of reaction from plant data ---------------- */
  const LADDER = [
    ['Cellular reactions, water treatment plants', 1e-4, 1.2e-3],
    ['Human at rest', 1.5e-3, 2.5e-3],
    ['Human working hard', 5e-3, 9e-3],
    ['Gases in porous catalyst particles', 3e-2, 5],
    ['Coal furnaces', 2, 25],
    ['Jet engines', 5e2, 3e3],
    ['Rocket engines', 6e3, 5e4],
    ['Bimolecular reactions, every collision counts (~1 atm, 400 °C)', 1e9, 1e10]
  ];
  const DAY = 86400;
  ST.tool({
    id: 'ratedef', tag: 'O-101', ch: 1, title: 'Rate of reaction from process data',
    lede: 'Express how fast something reacts per unit volume of fluid, mass of solid, interfacial surface, volume of solid or volume of reactor, and see where it sits on the book’s ladder of reaction rates.',
    keys: 'rate definition intensive rate per volume per mass per surface rocket engine living person figure 1.3',
    sweep: false,
    inputs: [
      { group: 'How much reacts' },
      { k: 'how', t: 'select', label: 'Amount of component i reacting', opts: [['mol', 'given directly, mol/s'], ['mass', 'from a mass flow of some product'], ['energy', 'from an energy release']], v: 'mass' },
      { k: 'dN', label: 'Moles of i reacting', unit: 'mol/s', v: 1, min: 1e-6, max: 1e9, scale: 'log', if: v => v.how === 'mol' },
      { k: 'mdot', label: 'Mass flow of measured species', unit: 'kg/s', v: 108, min: 1e-9, max: 1e6, scale: 'log', if: v => v.how === 'mass' },
      { k: 'mw', label: 'Its molar mass', unit: 'g/mol', v: 18, min: 1, max: 1000, scale: 'log', if: v => v.how === 'mass' },
      { k: 'E', label: 'Energy released', unit: 'kJ/day', v: 6000, min: 1e-3, max: 1e12, scale: 'log', if: v => v.how === 'energy' },
      { k: 'dH', label: 'Heat released per mole of fuel', unit: 'kJ/mol', v: 2816, min: 1e-3, max: 1e5, scale: 'log', if: v => v.how === 'energy' },
      { k: 'nu', label: 'Moles of i per mole measured', v: 1, min: 0.01, max: 100, scale: 'log', hint: 'Rocket: 1 mol H₂ per mol H₂O. Person: 6 mol O₂ per mol glucose.' },
      { group: 'Size of the thing' },
      { k: 'V', label: 'Fluid volume V', unit: 'm³', v: 0.2121, min: 1e-9, max: 1e9, scale: 'log' },
      { k: 'W', label: 'Mass of solid W (optional)', unit: 'kg', v: 0, min: 0, max: 1e9, slider: false },
      { k: 'S', label: 'Interfacial surface S (optional)', unit: 'm²', v: 0, min: 0, max: 1e12, slider: false },
      { k: 'Vs', label: 'Volume of solid (optional)', unit: 'm³', v: 0, min: 0, max: 1e9, slider: false },
      { k: 'Vr', label: 'Reactor volume (optional)', unit: 'm³', v: 0, min: 0, max: 1e9, slider: false }
    ],
    presets: [
      { name: 'Example 1.1', v: { how: 'mass', mdot: 108, mw: 18, nu: 1, V: Math.PI * 0.3 * 0.3 * 0.75 }, book: { r: 2.829e4 }, note: 'Rocket engine: 108 kg/s of steam from a chamber 75 cm long and 60 cm across. Use ν = 0.5 for oxygen (book: 1.415 × 10⁴).' },
      { name: 'Example 1.2', v: { how: 'energy', E: 6000, dH: 2816, nu: 6, V: 0.075 }, book: { r: { v: 0.002, tol: 0.06 } }, note: 'A 75 kg person eating 6000 kJ a day of glucose, density about 1000 kg/m³.' },
      { name: 'Problem 1.1', v: { how: 'mass', mdot: 32000 * 0.2 / DAY, mw: 32, nu: 1, V: 32000 * 8 / 24 }, note: '32 000 m³/day of waste water with a BOD of 200 mg O₂/L held 8 h.' }
    ],
    eqs: [
      ['r_i = \\frac{1}{V}\\frac{dN_i}{dt},\\quad r_i\' = \\frac{1}{W}\\frac{dN_i}{dt},\\quad r_i\'\' = \\frac{1}{S}\\frac{dN_i}{dt},\\quad r_i\'\'\' = \\frac{1}{V_s}\\frac{dN_i}{dt},\\quad r_i\'\'\'\' = \\frac{1}{V_r}\\frac{dN_i}{dt}', 'Eqs. 1.2–1.6'],
      ['V r_i = W r_i\' = S r_i\'\' = V_s r_i\'\'\' = V_r r_i\'\'\'\'', 'Eq. 1.7']
    ],
    run(v) {
      let dN;
      if (v.how === 'mol') dN = +v.dN * +v.nu;
      else if (v.how === 'mass') dN = (+v.mdot * 1000 / +v.mw) * +v.nu;
      else dN = (+v.E / +v.dH) * +v.nu / DAY;
      const r = dN / +v.V;
      const kpi = [{ k: 'r', l: '−r_i per volume of fluid', v: r, u: 'mol/(m³·s)', main: true }, { k: 'dN', l: 'Moles reacting', v: dN, u: 'mol/s' }];
      const rows = [['fluid volume V', +v.V, 'm³', r, 'mol/(m³ fluid·s)']];
      if (+v.W > 0) { kpi.push({ k: 'rW', l: '−r′ per mass of solid', v: dN / v.W, u: 'mol/(kg·s)' }); rows.push(['mass of solid W', +v.W, 'kg', dN / v.W, 'mol/(kg·s)']); }
      if (+v.S > 0) { kpi.push({ k: 'rS', l: '−r″ per surface', v: dN / v.S, u: 'mol/(m²·s)' }); rows.push(['surface S', +v.S, 'm²', dN / v.S, 'mol/(m²·s)']); }
      if (+v.Vs > 0) { kpi.push({ k: 'rVs', l: '−r‴ per solid volume', v: dN / v.Vs, u: 'mol/(m³ solid·s)' }); rows.push(['solid volume', +v.Vs, 'm³', dN / v.Vs, 'mol/(m³ solid·s)']); }
      if (+v.Vr > 0) { kpi.push({ k: 'rVr', l: '−r⁗ per reactor volume', v: dN / v.Vr, u: 'mol/(m³ reactor·s)' }); rows.push(['reactor volume', +v.Vr, 'm³', dN / v.Vr, 'mol/(m³ reactor·s)']); }
      const yy = LADDER.map((_, i) => LADDER.length - i);
      const series = LADDER.map((L, i) => ({ name: L[0], x: [L[1], L[2]], y: [yy[i], yy[i]], color: 8, w: 7, legend: false }));
      const marks = LADDER.map((L, i) => ({ t: 'text', x: Math.sqrt(L[1] * L[2]), y: yy[i] + 0.38, s: L[0].replace(/,.*/, ''), anchor: 'middle', color: 8, noRange: true }));
      marks.push({ t: 'v', x: r, label: 'yours: ' + N.fmt(r, 3), color: 2, dash: false });
      const info = [];
      const near = LADDER.filter(L => r >= L[1] / 3 && r <= L[2] * 3).map(L => L[0].toLowerCase());
      if (near.length) info.push('Comparable to ' + near.join('; ') + '.');
      return {
        kpi, info,
        plots: [{ title: 'Where this rate sits (Fig. 1.3, bands read from the figure)', wide: true, h: 330, legend: false, x: { label: 'moles of A disappearing per m³ of thing per second', log: true, min: 3e-5, max: 3e10 }, y: { label: '', min: 0.3, max: LADDER.length + 1, tight: true, ticks: false }, series, marks }],
        tables: [{ title: 'Rate on each basis', cols: ['basis', 'size', 'unit', 'rate', 'rate unit'], rows }]
      };
    }
  });

  /* ---------------- Ch 2: rate constants, units, stoichiometry ---------------- */
  const CU = { 'mol/L': 1, 'mol/m³': 1e-3, 'kmol/m³': 1, 'mol/cm³': 1000, 'mmol/L': 1e-3, 'lbmol/ft³': 16.0185 };
  const TU = { s: 1, min: 60, h: 3600, day: 86400 };
  const PU = { atm: 1, Pa: 1 / 101325, kPa: 1 / 101.325, bar: 1 / 1.01325, mmHg: 1 / 760, psi: 1 / 14.6959 };
  const opts = o => Object.keys(o).map(k => [k, k]);
  ST.tool({
    id: 'kunits', tag: 'H-201', ch: 2, sec: '2.1', title: 'Rate constants: units, pressure basis, stoichiometry',
    lede: 'Convert an nth-order rate constant between concentration and time units, between pressure and concentration bases for gases, and turn the rate of one component into the rates of all the others.',
    keys: 'rate constant units conversion kp kc pressure units order stoichiometric rates elementary',
    sweep: false,
    inputs: [
      { group: 'Rate constant' },
      { k: 'n', label: 'Overall order n', v: 2, min: -1, max: 4, step: 0.5 },
      { k: 'k', label: 'k as given', v: 0.005, min: 1e-12, max: 1e12, slider: false },
      { k: 'basis', t: 'select', label: 'Given on which basis', opts: [['C', '−r_A = k C_A^n (concentrations)'], ['p', '−r_A = k p_A^n (rate per volume, pressures)'], ['dp', '−dp_A/dt = k p_A^n (pressure rate)']], v: 'C' },
      { k: 'cu1', t: 'select', label: 'Concentration unit given', opts: opts(CU), v: 'mol/cm³', if: v => v.basis === 'C' || v.basis === 'p' },
      { k: 'pu1', t: 'select', label: 'Pressure unit given', opts: opts(PU), v: 'atm', if: v => v.basis !== 'C' },
      { k: 'tu1', t: 'select', label: 'Time unit given', opts: opts(TU), v: 'min' },
      { k: 'T', label: 'Temperature', unit: 'K', v: 400, min: 1, max: 3000, if: v => v.basis !== 'C' },
      { group: 'Wanted as −r_A = k C_A^n in' },
      { k: 'cu2', t: 'select', label: 'Concentration unit', opts: opts(CU), v: 'mol/L' },
      { k: 'tu2', t: 'select', label: 'Time unit', opts: opts(TU), v: 'h' },
      { group: 'Rates of the other components' },
      { k: 'stoich', t: 'text', label: 'Stoichiometric equation', v: '2NO2 + 0.5O2 -> N2O5', mono: true },
      { k: 'key', t: 'text', label: 'Component with a known rate', v: 'NO2', mono: true },
      { k: 'rkey', label: 'Its rate of disappearance', v: 1, min: -1e9, max: 1e9, slider: false }
    ],
    presets: [
      { name: 'Problem 2.6', v: { n: 2, k: 0.005, basis: 'C', cu1: 'mol/cm³', tu1: 'min', cu2: 'mol/L', tu2: 'h' }, note: '−r_A = 0.005 C_A² in mol/cm³ and minutes, wanted in mol/L and hours.' },
      { name: 'Problem 2.7', v: { n: 2, k: 3.66, basis: 'dp', pu1: 'atm', tu1: 'h', T: 400, cu2: 'mol/L', tu2: 'h' }, note: '−dp_A/dt = 3.66 p_A² atm/h at 400 K.' },
      { name: 'Example 3.4, 400 K', v: { n: 2, k: 2.3, basis: 'p', cu1: 'mol/m³', pu1: 'atm', tu1: 's', T: 400, cu2: 'mol/m³', tu2: 's' }, note: '−r_A = 2.3 p_A² mol/(m³·s·atm²) at 400 K. Repeat at 500 K to see why E differs with the basis.' },
      { name: 'Problem 2.2', v: { stoich: '2NO2 + 0.5O2 -> N2O5', key: 'NO2', rkey: 1 } }
    ],
    eqs: [
      ['k_2 = k_1\\left(\\frac{c_1}{c_2}\\right)^{1-n}\\frac{t_2}{t_1}', '', 'Units of k are (concentration)^{1−n}·time⁻¹.'],
      ['-r_A = k_pp_A^n = k_p(RT)^nC_A^n \\;\\Rightarrow\\; k_C = k_p(RT)^n', 'Eq. 2.6'],
      ['-\\frac{dp_A}{dt} = k\'_pp_A^n \\;\\Rightarrow\\; k_C = k\'_p(RT)^{n-1}', '', 'Constant volume ideal gas: $-r_A = -(1/RT)\\,dp_A/dt$.'],
      ['\\frac{-r_A}{a} = \\frac{-r_B}{b} = \\frac{r_R}{r} = \\frac{r_S}{s}', 'Eq. 2.1']
    ],
    run(v) {
      const n = +v.n, warn = [];
      let kC, cuIn;
      const Rlit = 0.082057;
      if (v.basis === 'C') { kC = +v.k; cuIn = CU[v.cu1]; }
      else {
        const patm = PU[v.pu1];
        const RT = Rlit * (+v.T);
        /* PU = atm per pressure unit; R in L·atm/(mol·K) so concentrations come out in mol/L */
        if (v.basis === 'p') kC = (+v.k) * CU[v.cu1] * Math.pow(RT / patm, n);
        else kC = (+v.k) * Math.pow(RT / patm, n - 1);
        cuIn = 1;
      }
      const k2 = kC * Math.pow(cuIn / CU[v.cu2], 1 - n) * (TU[v.tu2] / TU[v.tu1]);
      const unit = KIT.kUnit(n, v.cu2, v.tu2);
      const kpi = [{ k: 'k2', l: 'k', v: k2, u: unit, main: true, s: '−r_A = k C_A^n' }];
      if (v.basis !== 'C') kpi.push({ k: 'kCL', l: 'k in mol/L and ' + v.tu1, v: kC, u: KIT.kUnit(n, 'mol/L', v.tu1) });
      /* stoichiometric rates */
      let rows = [];
      try {
        const net = ST.parseNet(String(v.stoich).replace(/->|→|=/, '->') + ' ; k = 1');
        const r0 = net.rx[0];
        const nu = {};
        r0.reac.forEach(([i, c]) => { nu[net.species[i]] = -c; });
        r0.prod.forEach(([i, c]) => { nu[net.species[i]] = (nu[net.species[i]] || 0) + c; });
        const key = String(v.key).trim();
        if (!(key in nu)) throw new Error('“' + key + '” is not in the equation');
        const rate = -(+v.rkey) / nu[key];
        rows = Object.keys(nu).map(s => [s, nu[s], nu[s] * rate, nu[s] < 0 ? 'disappearing' : 'forming']);
        kpi.push({ k: 'rx', l: 'Rate of reaction (per unit of stoichiometric equation)', v: rate, s: 'r_i = ν_i × this' });
      } catch (e) { warn.push('Stoichiometry: ' + e.message); }
      return { kpi, warn, tables: [{ title: 'Rate of formation of each component (negative = disappearing)', cols: ['component', 'ν_i', 'r_i', ''], rows }] };
    }
  });

  /* ---------------- Ch 2/3: Arrhenius ---------------- */
  ST.tool({
    id: 'arrhenius', tag: 'H-202', ch: 2, sec: '2.2, 3.3', title: 'Arrhenius analysis and temperature sensitivity',
    lede: 'Activation energy from two points, from times at two temperatures, or by fitting rate data; rates at other temperatures; the temperature rise that doubles a rate; and curved Arrhenius plots from two competing steps.',
    keys: 'activation energy arrhenius frequency factor ln k 1/T temperature dependence pasteurization doubling table 2.1 collision theory transition state',
    inputs: [
      { group: 'Question' },
      { k: 'mode', t: 'select', label: 'What do you have', opts: [['fit', 'k (or a rate) at several temperatures'], ['two', 'k at two temperatures'], ['time', 'time for the same result at two temperatures'], ['pred', 'E and k at one temperature'], ['mech', 'two competing or consecutive steps']], v: 'fit' },
      { k: 'tunit', t: 'seg', label: 'Temperatures in', opts: [['K', 'K'], ['C', '°C']], v: 'K' },
      { k: 'data', t: 'table', label: 'Data', cols: ['T', 'k or rate'], v: [[781.15, 0.1059], [700.15, 0.0031], [666.15, 0.000588], [629.15, 80.9e-6], [556.15, 0.942e-6]], if: v => v.mode === 'fit' },
      { k: 'inv', t: 'check', label: 'Second column is a time (rate ∝ 1/time)', v: false, if: v => v.mode === 'fit' },
      { k: 'm', label: 'Exponent m in k = k₀Tᵐe^{−E/RT}', v: 0, min: -1, max: 3, step: 0.5, if: v => v.mode === 'fit', hint: '0 is Arrhenius; ½ collision theory; 1 transition state theory.' },
      { k: 'pconv', t: 'check', label: 'k was measured with partial pressures (convert to concentrations)', v: false, if: v => v.mode === 'fit' },
      { k: 'n', label: 'Order n for the conversion', v: 2, min: 0, max: 4, step: 0.5, if: v => v.mode === 'fit' && v.pconv },
      { k: 'T1', label: 'T₁', v: 336, min: 1, max: 3000, if: v => v.mode === 'two' || v.mode === 'time' || v.mode === 'pred' },
      { k: 'k1', label: 'k₁', v: 1, min: 1e-30, max: 1e30, slider: false, if: v => v.mode === 'two' || v.mode === 'pred' },
      { k: 't1', label: 'Time t₁', v: 1800, min: 1e-9, max: 1e12, slider: false, if: v => v.mode === 'time' },
      { k: 'T2', label: 'T₂', v: 347, min: 1, max: 3000, if: v => v.mode !== 'fit' && v.mode !== 'mech' },
      { k: 'k2', label: 'k₂', v: 2, min: 1e-30, max: 1e30, slider: false, if: v => v.mode === 'two' },
      { k: 't2', label: 'Time t₂', v: 15, min: 1e-9, max: 1e12, slider: false, if: v => v.mode === 'time' },
      { k: 'E', label: 'Activation energy E', unit: 'kJ/mol', v: 175, min: 1, max: 1000, if: v => v.mode === 'pred' },
      { group: 'Two steps', if: v => v.mode === 'mech' },
      { k: 'arr', t: 'seg', label: 'Arrangement', opts: [['par', 'parallel paths'], ['ser', 'steps in series']], v: 'par', if: v => v.mode === 'mech' },
      { k: 'Ea', label: 'E of step 1', unit: 'kJ/mol', v: 40, min: 1, max: 600, if: v => v.mode === 'mech' },
      { k: 'Eb', label: 'E of step 2', unit: 'kJ/mol', v: 160, min: 1, max: 600, if: v => v.mode === 'mech' },
      { k: 'Tx', label: 'Temperature where both steps are equal', unit: 'K', v: 600, min: 100, max: 3000, if: v => v.mode === 'mech' },
      { k: 'Tlo', label: 'Range from', unit: 'K', v: 400, min: 50, max: 3000, if: v => v.mode === 'mech' },
      { k: 'Thi', label: 'to', unit: 'K', v: 1000, min: 50, max: 4000, if: v => v.mode === 'mech' }
    ],
    presets: [
      { name: 'Example 2.3', v: { mode: 'time', tunit: 'K', T1: 336, t1: 1800, T2: 347, t2: 15 }, book: { E: 422 }, note: 'Milk pasteurized by 30 min at 63 °C or 15 s at 74 °C.' },
      { name: 'Problem 3.31', v: { mode: 'fit', tunit: 'C', data: [[508, 0.1059], [427, 0.00310], [393, 0.000588], [356, 80.9e-6], [283, 0.942e-6]], inv: false, m: 0, pconv: false }, note: 'Bodenstein’s data for the decomposition of hydrogen iodide, k in cm³/(mol·s).' },
      { name: 'Example 3.4', v: { mode: 'fit', tunit: 'K', data: [[400, 2.3], [500, 2.3]], inv: false, m: 0, pconv: true, n: 2 }, book: { E: { v: 7.394, tol: 0.01 } }, note: 'The same rate in pressure units at 400 K and 500 K gives E = 0 on a pressure basis but about 7.4 kJ/mol on a concentration basis.' },
      { name: 'Problem 2.9', v: { mode: 'pred', tunit: 'C', E: 300, T1: 500, k1: 1, T2: 650 }, note: 'Ethane pyrolysis, E ≈ 300 kJ/mol: how much faster at 650 °C than at 500 °C?' },
      { name: 'Problem 2.10', v: { mode: 'two', tunit: 'K', T1: 1000, k1: 1, T2: 1100, k2: 20 }, note: 'n-Nonane cracks 20 times faster at 1100 K than at 1000 K.' },
      { name: 'Problem 2.11', v: { mode: 'fit', tunit: 'C', data: [[13, 150], [16, 160], [22, 230], [24, 295], [28, 370]], inv: false, m: 0, pconv: false }, note: 'Running speed of Oregon ants (m/h) against temperature.' },
      { name: 'Problem 2.12', v: { mode: 'pred', tunit: 'K', E: 175, T1: 780, k1: 1, T2: 792 }, note: 'Raising the set point from 780 K to 792 K with E = 175 kJ/mol.' },
      { name: 'Problem 2.13', v: { mode: 'fit', tunit: 'C', data: [[22.0, 87], [23.4, 85], [26.3, 74], [24.3, 78], [21.1, 90], [22.7, 84]], inv: true, m: 0, pconv: false }, note: 'Watermelon growing days against mean temperature; growth rate taken as 1/days.' },
      { name: 'Problem 2.14', v: { mode: 'fit', tunit: 'K', data: [60, 65, 70, 75, 80].map(F => [(F - 32) / 1.8 + 273.15, F - 40]), inv: false, m: 0, pconv: false }, note: 'Dolbear’s cricket law, chirps per 15 s = T(°F) − 40, between 60 and 80 °F.' },
      { name: 'Figure 2.3', v: { mode: 'mech', arr: 'par', Ea: 40, Eb: 160, Tx: 600, Tlo: 400, Thi: 1000 }, note: 'Parallel paths: the apparent E rises with temperature.' }
    ],
    eqs: [
      ['k = k_0e^{-E/RT}', 'Eq. 2.34'],
      ['\\ln\\frac{k_2}{k_1} = \\ln\\frac{r_2}{r_1} = \\frac{E}{R}\\left(\\frac{1}{T_1}-\\frac{1}{T_2}\\right)', 'Eq. 2.35'],
      ['k \\propto T^me^{-E/RT},\\qquad 0\\le m\\le 1', 'Eq. 2.36', 'Collision and transition-state theory; the exponential dominates.'],
      ['\\Delta T_{\\times 2}:\\quad \\frac{E}{R}\\left(\\frac{1}{T}-\\frac{1}{T+\\Delta T}\\right)=\\ln 2', 'Table 2.1']
    ],
    method: ['Fits regress ln k − m ln T against 1/T by ordinary least squares; the 95% interval on E uses Student’s t with n − 2 degrees of freedom.', 'With the pressure option, each k measured as $-r_A = k_pp_A^n$ is converted to $k_C = k_p(RT)^n$ before fitting, which is the correction Example 3.4 makes.'],
    run(v) {
      const toK = T => (v.tunit === 'C' ? +T + 273.15 : +T);
      const kpi = [], warn = [], plots = [], tables = [];
      const dbl = (E, T) => { const a = 1 / T - Math.log(2) * RG / E; return 1 / a - T; };
      if (v.mode === 'fit') {
        let pts = (v.data || []).filter(r => isFinite(r[0]) && isFinite(r[1]) && r[1] > 0).map(r => ({ T: toK(r[0]), k: v.inv ? 1 / r[1] : r[1] }));
        if (pts.length < 2) throw new Error('Enter at least two temperatures.');
        if (v.pconv) pts = pts.map(p => ({ T: p.T, k: p.k * Math.pow(0.082057 * p.T, +v.n), raw: p.k }));
        const m = +v.m || 0;
        const lr = N.linreg(pts.map(p => 1 / p.T), pts.map(p => Math.log(p.k) - m * Math.log(p.T)));
        const E = -lr.b * RG, k0 = Math.exp(lr.a);
        const tq = lr.dof > 0 ? N.tinv(0.975, lr.dof) : NaN;
        kpi.push({ k: 'E', l: 'Activation energy E', v: E / 1000, u: 'kJ/mol', main: true, s: isFinite(tq) && isFinite(lr.sb) ? '± ' + N.fmt(tq * lr.sb * RG / 1000, 3) + ' (95%)' : (pts.length === 2 ? 'two points, no interval' : '') });
        kpi.push({ k: 'k0', l: 'Frequency factor k₀', v: k0, s: m ? 'with Tᵐ, m = ' + m : 'same units as k' });
        kpi.push({ k: 'ER', l: 'E/R', v: E / RG, u: 'K' });
        if (pts.length > 2) kpi.push({ k: 'r2', l: 'R²', v: lr.r2 });
        const Tm = N.mean(pts.map(p => p.T));
        kpi.push({ k: 'dbl', l: 'Rise to double the rate', v: dbl(E, Tm), u: 'K', s: 'at ' + N.fmt(Tm, 4) + ' K' });
        const xs = pts.map(p => 1000 / p.T), xr = N.linspace(Math.min(...xs) * 0.98, Math.max(...xs) * 1.02, 40);
        plots.push({ title: 'Arrhenius plot', x: { label: '1000/T', unit: '1/K' }, y: { label: 'k' + (v.pconv ? '_C' : ''), log: true }, series: [{ name: 'data', type: 'pts', x: xs, y: pts.map(p => p.k), color: 2 }, { name: 'fit', x: xr, y: xr.map(x => k0 * Math.pow(1000 / x, m) * Math.exp(-E / RG * x / 1000)), color: 1 }] });
        tables.push({ title: 'Data and fit', cols: ['T, K', 'k used', 'fitted k', 'residual in ln k'].concat(v.pconv ? ['k as measured'] : []), rows: pts.map(p => { const f = k0 * Math.pow(p.T, m) * Math.exp(-E / RG / p.T); return [p.T, p.k, f, Math.log(p.k / f)].concat(v.pconv ? [p.raw] : []); }) });
        if (v.pconv) warn.push('On the pressure basis the same data give E = ' + N.fmt(-N.linreg(pts.map(p => 1 / p.T), pts.map(p => Math.log(p.raw))).b * RG / 1000, 4) + ' kJ/mol. Use the concentration-based value.');
      } else if (v.mode === 'two' || v.mode === 'time') {
        const T1 = toK(v.T1), T2 = toK(v.T2);
        const ratio = v.mode === 'two' ? (+v.k2) / (+v.k1) : (+v.t1) / (+v.t2);
        const E = RG * Math.log(ratio) / (1 / T1 - 1 / T2);
        kpi.push({ k: 'E', l: 'Activation energy E', v: E / 1000, u: 'kJ/mol', main: true }, { k: 'ratio', l: 'Rate ratio r₂/r₁', v: ratio }, { k: 'dbl', l: 'Rise to double the rate', v: dbl(E, T1), u: 'K', s: 'at T₁' });
        if (v.mode === 'two') kpi.push({ k: 'k0', l: 'k₀', v: (+v.k1) * Math.exp(E / RG / T1) });
      } else if (v.mode === 'pred') {
        const T1 = toK(v.T1), T2 = toK(v.T2), E = +v.E * 1000;
        const ratio = Math.exp(E / RG * (1 / T1 - 1 / T2));
        kpi.push({ k: 'ratio', l: 'Rate ratio r₂/r₁', v: ratio, main: true }, { k: 'k2', l: 'k at T₂', v: (+v.k1) * ratio }, { k: 'dbl', l: 'Rise to double the rate', v: dbl(E, T1), u: 'K', s: 'at T₁' });
      } else {
        const Ea = +v.Ea * 1000, Eb = +v.Eb * 1000, Tx = +v.Tx;
        const ka = T => Math.exp(-Ea / RG * (1 / T - 1 / Tx)), kb = T => Math.exp(-Eb / RG * (1 / T - 1 / Tx));
        const kt = T => (v.arr === 'par' ? ka(T) + kb(T) : 1 / (1 / ka(T) + 1 / kb(T)));
        const Ts = N.linspace(+v.Tlo, +v.Thi, 120), x = Ts.map(T => 1000 / T);
        plots.push({ title: 'ln k against 1/T', x: { label: '1000/T', unit: '1/K' }, y: { label: 'k (relative)', log: true }, series: [{ name: 'overall', x, y: Ts.map(kt), color: 1, w: 2.6 }, { name: 'step 1', x, y: Ts.map(ka), color: 2, dash: '5 4' }, { name: 'step 2', x, y: Ts.map(kb), color: 3, dash: '5 4' }] });
        const Eapp = T => -RG * (Math.log(kt(T * 1.0005)) - Math.log(kt(T / 1.0005))) / (1 / (T * 1.0005) - 1 / (T / 1.0005)) / 1000;
        plots.push({ title: 'Apparent activation energy', x: { label: 'T', unit: 'K' }, y: { label: 'E_{app}', unit: 'kJ/mol' }, series: [{ name: 'E_app', x: Ts, y: Ts.map(Eapp), color: 1 }] });
        kpi.push({ k: 'Elo', l: 'E apparent at ' + v.Tlo + ' K', v: Eapp(+v.Tlo), u: 'kJ/mol' }, { k: 'Ehi', l: 'E apparent at ' + v.Thi + ' K', v: Eapp(+v.Thi), u: 'kJ/mol' });
        warn.push(v.arr === 'par' ? 'Parallel paths: the path with the larger E takes over at high temperature, so the apparent E rises with T (Fig. 2.3).' : 'Steps in series: the slower step controls; at high T it is the one with the smaller E, so the apparent E falls as T rises (Fig. 2.3).');
      }
      /* Table 2.1 regenerated */
      if (v.mode !== 'mech') {
        const Es = [40, 160, 280, 400], Tc = [0, 400, 1000, 2000];
        tables.push({ title: 'Temperature rise that doubles the rate (Table 2.1), K', cols: ['T, °C'].concat(Es.map(e => 'E = ' + e + ' kJ/mol')), rows: Tc.map(t => [t].concat(Es.map(e => dbl(e * 1000, t + 273.15)))), d: 3 });
      }
      return { kpi, warn, plots, tables };
    }
  });

  /* ---------------- Ch 2: mechanism tester ---------------- */
  ST.tool({
    id: 'mechanism', tag: 'H-203', ch: 2, sec: '2.1, 2.3', title: 'Mechanism tester: full kinetics against the steady-state rate',
    lede: 'Write a set of elementary steps, simulate them exactly in a batch, and compare the true rate of formation of a product with the rate law you derived using the steady-state or equilibrium assumption.',
    keys: 'mechanism elementary steps intermediate steady state approximation pssh michaelis menten lindemann chain reaction free radicals hbr ozone n2o5 apparent order',
    units: { C: 'mol/L', t: 's' }, slow: true,
    inputs: [
      { group: 'Mechanism' },
      { k: 'net', t: 'area', label: 'Elementary steps (one per line)', rows: 7, v: 'B + B <=> B2s ; k = 2, kr = 50\nA + B2s -> AB + B ; k = 400', hint: 'Format: A + B -> R ; k = 1. Reversible: <=> with kr. Global parameters: name = value.' },
      { k: 'init', t: 'text', label: 'Initial concentrations', v: 'A = 1, B = 1', mono: true },
      { k: 'tEnd', label: 'Run time', unit: '[t]', v: 20, min: 1e-6, max: 1e8, scale: 'log' },
      { group: 'Your rate law' },
      { k: 'prod', t: 'text', label: 'Product followed', v: 'AB', mono: true },
      { k: 'law', t: 'text', label: 'Proposed rate of formation', v: 'k1*k3*A*B^2/(k2 + k3*A)', mono: true, hint: 'Use species names and the constants defined below.' },
      { k: 'consts', t: 'text', label: 'Constants used in it', v: 'k1 = 2, k2 = 50, k3 = 400', mono: true },
      { k: 'logy', t: 'check', label: 'Log scale for concentrations', v: true }
    ],
    presets: [
      { name: 'Example 2.1, model 2', v: { net: 'B + B <=> B2s ; k = 2, kr = 50\nA + B2s -> AB + B ; k = 400', init: 'A = 1, B = 1', tEnd: 20, prod: 'AB', law: 'k1*k3*A*B^2/(k2 + k3*A)', consts: 'k1 = 2, k2 = 50, k3 = 400' }, note: 'The intermediate B₂* stays tiny, and when k₃[A] ≫ k₂ the rate reduces to k₁[B]², the observed second order in B.' },
      { name: 'Example 2.2, Michaelis–Menten', v: { net: 'A + E <=> X ; k = 100, kr = 10\nX -> R + E ; k = 5', init: 'A = 1, E = 0.002', tEnd: 1500, prod: 'R', law: 'k3*E0*A/(CM + A)', consts: 'k3 = 5, E0 = 0.002, CM = 0.15', logy: true }, note: 'Steady state on the complex gives $C_M$ = (k₂ + k₃)/k₁ = 0.15. Raise $E_0$ towards $C_A$ and watch the approximation fail.' },
      { name: 'Lindemann (azomethane)', v: { net: 'A + A <=> As + A ; k = 1, kr = 100\nAs -> R + S ; k = 50', init: 'A = 5', tEnd: 20, prod: 'R', law: 'k1*k3*A^2/(k2*A + k3)', consts: 'k1 = 1, k2 = 100, k3 = 50', logy: true }, note: 'First order while k₂[A] ≫ k₃ (high [A]); second order once [A] falls below k₃/k₂ = 0.5. Watch the apparent order shift.' },
      { name: 'H₂ + Br₂ chain', v: { net: 'Br2 <=> Br + Br ; k = 0.01, kr = 1000\nBr + H2 -> HBr + H ; k = 5\nH + Br2 -> HBr + Br ; k = 5000\nH + HBr -> H2 + Br ; k = 500', init: 'H2 = 1, Br2 = 1', tEnd: 200, prod: 'HBr', law: '2*k2*sqrt(k1/k5)*H2*sqrt(Br2)/(1 + (k4/k3)*HBr/Br2)', consts: 'k1 = 0.01, k5 = 1000, k2 = 5, k3 = 5000, k4 = 500', logy: true }, note: 'Illustrative constants. The steady-state result has the form of the observed rate quoted in §2.1.' },
      { name: 'Problem 2.20, ozone', v: { net: 'O3 <=> O2 + O ; k = 1, kr = 1000\nO + O3 -> O2 + O2 ; k = 10', init: 'O3 = 1, O2 = 0.5', tEnd: 300, prod: 'O2', law: '3*k1*k3*O3^2/(k2*O2)', consts: 'k1 = 1, k2 = 1000, k3 = 10', logy: true }, note: 'Two-step scheme giving −r(O₃) = k[O₃]²/[O₂] (illustrative constants; the factor 3 counts O₂ from both steps).' },
      { name: 'Problem 2.18, N₂O₅ (Ogg)', v: { net: 'N2O5 <=> NO2 + NO3 ; k = 1, kr = 50\nNO2 + NO3 -> NO + O2 + NO2 ; k = 0.5\nNO + NO3 -> NO2 + NO2 ; k = 1000', init: 'N2O5 = 1', tEnd: 200, prod: 'O2', law: 'k1*k2/(k_1 + 2*k2)*N2O5', consts: 'k1 = 1, k_1 = 50, k2 = 0.5', logy: true }, note: 'Ogg’s scheme gives first-order decomposition of N₂O₅.' }
    ],
    eqs: [
      ['\\frac{d[X]}{dt} \\approx 0 \\quad\\text{for each unmeasured intermediate X}', '§2.1', 'Type 1 intermediate: the steady-state approximation.'],
      ['[C_0] = [C] + [X]', '§2.1', 'Type 2: a catalyst present free and combined.'],
      ['-r_A = \\frac{k_3[E_0][A]}{[M]+[A]},\\qquad [M] = \\frac{k_2+k_3}{k_1}', 'Example 2.2']
    ],
    method: ['The mechanism is integrated without any approximation (Dormand–Prince, switching to the Rosenbrock ode23s method when the intermediates make the system stiff). The true rate of formation of the product is the sum over all steps.', 'Your proposed rate is evaluated on the simulated concentrations of the measured species. Where the two curves lie on top of each other after the brief start-up transient, the approximation is good.'],
    run(v) {
      const net = ST.parseNet(v.net);
      const init = ST.K.parseParams(v.init);
      const C0 = net.species.map(s => init[s] || 0);
      Object.keys(init).forEach(s => { if (net.species.indexOf(s) < 0) throw new Error('“' + s + '” is not in the mechanism'); });
      const ip = net.idx(String(v.prod).trim());
      if (ip < 0) throw new Error('Product “' + v.prod + '” is not in the mechanism');
      const consts = ST.K.parseParams(v.consts);
      const fn = ST.X.compile(v.law, net.species.concat(Object.keys(consts)));
      const cv = Object.values(consts);
      const tEnd = +v.tEnd;
      const sol = N.ode((t, C) => net.prod(C), 0, C0, tEnd, { rtol: 1e-8, atol: 1e-14, minPts: 100 });
      const ts = v.logy ? [0].concat(N.logspace(tEnd * 1e-5, tEnd, 160)) : N.linspace(0, tEnd, 200);
      const Cs = ts.map(sol.at);
      const trueR = Cs.map(C => net.prod(C)[ip]);
      const propR = Cs.map(C => fn(C.map(c => Math.max(c, 0)).concat(cv)));
      const series = net.species.map((s, i) => ({ name: s, x: ts, y: Cs.map(C => C[i]), color: i + 1 }));
      const rr = [{ name: 'true rate (full mechanism)', x: ts, y: trueR, color: 1, w: 2.6 }, { name: 'your rate law', x: ts, y: propR, color: 2, dash: '6 4' }];
      const peak = Math.max(...trueR.map(Math.abs));
      let dev = 0;
      ts.forEach((t, i) => { if (t > tEnd * 0.05 && Math.abs(trueR[i]) > 1e-3 * peak && isFinite(propR[i])) dev = Math.max(dev, Math.abs(propR[i] - trueR[i]) / Math.abs(trueR[i])); });
      const kpi = [
        { k: 'dev', l: 'Largest relative gap after start-up', v: dev * 100, u: '%', main: true, state: dev < 0.02 ? '' : dev < 0.1 ? 'warn' : 'bad', s: 'after 5% of the run, where the rate is above 0.1% of its peak' },
        { k: 'conv', l: 'Product at end', v: sol.last()[ip], u: v._cu }
      ];
      const reac = net.species.findIndex((s, i) => C0[i] > 0);
      if (reac >= 0) {
        const pairs = Cs.map((C, i) => [C[reac], trueR[i]]).filter(p => p[0] > 0 && p[1] > 0 && p[0] < C0[reac] * 0.95);
        if (pairs.length > 4) {
          const lr = N.linreg(pairs.map(p => Math.log(p[0])), pairs.map(p => Math.log(p[1])));
          kpi.push({ k: 'order', l: 'Apparent order in ' + net.species[reac], v: lr.b, s: 'slope of ln rate vs ln C over the run' });
        }
      }
      const pos = v.logy ? { log: true, min: tEnd * 1e-5 } : { min: 0 };
      return {
        kpi,
        plots: [
          { title: 'Concentrations in the batch', x: Object.assign({ label: 't', unit: v._tu }, pos), y: { label: 'C', unit: v._cu, log: !!v.logy }, series: series.map(s => (v.logy ? Object.assign(s, { x: s.x.slice(1), y: s.y.slice(1) }) : s)) },
          { title: 'Rate of formation of ' + v.prod, x: Object.assign({ label: 't', unit: v._tu }, pos), y: { label: 'r', unit: v._cu + '/' + v._tu }, series: rr.map(s => (v.logy ? Object.assign(s, { x: s.x.slice(1), y: s.y.slice(1) }) : s)) }
        ],
        tables: [{ title: 'Steps as read', cols: ['step', 'forward k', 'reverse k'], rows: net.rx.map(r => [r.text, r.custom ? 'custom' : r.k, r.rev ? r.kr : '—']) }]
      };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
