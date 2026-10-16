/* Appendix — Miscellany: units and conversion factors (A–L), the gas constant, orders of magnitude of
   diffusivity, concentration, conductivity, transfer coefficients and rates (K–P), and the
   dimensionless groups with the packed-bed correlations (N, O, Q). */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N;
  const fmt = (x, s = 4) => N.fmt(x, s);

  /* factors: how many of each unit make one SI unit (the appendix's own numbers where it gives them) */
  const U = {
    len: { t: 'Length', si: 'm', u: [['Å', 1e10], ['μm', 1e6], ['mm', 1e3], ['cm', 100], ['m', 1], ['in', 39.37], ['ft', 3.28083], ['mile', 0.0006214]] },
    vol: { t: 'Volume', si: 'm³', u: [['in³', 61023], ['fluid oz', 33817.6], ['L', 1000], ['US gal', 264.2], ['Imp. gal', 220.2], ['ft³', 35.318], ['bbl (oil, 42 US gal)', 6.29], ['drum (55 US gal)', 4.80], ['m³', 1]] },
    mass: { t: 'Mass', si: 'kg', u: [['oz (avoirdupois)', 35.27], ['lb', 2.205], ['g', 1000], ['kg', 1], ['short ton (2000 lb)', 0.0011025], ['metric ton', 0.001], ['long ton (2240 lb)', 0.0009842]] },
    force: { t: 'Force', si: 'N', u: [['N', 1], ['dyne', 1e5], ['kg_f', 0.101972], ['lb_f', 0.224809]] },
    pres: { t: 'Pressure', si: 'Pa', def: 'atm', u: [['Pa', 1], ['kPa', 1e-3], ['bar', 1e-5], ['atm', 1 / 101325], ['mm Hg', 760 / 101325], ['psi (lb_f/in²)', 14.696 / 101325], ['in Hg', 29.92 / 101325], ['ft H₂O', 33.93 / 101325], ['in H₂O', 1 / 248.86], ['dyne/cm²', 10]] },
    energy: { t: 'Work, energy and heat', si: 'J', def: 'kcal', u: [['erg', 1e7], ['J', 1], ['ft·lb_f', 0.737562], ['cal', 0.238846], ['kg_f·m', 0.101972], ['L·atm', 0.009869], ['Btu', 0.0009478], ['kcal', 0.00023885], ['hp·hr', 3.72506e-7], ['kW·hr', 2.77778e-7]] },
    R: { t: 'Gas constant R', si: 'J/mol·K', def: 'J/mol·K', u: [['J/mol·K', 1], ['cal/mol·K', 1.987 / 8.314], ['ft³·atm/lb mol·°R', 0.7302 / 8.314], ['L·atm/mol·K', 0.08206 / 8.314], ['Btu/lb mol·°R', 1.987 / 8.314], ['m³·atm/mol·K', 82.06e-6 / 8.314], ['Pa·m³/mol·K', 1], ['Pa·L/mol·K', 1000]] },
    mw: { t: 'Molecular weight', si: 'kg/mol', def: 'g/mol', u: [['kg/mol', 1], ['g/mol', 1000], ['lb/lb mol', 1000]] },
    visc: { t: 'Viscosity μ', si: 'kg/m·s', def: 'cP', u: [['Pl = kg/m·s = Pa·s', 1], ['poise', 10], ['cP', 1000], ['lb_m/ft·s', 0.672], ['lb_m/ft·hr', 2420]] },
    dens: { t: 'Density ρ', si: 'kg/m³', u: [['kg/m³', 1], ['g/cm³', 1e-3], ['lb_m/ft³', 0.062428]] },
    diff: { t: 'Diffusivity 𝒟', si: 'm²/s', def: 'cm²/s', u: [['m²/s', 1], ['cm²/s', 1e4], ['ft²/hr', 3.875e4]] },
    conc: { t: 'Concentration C_A', si: 'mol/m³', u: [['mol/m³', 1], ['mol/L', 1e-3], ['kmol/m³', 1e-3], ['lb mol/ft³', 6.24e-5]] },
    cond: { t: 'Thermal conductivity k', si: 'W/m·K', u: [['W/m·K', 1], ['cal/m·K·s', 0.239], ['Btu/hr·ft·°F', 0.578]] },
    htc: { t: 'Heat transfer coefficient h', si: 'W/m²·K', u: [['W/m²·K', 1], ['cal/m²·K·s', 0.239], ['Btu/hr·ft²·°F', 0.176]] },
    mtc: { t: 'Mass transfer coefficient k_g', si: 'm³ gas/m² surface·s', u: [['m/s', 1], ['cm/s', 100], ['ft/hr', 11811]] },
    rate: { t: 'Rate of reaction −r_A‴', si: 'mol/m³·s', u: [['mol/m³·s', 1], ['mol/L·s', 1e-3], ['mol/L·min', 0.06], ['kmol/m³·hr', 3.6], ['lb mol/ft³·hr', 6.24e-5 * 3600]] },
    temp: { t: 'Temperature', si: 'K', def: '°C', u: [['K'], ['°C'], ['°F'], ['°R']] }
  };
  const toK = { K: x => x, '°C': x => x + 273.15, '°F': x => (x - 32) / 1.8 + 273.15, '°R': x => x / 1.8 };
  const fromK = { K: x => x, '°C': x => x - 273.15, '°F': x => (x - 273.15) * 1.8 + 32, '°R': x => x * 1.8 };
  /* orders of magnitude, after the appendix's charts: [label, low, high] or [label, value] in SI */
  const SPEC = {
    diff: { t: 'Diffusivity (Appendix K)', r: [['in liquids', 1e-10, 1e-8], ['gas in porous media, Knudsen regime or small pores', 1e-8, 5e-6], ['commercial catalyst', 1e-7, 2e-6], ['gas in porous media, large pores or bulk diffusion', 5e-7, 1e-5], ['in gas at 1 atm', 5e-6, 1e-4]] },
    conc: { t: 'Concentration (Appendix L)', r: [['1% in gas at 1 atm, 1000 K', 0.01 * 101325 / (8.314 * 1000)], ['pure gas at 1 atm, 273 K', 101325 / (8.314 * 273.15)], ['1 N aqueous', 1000], ['pure gas at 1000 atm, 273 K (ideal)', 1000 * 101325 / (8.314 * 273.15)], ['pure water', 1000 / 0.018015]] },
    cond: { t: 'Thermal conductivity (Appendix M)', r: [['gases', 0.007, 0.05], ['insulators', 0.02, 0.07], ['liquids and H₂ gas', 0.08, 0.7], ['water', 0.6], ['nonmetals, porous structures, alumina, silica, activated carbon', 0.05, 5], ['metals', 10, 420]] },
    htc: { t: 'Heat transfer coefficient (Appendix N)', r: [['gas to particle', 8, 1200], ['gas to fine entrained particles (FCC, fast fluidized)', 1000, 1200], ['liquid to particle', 80, 1200]] },
    mtc: { t: 'Mass transfer coefficient (Appendix O)', r: [['liquid to particle', 2e-7, 2e-5], ['gas to particle', 0.02, 2]] },
    rate: { t: 'Rate of reaction (Appendix P)', r: [['cellular reactions, industrial water treatment', 3e-5, 3e-4], ['human being at rest', 2e-3], ['working hard', 6e-3], ['gases in porous catalyst particles', 0.1, 10], ['coal furnaces', 3, 60], ['jet engines', 100, 2000], ['rocket engines', 3000, 3e4], ['bimolecular reaction, every collision counting, 1 atm, 400 °C', 1e9, 3e10]] },
    visc: { t: 'Viscosity (Appendix I)', r: [['gases', 1e-5], ['air at 20 °C', 1.8e-5], ['water at 20 °C', 1e-3]] },
    dens: { t: 'Density (Appendix J)', r: [['air at 20 °C, 1 atm', 1.2], ['water', 1000]] }
  };

  /* unit names stay upright: an empty $$ switches the light-math lexer to explicit mode */
  const up = u => u + '$$';
  function spectrum(key, xSI, unit) {
    const S = SPEC[key];
    const rows = S.r.length;
    const lo = Math.min(...S.r.map(r => r[1]), xSI > 0 ? xSI : Infinity), hi = Math.max(...S.r.map(r => r[2] || r[1]), xSI > 0 ? xSI : 0);
    const x0 = Math.pow(10, Math.floor(Math.log10(lo)) - 0.3), x1 = Math.pow(10, Math.ceil(Math.log10(hi)) + 0.3);
    const series = S.r.map((r, i) => {
      const y = rows - i, c = (i % 7) + 2;
      return r.length === 3 ? { name: r[0], x: [r[1], r[2]], y: [y, y], color: c, w: 9, noRange: true } : { name: r[0], x: [r[1]], y: [y], type: 'pts', color: c, r: 6, noRange: true };
    });
    const marks = xSI > 0 ? [{ t: 'v', x: xSI, color: 1, label: 'your value', dash: false }] : [];
    return { title: S.t + ', ' + unit, x: { label: up(unit), log: true, min: x0, max: x1 }, y: { label: '', min: 0.3, max: rows + 0.7, ticks: false }, series, marks, hmax: 60 + 30 * rows, cap: 'Typical ranges read from the appendix charts: orders of magnitude, not data.' };
  }

  ST.tool({
    id: 'units', tag: 'A-01', ch: 31, sec: 'A–M, P', title: 'Units, conversion factors and orders of magnitude',
    lede: 'Every conversion table of the appendix in one place: length, volume, mass, force, pressure, energy, the gas constant in eight sets of units, molecular weight, viscosity, density, diffusivity, concentration, conductivity, heat and mass transfer coefficients, reaction rate and temperature. Type a value in any unit and see it in all the others, then see where it falls among the typical magnitudes the appendix charts.',
    keys: 'units conversion factors SI english angstrom micron inch foot mile liter gallon barrel drum pound ton pascal atm bar psi mm hg inch water joule calorie btu kcal kwh gas constant R viscosity poise centipoise poiseuille density diffusivity concentration thermal conductivity heat transfer coefficient mass transfer coefficient rate of reaction temperature kelvin rankine fahrenheit orders of magnitude appendix miscellany newton gc',
    inputs: [
      { k: 'qty', t: 'select', label: 'Quantity', opts: Object.keys(U).map(k => [k, U[k].t]), v: 'pres' },
      { k: 'val', label: 'Value', v: 1, slider: false }
    ].concat(Object.keys(U).map(k => ({ k: 'u_' + k, t: 'select', label: 'Given in', opts: U[k].u.map(r => [r[0], r[0]]), v: U[k].def || U[k].u[0][0], if: v => v.qty === k }))),
    presets: [
      { name: 'One atmosphere in every unit', v: { qty: 'pres', val: 1, u_pres: 'atm' } },
      { name: 'Air at 20 °C: ρ = 1.20 kg/m³ (Appendix J)', v: { qty: 'dens', val: 1.2, u_dens: 'kg/m³' } },
      { name: 'R in every unit (Appendix H)', v: { qty: 'R', val: 8.314, u_R: 'J/mol·K' } },
      { name: 'Diffusivity of a gas, 0.1 cm²/s', v: { qty: 'diff', val: 0.1, u_diff: 'cm²/s' } },
      { name: 'Rate in a catalyst pellet, 1 mol/L·min', v: { qty: 'rate', val: 1, u_rate: 'mol/L·min' } },
      { name: 'Normal body temperature, 98.6 °F', v: { qty: 'temp', val: 98.6, u_temp: '°F' } }
    ],
    eqs: [
      ['F = \\frac{ma}{g_c},\\qquad g_c = 1\\,\\mathrm{\\frac{kg\\cdot m}{s^2\\cdot N}}', 'Appendix A', 'At the earth’s surface a = g = 9.806 m/s².'],
      ['1\\ \\mathrm{atm} = 760\\ \\mathrm{mm\\,Hg} = 14.696\\ \\mathrm{lb_f/in^2} = 101\\,325\\ \\mathrm{Pa}', 'Appendix E', '1 bar = 10⁵ Pa, close to 1 atm.'],
      ['\\mathbf{R} = 8.314\\ \\mathrm{\\frac{J}{mol\\cdot K}} = 1.987\\ \\mathrm{\\frac{cal}{mol\\cdot K}} = 0.082\\,06\\ \\mathrm{\\frac{L\\cdot atm}{mol\\cdot K}}', 'Appendix H'],
      ['\\mathbf{R} = 0.7302\\ \\mathrm{\\frac{ft^3\\cdot atm}{lb\\,mol\\cdot{}^{\\circ}R}} = 1.987\\ \\mathrm{\\frac{Btu}{lb\\,mol\\cdot{}^{\\circ}R}}', 'Appendix H'],
      ['1\\ \\mathrm{Pl} = 1\\ \\mathrm{\\frac{kg}{m\\cdot s}} = 10\\ \\mathrm{poise} = 1000\\ \\mathrm{cP}', 'Appendix I']
    ],
    method: ['Each quantity is stored as the number of units that make one SI unit, with the appendix’s own factors where it prints them; a value goes to SI and back out to every unit. Temperatures convert through kelvin. The magnitude charts are bands read from the appendix figures.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const q = U[v.qty], from = v['u_' + v.qty], x = +v.val;
      if (!q) throw new Error('Pick a quantity.');
      if (v.qty === 'temp') {
        const K = toK[from](x);
        if (K < 0) warn.push('That is below absolute zero.');
        kpi.push({ k: 'K', l: 'Kelvin', v: K, u: 'K', main: true });
        ['°C', '°F', '°R'].forEach((u, i) => kpi.push({ k: 't' + i, l: up(u), v: fromK[u](K) }));
        tables.push({ title: 'Temperature', cols: ['unit', 'value'], rows: ['K', '°C', '°F', '°R'].map(u => [up(u), fromK[u](K)]) });
        info.push('Absolute temperatures (K or °R) go into the Arrhenius equation and the gas law; differences are the same in K and °C, and in °R and °F.');
        return { kpi, warn, info, plots, tables };
      }
      const f = q.u.find(r => r[0] === from);
      const si = x / f[1];
      kpi.push({ k: 'si', l: up('In SI, ' + q.si), v: si, main: true });
      q.u.filter(r => r[0] !== from && r[1] !== 1).slice(0, 5).forEach((r, i) => kpi.push({ k: 'c' + i, l: up(r[0]), v: si * r[1] }));
      tables.push({ title: up(q.t + ': ' + fmt(x, 6) + ' ' + from), cols: ['unit', 'value', up('one ' + from + ' is')], rows: q.u.map(r => [up(r[0]), si * r[1], r[1] / f[1]]) });
      if (SPEC[v.qty]) plots.push(spectrum(v.qty, si, q.si));
      if (v.qty === 'mw') info.push('In SI the molecular weight carries units: (mw)_{O₂} = 0.032 kg/mol, (mw)_{air} = 0.0289 kg/mol (Appendix G).');
      if (v.qty === 'diff') info.push('𝒟 does not depend on π in liquids; for bulk diffusion of gases 𝒟 ∝ T^{3/2}/π, and in the Knudsen regime 𝒟 ∝ T^{1/2}, independent of π (A-02 computes both).');
      if (v.qty === 'cond') info.push('k and k_eff do not depend on π. In gas or liquid k is per m of fluid; in porous structures k_eff is per m of structure.');
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* A-02  dimensionless groups and transfer coefficients                   */
  /* ===================================================================== */
  const GROUPS = [
    ['Re', 'Reynolds', 'duρ/μ', 'inertial / viscous effects'],
    ['Sc', 'Schmidt', 'μ/ρ𝒟', 'molecular momentum / mass transfer'],
    ['Pr', 'Prandtl', 'C_{p}μ/k', 'molecular momentum / heat transfer'],
    ['Nu', 'Nusselt', 'hd/k', 'total / molecular heat transfer'],
    ['Sh', 'Sherwood', 'k_{g}d/𝒟', 'total / molecular mass transfer'],
    ['Pe', 'Peclet', 'duρC_{p}/k = Re·Pr', 'total momentum / molecular heat transfer'],
    ['Bo', 'Bodenstein', 'ud/𝒟 = Re·Sc', 'total momentum / molecular mass transfer'],
    ['Dud', 'intensity of axial dispersion', 'D/ud', 'axial dispersion / bulk flow'],
    ['DuL', 'vessel dispersion number', 'D/uL', 'the same, over the vessel']
  ];
  ST.tool({
    id: 'groups', tag: 'A-02', ch: 31, sec: 'G–Q', title: 'Dimensionless groups, gas properties and transfer coefficients',
    lede: 'From a fluid’s properties and a flow, all the groups of the appendix — Re, Sc, Pr, Nu, Sh, Pe, Bodenstein and the dispersion numbers — with what each one measures; heat and mass transfer coefficients for particles in packed beds from Nu, Sh = 2 + 1.8Re^{1/2}(Pr, Sc)^{1/3}; the ideal-gas density; and whether a gas diffuses in a pore by bulk or Knudsen diffusion, with 𝒟 = 1.534d√(T/(mw)).',
    keys: 'dimensionless groups reynolds schmidt prandtl nusselt sherwood peclet bodenstein dispersion number intensity of axial dispersion vessel dispersion number packed bed heat transfer coefficient mass transfer coefficient particle correlation ideal gas density knudsen diffusion bulk diffusion pore diameter molecular weight appendix',
    inputs: [
      { group: 'Fluid' },
      { k: 'ideal', t: 'check', label: 'Gas: density from the ideal-gas law', v: true },
      { k: 'p', label: 'Pressure π', unit: 'Pa', v: 101325, min: 1e-3, max: 1e9, scale: 'log' },
      { k: 'T', label: 'Temperature', unit: 'K', v: 293, min: 1, max: 5000, scale: 'log' },
      { k: 'mw', label: 'Molecular weight (mw)', unit: 'kg/mol', v: 0.0289, min: 1e-4, max: 10, scale: 'log' },
      { k: 'rho', label: 'Density ρ', unit: 'kg/m³', v: 1000, min: 1e-6, max: 1e5, scale: 'log', if: v => !v.ideal },
      { k: 'mu', label: 'Viscosity μ', unit: 'kg/m·s', v: 1.8e-5, min: 1e-7, max: 1e3, scale: 'log' },
      { k: 'Dm', label: 'Molecular diffusivity 𝒟', unit: 'm²/s', v: 2e-5, min: 1e-14, max: 1e-2, scale: 'log' },
      { k: 'Cp', label: 'Heat capacity C_p', unit: 'J/kg·K', v: 1005, min: 1, max: 1e5, scale: 'log' },
      { k: 'kf', label: 'Thermal conductivity k', unit: 'W/m·K', v: 0.0257, min: 1e-4, max: 1e3, scale: 'log' },
      { group: 'Flow' },
      { k: 'd', label: 'Particle or tube diameter d', unit: 'm', v: 0.005, min: 1e-9, max: 10, scale: 'log' },
      { k: 'u', label: 'Velocity u', unit: 'm/s', v: 1, min: 1e-9, max: 1e3, scale: 'log' },
      { k: 'Dax', label: 'Axial dispersion coefficient D (0 if unknown)', unit: 'm²/s', v: 0, min: 0, max: 100, slider: false },
      { k: 'L', label: 'Vessel length L', unit: 'm', v: 1, min: 1e-6, max: 1e4, scale: 'log', if: v => +v.Dax > 0 },
      { group: 'Pores (Appendix K)' },
      { k: 'dp', label: 'Pore diameter', unit: 'm', v: 1e-8, min: 1e-11, max: 1e-2, scale: 'log' }
    ],
    presets: [
      { name: 'Air at 20 °C, 1 atm (Appendix J)', v: { ideal: true, p: 101325, T: 293, mw: 0.0289, mu: 1.8e-5, Dm: 2e-5, Cp: 1005, kf: 0.0257, d: 0.005, u: 1, Dax: 0, dp: 1e-8 }, book: { rho: { v: 1.2, tol: 0.005 } }, check: { Pr: { v: 0.70, tol: 0.08 } }, note: 'The book: ρ = (101 325)(0.0289)/(8.314)(293) = 1.20 kg/m³, and Pr = 0.66 to 0.75 for air and the common gases.' },
      { name: 'Typical gas: Sc ≅ 1 (Appendix Q)', v: { ideal: false, rho: 1, T: 293, mu: 1e-5, Dm: 1e-5, Cp: 1000, kf: 0.025, d: 0.005, u: 1, Dax: 0, dp: 1e-8 }, book: { Sc: { v: 1, tol: 1e-9 } } },
      { name: 'Typical liquid: Sc ≅ 10³ (Appendix Q)', v: { ideal: false, rho: 1000, T: 293, mu: 1e-3, Dm: 1e-9, Cp: 4180, kf: 0.6, d: 0.001, u: 0.01, Dax: 0, dp: 1e-8 }, book: { Sc: { v: 1000, tol: 1e-9 } }, check: { Pr: { v: 7, tol: 0.05 } } },
      { name: 'Knudsen diffusion in a 10 nm pore', v: { ideal: true, p: 101325, T: 600, mw: 0.028, mu: 3e-5, Dm: 6e-5, Cp: 1100, kf: 0.045, d: 0.003, u: 0.5, Dax: 0, dp: 1e-8 }, expect: { DK: { v: 1.534e-8 * Math.sqrt(600 / 0.028), tol: 1e-9 } } },
      { name: 'Packed bed with measured dispersion', v: { ideal: true, p: 101325, T: 500, mw: 0.03, mu: 2.6e-5, Dm: 4e-5, Cp: 1050, kf: 0.04, d: 0.004, u: 0.4, Dax: 0.002, L: 2, dp: 2e-7 } }
    ],
    eqs: [
      ['\\rho = \\frac{\\pi(mw)}{\\mathbf{R}T}', 'Appendix J'],
      ['\\mathscr{D}_K = 1.534\\,d\\sqrt{T/(mw)}', 'Appendix K', 'Knudsen regime when πd < 0.01 Pa·m; d in m, T in K, (mw) in kg/mol.'],
      ['Nu = \\frac{hd_p}{k} = 2 + 1.8\\,Re_p^{1/2}Pr^{1/3}', 'Appendix N', 'Packed beds, Re_p > 100.'],
      ['Sh = \\frac{k_gd_p}{\\mathscr{D}} = 2 + 1.8\\,Re_p^{1/2}Sc^{1/3}', 'Appendix O', 'Packed beds, Re_p > 80.'],
      ['Pe = Re\\,Pr,\\qquad Bo = Re\\,Sc = \\frac{ud}{\\mathscr{D}}', 'Appendix Q'],
      ['\\frac{\\mathbf{D}}{ud},\\qquad \\frac{\\mathbf{D}}{uL}', 'Appendix Q', 'Not the reciprocal of a Peclet or Bodenstein number: D is a dispersion coefficient, not a molecular diffusivity.']
    ],
    method: ['Straight evaluation of the appendix’s definitions and correlations. The correlations hold above the Reynolds numbers shown; below them the tool still reports the numbers but flags them.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const T = +v.T, mu = +v.mu, Dm = +v.Dm, Cp = +v.Cp, kf = +v.kf, d = +v.d, u = +v.u;
      const rho = v.ideal ? +v.p * +v.mw / (8.314 * T) : +v.rho;
      if (v.ideal) kpi.push({ k: 'rho', l: 'Gas density π(mw)/RT', v: rho, u: 'kg/m³' });
      const Re = d * u * rho / mu, Sc = mu / (rho * Dm), Pr = Cp * mu / kf;
      const Nu = 2 + 1.8 * Math.sqrt(Re) * Math.cbrt(Pr), Sh = 2 + 1.8 * Math.sqrt(Re) * Math.cbrt(Sc);
      const h = Nu * kf / d, kg = Sh * Dm / d;
      kpi.push({ k: 'Re', l: 'Reynolds number du ρ/μ', v: Re, main: true });
      kpi.push({ k: 'Sc', l: 'Schmidt number μ/ρ𝒟', v: Sc, s: Sc < 10 ? 'gas-like (≅ 1)' : Sc > 100 ? 'liquid-like (≅ 10³)' : '' });
      kpi.push({ k: 'Pr', l: 'Prandtl number C_pμ/k', v: Pr });
      kpi.push({ k: 'h', l: 'Particle h, packed bed', v: h, u: 'W/m²·K', s: 'Nu = ' + fmt(Nu, 4) + (Re < 100 ? ' (Re_p below 100)' : '') });
      kpi.push({ k: 'kg', l: 'Particle k_g, packed bed', v: kg, u: 'm/s', s: 'Sh = ' + fmt(Sh, 4) + (Re < 80 ? ' (Re_p below 80)' : '') });
      if (Re < 80) warn.push('Re_p = ' + fmt(Re, 3) + ': the packed-bed correlations are given for Re_p above 80–100; at low Re they tend to the single-sphere limit Nu = Sh = 2.');
      const Pe = Re * Pr, Bo = Re * Sc;
      const vals = { Re, Sc, Pr, Nu, Sh, Pe, Bo };
      if (+v.Dax > 0) { vals.Dud = +v.Dax / (u * d); vals.DuL = +v.Dax / (u * +v.L); kpi.push({ k: 'DuL', l: 'Vessel dispersion number D/uL', v: vals.DuL, s: vals.DuL < 0.01 ? 'small: close to plug flow' : vals.DuL > 0.2 ? 'large: far from plug flow' : '' }); }
      /* pores */
      const pd = +v.p * +v.dp;
      const DK = 1.534 * +v.dp * Math.sqrt(T / +v.mw);
      kpi.push({ k: 'DK', l: 'Knudsen diffusivity in the pore', v: DK, u: 'm²/s', s: pd < 0.01 ? 'πd = ' + fmt(pd, 3) + ' Pa·m < 0.01: Knudsen regime' : 'πd = ' + fmt(pd, 3) + ' Pa·m: bulk diffusion dominates' });
      tables.push({ title: 'Dimensionless groups (Appendix Q)', cols: ['group', 'name', 'definition', 'value', 'what it compares'], rows: GROUPS.filter(g => vals[g[0]] !== undefined).map(g => [g[0] === 'Dud' ? 'D/ud' : g[0] === 'DuL' ? 'D/uL' : g[0], g[1], g[2], vals[g[0]], g[3]]) });
      /* where h and kg fall */
      const re = N.logspace(1, 1e5, 120);
      plots.push({ title: 'Packed-bed correlations (Appendix N, O)', x: { label: 'Re_p', log: true }, y: { label: 'Nu or Sh', log: true },
        series: [{ name: 'Sh, gases (Sc = 1)', x: re, y: re.map(r => 2 + 1.8 * Math.sqrt(r)), color: 2 }, { name: 'Sh, liquids (Sc = 10³)', x: re, y: re.map(r => 2 + 18 * Math.sqrt(r)), color: 3 }, { name: 'Nu, this fluid', x: re, y: re.map(r => 2 + 1.8 * Math.sqrt(r) * Math.cbrt(Pr)), color: 4, dash: '4 4' }],
        marks: [{ t: 'pt', x: Re, y: Sh, color: 1, label: 'Sh here' }, { t: 'pt', x: Re, y: Nu, color: 4, label: 'Nu here', dy: 14 }, { t: 'band', x0: 10, x1: 80, color: 'var(--grid)' }],
        cap: 'Shaded: below the range the correlations were fitted to.' });
      plots.push(spectrum('htc', h, 'W/m²·K'));
      plots.push(spectrum('mtc', kg, 'm/s'));
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
