/* Chapter 25 — fluid–particle reactions: kinetics. The shrinking-core model for flat plates,
   cylinders and spheres (film, ash, reaction and all combined), shrinking particles with film and
   reaction in series, finding the controlling step from runs, and the combustion chart (Fig. 25.12). */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N;
  const fmt = (x, s = 4) => N.fmt(x, s);
  const FP = ST.FP = {};
  const Rg = 8.314462618;

  /* t/τ as a function of X_B for each controlling step (Table 25.1) */
  FP.g = function (mech, shape, X) {
    X = Math.min(Math.max(X, 0), 1);
    const y = 1 - X;
    switch (mech) {
      case 'film': return X;
      case 'ash': return shape === 'slab' ? X * X : shape === 'cyl' ? X + (y > 0 ? y * Math.log(y) : 0) : 1 - 3 * Math.cbrt(y * y) + 2 * y;
      case 'rxn': return shape === 'slab' ? X : shape === 'cyl' ? 1 - Math.sqrt(y) : 1 - Math.cbrt(y);
      case 'stokes': return 1 - Math.cbrt(y * y);
      case 'large': return 1 - Math.sqrt(y);
      case 'srxn': return 1 - Math.cbrt(y);
      default: return NaN;
    }
  };
  /* X_B from t/τ */
  FP.Xof = function (mech, shape, s) {
    if (!(s > 0)) return 0;
    if (s >= 1) return 1;
    switch (mech) {
      case 'film': return s;
      case 'rxn': return shape === 'slab' ? s : shape === 'cyl' ? 1 - (1 - s) * (1 - s) : 1 - Math.pow(1 - s, 3);
      case 'stokes': return 1 - Math.pow(1 - s, 1.5);
      case 'large': return 1 - (1 - s) * (1 - s);
      case 'srxn': return 1 - Math.pow(1 - s, 3);
      case 'ash':
        if (shape === 'slab') return Math.sqrt(s);
        if (shape === 'cyl') return N.brent(X => FP.g('ash', 'cyl', X) - s, 0, 1, { tol: 1e-14 });
        { const u = N.brent(u => 1 - 3 * u * u + 2 * u * u * u - s, 0, 1, { tol: 1e-14 }); return 1 - u * u * u; }
      default: return NaN;
    }
  };
  /* resistances in series, constant size: t(X) = Σ τ_m g_m(X) (Eq. 25.32), X(t) by Brent */
  FP.combo = function (taus, shape) {
    const tX = X => (taus.film ? taus.film * FP.g('film', shape, X) : 0) + (taus.ash ? taus.ash * FP.g('ash', shape, X) : 0) + (taus.rxn ? taus.rxn * FP.g('rxn', shape, X) : 0);
    const tau = tX(1);
    const Xt = t => (t <= 0 ? 0 : t >= tau ? 1 : N.brent(X => tX(X) - t, 0, 1, { tol: 1e-13 }));
    return { tX, Xt, tau };
  };
  FP.MECH = { film: 'gas film', ash: 'ash layer', rxn: 'reaction', stokes: 'gas film, Stokes regime (shrinking)', large: 'gas film, large particles (shrinking)', srxn: 'reaction (shrinking particle)' };
  FP.EXPO = { film: [1.5, 2], ash: [2, 2], rxn: [1, 1], stokes: [2, 2], large: [1.5, 1.5], srxn: [1, 1] };
  /* geometry factor in τ: plate 1, cylinder 2, sphere 3 (film); 2, 4, 6 (ash) */
  FP.F = { slab: 1, cyl: 2, sph: 3 };
  /* air at 1 atm, for the combustion chart */
  FP.air = T => { const mu = 1.716e-5 * Math.pow(T / 273.15, 1.5) * (273.15 + 110.4) / (T + 110.4); return { mu, rho: 101325 * 0.028965 / (Rg * T), D: 0.181e-4 * Math.pow(T / 273.15, 1.75) }; };
  /* Froessling, Eq. 25.24 */
  FP.kgF = (d, u, rho, mu, D, y = 1) => D / (d * y) * (2 + 0.6 * Math.cbrt(mu / (rho * D)) * Math.sqrt(Math.max(d * u * rho / mu, 0)));

  /* ===================================================================== */
  /* N-2501  single particles: shrinking core and shrinking sphere          */
  /* ===================================================================== */
  ST.tool({
    id: 'scm', tag: 'N-2501', ch: 25, sec: 'Ch. 25', title: 'Fluid–particle kinetics: shrinking core and shrinking particles',
    lede: 'How long a single particle takes to react with a surrounding fluid, by the shrinking-core model: flat plates, cylinders and spheres with gas film, ash layer and surface reaction acting alone or all together, and particles that shrink as they burn with film and reaction in series (Froessling’s film coefficient updated as they shrink). The tool gives every τ, the conversion–time curve, how the resistances shift as the core recedes — and from your runs at one or several sizes it finds which step controls, as Chapter 25 asks.',
    keys: 'fluid particle noncatalytic gas solid shrinking core model unreacted core progressive conversion ash layer diffusion gas film reaction control time for complete conversion tau shrinking particle stokes regime froessling sphere cylinder flat plate combination of resistances controlling step particle size temperature carbon combustion figure 25.9 25.10 25.11 25.12 table 25.1 problem 25.1 25.3 25.4 25.5 25.6 25.7 25.8 25.9 25.10 25.11 25.12 zinc blende roasting graphite burning uranium oxide coal pile',
    units: { t: 's' },
    inputs: [
      { group: 'Task' },
      { k: 'mode', t: 'seg', label: 'Do', opts: [['part', 'one particle from its properties'], ['data', 'controlling step from runs'], ['fig', 'carbon combustion chart (Fig. 25.12)']], v: 'part' },
      { group: 'Particle', if: v => v.mode !== 'fig' },
      { k: 'beh', t: 'seg', label: 'During reaction the particle', opts: [['const', 'keeps its size (ash or firm product)'], ['shrink', 'shrinks (no ash)']], v: 'const', if: v => v.mode !== 'fig' },
      { k: 'shape', t: 'seg', label: 'Shape', opts: [['slab', 'flat plate'], ['cyl', 'cylinder'], ['sph', 'sphere']], v: 'sph', if: v => v.mode !== 'fig' && v.beh === 'const' },
      { k: 'R', label: v => (v.shape === 'slab' && v.beh === 'const' ? 'Half thickness L' : 'Radius R'), unit: 'mm', v: 1, min: 1e-6, max: 1e6, scale: 'log', if: v => v.mode === 'part' },
      { k: 'rhoB', label: 'Molar density of B in the solid ρ_B', unit: 'mol/m³', v: 42500, min: 1e-6, max: 1e7, scale: 'log', if: v => v.mode === 'part' },
      { k: 'b', label: 'b, mol B per mol A', v: 2 / 3, min: 0.01, max: 20, if: v => v.mode === 'part' },
      { group: 'Fluid', if: v => v.mode === 'part' },
      { k: 'cgiv', t: 'seg', label: 'C_Ag from', opts: [['C', 'concentration'], ['yPT', 'mole fraction, π and T']], v: 'yPT', if: v => v.mode === 'part' },
      { k: 'CA', label: 'C_Ag', unit: 'mol/m³', v: 1, min: 1e-12, max: 1e6, if: v => v.mode === 'part' && v.cgiv === 'C' },
      { k: 'yA', label: 'Mole fraction of A', v: 0.08, min: 1e-9, max: 1, if: v => v.mode === 'part' && v.cgiv === 'yPT' },
      { k: 'P', label: 'Pressure π', unit: 'atm', v: 1, min: 1e-6, max: 1e4, if: v => v.mode === 'part' && v.cgiv === 'yPT' },
      { k: 'TC', label: 'Temperature', unit: '°C', v: 900, min: -270, max: 4000, if: v => v.mode === 'part' && v.cgiv === 'yPT' },
      { group: 'Resistances', if: v => v.mode === 'part' },
      { k: 'uf', t: 'check', label: 'Gas film', v: false, if: v => v.mode === 'part' },
      { k: 'kgs', t: 'seg', label: 'k_g', opts: [['given', 'given'], ['fro', 'Froessling (Eq. 25.24)']], v: 'given', if: v => v.mode === 'part' && v.uf },
      { k: 'kg', label: 'k_g', unit: 'cm/s', v: 10, min: 1e-9, max: 1e9, scale: 'log', if: v => v.mode === 'part' && v.uf && v.kgs === 'given' },
      { k: 'u', label: 'Relative velocity of gas and particle u', unit: 'cm/s', v: 10, min: 0, max: 1e6, if: v => v.mode === 'part' && v.uf && v.kgs === 'fro' },
      { k: 'Dg', label: 'Diffusivity of A in the gas D', unit: 'cm²/s', v: 2, min: 1e-9, max: 1e4, scale: 'log', if: v => v.mode === 'part' && v.uf && v.kgs === 'fro' },
      { k: 'nu', label: 'Kinematic viscosity of the gas μ/ρ', unit: 'cm²/s', v: 2, min: 1e-9, max: 1e4, scale: 'log', if: v => v.mode === 'part' && v.uf && v.kgs === 'fro' },
      { k: 'ua', t: 'check', label: 'Ash layer', v: true, if: v => v.mode === 'part' && v.beh === 'const' },
      { k: 'De', label: 'Effective diffusivity in the ash D_e', unit: 'cm²/s', v: 0.08, min: 1e-12, max: 1e4, scale: 'log', if: v => v.mode === 'part' && v.beh === 'const' && v.ua },
      { k: 'ur', t: 'check', label: 'Surface reaction', v: true, if: v => v.mode === 'part' },
      { k: 'kr', label: 'k″ (first order in A)', unit: 'cm/s', v: 2, min: 1e-12, max: 1e9, scale: 'log', if: v => v.mode === 'part' && v.ur },
      { k: 'tq', label: 'Time of interest', unit: '[t]', v: 1800, min: 0, max: 1e12, if: v => v.mode === 'part' },
      { group: 'Runs', if: v => v.mode === 'data' },
      { k: 'meas', t: 'seg', label: 'Measured', opts: [['X', 'conversion X_B'], ['r', 'core (or particle) radius r/R']], v: 'X', if: v => v.mode === 'data' },
      { k: 'runs', t: 'table', label: 'Runs: size, time, measurement', cols: ['size', 't', 'X_B or r/R'], v: [[1, 4, 1], [1.5, 6, 1]], if: v => v.mode === 'data' },
      { k: 'doneEq', t: 'check', label: 'A run with X_B = 1 gives the time for complete conversion (τ), not just an upper bound', v: true, if: v => v.mode === 'data' },
      { group: 'Chart', if: v => v.mode === 'fig' },
      { k: 'Ek', label: 'Activation energy of the surface reaction', unit: 'kJ/mol', v: 184, min: 0, max: 1000, if: v => v.mode === 'fig' },
      { k: 'Ak', label: 'k″ = A √T e^{−E/RT}: A', unit: 'cm/(s·K^{1/2})', v: 4.39e7, min: 1e-6, max: 1e20, scale: 'log', if: v => v.mode === 'fig', hint: 'The default puts k″ through Fig. 25.12 (1 cm/s near 1050 K).' },
      { k: 'yO2', label: 'Mole fraction of O₂', v: 0.21, min: 1e-6, max: 1, if: v => v.mode === 'fig' }
    ],
    presets: [
      { name: 'Problem 25.4 — roasting zinc blende, R = 1 mm', v: { mode: 'part', beh: 'const', shape: 'sph', R: 1, rhoB: 42500, b: 2 / 3, cgiv: 'yPT', yA: 0.08, P: 1, TC: 900, uf: false, ua: true, De: 0.08, ur: true, kr: 2, tq: 1800, _tu: 's' }, expect: { tau: { v: 0.0425 * 0.1 / ((2 / 3) * 2 * 0.08 / (82.057 * 1173.15)) + 0.0425 * 0.01 / (6 * (2 / 3) * 0.08 * 0.08 / (82.057 * 1173.15)), tol: 1e-4 } }, note: '2ZnS + 3O₂ → 2ZnO + 2SO₂ in 8% O₂ at 900 °C and 1 atm, ρ_B = 0.0425 mol/cm³, k″ = 2 cm/s, D_e = 0.08 cm²/s; film resistance negligible. τ_total = ρ_BR/(bk″C) + ρ_BR²/(6bD_eC).' },
      { name: 'Problem 25.5 — the same with R doubled is 25.6', v: { mode: 'part', beh: 'const', shape: 'sph', R: 2, rhoB: 42500, b: 2 / 3, cgiv: 'yPT', yA: 0.08, P: 1, TC: 900, uf: false, ua: true, De: 0.08, ur: true, kr: 2, tq: 1800, _tu: 's' }, note: 'Doubling the size: ash diffusion (∝ R²) gains on reaction (∝ R).' },
      { name: 'Problem 25.3 — burning graphite, R₀ = 5 mm', v: { mode: 'part', beh: 'shrink', R: 5, rhoB: 2200 / 12.011 * 1000, b: 1, cgiv: 'yPT', yA: 0.08, P: 1, TC: 900, uf: false, ur: true, kr: 20, tq: 1800, _tu: 's' }, expect: { tau: { v: (2200 / 12.011 * 1000) * 0.005 / (0.2 * 0.08 * 101325 / (8.314462618 * 1173.15)), tol: 1e-6 } }, note: 'Graphite, ρ = 2.2 g/cm³, k″ = 20 cm/s, 8% O₂ at 900 °C; the film offers no resistance at high gas velocity. τ = ρ_B R₀/(k″C_Ag).' },
      { name: 'Film and reaction in series, shrinking char particle', v: { mode: 'part', beh: 'shrink', R: 0.5, rhoB: 150000, b: 1, cgiv: 'yPT', yA: 0.21, P: 1, TC: 1200, uf: true, kgs: 'fro', u: 100, Dg: 2.5, nu: 2.3, ur: true, kr: 300, tq: 1, _tu: 's' }, note: 'As the particle shrinks Froessling’s k_g rises, so control passes from the film toward the reaction (Fig. 25.11).' },
      { name: 'Problem 25.1 — 7/8 in one hour, all in two', v: { mode: 'data', beh: 'const', shape: 'sph', meas: 'X', runs: [[1, 1, 0.875], [1, 2, 1]], doneEq: true, _tu: 'h' }, note: 'Film diffusion would give X proportional to t; reaction and ash both pass through t/τ = ½ at X = 7/8 — the data cannot separate them.' },
      { name: 'Problem 25.7 — complete in 4 and 6 min', v: { mode: 'data', beh: 'const', shape: 'sph', meas: 'X', runs: [[1, 4, 1], [1.5, 6, 1]], _tu: 'min' }, check: { verdict: 'reaction' } },
      { name: 'Problem 25.8 — 30% at 2 s, 75% at 5 s', v: { mode: 'data', beh: 'const', shape: 'sph', meas: 'X', runs: [[1, 2, 0.3], [1, 5, 0.75]], _tu: 's' }, check: { verdict: 'gas film' } },
      { name: 'Problem 25.9 — complete in 200 and 450 s', v: { mode: 'data', beh: 'const', shape: 'sph', meas: 'X', runs: [[1, 200, 1], [1.5, 450, 1]], doneEq: true, _tu: 's' }, check: { verdict: 'ash layer' }, note: 'τ ∝ R²: ash diffusion, or film diffusion at the top of its range; with a firm product the ash controls.' },
      { name: 'Problem 25.10 — 2 mm at 87.5%, 1 mm complete, both at 1 s', v: { mode: 'data', beh: 'const', shape: 'sph', meas: 'X', runs: [[2, 1, 0.875], [1, 1, 1]], doneEq: true, _tu: 's' }, check: { verdict: 'reaction' }, note: 'Read as: the 1 mm particles take exactly 1 s to convert. Only τ ∝ R then fits.' },
      { name: 'Problem 25.11 — reducing UO₃', v: { mode: 'data', beh: 'const', shape: 'sph', meas: 'X', runs: [[1, 0.18, 0.45], [1, 0.347, 0.68], [1, 0.453, 0.8], [1, 0.567, 0.95], [1, 0.733, 0.98]], _tu: 'h' }, check: { verdict: 'reaction' } },
      { name: 'Problem 25.12 — a burning coal pile', v: { mode: 'data', beh: 'shrink', meas: 'r', runs: [[1, 24, 0.95]], _tu: 'h' }, note: 'The silhouette shrinks 5% in a day. If reaction at the surface controls (the pile is far too big for film control), R falls linearly: 14% of the mass goes in that day and the fire burns out in about 480 hours.' },
      { name: 'Figure 25.12 — burning carbon', v: { mode: 'fig', Ek: 184, Ak: 4.39e7, yO2: 0.21 } }
    ],
    eqs: [
      ['\\frac{t}{\\tau} = X_B,\\qquad \\tau = \\frac{\\rho_BR}{3bk_gC_{Ag}}', 'Eqs. 25.10–25.11', 'Gas film controls (sphere).'],
      ['\\frac{t}{\\tau} = 1-3(1-X_B)^{2/3}+2(1-X_B),\\qquad \\tau = \\frac{\\rho_BR^2}{6b\\mathscr{D}_eC_{Ag}}', 'Eqs. 25.17–25.18', 'Ash diffusion controls (sphere).'],
      ['\\frac{t}{\\tau} = 1-(1-X_B)^{1/3},\\qquad \\tau = \\frac{\\rho_BR}{bk^{\\prime\\prime}C_{Ag}}', 'Eqs. 25.22–25.23', 'Reaction controls (sphere).'],
      ['\\frac{t}{\\tau} = 1-(1-X_B)^{2/3},\\qquad \\tau = \\frac{\\rho_BR_0^2}{2bC_{Ag}\\mathscr{D}}', 'Eqs. 25.29–25.30', 'Shrinking sphere, Stokes regime.'],
      ['\\frac{k_gd_py}{\\mathscr{D}} = 2+0.6\\,Sc^{1/3}Re^{1/2}', 'Eq. 25.24', 'Froessling.'],
      ['t_{total} = t_{film\\ alone} + t_{ash\\ alone} + t_{reaction\\ alone}', 'Eq. 25.32'],
      ['-\\frac{dr_c}{dt} = \\cfrac{bC_A/\\rho_B}{\\cfrac{r_c^2}{R^2k_g}+\\cfrac{(R-r_c)r_c}{R\\mathscr{D}_e}+\\cfrac{1}{k^{\\prime\\prime}}}', 'Eq. 25.33b'],
      ['t\\propto R^{1.5\\text{ to }2.0}\\ (\\text{film}),\\qquad t\\propto R^2\\ (\\text{ash}),\\qquad t\\propto R\\ (\\text{reaction})', 'Eqs. 25.36–25.38']
    ],
    method: ['Table 25.1 is implemented for every shape and inverted exactly (closed forms, or Brent’s method for ash diffusion in cylinders and spheres). With several resistances in a particle of constant size the times add exactly (Eq. 25.32), so t(X) is a sum and X(t) is found by Brent. A shrinking particle is integrated as dR/dt = −(bC_A/ρ_B)/(1/k_g + 1/k″) with k_g re-evaluated from Froessling’s equation at every size.', 'Runs: for each mechanism every run gives a τ (or, for a run that reached completion, an upper bound on τ). Values from one size must agree with one another, and across sizes τ must scale as R, R², or R^{1.5–2} (film). A mechanism is called consistent when the spread is below 6% and no completed run contradicts it.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 's', ts = { s: 1, min: 60, h: 3600, day: 86400 }[tu] || 1;
      const shp = v.beh === 'shrink' ? 'sph' : v.shape;
      const curveSet = v.beh === 'shrink' ? ['stokes', 'large', 'srxn'] : ['film', 'ash', 'rxn'];
      const fig = () => {
        const ss = N.linspace(0, 1, 160);
        const cols = { film: 2, ash: 3, rxn: 1, stokes: 2, large: 4, srxn: 1 };
        plots.push({ title: 'Conversion against time for each controlling step (Fig. 25.10)', x: { label: 't/τ', min: 0, max: 1 }, y: { label: '1 − X_B', min: 0, max: 1 }, series: curveSet.map(m => ({ name: FP.MECH[m], x: ss, y: ss.map(s => 1 - FP.Xof(m, shp, s)), color: cols[m], w: 1.8 })), cap: PEL_NAME(shp) + (v.beh === 'shrink' ? ', shrinking.' : ', constant size.') });
      };
      const PEL_NAME = s => (s === 'slab' ? 'Flat plate' : s === 'cyl' ? 'Cylinder' : 'Sphere');
      if (v.mode === 'fig') {
        const Tg = N.linspace(1000, 2200, 120), E = +v.Ek * 1000, A = +v.Ak;
        const kr = T => A * Math.sqrt(T) * Math.exp(-E / (Rg * T));
        const ser = [{ name: 'surface reaction k″ (Eq. 25.39 form)', x: Tg, y: Tg.map(kr), color: 1, w: 2.2 }];
        const ds = [0.1, 1, 10], us = [1, 10, 100, 1000];
        ds.forEach((dmm, i) => us.forEach((ucm, j) => {
          const kg = Tg.map(T => { const a = FP.air(T); return 100 * FP.kgF(dmm / 1000, ucm / 100, a.rho, a.mu, a.D); });
          ser.push({ name: j === 0 ? 'film k_g, d_p = ' + dmm + ' mm' : undefined, x: Tg, y: kg, color: i + 2, dash: '4 3', w: 1, legend: j === 0 });
          ser.push({ name: j === 0 ? 'actual k̄, d_p = ' + dmm + ' mm' : undefined, x: Tg, y: Tg.map((T, q) => 1 / (1 / kr(T) + 1 / kg[q])), color: i + 2, w: 1.8, legend: j === 0 });
        }));
        plots.push({ title: 'Combustion of pure carbon particles (after Fig. 25.12)', x: { label: 'T', unit: 'K', min: 1000, max: 2200 }, y: { label: 'k″, k_g, k̄', unit: 'cm/s', log: true, min: 1, max: 1000 }, series: ser, hmax: 520, cap: 'For each size the four lines are u = 1, 10, 100 and 1000 cm/s. At low T the steep reaction line controls; at high T the curves flatten onto the film lines, which barely depend on temperature but strongly on size and velocity.' });
        const T0 = 1400, a0 = FP.air(T0);
        kpi.push({ k: 'kr1400', l: 'k″ at 1400 K', v: kr(T0), u: 'cm/s' });
        kpi.push({ k: 'kg1400', l: 'k_g at 1400 K, 1 mm, 10 cm/s', v: 100 * FP.kgF(1e-3, 0.1, a0.rho, a0.mu, a0.D), u: 'cm/s', s: 'Sc = ' + fmt(a0.mu / (a0.rho * a0.D), 3) });
        info.push('Air properties from Sutherland’s law and D_O₂ = 0.181 (T/273)^{1.75} cm²/s; the film follows Froessling. The preexponential factor is set so the reaction line passes through the figure; enter your own k″(T) to replace it.');
        return { kpi, warn, info, plots, tables };
      }
      if (v.mode === 'data') {
        const rows = (v.runs || []).filter(r => r && r.slice(0, 3).every(z => z !== '' && isFinite(+z))).map(r => ({ R: +r[0], t: +r[1], m: +r[2] }));
        if (!rows.length) throw new Error('Enter at least one run.');
        rows.forEach(r => {
          if (v.meas === 'r') r.X = v.beh === 'shrink' ? 1 - Math.pow(r.m, 3) : shp === 'slab' ? 1 - r.m : shp === 'cyl' ? 1 - r.m * r.m : 1 - Math.pow(r.m, 3);
          else r.X = r.m;
          r.X = Math.min(Math.max(r.X, 0), 1);
        });
        const verdicts = [], tabRows = [];
        const eq = v.doneEq !== false;
        curveSet.forEach(mech => {
          const [p0, p1] = FP.EXPO[mech];
          const est = rows.filter(r => r.X > 0 && (r.X < 1 - 1e-9 || eq)).map(r => ({ R: r.R, tau: r.X >= 1 - 1e-9 ? r.t : r.t / FP.g(mech, shp, r.X) }));
          const bound = eq ? [] : rows.filter(r => r.X >= 1 - 1e-9);
          const sizes = [...new Set(rows.map(r => r.R))];
          let best = null;
          const ps = p0 === p1 ? [p0] : N.linspace(p0, p1, 11);
          ps.forEach(p => {
            let A = NaN, spread = 0, sig = 0, viol = 0;
            const src = est.length ? est : bound.map(r => ({ R: r.R, tau: r.t }));
            if (src.length) {
              const la = src.map(e => Math.log(e.tau) - p * Math.log(e.R)), m = N.mean(la);
              A = Math.exp(m);
              spread = src.length > 1 ? Math.max(...la) - Math.min(...la) : 0;
              sig = Math.sqrt(N.mean(la.map(q => (q - m) * (q - m))));
            }
            if (est.length) bound.forEach(r => { const tau = A * Math.pow(r.R, p); if (tau > r.t * 1.06) viol = Math.max(viol, Math.log(tau / r.t)); });
            const score = sig + viol;
            if (!best || score < best.score) best = { p, A, spread, sig, viol, score };
          });
          verdicts.push({ mech, ok: best.spread + best.viol < Math.log(1.06), best, sizes: sizes.length, n: est.length });
          rows.forEach(r => tabRows.push([FP.MECH[mech], r.R, r.t, r.X, r.X >= 1 - 1e-9 ? (eq ? r.t : '≤ ' + N.fmt(r.t, 4)) : r.X > 0 ? r.t / FP.g(mech, shp, r.X) : '—', best.A * Math.pow(r.R, best.p)]));
        });
        let okM = verdicts.filter(q => q.ok), how = 'within 6%';
        if (!okM.length && verdicts[0].n >= 3) {
          const srt = verdicts.slice().sort((x, y) => x.best.score - y.best.score);
          if (srt[0].best.sig < 0.15 && srt[0].best.score < 0.6 * srt[1].best.score) { okM = [srt[0]]; how = 'best of the three, scatter ' + fmt((Math.exp(srt[0].best.sig) - 1) * 100, 2) + '%'; }
        }
        if (v.beh === 'const' && okM.length > 1 && okM.some(q => q.mech === 'film') && okM.some(q => q.mech === 'ash')) {
          okM = okM.filter(q => q.mech !== 'film');
          info.push('Film and ash both fit, but with a firm, nonflaking product the ash layer’s resistance is usually far larger than the film’s, so film control is set aside (Ch. 25.5, “ash versus film resistance”).');
        }
        const ambiguous = okM.length > 1;
        kpi.push({ k: 'verdict', l: 'Controlling step', v: okM.length === 1 ? FP.MECH[okM[0].mech] : okM.length ? 'not decided' : 'none fits', main: true, s: okM.length > 1 ? 'consistent with ' + okM.map(q => FP.MECH[q.mech]).join(' and ') : okM.length ? how : 'the data fit no single controlling step' });
        verdicts.forEach(q => kpi.push({ k: 'tau_' + q.mech, l: 'τ if ' + FP.MECH[q.mech] + ' controls', v: q.best.A * Math.pow(rows[0].R, q.best.p), u: tu, s: (okM.indexOf(q) >= 0 ? 'consistent' : 'inconsistent') + (q.sizes > 1 ? '; τ ∝ R^' + fmt(q.best.p, 3) : '') + (q.best.spread > 0 ? '; spread ' + fmt((Math.exp(q.best.spread) - 1) * 100, 3) + '%' : ''), state: okM.indexOf(q) >= 0 ? 'ok' : 'bad' }));
        if (ambiguous) info.push('More than one controlling step is consistent with these runs. Temperature is the usual discriminator: the chemical step is far more temperature-sensitive than the physical ones (Ch. 25.5); so are runs at other sizes, or a look at whether a firm ash layer forms.');
        if (v.beh === 'const' && okM.some(q => q.mech === 'film') && okM.some(q => q.mech === 'ash')) info.push('When a hard, nonflaking ash forms, its resistance is usually far larger than the gas film’s, so film control can be ruled out (Ch. 25.5).');
        tables.push({ title: 'Each run read by each mechanism (Table 25.1)', cols: ['mechanism', 'size', 't', 'X_B', 'τ from this run', 'τ fitted'], rows: tabRows });
        /* runs on the curves */
        fig();
        const best = okM.length ? okM[0] : verdicts.reduce((a, b) => (b.best.score < a.best.score ? b : a));
        plots[0].series.push({ name: 'runs (t/τ of ' + FP.MECH[best.mech] + ')', x: rows.map(r => r.t / (best.best.A * Math.pow(r.R, best.best.p))), y: rows.map(r => 1 - r.X), type: 'pts', color: 8 });
        return { kpi, warn, info, plots, tables };
      }
      /* ---------- one particle from its properties ---------- */
      const R = +v.R / 1000, rhoB = +v.rhoB, b = +v.b;
      const CA = v.cgiv === 'C' ? +v.CA : +v.yA * +v.P * 101325 / (Rg * (+v.TC + 273.15));
      const F = FP.F[shp];
      const tq = +v.tq * ts;
      if (v.beh === 'const') {
        let kg = NaN;
        if (v.uf) kg = v.kgs === 'given' ? +v.kg / 100 : FP.kgF(2 * R, +v.u / 100, 1, +v.nu * 1e-4, +v.Dg * 1e-4);
        const taus = {
          film: v.uf ? rhoB * R / (F * b * kg * CA) : 0,
          ash: v.ua ? rhoB * R * R / (2 * F * b * (+v.De * 1e-4) * CA) : 0,
          rxn: v.ur ? rhoB * R / (b * (+v.kr / 100) * CA) : 0
        };
        if (!(taus.film || taus.ash || taus.rxn)) throw new Error('Switch on at least one resistance.');
        const C = FP.combo(taus, shp);
        kpi.push({ k: 'tau', l: 'Time for complete conversion τ', v: C.tau / ts, u: tu, main: true, s: fmt(C.tau / 3600, 4) + ' h' });
        ['film', 'ash', 'rxn'].forEach(m => { if (taus[m]) kpi.push({ k: 'tau_' + m, l: 'τ with ' + FP.MECH[m] + ' alone', v: taus[m] / ts, u: tu, s: fmt(taus[m] / C.tau * 100, 3) + '% of the total' }); });
        kpi.push({ k: 'Xq', l: 'X_B after ' + fmt(+v.tq, 4) + ' ' + tu, v: C.Xt(tq) });
        kpi.push({ k: 'CA', l: 'C_Ag', v: CA, u: 'mol/m³' });
        if (v.uf) kpi.push({ k: 'kg', l: 'k_g', v: kg * 100, u: 'cm/s' });
        /* resistance shares against conversion, Eq. 25.33b (sphere) generalized */
        const Xs = N.linspace(0, 0.999, 120);
        const share = X => {
          const tx = m => (taus[m] ? taus[m] * (FP.g(m, shp, Math.min(X + 1e-4, 1)) - FP.g(m, shp, Math.max(X - 1e-4, 0))) : 0);
          const a = tx('film'), bb = tx('ash'), c = tx('rxn'), s = a + bb + c;
          return [a / s, bb / s, c / s];
        };
        const sh = Xs.map(share);
        plots.push({ title: 'Which step holds the particle back as the core recedes', x: { label: 'X_B', min: 0, max: 1 }, y: { label: 'share of the instantaneous resistance', min: 0, max: 1 }, series: [{ name: 'gas film', x: Xs, y: sh.map(q => q[0]), color: 2 }, { name: 'ash layer', x: Xs, y: sh.map(q => q[1]), color: 3 }, { name: 'reaction', x: Xs, y: sh.map(q => q[2]), color: 1 }], cap: 'The ash builds up as the core shrinks, while the reacting surface shrinks: their roles change with conversion (Ch. 25.4).' });
        fig();
        const tg = N.linspace(0, C.tau, 140);
        plots[1].series.push({ name: 'this particle, all resistances', x: tg.map(t => t / C.tau), y: tg.map(t => 1 - C.Xt(t)), color: 8, w: 2.6 });
        tables.push({ title: 'Core radius and conversion against time', cols: ['t, ' + tu, 't/τ', 'X_B', 'r_c/R'], rows: [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1].map(f => { const X = C.Xt(f * C.tau); return [f * C.tau / ts, f, X, shp === 'slab' ? 1 - X : shp === 'cyl' ? Math.sqrt(1 - X) : Math.cbrt(1 - X)]; }) });
        return { kpi, warn, info, plots, tables };
      }
      /* shrinking sphere: film (Froessling or given, updated with size) and reaction in series */
      const kgOf = Rn => (!v.uf ? Infinity : v.kgs === 'given' ? +v.kg / 100 * (R / Rn) : FP.kgF(2 * Rn, +v.u / 100, 1, +v.nu * 1e-4, +v.Dg * 1e-4));
      const kr = v.ur ? +v.kr / 100 : Infinity;
      if (!v.uf && !v.ur) throw new Error('Switch on the film, the reaction, or both.');
      const dRdt = Rn => -(b * CA / rhoB) / (1 / kgOf(Math.max(Rn, 1e-12)) + 1 / kr);
      const tau = N.quad(Rn => -1 / dRdt(Rn), 0, R, { rtol: 1e-9 });
      kpi.push({ k: 'tau', l: 'Time for the particle to disappear τ', v: tau / ts, u: tu, main: true, s: fmt(tau / 3600, 4) + ' h' });
      const tR = Rn => N.quad(q => -1 / dRdt(q), Rn, R, { rtol: 1e-9 });
      const Rt = t => (t >= tau ? 0 : t <= 0 ? R : N.brent(Rn => tR(Rn) - t, 0, R, { tol: 1e-12 * R }));
      const Rq = Rt(tq);
      kpi.push({ k: 'Xq', l: 'X_B after ' + fmt(+v.tq, 4) + ' ' + tu, v: 1 - Math.pow(Rq / R, 3), s: 'R = ' + fmt(Rq * 1000, 4) + ' mm' });
      if (v.ur) kpi.push({ k: 'tau_r', l: 'τ with reaction alone', v: rhoB * R / (b * kr * CA) / ts, u: tu });
      if (v.uf && v.kgs === 'fro') {
        const DgS = +v.Dg * 1e-4;
        kpi.push({ k: 'tau_s', l: 'τ with the film alone, Stokes regime (Eq. 25.29)', v: rhoB * R * R / (2 * b * CA * DgS) / ts, u: tu, s: 'a lower bound on film control' });
        kpi.push({ k: 'Re0', l: 'Re at the start', v: 2 * R * (+v.u / 100) / (+v.nu * 1e-4) });
      }
      const tg = N.linspace(0, tau, 160), Rs = tg.map(Rt);
      const ser = [{ name: 'R/R₀', x: tg.map(t => t / ts), y: Rs.map(q => q / R), color: 1, w: 2.2 }, { name: '1 − X_B', x: tg.map(t => t / ts), y: Rs.map(q => Math.pow(q / R, 3)), color: 2 }];
      if (v.uf && v.ur) ser.push({ name: 'film share of resistance', x: tg.map(t => t / ts), y: Rs.map(q => { const kg = kgOf(Math.max(q, 1e-12)); return (1 / kg) / (1 / kg + 1 / kr); }), color: 3, dash: '5 4' });
      plots.push({ title: 'The shrinking particle', x: { label: 't', unit: tu, min: 0 }, y: { label: 'fraction', min: 0, max: 1 }, series: ser });
      fig();
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
