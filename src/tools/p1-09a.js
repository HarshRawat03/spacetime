/* Chapter 9 — Temperature and pressure effects: thermochemistry and the X–T chart */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, TH = ST.TH, D = ST.D;
  const C2K = TH.C2K, K2C = TH.K2C, R = TH.R;

  /* ---------- shared input block for one reaction on the X–T plane (H-902 … H-904, Ch 19) ---------- */
  const FORMS = {
    'A=R': [1, 0, 1, 0], 'A=2R': [1, 0, 2, 0], 'A=R+S': [1, 0, 1, 1], 'A+B=R': [1, 1, 1, 0],
    'A+B=R+S': [1, 1, 1, 1], 'A+B=2R': [1, 1, 2, 0], '2A=R+S': [2, 0, 1, 1], '2A=R': [2, 0, 1, 0]
  };
  const FORM_OPTS = [['A=R', 'A ⇌ R'], ['A=2R', 'A ⇌ 2R'], ['A=R+S', 'A ⇌ R + S'], ['A+B=R', 'A + B ⇌ R'], ['A+B=R+S', 'A + B ⇌ R + S'], ['A+B=2R', 'A + B ⇌ 2R'], ['2A=R+S', '2A ⇌ R + S'], ['2A=R', '2A ⇌ R']];
  const XTI = ST.XTI = { FORMS, FORM_OPTS };
  const fwdOrder = v => { const f = FORMS[v.form]; return v.rev ? f[0] + f[1] : (+v.oA || 0) + (f[1] ? (+v.oB || 0) : 0); };
  XTI.kUnit = v => { const n = fwdOrder(v), t = v._tu || 'min'; if (Math.abs(n - 1) < 1e-9) return '1/' + t; const e = n - 1; return '(L/mol)' + (Math.abs(e - 1) < 1e-9 ? '' : '^{' + N.fmt(e, 3) + '}') + '/' + t; };
  const dnOf = v => { const f = FORMS[v.form]; return f[2] + f[3] - f[0] - f[1]; };
  const kEqUnit = v => { const dn = dnOf(v); if (!dn) return ''; return (v.phase === 'gas' ? 'atm' : 'mol/L') + (dn === 1 ? '' : '^{' + dn + '}'); };

  XTI.inputs = function (d = {}) {
    const has = i => v => FORMS[v.form][i] > 0;
    return [
      { group: 'Reaction' },
      { k: 'form', t: 'select', label: 'Stoichiometry (A is the key reactant)', opts: FORM_OPTS, v: d.form || 'A=R' },
      { k: 'rev', t: 'check', label: 'Reversible, elementary both ways', v: d.rev !== undefined ? d.rev : true },
      { k: 'oA', label: 'Order in A', v: d.oA !== undefined ? d.oA : 1, min: 0, max: 4, step: 0.05, if: v => !v.rev },
      { k: 'oB', label: 'Order in B', v: d.oB !== undefined ? d.oB : 1, min: 0, max: 4, step: 0.05, if: v => !v.rev && FORMS[v.form][1] > 0 },
      { k: 'phase', t: 'seg', label: 'Phase', opts: [['liq', 'liquid, constant density'], ['gas', 'ideal gas, constant P']], v: d.phase || 'liq' },
      { group: 'Feed' },
      { k: 'CA0', label: '$C_{A0}$', unit: 'mol/L', v: d.CA0 || 1, min: 1e-6, max: 1e3, scale: 'log', if: v => v.phase === 'liq' },
      { k: 'P', label: 'Pressure', unit: 'atm', v: d.P || 1, min: 1e-3, max: 1e3, scale: 'log', if: v => v.phase === 'gas' },
      { k: 'yA0', label: 'Mole fraction of A in the feed', v: d.yA0 || 1, min: 1e-4, max: 1, if: v => v.phase === 'gas' },
      { k: 'MB', label: 'B fed per mol A', v: d.MB !== undefined ? d.MB : 1, min: 0, max: 100, if: has(1) },
      { k: 'MR', label: 'R fed per mol A', v: d.MR || 0, min: 0, max: 100 },
      { k: 'MS', label: 'S fed per mol A', v: d.MS || 0, min: 0, max: 100, if: has(3) },
      { group: 'Forward rate constant, k₁ = A₁ e^{−E₁/RT}' },
      { k: 'ksrc', t: 'seg', label: 'Given as', opts: [['arr', 'A₁ and E₁'], ['ref', 'k₁ at one T'], ['runs', 'two batch runs']], v: d.ksrc || 'arr' },
      { k: 'A1', label: '$A_1$', unit: XTI.kUnit, v: d.A1 || 3.39e7, min: 1e-12, max: 1e40, scale: 'log', if: v => v.ksrc === 'arr' },
      { k: 'k1r', label: '$k_1$ at $T_k$', unit: XTI.kUnit, v: d.k1r || 0.1, min: 1e-12, max: 1e12, scale: 'log', if: v => v.ksrc === 'ref' },
      { k: 'Tk', label: '$T_k$', unit: '°C', v: d.Tk !== undefined ? d.Tk : 25, min: -100, max: 1500, if: v => v.ksrc === 'ref' },
      { k: 'E1', label: '$E_1$', unit: 'J/mol', v: d.E1 !== undefined ? d.E1 : 48900, min: 0, max: 400000, if: v => v.ksrc !== 'runs' },
      { k: 'r1T', label: 'Run 1: temperature', unit: '°C', v: d.r1T !== undefined ? d.r1T : 65, min: -50, max: 500, if: v => v.ksrc === 'runs' },
      { k: 'r1t', label: 'Run 1: time', unit: '[t]', v: d.r1t || 1, min: 1e-6, max: 1e6, scale: 'log', if: v => v.ksrc === 'runs' },
      { k: 'r1X', label: 'Run 1: conversion reached', v: d.r1X || 0.581, min: 1e-6, max: 0.999999, if: v => v.ksrc === 'runs' },
      { k: 'r2T', label: 'Run 2: temperature', unit: '°C', v: d.r2T !== undefined ? d.r2T : 25, min: -50, max: 500, if: v => v.ksrc === 'runs' },
      { k: 'r2t', label: 'Run 2: time', unit: '[t]', v: d.r2t || 10, min: 1e-6, max: 1e6, scale: 'log', if: v => v.ksrc === 'runs' },
      { k: 'r2X', label: 'Run 2: conversion reached', v: d.r2X || 0.6, min: 1e-6, max: 0.999999, if: v => v.ksrc === 'runs' },
      { group: 'Equilibrium and heat of reaction' },
      { k: 'Ksrc', t: 'seg', label: 'K from', opts: [['thermo', 'ΔG° and ΔH° at 25 °C'], ['KT', 'K at one T, and ΔH°'], ['arr2', 'reverse A₂, E₂']], v: d.Ksrc || 'thermo', if: v => v.rev },
      { k: 'dG', label: 'ΔG° at 25 °C', unit: 'J/mol', v: d.dG !== undefined ? d.dG : -14130, min: -5e5, max: 5e5, if: v => v.rev && v.Ksrc === 'thermo', hint: 'Per mole of reaction as written. Standard states: 1 mol/L for liquids, 1 atm for gases.' },
      { k: 'KT', label: 'K at $T_K$', unit: kEqUnit, v: d.KT || 10, min: 1e-30, max: 1e30, scale: 'log', if: v => v.rev && v.Ksrc === 'KT' },
      { k: 'TK', label: '$T_K$', unit: '°C', v: d.TK !== undefined ? d.TK : 25, min: -100, max: 2000, if: v => v.rev && v.Ksrc === 'KT' },
      { k: 'A2', label: '$A_2$ (reverse)', v: d.A2 || 1.8e18, min: 1e-12, max: 1e40, scale: 'log', if: v => v.rev && v.Ksrc === 'arr2' },
      { k: 'E2', label: '$E_2$ (reverse)', unit: 'J/mol', v: d.E2 !== undefined ? d.E2 : 124200, min: 0, max: 600000, if: v => v.rev && v.Ksrc === 'arr2' },
      { k: 'dH', label: 'ΔH_r at 25 °C', unit: 'J/mol', v: d.dH !== undefined ? d.dH : -75300, min: -1e6, max: 1e6, if: v => !(v.rev && v.Ksrc === 'arr2'), hint: 'Per mole of reaction as written; negative for exothermic.' },
      { k: 'dCp', label: 'ΔC_p of reaction', unit: 'J/(mol·K)', v: d.dCp || 0, min: -1000, max: 1000, if: v => !(v.rev && v.Ksrc === 'arr2') },
      { group: 'Energy and temperature limits' },
      { k: 'Cp', label: "$C'_p$ of the feed, per mol A", unit: 'J/(mol A·K)', v: d.Cp || 1046, min: 1, max: 1e6, scale: 'log', hint: 'Includes solvent and inerts. Product stream: C″_p = C′_p + ΔC_p/a.' },
      { k: 'Tlo', label: 'Lowest allowable T', unit: '°C', v: d.Tlo !== undefined ? d.Tlo : 5, min: -100, max: 1500 },
      { k: 'Thi', label: 'Highest allowable T', unit: '°C', v: d.Thi !== undefined ? d.Thi : 95, min: -100, max: 2000 }
    ];
  };

  /* forward constant from two batch runs (constant density) with K(T) known */
  function fitRuns(p, v) {
    if (p.phase === 'gas') throw new Error('Batch runs are read as constant-density runs; switch to the liquid phase or give A₁ and E₁.');
    if (p.rev && v.Ksrc === 'arr2') throw new Error('Batch runs need K from thermodynamics (ΔG°, or K at one T), not from the reverse rate constant.');
    const unit = TH.model(Object.assign({}, p, { A1: 1, E1: 0 }));
    const one = (Tc, t, X) => {
      const T = C2K(Tc), xe = unit.Xe(T);
      if (!(X < xe)) throw new Error('A run at ' + N.fmt(Tc, 4) + ' °C reports X = ' + N.fmt(X, 4) + ', at or beyond equilibrium (' + N.fmt(xe, 4) + ').');
      const I = N.quad(x => 1 / unit.rate(x, T), 0, X, { rtol: 1e-10 });
      return unit.CA0f(T) * I / t;
    };
    const ka = one(+v.r1T, +v.r1t, +v.r1X), kb = one(+v.r2T, +v.r2t, +v.r2X);
    const Ta = C2K(+v.r1T), Tb = C2K(+v.r2T);
    if (Math.abs(Ta - Tb) < 1e-6) throw new Error('The two runs must be at different temperatures.');
    const E1 = R * Math.log(ka / kb) / (1 / Tb - 1 / Ta);
    return { ka, kb, Ta, Tb, E1, A1: ka * Math.exp(E1 / (R * Ta)) };
  }

  XTI.model = function (v) {
    const [a, b, r, s] = FORMS[v.form];
    const p = {
      a, b, r, s, rev: !!v.rev, oA: +v.oA, oB: +v.oB, phase: v.phase, CA0: +v.CA0, P: +v.P, yA0: +v.yA0,
      MB: b ? +v.MB : 0, MR: +v.MR || 0, MS: s ? +v.MS || 0 : 0, Cp: +v.Cp, Tlo: C2K(+v.Tlo), Thi: C2K(+v.Thi), tsec: TH.TSEC[v._tu || 'min'] || 60
    };
    if (!(p.Thi > p.Tlo)) throw new Error('The highest allowable temperature must be above the lowest.');
    if (p.phase === 'gas' && p.yA0 * (1 + p.MB + p.MR + p.MS) > 1 + 1e-9) throw new Error('Feed mole fractions add to more than 1: lower y_A0 or the feed ratios.');
    if (b && !(p.MB > 0)) throw new Error('Feed some B: the reaction needs it.');
    if (p.rev && !(r + s > 0)) throw new Error('A reversible reaction needs products.');
    const T0 = 298.15;
    let c = TH.coef({ T0, H0: +v.dH, da: +v.dCp, G0: +v.dG });
    if (p.rev && v.Ksrc === 'KT') {
      const c0 = TH.coef({ T0, H0: +v.dH, da: +v.dCp, lnK0: 0 });
      c = TH.coef({ T0, H0: +v.dH, da: +v.dCp, lnK0: Math.log(+v.KT) - TH.lnK(c0, C2K(+v.TK)) });
    }
    p.c = c; p.Ktype = p.rev ? v.Ksrc : 'thermo'; p.A2 = +v.A2; p.E2 = +v.E2;
    let fit = null;
    if (v.ksrc === 'ref') { p.E1 = +v.E1; p.A1 = +v.k1r * Math.exp(p.E1 / (R * C2K(+v.Tk))); }
    else if (v.ksrc === 'runs') { fit = fitRuns(p, v); p.A1 = fit.A1; p.E1 = fit.E1; }
    else { p.A1 = +v.A1; p.E1 = +v.E1; }
    const m = TH.model(p);
    m.fit = fit;
    m.label = FORM_OPTS.find(o => o[0] === v.form)[1].replace('⇌', p.rev ? '⇌' : '→');
    return m;
  };

  /* chart + overlays, shared with the design tools */
  XTI.chartSpec = function (m, o = {}) {
    const ch = TH.chart(m, o);
    return {
      title: o.title || 'Conversion–temperature chart, contours of −r_A (Fig. 9.3)',
      x: ch.x, y: ch.y, series: ch.series.concat(o.extra || []), marks: ch.marks.concat(o.marks || []), aspect: o.aspect || 0.5, hmax: 560, wide: true,
      cap: o.cap || 'Contours are rates of reaction of A, ' + (o.ru || 'mol/(L·min)') + '; labels sit on the locus of maximum rates.'
    };
  };

  /* ===================================================================== */
  /* H-901  heats of reaction, K(T) and equilibrium conversion              */
  /* ===================================================================== */
  const SP_WGS = '# water-gas shift; Cp/R = a + bT + cT² + d/T² (Smith, Van Ness & Abbott form)\n# name   nu   n0    dHf kJ/mol  dGf kJ/mol   a      b         c   d\nCO      -1    1    -110.53     -137.17     3.376  0.557e-3  0   -0.031e5\nH2O     -1    1    -241.818    -228.572    3.470  1.450e-3  0    0.121e5\nCO2      1    0    -393.509    -394.359    5.457  1.045e-3  0   -1.157e5\nH2       1    0       0           0        3.249  0.422e-3  0    0.083e5';
  const SP_NH3 = '# ammonia synthesis, ideal gas\n# name   nu    n0   dHf      dGf      a      b         c   d\nN2      -0.5  1    0        0        3.280  0.593e-3  0    0.040e5\nH2      -1.5  3    0        0        3.249  0.422e-3  0    0.083e5\nNH3      1    0   -46.11   -16.45    3.578  3.020e-3  0   -0.186e5';
  const SP_SO2 = '# SO2 oxidation in air (8% SO2, 13% O2, rest N2)\n# name   nu    n0    dHf       dGf       a      b         c   d\nSO2     -1    8    -296.83   -300.194   5.699  0.801e-3  0   -1.015e5\nO2      -0.5  13     0         0        3.639  0.506e-3  0   -0.227e5\nSO3      1    0    -395.72   -371.06    8.060  1.056e-3  0   -2.028e5\nN2       0    79     0         0        3.280  0.593e-3  0    0.040e5';

  function parseSp(text, cpR) {
    const rows = [];
    String(text).split('\n').forEach((raw, i) => {
      const line = raw.replace(/#.*$/, '').trim();
      if (!line) return;
      const t = line.split(/[\s,;]+/);
      if (t.length < 6) throw new Error('Line ' + (i + 1) + ': expected name, ν, n₀, ΔH_f, ΔG_f, then Cp coefficients.');
      const num = t.slice(1).map(x => parseFloat(x));
      if (num.slice(0, 5).some(x => !isFinite(x))) throw new Error('Line ' + (i + 1) + ': cannot read the numbers.');
      const f = cpR ? R : 1;
      rows.push({ name: t[0], nu: num[0], n0: num[1], dHf: num[2] * 1000, dGf: num[3] * 1000, a: num[4] * f, b: (num[5] || 0) * f, c: (num[6] || 0) * f, d: (num[7] || 0) * f });
    });
    if (rows.length < 2) throw new Error('List at least two species.');
    if (!rows.some(r => r.nu < 0) || !rows.some(r => r.nu > 0)) throw new Error('Give reactants (ν < 0) and products (ν > 0).');
    return rows;
  }

  ST.tool({
    id: 'thermo', tag: 'H-901', ch: 9, sec: 'Ch. 9.1', title: 'Heat of reaction, equilibrium constant and equilibrium conversion against temperature',
    lede: 'ΔH_r(T) from heat capacities, K(T) from the van ’t Hoff equation integrated with a temperature-dependent ΔC_p, and the equilibrium composition for gases or solutions, with pressure and inerts. Either enter the reaction totals or a table of species with formation data.',
    keys: 'heat of reaction temperature kirchhoff van t hoff equilibrium constant gibbs free energy equilibrium conversion pressure inerts le chatelier adiabatic equilibrium temperature heat capacity polynomial formation enthalpy',
    inputs: [
      { group: 'Data' },
      { k: 'mode', t: 'seg', label: 'Enter', opts: [['rxn', 'reaction totals'], ['sp', 'species table']], v: 'rxn' },
      { k: 'form', t: 'select', label: 'Stoichiometry', opts: FORM_OPTS, v: 'A=R', if: v => v.mode === 'rxn' },
      { k: 'dH', label: 'ΔH°_r at 25 °C', unit: 'J/mol', v: -75300, min: -1e7, max: 1e7, if: v => v.mode === 'rxn' },
      { k: 'Gsrc', t: 'seg', label: 'Equilibrium given by', opts: [['dG', 'ΔG° at 25 °C'], ['K', 'K at 25 °C']], v: 'dG', if: v => v.mode === 'rxn' },
      { k: 'dG', label: 'ΔG° at 25 °C', unit: 'J/mol', v: -14130, min: -1e7, max: 1e7, if: v => v.mode === 'rxn' && v.Gsrc === 'dG' },
      { k: 'K298', label: 'K at 25 °C', v: 300, min: 1e-40, max: 1e40, scale: 'log', if: v => v.mode === 'rxn' && v.Gsrc === 'K' },
      { k: 'cpmode', t: 'seg', label: 'Heat capacities', opts: [['none', 'ΔC_p = 0'], ['mean', 'mean C_p of each species'], ['delta', 'Δα, Δβ, Δγ']], v: 'none', if: v => v.mode === 'rxn' },
      { k: 'CpA', label: 'C_p of A', unit: 'J/(mol·K)', v: 35, min: 0, max: 1000, if: v => v.mode === 'rxn' && v.cpmode === 'mean' },
      { k: 'CpB', label: 'C_p of B', unit: 'J/(mol·K)', v: 45, min: 0, max: 1000, if: v => v.mode === 'rxn' && v.cpmode === 'mean' && FORMS[v.form][1] > 0 },
      { k: 'CpR', label: 'C_p of R', unit: 'J/(mol·K)', v: 70, min: 0, max: 1000, if: v => v.mode === 'rxn' && v.cpmode === 'mean' },
      { k: 'CpS', label: 'C_p of S', unit: 'J/(mol·K)', v: 40, min: 0, max: 1000, if: v => v.mode === 'rxn' && v.cpmode === 'mean' && FORMS[v.form][3] > 0 },
      { k: 'da', label: 'Δα', unit: 'J/(mol·K)', v: 0, min: -1000, max: 1000, if: v => v.mode === 'rxn' && v.cpmode === 'delta' },
      { k: 'db', label: 'Δβ', unit: 'J/(mol·K²)', v: 0, min: -1, max: 1, if: v => v.mode === 'rxn' && v.cpmode === 'delta' },
      { k: 'dc', label: 'Δγ', unit: 'J/(mol·K³)', v: 0, min: -1e-3, max: 1e-3, if: v => v.mode === 'rxn' && v.cpmode === 'delta' },
      { k: 'sp', t: 'area', label: 'Species: name, ν, feed moles, ΔH°_f and ΔG°_f (kJ/mol, 25 °C), C_p = a + bT + cT² + d/T²', v: SP_WGS, rows: 7, if: v => v.mode === 'sp' },
      { k: 'cpR', t: 'check', label: 'C_p coefficients are given as C_p/R', v: true, if: v => v.mode === 'sp' },
      { group: 'Equilibrium mixture' },
      { k: 'std', t: 'seg', label: 'Standard state', opts: [['gas', 'ideal gas'], ['conc', '1 mol/L solute'], ['x', 'ideal solution (x)']], v: 'conc' },
      { k: 'Pstd', t: 'seg', label: 'Standard pressure', opts: [['1', '1 atm'], ['0.986923', '1 bar']], v: '1', if: v => v.std === 'gas' },
      { k: 'P', label: 'Pressure', unit: 'atm', v: 1, min: 1e-4, max: 1e4, scale: 'log', if: v => v.std === 'gas' },
      { k: 'CA0', label: 'Concentration scale: $C_{A0}$', unit: 'mol/L', v: 1, min: 1e-8, max: 1e3, scale: 'log', if: v => v.std === 'conc' && v.mode === 'rxn', hint: 'With 1 mol/L standard states K carries units unless Δn = 0.' },
      { k: 'MB', label: 'B fed per mol A', v: 1, min: 0, max: 1000, if: v => v.mode === 'rxn' && FORMS[v.form][1] > 0 },
      { k: 'MR', label: 'R fed per mol A', v: 0, min: 0, max: 1000, if: v => v.mode === 'rxn' },
      { k: 'MS', label: 'S fed per mol A', v: 0, min: 0, max: 1000, if: v => v.mode === 'rxn' && FORMS[v.form][3] > 0 },
      { k: 'MI', label: 'Inerts fed per mol A', v: 0, min: 0, max: 1000, if: v => v.mode === 'rxn' && v.std !== 'conc' },
      { group: 'Temperatures' },
      { k: 'T2', label: 'Evaluate at', unit: '°C', v: 65, min: -200, max: 3000 },
      { k: 'Ta', label: 'Plot from', unit: '°C', v: 0, min: -200, max: 3000 },
      { k: 'Tb', label: 'Plot to', unit: '°C', v: 100, min: -200, max: 3500 },
      { k: 'Xt', label: 'Equilibrium conversion wanted', v: 0.75, min: 0.001, max: 0.999 },
      { group: 'Adiabatic reactor' },
      { k: 'adOn', t: 'check', label: 'Find the adiabatic equilibrium point', v: false },
      { k: 'Tf', label: 'Feed temperature', unit: '°C', v: 25, min: -200, max: 3000, if: v => v.adOn },
      { k: 'Cpf', label: "Feed $C'_p$ per mol of key reactant", unit: 'J/(mol·K)', v: 1046, min: 1, max: 1e6, scale: 'log', if: v => v.adOn }
    ],
    presets: [
      { name: 'Example 9.1', v: { mode: 'rxn', form: 'A+B=2R', dH: -50000, Gsrc: 'K', K298: 1e8, cpmode: 'mean', CpA: 35, CpB: 45, CpR: 70, std: 'gas', P: 1, MB: 1, T2: 1025, Ta: 25, Tb: 1100, adOn: false }, book: { dH2: 10000 }, note: 'Gas reaction A + B → 2R with ΔH_r = −50 kJ at 25 °C and mean heat capacities 35, 45 and 70 J/mol·K. At 1025 °C the reaction has turned endothermic: +10 kJ/mol.' },
      { name: 'Example 9.2', v: { mode: 'rxn', form: 'A=R', dH: -75300, Gsrc: 'dG', dG: -14130, cpmode: 'none', std: 'conc', CA0: 1, T2: 25, Ta: 0, Tb: 100, Xt: 0.75, adOn: false }, book: { K298: 300, Tlim: { v: 78, tol: 0.01 } }, note: 'Aqueous A ⇌ R with ΔG° = −14 130 and ΔH° = −75 300 J/mol. K = 300 at 25 °C, and X_{Ae} ≥ 0.75 needs T below 78 °C.' },
      { name: 'Example 9.2 fed adiabatically at 25 °C', v: { mode: 'rxn', form: 'A=R', dH: -75300, Gsrc: 'dG', dG: -14130, cpmode: 'none', std: 'conc', CA0: 4, T2: 65, Ta: 0, Tb: 120, Xt: 0.75, adOn: true, Tf: 25, Cpf: 1046 }, note: 'The 4 mol/L solution of Example 9.5 (C′_p = 250 cal per mol A·K): where the adiabatic line meets the equilibrium curve.' },
      { name: 'Water–gas shift', v: { mode: 'sp', sp: SP_WGS, cpR: true, std: 'gas', Pstd: '0.986923', P: 1, T2: 400, Ta: 150, Tb: 1100, Xt: 0.9, adOn: true, Tf: 200, Cpf: 75 }, note: 'CO + H₂O ⇌ CO₂ + H₂ from formation data and C_p polynomials. K falls through 1 near 820 °C; pressure has no effect because Δn = 0.' },
      { name: 'Ammonia synthesis', v: { mode: 'sp', sp: SP_NH3, cpR: true, std: 'gas', Pstd: '0.986923', P: 100, T2: 450, Ta: 200, Tb: 800, Xt: 0.3, adOn: false }, note: '½N₂ + 1.5H₂ ⇌ NH₃ treated as an ideal gas; the real mixture at 100 atm deviates noticeably, so read this as the ideal-gas estimate.' },
      { name: 'SO₂ oxidation in air', v: { mode: 'sp', sp: SP_SO2, cpR: true, std: 'gas', Pstd: '0.986923', P: 1, T2: 450, Ta: 300, Tb: 900, Xt: 0.97, adOn: true, Tf: 420, Cpf: 400 }, note: 'Converter feed with 8% SO₂: the equilibrium curve and the adiabatic first-bed limit.' }
    ],
    eqs: [
      ['\\Delta H_{r2} = \\Delta H_{r1} + \\int_{T_1}^{T_2} \\Delta C_p\\,dT,\\quad \\Delta C_p = \\Delta\\alpha + \\Delta\\beta T + \\Delta\\gamma T^2', 'Eqs. 4–7'],
      ['\\Delta G^\\circ = -RT\\ln K', 'Eq. 9'],
      ['\\frac{d\\ln K}{dT} = \\frac{\\Delta H_r}{RT^2}', 'Eq. 15'],
      ['\\ln\\frac{K_2}{K_1} = -\\frac{\\Delta H_r}{R}\\left(\\frac{1}{T_2}-\\frac{1}{T_1}\\right)', 'Eq. 16', 'Constant ΔH_r.'],
      ['R\\ln K = -\\frac{\\Delta H_0}{T} + \\Delta\\alpha\\ln T + \\frac{\\Delta\\beta}{2}T + \\frac{\\Delta\\gamma}{6}T^2 + M', 'Eq. 19', 'With ΔH_r = ΔH₀ + ΔαT + ΔβT²/2 + ΔγT³/3.'],
      ['K_p = \\prod_i \\left(\\frac{y_i P}{P^\\circ}\\right)^{\\nu_i} = K_y\\left(\\frac{P}{P^\\circ}\\right)^{\\Delta n}', 'Eq. 11–13'],
      ['X_A = \\frac{C\'_p\\,(T-T_1)}{-\\Delta H_r(T)}', 'Eq. 22', 'Adiabatic operating line.']
    ],
    method: ['ΔH_r(T) and ln K(T) are integrated in closed form for ΔC_p = Δa + ΔbT + ΔcT² + Δd/T². The equilibrium extent is found by bracketed root-finding in a logit coordinate between the extent bounds, which keeps full precision for conversions within 10⁻¹² of 0 or 1.', 'With a species table the reaction is assembled from formation data; species with ν = 0 are inerts. Ideal-gas activities are y_iP/P°, solution activities C_i/(1 mol/L) or x_i.'],
    run(v) {
      const T0 = 298.15, kpi = [], warn = [], info = [];
      let sp, c, names, keyIdx, dnTot, spCp = null;
      if (v.mode === 'sp') {
        const rows = parseSp(v.sp, v.cpR);
        sp = rows.map(r => ({ nu: r.nu, n0: r.n0 }));
        names = rows.map(r => r.name);
        const S = f => N.sum(rows.map(r => r.nu * r[f]));
        c = TH.coef({ T0, H0: S('dHf'), G0: S('dGf'), da: S('a'), db: S('b'), dc: S('c'), dd: S('d') });
        keyIdx = rows.findIndex(r => r.nu < 0);
        if (!(rows[keyIdx].n0 > 0)) throw new Error('Feed some of the first reactant, ' + names[keyIdx] + '.');
        spCp = rows;
      } else {
        const [a, b, r, s] = FORMS[v.form];
        names = ['A', 'B', 'R', 'S', 'I'];
        const nu = [-a, -b, r, s, 0], n0 = [1, b ? +v.MB : 0, +v.MR || 0, s ? +v.MS || 0 : 0, v.std !== 'conc' ? +v.MI || 0 : 0];
        sp = nu.map((x, i) => ({ nu: x, n0: n0[i] }));
        const keep = sp.map((x, i) => x.nu !== 0 || x.n0 > 0);
        sp = sp.filter((x, i) => keep[i]); names = names.filter((x, i) => keep[i]);
        let da = 0, db = 0, dc = 0;
        if (v.cpmode === 'mean') {
          da = r * +v.CpR + (s ? s * +v.CpS : 0) - a * +v.CpA - (b ? b * +v.CpB : 0);
          spCp = [[-a, +v.CpA], [-b, +v.CpB], [r, +v.CpR], [s, +v.CpS]].filter(q => q[0]).map(q => ({ nu: q[0], a: q[1], b: 0, c: 0, d: 0 }));
        }
        if (v.cpmode === 'delta') { da = +v.da; db = +v.db; dc = +v.dc; }
        const lnK0 = v.Gsrc === 'K' ? Math.log(+v.K298) : -(+v.dG) / (R * T0);
        c = TH.coef({ T0, H0: +v.dH, da, db, dc, lnK0 });
        keyIdx = 0;
        if (b && !(+v.MB > 0)) throw new Error('Feed some B.');
      }
      dnTot = N.sum(sp.map(x => x.nu));
      const nuKey = -sp[keyIdx].nu, keyName = names[keyIdx];
      const basis = v.std, P = +v.P, Pstd = +v.Pstd || 1;
      const eqo = { basis, P, Pstd, cs: v.mode === 'rxn' ? +v.CA0 : 1, key: keyIdx };
      const eqAt = (T, o2) => TH.equil(sp, TH.lnK(c, T), Object.assign({}, eqo, o2 || {}));
      const T2 = C2K(+v.T2), Ta = C2K(+v.Ta), Tb = C2K(+v.Tb);
      if (!(Tb > Ta)) throw new Error('“Plot to” must be above “plot from”.');
      if (Ta <= 0) throw new Error('Temperatures must be above absolute zero.');
      const dH2 = TH.dH(c, T2), dG2 = TH.dG(c, T2), lnK2 = TH.lnK(c, T2);
      const e2 = eqAt(T2);
      kpi.push({ k: 'dH2', l: 'ΔH_r at ' + N.fmt(+v.T2, 5) + ' °C', v: dH2, u: 'J/mol', main: true, s: dH2 < 0 ? 'exothermic' : 'endothermic' });
      kpi.push({ k: 'K2', l: 'K at ' + N.fmt(+v.T2, 5) + ' °C', v: Math.exp(lnK2), s: 'ΔG° = ' + N.fmt(dG2, 4) + ' J/mol' });
      kpi.push({ k: 'Xe2', l: 'X_e of ' + keyName + ' at ' + N.fmt(+v.T2, 5) + ' °C', v: e2 ? e2.X : NaN, main: true });
      kpi.push({ k: 'K298', l: 'K at 25 °C', v: Math.exp(c.lnK0) }, { k: 'dH298', l: 'ΔH°_r at 25 °C', v: c.H0, u: 'J/mol' });
      if (!c.flat) kpi.push({ k: 'dCp2', l: 'ΔC_p at ' + N.fmt(+v.T2, 5) + ' °C', v: TH.dCp(c, T2), u: 'J/(mol·K)' });
      /* temperature limit for a wanted conversion */
      const Ts = N.linspace(Ta, Tb, 400);
      const Xs = Ts.map(T => { const e = eqAt(T); return e ? e.X : NaN; });
      const Xt = +v.Xt;
      let Tlim = NaN, side = '';
      for (let i = 1; i < Ts.length; i++) if ((Xs[i - 1] - Xt) * (Xs[i] - Xt) <= 0 && isFinite(Xs[i]) && isFinite(Xs[i - 1])) {
        Tlim = N.brent(T => eqAt(T).X - Xt, Ts[i - 1], Ts[i], { tol: 1e-10 });
        side = Xs[i] < Xs[i - 1] ? 'below' : 'above';
        break;
      }
      kpi.push({ k: 'Tlim', l: 'X_e ≥ ' + N.fmt(Xt, 3) + ' needs T ' + (side || 'in a range not plotted'), v: isFinite(Tlim) ? K2C(Tlim) : NaN, u: '°C', main: true });
      if (!isFinite(Tlim)) info.push('The wanted equilibrium conversion is not crossed between ' + N.fmt(+v.Ta, 4) + ' and ' + N.fmt(+v.Tb, 4) + ' °C.');
      info.push('ΔH_r ' + (dH2 < 0 ? '< 0: equilibrium conversion falls as T rises.' : '> 0: equilibrium conversion rises with T.') + (dnTot ? ' Δn = ' + N.fmt(dnTot, 3) + (basis === 'gas' ? (dnTot < 0 ? ': higher pressure raises conversion.' : ': higher pressure lowers conversion.') : '') : ' Δn = 0: pressure and inerts do not move the equilibrium of a gas.'));
      /* adiabatic equilibrium */
      let adPt = null;
      if (v.adOn) {
        const Tf = C2K(+v.Tf), Cpf = +v.Cpf;
        const xline = T => Cpf * (T - Tf) / (-TH.dH(c, T) / nuKey);
        const g = T => xline(T) - eqAt(T).X;
        const dir = TH.dH(c, Tf) < 0 ? 1 : -1;
        const span = N.linspace(Tf, Math.max(60, Tf + dir * 2500), 600);
        let Tad = NaN;
        for (let i = 1; i < span.length; i++) { const g0 = g(span[i - 1]), g1 = g(span[i]); if (isFinite(g0) && isFinite(g1) && g0 * g1 <= 0) { Tad = N.brent(g, span[i - 1], span[i], { tol: 1e-9 }); break; } }
        if (isFinite(Tad)) {
          adPt = { T: Tad, X: eqAt(Tad).X };
          kpi.push({ k: 'Tad', l: 'Adiabatic equilibrium temperature', v: K2C(Tad), u: '°C', s: 'X = ' + N.fmt(adPt.X, 4) });
          kpi.push({ k: 'Xad', l: 'Adiabatic equilibrium conversion', v: adPt.X });
        } else warn.push('The adiabatic line from ' + N.fmt(+v.Tf, 4) + ' °C does not meet the equilibrium curve.');
      }
      /* plots */
      const TsC = Ts.map(K2C);
      const xeSeries = [{ name: 'X_e at P = ' + N.fmt(P, 3) + ' atm', x: TsC, y: Xs, color: 1, w: 2 }];
      if (basis === 'gas' && dnTot) [0.1, 10].forEach((f, i) => xeSeries.push({ name: 'P = ' + N.fmt(P * f, 3) + ' atm', x: TsC, y: Ts.map(T => { const e = eqAt(T, { P: P * f }); return e ? e.X : NaN; }), color: i ? 3 : 6, dash: '5 4' }));
      if (basis !== 'conc' && dnTot && v.mode === 'rxn') {
        const mi = +v.MI || 0, more = 4 * (1 + mi);
        xeSeries.push({ name: N.fmt(more, 3) + ' mol inerts per mol A', x: TsC, y: Ts.map(T => { const sp2 = sp.map((x, i) => (names[i] === 'I' ? { nu: 0, n0: more } : x)); const e = TH.equil(names.indexOf('I') >= 0 ? sp2 : sp2.concat([{ nu: 0, n0: more }]), TH.lnK(c, T), eqo); return e ? e.X : NaN; }), color: 5, dash: '2 3' });
      }
      const marks = [{ t: 'h', y: Xt, color: 8, label: 'X = ' + N.fmt(Xt, 3) }];
      if (isFinite(Tlim)) marks.push({ t: 'v', x: K2C(Tlim), color: 2, label: N.fmt(K2C(Tlim), 4) + ' °C' });
      if (adPt) {
        const Tf = C2K(+v.Tf), Cpf = +v.Cpf;
        const TT = N.linspace(Tf, adPt.T, 40);
        xeSeries.push({ name: 'adiabatic line', x: TT.map(K2C), y: TT.map(T => Cpf * (T - Tf) / (-TH.dH(c, T) / nuKey)), color: 4, w: 1.6 });
        marks.push({ t: 'pt', x: K2C(adPt.T), y: adPt.X, color: 4, label: 'adiabatic equilibrium' });
      }
      const plots = [{ title: 'Equilibrium conversion of ' + keyName + ' (Fig. E9.2)', x: { label: 'T', unit: '°C', min: +v.Ta, max: +v.Tb }, y: { label: 'X_{' + keyName + 'e}', min: 0, max: 1 }, series: xeSeries, marks }];
      plots.push({ title: 'Heat of reaction and free energy', x: { label: 'T', unit: '°C', min: +v.Ta, max: +v.Tb }, y: { label: 'energy', unit: 'kJ/mol' }, series: [{ name: 'ΔH_r', x: TsC, y: Ts.map(T => TH.dH(c, T) / 1000), color: 1 }, { name: 'ΔG°', x: TsC, y: Ts.map(T => TH.dG(c, T) / 1000), color: 2 }, { name: 'TΔS°', x: TsC, y: Ts.map(T => T * TH.dS(c, T) / 1000), color: 3, dash: '5 4' }], marks: [{ t: 'h', y: 0, color: 8 }, { t: 'v', x: +v.T2, color: 8, label: N.fmt(+v.T2, 4) + ' °C' }] });
      const inv = Ts.map(T => 1000 / T);
      plots.push({ title: 'van ’t Hoff plot', x: { label: '1000/T', unit: '1/K' }, y: { label: 'ln K' }, series: [{ name: 'ln K', x: inv, y: Ts.map(T => TH.lnK(c, T)), color: 1 }, { name: 'constant-ΔH_r line through 25 °C', x: inv, y: Ts.map(T => c.lnK0 - c.H0 / R * (1 / T - 1 / T0)), color: 8, dash: '4 4' }], cap: 'Slope = −ΔH_r/R. A curved line means ΔC_p matters.' });
      /* Table E9.2 */
      const tT = N.linspace(+v.Ta, +v.Tb, 11);
      const rows = tT.map(tc => { const T = C2K(tc), e = eqAt(T); return [tc, Math.exp(TH.lnK(c, T)), e ? e.X : NaN, TH.dH(c, T), TH.dG(c, T)]; });
      const tables = [{ title: 'K and X_e against temperature (Table E9.2)', cols: ['T, °C', 'K', 'X_{' + keyName + 'e}', 'ΔH_r, J/mol', 'ΔG°, J/mol'], rows }];
      if (e2) tables.push({ title: 'Equilibrium mixture at ' + N.fmt(+v.T2, 5) + ' °C', cols: ['species', 'ν', 'feed', 'at equilibrium', basis === 'gas' ? 'y_i' : basis === 'x' ? 'x_i' : 'C_i, mol/L'], rows: sp.map((x, i) => [names[i], x.nu, x.n0, e2.n[i], basis === 'conc' ? e2.n[i] * eqo.cs : e2.n[i] / e2.nt]) });
      /* enthalpy path (Fig. E9.1) */
      const intCp = (q, A, B) => q.a * (B - A) + q.b / 2 * (B * B - A * A) + q.c / 3 * (B ** 3 - A ** 3) - q.d * (1 / B - 1 / A);
      let sens = null;
      if (spCp) { let re = 0, pr = 0; spCp.forEach(x => { if (x.nu < 0) re += -x.nu * intCp(x, T2, T0); else if (x.nu > 0) pr += x.nu * intCp(x, T0, T2); }); sens = [re, pr]; }
      const d1 = c.H0, kj = x => N.fmt(x / 1000, 4) + ' kJ';
      const tc = N.fmt(+v.T2, 4) + ' °C';
      let ds = D.tag(95, 30, 'reactants, ' + tc) + D.tag(385, 30, 'products, ' + tc) + D.tag(95, 142, 'reactants, 25 °C') + D.tag(385, 142, 'products, 25 °C');
      ds += D.pipe([[95, 40], [95, 131]], { flow: false }) + D.pipe([[160, 142], [318, 142]], { flow: false }) + D.pipe([[385, 131], [385, 40]], { flow: false }) + D.pipe([[160, 30], [318, 30]], { flow: false });
      ds += D.text(239, 134, 'ΔH_r = ' + kj(d1), { fs: 10.5 }) + D.text(239, 22, 'ΔH_r = ' + kj(dH2), { fs: 10.5 });
      if (sens) ds += D.text(103, 90, 'ΔH₁ = ' + kj(sens[0]), { anchor: 'start', fs: 10 }) + D.text(377, 90, 'ΔH₃ = ' + kj(sens[1]), { anchor: 'end', fs: 10 });
      else ds += D.text(239, 90, '∫ΔC_p dT = ' + kj(dH2 - d1), { fs: 10 });
      const path = D.svg(480, 162, ds, 'enthalpy path: cool the reactants, react at 25 °C, heat the products');
      return { kpi, warn, info, plots, tables, diagram: path };
    }
  });

  /* ===================================================================== */
  /* H-902  the X–T chart                                                   */
  /* ===================================================================== */
  const EX93 = { form: 'A=R', rev: true, phase: 'liq', CA0: 1, MR: 0, ksrc: 'runs', r1T: 65, r1t: 1, r1X: 0.581, r2T: 25, r2t: 10, r2X: 0.6, Ksrc: 'thermo', dG: -14130, dH: -75300, dCp: 0, Cp: 1046, Tlo: 0, Thi: 100 };
  ST.tool({
    id: 'xtchart', tag: 'H-902', ch: 9, sec: 'Ch. 9.1', title: 'Conversion–temperature chart: rate contours, equilibrium, locus of maximum rates',
    lede: 'The working chart of nonisothermal design. Rate contours over the X–T plane for any single reaction, liquid or gas, reversible or not, with the equilibrium curve, the locus of maximum rates and any energy-balance line laid over it. Build the rate constants from two batch runs as in Example 9.3.',
    keys: 'conversion temperature chart rate contours equilibrium curve locus of maximum rates operating line energy balance reversible exothermic endothermic batch runs activation energy arrhenius example 9.3',
    units: { t: 'min' }, slow: true,
    inputs: XTI.inputs(EX93).concat([
      { group: 'Chart' },
      { k: 'Ta', label: 'Chart from', unit: '°C', v: 0, min: -150, max: 2000 },
      { k: 'Tb', label: 'Chart to', unit: '°C', v: 100, min: -150, max: 2500 },
      { k: 'levels', t: 'text', label: 'Contour levels (blank for automatic)', v: '', mono: true, hint: 'Rates of A, comma separated, e.g. 0.01, 0.1, 1, 10' },
      { k: 'pX', label: 'Probe: conversion', v: 0.5, min: 0, max: 0.9999 },
      { k: 'pT', label: 'Probe: temperature', unit: '°C', v: 65, min: -150, max: 2500 },
      { k: 'opOn', t: 'check', label: 'Draw an energy-balance line', v: true },
      { k: 'opT1', label: 'Feed temperature for that line', unit: '°C', v: 25, min: -150, max: 2500, if: v => v.opOn },
      { k: 'opQ', label: 'Heat added per mol A fed (negative to remove)', unit: 'J/mol', v: 0, min: -1e6, max: 1e6, if: v => v.opOn }
    ]),
    presets: [
      { name: 'Example 9.3 (C_A0 = 1)', v: Object.assign({}, EX93, { Ta: 0, Tb: 100, pX: 0.5, pT: 65, opOn: false }), book: { k1a: 0.942, k1b: { v: 0.0909, tol: 0.02 }, E1: 48900, Xep: 0.89 }, note: 'R-free A, 58.1% converted in 1 min at 65 °C and 60% in 10 min at 25 °C; reversible first-order kinetics with K from Example 9.2. The chart of Fig. E9.3.' },
      { name: 'Same, C_A0 = 4 with the book’s k₂', v: Object.assign({}, EX93, { CA0: 4, ksrc: 'arr', A1: Math.exp(17.34), E1: 48900, Ksrc: 'arr2', A2: Math.exp(42.04), E2: 124200, opOn: true, opT1: 16, opQ: 0, pX: 0.8, pT: 62 }), note: 'k₁ = exp(17.34 − 48 900/RT), k₂ = exp(42.04 − 124 200/RT) min⁻¹. The energy line from 16 °C is the best adiabatic plug-flow line of Example 9.6.' },
      { name: 'Endothermic reversible (gas)', v: { form: 'A=R+S', rev: true, phase: 'gas', P: 1, yA0: 0.5, MR: 0, MS: 0, ksrc: 'arr', A1: 1e10, E1: 120000, Ksrc: 'KT', KT: 0.05, TK: 400, dH: 90000, dCp: 0, Cp: 120, Tlo: 300, Thi: 700, Ta: 300, Tb: 700, pX: 0.5, pT: 550, opOn: true, opT1: 650, opQ: 0 }, note: 'Rates and equilibrium both rise with T, so the best path is the hottest allowed. The adiabatic line runs to the left.' },
      { name: 'Irreversible second order', v: { form: 'A+B=R', rev: false, oA: 1, oB: 1, phase: 'liq', CA0: 2, MB: 1.2, ksrc: 'ref', k1r: 0.05, Tk: 25, E1: 60000, dH: -60000, dCp: 0, Cp: 1500, Tlo: 0, Thi: 120, Ta: 0, Tb: 120, pX: 0.6, pT: 60, opOn: true, opT1: 20, opQ: 0 }, note: 'No equilibrium limit: rates always climb with temperature, so the locus of maximum rates sits on the upper limit.' }
    ],
    eqs: [
      ['-r_A = k_1 C_A - k_2 C_R = k_1 C_{A0}\\left[(1-X_A) - \\frac{M + X_A}{K}\\right]', 'Ex. 9.3'],
      ['-\\ln\\left(1-\\frac{X_A}{X_{Ae}}\\right) = \\frac{k_1 t}{X_{Ae}}', 'Eq. 3.54', 'Batch, reversible first order, R-free feed.'],
      ['k_1 = A_1 e^{-E_1/RT},\\quad K = k_1/k_2', 'Ex. 9.3'],
      ['\\left.\\frac{\\partial(-r_A)}{\\partial T}\\right|_{X_A} = 0', 'Fig. 9.5', 'Locus of maximum rates.'],
      ['X_A = \\frac{C\'_p\\,\\Delta T - Q}{-\\Delta H_{r2}}', 'Eq. 26', 'Energy-balance line with heat Q added per mole of A.']
    ],
    method: ['General elementary reversible rate written per mole of A, −r_A = a k₁[ΠC_reac^ν − ΠC_prod^ν/K_c], with K_c from thermodynamics (gas K_p converted by (RT)^−Δn) or from k₁/k₂. Ideal gases at constant pressure carry the density change with temperature and moles.', 'Contours are drawn by marching squares on a log-rate grid. Topmost points of each contour satisfy ∂r/∂T = 0, so the labels line up along the locus of maximum rates, as in the book’s figure.', 'From batch runs, k₁ at each temperature follows from t = C_A0∫dX/(−r_A) with K(T) known; two temperatures give E₁ and A₁.'],
    run(v) {
      const m = XTI.model(v), kpi = [], warn = [], info = [];
      const Ta = C2K(+v.Ta), Tb = C2K(+v.Tb);
      if (!(Tb > Ta)) throw new Error('“Chart to” must be above “chart from”.');
      const ru = 'mol/(L·' + (v._tu || 'min') + ')';
      if (m.fit) {
        kpi.push({ k: 'k1a', l: 'k₁ at ' + N.fmt(+v.r1T, 4) + ' °C (run 1)', v: m.fit.ka, u: XTI.kUnit(v), main: true });
        kpi.push({ k: 'k1b', l: 'k₁ at ' + N.fmt(+v.r2T, 4) + ' °C (run 2)', v: m.fit.kb, u: XTI.kUnit(v) });
        kpi.push({ k: 'E1', l: 'E₁ from the two runs', v: m.fit.E1, u: 'J/mol', main: true }, { k: 'lnA1', l: 'ln A₁', v: Math.log(m.fit.A1) });
      }
      const pT = C2K(+v.pT), pX = +v.pX;
      const xe = m.Xe(pT);
      kpi.push({ k: 'Xep', l: 'X_{Ae} at ' + N.fmt(+v.pT, 4) + ' °C', v: xe, main: !m.fit });
      if (m.rev) kpi.push({ k: 'Kp', l: 'K at ' + N.fmt(+v.pT, 4) + ' °C', v: m.Kc(pT), u: kEqUnit(Object.assign({}, v, { phase: 'liq' })) });
      kpi.push({ k: 'k1p', l: 'k₁ at ' + N.fmt(+v.pT, 4) + ' °C', v: m.k1(pT), u: XTI.kUnit(v) });
      const rp = m.rate(pX, pT);
      kpi.push({ k: 'rp', l: '−r_A at the probe', v: rp, u: ru, s: 'X = ' + N.fmt(pX, 3) + ', ' + N.fmt(+v.pT, 4) + ' °C', state: rp < 0 ? 'bad' : undefined });
      const To = m.Topt(pX, Ta, Tb), rmx = m.rate(pX, To);
      kpi.push({ k: 'Topt', l: 'Fastest temperature at X = ' + N.fmt(pX, 3), v: K2C(To), u: '°C', s: '−r_A = ' + N.fmt(rmx, 4) + ' ' + ru + (Math.abs(To - Tb) < 0.05 ? ' (at the chart edge)' : '') });
      kpi.push({ k: 'dHA', l: 'ΔH_r per mol A at the probe', v: m.dHA(pT), u: 'J/mol' });
      if (rp < 0) warn.push('The probe point lies beyond equilibrium: the reaction runs backwards there.');
      /* chart */
      const lv = String(v.levels || '').split(/[,;\s]+/).map(parseFloat).filter(x => x > 0);
      const extra = [], marks = [{ t: 'pt', x: +v.pT, y: pX, color: 8, label: 'probe' }];
      if (v.opOn) {
        const T1 = C2K(+v.opT1), Q = +v.opQ;
        const xs = N.linspace(0, Math.min(0.999, m.Xmax), 200);
        const pts = xs.map(X => [m.Tline(X, T1, Q), X]).filter(p => p[0] > 0);
        extra.push({ name: 'energy balance from ' + N.fmt(+v.opT1, 4) + ' °C', x: pts.map(p => K2C(p[0])), y: pts.map(p => p[1]), color: 6, w: 2 });
        const slope = m.Cp / (-m.dHA(T1));
        info.push('Energy-balance line: slope dX/dT = C′_p/(−ΔH_r) = 1/' + N.fmt(1 / slope, 4) + ' per K.');
        /* where it meets equilibrium */
        if (m.rev) {
          const g = X => m.rate(X, m.Tline(X, T1, Q));
          const Xc = N.root(g, 1e-9, m.Xmax * 0.999999, { n: 400 });
          if (isFinite(Xc)) { kpi.push({ k: 'Xlim', l: 'Line meets equilibrium at', v: Xc, s: N.fmt(K2C(m.Tline(Xc, T1, Q)), 4) + ' °C' }); marks.push({ t: 'pt', x: K2C(m.Tline(Xc, T1, Q)), y: Xc, color: 6, open: true }); }
        }
      }
      const chart = XTI.chartSpec(m, { Ta, Tb, levels: lv, extra, marks, ru });
      chart.x.min = +v.Ta; chart.x.max = +v.Tb;
      const plots = [chart];
      /* rate against T at fixed X (Fig. 9.2) */
      const Tg = N.linspace(Ta, Tb, 220), TgC = Tg.map(K2C);
      const xTop = m.rev ? m.Xe(Ta) : m.Xmax;
      const xl = [0, 0.2, 0.4, 0.6, 0.8, 0.9].filter(x => x < xTop * 0.995).slice(0, 6);
      plots.push({ title: 'Rate against temperature at fixed conversion (Fig. 9.2)', x: { label: 'T', unit: '°C', min: +v.Ta, max: +v.Tb }, y: { label: '−r_A', unit: ru, min: 0 }, series: xl.map((x, i) => ({ name: 'X_A = ' + x, x: TgC, y: Tg.map(T => { const r = m.rate(x, T); return r > 0 ? r : NaN; }), color: i + 1 })), cap: 'Reversible exothermic reactions show a maximum: the locus of maximum rates joins these peaks.' });
      const tl = N.linspace(Ta, Tb, 5);
      plots.push({ title: 'Rate against conversion at fixed temperature', x: { label: 'X_A', min: 0, max: Math.min(1, m.Xmax) }, y: { label: '−r_A', unit: ru, min: 0, log: false }, series: tl.map((T, i) => { const top = m.rev ? m.Xe(T) : m.Xmax * 0.999; const xs = N.linspace(0, top, 160); return { name: N.fmt(K2C(T), 4) + ' °C', x: xs, y: xs.map(x => Math.max(m.rate(x, T), 0)), color: i + 1 }; }) });
      const inv = Tg.map(T => 1000 / T);
      const arr = [{ name: 'ln k₁', x: inv, y: Tg.map(T => Math.log(m.k1(T))), color: 1 }];
      if (m.rev) arr.push({ name: 'ln k₂', x: inv, y: Tg.map(T => Math.log(m.k2(T))), color: 2 }, { name: 'ln K', x: inv, y: Tg.map(T => m.lnKc(T)), color: 3, dash: '5 4' });
      plots.push({ title: 'Arrhenius and van ’t Hoff lines', x: { label: '1000/T', unit: '1/K' }, y: { label: 'ln' }, series: arr });
      /* tables */
      const tT = N.linspace(+v.Ta, +v.Tb, 11);
      const t1 = { title: 'Rate constants and equilibrium', cols: ['T, °C', 'k₁', 'k₂', 'K', 'X_{Ae}', 'ΔH_r per mol A'], rows: tT.map(tc => { const T = C2K(tc); return [tc, m.k1(T), m.rev ? m.k2(T) : NaN, m.rev ? m.Kc(T) : NaN, m.Xe(T), m.dHA(T)]; }) };
      const tables = [t1];
      if (m.rev && m.exo) {
        const xs = N.linspace(0, Math.min(0.95, m.Xe(Ta) * 0.98), 11);
        tables.push({ title: 'Locus of maximum rates', cols: ['X_A', 'T, °C', '−r_A max', 'X_{Ae} at that T'], rows: xs.map(x => { const T = m.Topt(x, Ta - 300 > 100 ? Ta - 300 : 100, Tb + 300); return [x, K2C(T), m.rate(x, T), m.Xe(T)]; }) });
      }
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
