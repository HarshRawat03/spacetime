/* Chapter 22 — G/L reactions on solid catalysts: trickle beds, slurry reactors, bubble columns and
   three-phase fluidized beds. The general rate through gas film, liquid film, liquid–solid film and
   the porous catalyst, the two extremes, and the performance equations for every flow combination. */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, PEL = ST.PEL;
  const CAT = ST.CAT = ST.CAT || {};
  const fmt = (x, s = 4) => N.fmt(x, s);

  /* local rate per volume of reactor; returns the details of the resistance chain */
  function glRate(v, pA, CB) {
    const H = +v.HA, rhos = +v.rhos, fs = +v.fs, dp = +v.dp, De = +v.De, n = +v.n, m = +v.m, b = +v.b, kp = +v.kp;
    const ac = 6 * fs / dp, L = dp / 6;
    const CAsat = pA / H;
    /* gas + liquid film resistance on the liquid-concentration basis */
    const kgl = v.mt === 'sep' ? 1 / (1 / (H * +v.kAg) + 1 / +v.kAl) : +v.kgl;
    const kAc = +v.kAc;
    const CAbar = v.cbar === 'sat' ? CAsat : null;
    const chainA = CAm => {
      const kapp = kp * Math.pow(Math.max(CAm, 1e-300), n - 1) * Math.pow(Math.max(CB, 0), m) * rhos; /* per volume of catalyst, 1/s */
      const MT = L * Math.sqrt((n + 1) / 2 * kapp / De);
      const eta = v.eff === 'one' ? 1 : PEL.eta1('sph', MT);
      const R1 = 1 / kgl, R2 = 1 / (kAc * ac), R3 = 1 / (kapp * eta * fs);
      const r = CAsat / (R1 + R2 + R3);
      const CAs = CAsat - r * (R1 + R2);
      return { r, R: [R1, R2, R3], MT, eta, kapp, CAs };
    };
    let A = chainA(CAbar !== null ? CAbar : CAsat);
    if (CAbar === null) {
      /* self-consistent: pseudo-first-order constant at the concentration at the catalyst surface */
      let CAm = Math.max(A.CAs, CAsat * 1e-6);
      for (let i = 0; i < 60; i++) { const B = chainA(Math.max(CAm, CAsat * 1e-9)); const nx = Math.max(B.CAs, CAsat * 1e-9); if (Math.abs(nx - CAm) < 1e-10 * CAsat) { A = B; break; } CAm = 0.5 * (CAm + nx); A = B; }
    }
    /* extreme 2: B limiting (Eq. 22.5), A saturated in the liquid */
    const kBapp = b * kp * Math.pow(CAsat, n) * Math.pow(Math.max(CB, 1e-300), m - 1) * rhos;
    const MTB = L * Math.sqrt((m + 1) / 2 * kBapp / (+v.DeB || De));
    const etaB = v.eff === 'one' ? 1 : PEL.eta1('sph', MTB);
    const rB = CB / (1 / ((+v.kBc || kAc) * ac) + 1 / (kBapp * etaB * fs));
    return { rA: A.r, rB, A, CAsat, ac, kgl, MTB, etaB, limiting: A.r * b <= rB ? 'A' : 'B', r: Math.min(A.r, rB / b) };
  }

  ST.tool({
    id: 'glcat', tag: 'C-2201', ch: 22, sec: 'Ch. 22', title: 'G/L reactions on solid catalysts: trickle beds, slurries, three-phase beds',
    lede: 'A gas dissolves, crosses the liquid and the film around the catalyst, and reacts inside porous pellets with a liquid reactant. This tool assembles the resistances of Eq. 22.2 (with the pellet’s effectiveness factor for any order), tells which reactant really limits, and then sizes the contactor for every flow combination of the chapter: gas in plug or mixed flow or pure and recycled, liquid in plug, mixed or batch — conversion, or the batch time.',
    keys: 'three phase gas liquid solid catalyst trickle bed slurry reactor bubble column packed bubble column three phase fluidized bed hydrogenation henry law liquid film gas film liquid solid mass transfer effectiveness factor extreme excess of B excess of A example 22.1 22.2 problem 22.1 22.2 22.3 22.4 22.5 22.6 22.8 acetone butynediol glucose sorbitol aniline ethanol oxidation',
    inputs: [
      { group: 'Reaction A(g → l) + bB(l) → products, −r′_A = k′C_Aⁿ C_Bᵐ per kg catalyst' },
      { k: 'kp', label: 'k′', unit: 'm³…/(kg·s)', v: 2.35e-3, min: 1e-15, max: 1e6, scale: 'log' },
      { k: 'n', label: 'Order in A, n', v: 0.5, min: 0, max: 3, step: 0.1 },
      { k: 'm', label: 'Order in B, m', v: 0, min: 0, max: 3, step: 0.1 },
      { k: 'b', label: 'b, mol B per mol A', v: 1, min: 0.01, max: 20 },
      { group: 'Gas' },
      { k: 'gas', t: 'select', label: 'Gas', opts: [['pure', 'pure A, recycled: p_A constant'], ['plug', 'A with inerts, plug flow'], ['mixed', 'A with inerts, mixed flow']], v: 'pure' },
      { k: 'pi', label: 'Total pressure π', unit: 'Pa', v: 101325, min: 1, max: 1e9, scale: 'log' },
      { k: 'yA', label: 'Mole fraction of A entering', v: 1, min: 1e-6, max: 1, if: v => v.gas !== 'pure' },
      { k: 'vg', label: 'Gas flow v_g', unit: 'm³/s', v: 0.04, min: 1e-12, max: 1e6, if: v => v.gas !== 'pure' },
      { k: 'HA', label: 'Henry constant H_A = p_A/C_A', unit: 'Pa·m³/mol', v: 36845, min: 1e-6, max: 1e12, scale: 'log' },
      { group: 'Liquid' },
      { k: 'liq', t: 'select', label: 'Liquid', opts: [['mixed', 'continuous, mixed (tanks, slurry columns)'], ['plug', 'continuous, plug (trickle beds, towers)'], ['batch', 'batch']], v: 'mixed' },
      { k: 'vl', label: 'Liquid flow v_l', unit: 'm³/s', v: 1e-4, min: 1e-15, max: 1e6, if: v => v.liq !== 'batch' },
      { k: 'CB0', label: 'C_B0', unit: 'mol/m³', v: 1000, min: 1e-9, max: 1e6 },
      { k: 'Vl', label: 'Volume of the batch of liquid', unit: 'm³', v: 2, min: 1e-9, max: 1e6, if: v => v.liq === 'batch' },
      { k: 'XBt', label: 'Conversion of B wanted (batch)', v: 0.9, min: 0.001, max: 0.9999, if: v => v.liq === 'batch' },
      { group: 'Contactor and catalyst' },
      { k: 'Vr', label: 'Reactor volume V_r', unit: 'm³', v: 0.5, min: 1e-9, max: 1e6 },
      { k: 'fs', label: 'Solid loading f_s (m³ particles/m³ reactor)', v: 0.6, min: 1e-6, max: 0.7 },
      { k: 'dp', label: 'Particle diameter d_p', unit: 'm', v: 5e-3, min: 1e-8, max: 0.1, scale: 'log' },
      { k: 'rhos', label: 'Particle density ρ_s', unit: 'kg/m³', v: 4500, min: 1, max: 3e4 },
      { k: 'De', label: 'Effective diffusivity of A, D_e', unit: 'm³/(m cat·s)', v: 8e-10, min: 1e-20, max: 1, scale: 'log' },
      { k: 'eff', t: 'seg', label: 'Pellet effectiveness', opts: [['calc', 'from M_T (Eq. 18.28)'], ['one', 'take ℰ = 1']], v: 'calc' },
      { k: 'cbar', t: 'seg', label: 'C̄_A in the pseudo-first-order constant (Eq. 22.17)', opts: [['sat', 'saturation value (as in the book)'], ['self', 'self-consistent at the catalyst']], v: 'sat' },
      { group: 'Mass transfer' },
      { k: 'mt', t: 'seg', label: 'Gas–liquid coefficients', opts: [['comb', '(k_A a_i)_{g+l}'], ['sep', 'k_Ag a_i and k_Al a_i']], v: 'comb' },
      { k: 'kgl', label: '(k_A a_i)_{g+l}', unit: '1/s', v: 0.02, min: 1e-12, max: 1e6, if: v => v.mt === 'comb' },
      { k: 'kAg', label: 'k_Ag a_i', unit: 'mol/(m³·Pa·s)', v: 3e-4, min: 1e-15, max: 1e6, if: v => v.mt === 'sep' },
      { k: 'kAl', label: 'k_Al a_i', unit: '1/s', v: 0.02, min: 1e-12, max: 1e6, if: v => v.mt === 'sep' },
      { k: 'kAc', label: 'k_Ac (liquid to particle)', unit: 'm/s', v: 0.05 / (6 * 0.6 / 5e-3), min: 1e-15, max: 10, scale: 'log' }
    ],
    presets: [
      { name: 'Example 22.1 — acetone hydrogenation, packed bubble column', v: { kp: 2.35e-3, n: 0.5, m: 0, b: 1, gas: 'pure', pi: 101325, HA: 36845, liq: 'mixed', vl: 1e-4, CB0: 1000, Vr: 0.5, fs: 0.6, dp: 5e-3, rhos: 4500, De: 8e-10, eff: 'calc', cbar: 'sat', mt: 'comb', kgl: 0.02, kAc: 0.05 / (6 * 0.6 / 5e-3) }, book: { r: { v: 0.0317, tol: 0.01 }, XB: { v: 0.158, tol: 0.02 }, MT: { v: 64.4, tol: 0.01 }, eta: { v: 0.0155, tol: 0.01 }, sh1: { v: 0.58, tol: 0.02 } }, note: 'Aqueous acetone (1000 mol/m³, 10⁻⁴ m³/s) and pure H₂ at 1 atm (H_A = 36 845 Pa·m³/mol) through a 5 m, 0.1 m² column of 5 mm Raney nickel (f_s = 0.6). (k_A a_i)_{g+l} = 0.02 s⁻¹, k_Ac a_c = 0.05 s⁻¹. Book: M_T = 64.4, ℰ = 0.0155, −r‴ = 0.0317 mol/m³·s with 58%, 23% and 19% of the resistance in the three steps, X_B = 16%.' },
      { name: 'Example 22.2 — butynediol, batch slurry', v: { kp: 5e-5, n: 1, m: 1, b: 1, gas: 'pure', pi: 14.6 * 101325, HA: 148000, liq: 'batch', CB0: 2500, Vl: 2, XBt: 0.9, Vr: 2, fs: 0.0055, dp: 5e-5, rhos: 1450, De: 5e-10, eff: 'calc', cbar: 'sat', mt: 'comb', kgl: 0.277, kAc: 4.4e-4 }, book: { t: { v: 3460, tol: 0.03 } }, note: 'Palladium on 50 μm pellets (f_s = 0.0055) in a 2 m³ tank of butynediol (2500 mol/m³), pure H₂ recirculated. The text says 1.46 atm, but Fig. E22.2a and the book’s own arithmetic use 14.6 atm (C_A = 10 mol/m³), as here. 90% conversion takes 3460 s, 58 min.' },
      { name: 'Problem 22.1 — ethanol oxidation in a trickle bed', v: { kp: 1.77e-5, n: 1, m: 0, b: 1, gas: 'pure', pi: 10 * 101325, HA: 86000, liq: 'plug', vl: 2e-4, CB0: 400, Vr: 0.5, fs: 0.58, dp: 5e-3, rhos: 1800, De: 4.16e-10, eff: 'calc', cbar: 'sat', mt: 'sep', kAg: 3e-4, kAl: 0.02, kAc: 3.86e-4 }, note: 'Pure oxygen at 10 atm and 2–3% aqueous ethanol over palladium–alumina pellets.' },
      { name: 'Problem 22.2 — the same oxidation in a slurry column', v: { kp: 1.77e-5, n: 1, m: 0, b: 1, gas: 'pure', pi: 10 * 101325, HA: 86000, liq: 'mixed', vl: 2e-4, CB0: 400, Vr: 0.5, fs: 0.2, dp: 1e-4, rhos: 1800, De: 4.16e-10, eff: 'calc', cbar: 'sat', mt: 'comb', kgl: 0.052, kAc: 4e-4 } },
      { name: 'Problem 22.3 — glucose to sorbitol, slurry tank', v: { kp: 5.96e-6, n: 0.6, m: 1, b: 1, gas: 'pure', pi: 200 * 101325, HA: 277600, liq: 'mixed', vl: 0.01, CB0: 2000, Vr: 2, fs: 0.056, dp: 1e-5, rhos: 8900, De: 2e-9, eff: 'calc', cbar: 'sat', mt: 'comb', kgl: 0.05, kAc: 1e-4 } },
      { name: 'Problem 22.4 — multistage bubble column', v: { kp: 5.96e-6, n: 0.6, m: 1, b: 1, gas: 'pure', pi: 200 * 101325, HA: 277600, liq: 'plug', vl: 0.01, CB0: 2000, Vr: 2, fs: 0.4, dp: 1e-3, rhos: 8900, De: 2e-9, eff: 'calc', cbar: 'sat', mt: 'comb', kgl: 0.025, kAc: 1e-5 } },
      { name: 'Problem 22.5 — aniline in a three-phase fluidized bed', v: { kp: 0.05, n: 1, m: 0, b: 2 / 3, gas: 'pure', pi: 101325, HA: 28500, liq: 'batch', CB0: 1097, Vl: 0.65, XBt: 0.9, Vr: 1, fs: 0.25, dp: 3e-4, rhos: 750, De: 8.35e-10, eff: 'calc', cbar: 'sat', mt: 'comb', kgl: 0.04, kAc: 1e-5 }, note: '3H₂ + 2C₆H₅NH₂ → products: b = 2/3. Per m³ of bed, f_l = 0.65 m³ of liquid.' },
      { name: 'Problem 22.6 — aniline in a bubble column with an external loop', v: { kp: 0.05, n: 1, m: 0, b: 2 / 3, gas: 'pure', pi: 101325, HA: 28500, liq: 'batch', CB0: 1097, Vl: 1.0, XBt: 0.9, Vr: 1, fs: 0.4, dp: 3e-3, rhos: 750, De: 8.35e-10, eff: 'calc', cbar: 'sat', mt: 'comb', kgl: 0.02, kAc: 7e-5 }, note: 'The loop outside holds as much liquid as the column (f_l = 0.5), so the batch is twice the liquid in the reactor.' }
    ],
    eqs: [
      ['-r_A^{\\prime\\prime\\prime} = \\cfrac{p_{Ag}}{\\cfrac{1}{k_{Ag}a_i}+\\cfrac{H_A}{k_{Al}a_i}+\\cfrac{H_A}{k_{Ac}a_c}+\\cfrac{H_A}{(k_A^{\\prime\\prime\\prime}\\bar C_B)\\mathscr{E}_Af_s}}', 'Eq. 22.2'],
      ['-r_B^{\\prime\\prime\\prime} = \\cfrac{C_{Bl}}{\\cfrac{1}{k_{Bc}a_c}+\\cfrac{1}{(k_B^{\\prime\\prime\\prime}\\bar C_A)\\mathscr{E}_Bf_s}}', 'Eq. 22.3'],
      ['a_c = \\frac{6f_s}{d_p},\\qquad k^{\\prime\\prime\\prime} = \\rho_sk^\\prime,\\qquad H_A = \\frac{p_A}{C_A}', 'Ch. 22'],
      ['-r_A^{\\prime\\prime\\prime} = [(k_A^{\\prime\\prime\\prime}\\bar C_B^m)\\bar C_A^{n-1}]C_A', 'Eq. 22.17', 'Pseudo-first-order rate for general orders.'],
      ['F_{A0}X_A = \\frac{1}{b}F_{B0}X_B = (-r_A^{\\prime\\prime\\prime})V_r', 'Eq. 22.6', 'Mixed gas, any liquid flow (excess B).'],
      ['\\frac{V_r}{F_{A0}} = \\int_0^{X_A}\\frac{dX_A}{-r_A^{\\prime\\prime\\prime}},\\qquad 1-X_A = \\frac{p_A(\\pi-p_{A0})}{p_{A0}(\\pi-p_A)}', 'Eqs. 22.8, 22.10', 'Plug gas (excess B).'],
      ['t = \\frac{V_l}{bV_r}\\int_{C_{Bf}}^{C_{B0}}\\frac{dC_B}{-r_A^{\\prime\\prime\\prime}}', 'Eq. 22.12', 'Batch of liquid.'],
      ['\\frac{V_r}{F_{B0}} = \\int_0^{X_B}\\frac{dX_B}{-r_B^{\\prime\\prime\\prime}}\\ \\text{(plug L)},\\qquad \\frac{V_r}{F_{B0}} = \\frac{X_B}{-r_B^{\\prime\\prime\\prime}}\\ \\text{(mixed L)}', 'Eqs. 22.19–22.20', 'Excess of A.']
    ],
    method: ['At each point the rate through the A path (Eq. 22.2 with the pseudo-first-order constant of Eq. 22.17) and through the B path (Eq. 22.3) are both evaluated; the smaller one governs, which is the book’s test for the two extremes. The pellet effectiveness is the exact first-order sphere value at the generalized modulus M_T = L√((n+1)/2·k‴/D_e). C̄_A in the pseudo-first-order constant is taken at saturation as in the book, or iterated to the concentration at the catalyst surface.', 'Performance: with pure, recycled A the rate depends only on C_B, so mixed liquid is a root of F_B0 X_B = b(−r)V_r, plug liquid integrates along the reactor, and a batch integrates dC_B/dt = −b(−r)V_r/V_l. With A diluted in an inert gas the gas balance (Eqs. 22.6–22.10) is used with C_B held at its feed value — the excess-B case of the book.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const b = +v.b, Vr = +v.Vr, CB0 = +v.CB0;
      const pA0 = v.gas === 'pure' ? +v.pi : +v.pi * +v.yA;
      const R0 = glRate(v, pA0, CB0);
      const A = R0.A, tot = A.R[0] + A.R[1] + A.R[2];
      kpi.push({ k: 'CAsat', l: 'C_A in liquid saturated at the entering gas', v: R0.CAsat, u: 'mol/m³', s: CB0 > 3 * R0.CAsat ? 'C_B ≫ C_A: extreme 1, B in excess' : CB0 < R0.CAsat / 3 ? 'C_B ≪ C_A: extreme 2, A in excess' : 'neither extreme by concentrations' });
      kpi.push({ k: 'r', l: 'Rate at entry, −r‴_A (per m³ reactor)', v: R0.r, u: 'mol/(m³·s)', main: v.liq !== 'batch' && false, s: 'A path ' + fmt(R0.rA, 4) + '; B path ' + fmt(R0.rB / b, 4) + ' (as A) — ' + (R0.limiting === 'A' ? 'A limits (Eq. 22.4)' : 'B limits (Eq. 22.5)') });
      kpi.push({ k: 'MT', l: 'Thiele modulus for A (Eq. 18.28)', v: A.MT, s: CAT.regime ? CAT.regime(A.MT) : '' }, { k: 'eta', l: 'Effectiveness factor ℰ_A', v: A.eta });
      ['gas + liquid films', 'liquid to particle', 'reaction in pellets'].forEach((nm, i) => kpi.push({ k: 'sh' + (i + 1), l: 'Share of resistance: ' + nm, v: A.R[i] / tot }));
      const worst = A.R.indexOf(Math.max(...A.R));
      info.push(['Most of the resistance is in the gas–liquid films: favour contactors with large interfacial area (agitated slurries, bubble columns) or higher pressure.', 'The liquid-to-particle film dominates: raise the external area of solid — more solids or smaller particles.', 'The pellets dominate: ' + (A.eta < 0.3 ? 'pore diffusion is strong (ℰ = ' + fmt(A.eta, 3) + '), so smaller particles would pay' : 'more catalyst or a more active catalyst is the lever') + '.'][worst]);
      tables.push({ title: 'Resistances in series (Eq. 22.18 form, s)', cols: ['step', 'resistance', 'share'], rows: [['gas + liquid films 1/(k_A a_i)_{g+l}', A.R[0], A.R[0] / tot], ['liquid to particle 1/(k_Ac a_c)', A.R[1], A.R[1] / tot], ['reaction 1/(k‴C̄ⁿ⁻¹C̄_Bᵐ ℰ f_s)', A.R[2], A.R[2] / tot]] });
      const rAt = (pA, CB) => glRate(v, pA, CB).r;
      if (v.liq === 'batch') {
        if (v.gas !== 'pure') warn.push('For a batch of liquid the gas is treated as pure, recycled A (Example 22.2); the inert-gas cases use the same integral with p_A at the gas exit.');
        const Vl = +v.Vl, CBf = CB0 * (1 - +v.XBt);
        const t = Vl / (b * Vr) * N.quad(C => 1 / rAt(pA0, C), CBf, CB0, { rtol: 1e-8 });
        kpi.unshift({ k: 't', l: 'Time for X_B = ' + fmt(+v.XBt, 3), v: t, u: 's', main: true, s: fmt(t / 60, 4) + ' min' });
        const Cg = N.linspace(CBf, CB0, 120);
        plots.push({ title: 'Reaction time for the batch (Figs. 22.6, E22.2b)', x: { label: 'C_B', unit: 'mol/m³', min: 0 }, y: { label: '1/(−r‴_A)', unit: 'm³·s/mol', min: 0 }, series: [{ name: '1/(−r‴_A)', x: Cg, y: Cg.map(C => 1 / rAt(pA0, C)), color: 1, w: 2.2, type: 'area', line: true, fill: 0.16 }], cap: 'Shaded area × V_l/(bV_r) is the batch time.' });
        const rows = [CB0, CB0 * 0.4, CBf].map(C => { const q = glRate(v, pA0, C); return [C, q.A.MT, q.A.eta, q.r, 1 / q.r]; });
        tables.push({ title: 'Along the batch (as in Example 22.2)', cols: ['C_B', 'M_T', 'ℰ_A', '−r‴_A', '1/(−r‴_A)'], rows });
        return { kpi, warn, info, plots, tables };
      }
      const vl = +v.vl, FB0 = vl * CB0;
      if (v.gas === 'pure') {
        let XB;
        if (v.liq === 'mixed') {
          const g = X => FB0 * X - b * rAt(pA0, CB0 * (1 - X)) * Vr;
          XB = N.brent(g, 0, 1, { tol: 1e-14 });
        } else {
          const sol = N.ode((V, y) => [b * rAt(pA0, CB0 * (1 - Math.min(y[0], 1))) / FB0], 0, [0], Vr, { rtol: 1e-9, atol: 1e-14 });
          XB = sol.last()[0];
          const Vs = N.linspace(0, Vr, 80);
          plots.push({ title: 'Conversion of B along the contactor', x: { label: 'reactor volume passed', unit: 'm³', min: 0 }, y: { label: 'X_B', min: 0 }, series: [{ name: 'X_B', x: Vs, y: Vs.map(V => sol.at(V)[0]), color: 1, w: 2.2 }] });
        }
        kpi.unshift({ k: 'XB', l: 'Conversion of B', v: XB, main: true, s: (v.liq === 'mixed' ? 'liquid mixed' : 'liquid in plug flow') + ', pure A at constant p_A' });
        kpi.push({ k: 'FB0', l: 'F_B0', v: FB0, u: 'mol/s' });
        const Cg = N.linspace(CB0 * 0.01, CB0, 100);
        plots.push({ title: 'Rate against C_B at this gas pressure', x: { label: 'C_B', unit: 'mol/m³', min: 0 }, y: { label: '−r‴_A', unit: 'mol/(m³·s)', min: 0 }, series: [{ name: 'governing rate', x: Cg, y: Cg.map(C => rAt(pA0, C)), color: 1, w: 2.2 }, { name: 'A path (Eq. 22.4)', x: Cg, y: Cg.map(C => glRate(v, pA0, C).rA), color: 2, dash: '5 4' }, { name: 'B path ÷ b (Eq. 22.5)', x: Cg, y: Cg.map(C => glRate(v, pA0, C).rB / b), color: 3, dash: '2 3' }], cap: 'Where the dashed curves cross, the limiting reactant changes from B (left) to A (right).' });
        return { kpi, warn, info, plots, tables };
      }
      /* A with inerts: gas balance, B in excess (C_B ≈ C_B0) */
      const pi = +v.pi, FA0 = +v.vg * pA0 / (8.314462618 * 298.15);
      info.push('Gas molar flow is taken at 25 °C from v_g and the partial pressure; enter v_g at the conditions you mean.');
      const pOf = X => pi * pA0 * (1 - X) / (pi - pA0 * X);
      let XA;
      if (v.gas === 'mixed') XA = N.brent(X => FA0 * X - rAt(pOf(X), CB0) * Vr, 0, 1 - 1e-12, { tol: 1e-14 });
      else {
        const sol = N.ode((V, y) => [rAt(pOf(Math.min(y[0], 1 - 1e-12)), CB0) / FA0], 0, [0], Vr, { rtol: 1e-9, atol: 1e-14 });
        XA = sol.last()[0];
      }
      const XB = Math.min(1, FA0 * XA * b / FB0);
      kpi.unshift({ k: 'XA', l: 'Conversion of A in the gas', v: XA, main: true, s: 'gas in ' + v.gas + ' flow; B in excess' });
      kpi.push({ k: 'XB', l: 'Conversion of B', v: XB, s: XB > 0.2 ? 'not small: C_B is not constant, treat with care' : '' });
      if (XB > 0.2) warn.push('B is converted appreciably, so the excess-B assumption is stretched.');
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
