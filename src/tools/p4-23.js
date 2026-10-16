/* Chapter 23 — fluid–fluid reactions: kinetics. Gas film, liquid film and main body in series;
   the enhancement factor by van Krevelen–Hoftijzer (the book's Fig. 23.4), by Danckwerts surface
   renewal (DeCoursey), and by solving the film equations exactly; the eight regimes A–H. */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N;
  const fmt = (x, s = 4) => N.fmt(x, s);
  const FF = ST.FF = {};

  /* ---------- enhancement factors ---------- */
  FF.Evkh = function (M, Ei) {
    if (!(M > 0)) return 1;
    if (!isFinite(Ei)) return M > 300 ? M : M / Math.tanh(M);
    if (Ei <= 1 + 1e-12) return 1;
    if (!isFinite(M)) return Ei;
    const f = E => { const s = Math.sqrt(Math.max(Ei - E, 0) / (Ei - 1)), z = M * s; return E - (z < 1e-8 ? 1 : z > 300 ? z : z / Math.tanh(z)); };
    return N.brent(f, 1, Ei, { tol: 1e-12 * Ei });
  };
  FF.Edec = function (M, Ei) {
    if (!(M > 0)) return 1;
    if (!isFinite(Ei)) return Math.sqrt(1 + M * M);
    if (!isFinite(M)) return Ei;
    if (Ei <= 1 + 1e-12) return 1;
    const q = Ei - 1, M2 = M * M;
    return -M2 / (2 * q) + Math.sqrt(M2 * M2 / (4 * q * q) + Ei * M2 / q + 1);
  };

  /* exact film theory, A + bB, second order, no A in the bulk:
     a'' = M² a β,  β = 1 + (a − E(1 − x))/q,  a(0) = 1, a(1) = 0, E = −a'(0), q = E_i − 1 */
  function filmGrid(M, q) {
    const pts = [];
    const hmin = Math.min(1e-3, 0.02 / Math.max(M, 1));
    for (let x = hmin, h = hmin; x < 1; h *= 1.035, x += h) pts.push(x);
    for (let i = 1; i < 400; i++) pts.push(i / 400);
    if (isFinite(q) && M > 1) {
      const xr = 1 / (1 + q), w = Math.min(0.25, 8 * Math.cbrt(Math.max(q, 1e-6) / (M * M)) + 0.01);
      for (let i = 0; i <= 240; i++) { const x = xr - w + 2 * w * i / 240; if (x > 0 && x < 1) pts.push(x); }
    }
    pts.sort((a, b) => a - b);
    const x = [0];
    pts.forEach(p => { if (p - x[x.length - 1] > 1e-9 && 1 - p > 1e-9) x.push(p); });
    x.push(1);
    return x;
  }
  function filmSolveA(x, M, q, E, a) {
    const n = x.length - 1, M2 = M * M;
    for (let it = 0; it < 60; it++) {
      const lo = new Float64Array(n + 1), di = new Float64Array(n + 1), up = new Float64Array(n + 1), F = new Float64Array(n + 1);
      let nrm = 0;
      for (let i = 1; i < n; i++) {
        const h0 = x[i] - x[i - 1], h1 = x[i + 1] - x[i], s = 2 / (h0 + h1);
        const bRaw = 1 + (a[i] - E * (1 - x[i])) / q, b = bRaw > 0 ? bRaw : 0;
        F[i] = s * ((a[i + 1] - a[i]) / h1 - (a[i] - a[i - 1]) / h0) - M2 * a[i] * b;
        lo[i] = s / h0; up[i] = s / h1;
        di[i] = -s * (1 / h0 + 1 / h1) - M2 * (b + (bRaw > 0 ? a[i] / q : 0));
        nrm = Math.max(nrm, Math.abs(F[i]) / (Math.abs(di[i]) + 1e-300));
      }
      /* Thomas */
      const c = new Float64Array(n + 1), d = new Float64Array(n + 1);
      c[1] = up[1] / di[1]; d[1] = -F[1] / di[1];
      for (let i = 2; i < n; i++) { const m = di[i] - lo[i] * c[i - 1]; c[i] = up[i] / m; d[i] = (-F[i] - lo[i] * d[i - 1]) / m; }
      const dx = new Float64Array(n + 1);
      dx[n - 1] = d[n - 1];
      for (let i = n - 2; i >= 1; i--) dx[i] = d[i] - c[i] * dx[i + 1];
      let mx = 0;
      for (let i = 1; i < n; i++) { a[i] = Math.min(1, Math.max(0, a[i] + dx[i])); mx = Math.max(mx, Math.abs(dx[i])); }
      if (mx < 1e-12) break;
    }
    return a;
  }
  FF.filmExact = function (M, Ei, keep) {
    const q = Ei - 1;
    if (!(M > 0) || q <= 1e-12) return { E: 1 + (M > 0 ? 0 : 0) };
    if (!isFinite(q)) {
      const E = M > 300 ? M : M / Math.tanh(M);
      if (!keep) return { E };
      const xs = N.linspace(0, 1, 300);
      return { E, x: xs, a: xs.map(z => (M > 300 ? Math.exp(-M * z) : Math.sinh(M * (1 - z)) / Math.sinh(M))), b: xs.map(() => 1) };
    }
    if (!isFinite(M)) {
      const xr = 1 / Ei, xs = N.linspace(0, 1, 301).concat([xr]).sort((p, r) => p - r);
      return { E: Ei, x: xs, a: xs.map(z => Math.max(0, 1 - z / xr)), b: xs.map(z => (z < xr ? 0 : (z - xr) / (1 - xr))) };
    }
    const x = filmGrid(M, q), n = x.length - 1;
    let a = new Float64Array(n + 1);
    const xr = 1 / Ei;
    for (let i = 0; i <= n; i++) a[i] = Math.max(0, Math.min(1, M * xr > 3 ? Math.max(0, 1 - x[i] / xr) * Math.exp(-0 * x[i]) : Math.sinh(M * (1 - x[i])) / Math.sinh(M)));
    a[0] = 1; a[n] = 0;
    const slope = () => { const h0 = x[1] - x[0], h1 = x[2] - x[1]; return -(-(2 * h0 + h1) / (h0 * (h0 + h1)) * a[0] + (h0 + h1) / (h0 * h1) * a[1] - h0 / (h1 * (h0 + h1)) * a[2]); };
    const g = E => { filmSolveA(x, M, q, E, a); return E - slope(); };
    const Emax = Math.min(Ei, 1.02 * (M > 300 ? M : M / Math.tanh(M)) + 0.5);
    const E = N.brent(g, 1, Emax, { tol: 1e-11 * Emax });
    if (!keep) return { E };
    filmSolveA(x, M, q, E, a);
    const xs = Array.from(x), as = Array.from(a);
    return { E, x: xs, a: as, b: xs.map((z, i) => Math.max(0, 1 + (as[i] - E * (1 - z)) / q)) };
  };
  FF.E = function (theory, M, Ei) {
    if (theory === 'dec') return FF.Edec(M, Ei);
    if (theory === 'exact') return FF.filmExact(M, Ei).E;
    return FF.Evkh(M, Ei);
  };

  /* ---------- local rate (per volume of contactor), the general expression Eq. 23.5 ----------
     o: {kAga (∞ = no gas film), kAla, a, fl, H, DA, rB = D_B/D_A, b, k (second-order, ∞ = instantaneous), theory}
     CAb: dissolved A in the bulk (used for straight mass transfer, k = 0) */
  FF.local = function (o, pA, CB, CAb = 0) {
    const kAga = o.kAga, kAla = o.kAla, H = o.H, rB = o.rB, b = o.b, k2 = o.k, fl = o.fl;
    const kAl = kAla / o.a, renew = o.theory === 'dec';
    const res = { pA, CB };
    const gasF = isFinite(kAga);
    if (!(k2 > 0) || !(CB > 0)) {
      /* straight mass transfer */
      const R = (gasF ? 1 / kAga : 0) + H / kAla;
      res.r = Math.max(0, (pA - H * CAb) / R);
      res.E = 1; res.MH = 0; res.Ei = 1;
      res.pAi = gasF ? pA - res.r / kAga : pA; res.CAi = H > 0 ? res.pAi / H : NaN;
      res.R = [gasF ? 1 / kAga : 0, H / kAla, 0];
      res.share = res.R.map(q => q / R);
      res.case = 'physical';
      return res;
    }
    const MH = isFinite(k2) ? Math.sqrt(o.DA * k2 * CB) / kAl : Infinity;
    res.MH = MH;
    const EiOf = CAi => (renew ? Math.sqrt(1 / rB) + Math.sqrt(rB) * CB / (b * CAi) : 1 + rB * CB / (b * CAi));
    const ECof = CAi => {
      /* E·C_Ai, finite as C_Ai → 0 */
      if (!isFinite(MH)) return renew ? Math.sqrt(1 / rB) * CAi + Math.sqrt(rB) * CB / b : CAi + rB * CB / b;
      if (CAi <= 0) return 0;
      return FF.E(o.theory, MH, EiOf(CAi)) * CAi;
    };
    const bulkR = isFinite(k2) ? 1 / (k2 * CB * fl) : 0;
    const liq = CAi => { const ec = ECof(CAi); if (!(ec > 0)) return 0; const E = ec / Math.max(CAi, 1e-300); return CAi > 0 ? CAi / (1 / (kAla * E) + bulkR) : kAla * ec; };
    let CAi;
    if (H === 0) {
      res.r = gasF ? kAga * pA : Infinity; res.pAi = 0; res.CAi = NaN; res.E = NaN; res.Ei = Infinity;
      res.R = [1 / kAga, 0, 0]; res.share = [1, 0, 0]; res.case = 'gas film (very soluble gas)';
      return res;
    }
    if (!gasF) CAi = pA / H;
    else {
      const g = c => kAga * (pA - H * c) - liq(c);
      const g0 = g(0);
      if (g0 <= 0) CAi = 0;
      else CAi = N.brent(g, 0, pA / H, { fa: g0, tol: 1e-14 * pA / H + 1e-300 });
    }
    res.CAi = CAi; res.pAi = H * CAi;
    res.r = gasF ? kAga * (pA - H * CAi) : liq(CAi);
    res.Ei = CAi > 0 ? EiOf(CAi) : Infinity;
    res.E = CAi > 0 ? ECof(CAi) / CAi : Infinity;
    const Rg = gasF ? 1 / kAga : 0, Rl = CAi > 0 ? H / (kAla * res.E) : 0, Rb = H * bulkR;
    const tot = Rg + Rl + Rb || 1;
    res.R = [Rg, Rl, Rb]; res.share = [Rg / tot, Rl / tot, Rb / tot];
    /* regime, Fig. 23.3 */
    const kBla = kAla * (renew ? Math.sqrt(rB) : rB);
    if (!isFinite(MH)) res.case = CAi <= 0 ? 'B' : 'A';
    else if (MH > 2) res.case = CAi <= 0 ? 'B' : res.E > 0.9 * res.Ei ? 'A' : res.Ei > 5 * MH ? 'D' : 'C';
    else if (MH > 0.02) res.case = res.Ei > 5 * MH ? 'F' : 'E';
    else res.case = res.share[2] > 0.9 ? 'H' : 'G';
    res.kBla = kBla;
    return res;
  };
  FF.CASES = {
    A: 'instantaneous reaction on a plane in the liquid film (low C_B)',
    B: 'instantaneous reaction at the interface (high C_B): gas film controls',
    C: 'fast reaction in the liquid film, B depleted near the interface',
    D: 'fast pseudo-first-order reaction in the liquid film (high C_B)',
    E: 'intermediate rate: reaction in the film and the main body, B depleted',
    F: 'intermediate rate: reaction in the film and the main body, C_B high',
    G: 'slow reaction in the main body, film resistance still felt',
    H: 'very slow reaction in the main body, no mass-transfer resistance',
    physical: 'straight mass transfer (no reaction)'
  };

  /* common inputs for the fluid–fluid tools */
  FF.inputs = function (d = {}) {
    return [
      { group: 'Reaction A(g → l) + bB(l) → products, −r_A = k C_A C_B' },
      { k: 'inst', t: 'check', label: 'Instantaneous reaction (k → ∞)', v: !!d.inst },
      { k: 'k', label: 'Rate constant k', unit: 'm³/(mol·[t])', v: d.k !== undefined ? d.k : 1e6, min: 0, max: 1e20, if: v => !v.inst, hint: 'Zero for straight mass transfer.' },
      { k: 'b', label: 'b, mol B per mol A', v: d.b || 1, min: 0.01, max: 20 },
      { group: 'Mass transfer and solubility' },
      { k: 'gf', t: 'check', label: 'Gas-film resistance present (uncheck for pure A)', v: d.gf !== undefined ? d.gf : true },
      { k: 'kAga', label: 'k_Ag a', unit: 'mol/(m³·[t]·Pa)', v: d.kAga || 0.01, min: 1e-15, max: 1e12, scale: 'log', if: v => v.gf },
      { k: 'kAla', label: 'k_Al a', unit: '1/[t]', v: d.kAla || 20, min: 1e-15, max: 1e12, scale: 'log' },
      { k: 'a', label: 'Interfacial area a', unit: 'm²/m³', v: d.a || 20, min: 1e-6, max: 1e6, scale: 'log' },
      { k: 'fl', label: 'Liquid fraction f_l', v: d.fl || 0.98, min: 1e-4, max: 1 },
      { k: 'H', label: 'Henry constant H_A = p_A/C_A', unit: 'Pa·m³/mol', v: d.H !== undefined ? d.H : 1e5, min: 0, max: 1e12 },
      { k: 'DA', label: 'Diffusivity of A in the liquid D_A', unit: 'm²/[t]', v: d.DA || 1e-6, min: 1e-20, max: 1, scale: 'log' },
      { k: 'rB', label: 'D_B/D_A', v: d.rB || 1, min: 0.01, max: 100 },
      { k: 'theory', t: 'select', label: 'Enhancement factor from', opts: [['vkh', 'film theory, van Krevelen–Hoftijzer (Fig. 23.4)'], ['dec', 'surface renewal, Danckwerts–DeCoursey'], ['exact', 'film theory, exact numerical solution']], v: d.theory || 'vkh' }
    ];
  };
  FF.params = v => ({ kAga: v.gf ? +v.kAga : Infinity, kAla: +v.kAla, a: +v.a, fl: +v.fl, H: +v.H, DA: +v.DA, rB: +v.rB, b: +v.b, k: v.inst ? Infinity : +v.k, theory: v.theory });

  /* ===================================================================== */
  /* N-2301  rate of a G/L reaction                                         */
  /* ===================================================================== */
  ST.tool({
    id: 'glkin', tag: 'N-2301', ch: 23, sec: 'Ch. 23', title: 'Gas–liquid reaction rate: Hatta number, enhancement factor and regimes A–H',
    lede: 'The overall rate of absorption with reaction at any point of a contactor: gas film, liquid film with enhancement, and the main body of liquid in series (Eq. 23.5), with the interface conditions found self-consistently rather than guessed. Hatta modulus, instantaneous enhancement factor and E from Fig. 23.4 — or from surface-renewal theory, or from an exact numerical solution of the film equations — which of the eight regimes applies, where the resistance lies, how the concentrations run through both films, and how much reaction speeds up absorption.',
    keys: 'gas liquid absorption with chemical reaction two film theory hatta number enhancement factor instantaneous van krevelen hoftijzer danckwerts surface renewal decoursey penetration regimes A B C D E F G H reaction plane pseudo first order henry law gas film liquid film resistance example 23.1 problem 23.1 23.2 23.3 23.4 23.5 23.6 23.7 23.8 23.9 CO2 NaOH H2S MEA',
    units: { t: 'h' },
    inputs: [
      { group: 'Where in the contactor' },
      { k: 'pA', label: 'Partial pressure of A in the gas p_A', unit: 'Pa', v: 5000, min: 0, max: 1e9 },
      { k: 'CB', label: 'Concentration of B in the liquid C_B', unit: 'mol/m³', v: 100, min: 0, max: 1e7 },
      { k: 'mB', label: 'Order in B, m (−r = k C_A C_Bᵐ; book: treat as k C_Bᵐ⁻¹ · C_A C_B)', v: 1, min: 0, max: 4, step: 0.5 }
    ].concat(FF.inputs({ k: 1e6, b: 2, kAga: 0.01, kAla: 20, a: 20, fl: 0.98, H: 1e5, DA: 1e-6, rB: 1 })),
    presets: [
      { name: 'Example 23.1 — third-order reaction, air bubbling through B', v: { pA: 5000, CB: 100, mB: 2, inst: false, k: 1e6, b: 2, gf: true, kAga: 0.01, kAla: 20, a: 20, fl: 0.98, H: 1e5, DA: 1e-6, rB: 1, theory: 'vkh', _tu: 'h' }, book: { r: { v: 33, tol: 0.02 }, MH: { v: 100, tol: 0.001 }, shG: { v: 2 / 3, tol: 0.02 }, regime: 'D' }, note: '−r_A = k C_A C_B², k = 10⁶ m⁶/mol²·hr, treated as second order with k C_B. Book: M_H = 100, E = M_H, two-thirds of the resistance in the gas film, the reaction is pseudo-first order in the film at the interface, and −r‴ = 33 mol/hr·m³.' },
      { name: 'Problem 23.1 — k = 10, H = 10⁵', v: { pA: 100, CB: 100, mB: 1, inst: false, k: 10, b: 1, gf: true, kAga: 0.1, kAla: 100, a: 100, fl: 0.01, H: 1e5, DA: 1e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 23.2 — k = 10⁶, H = 10⁴', v: { pA: 100, CB: 100, mB: 1, inst: false, k: 1e6, b: 1, gf: true, kAga: 0.1, kAla: 100, a: 100, fl: 0.01, H: 1e4, DA: 1e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 23.3 — k = 10, H = 10³', v: { pA: 100, CB: 100, mB: 1, inst: false, k: 10, b: 1, gf: true, kAga: 0.1, kAla: 100, a: 100, fl: 0.01, H: 1e3, DA: 1e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 23.4 — k = 10⁻⁴, H = 1', v: { pA: 100, CB: 100, mB: 1, inst: false, k: 1e-4, b: 1, gf: true, kAga: 0.1, kAla: 100, a: 100, fl: 0.01, H: 1, DA: 1e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 23.5 — k = 10⁻², H = 1', v: { pA: 100, CB: 100, mB: 1, inst: false, k: 1e-2, b: 1, gf: true, kAga: 0.1, kAla: 100, a: 100, fl: 0.01, H: 1, DA: 1e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 23.6 — k = 10⁸, H = 1', v: { pA: 100, CB: 100, mB: 1, inst: false, k: 1e8, b: 1, gf: true, kAga: 0.1, kAla: 100, a: 100, fl: 0.01, H: 1, DA: 1e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 23.7 — Example 23.1 with C_B = 1', v: { pA: 5000, CB: 1, mB: 2, inst: false, k: 1e6, b: 2, gf: true, kAga: 0.01, kAla: 20, a: 20, fl: 0.98, H: 1e5, DA: 1e-6, rB: 1, theory: 'vkh', _tu: 'h' } },
      { name: 'Problem 23.8 — CO₂ into caustic soda', v: { pA: 1e5, CB: 500, mB: 1, inst: false, k: 10, b: 2, gf: true, kAga: 1e-2, kAla: 1e-2, a: 100, fl: 0.1, H: 25000, DA: 1.8e-9, rB: 3.06 / 1.8, theory: 'vkh', _tu: 's' }, note: 'Danckwerts’ data: k_Ag = 10⁻⁴ mol/m²·s·Pa and k_Al = 10⁻⁴ m/s with a = 100 m⁻¹ (the problem prints k_Ag a; the per-area reading is used), H_A = 25 000 Pa·m³/mol, D_A = 1.8 × 10⁻⁹, D_B = 3.06 × 10⁻⁹ m²/s, k = 10 m³/mol·s, f_l = 0.1, b = 2.' },
      { name: 'Problem 23.9 — H₂S into methanolamine', v: { pA: 0.02 * 101325, CB: 250, mB: 1, inst: true, b: 1, gf: true, kAga: 60 / 101325, kAla: 0.03, a: 100, fl: 0.1, H: 1e-4 * 101325, DA: 1e-9, rB: 0.64, theory: 'vkh', _tu: 's' }, note: '20 atm gas with 0.1% H₂S, 250 mol/m³ of free MEA with D_MEA = 0.64 D_H₂S, instantaneous reaction, k_Al a = 0.03 s⁻¹, k_Ag a = 60 mol/m³·s·atm, H_A = 10⁻⁴ m³·atm/mol. The comparison with pure water is the “versus straight mass transfer” line.' }
    ],
    eqs: [
      ['-r_A^{\\prime\\prime\\prime} = \\cfrac{p_A}{\\cfrac{1}{k_{Ag}a}+\\cfrac{H_A}{k_{Al}aE}+\\cfrac{H_A}{kC_Bf_l}}', 'Eq. 23.5', 'Gas film, liquid film, main body of liquid in series.'],
      ['M_H = \\frac{\\sqrt{\\mathscr{D}_{Al}kC_B}}{k_{Al}},\\qquad E_i = 1 + \\frac{\\mathscr{D}_BC_BH_A}{b\\mathscr{D}_Ap_{Ai}}', 'Eqs. 23.7–23.8, Fig. 23.4'],
      ['E = \\frac{M_H\\sqrt{(E_i-E)/(E_i-1)}}{\\tanh\\left(M_H\\sqrt{(E_i-E)/(E_i-1)}\\right)}', 'van Krevelen–Hoftijzer', 'The curves of Fig. 23.4.'],
      ['E = -\\frac{M_H^2}{2(E_i-1)} + \\sqrt{\\frac{M_H^4}{4(E_i-1)^2}+\\frac{E_iM_H^2}{E_i-1}+1}', 'DeCoursey', 'Surface renewal; E_i = √(D_A/D_B) + √(D_B/D_A) C_B/(bC_Ai).'],
      ['-r_A^{\\prime\\prime} = \\frac{\\frac{\\mathscr{D}_B}{\\mathscr{D}_A}\\frac{C_B}{b}+\\frac{p_A}{H_A}}{\\frac{1}{H_Ak_{Ag}}+\\frac{1}{k_{Al}}}\\ (\\text{A}),\\qquad -r_A^{\\prime\\prime} = k_{Ag}p_A\\ (\\text{B})', 'Eqs. 23.13, 23.16', 'Instantaneous: B when k_Ag p_A ≤ k_Bl C_B/b (Eq. 23.17).'],
      ['-r_A^{\\prime\\prime\\prime} = \\cfrac{p_A}{\\cfrac{1}{k_{Ag}a}+\\cfrac{H_A}{a\\sqrt{\\mathscr{D}_AkC_B}}}', 'Eq. 23.19', 'Case D, pseudo first order.'],
      ['-r_A^{\\prime\\prime\\prime} = \\frac{kf_l}{H_A}p_AC_B', 'Eq. 23.21', 'Case H.']
    ],
    method: ['The interface concentration is found by Brent’s method from k_Ag a(p_A − H_A C_Ai) = C_Ai/[1/(k_Al aE) + 1/(kC_Bf_l)], with E evaluated at the trial C_Ai because E_i depends on it — the guess-and-check of Example 23.1 done exactly. When the gas film can deliver no more than k_Bl C_B/b, the reaction plane sits at the interface (Case B).', 'The exact film solution uses the invariant D_B C_B − bD_A C_A = linear in x to reduce the two film equations to one, a″ = M² a[1 + (a − E(1 − x))/(E_i − 1)], solved by Newton’s method on a grid crowded at the interface and around the reaction plane, with E found by Brent so that E = −a′(0). It reproduces M/tanh M and E_i in its two limits; van Krevelen–Hoftijzer’s formula is within a few per cent of it.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'h';
      const o = FF.params(v);
      const mB = +v.mB, CB = +v.CB, pA = +v.pA;
      if (!v.inst && mB !== 1 && CB > 0) o.k = +v.k * Math.pow(CB, mB - 1);
      const R = FF.local(o, pA, CB, 0);
      const ru = 'mol/(m³·' + tu + ')';
      kpi.push({ k: 'r', l: 'Rate of reaction −r‴_A', v: R.r, u: ru, main: true });
      kpi.push({ k: 'regime', l: 'Regime (Fig. 23.3)', v: R.case, s: FF.CASES[R.case] || '' });
      kpi.push({ k: 'MH', l: 'Hatta modulus M_H', v: R.MH, s: R.MH > 2 ? 'reaction in the film' : R.MH > 0.02 ? 'film and main body' : 'main body only' });
      kpi.push({ k: 'Ei', l: 'Enhancement factor for instantaneous reaction E_i', v: R.Ei });
      kpi.push({ k: 'E', l: 'Enhancement factor E', v: R.E });
      kpi.push({ k: 'pAi', l: 'p_A at the interface', v: R.pAi, u: 'Pa', s: 'C_Ai = ' + fmt(R.CAi, 4) + ' mol/m³' });
      ['gas film', 'liquid film', 'main body of liquid'].forEach((nm, i) => kpi.push({ k: ['shG', 'shL', 'shB'][i], l: 'Resistance in the ' + nm, v: R.share[i], u: '' }));
      /* physical absorption for comparison */
      const P = FF.local(Object.assign({}, o, { k: 0 }), pA, CB, 0);
      kpi.push({ k: 'boost', l: 'Rate relative to physical absorption into fresh, A-free liquid', v: R.r / P.r, s: 'physical absorption: ' + fmt(P.r, 4) + ' ' + ru + (R.r < P.r ? '; below 1: the main body cannot use A as fast as fresh liquid would take it up' : '') });
      /* Case A/B test, Eq. 23.17 */
      if (v.gf) {
        const lhs = o.kAga * pA, rhs = R.kBla ? R.kBla * CB / o.b : NaN;
        if (isFinite(rhs)) info.push('Test of Eq. 23.17: k_Ag a p_A = ' + fmt(lhs, 4) + ' against k_Bl a C_B/b = ' + fmt(rhs, 4) + ' — for an instantaneous reaction this is Case ' + (lhs > rhs ? 'A (plane inside the film)' : 'B (plane at the interface)') + '.');
      }
      if (R.case === 'D' || R.case === 'B') info.push('Reaction is confined to the film next to the interface: a contactor with large interfacial area (spray or plate column, agitated tank) is favoured; liquid holdup matters little.');
      if (R.case === 'G' || R.case === 'H') info.push('Reaction happens in the main body of liquid: a contactor with a large liquid holdup (bubble column, tank) is favoured; interfacial area matters little.');
      /* compare theories */
      if (isFinite(R.MH) && isFinite(R.Ei) && R.MH > 0) {
        const th = [['film, van Krevelen–Hoftijzer', 'vkh'], ['surface renewal, DeCoursey', 'dec'], ['film, exact', 'exact']];
        tables.push({ title: 'The enhancement factor at this point by each theory', cols: ['theory', 'E (at this E_i)', 'rate, self-consistent'], rows: th.map(([nm, t]) => { const Rt = FF.local(Object.assign({}, o, { theory: t }), pA, CB, 0); return [nm, t === 'dec' ? FF.Edec(R.MH, R.Ei) : t === 'exact' ? FF.filmExact(R.MH, R.Ei).E : FF.Evkh(R.MH, R.Ei), Rt.r]; }) });
      }
      /* Fig. 23.4 */
      const Ms = N.logspace(0.1, 1000, 120), ser = [];
      [2, 5, 10, 20, 50, 100, 200, 500, 1000].forEach((Ei, i) => ser.push({ name: 'E_i = ' + Ei, x: Ms, y: Ms.map(M => FF.Evkh(M, Ei)), color: i % 2 ? 'var(--muted)' : 3, w: 1, legend: i === 0 || i === 8 ? true : false }));
      ser.push({ name: 'E_i = ∞ (M_H/tanh M_H)', x: Ms, y: Ms.map(M => M / Math.tanh(M)), color: 2, w: 1.6 });
      if (isFinite(R.Ei) && isFinite(R.MH)) {
        ser.push({ name: 'this E_i, van Krevelen–Hoftijzer', x: Ms, y: Ms.map(M => FF.Evkh(M, R.Ei)), color: 1, w: 2.4 });
        ser.push({ name: 'this E_i, surface renewal', x: Ms, y: Ms.map(M => FF.Edec(M, R.Ei)), color: 5, dash: '5 4', w: 1.6 });
        const Mx = N.logspace(0.1, 1000, 26);
        ser.push({ name: 'this E_i, exact film theory', x: Mx, y: Mx.map(M => FF.filmExact(M, R.Ei).E), type: 'pts', color: 4, r: 3 });
      }
      const marks = [{ t: 'v', x: 2, color: 'var(--muted)', label: 'M_H = 2' }, { t: 'v', x: 0.3, color: 'var(--muted)', label: '' }];
      if (isFinite(R.MH) && R.MH >= 0.1 && R.MH <= 1000 && isFinite(R.E)) marks.push({ t: 'pt', x: R.MH, y: R.E, color: 8, label: 'here' });
      plots.push({ title: 'Enhancement factor against Hatta modulus (Fig. 23.4)', x: { label: 'M_H', log: true, min: 0.1, max: 1000 }, y: { label: 'E', log: true, min: 1, max: 1000 }, series: ser, marks, hmax: 470, cap: 'Below M_H ≈ 0.3 reaction in the film is negligible (E → 1); for E_i > 5M_H the reaction is pseudo first order (E → M_H); for M_H > 5E_i it is instantaneous (E → E_i).' });
      /* concentration profiles through both films */
      if (R.CAi >= 0 && isFinite(o.H) && o.H > 0) {
        const CAstar = pA / o.H;
        const xs = [], ya = [], yb = [];
        /* gas film: linear from p_A to p_Ai, drawn on x ∈ [−1, 0] */
        [-1.3, -1].forEach(x => { xs.push(x); ya.push(1); yb.push(NaN); });
        xs.push(0); ya.push(R.CAi / CAstar); yb.push(NaN);
        if (isFinite(R.MH) && isFinite(R.Ei) && R.MH > 0) {
          const F = o.theory === 'exact' || true ? FF.filmExact(R.MH, R.Ei, true) : null;
          if (F && F.x) F.x.forEach((x, i) => { if (i % Math.max(1, Math.floor(F.x.length / 300)) === 0 || i === F.x.length - 1) { xs.push(x); ya.push(F.a[i] * R.CAi / CAstar); yb.push(F.b[i]); } });
        } else if (!isFinite(R.MH) && R.CAi > 0) {
          const xr = 1 / R.Ei;
          [[0, R.CAi / CAstar, 0], [xr, 0, 0], [1, 0, 1]].forEach(p => { xs.push(p[0]); ya.push(p[1]); yb.push(p[2]); });
        }
        xs.push(1.3); ya.push(0); yb.push(1);
        const Ab = xs.map((x, i) => (x <= 0 ? ya[i] : NaN)), Al = xs.map((x, i) => (x >= 0 ? ya[i] : NaN));
        plots.push({ title: 'Through the films: A (as C_A/C_A*) and B (as C_B/C_B,bulk)', x: { label: 'gas film ← interface → liquid film (x/x₀)', min: -1.3, max: 1.3 }, y: { label: 'relative concentration', min: 0, max: 1.05 }, series: [{ name: 'A, gas side (p_A/p_A,bulk)', x: xs, y: Ab, color: 2, w: 2.2 }, { name: 'A in the liquid film', x: xs, y: Al, color: 1, w: 2.2 }, { name: 'B in the liquid film', x: xs, y: yb, color: 3, w: 2.2 }], marks: [{ t: 'band', x0: -1, x1: 0, color: 'var(--grid)' }, { t: 'band', x0: 0, x1: 1, color: 'var(--grid)' }, { t: 'v', x: 0, color: 'var(--muted)', label: 'interface' }], cap: 'Exact film-theory profiles at this M_H and E_i. C_A* = p_A/H_A, the liquid concentration in equilibrium with the bulk gas.' });
      }
      /* rate against C_B */
      const Cg = N.logspace(Math.max(CB, 1e-6) * 1e-3, Math.max(CB, 1e-6) * 1e3, 90);
      const rk = Cg.map(c => { const oo = Object.assign({}, o); if (!v.inst && mB !== 1) oo.k = +v.k * Math.pow(c, mB - 1); return FF.local(oo, pA, c, 0); });
      plots.push({ title: 'Rate against the concentration of B', x: { label: 'C_B', unit: 'mol/m³', log: true }, y: { label: '−r‴_A', unit: ru, log: true }, series: [{ name: 'with reaction', x: Cg, y: rk.map(q => q.r), color: 1, w: 2.2 }, { name: 'straight mass transfer', x: Cg, y: Cg.map(() => P.r), color: 'var(--muted)', dash: '5 4' }, { name: 'gas film alone (k_Ag a p_A)', x: Cg, y: Cg.map(() => (v.gf ? o.kAga * pA : NaN)), color: 2, dash: '2 3' }], marks: [{ t: 'pt', x: CB, y: R.r, color: 8, label: R.case }], cap: 'Raising C_B moves the reaction toward the interface until the gas film alone limits the rate.' });
      plots.push({ title: 'Where the resistance lies', x: { cats: ['gas film', 'liquid film', 'main body'] }, y: { label: 'share', min: 0, max: 1.1 }, series: [{ name: 'share', x: [0, 1, 2], y: R.share, type: 'bar', color: 2, legend: false }], marks: R.share.map((q, i) => ({ t: 'text', x: i, y: q + 0.035, s: N.fmt(q, 3), anchor: 'middle' })) });
      const casesRows = 'ABCDEFGH'.split('').map(c => [c === R.case ? '▶ ' + c : c, FF.CASES[c]]);
      tables.push({ title: 'The eight regimes', cols: ['case', 'behaviour'], rows: casesRows });
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
