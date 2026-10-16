/* Chapter 20 — fluidized and other suspended-solids reactors: regimes and velocities, the
   Kunii–Levenspiel bubbling bed (single, reversible, series and parallel reactions), and the
   solids distribution and conversion of turbulent, fast-fluidized and pneumatic-conveying beds */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N;
  const fmt = (x, s = 4) => N.fmt(x, s);
  const g0 = 9.80665;
  const FZ = ST.FZ = {};

  /* dimensionless size and velocity (Eqs. 20.1–20.2) */
  FZ.dstar = (dp, rg, rs, mu) => dp * Math.cbrt(rg * (rs - rg) * g0 / (mu * mu));
  FZ.ustar = (u, rg, rs, mu) => u * Math.cbrt(rg * rg / (mu * (rs - rg) * g0));
  FZ.uOf = (us, rg, rs, mu) => us / Math.cbrt(rg * rg / (mu * (rs - rg) * g0));
  /* minimum fluidization, Ergun (Eq. 20.3): 150(1−ε)u* + 1.75 u*² d* = ε³ d*², in u* */
  FZ.umfStar = (ds, eps) => { const a = 1.75 * ds, b = 150 * (1 - eps), c = -eps * eps * eps * ds * ds; return (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a); };
  /* Wen and Yu: Re_mf = sqrt(33.7² + 0.0408 Ar) − 33.7, Ar = d*³ */
  FZ.umfStarWY = ds => (Math.sqrt(33.7 * 33.7 + 0.0408 * ds * ds * ds) - 33.7) / ds;
  /* terminal velocity, Haider and Levenspiel (Eqs. 20.4–20.5) */
  FZ.utStar = (ds, phi) => (phi >= 0.9999 ? 1 / (18 / (ds * ds) + 0.591 / Math.sqrt(ds)) : 1 / (18 / (ds * ds) + (2.335 - 1.744 * phi) / Math.sqrt(ds)));

  /* ===================================================================== */
  /* C-2001  fluidization regimes, velocities, Geldart group                */
  /* ===================================================================== */
  ST.tool({
    id: 'fluidize', tag: 'C-2001', ch: 20, sec: 'Ch. 20.1', title: 'Fluidization regimes: minimum fluidization, terminal velocity, Geldart group',
    lede: 'Where a gas–solid system sits among fixed, bubbling, turbulent, fast fluidized and pneumatic-conveying contacting: d*_p and u*, the minimum fluidizing velocity by Ergun, the terminal velocity of spheres and irregular particles, the Geldart class, and the transitions to fast fluidization and to pneumatic transport, on the general regime chart and the Geldart chart.',
    keys: 'fluidization regime minimum fluidizing velocity ergun terminal velocity haider levenspiel sphericity geldart classification A B C D bubbling turbulent fast fluidized pneumatic conveying choking dimensionless particle size velocity grace chart circulating fluidized bed',
    inputs: [
      { group: 'Particles and gas' },
      { k: 'dp', label: 'Particle size d_p', unit: 'μm', v: 100, min: 1, max: 1e5, scale: 'log' },
      { k: 'phi', label: 'Sphericity φ_s', v: 1, min: 0.1, max: 1 },
      { k: 'rhos', label: 'Particle density ρ_s', unit: 'kg/m³', v: 2000, min: 10, max: 2e4, scale: 'log' },
      { k: 'rhog', label: 'Gas density ρ_g', unit: 'kg/m³', v: 1.2, min: 0.01, max: 300, scale: 'log' },
      { k: 'mu', label: 'Gas viscosity μ', unit: 'Pa·s', v: 1.8e-5, min: 1e-7, max: 1e-2, scale: 'log' },
      { k: 'emf', label: 'Voidage at minimum fluidization ε_mf', v: 0.5, min: 0.3, max: 0.8 },
      { k: 'umfM', t: 'seg', label: 'u_mf from', opts: [['ergun', 'Ergun with ε_mf (Eq. 20.3)'], ['wy', 'Wen and Yu']], v: 'ergun' },
      { group: 'Operation' },
      { k: 'u0', label: 'Superficial gas velocity u₀', unit: 'm/s', v: 0.3, min: 1e-5, max: 100, scale: 'log' },
      { k: 'Gs', label: 'Solids circulation flux G_s', unit: 'kg/(m²·s)', v: 50, min: 0, max: 5000 }
    ],
    presets: [
      { name: 'FCC catalyst (Geldart A)', v: { dp: 60, phi: 0.9, rhos: 1500, rhog: 1.2, mu: 1.8e-5, emf: 0.45, u0: 0.3 } },
      { name: 'Sand in air (Geldart B)', v: { dp: 300, phi: 0.86, rhos: 2600, rhog: 1.2, mu: 1.8e-5, emf: 0.45, u0: 0.5 } },
      { name: 'Coarse grain (Geldart D, spoutable)', v: { dp: 2000, phi: 0.9, rhos: 1200, rhog: 1.2, mu: 1.8e-5, emf: 0.42, u0: 1.5 } },
      { name: 'Example 20.1 solids at 0.3 m/s', v: { dp: 68, phi: 1, rhos: 2000, rhog: 1.2, mu: 1.8e-5, emf: 0.5, u0: 0.3 } },
      { name: 'Riser of a cat cracker', v: { dp: 70, phi: 0.9, rhos: 1500, rhog: 1.5, mu: 2.5e-5, emf: 0.45, u0: 8, Gs: 400 } }
    ],
    eqs: [
      ['d_p^* = d_p\\left[\\frac{\\rho_g(\\rho_s-\\rho_g)g}{\\mu^2}\\right]^{1/3},\\qquad u^* = u\\left[\\frac{\\rho_g^2}{\\mu(\\rho_s-\\rho_g)g}\\right]^{1/3}', 'Eqs. 20.1–20.2'],
      ['150(1-\\varepsilon_{mf})u_{mf}^* + 1.75(u_{mf}^*)^2d_p^* = \\varepsilon_{mf}^3(d_p^*)^2', 'Eq. 20.3', 'Ergun at minimum fluidization.'],
      ['u_t^* = \\left[\\frac{18}{(d_p^*)^2}+\\frac{0.591}{(d_p^*)^{1/2}}\\right]^{-1}\\ \\text{(spheres)}', 'Eq. 20.4'],
      ['u_t^* = \\left[\\frac{18}{(d_p^*)^2}+\\frac{2.335-1.744\\phi_s}{(d_p^*)^{1/2}}\\right]^{-1}', 'Eq. 20.5', 'Irregular particles, 0.5 < φ_s < 1.'],
      ['u_{TB\\text{-}FF} = 1.53\\sqrt{\\frac{(\\rho_s-\\rho_g)gd_p}{\\rho_g}}', 'Bi et al.'],
      ['\\frac{u_{FF\\text{-}PC}}{\\sqrt{gd_p}} = 21.6\\left(\\frac{G_s}{\\rho_g u_{FF\\text{-}PC}}\\right)^{0.542}(d_p^*)^{0.315}', 'Bi and Fan'],
      ['d_{p,A/B}^* = 101\\left(\\frac{\\rho_g}{\\rho_s-\\rho_g}\\right)^{0.425},\\qquad d_{p,C/A}^* \\approx 0.68\\text{–}1.1', 'Fig. 20.5']
    ],
    method: ['All velocities are computed in the dimensionless form of Eqs. 20.1–20.5 and converted back. The fast-fluidization–pneumatic transition is implicit in u and is solved by Brent’s method. The onset of turbulent fluidization uses Bi and Grace’s Re_c = 0.565 Ar^{0.461}, one of the correlations behind the band of Fig. 20.5; the A/B boundary is the line of Fig. 20.5, C/A the middle of its band, and B/D Geldart’s (ρ_s − ρ_g)d_p² ≈ 10⁶ (g/cm³)(μm²).'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const dp = +v.dp * 1e-6, rs = +v.rhos, rg = +v.rhog, mu = +v.mu, phi = +v.phi, emf = +v.emf, u0 = +v.u0, Gs = +v.Gs;
      if (!(rs > rg)) throw new Error('The particles must be denser than the gas.');
      const ds = FZ.dstar(dp, rg, rs, mu), us0 = FZ.ustar(u0, rg, rs, mu);
      const umfS = v.umfM === 'wy' ? FZ.umfStarWY(ds) : FZ.umfStar(ds, emf), umf = FZ.uOf(umfS, rg, rs, mu);
      const utS = FZ.utStar(ds, phi), ut = FZ.uOf(utS, rg, rs, mu);
      const uTF = 1.53 * Math.sqrt((rs - rg) * g0 * dp / rg);
      const Ar = ds * ds * ds, uc = 0.565 * Math.pow(Ar, 0.461) * mu / (rg * dp);
      const fPC = u => u / Math.sqrt(g0 * dp) - 21.6 * Math.pow(Gs / (rg * u), 0.542) * Math.pow(ds, 0.315);
      const uPC = Gs > 0 ? N.brent(fPC, 1e-4, 1e3, { tol: 1e-12 }) : NaN;
      const dsAB = 101 * Math.pow(rg / (rs - rg), 0.425), dsCA = 0.9;
      const dpum = +v.dp, drho = (rs - rg) / 1000;
      const geld = ds < dsCA ? 'C (cohesive)' : ds < dsAB ? 'A (aeratable)' : drho * dpum * dpum < 1e6 ? 'B (sand-like)' : 'D (spoutable)';
      let regime;
      if (u0 < umf) regime = 'fixed bed';
      else if (geld[0] === 'A' && u0 < 3 * umf) regime = 'expanded bed (bubble-free, Geldart A)';
      else if (u0 < uc) regime = 'bubbling fluidized bed (BFB)';
      else if (u0 < uTF) regime = 'turbulent fluidized bed (TB)';
      else if (!(Gs > 0) || u0 < uPC) regime = 'fast fluidized bed (FF)';
      else regime = 'pneumatic conveying (PC)';
      kpi.push({ k: 'regime', l: 'Contacting regime', v: regime, main: true, s: 'u₀ = ' + fmt(u0 / umf, 3) + ' u_mf' });
      kpi.push({ k: 'geldart', l: 'Geldart group', v: geld });
      kpi.push({ k: 'umf', l: 'Minimum fluidizing velocity u_mf', v: umf, u: 'm/s', s: 'u*_mf = ' + fmt(umfS, 4) });
      kpi.push({ k: 'ut', l: 'Terminal velocity of one particle u_t', v: ut, u: 'm/s', s: 'u_t/u_mf = ' + fmt(ut / umf, 3) });
      kpi.push({ k: 'uc', l: 'Onset of turbulent fluidization u_c', v: uc, u: 'm/s', s: 'Bi and Grace' });
      kpi.push({ k: 'uTF', l: 'Turbulent → fast fluidized u_TB-FF', v: uTF, u: 'm/s' });
      if (Gs > 0) kpi.push({ k: 'uPC', l: 'Fast fluidized → pneumatic at this G_s', v: uPC, u: 'm/s' });
      kpi.push({ k: 'ds', l: 'd*_p and u*', v: ds, s: 'u* = ' + fmt(us0, 4) + ', Re_p = ' + fmt(us0 * ds, 4) });
      if (geld[0] === 'C') warn.push('Group C powders are cohesive: they channel and lift as plugs rather than fluidize.');
      if (geld[0] === 'D') info.push('Group D solids form large exploding bubbles and can spout; a spouted bed may suit them better.');
      /* regime chart (Fig. 20.5) */
      const dg = N.logspace(0.5, 300, 160);
      const ser = [
        { name: 'u*_mf (Eq. 20.3)', x: dg, y: dg.map(d => (v.umfM === 'wy' ? FZ.umfStarWY(d) : FZ.umfStar(d, emf))), color: 2, w: 2 },
        { name: 'u*_t, spheres (Eq. 20.4)', x: dg, y: dg.map(d => FZ.utStar(d, 1)), color: 1, w: 2 },
        { name: 'onset of turbulence (Bi and Grace)', x: dg, y: dg.map(d => 0.565 * Math.pow(d * d * d, 0.461) / d), color: 3, dash: '5 4' },
        { name: 'u*_TB-FF', x: dg, y: dg.map(d => FZ.ustar(1.53 * Math.sqrt((rs - rg) * g0 * (d / Math.cbrt(rg * (rs - rg) * g0 / (mu * mu))) / rg), rg, rs, mu)), color: 4, dash: '2 3' }
      ];
      if (phi < 0.9999) ser.splice(2, 0, { name: 'u*_t, φ_s = ' + fmt(phi, 3), x: dg, y: dg.map(d => FZ.utStar(d, phi)), color: 1, dash: '6 3' });
      plots.push({ title: 'General flow regime diagram (Fig. 20.5)', x: { label: 'd*_p', log: true, min: 0.5, max: 300 }, y: { label: 'u*', log: true, min: 1e-3, max: 100 }, series: ser, marks: [{ t: 'v', x: dsAB, color: 'var(--muted)', label: 'A | B' }, { t: 'v', x: dsCA, color: 'var(--muted)', label: 'C | A' }, { t: 'pt', x: ds, y: us0, color: 8, label: 'operating point' }], cap: 'Below u*_mf: fixed bed. Between u*_mf and the onset of turbulence: bubbling. Above u*_t particles blow out singly; whole beds blow out above u*_TB-FF.' });
      /* Geldart chart (Fig. 20.2) */
      const dpg = N.logspace(10, 2000, 120);
      const dAB = N.logspace(10, 2000, 120).map(d => d);
      const ABline = dAB.map(d => { const dr = 1; return 225 / d; });
      plots.push({ title: 'Geldart classification (Fig. 20.2)', x: { label: 'd_p', unit: 'μm', log: true, min: 10, max: 2000 }, y: { label: 'ρ_s − ρ_g', unit: 'g/cm³', log: true, min: 0.1, max: 10 }, series: [{ name: 'A | B: (ρ_s − ρ_g) d_p ≈ 225', x: dAB, y: ABline, color: 3, w: 2 }, { name: 'B | D: (ρ_s − ρ_g) d_p² ≈ 10⁶', x: dpg, y: dpg.map(d => 1e6 / (d * d)), color: 4, w: 2 }, { name: 'C | A (approximate)', x: dpg, y: dpg.map(d => Math.pow(10, 1.6 - 2.2 * Math.log10(d / 20))), color: 2, dash: '5 4' }], marks: [{ t: 'pt', x: dpum, y: drho, color: 8, label: 'these solids' }, { t: 'text', x: 18, y: 0.2, s: 'C' }, { t: 'text', x: 60, y: 0.6, s: 'A' }, { t: 'text', x: 400, y: 1.4, s: 'B' }, { t: 'text', x: 1200, y: 4, s: 'D' }], cap: 'The A | B line here is Geldart’s; the C | A boundary is only indicative.' });
      tables.push({ title: 'Velocities', cols: ['quantity', 'm/s', 'u*'], rows: [['u_mf', umf, umfS], ['u_c (turbulence)', uc, FZ.ustar(uc, rg, rs, mu)], ['u_t (one particle)', ut, utS], ['u_TB-FF', uTF, FZ.ustar(uTF, rg, rs, mu)]].concat(Gs > 0 ? [['u_FF-PC at G_s', uPC, FZ.ustar(uPC, rg, rs, mu)]] : []).concat([['u₀ (operating)', u0, us0]]) });
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* C-2002  the Kunii–Levenspiel bubbling bed                              */
  /* ===================================================================== */
  /* bubble-phase quantities for a K–L bed */
  FZ.kl = function (p) {
    const ubr = 0.711 * Math.sqrt(g0 * p.db);
    const ub = p.u0 - p.umf + ubr;
    const delta = p.deltaExact ? (p.u0 - p.umf) / (ub - p.umf) : p.u0 / ub;
    const ef = 1 - (1 - p.emf) * (1 - delta);
    const Kbc = 4.5 * p.umf / p.db + 5.85 * Math.sqrt(p.D) * Math.pow(g0, 0.25) / Math.pow(p.db, 1.25);
    const Kce = 6.77 * Math.sqrt(p.emf * p.D * ubr / Math.pow(p.db, 3));
    const fb = p.fb;
    const fc = delta * (1 - p.emf) * (3 * p.umf / p.emf / (ubr - p.umf / p.emf) + p.alpha);
    const fe = (1 - p.emf) * (1 - delta) - fc - fb;
    return { ubr, ub, delta, ef, Kbc, Kce, fb, fc, fe, ftot: 1 - ef, fast: ub / (p.umf / p.emf) };
  };
  /* first order: bracket of Eq. 20.21 (per volume of bed), and the five resistances */
  FZ.kf = (q, k) => q.fb * k + 1 / (1 / (q.delta * q.Kbc) + 1 / (q.fc * k + 1 / (1 / (q.delta * q.Kce) + 1 / (q.fe * k))));
  /* linear kinetics with a rate matrix Km (consumption = Km·C): bubble ODE matrix */
  const mat = {
    add: (A, B) => A.map((r, i) => r.map((x, j) => x + B[i][j])),
    sc: (A, s) => A.map(r => r.map(x => x * s)),
    eye: n => N.range(n).map(i => N.range(n).map(j => (i === j ? 1 : 0))),
    mul: (A, B) => A.map(r => B[0].map((_, j) => r.reduce((t, x, k) => t + x * B[k][j], 0))),
    inv: A => N.inv(A)
  };
  FZ.bubbleMatrix = function (q, Km) {
    const n = Km.length, I = mat.eye(n);
    const Me = mat.mul(mat.inv(mat.add(mat.sc(I, q.delta * q.Kce), mat.sc(Km, q.fe))), mat.sc(I, q.delta * q.Kce));
    const Lc = mat.add(mat.add(mat.sc(I, q.delta * q.Kbc), mat.sc(Km, q.fc)), mat.sc(mat.add(I, mat.sc(Me, -1)), q.delta * q.Kce));
    const Mc = mat.mul(mat.inv(Lc), mat.sc(I, q.delta * q.Kbc));
    const A = mat.add(mat.sc(Km, q.fb), mat.sc(mat.add(I, mat.sc(Mc, -1)), q.delta * q.Kbc));
    return { A, Mc, Me: mat.mul(Me, Mc) };
  };
  /* C(H) for dC/dz = −(A/u0) C by small-step matrix exponential */
  FZ.propagate = function (A, u0, C0, H) {
    const n = A.length, Z = mat.sc(A, -H / u0);
    let s = 0; const nrm = Math.max(...Z.map(r => r.reduce((t, x) => t + Math.abs(x), 0)));
    while (nrm / Math.pow(2, s) > 0.5) s++;
    const B = mat.sc(Z, 1 / Math.pow(2, s));
    let E = mat.eye(n), term = mat.eye(n);
    for (let k = 1; k < 20; k++) { term = mat.sc(mat.mul(term, B), 1 / k); E = mat.add(E, term); }
    for (let i = 0; i < s; i++) E = mat.mul(E, E);
    return E.map(r => r.reduce((t, x, j) => t + x * C0[j], 0));
  };

  ST.tool({
    id: 'klbfb', tag: 'C-2002', ch: 20, sec: 'Ch. 20.2–20.3', title: 'The Kunii–Levenspiel bubbling fluidized bed',
    lede: 'The bubbling bed as bubble, cloud and emulsion with one adjustable parameter, the bubble size. All flow quantities, the five resistances in series–parallel, the conversion, the concentration the solids actually see, the catalyst needed for a target, and the comparison with plug and mixed flow — for first-order, reversible, series A → R → S and parallel reactions, with bubble size as the lever.',
    keys: 'kunii levenspiel bubbling fluidized bed bfb bubble cloud emulsion wake bubble size rise velocity davidson interchange coefficient Kbc Kce fast bubble first order conversion series reaction intermediate maximum fig 20.10 20.11 20.12 example 20.1 problem 20.1 20.2 20.3 20.4 20.5 20.6 20.7 20.8 20.9 20.10 20.11 baffles catalyst mean concentration seen by solids',
    inputs: [
      { group: 'Bed and gas' },
      { k: 'amt', t: 'seg', label: 'Amount of solids given as', opts: [['W', 'mass W'], ['Hm', 'static height H_m']], v: 'W' },
      { k: 'W', label: 'Solids in the bed W', unit: 'kg', v: 7000, min: 1e-6, max: 1e9, scale: 'log', if: v => v.amt === 'W' },
      { k: 'Hm', label: 'Static (packed) bed height H_m', unit: 'm', v: 1.4, min: 1e-4, max: 100, if: v => v.amt === 'Hm' },
      { k: 'em', label: 'Voidage of the static bed ε_m', v: 0.4, min: 0.2, max: 0.8, if: v => v.amt === 'Hm' },
      { k: 'dt', label: 'Bed diameter', unit: 'm', v: 2, min: 1e-3, max: 50 },
      { k: 'rhos', label: 'Solid density ρ_s', unit: 'kg/m³', v: 2000, min: 10, max: 2e4 },
      { k: 'inert', label: 'Inert solids mixed in, per unit of catalyst', v: 0, min: 0, max: 100 },
      { k: 'u0', label: 'Superficial gas velocity u₀', unit: 'm/s', v: 0.3, min: 1e-4, max: 10 },
      { k: 'umf', label: 'u_mf', unit: 'm/s', v: 0.03, min: 1e-5, max: 5 },
      { k: 'emf', label: 'ε_mf', v: 0.5, min: 0.3, max: 0.8 },
      { k: 'db', label: 'Effective bubble size d_b', unit: 'm', v: 0.32, min: 1e-4, max: 5 },
      { k: 'D', label: 'Gas diffusivity D', unit: 'm²/s', v: 20e-6, min: 1e-10, max: 1e-2, scale: 'log' },
      { k: 'alpha', label: 'Wake fraction α', v: 0.33, min: 0, max: 2 },
      { k: 'fb', label: 'Solids in bubbles f_b', v: 0.001, min: 0, max: 0.01 },
      { k: 'dex', t: 'check', label: 'Use δ = (u₀ − u_mf)/(u_b − u_mf) instead of u₀/u_b', v: false },
      { group: 'Reaction (rate constants per volume of catalyst)' },
      { k: 'sch', t: 'select', label: 'Reaction', opts: [['first', 'A → R, first order'], ['rev', 'A ⇌ R, first order, reversible'], ['ser', 'A → R → S, first order'], ['par', 'A → R and A → S, first order']], v: 'first' },
      { k: 'k1', label: 'k‴ (or k‴₁)', unit: 'm³/(m³ cat·s)', v: 0.8, min: 1e-9, max: 1e6, scale: 'log' },
      { k: 'k2', label: v => (v.sch === 'ser' ? 'k‴₃ for R → S' : 'k‴₂ for A → S'), unit: 'm³/(m³ cat·s)', v: 0.05, min: 1e-9, max: 1e6, scale: 'log', if: v => v.sch === 'ser' || v.sch === 'par' },
      { k: 'Xe', label: 'Equilibrium conversion X_Ae', v: 0.94, min: 0.01, max: 0.9999, if: v => v.sch === 'rev' },
      { k: 'CA0', label: 'C_A0', unit: 'mol/m³', v: 100, min: 1e-9, max: 1e7 },
      { group: 'Design target (optional)' },
      { k: 'target', t: 'check', label: 'Find the solids needed for a conversion', v: false, if: v => v.sch !== 'ser' },
      { k: 'Xt', label: 'Conversion wanted', v: 0.9, min: 0.01, max: 0.9999, if: v => v.target && v.sch !== 'ser' }
    ],
    presets: [
      { name: 'Example 20.1 — 7 tons, 32 cm bubbles', v: { amt: 'W', W: 7000, dt: 2, rhos: 2000, inert: 0, u0: 0.3, umf: 0.03, emf: 0.5, db: 0.32, D: 20e-6, alpha: 0.33, fb: 0.001, dex: false, sch: 'first', k1: 0.8, CA0: 100, target: false }, book: { delta: { v: 0.196, tol: 0.01 }, ef: { v: 0.6, tol: 0.01 }, Kbc: { v: 0.614, tol: 0.01 }, Kce: { v: 0.133, tol: 0.01 }, fc: { v: 0.047, tol: 0.03 }, fe: { v: 0.352, tol: 0.01 }, H: { v: 2.785, tol: 0.01 }, X: { v: 0.32, tol: 0.02 }, Cbar: { v: 11, tol: 0.03 }, Xplug: { v: 0.95, tol: 0.01 } }, note: 'u₀ = 0.3 m/s through a 2 m bed of 7 tons of catalyst (u_mf = 0.03 m/s, ε_mf = 0.5, ρ_s = 2000 kg/m³), k‴ = 0.8 s⁻¹, D = 20 × 10⁻⁶ m²/s, α = 0.33, d_b = 0.32 m. Book: X_A = 32%, the solids see C̄_A = 11 mol/m³, a packed bed would give 95%.' },
      { name: 'Problem 20.1 — 14 tons', v: { amt: 'W', W: 14000, dt: 2, db: 0.32, u0: 0.3, sch: 'first', k1: 0.8, inert: 0, target: false } },
      { name: 'Problem 20.2 — 21 tons', v: { amt: 'W', W: 21000, dt: 2, db: 0.32, u0: 0.3, sch: 'first', k1: 0.8, inert: 0, target: false } },
      { name: 'Problem 20.3 — baffles, 16 cm bubbles', v: { amt: 'W', W: 7000, dt: 2, db: 0.16, u0: 0.3, sch: 'first', k1: 0.8, inert: 0, target: false }, book: { X: { v: 0.59, tol: 0.02 } }, note: 'Fig. 20.11: 59%.' },
      { name: 'Problem 20.4 — baffles, 8 cm bubbles', v: { amt: 'W', W: 7000, dt: 2, db: 0.08, u0: 0.3, sch: 'first', k1: 0.8, inert: 0, target: false }, book: { X: { v: 0.84, tol: 0.02 } }, note: 'Fig. 20.11: 84%.' },
      { name: 'Problem 20.5 — wider, shallower bed', v: { amt: 'W', W: 7000, dt: 2 * Math.SQRT2, db: 0.32, u0: 0.15, sch: 'first', k1: 0.8, inert: 0, target: false } },
      { name: 'Problem 20.7 — narrower, taller bed', v: { amt: 'W', W: 7000, dt: 2 / Math.SQRT2, db: 0.32, u0: 0.6, sch: 'first', k1: 0.8, inert: 0, target: false } },
      { name: 'Problem 20.8 — catalyst diluted 1 : 1 with inerts', v: { amt: 'W', W: 14000, dt: 2, db: 0.32, u0: 0.3, sch: 'first', k1: 0.8, inert: 1, target: false } },
      { name: 'Problem 20.9 — cumene, reversible (X_e = 94%)', v: { amt: 'Hm', Hm: 0.0762, em: 0.4, dt: 1, rhos: 2000, inert: 0, u0: 0.064, umf: 0.0061, emf: 0.5, db: 0.0135, D: 2e-5, alpha: 0.33, fb: 0.001, sch: 'rev', Xe: 0.94, k1: 1.4235, CA0: 1, target: false }, note: 'From the packed-bed run (60% conversion, H_m = 76.2 mm, u₀ = 64 mm/s, ε_m = 0.4): k‴ = 1.42 s⁻¹ for the reversible form −r‴ = k‴(C_A − C_Ae).' },
      { name: 'Problem 20.10 — nitrobenzene hydrogenation', v: { amt: 'Hm', Hm: 1.4, em: 0.4071, dt: 3.55, rhos: 2200, inert: 0, u0: 0.3, umf: 0.02, emf: 0.6, db: 0.1, D: 0.9e-4, alpha: 0.33, fb: 0.001, sch: 'first', k1: 1.2, CA0: 1, target: false } },
      { name: 'Problem 20.11 — pilot bed from a lab packed bed', v: { amt: 'Hm', Hm: 1, em: 0.5, dt: 1, rhos: 2000, inert: 0, u0: 0.2, umf: 0.032, emf: 0.5, db: 0.08, D: 3e-5, alpha: 0.34, fb: 0.001, sch: 'first', k1: 1.4026, CA0: 1, target: false }, note: 'The lab bed (10 cm, 2 cm/s, ε = 0.5, 97%) gives k‴ = 1.40 s⁻¹. The printed D = 3 × 10⁻⁵ is taken in m²/s, the usual size for a gas.' },
      { name: 'Figure 20.12 — A → R → S, plug flow limit', v: { amt: 'W', W: 7000, dt: 2, db: 0.32, u0: 0.3, sch: 'ser', k1: 0.8, k2: 0.05, inert: 0, CA0: 100 }, book: { CRplug: { v: 83, tol: 0.01 } }, note: 'k‴₁ = 0.8, k‴₃ = 0.05 s⁻¹: a fixed bed of 7 tons reaches C_R,max = 83. Fig. 20.12 gives 44 with 32 cm bubbles at 32 tons; the exact K–L balances solved here give about 47 at 32 tons. The figure’s numbers come from the closed-form multiple-reaction results of Levenspiel (1996), which are not reproduced here.' },
      { name: 'Figure 20.12 — A → R → S, 8 cm bubbles', v: { amt: 'W', W: 9500, dt: 2, db: 0.08, u0: 0.3, sch: 'ser', k1: 0.8, k2: 0.05, inert: 0, CA0: 100 }, book: { CRmax: { v: 73, tol: 0.03 } }, note: 'Fig. 20.12: 73 with 9.5 tons; the exact balances need about 10 tons for their maximum of 74.' },
      { name: 'Figure 20.12 — A → R → S, 16 cm bubbles', v: { amt: 'W', W: 15000, dt: 2, db: 0.16, u0: 0.3, sch: 'ser', k1: 0.8, k2: 0.05, inert: 0, CA0: 100 }, book: { CRmax: { v: 60, tol: 0.03 } }, note: 'Fig. 20.12: 60 with 15 tons; the exact balances reach 61 with about 17 tons.' },
      { name: 'Catalyst for 90% with 32 cm bubbles', v: { amt: 'W', W: 7000, dt: 2, db: 0.32, u0: 0.3, sch: 'first', k1: 0.8, inert: 0, target: true, Xt: 0.9 } }
    ],
    eqs: [
      ['u_{br} = 0.711(gd_b)^{1/2},\\qquad u_b = u_0 - u_{mf} + u_{br}', 'Eqs. 20.10–20.11'],
      ['\\delta = \\frac{u_0-u_{mf}}{u_b-u_{mf}} \\cong \\frac{u_0}{u_b},\\qquad 1-\\varepsilon_f = (1-\\varepsilon_{mf})(1-\\delta)', 'Eq. 20.12'],
      ['K_{bc} = 4.50\\frac{u_{mf}}{d_b} + 5.85\\frac{\\mathscr{D}^{1/2}g^{1/4}}{d_b^{5/4}}', 'Eq. 20.13'],
      ['K_{ce} = 6.77\\left(\\frac{\\varepsilon_{mf}\\mathscr{D}u_{br}}{d_b^3}\\right)^{1/2}', 'Eq. 20.14'],
      ['f_c = \\delta(1-\\varepsilon_{mf})\\left[\\frac{3u_{mf}/\\varepsilon_{mf}}{u_{br}-u_{mf}/\\varepsilon_{mf}}+\\alpha\\right],\\quad f_e = (1-\\varepsilon_{mf})(1-\\delta)-f_c-f_b', 'Eqs. 20.16–20.17'],
      ['\\ln\\frac{C_{A0}}{C_A} = \\left[f_bk^{\\prime\\prime\\prime} + \\cfrac{1}{\\cfrac{1}{\\delta K_{bc}} + \\cfrac{1}{f_ck^{\\prime\\prime\\prime} + \\cfrac{1}{\\cfrac{1}{\\delta K_{ce}}+\\cfrac{1}{f_ek^{\\prime\\prime\\prime}}}}}\\right]\\frac{H_{BFB}}{u_0}', 'Eq. 20.21'],
      ['\\bar C_A = \\frac{C_{A0}-C_A}{k^{\\prime\\prime\\prime}\\tau^{\\prime\\prime\\prime}},\\qquad H_{BFB} = \\frac{W}{\\rho_sA(1-\\varepsilon_f)}', 'Eqs. 20.19, 20.22'],
      ['\\ln\\frac{C_{A0}}{C_{Ap}} = k^{\\prime\\prime\\prime}\\tau^{\\prime\\prime\\prime} = \\frac{k^{\\prime\\prime\\prime}H_p(1-\\varepsilon_p)}{u_0}', 'Eq. 20.23', 'The same catalyst as a packed bed in plug flow.']
    ],
    method: ['Bubble, cloud and emulsion concentrations are linked by quasi-steady balances: for linear kinetics with rate matrix K, the cloud and emulsion concentrations are eliminated by solving two small linear systems at each level, leaving a bubble-phase equation u₀ dC_b/dz = −A C_b that is integrated exactly with a matrix exponential. For a single first-order reaction this is Eq. 20.21 term for term. Reversible reactions are written for C_A − C_Ae; series and parallel reactions use the 2 × 2 rate matrix, and the best catalyst load for the intermediate is found by golden-section search on the height.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const A = Math.PI * +v.dt * +v.dt / 4, rs = +v.rhos;
      const p = { u0: +v.u0, umf: +v.umf, emf: +v.emf, db: +v.db, D: +v.D, alpha: +v.alpha, fb: +v.fb, deltaExact: !!v.dex };
      if (!(p.u0 > p.umf)) throw new Error('u₀ must exceed u_mf for a bubbling bed.');
      const q = FZ.kl(p);
      const W = v.amt === 'W' ? +v.W : +v.Hm * (1 - +v.em) * A * rs;
      const H = W / (rs * A * (1 - q.ef));
      const cat = 1 / (1 + (+v.inert || 0));
      const k1 = +v.k1 * cat, k2 = +v.k2 * cat;
      const C0 = +v.CA0;
      kpi.push({ k: 'ubr', l: 'Rise velocity of a single bubble u_br', v: q.ubr, u: 'm/s', s: 'u_b = ' + fmt(q.ub, 4) + ' m/s' });
      kpi.push({ k: 'delta', l: 'Bed fraction in bubbles δ', v: q.delta });
      kpi.push({ k: 'ef', l: 'Voidage of the bubbling bed ε_f', v: q.ef });
      kpi.push({ k: 'Kbc', l: 'K_bc (bubble–cloud)', v: q.Kbc, u: '1/s' }, { k: 'Kce', l: 'K_ce (cloud–emulsion)', v: q.Kce, u: '1/s' });
      kpi.push({ k: 'fc', l: 'Solids in cloud and wake f_c', v: q.fc, s: 'f_b = ' + fmt(q.fb, 3) }, { k: 'fe', l: 'Solids in the emulsion f_e', v: q.fe });
      kpi.push({ k: 'H', l: 'Height of the bubbling bed', v: H, u: 'm', s: 'H_mf = ' + fmt(W / (rs * A * (1 - p.emf)), 4) + ' m' });
      if (q.fast < 5) warn.push('u_b/(u_mf/ε_mf) = ' + fmt(q.fast, 3) + ': bubbles are not fast, clouds are thick, and the K–L model is stretched.');
      if (p.db > 0.6 * +v.dt) warn.push('Bubbles this large in this vessel would slug; the model assumes freely bubbling beds.');
      if (q.fe < 0) throw new Error('These values leave no solids for the emulsion (f_e < 0): check α, δ and ε_mf.');
      if (v.sch === 'first' || v.sch === 'rev') {
        const Kf = FZ.kf(q, k1);
        const g = Kf * H / p.u0;
        const Xe = v.sch === 'rev' ? +v.Xe : 1;
        const X = Xe * (1 - Math.exp(-g));
        const kt = +v.k1 * W / rs / (p.u0 * A) * cat;
        const XplugB = Xe * (1 - Math.exp(-kt)), Xmix = Xe * kt / (1 + kt);
        kpi.unshift({ k: 'X', l: 'Conversion in the bubbling bed', v: X, main: true, s: 'packed bed ' + fmt(XplugB, 3) + ', mixed flow ' + fmt(Xmix, 3) });
        kpi.push({ k: 'Xplug', l: 'Same catalyst as a packed bed, plug flow (Eq. 20.23)', v: XplugB });
        kpi.push({ k: 'Kf', l: 'Effective rate constant K‴ (Eq. 20.21)', v: Kf / ((1 - q.ef) * cat), u: '1/s', s: 'k‴ = ' + fmt(k1 / cat, 4) + ' s⁻¹; ratio ' + fmt(Kf / ((1 - q.ef) * k1), 3) });
        const CA = C0 * (1 - X);
        const Cbar = v.sch === 'rev' ? NaN : (C0 - CA) * p.u0 * A / (k1 * W / rs);
        if (isFinite(Cbar)) kpi.push({ k: 'Cbar', l: 'Mean C_A seen by the solids (Eq. 20.22)', v: Cbar, u: 'mol/m³', s: 'gas leaves at ' + fmt(CA, 4) + ' mol/m³' });
        /* resistances */
        const R = [['bubble → cloud', 1 / (q.delta * q.Kbc)], ['cloud → emulsion', 1 / (q.delta * q.Kce)], ['reaction in bubble', 1 / (q.fb * k1)], ['reaction in cloud', 1 / (q.fc * k1)], ['reaction in emulsion', 1 / (q.fe * k1)]];
        tables.push({ title: 'The five resistances of Eq. 20.21 (per volume of bed, s)', cols: ['step', '1/(δK) or 1/(fk‴)', 'compared with 1/K_f'], rows: R.map(r => [r[0], r[1], r[1] * Kf]) });
        if (v.target) {
          const Xt = +v.Xt;
          if (!(Xt < Xe)) throw new Error('The conversion wanted is at or beyond equilibrium.');
          const Hn = p.u0 * Math.log(1 / (1 - Xt / Xe)) / Kf, Wn = Hn * rs * A * (1 - q.ef);
          kpi.push({ k: 'Wneed', l: 'Solids needed for X_A = ' + fmt(Xt, 3), v: Wn, u: 'kg', s: 'bed height ' + fmt(Hn, 4) + ' m; packed bed would need ' + fmt(Math.log(1 / (1 - Xt / Xe)) * p.u0 * A * rs / (+v.k1 * cat), 4) + ' kg' });
        }
        /* Fig. 20.10: C/C0 against k‴τ‴ for several bubble sizes */
        const taus = N.linspace(0, 8, 120);
        const curve = db => { const qq = FZ.kl(Object.assign({}, p, { db })); const r = FZ.kf(qq, k1) / ((1 - qq.ef) * k1); return taus.map(t => Math.exp(-r * t)); };
        const ser = [{ name: 'plug flow', x: taus, y: taus.map(t => Math.exp(-t)), color: 'var(--muted)', dash: '5 4' }, { name: 'mixed flow', x: taus, y: taus.map(t => 1 / (1 + t)), color: 'var(--muted)', dash: '2 3' }];
        [0.08, 0.16, 0.32].forEach((db, i) => ser.push({ name: 'd_b = ' + fmt(db * 100, 3) + ' cm', x: taus, y: curve(db), color: i + 3 }));
        ser.push({ name: 'this bed, d_b = ' + fmt(p.db * 100, 3) + ' cm', x: taus, y: curve(p.db), color: 1, w: 2.4 });
        const tNow = k1 * (1 - q.ef) * H / p.u0;
        plots.push({ title: 'Bubble size and performance (Fig. 20.10)', x: { label: 'f_total H k‴/u₀', min: 0, max: 8 }, y: { label: v.sch === 'rev' ? '1 − X_A/X_Ae' : 'C_A/C_A0', log: true, min: 0.01, max: 1 }, series: ser, marks: [{ t: 'pt', x: tNow, y: 1 - X / Xe, color: 8, label: 'this bed' }], cap: 'Large bubbles bypass the solids: the bed can fall far below mixed flow.' });
        /* X against bubble size */
        const dbs = N.logspace(0.01, 1, 80);
        plots.push({ title: 'Conversion against bubble size at this catalyst load', x: { label: 'd_b', unit: 'm', log: true }, y: { label: 'X_A', min: 0, max: 1 }, series: [{ name: 'bubbling bed', x: dbs, y: dbs.map(db => { const qq = FZ.kl(Object.assign({}, p, { db })); const HH = W / (rs * A * (1 - qq.ef)); return Xe * (1 - Math.exp(-FZ.kf(qq, k1) * HH / p.u0)); }), color: 1, w: 2.2 }, { name: 'packed bed', x: [0.01, 1], y: [XplugB, XplugB], color: 'var(--muted)', dash: '5 4' }, { name: 'mixed flow', x: [0.01, 1], y: [Xmix, Xmix], color: 'var(--muted)', dash: '2 3' }], marks: [{ t: 'v', x: p.db, color: 8, label: 'd_b' }], cap: 'Baffles or internals that keep bubbles small are the main lever on a bubbling bed (Problems 20.3–20.4).' });
        /* profiles up the bed: bubble, cloud, emulsion */
        const M = FZ.bubbleMatrix(q, [[k1]]);
        const zs = N.linspace(0, H, 80);
        const Cb = zs.map(z => Math.exp(-M.A[0][0] * z / p.u0));
        plots.push({ title: 'Reactant up the bed: bubble, cloud and emulsion', x: { label: 'height z', unit: 'm', min: 0, max: H }, y: { label: v.sch === 'rev' ? '(C − C_e)/(C_A0 − C_e)' : 'C_A/C_A0', min: 0, max: 1 }, series: [{ name: 'bubble gas', x: zs, y: Cb, color: 1, w: 2.2 }, { name: 'cloud', x: zs, y: Cb.map(c => c * M.Mc[0][0]), color: 2 }, { name: 'emulsion', x: zs, y: Cb.map(c => c * M.Me[0][0]), color: 3 }], cap: 'The solids in the emulsion see far less reactant than the gas leaving the bed carries.' });
        return { kpi, warn, info, plots, tables };
      }
      /* series and parallel: 2 × 2 rate matrix on [C_A, C_R] */
      const Km = v.sch === 'ser' ? [[k1, 0], [-k1, k2]] : [[k1 + k2, 0], [-k1, 0]];
      const M = FZ.bubbleMatrix(q, Km);
      const out = HH => FZ.propagate(M.A, p.u0, [C0, 0], HH);
      const C = out(H);
      if (v.sch === 'par') {
        const phi = k1 / (k1 + k2);
        kpi.unshift({ k: 'X', l: 'Conversion of A', v: 1 - C[0] / C0, main: true });
        kpi.push({ k: 'CR', l: 'C_R leaving', v: C[1], u: 'mol/m³', s: 'fraction of reacted A to R = ' + fmt(C[1] / (C0 - C[0]), 4) + ' (= k₁/(k₁ + k₂) = ' + fmt(phi, 4) + ', same as in a fixed bed)' });
        info.push('Side-by-side reactions of the same order keep their product distribution in a bubbling bed; only the conversion suffers.');
        return { kpi, warn, info, plots, tables };
      }
      /* series: best height for R */
      const Hof = Wq => Wq / (rs * A * (1 - q.ef));
      const fR = Wq => -out(Hof(Wq))[1];
      let lo = W * 1e-3, hi = W * 1e3;
      const Ws = N.logspace(lo, hi, 60), CRs = Ws.map(w => -fR(w));
      let bi = 0; CRs.forEach((c, i) => { if (c > CRs[bi]) bi = i; });
      const best = N.fmin(w => fR(Math.exp(w)), Math.log(Ws[Math.max(0, bi - 1)]), Math.log(Ws[Math.min(59, bi + 1)]));
      const Wb = Math.exp(best.x), CRmax = -best.f, Cb = out(Hof(Wb));
      const kap = k2 / k1, CRplug = C0 * Math.pow(kap, kap / (1 - kap));
      kpi.unshift({ k: 'CRmax', l: 'Best C_R in this bubbling bed', v: CRmax, u: 'mol/m³', main: true, s: 'X_A = ' + fmt(1 - Cb[0] / C0, 3) + ' at that load' });
      kpi.push({ k: 'Wmax', l: 'Catalyst for the best C_R', v: Wb, u: 'kg', s: 'bed height ' + fmt(Hof(Wb), 4) + ' m' });
      kpi.push({ k: 'CRplug', l: 'Best C_R in a fixed bed (plug flow)', v: CRplug, u: 'mol/m³', s: 'needs ' + fmt(Math.log(1 / kap) / (1 - kap) / (+v.k1 * cat) * p.u0 * A * rs, 4) + ' kg' });
      kpi.push({ k: 'CRnow', l: 'C_R leaving with the solids entered', v: C[1], u: 'mol/m³', s: 'X_A = ' + fmt(1 - C[0] / C0, 3) });
      const Wg = N.logspace(W * 0.02, Math.max(W, Wb) * 6, 140);
      const plugR = w => { const t = +v.k1 * cat * w / rs / (p.u0 * A); return C0 * (Math.abs(kap - 1) < 1e-9 ? t * Math.exp(-t) : (Math.exp(-t) - Math.exp(-kap * t)) / (kap - 1)); };
      plots.push({ title: 'Intermediate R against catalyst (Fig. 20.12)', x: { label: 'catalyst W', unit: 'kg', log: true }, y: { label: 'C_R', unit: 'mol/m³', min: 0 }, series: [{ name: 'bubbling bed, d_b = ' + fmt(p.db * 100, 3) + ' cm', x: Wg, y: Wg.map(w => out(Hof(w))[1]), color: 1, w: 2.2 }, { name: 'fixed bed, plug flow', x: Wg, y: Wg.map(plugR), color: 'var(--muted)', dash: '5 4' }], marks: [{ t: 'pt', x: Wb, y: CRmax, color: 8, label: 'best' }] });
      const rows = [0.04, 0.08, 0.16, 0.32, 0.64].map(db => {
        const qq = FZ.kl(Object.assign({}, p, { db })), MM = FZ.bubbleMatrix(qq, Km);
        const Hq = w => w / (rs * A * (1 - qq.ef));
        const f = lw => -FZ.propagate(MM.A, p.u0, [C0, 0], Hq(Math.exp(lw)))[1];
        const sc = N.fminScan(f, Math.log(W * 1e-3), Math.log(W * 1e3), 50);
        return [db * 100, -sc.f, Math.exp(sc.x), 1 - FZ.propagate(MM.A, p.u0, [C0, 0], Hq(Math.exp(sc.x)))[0] / C0];
      });
      tables.push({ title: 'Bubble size and the best yield of R (Fig. 20.12)', cols: ['d_b, cm', 'C_R,max', 'W needed, kg', 'X_A there'], rows: [['fixed bed', CRplug, Math.log(1 / kap) / (1 - kap) / (+v.k1 * cat) * p.u0 * A * rs, 1 - Math.pow(kap, 1 / (1 - kap))]].concat(rows) });
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* C-2003  turbulent, fast fluidized and pneumatic beds: f(z) and conversion */
  /* ===================================================================== */
  ST.tool({
    id: 'cfb', tag: 'C-2003', ch: 20, sec: 'Ch. 20.4', title: 'Turbulent, fast fluidized and pneumatic beds: solids profile and conversion',
    lede: 'The circulating-bed regimes with the book’s solids distribution f = f* + (f_d − f*) e^{−az}: how the catalyst inventory splits between the dense bottom and the lean freeboard, how much catalyst the gas actually meets, and the first-order conversion it can give — the plug-flow upper bound of Eq. 18.50, beside a pessimistic estimate with the dense zone well mixed.',
    keys: 'circulating fluidized bed cfb turbulent bed fast fluidized pneumatic conveying riser downer solids fraction exponential decay decay constant u0 a dense region lean region inventory conversion first order suspended solids eq 18.50 fcc riser downflow',
    inputs: [
      { group: 'Regime and solids' },
      { k: 'reg', t: 'select', label: 'Regime', opts: [['TB', 'turbulent bed'], ['FF', 'fast fluidized'], ['PC', 'pneumatic conveying']], v: 'FF' },
      { k: 'gcls', t: 'seg', label: 'Geldart class', opts: [['A', 'A (u₀a ≈ 2–4)'], ['AB', 'AB (≈ 5)'], ['B', 'B (≈ 7)']], v: 'A' },
      { k: 'u0a', label: 'u₀ a', unit: '1/s', v: 3, min: 0.1, max: 30 },
      { k: 'u0', label: 'Superficial gas velocity u₀', unit: 'm/s', v: 4, min: 0.05, max: 50 },
      { k: 'fs', label: 'f* (lean limit)', v: 0.015, min: 0.001, max: 0.1 },
      { k: 'fd', label: 'f_d (dense region)', v: 0.12, min: 0.01, max: 0.6 },
      { k: 'Ht', label: 'Height of the vessel', unit: 'm', v: 20, min: 0.1, max: 200 },
      { k: 'inv', t: 'seg', label: 'Size the dense region from', opts: [['Hd', 'its height'], ['W', 'the solids inventory']], v: 'Hd' },
      { k: 'Hd', label: 'Height of the dense region', unit: 'm', v: 2, min: 0, max: 200, if: v => v.inv === 'Hd' },
      { k: 'W', label: 'Solids inventory', unit: 'kg', v: 20000, min: 1, max: 1e9, scale: 'log', if: v => v.inv === 'W' },
      { k: 'dt', label: 'Vessel diameter', unit: 'm', v: 1, min: 0.01, max: 50 },
      { k: 'rhos', label: 'Solid density ρ_s', unit: 'kg/m³', v: 1500, min: 10, max: 2e4 },
      { group: 'Reaction' },
      { k: 'k', label: 'k‴, first order', unit: 'm³/(m³ cat·s)', v: 5, min: 1e-9, max: 1e6, scale: 'log' },
      { k: 'down', t: 'check', label: 'Downflow reactor (no dense bottom, near plug flow)', v: false }
    ],
    presets: [
      { name: 'Fast fluidized riser, Geldart A', v: { reg: 'FF', gcls: 'A', u0a: 3, u0: 4, fs: 0.015, fd: 0.12, Ht: 20, inv: 'Hd', Hd: 2, dt: 1, rhos: 1500, k: 5, down: false } },
      { name: 'Turbulent bed', v: { reg: 'TB', gcls: 'AB', u0a: 5, u0: 1.2, fs: 0.015, fd: 0.3, Ht: 10, inv: 'Hd', Hd: 3, dt: 2, rhos: 1500, k: 2, down: false } },
      { name: 'Pneumatic conveying', v: { reg: 'PC', gcls: 'B', u0a: 7, u0: 10, fs: 0.01, fd: 0.06, Ht: 25, inv: 'Hd', Hd: 0, dt: 1, rhos: 1500, k: 5, down: false } },
      { name: 'Downflow reactor (Fig. 20.17)', v: { reg: 'PC', gcls: 'A', u0a: 3, u0: 10, fs: 0.015, fd: 0.06, Ht: 20, inv: 'Hd', Hd: 0, dt: 1, rhos: 1500, k: 5, down: true } }
    ],
    eqs: [
      ['f = f^* + (f_d - f^*)e^{-az}', 'Figs. 20.13–20.17', 'Lean region above the dense one.'],
      ['u_0\\mathbf{a} = 2\\text{–}4\\ \\text{s}^{-1}\\ (\\text{A}),\\quad 5\\ (\\text{AB}),\\quad 7\\ (\\text{B})', 'Ch. 20.4'],
      ['\\ln\\frac{C_{A0}}{C_A} = \\frac{k^{\\prime\\prime\\prime}}{u_0}\\int_0^H f\\,dz', 'Eq. 18.51', 'Plug flow of gas, first order.']
    ],
    method: ['The solids profile is uniform at f_d up to the dense height and decays exponentially above it toward f*. When the inventory is given, the dense height is found by Brent’s method from W = ρ_s A ∫f dz. Conversion uses plug flow of gas through the whole profile (Eq. 18.51), the best this contacting can do; the pessimistic case treats the dense zone as one mixed region followed by plug flow in the lean zone. The book stresses that reliable flow models for these regimes do not yet exist.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const A = Math.PI * +v.dt * +v.dt / 4, rs = +v.rhos, u0 = +v.u0, a = +v.u0a / u0, fs = +v.fs, fd = +v.fd, Ht = +v.Ht, k = +v.k;
      if (!(fd >= fs)) throw new Error('f_d must be at least f*.');
      const integ = Hd => { const H1 = Math.min(Hd, Ht), L = Ht - H1; return fd * H1 + fs * L + (fd - fs) * (1 - Math.exp(-a * L)) / a; };
      let Hd = v.down ? 0 : +v.Hd;
      if (!v.down && v.inv === 'W') {
        const target = +v.W / (rs * A);
        if (target > fd * Ht) throw new Error('That inventory does not fit even with the whole vessel dense.');
        if (target < integ(0)) { Hd = 0; warn.push('The inventory is less than the lean profile alone holds; the dense region vanishes.'); }
        else Hd = N.brent(h => integ(h) - target, 0, Ht, { tol: 1e-10 });
      }
      const I = integ(Hd), W = I * rs * A;
      const fz = z => (z <= Hd ? fd : fs + (fd - fs) * Math.exp(-a * (z - Hd)));
      const Xplug = 1 - Math.exp(-k * I / u0);
      const dense = k * fd * Hd / u0, lean = k * (I - fd * Hd) / u0;
      const Xpess = 1 - Math.exp(-lean) / (1 + dense);
      kpi.push({ k: 'X', l: 'Conversion, plug flow of gas (upper bound)', v: Xplug, main: true, s: 'dense zone mixed instead: ' + fmt(Xpess, 4) });
      kpi.push({ k: 'If', l: '∫₀ᴴ f dz', v: I, u: 'm' }, { k: 'W', l: 'Solids in the vessel', v: W, u: 'kg', s: 'dense ' + fmt(fd * Hd * rs * A, 4) + ' kg, lean ' + fmt((I - fd * Hd) * rs * A, 4) + ' kg' });
      kpi.push({ k: 'Hd', l: 'Height of the dense region', v: Hd, u: 'm' }, { k: 'a', l: 'Decay constant a', v: a, u: '1/m', s: 'u₀a = ' + fmt(+v.u0a, 3) + ' s⁻¹' });
      kpi.push({ k: 'tg', l: 'Gas residence time H/u₀', v: Ht / u0, u: 's' });
      if (v.reg === 'FF' && !v.down) info.push('In the fast fluidized bed a dense, downflowing wall zone carries part of the solids and of the gas in near stagnant contact (Fig. 20.16); its exchange coefficient K_cw is not yet known, so this estimate leaves it out.');
      if (v.down) info.push('Downflow reactors keep gas and solids near plug flow and avoid the backmixed wall zone of upflow risers — the reason for interest in them for catalytic cracking.');
      const zs = N.linspace(0, Ht, 200);
      const cum = N.cumtrapz(zs, zs.map(fz));
      plots.push({ title: 'Solids up the vessel and conversion', x: { label: 'height z', unit: 'm', min: 0, max: Ht }, y: { label: 'f, X_A', min: 0, max: 1 }, series: [{ name: 'solids fraction f', x: zs, y: zs.map(fz), color: 3, type: 'area', line: true }, { name: 'X_A, plug flow of gas', x: zs, y: cum.map(c => 1 - Math.exp(-k * c / u0)), color: 1, w: 2.2 }], cap: 'The dense bottom holds most of the solids; the lean region still adds catalyst contact all the way up.' });
      tables.push({ title: 'Typical values (Figs. 20.13–20.17)', cols: ['regime', 'f*', 'f_d', 'u₀a, s⁻¹'], rows: [['turbulent', '0.01–0.02', '0.2–0.4', '2–7'], ['fast fluidized', '0.01–0.02', '0.06–0.2', '2–7'], ['pneumatic', '0.01–0.02', '≈ 0.06 at entry', '2–7']] });
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
