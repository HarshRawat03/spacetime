/* Chapters 17–18 (I) — rate steps in series and parallel; effectiveness factors for any kinetics;
   diagnosing film, pore and heat resistances from measurements */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, PEL = ST.PEL;
  const CAT = ST.CAT = ST.CAT || {};

  /* characteristic length L = V/S from the size the user knows */
  CAT.Lof = (shape, d) => (shape === 'slab' ? d / 2 : shape === 'cyl' ? d / 4 : d / 6);
  CAT.sizeWord = shape => (shape === 'slab' ? 'thickness' : 'diameter');
  CAT.regime = M => (M < 0.4 ? 'diffusion-free' : M > 4 ? 'strong pore diffusion' : 'intermediate');
  CAT.regimeW = MW => (MW < 0.15 ? 'diffusion-free' : MW > 4 ? 'strong pore diffusion' : 'intermediate');
  CAT.state = M => (M < 0.4 ? 'ok' : M > 4 ? 'bad' : 'warn');
  const fmt = (x, s = 4) => N.fmt(x, s);

  /* ===================================================================== */
  /* C-1701  rate steps in series and in parallel                           */
  /* ===================================================================== */
  const kU = v => (v.basis === 'vol' ? '1/[t]' : 'm/[t]');
  const rU = v => (v.basis === 'vol' ? 'mol/(m³·[t])' : 'mol/(m²·[t])');
  const krU = v => {
    const n = v.law === 'nth' ? +v.n : 1, base = v.basis === 'vol' ? '1/[t]' : 'm/[t]';
    if (v.law === 'lh' || n === 1) return base;
    const p = 1 - n;
    return 'mol^{' + fmt(p, 3) + '}·m^{' + fmt(-3 * p, 3) + '}·' + base;
  };

  function chain(v) {
    const law = v.law, kr = +v.kr * (+v.eta || 1), n = +v.n, Ce = +v.Ce || 0, K = +v.K || 0;
    const rxn = law === 'rev' ? (c => kr * (c - Ce)) : law === 'lh' ? (c => (c <= 0 ? 0 : kr * c / (1 + K * c))) : (c => (c <= 0 ? 0 : kr * Math.pow(c, n)));
    const H = v.i ? +v.H : 1;
    const kg = v.g ? +v.kg : Infinity, kl = v.l ? +v.kl : Infinity, ks = v.s ? +v.ks : Infinity;
    const rb = v.b ? (c => (c <= 0 ? 0 : +v.kb * Math.pow(c, +v.nb))) : (() => 0);
    /* march backward from a guessed surface concentration to the main fluid */
    const back = Cs => {
      const r = rxn(Cs);
      const Clb = Cs + r / ks;
      const rB = rb(Clb), J = r + rB;
      const Cli = Clb + J / kl;
      const Cgi = H * Cli;
      const C0 = Cgi + J / kg;
      return { Cs, r, Clb, rB, J, Cli, Cgi, C0 };
    };
    return { back, rxn, H, kg, kl, ks, rb, Ce: law === 'rev' ? Ce : 0 };
  }
  function solveChain(v, C0) {
    const ch = chain(v);
    const lo = ch.Ce, hi = C0 / ch.H;
    if (!(hi > lo)) throw new Error(v.law === 'rev' ? 'The main-fluid concentration must lie above equilibrium (C₀/H > C_e).' : 'The concentration must be positive.');
    const g = Cs => ch.back(Cs).C0 - C0;
    const Cs = N.brent(g, lo, hi, { tol: 1e-15 * Math.max(1, hi) });
    if (!isFinite(Cs)) throw new Error('Could not balance the steps; check the coefficients.');
    return Object.assign(ch.back(Cs), { ch });
  }

  ST.tool({
    id: 'resist', tag: 'C-1701', ch: 17, sec: 'Ch. 17', title: 'Rate steps in series and in parallel: film, interface, bulk and surface',
    lede: 'Build the chain a reactant travels — gas film, gas–liquid interface, liquid film, the main body of liquid (where it may also react), the film around a solid, and the reaction at the surface with any rate form. Every intermediate concentration is eliminated exactly, so nonlinear steps need no first-order approximation; the tool shows where the driving force is spent and which step is worth improving.',
    keys: 'heterogeneous resistances in series parallel film mass transfer surface reaction overall rate example 17.1 17.2 17.3 17.4 problem 17.1 17.3 slurry burning carbon fermentation langmuir hinshelwood henry partition controlling step',
    units: { t: 's' },
    inputs: [
      { group: 'Basis' },
      { k: 'basis', t: 'seg', label: 'Rates and coefficients per unit of', opts: [['area', 'surface: k, rates per m²'], ['vol', 'volume: k·a, rates per m³']], v: 'area' },
      { k: 'C0', label: 'Concentration of A in the main fluid', unit: 'mol/m³', v: 1, min: 1e-9, max: 1e7, scale: 'log' },
      { group: 'Steps in series, from the main fluid inward' },
      { k: 'g', t: 'check', label: 'Film on the fluid side (gas film)', v: true },
      { k: 'kg', label: 'Film coefficient k_g', unit: kU, v: 0.05, min: 1e-12, max: 1e6, scale: 'log', if: v => v.g },
      { k: 'i', t: 'check', label: 'Interface with phase partition (gas → liquid)', v: false },
      { k: 'H', label: 'Partition H = C_gas / C_liquid at equilibrium', v: 1, min: 1e-9, max: 1e9, scale: 'log', if: v => v.i },
      { k: 'l', t: 'check', label: 'Liquid film at the interface', v: false },
      { k: 'kl', label: 'Liquid-film coefficient k_l', unit: kU, v: 0.01, min: 1e-12, max: 1e6, scale: 'log', if: v => v.l },
      { k: 'b', t: 'check', label: 'A also reacts in the main body of liquid (parallel step)', v: false },
      { k: 'kb', label: 'Bulk rate constant k_b', unit: v => (+v.nb === 1 ? kU(v) : '(m³/mol)^{' + fmt(+v.nb - 1, 3) + '}·' + kU(v)), v: 0.01, min: 0, max: 1e9, if: v => v.b },
      { k: 'nb', label: 'Order of the bulk reaction', v: 1, min: 0, max: 4, step: 0.5, if: v => v.b },
      { k: 's', t: 'check', label: 'Film around the solid (catalyst or reacting particle)', v: false },
      { k: 'ks', label: 'Film coefficient to the solid k_s', unit: kU, v: 0.02, min: 1e-12, max: 1e6, scale: 'log', if: v => v.s },
      { group: 'Last step: reaction at the surface' },
      { k: 'law', t: 'select', label: 'Rate form', opts: [['nth', 'k″ C^n'], ['rev', 'k″ (C − C_e), reversible first order'], ['lh', 'k″ C / (1 + K C), Langmuir–Hinshelwood']], v: 'nth' },
      { k: 'kr', label: 'Rate constant k″', unit: krU, v: 0.02, min: 1e-15, max: 1e9, scale: 'log' },
      { k: 'n', label: 'Order n', v: 1, min: 0, max: 4, step: 0.5, if: v => v.law === 'nth' },
      { k: 'Ce', label: 'Equilibrium concentration C_e at the surface', unit: 'mol/m³', v: 0.2, min: 0, max: 1e7, if: v => v.law === 'rev' },
      { k: 'K', label: 'Adsorption constant K', unit: 'm³/mol', v: 5, min: 0, max: 1e9, if: v => v.law === 'lh' },
      { k: 'eta', label: 'Effectiveness factor of the solid (1 if nonporous)', v: 1, min: 1e-6, max: 1, hint: 'Porous solids: take η from C-1801.' }
    ],
    presets: [
      { name: 'Example 17.3 — linear: film then first-order surface', v: { basis: 'area', C0: 1, g: true, kg: 0.05, i: false, l: false, b: false, s: false, law: 'nth', kr: 0.02, n: 1, eta: 1 }, expect: { r: { v: 1 / 70, tol: 1e-6 } }, note: '1/k_g and 1/k″ add: −r″ = C/(1/k_g + 1/k″). With k_g = 0.05 and k″ = 0.02 m/s and C = 1 mol/m³ the overall rate is 1/70 mol/m²·s.' },
      { name: 'Example 17.4 — nonlinear: film then second order', v: { basis: 'area', C0: 5, g: true, kg: 0.05, i: false, l: false, b: false, s: false, law: 'nth', kr: 0.02, n: 2, eta: 1 }, expect: { r: { v: 0.125, tol: 1e-6 } }, note: 'The second-order surface step cannot be added as a resistance. The book’s closed form −r″ = k_g C − (k_g²/2k″)(√(1 + 4k″C/k_g) − 1) gives 0.125 mol/m²·s here; the chain solver gives the same without the algebra.' },
      { name: 'Problem 17.1 — reversible surface step', v: { basis: 'area', C0: 1, g: true, kg: 0.05, i: false, l: false, b: false, s: false, law: 'rev', kr: 0.02, Ce: 0.2, eta: 1 }, expect: { r: { v: 0.8 / 70, tol: 1e-6 } }, note: 'Gas film then a reversible first-order surface reaction: −r″ = (C − C_e)/(1/k_g + 1/k″).' },
      { name: 'Problem 17.3 — slurry reactor', v: { basis: 'vol', C0: 40, g: false, i: true, H: 2.5, l: true, kl: 0.05, b: false, s: true, ks: 0.2, law: 'nth', kr: 0.1, n: 1, eta: 1 }, expect: { r: { v: 16 / (20 + 5 + 10), tol: 1e-6 } }, note: 'Pure gas bubbles; the reactant dissolves (C* = p/H), crosses the liquid film at the bubble (k_l a_i), the film around the catalyst (k_c a_c) and reacts first order on the catalyst (k a_c): −r = C*/(1/k_l a_i + 1/k_c a_c + 1/k a_c).' },
      { name: 'Example 17.1 — burning carbon particle (illustrative values)', v: { basis: 'area', C0: 2.1, g: true, kg: 0.6, i: false, l: false, b: false, s: false, law: 'nth', kr: 0.9, n: 1, eta: 1 }, note: 'Two steps in series: oxygen diffuses through the gas film and reacts at the carbon surface. Values are illustrative (air at about 1200 K).' },
      { name: 'Example 17.2 — aerobic fermentation (illustrative values)', v: { basis: 'vol', C0: 8, g: true, kg: 2, i: true, H: 30, l: true, kl: 0.08, b: false, s: true, ks: 0.5, law: 'lh', kr: 0.4, K: 15, eta: 0.6 }, note: 'Oxygen from the bubble: gas film, interface, liquid film, main liquid, film around the cell clump, diffusion into the clump (η) and uptake by the cells (saturation kinetics). Values are illustrative; the point is to see which steps matter.' },
      { name: 'Saturation kinetics behind a film: order shifts with concentration', v: { basis: 'area', C0: 1, g: true, kg: 0.01, i: false, l: false, b: false, s: false, law: 'lh', kr: 0.05, K: 2, eta: 1 }, note: 'At low C the rate is first order and film-limited; at high C the surface saturates and the order falls toward zero. The apparent order is read off the rate curve.' }
    ],
    eqs: [
      ['-r_A^{\\prime\\prime} = k_g\\,(C_{Ag}-C_{As}) = k^{\\prime\\prime} C_{As} \\;\\Rightarrow\\; -r_A^{\\prime\\prime} = \\frac{C_{Ag}}{1/k_g + 1/k^{\\prime\\prime}}', 'Ex. 17.3', 'Linear steps in series: resistances add.'],
      ['-r_A^{\\prime\\prime} = k_g C_{Ag} - \\frac{k_g^2}{2k^{\\prime\\prime}}\\left(\\sqrt{1+\\frac{4k^{\\prime\\prime}C_{Ag}}{k_g}}-1\\right)', 'Ex. 17.4', 'Film then second order: no resistance sum.'],
      ['(-r_A)V = (-r_A^\\prime)W = (-r_A^{\\prime\\prime})S', 'Eq. 17.4', 'Put every step on the same basis before combining.'],
      ['\\text{series: } r = r_1 = r_2 = r_3,\\qquad \\text{parallel: } r = r_1 + r_2 + r_3', 'Ch. 17'],
      ['\\frac{C_{Ag}/H - C_e}{-r} = \\frac{1}{Hk_g} + \\frac{1}{k_l} + \\frac{1}{k_s} + \\frac{1}{\\eta k^{\\prime\\prime}}', 'Probs. 17.1, 17.3', 'All-linear chain written in liquid-side units.']
    ],
    method: ['The unknown intermediate concentrations are eliminated by marching backward from a trial surface concentration: the surface rate fixes the flux through the solid film, which fixes the main-liquid concentration, where any bulk reaction adds to the flux through the liquid film, and so on out to the main fluid. The trial value is then adjusted by Brent’s method until the main-fluid concentration matches. This is exact for any rate forms and needs no linearization.', 'Elasticities ∂ln r/∂ln k are central differences at ±1%. For an all-linear series chain they are the shares of the total resistance and add to one; the largest marks the step worth improving.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const C0 = +v.C0;
      const S = solveChain(v, C0), ch = S.ch, H = ch.H;
      const ru = rU(v);
      kpi.push({ k: 'r', l: 'Overall rate −r_A', v: S.J, u: ru, main: true });
      const rxR = (S.r > 0) ? (S.Cs - ch.Ce) / S.r : NaN;
      /* KPIs for each intermediate concentration */
      if (v.g) kpi.push({ k: 'Cgi', l: v.i ? 'C at the gas side of the interface' : 'C beyond the film', v: S.Cgi, u: 'mol/m³' });
      if (v.i) kpi.push({ k: 'Cli', l: 'C on the liquid side of the interface (= C_gi/H)', v: S.Cli, u: 'mol/m³' });
      if (v.l) kpi.push({ k: 'Clb', l: 'C in the main body of liquid', v: S.Clb, u: 'mol/m³' });
      kpi.push({ k: 'Cs', l: 'C at the reacting surface', v: S.Cs, u: 'mol/m³' });
      if (v.b) {
        kpi.push({ k: 'fb', l: 'Share of A reacting in the bulk liquid', v: S.rB / S.J, s: N.fmt(S.rB, 4) + ' ' + ru + ' in the bulk' });
      }
      /* elasticities */
      const base = S.J, el = [];
      const elas = (key, label) => {
        const vp = Object.assign({}, v), vm = Object.assign({}, v);
        vp[key] = +v[key] * 1.01; vm[key] = +v[key] / 1.01;
        const rp = solveChain(vp, C0).J, rm = solveChain(vm, C0).J;
        el.push({ label, e: Math.log(rp / rm) / (2 * Math.log(1.01)) });
      };
      if (v.g) elas('kg', 'gas film k_g');
      if (v.l) elas('kl', 'liquid film k_l');
      if (v.b) elas('kb', 'bulk reaction k_b');
      if (v.s) elas('ks', 'film to solid k_s');
      elas('kr', 'surface reaction k″');
      if (v.i) elas('H', 'partition H');
      const ctrl = el.filter(q => q.label !== 'partition H').reduce((a, b) => (Math.abs(b.e) > Math.abs(a.e) ? b : a), el[0]);
      kpi.push({ k: 'ctrl', l: 'Step most worth improving', v: ctrl.label, s: 'elasticity ' + N.fmt(ctrl.e, 3) + ': 1% better gives ' + N.fmt(ctrl.e, 3) + '% more rate' });
      /* apparent order from the rate curve */
      const Jp = solveChain(v, C0 * 1.01).J, Jm = solveChain(v, C0 / 1.01).J;
      const nApp = Math.log(Jp / Jm) / (2 * Math.log(1.01));
      kpi.push({ k: 'napp', l: 'Apparent order d ln r / d ln C', v: nApp, s: v.law === 'nth' ? 'true surface order ' + fmt(+v.n, 3) : '' });
      kpi.push({ k: 'kov', l: 'Overall coefficient −r / (C₀/H − C_e)', v: S.J / (C0 / H - ch.Ce), u: kU(v) });
      /* closed forms from the book where they exist */
      const onlyFilm = v.g && !v.i && !v.l && !v.b && !v.s;
      const kr = +v.kr * (+v.eta || 1);
      if (onlyFilm && v.law === 'nth' && +v.n === 2) {
        const kg = +v.kg, cf = kg * C0 - kg * kg / (2 * kr) * (Math.sqrt(1 + 4 * kr * C0 / kg) - 1);
        kpi.push({ k: 'closed', l: 'Closed form of Example 17.4', v: cf, u: ru });
      }
      const linear = (v.law === 'rev' || (v.law === 'nth' && +v.n === 1)) && (!v.b || +v.nb === 1);
      if (linear && !v.b) {
        const sumR = (v.g ? 1 / (H * +v.kg) : 0) + (v.l ? 1 / +v.kl : 0) + (v.s ? 1 / +v.ks : 0) + 1 / kr;
        kpi.push({ k: 'closed', l: 'Sum of resistances (all steps linear)', v: (C0 / H - ch.Ce) / sumR, u: ru, s: '1/k_overall = ' + N.fmt(sumR, 4) });
        info.push('All steps are linear, so the resistances 1/(Hk_g), 1/k_l, 1/k_s and 1/(ηk″) simply add (Example 17.3, Problems 17.1 and 17.3).');
      } else if (!linear) info.push('At least one step is nonlinear: resistances cannot be added (Example 17.4). The chain is solved exactly instead; the book’s suggestion of replacing the nonlinear step by an equivalent first-order one is shown in the rate chart as the chord through the operating point.');
      /* table of steps */
      const rows = [];
      let prev = C0 / H;
      if (v.g) { const c = S.Cgi / H; rows.push(['fluid (gas) film', +v.kg, 1 / (H * +v.kg), prev - c, (prev - c) / (C0 / H - ch.Ce)]); prev = c; }
      if (v.l) { const c = S.Clb; rows.push(['liquid film', +v.kl, 1 / +v.kl, prev - c, (prev - c) / (C0 / H - ch.Ce)]); prev = c; }
      if (v.s) { const c = S.Cs; rows.push(['film to the solid', +v.ks, 1 / +v.ks, prev - c, (prev - c) / (C0 / H - ch.Ce)]); prev = c; }
      rows.push(['surface reaction' + (v.law === 'nth' && +v.n !== 1 ? ' (nonlinear)' : ''), kr, linear ? 1 / kr : rxR, S.Cs - ch.Ce, (S.Cs - ch.Ce) / (C0 / H - ch.Ce)]);
      tables.push({ title: 'Where the driving force goes (liquid-side units)', cols: ['step', 'coefficient', 'resistance (drop/flux)', 'concentration drop', 'share of the total drop'], rows });
      /* plot 1: staircase across the films, nodes as categories */
      const U = u => u.replace(/\[t\]/g, v._tu || 's');
      const cats = ['main fluid'], xs = [0], ys = [C0 / H], marks = [];
      let x = 0;
      const film = (cIn, cOut, name) => {
        marks.push({ t: 'band', x0: x + 0.18, x1: x + 0.82, color: 'var(--grid)' });
        marks.push({ t: 'text', x: x + 0.5, y: Math.max(cIn, cOut) + 0.05 * C0 / H, s: name, anchor: 'middle' });
        xs.push(x + 0.18, x + 0.82); ys.push(cIn, cOut); x += 1;
      };
      if (v.g) { film(C0 / H, S.Cgi / H, 'gas film'); cats.push(v.i ? 'interface' : (v.l || v.s ? 'film edge' : 'surface')); xs.push(x); ys.push(S.Cgi / H); }
      if (v.l) { film(S.Cli, S.Clb, 'liquid film'); cats.push('bulk liquid'); xs.push(x); ys.push(S.Clb); }
      if (v.s) { film(S.Clb, S.Cs, 'film to solid'); cats.push('solid surface'); xs.push(x); ys.push(S.Cs); }
      if (cats.length === 1) cats[0] = 'main fluid = surface';
      if (ch.Ce > 0) marks.push({ t: 'h', y: ch.Ce, color: 5, label: 'C_e' });
      plots.push({ title: 'Concentration from the main fluid to the reacting surface', x: { cats }, y: { label: 'C, liquid-side units', unit: 'mol/m³', min: 0, max: C0 / H * 1.15 }, series: [{ name: 'profile', x: xs, y: ys, color: 1, w: 2.4 }], marks, cap: 'Within each film the profile is linear (film theory). Gas-side values are divided by H, which removes the jump at the interface.' });
      /* plot 2: rate vs main concentration */
      const Cs = N.logspace(C0 * 1e-3, C0 * 1e3, 120).filter(c => c / H > ch.Ce * 1.0001);
      const Rv = Cs.map(c => { try { return solveChain(v, c).J; } catch (e) { return NaN; } });
      const ser = [{ name: 'overall rate', x: Cs, y: Rv, color: 1, w: 2.2 }];
      const trans = (v.g ? 1 / (H * +v.kg) : 0) + (v.l ? 1 / +v.kl : 0) + (v.s ? 1 / +v.ks : 0);
      if (trans > 0) ser.push({ name: 'transport alone (surface step infinitely fast)', x: Cs, y: Cs.map(c => (c / H - ch.Ce) / trans), color: 2, dash: '5 4' });
      ser.push({ name: 'surface step alone (no transport resistance)', x: Cs, y: Cs.map(c => ch.rxn(c / H) + (v.b ? ch.rb(c / H) : 0)), color: 3, dash: '2 3' });
      if (!linear) ser.push({ name: 'equivalent first order through this point', x: Cs, y: Cs.map(c => S.J * (c / H - ch.Ce) / (C0 / H - ch.Ce)), color: 4, dash: '8 3 2 3' });
      plots.push({ title: 'Overall rate against the main-fluid concentration', x: { label: 'C in the main fluid', unit: 'mol/m³', log: true }, y: { label: '−r_A', unit: U(ru), log: true }, series: ser, marks: [{ t: 'pt', x: C0, y: S.J, color: 8, label: 'here' }], cap: 'The overall curve follows whichever step is slower; where it bends the apparent order changes.' });
      /* plot 3: elasticities */
      plots.push({ title: 'Payoff from improving each step: ∂ln r / ∂ln k', x: { cats: el.map(q => q.label) }, y: { label: 'elasticity', min: Math.min(0, ...el.map(q => q.e)) * 1.1, max: Math.max(1, ...el.map(q => q.e)) * 1.1 }, series: [{ name: 'elasticity', x: el.map((q, i) => i), y: el.map(q => q.e), type: 'bar', color: 2, legend: false }], marks: el.map((q, i) => ({ t: 'text', x: i, y: Math.max(q.e, 0) + 0.035, s: N.fmt(q.e, 3), anchor: 'middle' })), cap: 'For a chain of linear steps the bars add to one, and each bar is that step’s share of the total resistance.' });
      if (S.Cs < 1e-3 * C0 / H && ch.Ce === 0) info.push('The surface concentration is almost zero: the process is transport-controlled, and improving the surface reaction would not help.');
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* C-1801  effectiveness factor: any kinetics, any shape, heat effects    */
  /* ===================================================================== */
  const familyCache = new Map();
  const traceC = (o, opt) => {
    const key = JSON.stringify([o, opt || {}]);
    if (familyCache.has(key)) return familyCache.get(key);
    const tr = PEL.trace(o, opt);
    if (familyCache.size > 60) familyCache.clear();
    familyCache.set(key, tr);
    return tr;
  };
  CAT.trace = traceC;
  function kinOf(v) {
    const o = { shape: v.shape, kind: v.kind === 'nth' ? 'nth' : v.kind, n: v.kind === 'nth' ? +v.n : 1 };
    if (v.kind === 'rev') o.ce = +v.ce;
    if (v.kind === 'lh' || v.kind === 'lh2') o.K = +v.KC;
    if (v.noniso) { o.beta = +v.beta; o.gamma = +v.gamma; }
    if (v.vol && v.kind === 'nth' && +v.n === 1 && !v.noniso) o.a = (+v.nu - 1) * +v.yA;
    return o;
  }
  const kinLabel = v => (v.kind === 'nth' ? 'order ' + fmt(+v.n, 3) : v.kind === 'rev' ? 'reversible first order' : v.kind === 'lh' ? 'kC/(1+KC), KC_s = ' + fmt(+v.KC, 3) : 'kC/(1+KC)², KC_s = ' + fmt(+v.KC, 3)) + (v.noniso ? ', β = ' + fmt(+v.beta, 3) + ', γ = ' + fmt(+v.gamma, 3) : '');

  ST.tool({
    id: 'thiele', tag: 'C-1801', ch: 18, sec: 'Ch. 18.2–18.4', title: 'Effectiveness factor for any kinetics, shape and heat effect',
    lede: 'The effectiveness factor of a porous pellet — flat plate, cylinder or sphere — for nth-order, reversible or saturation kinetics, found by solving the pellet equation itself rather than reading a chart. It traces the S-shaped curves of hot exothermic pellets with all three steady states, the dead core of orders below one, molar change in the pores, an external film, egg-shell pellets with an inert core, and mixtures of sizes. From physical data or from an observed rate it tells you the regime, the largest pellet still free of pore resistance and how the order and activation energy you measure are falsified.',
    keys: 'effectiveness factor thiele modulus weisz wagner modulus pore diffusion catalyst pellet slab cylinder sphere generalized modulus dead core zero order half order second order reversible langmuir hinshelwood nonisothermal prater weisz hicks multiple steady states volume change egg shell biot film particle size mixture apparent activation energy falsified kinetics problem 18.20 18.21 18.22 18.36 18.38 figure 18.4 18.5 18.6 18.7 18.8',
    units: { t: 's' },
    inputs: [
      { group: 'Pellet and what you know' },
      { k: 'shape', t: 'seg', label: 'Shape', opts: [['slab', 'flat plate'], ['cyl', 'cylinder'], ['sph', 'sphere']], v: 'sph' },
      { k: 'spec', t: 'seg', label: 'Specify by', opts: [['phys', 'size, k, D_e'], ['MT', 'Thiele modulus M_T'], ['MW', 'observed rate (M_W)']], v: 'MT' },
      { k: 'dp', label: v => 'Pellet ' + CAT.sizeWord(v.shape), unit: 'mm', v: 6, min: 1e-6, max: 1e5, scale: 'log', if: v => v.spec === 'phys' },
      { k: 'k3', label: 'Rate constant k‴ (per volume of pellet)', unit: v => (v.kind === 'nth' && +v.n !== 1 ? '(m³/mol)^{' + fmt(+v.n - 1, 3) + '}/[t]' : '1/[t]'), v: 0.1, min: 1e-15, max: 1e12, scale: 'log', if: v => v.spec === 'phys' },
      { k: 'De', label: 'Effective diffusivity D_e', unit: 'm³/(m cat·[t])', v: 1e-6, min: 1e-20, max: 1, scale: 'log', if: v => v.spec === 'phys' },
      { k: 'Cs', label: 'Concentration of A at the pellet surface', unit: 'mol/m³', v: 100, min: 1e-12, max: 1e7, scale: 'log', if: v => v.spec === 'phys' },
      { k: 'MTin', label: v => (v.axis === 'surf' ? 'Thiele modulus on the surface rate, L√(k‴C_sⁿ⁻¹/D_e)' : 'Generalized Thiele modulus M_T (Eq. 18.26)'), v: 1, min: 1e-4, max: 1e4, scale: 'log', if: v => v.spec === 'MT' },
      { k: 'MWin', label: 'Wagner–Weisz–Wheeler modulus M_W = L²(−r‴_obs)/(C_As D_e)', v: 1, min: 1e-6, max: 1e6, scale: 'log', if: v => v.spec === 'MW' },
      { group: 'Kinetics at the surface temperature' },
      { k: 'kind', t: 'select', label: 'Rate form', opts: [['nth', 'nth order: k‴ Cⁿ'], ['rev', 'reversible first order: k‴ (C − C_e)'], ['lh', 'saturation: k‴ C / (1 + K C)'], ['lh2', 'dual site: k‴ C / (1 + K C)²']], v: 'nth' },
      { k: 'n', label: 'Order n', v: 1, min: 0, max: 5, step: 0.25, if: v => v.kind === 'nth' },
      { k: 'ce', label: 'C_e / C_As (equilibrium over surface value)', v: 0.3, min: 0, max: 0.999, if: v => v.kind === 'rev' },
      { k: 'KC', label: 'K C_As (adsorption strength at the surface)', v: 5, min: 0, max: 1e4, scale: 'log', if: v => v.kind === 'lh' || v.kind === 'lh2' },
      { group: 'Heat effects inside the pellet (Prater, Weisz–Hicks)' },
      { k: 'noniso', t: 'check', label: 'Temperature varies inside the pellet', v: false },
      { k: 'beta', label: 'Prater number β = D_e C_As(−ΔH)/(k_eff T_s)', v: 0.4, min: -0.9, max: 2, step: 0.05, if: v => v.noniso, hint: 'Positive for exothermic. ΔT_max inside the pellet = β T_s (Eq. 18.38).' },
      { k: 'gamma', label: 'Arrhenius number γ = E/(R T_s)', v: 20, min: 0, max: 80, if: v => v.noniso },
      { group: 'More physics', if: v => true },
      { k: 'vol', t: 'check', label: 'Molar change in the pores, A → νR (first order, isothermal)', v: false },
      { k: 'nu', label: 'ν (moles of product per mole of A)', v: 2, min: 0.05, max: 10, if: v => v.vol },
      { k: 'yA', label: 'Mole fraction of A at the pellet surface', v: 1, min: 0.001, max: 1, if: v => v.vol },
      { k: 'film', t: 'check', label: 'Film around the pellet (isothermal nth order)', v: false },
      { k: 'Bi', label: 'Mass Biot number k_g L / D_e', v: 10, min: 1e-4, max: 1e6, scale: 'log', if: v => v.film },
      { k: 'shell', t: 'check', label: 'Egg-shell pellet: only an outer layer is active', v: false },
      { k: 'delta', label: 'Active layer thickness / R (half thickness for a plate)', v: 0.1, min: 0.001, max: 1, if: v => v.shell },
      { k: 'mix', t: 'check', label: 'Mixture of pellet sizes (Aris)', v: false },
      { k: 'mixT', t: 'table', label: 'Size relative to the pellet above, volume fraction', cols: ['d / d_ref', 'volume fraction'], v: [[0.5, 0.2], [1, 0.5], [2, 0.3]], if: v => v.mix },
      { group: 'Falsified kinetics under strong pore resistance' },
      { k: 'Et', label: 'True activation energy E', unit: 'kJ/mol', v: 120, min: 0, max: 1000 },
      { k: 'Ed', label: 'Activation energy of diffusion E_diff', unit: 'kJ/mol', v: 10, min: 0, max: 200 },
      { group: 'Chart' },
      { k: 'axis', t: 'seg', label: 'Modulus on the axis (and as input)', opts: [['gen', 'generalized M_T (Eq. 18.26)'], ['surf', 'L√(k‴C_sⁿ⁻¹/D_e) at the surface']], v: 'gen', hint: 'The generalized modulus makes every curve end on ℰ = 1/M_T; the surface modulus is the one Weisz and Hicks used (their φ = 3× this for a sphere).' },
      { k: 'fam', t: 'select', label: 'Family of curves to draw', opts: [['shapes', 'three shapes, first order (Fig. 18.6)'], ['orders', 'orders 0, ½, 1, 2, 3 (generalized modulus)'], ['beta', 'heat effects: β from −0.4 to 0.8 (Fig. 18.8)'], ['vol', 'molar change: ν = ½, 1, 2 (Fig. 18.5)'], ['none', 'current kinetics only']], v: 'shapes' }
    ],
    presets: [
      { name: 'Figure 18.4 — single pore, mL = 2', v: { shape: 'slab', spec: 'MT', MTin: 2, kind: 'nth', n: 1, noniso: false, vol: false, film: false, shell: false, mix: false, fam: 'shapes' }, book: { eta: { v: 0.48, tol: 0.01 } }, note: 'A single pore with a sealed end is a flat plate: η = tanh(mL)/mL. The book reads 0.92, 0.76, 0.48, 0.20 and 0.10 for mL = 0.5, 1, 2, 5 and 10.' },
      { name: 'Figure 18.6 — first order sphere at M_T = 1', v: { shape: 'sph', spec: 'MT', MTin: 1, kind: 'nth', n: 1, noniso: false, vol: false, film: false, shell: false, mix: false, fam: 'shapes' }, book: { eta: { v: 0.6716, tol: 0.001 } }, note: 'Eq. 18.22: η = (1/M_T)(1/tanh 3M_T − 1/3M_T).' },
      { name: 'Problem 18.21 — second order, 6 mm sphere', v: { shape: 'sph', spec: 'phys', dp: 6, k3: 0.1, De: 1e-6, Cs: 100, kind: 'nth', n: 2, noniso: false, vol: false, film: false, shell: false, mix: false, fam: 'orders' }, note: 'd_p = 6 mm, D_e = 10⁻⁶ m³/m cat·s, C_A = 100 mol/m³, −r‴ = 0.1 C² mol/m³ cat·s. M_T by Eq. 18.28 is 3.87; the exact second-order sphere gives the slowdown.' },
      { name: 'Problem 18.20 — which ball size to order', v: { shape: 'sph', spec: 'phys', dp: 6, k3: 0.88, De: 2e-6, Cs: 200, kind: 'nth', n: 1, noniso: false, vol: false, film: false, shell: false, mix: false, fam: 'shapes' }, note: '−r‴ = 0.88 C_A mol/m³ cat·s with D_e = 2 × 10⁻⁶ m³/m cat·s: the largest balls still free of pore resistance (M_T = 0.4) are read from “largest diffusion-free size”.' },
      { name: 'Problem 18.22 — size that keeps pores out of it', v: { shape: 'sph', spec: 'phys', dp: 5, k3: 0.1, De: 1e-7, Cs: 10, kind: 'nth', n: 1, noniso: false, vol: false, film: false, shell: false, mix: false, fam: 'shapes' }, note: '−r‴ = 10⁻⁶ mol/cm³ cat·s at C_A = 10⁻⁵ mol/cm³, so k‴ = 0.1 s⁻¹; D_e = 10⁻³ cm³/cm cat·s = 10⁻⁷ m³/m cat·s.' },
      { name: 'Half order with a dead core', v: { shape: 'sph', spec: 'MT', MTin: 4, kind: 'nth', n: 0.5, noniso: false, vol: false, film: false, shell: false, mix: false, fam: 'orders' }, note: 'Below first order the reactant runs out before the centre: the inner part of the pellet does nothing.' },
      { name: 'Hot pellet with three steady states (Weisz–Hicks)', v: { shape: 'sph', spec: 'MT', axis: 'surf', MTin: 0.2189, kind: 'nth', n: 1, noniso: true, beta: 0.4, gamma: 20, vol: false, film: false, shell: false, mix: false, fam: 'beta' }, note: 'β = 0.4, γ = 20 and Weisz–Hicks φ = 3 × 0.2189 = 0.657: in a narrow band of modulus the exothermic pellet can be cold, ignited, or balanced on an unstable middle state (the hysteresis of Fig. 18.8). The three effectiveness factors are 1.40, 3.64 and 9.98.' },
      { name: 'Endothermic pellet', v: { shape: 'sph', spec: 'MT', axis: 'gen', MTin: 1, kind: 'nth', n: 1, noniso: true, beta: -0.2, gamma: 20, vol: false, film: false, shell: false, mix: false, fam: 'beta' }, note: 'The pellet interior cools, and η falls below the isothermal value.' },
      { name: 'Figure 18.5 — volume doubles in the pores', v: { shape: 'slab', spec: 'MT', MTin: 2, kind: 'nth', n: 1, noniso: false, vol: true, nu: 2, yA: 1, film: false, shell: false, mix: false, fam: 'vol' }, note: 'Product molecules crowd out of the pores and hinder reactant getting in, so η drops; contraction does the opposite (Thiele 1939).' },
      { name: 'Problem 18.38 — egg-shell platinum (0.3 mm on 6 mm)', v: { shape: 'sph', spec: 'MT', MTin: 10, kind: 'nth', n: 1, noniso: false, vol: false, film: false, shell: true, delta: 0.1, mix: false, fam: 'shapes' }, note: 'Under strong pore resistance only the outside of the pellet works anyway: impregnating just the outer 0.3 mm of a 6 mm ball saves most of the platinum at almost the same rate.' },
      { name: 'Film plus pores', v: { shape: 'sph', spec: 'MT', MTin: 3, kind: 'nth', n: 1, noniso: false, vol: false, film: true, Bi: 5, shell: false, mix: false, fam: 'shapes' }, note: 'With an external film the surface concentration falls below the bulk value and the overall effectiveness drops further.' },
      { name: 'Reversible first order (Eq. 18.27)', v: { shape: 'sph', spec: 'MT', MTin: 2, kind: 'rev', ce: 0.5, noniso: false, vol: false, film: false, shell: false, mix: false, fam: 'none' }, note: 'With the generalized modulus the reversible reaction falls on the first-order curve.' },
      { name: 'Strongly adsorbed reactant', v: { shape: 'slab', spec: 'MT', MTin: 2, kind: 'lh', KC: 50, noniso: false, vol: false, film: false, shell: false, mix: false, fam: 'orders' }, note: 'Saturation kinetics behave nearly zero order near the surface: compare with the order-0 curve.' },
      { name: 'Example 18.1 pellet from its Weisz modulus', v: { shape: 'sph', spec: 'MW', MWin: 16, kind: 'nth', n: 1, noniso: false, vol: false, film: false, shell: false, mix: false, fam: 'shapes' }, note: 'M_W = 16 from the observed rate: strong pore diffusion, M_T ≈ 16.3 and η ≈ 0.06 for a first-order sphere.' }
    ],
    eqs: [
      ['\\mathscr{E}_{plate} = \\frac{\\tanh M_T}{M_T},\\qquad \\mathscr{E}_{cyl} = \\frac{1}{M_T}\\frac{I_1(2M_T)}{I_0(2M_T)}', 'Eqs. 18.20–18.21'],
      ['\\mathscr{E}_{sphere} = \\frac{1}{M_T}\\left(\\frac{1}{\\tanh 3M_T}-\\frac{1}{3M_T}\\right)', 'Eq. 18.22'],
      ['M_T = L\\sqrt{k^{\\prime\\prime\\prime}/\\mathscr{D}_e},\\qquad L = \\frac{V_p}{S_{ex}} = \\frac{\\text{thickness}}{2},\\ \\frac{R}{2},\\ \\frac{R}{3}', 'Eqs. 18.13, 18.23'],
      ['M_T = \\frac{(-r_{As}^{\\prime\\prime\\prime})L}{\\left[2\\mathscr{D}_e\\int_{C_{Ae}}^{C_{As}}(-r_A^{\\prime\\prime\\prime})\\,dC_A\\right]^{1/2}}', 'Eq. 18.26', 'Generalized modulus: every rate form follows the first-order curve closely.'],
      ['M_T = L\\sqrt{\\frac{(n+1)k^{\\prime\\prime\\prime}C_{As}^{n-1}}{2\\mathscr{D}_e}}\\ \\text{(nth order)},\\qquad M_T = L\\sqrt{\\frac{k^{\\prime\\prime\\prime}}{\\mathscr{D}_e X_{Ae}}}\\ \\text{(reversible)}', 'Eqs. 18.27–18.28'],
      ['M_W = M_T^2\\,\\mathscr{E} = \\frac{L^2(-r_A^{\\prime\\prime\\prime})_{obs}}{C_{As}\\mathscr{D}_e}', 'Eq. 18.24', 'Diffusion-free below 0.15, strong pore resistance above 4.'],
      ['\\frac{1}{\\zeta^s}\\frac{d}{d\\zeta}\\left(\\zeta^s\\frac{dc}{d\\zeta}\\right) = f(c)\\,e^{\\gamma\\beta(1-c)/(1+\\beta(1-c))}', 'Pellet equation', 'Slab, cylinder, sphere: s = 0, 1, 2; c = C/C_As.'],
      ['\\frac{T}{T_s} = 1+\\beta(1-c),\\qquad \\beta = \\frac{\\mathscr{D}_e C_{As}(-\\Delta H_r)}{k_{eff}T_s}', 'Eqs. 18.37–18.38', 'Prater: temperature and concentration profiles have the same shape.'],
      ['n_{obs} = \\frac{n+1}{2},\\qquad \\mathbf{E}_{obs} = \\frac{\\mathbf{E}_{true}+\\mathbf{E}_{diff}}{2}', 'Eqs. 18.30, 18.32', 'Strong pore resistance.'],
      ['\\frac{r_1}{r_2} = \\frac{R_2}{R_1}\\ \\text{(strong)},\\quad 1\\ \\text{(diffusion-free)};\\qquad \\bar{\\mathscr{E}} = \\mathscr{E}_1 f_1 + \\mathscr{E}_2 f_2 + \\cdots', 'Eq. 18.25, Aris']
    ],
    method: ['The pellet equation is integrated outward from the centre with the centre concentration as the parameter, until c = 1 marks the surface; the distance travelled is the modulus and the slope there is η. Sweeping the centre value traces the whole η–M_T curve, including the folds of nonisothermal pellets, and adaptive refinement keeps neighbouring points within 0.035 decades. Values at a requested modulus are then polished by Brent’s method on the centre value, so they are exact to the integrator tolerance (10⁻¹⁰) — for example 0.671636 for the first-order sphere at M_T = 1.', 'Orders below one are continued past the onset of a dead core by starting from its edge with the exact local power-law behaviour. Beyond the traced range η is extended along η·M_T → 1 with a 1/M_T correction fitted to the last point.', 'Egg-shell pellets are solved at fixed modulus with an inert core of uniform concentration. The film is coupled by solving the flux balance k_g(C_b − C_s) = L k C_sⁿ η(M_s) for C_s. Apparent order and activation energy come from the local slope m = d ln η / d ln M_T: E_obs = E(1 + m/2) − E_diff m/2, which tends to Eq. 18.32 as m → −1.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const o = kinOf(v);
      if (v.vol && !(v.kind === 'nth' && +v.n === 1 && !v.noniso)) warn.push('Molar change is handled for isothermal first-order reactions only; it is ignored here.');
      const tr = traceC(o);
      if (!tr.pts.length) throw new Error('Could not trace the pellet equation for these kinetics.');
      const g = tr.g, s = tr.s, sh = v.shape, tu = v._tu || 's';
      /* the modulus on the surface rate (MT0 = L sqrt(r_s/(D C_s))) and the book's generalized one */
      let MT0 = NaN, L = NaN, rs = NaN;
      if (v.spec === 'phys') {
        L = CAT.Lof(sh, +v.dp / 1000);
        const Cs = +v.Cs, k = +v.k3;
        rs = v.kind === 'nth' ? k * Math.pow(Cs, +v.n) : v.kind === 'rev' ? k * Cs * (1 - +v.ce) : v.kind === 'lh' ? k * Cs / (1 + +v.KC) : k * Cs / Math.pow(1 + +v.KC, 2);
        MT0 = L * Math.sqrt(rs / (+v.De * Cs));
      } else if (v.spec === 'MT') MT0 = v.axis === 'surf' ? +v.MTin : +v.MTin / g;
      let sols = [];
      if (v.spec === 'MW') {
        const MW = +v.MWin;
        const fromMW = PEL.findAll(tr, p => Math.log(p.eta * p.MT * p.MT / MW));
        if (!fromMW.length) {
          const last = tr.pts[tr.pts.length - 1];
          /* beyond the traced range: η = 1/(g M), M_W = η M² ⇒ M = g M_W */
          if (MW > last.eta * last.MT * last.MT) sols = [{ MT: g * MW, eta: 1 / (g * g * MW), asym: true }];
          else sols = [{ eta: 1, MT: Math.sqrt(MW) }];
        } else sols = fromMW;
        MT0 = sols[0].MT;
        kpi.push({ k: 'MW', l: 'M_W (given)', v: MW, s: CAT.regimeW(MW), state: CAT.state(MW < 0.15 ? 0 : MW > 4 ? 9 : 1) });
      } else {
        sols = PEL.etaAt(tr, MT0, true).map(e => Object.assign(e, { MT: MT0 }));
      }
      const Mgen = MT0 * g;
      sols.sort((a, b) => a.eta - b.eta);
      const main = sols[0];
      const eta = main.eta;
      kpi.push({ k: 'eta', l: sols.length > 1 ? 'Effectiveness factor ℰ (cold branch)' : 'Effectiveness factor ℰ', v: eta, main: true, s: CAT.regime(Mgen), state: CAT.state(Mgen) });
      if (sols.length > 1) {
        kpi.push({ k: 'etaMid', l: 'Middle steady state (unstable)', v: sols[1].eta, state: 'warn' });
        kpi.push({ k: 'etaHi', l: 'Ignited steady state', v: sols[sols.length - 1].eta, state: 'bad', s: 'centre ' + fmt(1 + (+v.beta) * (1 - (sols[sols.length - 1].c0 || 0)), 4) + ' T_s' });
        warn.push('Three steady states at this modulus: which one the pellet sits on depends on its history (Fig. 18.8 hysteresis).');
      }
      kpi.push({ k: 'MT', l: 'Thiele modulus M_T (generalized, Eq. 18.26)', v: Mgen, s: g !== 1 ? 'L√(k/D_e) on the surface rate = ' + fmt(MT0, 4) : '' });
      kpi.push({ k: 'MWo', l: 'Weisz modulus M_W = M_T² ℰ', v: eta * MT0 * MT0, s: CAT.regimeW(eta * MT0 * MT0) });
      if (main.core !== undefined && main.core > 0) kpi.push({ k: 'core', l: 'Dead core (fraction of the radius)', v: main.core, s: N.fmt(Math.pow(main.core, s + 1) * 100, 3) + '% of the pellet volume idle' });
      else if (main.c0 !== undefined) kpi.push({ k: 'c0', l: 'Concentration at the centre, C_center/C_As', v: main.c0 });
      if (v.noniso) kpi.push({ k: 'Tc', l: 'Centre temperature T_center/T_s (Prater)', v: 1 + (+v.beta) * (1 - (main.c0 || 0)) });
      if (v.spec === 'phys') {
        const robs = eta * rs;
        kpi.push({ k: 'robs', l: 'Rate per pellet volume −r‴_A', v: robs, u: 'mol/(m³ cat·' + tu + ')', s: 'without pore resistance ' + fmt(rs, 4) });
        const dFree = +v.dp * 0.4 / Mgen, dStrong = +v.dp * 4 / Mgen;
        kpi.push({ k: 'dpFree', l: 'Largest diffusion-free ' + CAT.sizeWord(sh) + ' (M_T = 0.4)', v: dFree, u: 'mm' });
        kpi.push({ k: 'dpStrong', l: CAT.sizeWord(sh)[0].toUpperCase() + CAT.sizeWord(sh).slice(1) + ' where pore resistance becomes strong (M_T = 4)', v: dStrong, u: 'mm' });
      }
      /* local slope, apparent order and activation energy */
      const etaAtGen = Mg => PEL.etaAt(tr, Mg / g)[0].eta;
      const h = 0.01, m = Math.log(etaAtGen(Mgen * (1 + h)) / etaAtGen(Mgen / (1 + h))) / (2 * Math.log(1 + h));
      const Et = +v.Et, Ed = +v.Ed;
      const Eobs = Et * (1 + m / 2) - Ed * m / 2;
      kpi.push({ k: 'Eobs', l: 'Activation energy you would measure', v: Eobs, u: 'kJ/mol', s: 'd ln ℰ/d ln M_T = ' + fmt(m, 3) + ' (strong limit: (E + E_diff)/2 = ' + fmt((Et + Ed) / 2, 4) + ')' });
      if (v.kind === 'nth' && !v.noniso) kpi.push({ k: 'nobs', l: 'Order you would measure', v: +v.n + m * (+v.n - 1) / 2, s: 'strong limit (n + 1)/2 = ' + fmt((+v.n + 1) / 2, 3) });
      /* film */
      if (v.film) {
        if (v.kind !== 'nth' || v.noniso) warn.push('The film coupling is for isothermal nth-order kinetics; it is skipped here.');
        else {
          const F = PEL.withFilm(tr, +v.n, MT0, +v.Bi);
          kpi.push({ k: 'etaO', l: 'Overall effectiveness (film + pores), on the bulk concentration', v: F.etaO, s: 'C_s/C_bulk = ' + fmt(F.cs, 4) + ', pores alone at C_s: ' + fmt(F.etaIn, 4) });
          if (F.cs < 0.9) info.push('The film takes ' + fmt((1 - F.cs) * 100, 3) + '% of the concentration before the reactant reaches the pellet.');
        }
      }
      /* egg shell */
      if (v.shell) {
        if (v.noniso || !(v.kind === 'nth' || v.kind === 'rev' || v.kind === 'lh' || v.kind === 'lh2')) warn.push('Egg-shell pellets are solved for isothermal kinetics only.');
        else {
          const rho = 1 - +v.delta;
          const E = PEL.solvePhi(o, MT0, { rhoC: rho });
          const frac = 1 - Math.pow(rho, s + 1);
          kpi.push({ k: 'shellRate', l: 'Rate of the egg-shell pellet / uniformly loaded pellet', v: E.eta / eta, s: 'same k‴ in the active layer; ' + fmt(frac * 100, 3) + '% of the volume loaded' });
          kpi.push({ k: 'saved', l: 'Active material saved', v: 1 - frac, u: '', s: 'ℰ on the active layer: ' + fmt(E.etaActive, 4) });
        }
      }
      /* mixture of sizes */
      if (v.mix) {
        const rows = (v.mixT || []).filter(r => r && r.length >= 2 && +r[0] > 0 && +r[1] > 0);
        if (rows.length) {
          const fs = rows.map(r => +r[1]), ftot = N.sum(fs);
          const etas = rows.map(r => etaAtGen(Mgen * +r[0]));
          const ebar = N.sum(etas.map((e, i) => e * fs[i])) / ftot;
          kpi.push({ k: 'etaMix', l: 'Mean effectiveness of the mixture ℰ̄ (Aris)', v: ebar });
          tables.push({ title: 'Mixture of sizes', cols: ['d / d_ref', 'volume fraction', 'M_T', 'ℰ'], rows: rows.map((r, i) => [+r[0], fs[i] / ftot, Mgen * +r[0], etas[i]]) });
        }
      }
      /* three shapes at the same modulus (first order) */
      tables.push({ title: 'Same modulus, three shapes (first order, exact)', cols: ['shape', 'M_T', 'ℰ', 'regime'], rows: ['slab', 'cyl', 'sph'].map(q => [PEL.NAME[q], Mgen, PEL.eta1(q, Mgen), CAT.regime(Mgen)]) });
      /* ----- chart 1: η vs M_T family ----- */
      const ser = [];
      const grid = N.logspace(0.05, 50, 160);
      const ax = v.axis === 'surf' ? (p, trq) => p.MT : (p, trq) => p.MT * trq.g;
      const axCur = M => (v.axis === 'surf' ? M : M * g);
      const add = (name, trq, color, dash) => { const P = trq.pts; ser.push({ name, x: P.map(p => ax(p, trq)), y: P.map(p => p.eta), color, dash, w: 1.6 }); };
      if (v.fam === 'shapes') ['slab', 'cyl', 'sph'].forEach((q, i) => ser.push({ name: PEL.NAME[q] + ', first order', x: grid, y: grid.map(M => PEL.eta1(q, M)), color: i + 2, w: 1.5 }));
      else if (v.fam === 'orders') [0, 0.5, 1, 2, 3].forEach((n, i) => add('n = ' + fmt(n, 2), traceC({ shape: sh, kind: 'nth', n }, { maxPts: 140, dl: 0.07 }), i + 2));
      else if (v.fam === 'beta') [-0.4, -0.2, 0.2, 0.4, 0.6, 0.8].forEach((b, i) => add('β = ' + fmt(b, 2), traceC({ shape: sh, kind: 'nth', n: 1, beta: b, gamma: v.noniso ? +v.gamma : 20 }, { maxPts: 150, dl: 0.07 }), (i % 7) + 2));
      else if (v.fam === 'vol') [0.5, 1, 2].forEach((nu, i) => add('V_out/V_in = ' + fmt(nu, 2), traceC({ shape: sh, kind: 'nth', n: 1, a: nu - 1 }, { maxPts: 140, dl: 0.07 }), i + 2));
      if (v.fam === 'beta' || v.fam === 'shapes') { const gq = N.logspace(1e-3, 1e3, 200); if (v.fam === 'beta') ser.push({ name: 'β = 0 (isothermal)', x: gq, y: gq.map(M => PEL.eta1(sh, M)), color: 'var(--muted)', w: 1.4 }); }
      add('this pellet: ' + kinLabel(v), tr, 1);
      ser[ser.length - 1].w = 2.6;
      ser.push({ name: 'ℰ = 1/M_T', x: [1e-3, 1e3], y: [1e3, 1e-3], color: 'var(--muted)', dash: '4 4', legend: false });
      const marks = [{ t: 'v', x: 0.4, label: '0.4', color: 'var(--muted)' }, { t: 'v', x: 4, label: '4', color: 'var(--muted)' }];
      sols.forEach((q, i) => marks.push({ t: 'pt', x: axCur(q.MT || MT0), y: q.eta, color: 8, label: sols.length > 1 ? ['cold', 'unstable', 'ignited'][Math.min(i, 2)] : 'here', open: sols.length > 1 && i === 1 }));
      const xc = axCur(MT0), allX = ser.filter(q => q.name && q.name.indexOf('β') === 0).map(q => Math.min(...q.x.filter(x => x > 0)));
      const xmin = Math.max(1e-3, Math.min(0.05, xc / 5, ...allX.map(x => x * 2))), xmax = Math.max(50, xc * 5);
      const ymax = Math.max(1.6, ...sols.map(q => q.eta * 2), ...ser.filter(q => q.name && q.name.indexOf('β') === 0).map(q => Math.max(...q.y.filter(isFinite)) * 1.3));
      plots.push({ title: 'Effectiveness factor against the ' + (v.axis === 'surf' ? 'Thiele modulus on the surface rate' : 'generalized Thiele modulus'), x: { label: v.axis === 'surf' ? 'L√(k‴C_sⁿ⁻¹/D_e)' : 'M_T (Eq. 18.26)', log: true, min: xmin, max: xmax }, y: { label: 'ℰ', log: true, min: 0.02, max: Math.min(ymax, 200) }, series: ser, marks, hmax: 470, cap: 'Left of 0.4 the pellet is free of pore resistance; right of 4 every curve falls as 1/M_T and only the outer skin works.' });
      /* ----- chart 2: profile ----- */
      const pser = [];
      sols.forEach((q, i) => {
        const pr = v.spec === 'MW' && q.q === undefined && q.zc === undefined ? null : PEL.profile(tr, q, q.MT || MT0);
        if (pr) pser.push({ name: (sols.length > 1 ? ['cold', 'unstable', 'ignited'][Math.min(i, 2)] + ': ' : '') + 'C/C_As', x: pr.x, y: pr.c, color: i + 1, dash: sols.length > 1 && i === 1 ? '5 4' : undefined, w: 2 });
        if (pr && v.noniso) pser.push({ name: (sols.length > 1 ? ['cold', 'unstable', 'ignited'][Math.min(i, 2)] + ': ' : '') + '(T − T_s)/(β T_s)', x: pr.x, y: pr.c.map(c => 1 - c), color: i + 4, dash: '2 3' });
      });
      if (pser.length) plots.push({ title: 'Inside the pellet', x: { label: v.shape === 'slab' ? 'x / half-thickness' : 'r / R', min: 0, max: 1 }, y: { label: 'C/C_As', min: 0, max: 1.02 }, series: pser, marks: main.core > 0 ? [{ t: 'band', x0: 0, x1: main.core, color: 'var(--grid)' }, { t: 'text', x: main.core / 2, y: 0.5, s: 'dead core', anchor: 'middle' }] : [], cap: v.noniso ? 'Prater: the temperature profile is the mirror image of the concentration profile, T − T_s = β T_s (1 − C/C_As).' : 'Centre on the left, surface on the right.' });
      /* ----- chart 3: η against M_W ----- */
      const P = tr.pts;
      plots.push({ title: 'Effectiveness against the observable modulus M_W (Fig. 18.7)', x: { label: 'M_W = M_T² ℰ', log: true, min: 0.01, max: 1e3 }, y: { label: 'ℰ', log: true, min: 0.01, max: Math.min(ymax, 200) }, series: [{ name: kinLabel(v), x: P.map(p => p.eta * p.MT * p.MT), y: P.map(p => p.eta), color: 1, w: 2.2 }], marks: [{ t: 'v', x: 0.15, label: '0.15', color: 'var(--muted)' }, { t: 'v', x: 4, label: '4', color: 'var(--muted)' }].concat(sols.map(q => ({ t: 'pt', x: q.eta * Math.pow(q.MT || MT0, 2), y: q.eta, color: 8 }))), cap: 'Use this when you have measured a rate but do not know k‴: every quantity in M_W is observable.' });
      /* ----- chart 4: rate vs particle size ----- */
      const xsz = P.map(p => p.MT / MT0), ysz = P.map(p => p.eta);
      plots.push({ title: 'Rate per unit catalyst against pellet size', x: { label: 'size / present size', log: true, min: 0.01, max: 100 }, y: { label: 'rate / diffusion-free rate', log: true, min: 0.005, max: Math.min(ymax, 200) }, series: [{ name: 'ℰ at that size', x: xsz, y: ysz, color: 2, w: 2.2 }], marks: [{ t: 'pt', x: 1, y: eta, color: 8, label: 'present' }], cap: 'Diffusion-free: rate independent of size. Strong resistance: rate ∝ 1/size (Eq. 18.25).' });
      if (Mgen > 4) info.push('Strong pore resistance: the reaction is confined to an outer skin; smaller pellets (or an egg-shell loading) would use the catalyst better.');
      if (Mgen < 0.4) info.push('Free of pore resistance: larger pellets would cut the pressure drop without losing rate, up to about ' + (v.spec === 'phys' ? fmt(+v.dp * 0.4 / Mgen, 3) + ' mm.' : '0.4/M_T times the present size.'));
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* C-1802  which resistance controls? single rate or a set of runs        */
  /* ===================================================================== */
  ST.tool({
    id: 'diagnose', tag: 'C-1802', ch: 18, sec: 'Ch. 18.6', title: 'Which resistance controls: film, pores or heat',
    lede: 'Interpret catalytic rate measurements the way Chapter 18 does. From one observed rate and the pellet’s properties: the film-to-reaction ratio, the Weisz modulus, the true rate constant hidden behind the pores, and the temperature rise across the film and inside the pellet. From a set of runs at different pellet sizes, stirring or recycle rates and temperatures: whether the film matters, where each pellet sits on the effectiveness curve, and the true activation energy.',
    keys: 'controlling resistance diagnosis film mass transfer pore diffusion weisz modulus heat effects temperature rise prater basket reactor spinning rate recycle rate particle size activation energy example 18.1 18.7 problem 18.23 18.24 18.25 18.26 18.27 mears froessling ranz',
    units: { t: 'h' },
    inputs: [
      { group: 'What you have' },
      { k: 'mode', t: 'seg', label: 'Data', opts: [['one', 'one rate and the pellet properties'], ['runs', 'runs at different sizes, stirring, T']], v: 'one' },
      { k: 'shape', t: 'seg', label: 'Pellet shape', opts: [['slab', 'flat plate'], ['cyl', 'cylinder'], ['sph', 'sphere']], v: 'sph' },
      { k: 'n', label: 'Order of the surface reaction', v: 1, min: 0, max: 4, step: 0.5 },
      { group: 'Observed rate and properties', if: v => v.mode === 'one' },
      { k: 'dp', label: v => 'Pellet ' + CAT.sizeWord(v.shape), unit: 'mm', v: 2.4, min: 1e-6, max: 1e5, scale: 'log', if: v => v.mode === 'one' },
      { k: 'robs', label: 'Observed rate −r‴_A,obs (per pellet volume)', unit: 'mol/(m³ cat·[t])', v: 1e5, min: 0, max: 1e15, if: v => v.mode === 'one' },
      { k: 'Cg', label: 'Concentration of A in the gas, C_Ag', unit: 'mol/m³', v: 20, min: 1e-12, max: 1e7, if: v => v.mode === 'one' },
      { k: 'De', label: 'Effective diffusivity D_e', unit: 'm³/(m cat·[t])', v: 5e-5, min: 1e-20, max: 1e3, scale: 'log', if: v => v.mode === 'one' },
      { k: 'kg', label: 'Film mass transfer coefficient k_g', unit: 'm³/(m² cat·[t])', v: 300, min: 1e-12, max: 1e9, scale: 'log', if: v => v.mode === 'one' },
      { k: 'h', label: 'Film heat transfer coefficient h', unit: 'kJ/(m² cat·[t]·K)', v: 160, min: 1e-12, max: 1e9, scale: 'log', if: v => v.mode === 'one' },
      { k: 'keff', label: 'Effective conductivity of the pellet k_eff', unit: 'kJ/(m cat·[t]·K)', v: 1.6, min: 1e-12, max: 1e9, scale: 'log', if: v => v.mode === 'one' },
      { k: 'dH', label: 'Heat of reaction ΔH_r', unit: 'kJ/mol', v: -160, min: -1e4, max: 1e4, if: v => v.mode === 'one' },
      { k: 'Ts', label: 'Gas temperature', unit: '°C', v: 336, min: -273, max: 3000, if: v => v.mode === 'one' },
      { k: 'E', label: 'Activation energy (for the Arrhenius number)', unit: 'kJ/mol', v: 100, min: 0, max: 1000, if: v => v.mode === 'one' },
      { group: 'Runs', if: v => v.mode === 'runs' },
      { k: 'dtype', t: 'seg', label: 'Each run gives', opts: [['rate', 'a rate at a concentration'], ['mfr', 'W, C_A0, v, X_A in mixed flow']], v: 'rate', if: v => v.mode === 'runs' },
      { k: 'runsR', t: 'table', label: 'Runs: size, agitation (1 low … 3 high), T, C_A, −r_A', cols: ['d_p', 'agitation', 'T, K', 'C_A', '−r_A'], v: [[1, 3, 500, 1, 3], [3, 1, 500, 1, 1], [3, 3, 500, 1, 1]], if: v => v.mode === 'runs' && v.dtype === 'rate' },
      { k: 'runsM', t: 'table', label: 'Runs: size, agitation, T, W, C_A0, v, X_A', cols: ['d_p', 'agitation', 'T, K', 'W', 'C_A0', 'v', 'X_A'], v: [[3, 3, 500, 1, 100, 9, 0.4], [12, 3, 500, 4, 300, 8, 0.6]], if: v => v.mode === 'runs' && v.dtype === 'mfr' },
      { k: 'eps', label: 'ε_A (expansion) for the mixed-flow runs', v: 0, min: -0.99, max: 10, if: v => v.mode === 'runs' && v.dtype === 'mfr' },
      { k: 'Ed', label: 'Activation energy of diffusion (for correcting E)', unit: 'kJ/mol', v: 0, min: 0, max: 100, if: v => v.mode === 'runs' }
    ],
    presets: [
      { name: 'Example 18.1 — search for the controlling mechanism', v: { mode: 'one', shape: 'sph', n: 1, dp: 2.4, robs: 1e5, Cg: 20, De: 5e-5, kg: 300, h: 160, keff: 1.6, dH: -160, Ts: 336, E: 100, _tu: 'h' }, book: { film: { v: 1 / 150, tol: 0.01 }, MW: { v: 16, tol: 0.01 }, dTp: { v: 0.1, tol: 0.01 }, dTf: { v: 40, tol: 0.01 } }, note: 'd_p = 2.4 mm, D_e = 5 × 10⁻⁵ m³/hr·m cat, k_eff = 1.6 kJ/hr·m·K, h = 160 kJ/hr·m²·K, k_g = 300 m³/hr·m², ΔH = −160 kJ/mol, C_Ag = 20 mol/m³, −r‴ = 10⁵ mol/hr·m³ cat. Book: film ratio 1/150, M_W = 16, ΔT pellet 0.1 °C, ΔT film 40 °C.' },
      { name: 'Example 18.7 — basket reactor runs', v: { mode: 'runs', shape: 'sph', n: 1, dtype: 'rate', runsR: [[1, 3, 500, 1, 3], [3, 1, 500, 1, 1], [3, 3, 500, 1, 1]], Ed: 0 }, book: { film: 'negligible', pore: 'strong' }, note: 'Pellets of size 1 and 3, the larger ones at two spinning rates, same leaving concentration. Book: negligible film resistance, strong pore diffusion.' },
      { name: 'Problem 18.23 — mixed flow runs', v: { mode: 'runs', shape: 'sph', n: 1, dtype: 'mfr', runsM: [[3, 3, 500, 1, 100, 9, 0.4], [12, 3, 500, 4, 300, 8, 0.6]], eps: 0, Ed: 0 }, note: 'd_p = 3 and 12, first-order decomposition A → R. For first order the feed concentration cancels: k′ = v X/(W(1 − X)).' },
      { name: 'Problem 18.24 — mixed flow runs', v: { mode: 'runs', shape: 'sph', n: 1, dtype: 'mfr', runsM: [[4, 3, 500, 1, 300, 60, 0.8], [8, 3, 500, 3, 100, 160, 0.6]], eps: 0, Ed: 0 }, note: 'As printed, the smaller pellet is three times as active for half the size — more than pore diffusion alone can explain.' },
      { name: 'Problem 18.25 — mixed flow runs', v: { mode: 'runs', shape: 'sph', n: 1, dtype: 'mfr', runsM: [[2, 3, 500, 4, 75, 10, 0.2], [1, 3, 500, 6, 100, 5, 0.6]], eps: 0, Ed: 0 }, check: { pore: 'strong' }, note: 'Halving the size doubles k′: strong pore diffusion.' },
      { name: 'Problem 18.26 — activation energy from runs', v: { mode: 'runs', shape: 'sph', n: 1, dtype: 'rate', runsR: [[1, 3, 480, 20, 1], [2, 3, 480, 40, 1], [2, 3, 500, 40, 3]], Ed: 0 }, note: 'Sizes 1 and 2 at 480 K show strong pore resistance; the 500 K run gives the observed E, and the true E is about twice that.' },
      { name: 'Problem 18.27 — recycle reactor runs', v: { mode: 'runs', shape: 'sph', n: 1, dtype: 'rate', runsR: [[1, 2, 500, 1, 4], [1, 3, 500, 1, 4], [2, 3, 500, 1, 3], [2, 2, 500, 1, 3]], Ed: 0 }, note: 'Same leaving composition in every run; recycle “high” (2) and “higher still” (3). Rate unchanged by recycle; doubling the size cuts it from 4 to 3.' }
    ],
    eqs: [
      ['\\frac{\\text{observed rate}}{\\text{film-limited rate}} = \\frac{k^{\\prime\\prime\\prime}_{obs}V_p}{k_g S_{ex}} = \\frac{(-r^{\\prime\\prime\\prime}_{A,obs})\\,L}{k_g C_{Ag}}', 'Eq. 18.58, Ex. 18.1'],
      ['M_W = \\frac{(-r_A^{\\prime\\prime\\prime})_{obs}L^2}{\\mathscr{D}_e C_{Ag}}', 'Eq. 18.24'],
      ['\\Delta T_{film} = \\frac{L(-r^{\\prime\\prime\\prime}_{A,obs})(-\\Delta H_r)}{h},\\qquad \\Delta T_{particle} = \\frac{\\mathscr{D}_e(C_{As}-C_{A,center})(-\\Delta H_r)}{k_{eff}}', 'Eqs. 18.36, 18.38'],
      ['\\frac{k_g d_p}{\\mathscr{D}} = 2 + 0.6\\,Re^{1/2}Sc^{1/3}\\ \\text{(single particle)},\\qquad 2 + 1.8\\,Re^{1/2}Sc^{1/3}\\ \\text{(packed bed, }Re>80)', 'Froessling, Ranz'],
      ['\\frac{r_1}{r_2} = \\frac{\\mathscr{E}(M_{T1})}{\\mathscr{E}(M_{T1}R_2/R_1)}', 'Eq. 18.25', 'Solved for M_T1: places each size on the effectiveness curve.'],
      ['\\mathbf{E}_{obs} = \\mathbf{E}\\left(1+\\tfrac{m}{2}\\right) - \\mathbf{E}_{diff}\\tfrac{m}{2},\\quad m = \\frac{d\\ln\\mathscr{E}}{d\\ln M_T}\\ \\to\\ \\frac{\\mathbf{E}+\\mathbf{E}_{diff}}{2}', 'Eqs. 18.31–18.33']
    ],
    method: ['One rate: the film ratio, M_W, and the temperature rises are the book’s expressions; the Weisz modulus is inverted on the exact effectiveness curve of the chosen shape and order (C-1801) to give M_T, ℰ, and the true k‴ = r_obs/(ℰ C_Asⁿ). The surface concentration is corrected for the film before M_W is evaluated; the book uses C_Ag, which differs only when the film matters.', 'Runs: runs at the same size and temperature but different agitation reveal film resistance. For two sizes at the same temperature, the ratio of rate constants is set equal to ℰ(M_T1)/ℰ(M_T1·d₂/d₁) and solved for M_T1 on the exact curve. For two temperatures at the same size, the observed E comes from the Arrhenius slope and is corrected with the local slope m of the effectiveness curve.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const sh = v.shape, n = +v.n, tu = v._tu || 'h';
      const tr = traceC({ shape: sh, kind: 'nth', n });
      const g = tr.g;
      if (v.mode === 'one') {
        const L = CAT.Lof(sh, +v.dp / 1000), r = +v.robs, Cg = +v.Cg, De = +v.De, kg = +v.kg;
        const film = r * L / (kg * Cg);
        const Cs = Cg * (1 - film);
        kpi.push({ k: 'film', l: 'Observed rate / film-limited rate (Eq. 18.58)', v: film, s: film < 0.05 ? 'film resistance negligible' : film > 0.5 ? 'film resistance dominant' : 'film resistance significant', state: film < 0.05 ? 'ok' : film > 0.5 ? 'bad' : 'warn' });
        const MW = r * L * L / (De * Cg);
        kpi.push({ k: 'MW', l: 'Weisz modulus M_W (on C_Ag, as in the book)', v: MW, main: true, s: CAT.regimeW(MW), state: CAT.state(MW < 0.15 ? 0 : MW > 4 ? 9 : 1) });
        const MWs = r * L * L / (De * Cs);
        const sol = PEL.findAll(tr, p => Math.log(p.eta * p.MT * p.MT / MWs));
        let eta, MT0;
        if (sol.length) { eta = sol[0].eta; MT0 = sol[0].MT; }
        else { const last = tr.pts[tr.pts.length - 1]; if (MWs > last.eta * last.MT * last.MT) { MT0 = MWs * g; eta = 1 / (g * MT0); } else { MT0 = Math.sqrt(MWs); eta = 1; } }
        kpi.push({ k: 'MT', l: 'Thiele modulus M_T behind it', v: MT0 * g, s: CAT.regime(MT0 * g) });
        kpi.push({ k: 'eta', l: 'Effectiveness factor ℰ', v: eta });
        const ktrue = r / (eta * Math.pow(Cs, n));
        kpi.push({ k: 'ktrue', l: 'True rate constant k‴ (free of pores and film)', v: ktrue, u: n === 1 ? '1/' + tu : '(m³/mol)^{' + fmt(n - 1, 3) + '}/' + tu });
        const dH = +v.dH, dTp = De * Cs * Math.abs(dH) / +v.keff, dTf = L * r * Math.abs(dH) / +v.h;
        kpi.push({ k: 'dTp', l: 'Largest ΔT inside the pellet (Eq. 18.38)', v: dTp, u: 'K', s: dTp < 1 ? 'pellet practically isothermal' : 'gradient inside the pellet', state: dTp < 1 ? 'ok' : dTp > 10 ? 'bad' : 'warn' });
        kpi.push({ k: 'dTf', l: 'Largest ΔT across the gas film (Eq. 18.36)', v: dTf, u: 'K', s: dH < 0 ? 'pellet hotter than the gas' : 'pellet colder than the gas', state: dTf < 1 ? 'ok' : dTf > 10 ? 'bad' : 'warn' });
        const TsK = +v.Ts + 273.15, beta = (dH < 0 ? 1 : -1) * dTp / TsK, gam = +v.E * 1000 / (8.314462618 * TsK);
        kpi.push({ k: 'beta', l: 'Prater number β and Arrhenius number γ', v: beta, s: 'γ = ' + fmt(gam, 4) + '; βγ = ' + fmt(beta * gam, 3) + (Math.abs(beta * gam) < 0.05 ? ' — heat effects inside the pellet negligible' : '') });
        if (Math.abs(beta * gam) > 0.05) {
          const trN = traceC({ shape: sh, kind: 'nth', n, beta, gamma: gam });
          const etN = PEL.etaAt(trN, MT0, true);
          info.push('With the Prater rise inside the pellet the effectiveness at this M_T becomes ' + etN.map(e => fmt(e.eta, 4)).join(' / ') + ' (C-1801 has the full curves).');
        }
        const shares = [film, 1 - film];
        tables.push({ title: 'Verdict', cols: ['resistance', 'measure', 'value', 'conclusion'], rows: [['gas film (mass)', 'observed / film-limited rate', film, film < 0.05 ? 'negligible' : 'matters'], ['pores', 'M_W', MW, CAT.regimeW(MW)], ['pellet interior (heat)', 'ΔT_max, K', dTp, dTp < 1 ? 'isothermal pellet' : 'nonisothermal'], ['gas film (heat)', 'ΔT_max, K', dTf, dTf < 1 ? 'negligible' : 'pellet ' + (dH < 0 ? 'hotter' : 'colder') + ' than the gas']] });
        plots.push({ title: 'Where the concentration goes: gas → surface → pellet centre', x: { cats: ['lost across the film', 'left at the pellet surface'] }, y: { label: 'share of C_Ag', min: 0, max: 1.1 }, series: [{ name: 'share', x: [0, 1], y: shares, type: 'bar', color: 2, legend: false }], marks: shares.map((q, i) => ({ t: 'text', x: i, y: q + 0.035, s: N.fmt(q, 3), anchor: 'middle' })), cap: 'Of the available concentration the pellet then uses a fraction ℰ = ' + fmt(eta, 3) + ' of its interior.' });
        const P = tr.pts;
        plots.push({ title: 'This pellet on the M_W chart', x: { label: 'M_W', log: true, min: 0.01, max: 1e3 }, y: { label: 'ℰ', log: true, min: 0.01, max: 1.5 }, series: [{ name: PEL.NAME[sh] + ', order ' + fmt(n, 3), x: P.map(p => p.eta * p.MT * p.MT), y: P.map(p => p.eta), color: 1, w: 2.2 }], marks: [{ t: 'pt', x: MWs, y: eta, color: 8, label: 'measured' }, { t: 'v', x: 0.15, color: 'var(--muted)', label: '0.15' }, { t: 'v', x: 4, color: 'var(--muted)', label: '4' }] });
        if (film < 0.05 && MW > 4) info.push('As in Example 18.1: the film is not limiting, but the pores are; the pellet is nearly isothermal while the film may still carry a temperature difference.');
        return { kpi, warn, info, plots, tables };
      }
      /* ---------- runs ---------- */
      const raw = v.dtype === 'mfr' ? (v.runsM || []) : (v.runsR || []);
      const runs = raw.filter(r => r && r.length >= 5 && r.slice(0, v.dtype === 'mfr' ? 7 : 5).every(x => x !== null && x !== '' && isFinite(+x))).map(r => {
        const d = +r[0], ag = +r[1], T = +r[2];
        if (v.dtype === 'mfr') {
          const W = +r[3], CA0 = +r[4], vv = +r[5], X = +r[6], eps = +v.eps;
          const CA = CA0 * (1 - X) / (1 + eps * X), rate = vv * CA0 * X / W;
          return { d, ag, T, C: CA, rate, k: rate / Math.pow(CA, n) };
        }
        const C = +r[3], rate = +r[4];
        return { d, ag, T, C, rate, k: rate / Math.pow(C, n) };
      });
      if (runs.length < 2) throw new Error('Enter at least two runs.');
      tables.push({ title: 'Runs and rate constants', cols: ['run', 'd_p', 'agitation', 'T, K', 'C_A', '−r_A', 'k = −r_A/C_Aⁿ'], rows: runs.map((q, i) => [i + 1, q.d, q.ag, q.T, q.C, q.rate, q.k]) });
      /* film: same size and T, different agitation */
      let filmVerdict = 'not tested', filmMax = 0;
      const key = q => q.d + '|' + q.T;
      const groups = {};
      runs.forEach(q => { (groups[key(q)] = groups[key(q)] || []).push(q); });
      Object.values(groups).forEach(gr => {
        if (gr.length < 2) return;
        const ags = new Set(gr.map(q => q.ag));
        if (ags.size < 2) return;
        const ks = gr.map(q => q.k), dev = (Math.max(...ks) - Math.min(...ks)) / Math.max(...ks);
        filmMax = Math.max(filmMax, dev);
        filmVerdict = filmMax < 0.05 ? 'negligible' : 'present';
      });
      kpi.push({ k: 'film', l: 'Film resistance', v: filmVerdict, s: filmVerdict === 'not tested' ? 'no pair of runs differs only in agitation' : 'largest change of k with agitation ' + fmt(filmMax * 100, 3) + '%', state: filmVerdict === 'present' ? 'bad' : filmVerdict === 'negligible' ? 'ok' : 'warn', main: true });
      if (filmVerdict === 'negligible') info.push('Changing the stirring (or recycle) rate does not change the rate, so the film is not limiting; by Eq. 18.57 it then does not limit smaller pellets either.');
      /* pores: two sizes at the same T, best agitation */
      const byT = {};
      runs.forEach(q => { const cur = byT[q.T] && byT[q.T][q.d]; byT[q.T] = byT[q.T] || {}; if (!cur || q.ag > cur.ag) byT[q.T][q.d] = q; });
      let poreVerdict = 'not tested', poreRows = [];
      const Ts = Object.keys(byT).map(Number).sort((a, b) => a - b);
      const etaG = Mg => PEL.etaAt(tr, Mg / g)[0].eta;
      const placed = {};
      Ts.forEach(T => {
        const ds = Object.keys(byT[T]).map(Number).sort((a, b) => a - b);
        if (ds.length < 2) return;
        const a = byT[T][ds[0]], b = byT[T][ds[ds.length - 1]];
        const ratio = a.k / b.k, sr = b.d / a.d;
        let verdict, M1 = NaN;
        if (ratio > sr * 1.03) verdict = 'inconsistent';
        else if (ratio >= sr * 0.97) verdict = 'strong';
        else if (ratio <= 1.03) verdict = 'diffusion-free';
        else {
          verdict = 'intermediate';
          const f = lm => Math.log(etaG(Math.exp(lm)) / etaG(Math.exp(lm) * sr)) - Math.log(ratio);
          const lm = N.brent(f, Math.log(1e-3), Math.log(1e3), { tol: 1e-12 });
          M1 = Math.exp(lm);
        }
        if (verdict === 'strong') M1 = 8; /* anywhere right of 4 fits; drawn at 8 */
        poreRows.push([T, a.d, b.d, sr, ratio, verdict, M1, isFinite(M1) ? etaG(M1) : NaN, isFinite(M1) ? etaG(M1 * sr) : NaN]);
        if (isFinite(M1)) { placed[T] = { d: a.d, M: M1, verdict }; }
        poreVerdict = poreVerdict === 'not tested' ? verdict : poreVerdict === verdict ? verdict : 'mixed';
      });
      kpi.push({ k: 'pore', l: 'Pore diffusion', v: poreVerdict === 'strong' ? 'strong' : poreVerdict, s: poreRows.length ? 'k ratio ' + fmt(poreRows[0][4], 4) + ' for a size ratio ' + fmt(poreRows[0][3], 4) : 'no two sizes at the same temperature', state: poreVerdict === 'diffusion-free' ? 'ok' : poreVerdict === 'intermediate' ? 'warn' : 'bad' });
      if (poreRows.length) {
        tables.push({ title: 'Placing the pellets on the effectiveness curve (Eq. 18.25)', cols: ['T, K', 'small d', 'large d', 'size ratio', 'k ratio', 'regime', 'M_T (small)', 'ℰ small', 'ℰ large'], rows: poreRows });
        const r0 = poreRows[0];
        if (r0[5] === 'intermediate') {
          kpi.push({ k: 'MTs', l: 'M_T of the smaller pellets', v: r0[6], s: 'ℰ = ' + fmt(r0[7], 4) + '; larger: M_T = ' + fmt(r0[6] * r0[3], 4) + ', ℰ = ' + fmt(r0[8], 4) });
          const small = byT[r0[0]][r0[1]];
          kpi.push({ k: 'ktrue', l: 'True rate constant (k of the small pellets / ℰ)', v: small.k / r0[7] });
        }
        if (r0[5] === 'inconsistent') warn.push('The smaller pellets are more active than the size ratio allows (at most ' + fmt(r0[3], 3) + '×). Pore diffusion alone cannot do this: suspect film resistance on the larger pellets, different catalyst, or the data.');
        if (r0[5] === 'strong') info.push('The rate constant scales as 1/d_p: strong pore resistance (Eq. 18.25). The true k is hidden; only k/M_T is measured.');
      }
      /* temperature: same size, different T */
      const sizes = [...new Set(runs.map(q => q.d))];
      sizes.forEach(d => {
        const atD = Ts.map(T => byT[T][d]).filter(Boolean);
        if (atD.length < 2) return;
        const a = atD[0], b = atD[atD.length - 1];
        const Eobs = 8.314462618e-3 * Math.log(b.k / a.k) / (1 / a.T - 1 / b.T);
        let m = NaN, note = '';
        const pl = placed[a.T];
        if (pl) {
          const Ma = pl.M * (d / pl.d);
          const hh = 0.01;
          m = Math.log(etaG(Ma * (1 + hh)) / etaG(Ma / (1 + hh))) / (2 * Math.log(1 + hh));
          if (pl.verdict === 'strong') m = -1;
        }
        const Ed = +v.Ed;
        const Etrue = isFinite(m) ? (Eobs + Ed * m / 2) / (1 + m / 2) : NaN;
        kpi.push({ k: 'Eobs', l: 'Observed activation energy (size ' + fmt(d, 3) + ')', v: Eobs, u: 'kJ/mol', s: a.T + ' → ' + b.T + ' K' });
        if (isFinite(Etrue)) kpi.push({ k: 'Etrue', l: 'True activation energy', v: Etrue, u: 'kJ/mol', s: m === -1 ? 'strong pore resistance: E = 2E_obs − E_diff' : 'corrected with d ln ℰ/d ln M_T = ' + fmt(m, 3) });
        else info.push('To correct the observed E for pore resistance, include two sizes at one of these temperatures.');
      });
      /* chart: runs on the effectiveness curve */
      const P = tr.pts, marks = [];
      poreRows.forEach(r => { if (isFinite(r[6])) { marks.push({ t: 'pt', x: r[6], y: r[7], color: 2, label: 'd = ' + fmt(r[1], 3) }); marks.push({ t: 'pt', x: r[6] * r[3], y: r[8], color: 3, label: 'd = ' + fmt(r[2], 3) }); } });
      plots.push({ title: 'The runs placed on the effectiveness curve', x: { label: 'M_T', log: true, min: 0.05, max: 100 }, y: { label: 'ℰ', log: true, min: 0.01, max: 1.5 }, series: [{ name: PEL.NAME[sh] + ', order ' + fmt(n, 3), x: P.map(p => p.MT * g), y: P.map(p => p.eta), color: 1, w: 2.2 }], marks: marks.concat([{ t: 'v', x: 0.4, color: 'var(--muted)', label: '0.4' }, { t: 'v', x: 4, color: 'var(--muted)', label: '4' }]), cap: 'Two sizes at the same temperature fix both points: their k ratio equals the ratio of their ℰ values, a fixed distance apart along the M_T axis.' });
      const ds = runs.map(q => q.d), ks = runs.map(q => q.k);
      plots.push({ title: 'Rate constant against pellet size', x: { label: 'd_p', log: true }, y: { label: 'k', log: true }, series: [{ name: 'runs', x: ds, y: ks, type: 'pts', color: 2 }], cap: 'A slope of −1 on these axes is the signature of strong pore resistance; a slope of 0 means no pore effect.' });
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
