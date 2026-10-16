/* Chapter 26 — fluid–particle reactors: design. Mixtures of sizes, plug, mixed and arbitrary
   residence-time distributions of solids, elutriation, and the gas composition coupled to the
   solids conversion (mixed or plug gas, any feed ratio). */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, FP = ST.FP;
  const fmt = (x, s = 4) => N.fmt(x, s);

  /* single-particle conversion against time for a size, given τ components at the reference size */
  function particle(v, R, cScale) {
    const Rr = +v.Rref, shp = v.beh === 'shrink' ? 'sph' : v.shape;
    const sc = p => Math.pow(R / Rr, p) * cScale;
    if (v.mech === 'comb') {
      const taus = { film: +v.tf * sc(+v.pf), ash: +v.ta * sc(2), rxn: +v.tr * sc(1) };
      const C = FP.combo(taus, shp);
      return { tau: C.tau, Xt: C.Xt };
    }
    const p = v.mech === 'film' || v.mech === 'stokes' || v.mech === 'large' ? +v.pf : v.mech === 'ash' ? 2 : 1;
    const tau = +v.tau * sc(+v.pexp > 0 ? +v.pexp : p);
    return { tau, Xt: t => FP.Xof(v.mech, shp, t / tau) };
  }
  /* mean unconverted fraction for one size under a solids RTD */
  function unconv(v, P, tm) {
    const tau = P.tau;
    if (v.flow === 'plug') return 1 - P.Xt(tm.tp);
    if (v.flow === 'mixed' || v.flow === 'elut') {
      const tb = tm.tb;
      /* in w = t/t̄, cut where e^{−w} is negligible so the quadrature cannot miss a narrow E */
      return N.quad(w => (1 - P.Xt(w * tb)) * Math.exp(-w), 0, Math.min(tau / tb, 60), { rtol: 1e-10 });
    }
    if (v.flow === 'tis') {
      const n = Math.max(1, Math.round(+v.Nt)), tb = tm.tb, lg = N.lgamma(n);
      const E = t => (t <= 0 ? (n === 1 ? 1 / tb : 0) : Math.exp((n - 1) * Math.log(t) + n * Math.log(n / tb) - n * t / tb - lg));
      return N.quad(t => (1 - P.Xt(t)) * E(t), 0, Math.min(tau, tb * (1 + 15 / Math.sqrt(n) + 40 / n)), { rtol: 1e-9 });
    }
    if (v.flow === 'rect') {
      const a = +v.t1, b = +v.t2;
      if (!(b > a)) return 1 - P.Xt(a);
      return N.quad(t => 1 - P.Xt(t), a, Math.min(b, Math.max(a, tau)), { rtol: 1e-10 }) / (b - a);
    }
    if (v.flow === 'table') {
      const rows = tm.Etab;
      let s = 0;
      for (let i = 1; i < rows.length; i++) {
        const [t0, e0] = rows[i - 1], [t1, e1] = rows[i];
        if (t0 >= tau) break;
        const hi = Math.min(t1, tau);
        s += N.quad(t => (1 - P.Xt(t)) * (e0 + (e1 - e0) * (t - t0) / (t1 - t0)), t0, hi, { rtol: 1e-8 });
      }
      return s;
    }
    return NaN;
  }
  /* book series (Eqs. 26.8b, 26.10b, 26.11) for mixed flow, single mechanism */
  const SERIES = {
    film: [1 / 2, -1 / 6, 1 / 24],
    rxn: [1 / 4, -1 / 20, 1 / 120],
    ash: [1 / 5, -19 / 420, 41 / 4620, -0.00149]
  };
  const seriesU = (mech, s, nterm) => SERIES[mech].slice(0, nterm || SERIES[mech].length).reduce((t, c, i) => t + c * Math.pow(s, i + 1), 0);

  ST.tool({
    id: 'fpdesign', tag: 'N-2601', ch: 26, sec: 'Ch. 26', title: 'Fluid–particle reactors: size mixtures, solids flow and gas coupling',
    lede: 'The mean conversion of solids leaving a fluid–particle reactor — rotary kiln, moving grate, fluidized bed, moving bed — from the conversion–time behaviour of single particles (any controlling step, or all of them together, constant size or shrinking) and the residence-time distribution of the solids: plug, mixed, tanks in series, a rectangular window or any measured E curve; with a mixture of sizes, elutriation of fines, and a gas whose composition follows from its own feed rate in mixed or plug flow. It also runs backwards: the residence time, or the solids hold-up, needed for a target conversion.',
    keys: 'fluid particle reactor design solids conversion size distribution mixture plug flow solids mixed flow fluidized bed rotary kiln moving grate elutriation residence time distribution of solids exact integral series expansion gas composition stoichiometric feed mixed gas plug gas log mean weight of bed example 26.1 26.2 26.3 26.4 problem 26.1 26.2 26.3 26.4 26.5 26.6 26.7 26.8 26.9 26.11 26.12 26.13 26.14 figure 26.4 26.5',
    units: { t: 'min' },
    inputs: [
      { group: 'Single particles (at the reference size, in the feed gas C_A0)' },
      { k: 'beh', t: 'seg', label: 'Particles', opts: [['const', 'keep their size'], ['shrink', 'shrink (no ash)']], v: 'const' },
      { k: 'shape', t: 'seg', label: 'Shape', opts: [['slab', 'flat plate'], ['cyl', 'cylinder'], ['sph', 'sphere']], v: 'sph', if: v => v.beh === 'const' },
      { k: 'mech', t: 'select', label: 'Controlling step', opts: [['rxn', 'reaction'], ['ash', 'ash diffusion'], ['film', 'gas film'], ['comb', 'all three together (times add, Eq. 25.32)'], ['stokes', 'shrinking: film, Stokes regime'], ['large', 'shrinking: film, large particles'], ['srxn', 'shrinking: reaction']], v: 'rxn' },
      { k: 'Rref', label: 'Reference size', unit: 'μm', v: 100, min: 1e-6, max: 1e9, scale: 'log' },
      { k: 'tau', label: 'τ at the reference size', unit: '[t]', v: 10, min: 1e-12, max: 1e12, scale: 'log', if: v => v.mech !== 'comb' },
      { k: 'pexp', label: 'Size exponent of τ (0 = the mechanism’s own)', v: 0, min: 0, max: 3, step: 0.1, if: v => v.mech !== 'comb' },
      { k: 'tf', label: 'τ_film at the reference size', unit: '[t]', v: 2, min: 0, max: 1e12, if: v => v.mech === 'comb' },
      { k: 'ta', label: 'τ_ash at the reference size', unit: '[t]', v: 3, min: 0, max: 1e12, if: v => v.mech === 'comb' },
      { k: 'tr', label: 'τ_reaction at the reference size', unit: '[t]', v: 5, min: 0, max: 1e12, if: v => v.mech === 'comb' },
      { k: 'pf', label: 'Size exponent of film control (1.5–2)', v: 2, min: 1, max: 2.5, step: 0.1, if: v => v.mech === 'film' || v.mech === 'comb' || v.mech === 'stokes' || v.mech === 'large' },
      { k: 'sizes', t: 'table', label: 'Feed: particle size and mass fraction', cols: ['R, μm', 'fraction', 'κ, 1/[t] (elutriation)'], v: [[50, 0.3, 0], [100, 0.4, 0], [200, 0.3, 0]] },
      { group: 'Flow of solids' },
      { k: 'flow', t: 'select', label: 'Solids residence times', opts: [['plug', 'plug flow (kiln, grate, moving bed)'], ['mixed', 'mixed flow (fluidized bed)'], ['elut', 'mixed flow with elutriation of fines'], ['tis', 'N mixed stages in series'], ['rect', 'uniform between t₁ and t₂'], ['table', 'measured E curve']], v: 'plug' },
      { k: 'find', t: 'seg', label: 'Find', opts: [['X', 'mean conversion'], ['t', 'time (and hold-up) for a target']], v: 'X', if: v => ['plug', 'mixed', 'elut', 'tis'].indexOf(v.flow) >= 0 },
      { k: 'Xt', label: 'Target mean conversion X̄_B', v: 0.9, min: 0.001, max: 1, if: v => v.find === 't' && ['plug', 'mixed', 'elut', 'tis'].indexOf(v.flow) >= 0 },
      { k: 'tp', label: 'Time in the reactor t_p', unit: '[t]', v: 8, min: 0, max: 1e12, if: v => v.flow === 'plug' && v.find !== 't' },
      { k: 'hold', t: 'seg', label: 'Mean residence time from', opts: [['tb', 't̄ directly'], ['WF', 'hold-up W and feed rate F']], v: 'tb', if: v => ['mixed', 'elut', 'tis'].indexOf(v.flow) >= 0 && (v.find !== 't' || v.flow === 'elut') },
      { k: 'tb', label: 'Mean residence time t̄', unit: '[t]', v: 10, min: 1e-12, max: 1e12, if: v => ['mixed', 'tis'].indexOf(v.flow) >= 0 && v.hold === 'tb' && v.find !== 't' },
      { k: 'W', label: 'Solids in the bed W', unit: 'kg', v: 10, min: 1e-12, max: 1e12, if: v => (['mixed', 'tis'].indexOf(v.flow) >= 0 && v.hold === 'WF' && v.find !== 't') || (v.flow === 'elut' && v.find !== 't') },
      { k: 'F', label: 'Solids feed rate F', unit: 'kg/[t]', v: 1, min: 1e-12, max: 1e12, if: v => ['mixed', 'elut', 'tis'].indexOf(v.flow) >= 0 && (v.hold === 'WF' || v.flow === 'elut' || v.find === 't') },
      { k: 'Nt', label: 'Number of stages N', v: 2, min: 1, max: 50, step: 1, if: v => v.flow === 'tis' },
      { k: 't1', label: 't₁', unit: '[t]', v: 15, min: 0, max: 1e12, if: v => v.flow === 'rect' },
      { k: 't2', label: 't₂', unit: '[t]', v: 45, min: 0, max: 1e12, if: v => v.flow === 'rect' },
      { k: 'Etab', t: 'table', label: 'E curve of the solids', cols: ['t', 'E'], units: ['[t]', '1/[t]'], v: [[0, 0], [10, 0.05], [20, 0.04], [30, 0.01], [40, 0]], if: v => v.flow === 'table' },
      { group: 'Gas' },
      { k: 'gas', t: 'select', label: 'Gas composition', opts: [['uni', 'uniform at C_A0 (gas in large excess)'], ['mix', 'gas in mixed flow, set by its feed rate'], ['plug', 'gas in plug flow (solids mixed)']], v: 'uni' },
      { k: 'phi', label: 'Gas fed / gas needed for complete conversion of B', v: 1, min: 0.01, max: 100, if: v => v.gas !== 'uni' }
    ],
    presets: [
      { name: 'Example 26.1 — size mixture on a moving grate', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 100, tau: 10, pexp: 0, sizes: [[50, 0.3, 0], [100, 0.4, 0], [200, 0.3, 0]], flow: 'plug', find: 'X', tp: 8, gas: 'uni', _tu: 'min' }, book: { X: { v: 0.932, tol: 0.002 } }, note: '30% of 50 μm, 40% of 100 μm, 30% of 200 μm with τ = 5, 10, 20 min (reaction controls) for 8 min on a grate in a crossflow of gas. Book: 93.2% converted.' },
      { name: 'Example 26.2 — pyrrhotite, mixed flow, reaction bound', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 100, tau: 20, pexp: 0, sizes: [[100, 1, 0]], flow: 'mixed', find: 'X', hold: 'tb', tb: 60, gas: 'uni', _tu: 'min' }, book: { U: { v: 0.078, tol: 0.01 } }, note: 'τ ∝ R^{1.5} lies between reaction and ash control; with τ = 20 min and t̄ = 60 min, reaction control leaves 7.8% unconverted …' },
      { name: 'Example 26.2 — pyrrhotite, mixed flow, ash bound', v: { beh: 'const', shape: 'sph', mech: 'ash', Rref: 100, tau: 20, pexp: 0, sizes: [[100, 1, 0]], flow: 'mixed', find: 'X', hold: 'tb', tb: 60, gas: 'uni', _tu: 'min' }, book: { U: { v: 0.062, tol: 0.01 } }, note: '… and ash control 6.2%: about 7% of the sulfide remains.' },
      { name: 'Example 26.3 — size mixture in a fluidized bed', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 100, tau: 10, pexp: 0, sizes: [[50, 0.3, 0], [100, 0.4, 0], [200, 0.3, 0]], flow: 'mixed', find: 'X', hold: 'WF', W: 10, F: 1, gas: 'uni', _tu: 'min' }, book: { X: { v: 0.778, tol: 0.002 } }, note: 'The same feed in a fluidized bed holding 10 kg with 1 kg/min of feed: t̄ = 10 min. Book: 77.8%.' },
      { name: 'Example 26.4 — bed weight for 90%, stoichiometric gas', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 100, tau: 60, pexp: 0, sizes: [[100, 1, 0]], flow: 'mixed', find: 't', Xt: 0.9, F: 1000 / 60, gas: 'mix', phi: 1, _tu: 'min' }, book: { Wseries: { v: 23000, tol: 0.01 }, W: { v: 23000, tol: 0.01 } }, note: '1 ton/hr of B to 90% with exactly stoichiometric A in mixed flow: the solids see C_A = 0.1 C_A0, so τ = 10 hr. Book, from the first two terms of Eq. 26.10b: τ/t̄ = 0.435, t̄ = 23 hr, W = 23 tons; the exact integral agrees.' },
      { name: 'Problem 26.8 — twice the stoichiometric gas', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 100, tau: 60, pexp: 0, sizes: [[100, 1, 0]], flow: 'mixed', find: 't', Xt: 0.9, F: 1000 / 60, gas: 'mix', phi: 2, _tu: 'min' } },
      { name: 'Problem 26.9 — gas in plug flow', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 100, tau: 60, pexp: 0, sizes: [[100, 1, 0]], flow: 'mixed', find: 't', Xt: 0.9, F: 1000 / 60, gas: 'plug', phi: 1, _tu: 'min' } },
      { name: 'Problems 26.1–26.2 — reactor twice as big (ash, 80%)', v: { beh: 'const', shape: 'sph', mech: 'ash', Rref: 100, tau: 1, pexp: 0, sizes: [[100, 1, 0]], flow: 'plug', find: 'X', tp: 2 * FP.g('ash', 'sph', 0.8), gas: 'uni', _tu: 'min' }, note: 'With τ = 1, 80% conversion in plug flow needs t = 0.374; twice the reactor gives 96.5%. Switch the solids to mixed flow for Problem 26.2.' },
      { name: 'Problem 26.3 — kiln, time for complete conversion', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 2000, tau: 240, pexp: 0, sizes: [[500, 0.2, 0], [1000, 0.3, 0], [2000, 0.5, 0]], flow: 'plug', find: 't', Xt: 1, gas: 'uni', _tu: 'min' }, check: { T: { v: 240, tol: 1e-9 } }, note: '20% of 1 mm and smaller, 30% of 2 mm, 50% of 4 mm particles (diameters; radii entered), τ = 4 h for 4 mm: complete conversion takes 4 h.' },
      { name: 'Problem 26.4 — kiln, 15 minutes', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 2000, tau: 240, pexp: 0, sizes: [[500, 0.2, 0], [1000, 0.3, 0], [2000, 0.5, 0]], flow: 'plug', find: 'X', tp: 15, gas: 'uni', _tu: 'min' }, check: { X: { v: 0.2 * (1 - Math.pow(0.75, 3)) + 0.3 * (1 - Math.pow(0.875, 3)) + 0.5 * (1 - Math.pow(0.9375, 3)), tol: 1e-6 } } },
      { name: 'Problem 26.6 — scale-up of a fluidized roaster', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 300, tau: 126.2 / 60, pexp: 0, sizes: [[300, 1, 0]], flow: 'mixed', find: 't', Xt: 0.98, F: 4000 / 60, gas: 'uni', _tu: 'min' }, note: 'The bench bed (10 g/s through 1000 g, so t̄ = 100 s, and 75% converted) fixes τ = 126 s for these 0.3 mm particles under reaction control — check it with τ = 2.103 min, t̄ = 1.667 min. The commercial bed for 4 t/hr to 98% then needs t̄ ≈ 26 min and about 1.7 t of solids.' },
      { name: 'Problems 26.11–26.13 — moving bed, uniform E', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 100, tau: 60, pexp: 0, sizes: [[100, 1, 0]], flow: 'rect', t1: 15, t2: 45, gas: 'uni', _tu: 'min' }, note: 'Fe₂O₃ → FeS with τ = 1 h; E uniform from 15 to 45 min (P26.12). Try 0–90 and 30–90 for P26.11 and P26.13.' },
      { name: 'Problem 26.14 — plug flow, 45 min', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 100, tau: 60, pexp: 0, sizes: [[100, 1, 0]], flow: 'plug', find: 'X', tp: 45, gas: 'uni', _tu: 'min' }, check: { X: { v: 1 - Math.pow(0.25, 3), tol: 1e-9 } } },
      { name: 'Elutriation of fines from a fluidized bed', v: { beh: 'const', shape: 'sph', mech: 'rxn', Rref: 100, tau: 10, pexp: 0, sizes: [[25, 0.2, 0.5], [50, 0.3, 0.08], [100, 0.3, 0], [200, 0.2, 0]], flow: 'elut', find: 'X', W: 10, F: 1, gas: 'uni', _tu: 'min' }, note: 'Fines blown out of the bed (κ, the fraction of their bed inventory elutriated per minute) leave sooner and less converted: here the smallest size spends only about 2 minutes in the bed.' }
    ],
    eqs: [
      ['1-\\bar X_B = \\sum_{R(t_p=\\tau)}^{R_m}[1-X_B(R_i)]\\frac{F(R_i)}{F}', 'Eq. 26.2', 'Plug flow of solids, size mixture.'],
      ['1-\\bar X_B = \\int_0^{\\tau}(1-X_B)\\frac{e^{-t/\\bar t}}{\\bar t}\\,dt', 'Eq. 26.6', 'Mixed flow, one size.'],
      ['\\bar X_B = 3\\frac{\\bar t}{\\tau}-6\\left(\\frac{\\bar t}{\\tau}\\right)^2+6\\left(\\frac{\\bar t}{\\tau}\\right)^3\\left(1-e^{-\\tau/\\bar t}\\right)', 'Eq. 26.10a', 'Reaction control, mixed flow.'],
      ['1-\\bar X_B = \\frac{1}{5}\\frac{\\tau}{\\bar t}-\\frac{19}{420}\\left(\\frac{\\tau}{\\bar t}\\right)^2+\\frac{41}{4620}\\left(\\frac{\\tau}{\\bar t}\\right)^3-0.00149\\left(\\frac{\\tau}{\\bar t}\\right)^4+\\cdots', 'Eq. 26.11', 'Ash control, mixed flow.'],
      ['\\bar t = \\frac{W}{F},\\qquad 1-\\bar X_B = \\sum[1-\\bar X_B(R_i)]\\frac{F(R_i)}{F}', 'Eqs. 26.13–26.15'],
      ['1-\\bar X_B = \\int_0^{\\tau}(1-X_B)_{single}\\,\\mathbf{E}_{solid}\\,dt', 'Eq. 26.19', 'Any solids RTD.'],
      ['\\tau \\propto \\frac{1}{C_A}', 'Ex. 26.4', 'Every step of the shrinking-core model is first order in the gas.']
    ],
    method: ['Every size is treated as a macrofluid: its unconverted fraction is the integral of the single-particle 1 − X_B(t) over the residence-time distribution of the solids, done by adaptive Gauss–Kronrod quadrature up to τ, so no series is truncated (the book’s series values are shown beside the exact ones). Sizes are weighted by their mass fraction in the stream leaving, which equals the feed for a single exit stream.', 'Elutriation: each size leaves by overflow at rate W(R)·F₁/W and by entrainment at κ(R)W(R); the overflow rate is found by Brent’s method so that the bed holds W, which gives each size its own mean residence time 1/(F₁/W + κ). Gas coupling: τ scales as C_A0/C_A. With mixed gas the solids see the exit gas, C_A = C_A0(1 − X̄_B/φ); with plug gas through well-mixed solids they see the log-mean of inlet and outlet. The conversion and the gas composition are solved together by Brent’s method.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'min';
      const rows = (v.sizes || []).filter(r => r && +r[0] > 0 && +r[1] > 0).map(r => ({ R: +r[0], f: +r[1], kap: +r[2] || 0 }));
      if (!rows.length) throw new Error('Give at least one particle size with a positive fraction.');
      const ft = N.sum(rows.map(r => r.f));
      rows.forEach(r => { r.f /= ft; });
      if (v.gas === 'plug' && v.flow !== 'mixed' && v.flow !== 'elut') throw new Error('Plug flow of gas is handled for well-mixed solids, where every particle sees the bed-average gas.');
      const Etab = (v.Etab || []).filter(r => r && isFinite(+r[0]) && isFinite(+r[1]) && r[0] !== '' && r[1] !== '').map(r => [+r[0], Math.max(+r[1], 0)]).sort((a, b) => a[0] - b[0]);
      if (v.flow === 'table') {
        if (Etab.length < 2) throw new Error('Give at least two points of the E curve.');
        const area = N.trapz(Etab.map(r => r[0]), Etab.map(r => r[1]));
        if (!(area > 0)) throw new Error('The E curve has no area.');
        Etab.forEach(r => { r[1] /= area; });
      }
      /* residence time of each size, for a given overall time parameter T (t_p or t̄) */
      const tmFor = T => {
        if (v.flow === 'elut') {
          /* T = W/F; solve for λ = F1/W */
          const lam = rows.every(r => r.kap === 0) ? 1 / T : (() => {
            const g = l => N.sum(rows.map(r => r.f / (l + r.kap))) - T;
            const lo = 1e-12 / T, hi = 1e6 / T;
            if (!(g(lo) > 0)) return NaN;
            return N.brent(g, lo, hi, { tol: 1e-14 / T });
          })();
          return rows.map(r => ({ tb: 1 / (lam + r.kap), lam }));
        }
        return rows.map(() => ({ tp: T, tb: T, Etab }));
      };
      const evalAt = (T, cScale) => {
        const tms = tmFor(T);
        if (['mixed', 'elut', 'tis'].indexOf(v.flow) >= 0 && tms.some(q => !isFinite(q.tb))) return { U: NaN };
        const per = rows.map((r, i) => { const P = particle(v, r.R, cScale); const u = Math.max(0, unconv(v, P, tms[i])); return { R: r.R, f: r.f, tau: P.tau, U: u, tb: tms[i].tb }; });
        return { U: N.sum(per.map(q => q.f * q.U)), per, tms };
      };
      /* gas coupling: returns concentration scale C_A0/C_A for a given X̄ */
      const phi = +v.phi;
      const cOf = X => {
        if (v.gas === 'uni') return 1;
        const XA = X / phi;
        if (XA >= 1) return Infinity;
        if (v.gas === 'mix') return 1 / (1 - XA);
        const r = 1 - XA;
        return r >= 1 - 1e-12 ? 1 : Math.log(1 / r) / (1 - r);
      };
      const solveX = T => {
        if (v.gas === 'uni') { const e = evalAt(T, 1); return Object.assign(e, { X: 1 - e.U, c: 1 }); }
        const top = Math.min(1, phi) * (1 - 1e-12);
        const g = X => { const e = evalAt(T, cOf(X)); return (1 - e.U) - X; };
        const X = N.brent(g, 0, top, { tol: 1e-12 });
        const e = evalAt(T, cOf(X));
        return Object.assign(e, { X, c: cOf(X) });
      };
      const Tlabel = v.flow === 'plug' ? 't_p' : 't̄';
      let T, S;
      if (v.find === 't' && ['plug', 'mixed', 'elut', 'tis'].indexOf(v.flow) >= 0) {
        const Xt = +v.Xt;
        if (v.gas !== 'uni' && Xt >= phi) throw new Error('With φ = ' + fmt(phi, 3) + ' the gas cannot convert more than ' + fmt(phi * 100, 3) + '% of B.');
        const cs = cOf(Xt);
        const taumax = Math.max(...rows.map(r => particle(v, r.R, cs).tau));
        if (v.flow !== 'plug' && Xt >= 1 - 1e-9) throw new Error('Mixed flow always lets some solids out early: complete conversion needs plug flow.');
        if (v.flow === 'plug') T = Xt >= 1 - 1e-9 ? taumax : N.brent(q => (1 - evalAt(q, cs).U) - Xt, 0, taumax, { tol: 1e-12 * taumax });
        else {
          let hi = taumax;
          while (1 - evalAt(hi, cs).U < Xt && hi < 1e9 * taumax) hi *= 2;
          T = N.brent(q => (1 - evalAt(q, cs).U) - Xt, taumax * 1e-9, hi, { tol: 1e-12 * hi });
        }
        S = Object.assign(evalAt(T, cs), { X: Xt, c: cs });
        kpi.push({ k: 'T', l: (v.flow === 'plug' ? 'Time in the reactor' : 'Mean residence time') + ' for X̄_B = ' + fmt(Xt, 4), v: T, u: tu, main: true });
        if (v.flow !== 'plug' && +v.F > 0) kpi.push({ k: 'W', l: 'Solids hold-up W = t̄F', v: T * +v.F, u: 'kg', s: 'feed ' + fmt(+v.F, 4) + ' kg/' + tu });
        /* book series check for one size, single mechanism, mixed */
        if (v.flow === 'mixed' && rows.length === 1 && SERIES[v.mech]) {
          const tau = particle(v, rows[0].R, cs).tau;
          const nt = 2;
          const sOf = N.brent(s => seriesU(v.mech, s, nt) - (1 - Xt), 1e-9, 2, { tol: 1e-14 });
          if (isFinite(sOf)) {
            const Ts = tau / sOf;
            kpi.push({ k: 'Tseries', l: 'Same from the first two terms of the book’s series', v: Ts, u: tu, s: 'τ/t̄ = ' + fmt(sOf, 4) });
            if (+v.F > 0) kpi.push({ k: 'Wseries', l: 'Hold-up from the two-term series', v: Ts * +v.F, u: 'kg' });
          }
        }
      } else {
        T = v.flow === 'plug' ? +v.tp : (v.flow === 'mixed' || v.flow === 'tis') ? (v.hold === 'WF' ? +v.W / +v.F : +v.tb) : v.flow === 'elut' ? +v.W / +v.F : NaN;
        S = solveX(T);
        if (!isFinite(S.X)) throw new Error('No consistent solution: check the hold-up against the elutriation rates.');
        kpi.push({ k: 'X', l: 'Mean conversion of solids X̄_B', v: S.X, main: true });
        kpi.push({ k: 'U', l: 'Unconverted fraction 1 − X̄_B', v: 1 - S.X });
        if (isFinite(T)) kpi.push({ k: 'T', l: v.flow === 'plug' ? 'Time in the reactor t_p' : v.flow === 'elut' ? 'W/F' : 'Mean residence time t̄', v: T, u: tu });
        if (v.flow === 'mixed' && rows.length === 1 && SERIES[v.mech] && v.gas === 'uni') {
          const s = S.per[0].tau / T;
          kpi.push({ k: 'Useries', l: '1 − X̄_B from the book’s series (Eqs. 26.8b, 26.10b, 26.11)', v: seriesU(v.mech, s), s: s > 1 ? 'τ/t̄ > 1: the series is no longer reliable' : '' });
        }
      }
      if (v.gas !== 'uni') {
        kpi.push({ k: 'CA', l: 'C_A the solids see / C_A0', v: 1 / S.c, s: v.gas === 'plug' ? 'log mean of inlet and outlet' : 'the exit gas' });
        kpi.push({ k: 'XA', l: 'Conversion of the gas X_A', v: S.X / phi });
      }
      if (v.flow === 'elut') {
        const lam = S.tms[0].lam;
        kpi.push({ k: 'over', l: 'Solids leaving by overflow / feed', v: lam * T, s: 'the rest, ' + fmt((1 - lam * T) * 100, 3) + '%, is elutriated' });
      }
      /* per-size table */
      tables.push({ title: 'Each size', cols: ['R, μm', 'fraction', 'τ(R), ' + tu, v.flow === 'plug' ? 'X_B' : 't̄(R), ' + tu, 'unconverted 1 − X̄', 'contribution to 1 − X̄'], rows: S.per.map(q => [q.R, q.f, q.tau, v.flow === 'plug' ? 1 - q.U : q.tb, q.U, q.f * q.U]) });
      if (v.flow === 'plug') { const done = S.per.filter(q => q.U <= 1e-12).map(q => q.R); if (done.length) info.push('Sizes up to ' + fmt(Math.max(...done), 4) + ' μm are converted completely and drop out of Eq. 26.2.'); }
      /* plots: single-particle curves with the RTD */
      const tmax = Math.max(...S.per.map(q => q.tau)) * 1.05;
      const tg = N.linspace(0, tmax, 200);
      const ser = rows.map((r, i) => { const P = particle(v, r.R, S.c || 1); return { name: '1 − X_B, R = ' + fmt(r.R, 4) + ' μm', x: tg, y: tg.map(t => 1 - P.Xt(t)), color: (i % 7) + 1, w: 1.8 }; });
      const marks = [];
      if (v.flow === 'plug') marks.push({ t: 'v', x: T, color: 8, label: 't_p' });
      else {
        let Ef = null;
        if (v.flow === 'mixed') Ef = t => Math.exp(-t / T) / T;
        else if (v.flow === 'tis') { const n = Math.max(1, Math.round(+v.Nt)), lg = N.lgamma(n); Ef = t => (t <= 0 ? (n === 1 ? 1 / T : 0) : Math.exp((n - 1) * Math.log(t) + n * Math.log(n / T) - n * t / T - lg)); }
        else if (v.flow === 'rect') Ef = t => (t >= +v.t1 && t <= +v.t2 ? 1 / (+v.t2 - +v.t1) : 0);
        else if (v.flow === 'table') Ef = t => N.interpClamp(Etab.map(r => r[0]), Etab.map(r => r[1]), t);
        if (Ef) { const Ev = tg.map(Ef), mx = Math.max(...Ev) || 1; ser.push({ name: 'E of the solids (scaled)', x: tg, y: Ev.map(q => q / mx), color: 'var(--muted)', type: 'area', fill: 0.18, line: true }); }
      }
      plots.push({ title: 'Single particles and the time the solids spend in the reactor', x: { label: 't', unit: tu, min: 0, max: tmax }, y: { label: '1 − X_B', min: 0, max: 1 }, series: ser, marks, cap: 'Each size reacts along its own curve; the mean conversion weights these curves by the residence-time distribution of the solids (Eq. 26.19).' });
      /* Figs 26.4–26.5 for a single size */
      const us = N.logspace(0.01, 0.97, 32), shp = v.beh === 'shrink' ? 'sph' : v.shape;
      const gp = (m, X) => {
        const y = Math.max(1 - X, 1e-300);
        switch (m) {
          case 'film': return 1;
          case 'rxn': return shp === 'slab' ? 1 : shp === 'cyl' ? 0.5 / Math.sqrt(y) : Math.pow(y, -2 / 3) / 3;
          case 'ash': return shp === 'slab' ? 2 * X : shp === 'cyl' ? -Math.log(y) : 2 * (Math.pow(y, -1 / 3) - 1);
          case 'stokes': return 2 / 3 * Math.pow(y, -1 / 3);
          case 'large': return 0.5 / Math.sqrt(y);
          case 'srxn': return Math.pow(y, -2 / 3) / 3;
        }
        return NaN;
      };
      /* 1 − X̄ for one size in mixed flow at s = τ/t̄, integrating over X instead of t (no inversion) */
      const Umix = (m, s) => { const Xc = s * 1 > 60 ? FP.Xof(m, shp, 60 / s) : 1; return N.quad(X => (1 - X) * s * Math.exp(-s * FP.g(m, shp, X)) * gp(m, X), 0, Xc, { rtol: 1e-8 }); };
      const mixedS = (m, U) => 1 / N.brent(s => Umix(m, s) - U, 1e-6, 1e4, { tol: 1e-10 });
      const mechs = v.beh === 'shrink' ? [['stokes', 2], ['large', 4], ['srxn', 1]] : [['film', 2], ['rxn', 1], ['ash', 3]];
      const curves = mechs.map(([m, c]) => ({ m, c, y: us.map(u => mixedS(m, u)) }));
      plots.push({ title: 'Mixed flow of one size: t̄/τ for a given conversion (Fig. 26.4)', x: { label: '1 − X̄_B', log: true, min: 0.01, max: 1 }, y: { label: 't̄/τ', log: true, min: 0.1, max: 100 }, series: curves.map(q => ({ name: FP.MECH[q.m], x: us, y: q.y, color: q.c, w: 1.8 })) });
      plots.push({ title: 'Mixed against plug flow of solids: hold-up ratio (Fig. 26.5)', x: { label: '1 − X̄_B', log: true, min: 0.01, max: 1 }, y: { label: 't̄_mixed / t_plug', log: true, min: 1, max: 50 }, series: curves.map(q => ({ name: FP.MECH[q.m], x: us, y: us.map((u, i) => q.y[i] / FP.g(q.m, shp, 1 - u)), color: q.c, w: 1.8 })), cap: 'At high conversion mixed flow of solids needs a far larger hold-up than plug flow.' });
      if (rows.length > 1) plots.push({ title: 'Feed and unconverted solids by size', x: { cats: rows.map(r => fmt(r.R, 4) + ' μm') }, y: { label: 'fraction of the feed', min: 0, max: Math.max(...rows.map(r => r.f)) * 1.2 }, series: [{ name: 'feed', x: rows.map((r, i) => i - 0.18), y: rows.map(r => r.f), type: 'bar', color: 3, bw: 26 }, { name: 'leaving unconverted', x: rows.map((r, i) => i + 0.18), y: S.per.map(q => q.f * q.U), type: 'bar', color: 1, bw: 26 }] });
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
