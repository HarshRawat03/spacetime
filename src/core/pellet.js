/* Catalyst pellets: effectiveness factors for slab, cylinder and sphere with any rate form,
   isothermal or with Prater's temperature coupling (Weisz–Hicks), dead cores for orders below one,
   inert cores (egg-shell catalysts), molar change in the pores (Thiele), external films,
   parallel reactions inside the pellet and the linear series scheme A → R → S (Wheeler). */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N;
  const PEL = ST.PEL = {};
  PEL.S = { slab: 0, cyl: 1, sph: 2 };
  const sOf = sh => (typeof sh === 'number' ? sh : PEL.S[sh] !== undefined ? PEL.S[sh] : 2);
  PEL.sOf = sOf;
  PEL.NAME = { slab: 'flat plate', cyl: 'cylinder', sph: 'sphere' };

  /* exponentially scaled modified Bessel functions e^{-x} I_ν(x), ν = 0, 1 */
  N.besselIe = function (nu, x) {
    x = Math.abs(x);
    if (x <= 17) {
      const q = x * x / 4;
      let term = nu ? x / 2 : 1, sum = term;
      for (let k = 1; k < 300; k++) { term *= q / (k * (k + nu)); sum += term; if (term < 1e-17 * sum) break; }
      return sum * Math.exp(-x);
    }
    const mu = 4 * nu * nu;
    let term = 1, sum = 1, prev = Infinity;
    for (let k = 1; k < 60; k++) {
      const t = term * -(mu - (2 * k - 1) * (2 * k - 1)) / (k * 8 * x);
      if (Math.abs(t) > prev) break;
      prev = Math.abs(t); term = t; sum += term;
      if (Math.abs(term) < 1e-17) break;
    }
    return sum / Math.sqrt(2 * Math.PI * x);
  };

  /* first-order effectiveness factor, exact, M_T = L sqrt(k/D) with L = V/S (Eqs. 18.20–18.22) */
  PEL.eta1 = function (shape, M) {
    const s = sOf(shape);
    if (!(M > 0)) return 1;
    if (s === 0) return M < 1e-4 ? 1 - M * M / 3 : Math.tanh(M) / M;
    if (s === 1) return M < 1e-4 ? 1 - M * M / 2 : N.besselIe(1, 2 * M) / (M * N.besselIe(0, 2 * M));
    if (M < 1e-4) return 1 - 0.6 * M * M;
    const p = 3 * M;
    return (1 / Math.tanh(p) - 1 / p) / M;
  };

  /* dimensionless rate f(c), normalized so that f(1) = 1 at the pellet surface.
     o: {kind:'nth'|'rev'|'lh'|'lh2', n, ce, K, beta, gamma} */
  PEL.rate = function (o) {
    const n = +o.n || 0, ce = +o.ce || 0, K = +o.K || 0, b = +o.beta || 0, g = +o.gamma || 0;
    let base;
    if (o.kind === 'rev') base = c => (c - ce) / (1 - ce);
    else if (o.kind === 'lh') base = c => (c <= 0 ? 0 : c * (1 + K) / (1 + K * c));
    else if (o.kind === 'lh2') base = c => (c <= 0 ? 0 : c * Math.pow((1 + K) / (1 + K * c), 2));
    else base = n === 0 ? (c => (c > 0 ? 1 : 0)) : (c => (c <= 0 ? 0 : Math.pow(c, n)));
    if (!b || !g) return base;
    return c => { const u = b * (1 - c); return base(c) * Math.exp(g * u / (1 + u)); };
  };
  /* lower limit of c (equilibrium for reversible) */
  PEL.cmin = o => (o.kind === 'rev' ? +o.ce || 0 : 0);
  /* generalized-modulus factor (Eq. 18.26): M_gen = M_T / sqrt(2 ∫_{cmin}^{1} f dc) */
  PEL.genFactor = (f, cmin) => 1 / Math.sqrt(2 * N.quad(f, cmin, 1, { rtol: 1e-10 }));

  /* integrate outward from the centre (or from the edge of a dead core, or from the edge of an
     inert core) until c = 1 at the surface. Returns {phi = R sqrt(k/D) in ζ units, eta, P[]}.
     o: {a: molar-change factor (ν−1)y_s, zc/deadN: dead core, zIn: inert core radius, extra: [w_j(c)], full} */
  function shoot(f, s, c0, o = {}) {
    const a = o.a || 0, ex = o.extra || [], ne = ex.length;
    let z0, y0;
    if (o.zc !== undefined) {
      const n = o.deadN, p = 2 / (1 - n), A = Math.pow((1 - n) * (1 - n) / (2 * (1 + n)), 1 / (1 - n));
      const del = 1e-4 * Math.max(1, o.zc) + 1e-6;
      z0 = o.zc + del; y0 = [A * Math.pow(del, p), A * p * Math.pow(del, p - 1)];
    } else if (o.zIn > 0) { z0 = o.zIn; y0 = [c0, 0]; }
    else {
      const fc = f(c0) * (1 + a * c0);
      z0 = 1e-7; y0 = [c0 + fc * z0 * z0 / (2 * (s + 1)), fc * z0 / (s + 1)];
    }
    for (let j = 0; j < ne; j++) y0.push(0);
    const rhs = (z, y) => {
      const c = Math.max(y[0], 0), fc = f(c);
      const d = [y[1], (a ? fc * (1 + a * c) + a * y[1] * y[1] / (1 + a * c) : fc) - (s ? s * y[1] / z : 0)];
      for (let j = 0; j < ne; j++) d.push((s ? Math.pow(z, s) : 1) * ex[j](c));
      return d;
    };
    const atol = Math.max(y0[0], 1e-300) * 1e-9;
    const opt = { rtol: 1e-10, atol, event: (z, y) => y[0] - 1, minPts: o.full ? 120 : 16, stop: (z, y) => !isFinite(y[0]) };
    let span = 4, za = z0, ya = y0, sol = null;
    const T = [], Y = [];
    for (let k = 0; k < 30; k++) {
      sol = N.ode(rhs, za, ya, za + span, opt);
      if (o.full) { const st = T.length ? 1 : 0; for (let i = st; i < sol.t.length; i++) { T.push(sol.t[i]); Y.push(sol.y[i]); } }
      if (sol.status === 'event') break;
      if (sol.status !== 'ok') return null;
      za = sol.t[sol.t.length - 1]; ya = sol.last();
      if (ya[1] <= 0 && f(Math.max(ya[0], 0)) <= 0) return null;
      span *= 2;
      if (za > 2e5) return null;
    }
    if (!sol || sol.status !== 'event') return null;
    const Z = sol.t[sol.t.length - 1], y = sol.last();
    const r = { phi: Z, eta: (s + 1) * y[1] / ((1 + a) * Z), flux: y[1] / (1 + a), P: ex.map((_, j) => y[2 + j]) };
    if (o.full) { r.z = T; r.c = Y.map(q => q[0]); r.z0 = z0; }
    return r;
  }
  PEL.shoot = shoot;

  /* adaptive tracing in a parameter q: refine wherever neighbours are far apart on log axes */
  function adaptive(at, qhi, qlo, n0, maxPts, dl) {
    let P = N.linspace(qhi, qlo, n0).map(at);
    for (let pass = 0; pass < 12 && P.length < maxPts; pass++) {
      const nxt = [P[0]];
      let added = 0;
      for (let i = 1; i < P.length; i++) {
        const a = P[i - 1], b = P[i];
        if (!a.bad && !b.bad && P.length + added < maxPts) {
          const gap = Math.max(Math.abs(Math.log10(b.phi / a.phi)), Math.abs(Math.log10(b.eta / a.eta)));
          if (gap > dl && Math.abs(a.q - b.q) > 1e-6) { nxt.push(at(0.5 * (a.q + b.q))); added++; }
        }
        nxt.push(b);
      }
      P = nxt;
      if (!added) break;
    }
    return P.filter(p => !p.bad);
  }

  /* trace the η–M_T curve parametrically in the centre concentration c0; works through S-shaped,
     multiple-solution curves. o adds {shape, a (molar change)} to the kinetics. */
  PEL.trace = function (o, opt = {}) {
    const s = sOf(o.shape), a = +o.a || 0;
    const f = PEL.rate(o), cmin = PEL.cmin(o), span = 1 - cmin;
    const toC = q => cmin + span / (1 + Math.exp(-q));
    const at = q => { const c0 = toC(q), r = shoot(f, s, c0, { a }); return r && isFinite(r.eta) && r.eta > 0 ? { q, c0, phi: r.phi, eta: r.eta } : { q, c0, bad: true }; };
    const pts = adaptive(at, 23, opt.qmin || -69, opt.n0 || 47, opt.maxPts || 360, opt.dl || 0.035).map(p => ({ phi: p.phi, eta: p.eta, c0: p.c0, q: p.q }));
    const n = +o.n;
    if ((o.kind === 'nth' || !o.kind) && n < 1 && pts.length && !(+o.beta && +o.gamma) && !a) {
      const phiC = pts[pts.length - 1].phi;
      N.logspace(1e-4 * phiC, 3e3 * phiC, 80).forEach(zc => { const r = shoot(f, s, 0, { zc, deadN: n }); if (r && isFinite(r.eta)) pts.push({ phi: r.phi, eta: r.eta, c0: 0, core: zc / r.phi, zc }); });
    }
    const g = a ? 1 : PEL.genFactor(f, cmin);
    pts.forEach(p => { p.MT = p.phi / (s + 1); p.Mgen = p.MT * g; });
    return { pts, s, f, g, cmin, toC, o, a };
  };
  /* all effectiveness factors at a given M_T (L sqrt(k/D) on the surface rate); refine = polish by Brent */
  PEL.etaAt = function (tr, MT, refine) {
    const out = [], P = tr.pts;
    if (!P.length) return [{ eta: NaN }];
    for (let i = 1; i < P.length; i++) {
      const a = P[i - 1], b = P[i];
      if ((a.MT - MT) * (b.MT - MT) <= 0 && a.MT !== b.MT) {
        const t = (Math.log(MT) - Math.log(a.MT)) / (Math.log(b.MT) - Math.log(a.MT));
        let e = { eta: Math.exp(Math.log(a.eta) + t * (Math.log(b.eta) - Math.log(a.eta))), c0: a.c0 + t * (b.c0 - a.c0), core: a.core !== undefined && b.core !== undefined ? a.core + t * (b.core - a.core) : undefined };
        if (e.core !== undefined) e.zc = a.zc + t * (b.zc - a.zc);
        if (refine && a.core === undefined && b.core === undefined) {
          const phiT = MT * (tr.s + 1);
          const g = q => { const r = shoot(tr.f, tr.s, tr.toC(q), { a: tr.a }); return r ? Math.log(r.phi / phiT) : NaN; };
          const q = N.brent(g, a.q, b.q, { tol: 1e-12 });
          if (isFinite(q)) { const r = shoot(tr.f, tr.s, tr.toC(q), { a: tr.a }); if (r) e = { eta: r.eta, c0: tr.toC(q), q }; }
        } else if (refine && a.core !== undefined && b.core !== undefined) {
          const phiT = MT * (tr.s + 1), n = +tr.o.n;
          const g = lz => { const r = shoot(tr.f, tr.s, 0, { zc: Math.exp(lz), deadN: n }); return r ? Math.log(r.phi / phiT) : NaN; };
          const lz = N.brent(g, Math.log(a.zc), Math.log(b.zc), { tol: 1e-12 });
          if (isFinite(lz)) { const r = shoot(tr.f, tr.s, 0, { zc: Math.exp(lz), deadN: n }); if (r) e = { eta: r.eta, c0: 0, core: Math.exp(lz) / r.phi, zc: Math.exp(lz) }; }
        }
        out.push(e);
      }
    }
    if (!out.length) {
      if (MT < P[0].MT) out.push({ eta: P[0].eta, c0: P[0].c0 });
      else if (MT > P[P.length - 1].MT) {
        /* beyond the traced range: η·M_gen → 1 with a deficit that falls off as 1/M_T */
        const L = P[P.length - 1], a = (1 - L.eta * L.Mgen) * L.MT;
        out.push({ eta: (1 - a / MT) / (tr.g * MT), c0: 0, core: L.core !== undefined ? 1 - (1 - L.core) * L.MT / MT : undefined, asym: true });
      }
    }
    return out;
  };
  /* solve target(pt) = 0 on the traced curve between neighbouring points a and b (exactly, by shooting) */
  PEL.curveSolve = function (tr, a, b, target) {
    if (a.core === undefined && b.core === undefined) {
      const at = q => { const r = shoot(tr.f, tr.s, tr.toC(q), { a: tr.a }); return r ? { phi: r.phi, eta: r.eta, MT: r.phi / (tr.s + 1), c0: tr.toC(q), q } : null; };
      const g = q => { const p = at(q); return p ? target(p) : NaN; };
      const q = N.brent(g, a.q, b.q, { tol: 1e-12 });
      return isFinite(q) ? at(q) : null;
    }
    if (a.core !== undefined && b.core !== undefined) {
      const n = +tr.o.n;
      const at = lz => { const r = shoot(tr.f, tr.s, 0, { zc: Math.exp(lz), deadN: n }); return r ? { phi: r.phi, eta: r.eta, MT: r.phi / (tr.s + 1), c0: 0, core: Math.exp(lz) / r.phi, zc: Math.exp(lz) } : null; };
      const g = lz => { const p = at(lz); return p ? target(p) : NaN; };
      const lz = N.brent(g, Math.log(a.zc), Math.log(b.zc), { tol: 1e-12 });
      return isFinite(lz) ? at(lz) : null;
    }
    return null;
  };
  /* all points on the curve where target changes sign (e.g. a given Weisz modulus) */
  PEL.findAll = function (tr, target) {
    const P = tr.pts, out = [];
    for (let i = 1; i < P.length; i++) {
      const fa = target(P[i - 1]), fb = target(P[i]);
      if (isFinite(fa) && isFinite(fb) && fa * fb <= 0 && !(fa === 0 && i > 1)) {
        const p = PEL.curveSolve(tr, P[i - 1], P[i], target);
        if (p) { p.Mgen = p.MT * tr.g; out.push(p); }
        else { const t = fa / (fa - fb), a = P[i - 1], b = P[i]; out.push({ MT: a.MT + t * (b.MT - a.MT), eta: a.eta + t * (b.eta - a.eta), c0: a.c0 + t * (b.c0 - a.c0), Mgen: (a.MT + t * (b.MT - a.MT)) * tr.g }); }
      }
    }
    return out;
  };

  /* concentration profile c(r/R) for one solution returned by etaAt */
  PEL.profile = function (tr, e, MT) {
    const phiT = MT * (tr.s + 1);
    let r = null;
    if (e.asym) {
      /* thin reaction zone: exponential penetration with the generalized modulus, as a picture */
      const xs = N.linspace(0, 1, 240), lam = phiT / (tr.g * Math.SQRT2);
      return { x: xs, c: xs.map(x => Math.exp(-Math.min(700, lam * (1 - x)))), approx: true };
    }
    if (e.zc !== undefined) r = shoot(tr.f, tr.s, 0, { zc: e.zc, deadN: +tr.o.n, full: true });
    else r = shoot(tr.f, tr.s, e.c0, { a: tr.a, full: true });
    if (!r) return null;
    const Z = r.phi, x = [], c = [];
    if (e.zc !== undefined) { x.push(0, e.zc / Z); c.push(0, 0); }
    else { x.push(0); c.push(e.c0); }
    r.z.forEach((z, i) => { x.push(z / Z); c.push(Math.min(Math.max(r.c[i], 0), 1)); });
    return { x, c };
  };

  /* isothermal pellet at a fixed modulus by Brent on the centre value (unique for monotone kinetics).
     o: kinetics + shape; opt: {rhoC: inert core fraction of the radius (egg-shell), extra: [w_j], a, full} */
  PEL.solvePhi = function (o, MT, opt = {}) {
    const s = sOf(o.shape), f = o.f || PEL.rate(o), cmin = PEL.cmin(o), span = 1 - cmin, a = opt.a || 0;
    const rho = opt.rhoC || 0, phiT = MT * (s + 1);
    const toC = q => cmin + span / (1 + Math.exp(-q));
    const run = (q, full) => shoot(f, s, toC(q), { a, zIn: rho > 0 ? rho * phiT : 0, extra: opt.extra, full });
    const g = q => { const r = run(q); return r ? Math.log(r.phi / phiT) : NaN; };
    let lo = -60, hi = 30, glo = g(lo);
    const ghi = g(hi);
    for (let k = 0; k < 8 && !(glo > 0); k++) { lo -= 60; glo = g(lo); }
    const shellFrac = rho > 0 ? 1 - Math.pow(rho, s + 1) : 1;
    if (!(glo > 0) || !(ghi < 0)) {
      const gg = PEL.genFactor(f, cmin);
      return { eta: 1 / (gg * MT), etaActive: 1 / (gg * MT * shellFrac), c0: cmin, asym: true, s, phi: phiT };
    }
    const q = N.brent(g, lo, hi, { fa: glo, fb: ghi, tol: 1e-12 });
    const res = run(q, opt.full);
    if (!res) return { eta: NaN, etaActive: NaN, s };
    res.c0 = toC(q);
    res.etaActive = res.eta / shellFrac;
    res.s = s;
    return res;
  };

  /* overall effectiveness with an external film, isothermal nth order, L = V/S:
     k_g (C_b − C_s) = L k C_s^n η(M_s)  ⇒  Bi (1 − c_s) = M_b² c_s^n η(M_s),  Bi = k_g L / D_e */
  PEL.withFilm = function (tr, n, MTb, Bi) {
    const g = cs => {
      const Ms = MTb * Math.pow(Math.max(cs, 1e-300), (n - 1) / 2);
      const e = PEL.etaAt(tr, Ms)[0].eta;
      return (1 - cs) * Bi - e * MTb * MTb * Math.pow(cs, n);
    };
    const cs = N.brent(g, 1e-14, 1, { tol: 1e-14 });
    const Ms = MTb * Math.pow(cs, (n - 1) / 2), e = PEL.etaAt(tr, Ms, true)[0].eta;
    return { cs, etaIn: e, etaO: e * Math.pow(cs, n), Ms };
  };

  /* parallel reactions in the pellet: rates k_j C^{a_j}; returns η (on the total surface rate) and
     the fraction of A going to the first path, observed and at the surface (true). */
  PEL.parallel = function (shape, Cs, ks, as, D, Lc) {
    const r1 = ks[0] * Math.pow(Cs, as[0]), r2 = ks[1] * Math.pow(Cs, as[1]), rs = r1 + r2;
    const w1 = c => (c <= 0 ? 0 : r1 * Math.pow(c, as[0]) / rs), w2 = c => (c <= 0 ? 0 : r2 * Math.pow(c, as[1]) / rs);
    const f = c => w1(c) + w2(c), s = sOf(shape);
    const MT = Lc * Math.sqrt(rs / (D * Cs)), phiT = MT * (s + 1), phiTrue = r1 / rs;
    const g = q => { const r = shoot(f, s, 1 / (1 + Math.exp(-q)), { extra: [w1] }); return r ? Math.log(r.phi / phiT) : NaN; };
    let lo = -60, glo = g(lo);
    for (let k = 0; k < 8 && !(glo > 0); k++) { lo -= 60; glo = g(lo); }
    const ghi = g(30);
    if (!(glo > 0) || !(ghi < 0)) {
      /* thin-zone limit: along the profile dζ = dc / sqrt(2F(c)), F(c) = ∫_0^c f */
      return Object.assign(PEL.parallelThin(w1, f), { MT, phiTrue, asym: true, eta: Math.sqrt(2 * N.quad(f, 0, 1)) / MT });
    }
    const q = N.brent(g, lo, 30, { fa: glo, fb: ghi, tol: 1e-12 });
    const r = shoot(f, s, 1 / (1 + Math.exp(-q)), { extra: [w1] });
    const Z = r.phi, tot = (s ? Math.pow(Z, s) : 1) * r.flux;
    return { MT, eta: r.eta, phiObs: r.P[0] / tot, phiTrue, c0: 1 / (1 + Math.exp(-q)) };
  };

  /* fraction of A taken by path 1 when the reaction zone is a thin outer skin (exact strong-diffusion limit) */
  PEL.parallelThin = function (w1, f) {
    const F = c => N.quad(f, 0, c, { rtol: 1e-9 });
    const tot = Math.sqrt(2 * F(1));
    const num = N.quad(c => { const Fc = F(c); return Fc > 0 ? w1(c) / Math.sqrt(2 * Fc) : 0; }, 0, 1, { rtol: 1e-7 });
    return { phiObs: num / tot };
  };

  /* linear series A → R → S inside a pellet (first order both, diffusivities D_A, D_R):
     net rates per pellet volume, linear in the surface concentrations:
       −r_A = k1 η1 C_As;   r_R = a C_As − k2 η2 C_Rs  */
  PEL.series = function (shape, L, k1, k2, DA, DR) {
    const M1 = L * Math.sqrt(k1 / DA), M2 = L * Math.sqrt(k2 / DR);
    const e1 = PEL.eta1(shape, M1), e2 = PEL.eta1(shape, M2);
    const den = k2 - k1 * DR / DA;
    let aR;
    if (Math.abs(den) > 1e-7 * k2) {
      const alpha = k1 / den;
      aR = k1 * e1 - k2 * alpha * (e1 - e2);
    } else {
      const h = 1e-4;
      aR = 0.5 * (PEL.series(shape, L, k1, k2 * (1 + h), DA, DR).a + PEL.series(shape, L, k1, k2 * (1 - h), DA, DR).a);
    }
    return { M1, M2, e1, e2, a: aR, kA: k1 * e1, kR: k2 * e2 };
  };

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
