/* One reaction off the isothermal line: ΔH(T) and K(T) from heat capacities, equilibrium
   extent, rates on the conversion–temperature plane, the locus of maximum rates,
   energy-balance (operating) lines, and the classical optimum for adiabatic staging. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N;
  const TH = ST.TH = {};
  const R = 8.314462618, RL = 0.0820573661, K0 = 273.15;
  Object.assign(TH, { R, RL, K0 });
  TH.C2K = c => c + K0;
  TH.K2C = k => k - K0;
  TH.CAL = 4.184;
  TH.TSEC = { s: 1, min: 60, h: 3600, day: 86400 };

  /* ---------- thermochemistry; ΔCp = da + db·T + dc·T² + dd/T², J/(mol·K), T in K ---------- */
  TH.coef = function (o) {
    const c = { T0: o.T0 || 298.15, H0: +o.H0 || 0, da: +o.da || 0, db: +o.db || 0, dc: +o.dc || 0, dd: +o.dd || 0 };
    c.lnK0 = o.lnK0 !== undefined ? +o.lnK0 : o.G0 !== undefined ? -o.G0 / (R * c.T0) : 0;
    c.J = c.H0 - (c.da * c.T0 + c.db / 2 * c.T0 * c.T0 + c.dc / 3 * c.T0 ** 3 - c.dd / c.T0);
    c.flat = !c.da && !c.db && !c.dc && !c.dd;
    return c;
  };
  TH.dH = (c, T) => c.J + c.da * T + c.db / 2 * T * T + c.dc / 3 * T ** 3 - c.dd / T;
  TH.dCp = (c, T) => c.da + c.db * T + c.dc * T * T + c.dd / (T * T);
  const Fk = (c, T) => -c.J / T + c.da * Math.log(T) + c.db / 2 * T + c.dc / 6 * T * T + c.dd / (2 * T * T);
  TH.lnK = (c, T) => c.lnK0 + (Fk(c, T) - Fk(c, c.T0)) / R;
  TH.dG = (c, T) => -R * T * TH.lnK(c, T);
  TH.dS = (c, T) => (TH.dH(c, T) - TH.dG(c, T)) / T;

  /* ---------- equilibrium extent of one reaction ----------
     sp: [{nu, n0}] (nu < 0 for reactants); lnK thermodynamic;
     o.basis: 'gas' (ideal, P/Pstd), 'conc' (C_i = n_i·cs relative to 1 mol/L), 'x' (ideal solution, mole fractions).
     Works in a logit coordinate between the extent bounds so that conversions near 0 or 1 keep their digits. */
  TH.equil = function (sp, lnK, o = {}) {
    const basis = o.basis || 'gas', P = o.P || 1, Pstd = o.Pstd || 1, cs = o.cs || 1, nI = o.nI || 0;
    let lo = -Infinity, hi = Infinity;
    sp.forEach(s => {
      if (s.nu < 0) hi = Math.min(hi, s.n0 / -s.nu);
      else if (s.nu > 0) lo = Math.max(lo, -s.n0 / s.nu);
    });
    if (!isFinite(hi) || !isFinite(lo) || !(hi > lo)) return null;
    const span = hi - lo;
    const at = u => {
      const s = 1 / (1 + Math.exp(-u)), sc = 1 / (1 + Math.exp(u));
      const dLo = span * s, dHi = span * sc;
      const n = sp.map(x => (x.nu < 0 ? (x.n0 + x.nu * hi) - x.nu * dHi : x.nu > 0 ? (x.n0 + x.nu * lo) + x.nu * dLo : x.n0));
      return { n, xi: lo + dLo };
    };
    const g = u => {
      const a = at(u);
      let nt = nI;
      a.n.forEach(x => { nt += x; });
      let f = -lnK;
      sp.forEach((x, i) => {
        if (!x.nu) return;
        const ni = Math.max(a.n[i], 1e-320);
        const act = basis === 'gas' ? ni / nt * P / Pstd : basis === 'x' ? ni / nt : ni * cs;
        f += x.nu * Math.log(act);
      });
      return f;
    };
    let u;
    const ga = g(-740), gb = g(740);
    if (ga >= 0) u = -740; else if (gb <= 0) u = 740; else u = N.brent(g, -740, 740, { fa: ga, fb: gb, tol: 1e-12 });
    const a = at(u);
    const key = o.key !== undefined ? o.key : sp.findIndex(x => x.nu < 0);
    const nt = nI + N.sum(a.n);
    return { xi: a.xi, n: a.n, nt, X: (sp[key].n0 - a.n[key]) / sp[key].n0, Xc: a.n[key] / sp[key].n0, lo, hi };
  };

  /* ---------- one reaction on the X–T plane ----------
     p: {a,b,r,s (stoichiometry, A is the key reactant), rev, oA, oB (orders if irreversible),
         phase 'liq'|'gas', CA0 (liq), P atm, yA0 (gas), MB, MR, MS (feed moles per mol A),
         A1, E1 (k1 = A1·exp(−E1/RT)), Ktype 'thermo'|'arr2', c (TH.coef for the reaction as written,
         standard state 1 mol/L for liquids, 1 atm for gases), A2, E2,
         Cp (C'p, J per mol A fed per K), Tlo, Thi (K, allowable range), tsec (seconds per time unit)} */
  TH.model = function (p) {
    const a = p.a || 1, b = p.b || 0, r = p.r || 0, s = p.s || 0, dn = r + s - a - b;
    const gas = p.phase === 'gas', rev = !!p.rev;
    const eps = gas ? p.yA0 * dn / a : 0;
    const MB = +p.MB || 0, MR = +p.MR || 0, MS = +p.MS || 0;
    const CA0f = gas ? (T => p.yA0 * p.P / (RL * T)) : (() => +p.CA0);
    const Xmax = b > 0 ? Math.min(1, MB * a / b) : 1;
    const k1 = T => p.A1 * Math.exp(-p.E1 / (R * T));
    const lnKc = !rev ? () => Infinity : p.Ktype === 'arr2' ? (T => Math.log(p.A1 / p.A2) - (p.E1 - p.E2) / (R * T)) : (T => TH.lnK(p.c, T) - (gas ? dn * Math.log(RL * T) : 0));
    const k2 = T => (rev ? k1(T) * Math.exp(-lnKc(T)) : 0);
    const dHfix = p.dHfix !== undefined && isFinite(+p.dHfix) ? +p.dHfix : NaN;
    const dHrx = isFinite(dHfix) ? (() => dHfix) : rev && p.Ktype === 'arr2' ? (T => (p.E1 - p.E2) + (gas ? dn * R * T : 0)) : (T => TH.dH(p.c, T));
    const dHA = T => dHrx(T) / a;
    const flatH = isFinite(dHfix) ? true : rev && p.Ktype === 'arr2' ? !gas || dn === 0 : p.c.flat;
    const conc = (X, T) => {
      const c0 = CA0f(T), den = 1 + eps * X;
      return [c0 * (1 - X) / den, c0 * (MB - b / a * X) / den, c0 * (MR + r / a * X) / den, c0 * (MS + s / a * X) / den];
    };
    const pw = (x, e) => (e === 0 ? 1 : e === 1 ? Math.max(x, 0) : Math.pow(Math.max(x, 0), e));
    const oA = rev ? a : p.oA !== undefined ? +p.oA : a, oB = rev ? b : p.oB !== undefined ? +p.oB : b;
    const rate = (X, T) => {
      const C = conc(X, T);
      const f = pw(C[0], oA) * (b ? pw(C[1], oB) : 1);
      if (!rev) return a * k1(T) * f;
      const back = (r ? pw(C[2], r) : 1) * (s ? pw(C[3], s) : 1);
      return a * k1(T) * (f - back * Math.exp(-lnKc(T)));
    };
    const L = x => Math.log(Math.max(x, 1e-300));
    const gEq = (X, T) => {
      const C = conc(X, T);
      return a * L(C[0]) + (b ? b * L(C[1]) : 0) - (r ? r * L(C[2]) : 0) - (s ? s * L(C[3]) : 0) + lnKc(T);
    };
    const Xe = T => {
      if (!rev) return Xmax;
      const top = Xmax * (1 - 1e-14), g0 = gEq(0, T);
      if (!(g0 > 0)) return 0;
      const gt = gEq(top, T);
      if (gt >= 0) return top;
      return N.brent(X => gEq(X, T), 0, top, { fa: g0, fb: gt, tol: 1e-15 });
    };
    const Tlo = p.Tlo, Thi = p.Thi;
    const ToptIn = (X, lo, hi) => {
      const f = T => -rate(X, T);
      const sc = N.fminScan(f, lo, hi, 28);
      return sc.x;
    };
    const Topt = (X, lo = Tlo, hi = Thi) => ToptIn(X, lo, hi);
    const rmax = (X, lo, hi) => rate(X, Topt(X, lo, hi));
    const Cp = +p.Cp;
    /* energy balance through feed temperature T1 with heat Q added per mol A fed */
    const Xline = (T, T1, Q = 0) => (Cp * (T - T1) - Q) / (-dHA(T));
    const Tline = (X, T1, Q = 0) => {
      if (flatH) return T1 + (Q - dHA(T1) * X) / Cp;
      const h = T => Cp * (T - T1) - Q + dHA(T) * X;
      return N.brent(h, 1, 1e5, { tol: 1e-10 });
    };
    /* feed temperature whose adiabatic line passes through (X, T) */
    const T1for = (X, T, Q = 0) => T - (Q - dHA(T) * X) / Cp;
    const exo = dHA(0.5 * (Tlo + Thi)) < 0;
    return {
      p, a, b, r, s, dn, gas, rev, eps, Xmax, CA0f, k1, k2, lnKc, Kc: T => Math.exp(lnKc(T)), dHrx, dHA, flatH, exo,
      conc, rate, Xe, Topt, rmax, Cp, Xline, Tline, T1for, Tlo, Thi, tsec: p.tsec || 60
    };
  };

  /* V/F_A0 along a path T(X) */
  TH.vf = (m, Tof, Xa, Xb, rtol = 1e-7) => (Xb <= Xa ? 0 : N.quad(X => { const r = m.rate(X, Tof(X)); return r > 0 ? 1 / r : 1e300; }, Xa, Xb, { rtol }));

  /* optimum temperature progression in plug flow, with the profile along the reactor */
  TH.otp = function (m, Xf, o = {}) {
    const lo = o.lo || m.Tlo, hi = o.hi || m.Thi;
    if (!(Xf < m.Xe(lo))) return null;
    const vf = N.quad(X => 1 / m.rmax(X, lo, hi), 0, Xf, { rtol: 1e-7 });
    const xs = N.linspace(0, Xf, o.n || 161);
    const Ts = xs.map(X => m.Topt(X, lo, hi)), inv = xs.map((X, i) => 1 / m.rate(X, Ts[i]));
    const cum = N.cumtrapz(xs, inv), tot = cum[cum.length - 1];
    const frac = cum.map(c => c / tot);
    /* where the path leaves the upper limit */
    let XB = NaN, fB = NaN;
    for (let i = 1; i < xs.length; i++) if (Ts[i - 1] >= hi - 0.02 && Ts[i] < hi - 0.02) {
      XB = N.brent(X => m.Topt(X, lo, hi) - (hi - 0.02), xs[i - 1], xs[i], { tol: 1e-9 });
      if (!isFinite(XB)) XB = xs[i];
      fB = N.quad(X => 1 / m.rmax(X, lo, hi), 0, XB, { rtol: 1e-7 }) / vf;
      break;
    }
    return { vf, X: xs, T: Ts, inv, frac, XB, fB };
  };

  /* highest feed temperature whose adiabatic line still reaches Xf short of equilibrium (and of Tcap) */
  TH.adiabRange = function (m, Xf, Q = 0, o = {}) {
    const Tcap = o.Tcap || Infinity;
    const okAt = T1 => { const T = m.Tline(Xf, T1, Q); return m.rate(Xf, T) > 0 && (m.exo ? T <= Tcap : T1 <= Tcap); };
    let hiT = o.hiT || m.Thi + 50, loT = o.loT || m.Tlo - 150;
    loT = Math.max(loT, 120);
    if (!okAt(loT)) return null;
    if (okAt(hiT)) return { lo: loT, hi: hiT };
    let a = loT, b = hiT;
    for (let i = 0; i < 60; i++) { const c = 0.5 * (a + b); if (okAt(c)) a = c; else b = c; }
    return { lo: loT, hi: a };
  };

  /* best adiabatic plug flow: feed temperature minimizing V/FA0 */
  TH.adiabBest = function (m, Xf, Q = 0, o = {}) {
    const rg = TH.adiabRange(m, Xf, Q, o);
    if (!rg) return null;
    const f = T1 => { const v = TH.vf(m, X => m.Tline(X, T1, Q), 0, Xf, 1e-6); return isFinite(v) && v < 1e290 ? Math.log(v) : 700; };
    const best = N.fminScan(f, rg.lo, rg.hi * (1 - 1e-9), o.n || 40);
    const T1 = best.x, vf = TH.vf(m, X => m.Tline(X, T1, Q), 0, Xf);
    return { T1, vf, Tout: m.Tline(Xf, T1, Q), range: rg };
  };

  /* adiabatic plug flow with recycle on the line through T1: best recycle ratio */
  TH.recycleBest = function (m, Xf, T1, Q = 0) {
    const Tof = X => m.Tline(X, T1, Q);
    const vfR = R => (R <= 0 ? TH.vf(m, Tof, 0, Xf, 1e-6) : (1 + R) * TH.vf(m, Tof, R / (1 + R) * Xf, Xf, 1e-6));
    const f = w => { const v = vfR(Math.pow(10, w)); return v > 0 && v < 1e290 ? Math.log(v) : 700; };
    const sc = N.fminScan(f, -4, 3.5, 46);
    const R = Math.pow(10, sc.x), v = vfR(R), v0 = vfR(0);
    if (v0 <= v) return { R: 0, vf: v0, Xin: 0, Tin: T1 };
    return { R, vf: v, Xin: R / (1 + R) * Xf, Tin: Tof(R / (1 + R) * Xf) };
  };

  /* N mixed reactors in series, each at the temperature of maximum rate for its outlet */
  TH.mfrSeries = function (m, Xf, Nst, o = {}) {
    const lo = o.lo || m.Tlo, hi = o.hi || m.Thi;
    const cache = new Map();
    const rm = X => { const k = X.toPrecision(12); if (!cache.has(k)) cache.set(k, m.rmax(X, lo, hi)); return cache.get(k); };
    if (Nst === 1) return { X: [0, Xf], vf: Xf / rm(Xf), each: [Xf / rm(Xf)] };
    const xsOf = z => { const e = z.map(Math.exp).concat([1]), tot = N.sum(e); let c = 0; return [0].concat(e.map(v => { c += v; return Xf * c / tot; })); };
    const cost = z => { const xs = xsOf(z); let t = 0; for (let i = 1; i < xs.length; i++) { const r = rm(xs[i]); if (!(r > 0)) return 1e300; t += (xs[i] - xs[i - 1]) / r; } return t; };
    let z = new Array(Nst - 1).fill(0);
    let best = N.nelderMead(cost, z, { tol: 1e-12, maxIt: 900 * Nst, step: z.map(() => 0.6) });
    best = N.nelderMead(cost, best.x, { tol: 1e-14, maxIt: 900 * Nst, step: z.map(() => 0.15) });
    const xs = xsOf(best.x);
    const each = xs.slice(1).map((x, i) => (x - xs[i]) / rm(x));
    return { X: xs, vf: N.sum(each), each, T: xs.slice(1).map(x => m.Topt(x, lo, hi)) };
  };

  /* cooled plug flow with constant coolant temperature; Ua in W/(L·K) */
  TH.cooled = function (m, Xf, T1, Ua, Tc, o = {}) {
    const ua = Ua * m.tsec;
    const f = (v, y) => {
      const r = m.rate(y[0], y[1]);
      return [r, ((-m.dHA(y[1])) * r - ua * (y[1] - Tc)) / m.Cp];
    };
    const vmax = o.vmax || 1e6;
    const sol = N.ode(f, 0, [0, T1], vmax, { rtol: 1e-8, atol: 1e-10, event: (v, y) => y[0] - Xf, stop: (v, y) => !(m.rate(y[0], y[1]) > 1e-14 * m.rate(0, T1)) || y[1] > (o.Trun || 3000), minPts: 200 });
    const ok = sol.status === 'event';
    const vEnd = sol.t[sol.t.length - 1];
    let Thot = -Infinity, vHot = 0;
    sol.t.forEach((v, i) => { if (sol.y[i][1] > Thot) { Thot = sol.y[i][1]; vHot = v; } });
    return { ok, vf: ok ? vEnd : NaN, sol, Thot, vHot, Xend: sol.last()[0], Tend: sol.last()[1] };
  };

  /* Optimal adiabatic stages with interstage cooling (reversible exothermic reactions).
     Conditions: within each stage ∫ ∂(1/r)/∂T dX = 0 along its adiabatic line; between stages the
     rate leaving one stage equals the rate entering the next. Shooting on the first inlet temperature. */
  TH.stagesAd = function (m, Xf, Nst, o = {}) {
    const lo = o.lo || 200, hi = o.hi || m.Thi + 400;
    const dT = 0.01;
    const dinv = (X, T) => { const a = m.rate(X, T + dT), b = m.rate(X, T - dT); return -(a - b) / (2 * dT) / Math.pow(m.rate(X, T), 2); };
    const lineFrom = (Xs, Ts) => { const T1 = m.T1for(Xs, Ts); return X => m.Tline(X, T1); };
    const march = (Tin1, keep) => {
      let Xs = 0, Tin = Tin1;
      const st = [];
      for (let i = 0; i < Nst; i++) {
        const Tof = lineFrom(Xs, Tin);
        if (!(m.rate(Xs, Tin) > 0) || dinv(Xs, Tin) >= 0) return { fail: 'hot', i, st };
        const sol = N.ode((X, y) => [dinv(X, Tof(X))], Xs, [0], m.Xe(m.Tlo) * 0.999999, { rtol: 1e-7, atol: 1e-14, event: (X, y) => y[0], stop: (X) => !(m.rate(X, Tof(X)) > 0), minPts: 30 });
        if (sol.status !== 'event') return { fail: 'eq', i, st };
        const Xo = sol.t[sol.t.length - 1], To = Tof(Xo);
        st.push({ Xin: Xs, Tin, Xout: Xo, Tout: To, Tof });
        if (i === Nst - 1) break;
        const rt = m.rate(Xo, To), Tm = m.Topt(Xo, lo, hi);
        const g = T => m.rate(Xo, T) - rt;
        let a = Math.max(lo, Tm - 600);
        if (g(a) > 0) return { fail: 'cold', i, st };
        const Tn = N.brent(g, a, Tm, { tol: 1e-10 });
        Xs = Xo; Tin = Tn;
      }
      return { st, Xend: st[st.length - 1].Xout };
    };
    /* bracket the first inlet temperature */
    const Ts = N.linspace(Math.max(lo, m.Tlo - 120), m.Thi + 60, o.scan || 34);
    const res = Ts.map(T => { const r = march(T); return r.Xend !== undefined ? r.Xend - Xf : NaN; });
    let Tin1 = NaN;
    for (let i = 1; i < Ts.length; i++) {
      if (isFinite(res[i - 1]) && isFinite(res[i]) && res[i - 1] * res[i] <= 0) {
        Tin1 = N.brent(T => { const r = march(T); return r.Xend !== undefined ? r.Xend - Xf : (r.fail === 'hot' ? -1 : 1); }, Ts[i - 1], Ts[i], { tol: 1e-9 });
        break;
      }
    }
    if (!isFinite(Tin1)) return null;
    const r = march(Tin1);
    if (!r.st || r.st.length !== Nst) return null;
    r.st.forEach(s => { s.vf = TH.vf(m, s.Tof, s.Xin, s.Xout); });
    return { Tin1, stages: r.st, vf: N.sum(r.st.map(s => s.vf)) };
  };

  /* ---------- the X–T chart: rate contours, equilibrium, locus of maximum rates ---------- */
  TH.levels = function (zmin, zmax, n) {
    const out = [];
    const lo = Math.max(zmin, zmax - (n || 4));
    for (let e = Math.floor(lo) - 1; e <= Math.ceil(zmax); e++) [1, 2, 5].forEach(mm => { const z = e + Math.log10(mm); if (z >= lo && z <= zmax) out.push(z); });
    return out;
  };
  TH.chart = function (m, o = {}) {
    const tx = o.kelvin ? (T => T) : TH.K2C;
    const Ta = o.Ta || m.Tlo, Tb = o.Tb || m.Thi, nT = o.nT || 110, nX = o.nX || 90;
    const Xtop = Math.min(m.Xmax, o.Xtop || 1);
    const Ts = N.linspace(Ta, Tb, nT), Xs = N.linspace(0, Xtop * 0.998, nX);
    const scale = o.scale || 1;
    let zmin = Infinity, zmax = -Infinity;
    const Z = Xs.map(X => Ts.map(T => {
      const r = m.rate(X, T) * scale;
      if (!(r > 0)) return NaN;
      const z = Math.log10(r);
      if (z < zmin) zmin = z;
      if (z > zmax) zmax = z;
      return z;
    }));
    const lv = o.levels && o.levels.length ? o.levels.filter(x => x > 0).map(Math.log10) : TH.levels(zmin, zmax, o.decades || 4);
    const C = ST.Plot.contour(Ts.map(tx), Xs, Z, lv);
    const series = [], marks = [], placed = [];
    const xs0 = tx(Ta), xspan = Tb - Ta;
    /* label each contour at its top (which lies on the locus of maximum rates), skipping crowded ones */
    C.slice().reverse().forEach(c => {
      if (!c.x.length) return;
      series.push({ x: c.x, y: c.y, color: 'var(--muted)', w: 0.9, hover: false, legend: false, noRange: true });
      let bi = -1, by = -1;
      for (let j = 0; j < c.y.length; j++) if (isFinite(c.y[j]) && c.y[j] > by) { by = c.y[j]; bi = j; }
      if (bi < 0 || by > Xtop * 0.95) return;
      const nx = (c.x[bi] - xs0) / xspan, ny = by / Xtop;
      if (nx < 0.03 || nx > 0.97) return;
      if (placed.some(p => Math.hypot((p[0] - nx) * 1.6, p[1] - ny) < 0.075)) return;
      placed.push([nx, ny]);
      marks.push({ t: 'text', x: c.x[bi], y: by + 0.014 * Xtop, s: N.fmt(Math.pow(10, c.level), 2), anchor: 'middle', color: 'var(--muted)', noRange: true });
    });
    if (m.rev) {
      const Tg = N.linspace(Ta, Tb, 160);
      series.push({ name: 'equilibrium X_{Ae}', x: Tg.map(tx), y: Tg.map(m.Xe), color: 2, w: 2.2 });
    }
    /* locus of maximum rates (unconstrained within a wide window, drawn where it falls inside the chart) */
    if (m.rev && m.exo) {
      const Xl = N.linspace(0, Math.min(m.Xe(Ta), Xtop) * 0.995, 70);
      const Tl = Xl.map(X => m.Topt(X, Ta - 400 > 100 ? Ta - 400 : 100, Tb + 400));
      series.push({ name: 'locus of max rates', x: Tl.map(tx), y: Xl, color: 4, dash: '6 4', w: 1.8 });
    }
    return { series, marks, levels: lv.map(z => Math.pow(10, z)), x: { label: 'T', unit: o.kelvin ? 'K' : '°C', min: tx(Ta), max: tx(Tb) }, y: { label: 'X_A', min: 0, max: Xtop } };
  };

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
