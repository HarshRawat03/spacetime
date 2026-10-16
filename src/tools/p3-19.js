/* Chapter 19 — the packed bed catalytic reactor: staged adiabatic beds with intercooling, recycle,
   mixed stages on the optimum, cold-shot cooling, and the choice between them (Fig. 19.11 kinetics) */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, TH = ST.TH, XTI = ST.XTI;
  const fmt = (x, s = 4) => N.fmt(x, s);
  const MW = q => fmt(q / 1e6, 3) + ' MW';

  function pbModel(v) {
    const yA0 = +v.yA0;
    if (!(yA0 > 0 && yA0 <= 1)) throw new Error('The mole fraction of A in the feed must lie in (0, 1].');
    if (!(+v.Tmax > +v.Tmin)) throw new Error('The highest allowable temperature must be above the lowest.');
    return TH.model({
      a: 1, b: 0, r: 1, s: 0, rev: true, phase: 'gas', yA0, P: +v.P * 1e6 / 101325, MR: 0,
      A1: +v.A1 * 1000, E1: +v.E1, A2: +v.A2 * 1000, E2: +v.E2, Ktype: 'arr2', dHfix: +v.dH,
      Cp: +v.Cp / yA0, Tlo: +v.Tmin, Thi: +v.Tmax, tsec: 1
    });
  }
  /* W/F_A0 along an adiabatic line through (Xa, Ta) between Xa and Xb */
  const lineThrough = (m, Xa, Ta) => { const T1 = m.T1for(Xa, Ta); return X => m.Tline(X, T1); };

  /* ---- schemes; each returns {key, name, W, beds[], duties[], paths[], marks[], items (flowsheet), fail} ---- */
  function sOTP(m, c) {
    const r = TH.otp(m, c.Xf, { lo: c.Tmin, hi: c.Tmax });
    if (!r) return { key: 'otp', name: 'optimum temperature progression', fail: 'target beyond equilibrium' };
    const W = c.FA0 * r.vf, T0 = r.T[0], Te = r.T[r.T.length - 1];
    const d = [['preheat to ' + fmt(T0, 4) + ' K', c.FA0 * m.Cp * (T0 - c.Tfeed)], ['along the bed', c.FA0 * (m.Cp * (Te - T0) + m.dHA(T0) * c.Xf)], ['product cooler', c.FA0 * m.Cp * (c.Tprod - Te)]];
    return { key: 'otp', name: 'optimum temperature progression (ideal)', W, beds: [{ type: 'plug', W, Xin: 0, Xout: c.Xf, Tin: T0, Tout: Te }], duties: d, paths: [{ x: r.T, y: r.X }], items: [{ hx: MW(d[0][1]), after: fmt(T0, 4) + ' K' }, { tube: { prof: r.X.filter((x, i) => i % 16 === 0), tag: 'R-1901', label: 'wall-cooled', jacket: true, w: 150 }, after: fmt(Te, 4) + ' K' }, { hx: MW(d[2][1]), after: fmt(c.Tprod, 4) + ' K' }] };
  }
  function sAD(m, c, T1given) {
    let T1;
    if (isFinite(T1given)) T1 = T1given;
    else { const r = TH.adiabBest(m, c.Xf, 0, { Tcap: c.Tmax, n: 30 }); if (!r) return { key: 'ad', name: 'single adiabatic bed', fail: 'no adiabatic line reaches the target below T_max' }; T1 = r.T1; }
    const Tof = X => m.Tline(X, T1), vf = TH.vf(m, Tof, 0, c.Xf);
    if (!(vf < 1e250)) return { key: 'ad', name: 'single adiabatic bed', fail: 'the adiabatic line meets equilibrium first' };
    const W = c.FA0 * vf, To = Tof(c.Xf);
    const xs = N.linspace(0, c.Xf, 60);
    const d = [['preheat to ' + fmt(T1, 4) + ' K', c.FA0 * m.Cp * (T1 - c.Tfeed)], ['product cooler', c.FA0 * m.Cp * (c.Tprod - To)]];
    return { key: 'ad', name: 'single adiabatic bed, plug flow', W, Tmaxed: To, beds: [{ type: 'plug', W, Xin: 0, Xout: c.Xf, Tin: T1, Tout: To }], duties: d, paths: [{ x: xs.map(Tof), y: xs }], items: [{ hx: MW(d[0][1]), after: fmt(T1, 4) + ' K' }, { tube: { prof: xs.filter((x, i) => i % 6 === 0), tag: 'R-1901', label: 'adiabatic' }, after: fmt(To, 4) + ' K' }, { hx: MW(d[1][1]), after: fmt(c.Tprod, 4) + ' K' }] };
  }
  const vfRec = (m, Xf, T1, R) => { const Tof = X => m.Tline(X, T1); return R <= 0 ? TH.vf(m, Tof, 0, Xf, 1e-6) : (1 + R) * TH.vf(m, Tof, R / (1 + R) * Xf, Xf, 1e-6); };
  function sREC(m, c, T1g, Rg) {
    let T1 = T1g, R = Rg;
    if (!isFinite(T1) || !isFinite(R)) {
      const rg = TH.adiabRange(m, c.Xf, 0, { Tcap: c.Tmax });
      if (!rg) return { key: 'rec', name: 'single bed with recycle', fail: 'no line reaches the target below T_max' };
      if (!isFinite(T1)) {
        const f = t => { const b = isFinite(R) ? { vf: vfRec(m, c.Xf, t, R) } : TH.recycleBest(m, c.Xf, t); return b.vf > 0 && b.vf < 1e250 ? Math.log(b.vf) : 700; };
        T1 = N.fminScan(f, rg.lo, rg.hi * (1 - 1e-9), 18).x;
      }
      if (!isFinite(R)) R = TH.recycleBest(m, c.Xf, T1).R;
    }
    const vf = vfRec(m, c.Xf, T1, R);
    if (!(vf < 1e250)) return { key: 'rec', name: 'single bed with recycle', fail: 'the line meets equilibrium first' };
    const Tof = X => m.Tline(X, T1), Xi = R / (1 + R) * c.Xf, To = Tof(c.Xf), W = c.FA0 * vf;
    const xs = N.linspace(Xi, c.Xf, 40);
    const d = [['preheat fresh feed to ' + fmt(T1, 4) + ' K', c.FA0 * m.Cp * (T1 - c.Tfeed)], ['product cooler', c.FA0 * m.Cp * (c.Tprod - To)]];
    return { key: 'rec', name: 'single bed with recycle R = ' + fmt(R, 3), R, T1, W, beds: [{ type: 'rec', R, W, Xin: Xi, Xout: c.Xf, Tin: Tof(Xi), Tout: To }], duties: d, paths: [{ x: [T1, Tof(Xi)], y: [0, Xi], dash: '2 3' }, { x: xs.map(Tof), y: xs }], items: [{ hx: MW(d[0][1]), after: fmt(T1, 4) + ' K' }, { tube: { prof: xs.filter((x, i) => i % 5 === 0), tag: 'R-1901' }, after: fmt(To, 4) + ' K' }, { hx: MW(d[1][1]), after: fmt(c.Tprod, 4) + ' K' }], recycle: [1, 1], recLabel: 'R = ' + fmt(R, 3) };
  }
  function sMFR(m, c, Nst) {
    const r = TH.mfrSeries(m, c.Xf, Nst, { lo: c.Tmin, hi: c.Tmax });
    if (!r || !isFinite(r.vf)) return { key: 'mfr' + Nst, name: Nst + ' mixed', fail: 'no solution' };
    if (!r.T) r.T = r.X.slice(1).map(x => m.Topt(x, c.Tmin, c.Tmax));
    const beds = [], duties = [], items = [];
    let Tp = c.Tfeed;
    r.T.forEach((T, i) => {
      const Xa = r.X[i], Xb = r.X[i + 1], W = c.FA0 * r.each[i];
      const Q = c.FA0 * (m.Cp * (T - Tp) + m.dHA(T) * (Xb - Xa));
      beds.push({ type: 'mixed', W, Xin: Xa, Xout: Xb, Tin: Tp, Tout: T });
      duties.push(['exchanger in loop ' + (i + 1), Q]);
      items.push({ tank: { X: Xb, tag: 'R-19' + (i + 1), label: fmt(T, 4) + ' K', jacket: true }, after: 'X ' + fmt(Xb, 3), gap: 40 });
      Tp = T;
    });
    duties.push(['product cooler', c.FA0 * m.Cp * (c.Tprod - Tp)]);
    items.push({ hx: MW(duties[duties.length - 1][1]), after: fmt(c.Tprod, 4) + ' K' });
    const marks = r.T.map((T, i) => ({ t: 'pt', x: T, y: r.X[i + 1], color: 8, label: 'M' + (i + 1) }));
    return { key: 'mfr' + Nst, name: Nst === 1 ? 'one mixed bed on the optimum (R → ∞)' : Nst + ' mixed beds on the optimum (R → ∞)', W: c.FA0 * r.vf, beds, duties, paths: [], marks, items };
  }
  function sSTG(m, c, Nst) {
    const r = TH.stagesAd(m, c.Xf, Nst, { lo: Math.max(150, c.Tmin - 200), hi: c.Tmax + 400 });
    if (!r) return { key: 'stg' + Nst, name: Nst + ' adiabatic beds with intercooling', fail: 'no optimal staging found' };
    const beds = [], duties = [], items = [], paths = [];
    duties.push(['preheat to ' + fmt(r.stages[0].Tin, 4) + ' K', c.FA0 * m.Cp * (r.stages[0].Tin - c.Tfeed)]);
    items.push({ hx: MW(duties[0][1]), after: fmt(r.stages[0].Tin, 4) + ' K' });
    r.stages.forEach((s, i) => {
      const W = c.FA0 * s.vf;
      beds.push({ type: 'plug', W, Xin: s.Xin, Xout: s.Xout, Tin: s.Tin, Tout: s.Tout });
      const xs = N.linspace(s.Xin, s.Xout, 30);
      paths.push({ x: xs.map(s.Tof), y: xs });
      items.push({ tube: { prof: N.linspace(s.Xin, s.Xout, 8), tag: 'B-' + (i + 1), w: 72 }, after: fmt(s.Tout, 4) + ' K', gap: 44 });
      if (i < r.stages.length - 1) {
        const n = r.stages[i + 1];
        paths.push({ x: [s.Tout, n.Tin], y: [s.Xout, s.Xout], dash: '2 3' });
        const Q = c.FA0 * m.Cp * (n.Tin - s.Tout);
        duties.push(['intercooler ' + (i + 1), Q]);
        items.push({ hx: MW(Q), after: fmt(n.Tin, 4) + ' K' });
      }
    });
    const Te = r.stages[r.stages.length - 1].Tout;
    duties.push(['product cooler', c.FA0 * m.Cp * (c.Tprod - Te)]);
    items.push({ hx: MW(duties[duties.length - 1][1]), after: fmt(c.Tprod, 4) + ' K' });
    return { key: 'stg' + Nst, name: Nst + ' adiabatic beds with intercooling (Konoki optimum)', W: c.FA0 * r.vf, beds, duties, paths, items, Tmaxed: Math.max(...r.stages.map(s => s.Tout)) };
  }
  /* largest conversion reachable along a line before it meets equilibrium */
  const lineMax = (m, Tof, Xa) => {
    const g = X => m.rate(X, Tof(X));
    if (g(0.999999) > 0) return 0.999999;
    if (!(g(Xa) > 0)) return Xa;
    return N.brent(g, Xa, 0.999999, { tol: 1e-12 });
  };
  /* catalyst along a line with a smooth penalty when the line would cross equilibrium */
  const legCost = (m, Tof, Xa, Xb) => {
    const Xm = lineMax(m, Tof, Xa);
    if (Xb >= Xm - 1e-9) return { vf: NaN, pen: 30 + 300 * (Xb - Xm + 1e-9) };
    const vf = TH.vf(m, Tof, Xa, Xb, 1e-6);
    return vf < 1e250 ? { vf, pen: 0 } : { vf: NaN, pen: 30 };
  };
  const konokiStart = (m, c, Nst) => {
    const k = TH.stagesAd(m, c.Xf, Nst, { lo: Math.max(150, c.Tmin - 200), hi: c.Tmax + 400 });
    if (k) return k.stages;
    const ad = TH.adiabBest(m, c.Xf / Nst, 0, { Tcap: c.Tmax, n: 16 });
    const T = ad ? ad.T1 : 0.5 * (c.Tmin + c.Tmax);
    return N.range(Nst).map(i => ({ Xin: c.Xf * i / Nst, Xout: c.Xf * (i + 1) / Nst, Tin: T, Tout: T }));
  };
  /* staged beds each with recycle ratio R, intercooled; optimized by Nelder–Mead from Konoki's design */
  function sRST(m, c, Nst, R) {
    const Xf = c.Xf;
    const unpack = z => {
      const e = z.slice(0, Nst - 1).map(Math.exp).concat([1]), tot = N.sum(e);
      let cum = 0; const X = [0].concat(e.map(q => { cum += q; return Xf * cum / tot; }));
      return { X, Tin: z.slice(Nst - 1) };
    };
    const evalz = (z, keep) => {
      const { X, Tin } = unpack(z);
      let tot = 0, pen = 0; const st = [];
      for (let i = 0; i < Nst; i++) {
        const Xa = X[i], Xb = X[i + 1], Tof = lineThrough(m, Xa, Tin[i]);
        const Xm = (Xa + R * Xb) / (1 + R), To = Tof(Xb);
        if (To > c.Tmax) pen += (To - c.Tmax) * 0.02;
        const L = legCost(m, Tof, Xm, Xb);
        pen += L.pen;
        if (isFinite(L.vf)) tot += (1 + R) * L.vf;
        if (keep) st.push({ Xa, Xb, Xm, Tin: Tin[i], Tm: Tof(Xm), To, Tof, vf: (1 + R) * L.vf });
      }
      if (keep) return pen > 0 ? null : { st, vf: tot };
      return Math.log(Math.max(tot, 1e-300)) + pen;
    };
    const k0 = konokiStart(m, c, Nst);
    const w = k0.map(s => Math.max(s.Xout - s.Xin, 1e-6) / Xf), lw = w.map(q => Math.log(q / w[w.length - 1]));
    const z0 = lw.slice(0, Nst - 1).concat(k0.map(s => s.Tin));
    let best = N.nelderMead(z => evalz(z), z0, { tol: 1e-10, maxIt: 500 * Nst, step: z0.map((q, i) => (i < Nst - 1 ? 0.4 : 20)) });
    best = N.nelderMead(z => evalz(z), best.x, { tol: 1e-12, maxIt: 400 * Nst, step: z0.map((q, i) => (i < Nst - 1 ? 0.08 : 4)) });
    const r = evalz(best.x, true);
    if (!r) return { key: 'rst' + Nst, name: Nst + ' beds, each with recycle R = ' + fmt(R, 3), fail: 'no feasible design within T_max' };
    const beds = [], duties = [], items = [], paths = [];
    let Tp = c.Tfeed;
    r.st.forEach((s, i) => {
      const W = c.FA0 * s.vf;
      beds.push({ type: 'rec', R, W, Xin: s.Xm, Xout: s.Xb, Tin: s.Tm, Tout: s.To });
      const Q = c.FA0 * m.Cp * (s.Tin - Tp);
      duties.push([i ? 'intercooler ' + i : 'preheat to ' + fmt(s.Tin, 4) + ' K', Q]);
      items.push({ hx: MW(Q), after: fmt(s.Tin, 4) + ' K' }, { tube: { prof: N.linspace(s.Xm, s.Xb, 8), tag: 'B-' + (i + 1), label: 'R = ' + fmt(R, 3), w: 76 }, after: fmt(s.To, 4) + ' K', gap: 44 });
      const xs = N.linspace(s.Xm, s.Xb, 24);
      paths.push({ x: [s.Tin, s.Tm], y: [s.Xa, s.Xm], dash: '2 3' }, { x: xs.map(s.Tof), y: xs });
      if (i < r.st.length - 1) paths.push({ x: [s.To, r.st[i + 1].Tin], y: [s.Xb, s.Xb], dash: '2 3' });
      Tp = s.To;
    });
    duties.push(['product cooler', c.FA0 * m.Cp * (c.Tprod - Tp)]);
    items.push({ hx: MW(duties[duties.length - 1][1]), after: fmt(c.Tprod, 4) + ' K' });
    return { key: 'rst' + Nst, name: Nst + ' beds, each with recycle R = ' + fmt(R, 3), W: c.FA0 * r.vf, beds, duties, paths, items, Tmaxed: Math.max(...r.st.map(s => s.To)) };
  }
  /* cold-shot cooling: fresh cold feed split between stages; started from Konoki's intercooled design */
  function sCOLD(m, c, Nst) {
    const Xf = c.Xf, dH = -m.dHA(600), Cp = m.Cp;
    const unpack = z => {
      const fe = z.slice(0, Nst).map(Math.exp), ft = N.sum(fe), phi = fe.map(q => q / ft);
      return { phi, T1: z[Nst], xo: z.slice(Nst + 1).map(q => 1 / (1 + Math.exp(-q))) };
    };
    const run = (z, keep) => {
      const { phi, T1, xo } = unpack(z);
      let Fc = 0, Xp = 0, Tp = T1, tot = 0, pen = 0; const st = [];
      for (let i = 0; i < Nst; i++) {
        const Fn = Fc + phi[i];
        const Xin = i ? Xp * Fc / Fn : 0, Tin = i ? (Fc * Tp + phi[i] * c.Tfeed) / Fn : T1;
        const Tof = X => Tin + dH * (X - Xin) / Cp;
        const Xm = lineMax(m, Tof, Xin);
        const Xout = i === Nst - 1 ? Xf : Xin + (Math.min(Xf, Xm) - Xin) * xo[i];
        const To = Tof(Xout);
        if (To > c.Tmax) pen += (To - c.Tmax) * 0.02;
        const L = legCost(m, Tof, Xin, Xout);
        pen += L.pen;
        if (isFinite(L.vf)) tot += Fn * L.vf;
        if (keep) st.push({ phi: phi[i], Fn, Xin, Tin, Xout, To, Tof, vf: Fn * L.vf });
        Fc = Fn; Xp = Xout; Tp = To;
      }
      if (keep) return pen > 0 ? null : { st, vf: tot, T1 };
      return Math.log(Math.max(tot, 1e-300)) + pen;
    };
    /* start: same stage temperatures as the intercooled optimum, cooling done by cold feed instead */
    const k0 = konokiStart(m, c, Nst);
    const phi = [1];
    for (let i = 1; i < Nst; i++) {
      const Th = k0[i - 1].Tout, Tc = k0[i].Tin, f = Math.max(0.02, Math.min(0.8, (Th - Tc) / Math.max(Th - c.Tfeed, 1)));
      phi.push(f / (1 - f) * N.sum(phi));
    }
    const pt = N.sum(phi);
    const xo0 = k0.slice(0, Nst - 1).map(() => 0);
    const starts = [phi.map(q => Math.log(q / pt)).concat([k0[0].Tin]).concat(xo0)];
    starts.push(phi.map((q, i) => Math.log(i ? q / pt * 0.6 : q / pt)).concat([k0[0].Tin + 20]).concat(xo0));
    let bestv = Infinity, bestz = null;
    starts.forEach(z0 => {
      let b = N.nelderMead(z => run(z), z0, { tol: 1e-10, maxIt: 500 * Nst, step: z0.map((q, i) => (i === Nst ? 20 : 0.4)) });
      b = N.nelderMead(z => run(z), b.x, { tol: 1e-12, maxIt: 400 * Nst, step: z0.map((q, i) => (i === Nst ? 4 : 0.08)) });
      if (b.f < bestv) { bestv = b.f; bestz = b.x; }
    });
    const r = bestz ? run(bestz, true) : null;
    if (!r) return { key: 'cold' + Nst, name: Nst + ' beds with cold-shot cooling', fail: 'no feasible split within T_max and equilibrium' };
    const beds = [], duties = [], items = [], paths = [];
    duties.push(['preheat ' + fmt(r.st[0].phi * 100, 3) + '% of the feed to ' + fmt(r.T1, 4) + ' K', c.FA0 * r.st[0].phi * Cp * (r.T1 - c.Tfeed)]);
    items.push({ hx: MW(duties[0][1]), after: fmt(r.T1, 4) + ' K' });
    r.st.forEach((s, i) => {
      beds.push({ type: 'plug', W: c.FA0 * s.vf, Xin: s.Xin, Xout: s.Xout, Tin: s.Tin, Tout: s.To, phi: s.phi });
      const xs = N.linspace(s.Xin, s.Xout, 24);
      paths.push({ x: xs.map(s.Tof), y: xs });
      if (i < r.st.length - 1) { const n = r.st[i + 1]; paths.push({ x: [s.To, n.Tin], y: [s.Xout, n.Xin], dash: '2 3' }); }
      items.push({ tube: { prof: N.linspace(s.Xin, s.Xout, 8), tag: 'B-' + (i + 1), w: 70 }, after: (i < r.st.length - 1 ? '+' + fmt(r.st[i + 1].phi * 100, 3) + '% cold' : fmt(s.To, 4) + ' K'), gap: 52 });
    });
    const Te = r.st[r.st.length - 1].To;
    duties.push(['product cooler', c.FA0 * Cp * (c.Tprod - Te)]);
    items.push({ hx: MW(duties[duties.length - 1][1]), after: fmt(c.Tprod, 4) + ' K' });
    return { key: 'cold' + Nst, name: Nst + ' beds with cold-shot cooling', W: c.FA0 * r.vf, beds, duties, paths, items, Tmaxed: Math.max(...r.st.map(s => s.To)) };
  }

  const FIG1911 = { A1: Math.exp(7), E1: 83700, A2: Math.exp(18), E2: 167400, dH: -80000, Cp: 40, P: 2.494 };

  ST.tool({
    id: 'packedbed', tag: 'C-1901', ch: 19, sec: 'Ch. 19', title: 'Packed bed catalytic reactors: staging, intercooling, recycle and cold shot',
    lede: 'Design the adiabatic packed-bed system for an exothermic reversible catalytic reaction on the conversion–temperature chart of Fig. 19.11: one bed with the best feed temperature, one bed with the best recycle, mixed beds on the optimum with the rectangle rule, N beds with optimal intercooling (Konoki), N beds each with recycle, or N beds with cold-shot feed injection — each with its catalyst, heat duties and flowsheet, or all of them compared for one duty.',
    keys: 'packed bed catalytic reactor staged adiabatic intercooling interstage cooling cold shot recycle konoki horn optimum temperature progression mixed flow rectangles heat duty preheater cooler adiabatic slope choice of contacting fig 19.11 example 19.1 19.2 problem 19.13 19.14 19.15 19.16 19.17 19.18 19.19 19.20 19.21',
    slow: true,
    inputs: [
      { group: 'Feed and duty' },
      { k: 'FA0', label: 'F_A0', unit: 'mol/s', v: 100, min: 1e-6, max: 1e9, scale: 'log' },
      { k: 'yA0', label: 'Mole fraction of A in the feed (rest inert)', v: 0.125, min: 0.001, max: 1 },
      { k: 'Xf', label: 'Conversion wanted', v: 0.8, min: 0.01, max: 0.999 },
      { k: 'Tfeed', label: 'Feed available at', unit: 'K', v: 300, min: 100, max: 2000 },
      { k: 'Tprod', label: 'Product wanted at', unit: 'K', v: 300, min: 100, max: 2000 },
      { k: 'Tmax', label: 'Highest allowable temperature', unit: 'K', v: 900, min: 200, max: 3000 },
      { k: 'Tmin', label: 'Lowest temperature on the chart', unit: 'K', v: 300, min: 100, max: 2000 },
      { group: 'Arrangement' },
      { k: 'sch', t: 'select', label: 'Scheme', opts: [['all', 'compare every arrangement'], ['best1', 'best single bed (plug, recycle or mixed)'], ['otp', 'optimum temperature progression (ideal)'], ['ad', 'one adiabatic bed'], ['rec', 'one bed with recycle'], ['mfr', 'N mixed beds on the optimum (R → ∞)'], ['stg', 'N adiabatic beds with intercooling'], ['rst', 'N beds, each with recycle R'], ['cold', 'N beds with cold-shot cooling']], v: 'rec' },
      { k: 'Nst', label: 'Number of beds N', v: 2, min: 1, max: 5, step: 1, if: v => ['all', 'mfr', 'stg', 'rst', 'cold'].indexOf(v.sch) >= 0 },
      { k: 'adMode', t: 'seg', label: 'Adiabatic line', opts: [['best', 'best feed temperature'], ['T1', 'given feed temperature'], ['out', 'given outlet temperature']], v: 'best', if: v => v.sch === 'ad' || v.sch === 'rec' },
      { k: 'T1g', label: 'Feed temperature of the line', unit: 'K', v: 600, min: 100, max: 3000, if: v => (v.sch === 'ad' || v.sch === 'rec') && v.adMode === 'T1' },
      { k: 'Toutg', label: 'Outlet temperature', unit: 'K', v: 825, min: 100, max: 3000, if: v => (v.sch === 'ad' || v.sch === 'rec') && v.adMode === 'out' },
      { k: 'recMode', t: 'seg', label: 'Recycle ratio', opts: [['best', 'best'], ['given', 'given']], v: 'given', if: v => v.sch === 'rec' },
      { k: 'R', label: 'Recycle ratio R', v: 1, min: 0, max: 100, if: v => (v.sch === 'rec' && v.recMode === 'given') || v.sch === 'rst' },
      { group: 'Kinetics and thermodynamics (Fig. 19.11)' },
      { k: 'A1', label: 'k′₁ = A₁ e^{−E₁/RT}: A₁', unit: 'm³/(kg·s)', v: FIG1911.A1, min: 1e-20, max: 1e30, scale: 'log' },
      { k: 'E1', label: 'E₁', unit: 'J/mol', v: FIG1911.E1, min: 0, max: 1e6 },
      { k: 'A2', label: 'k′₂ = A₂ e^{−E₂/RT}: A₂', unit: 'm³/(kg·s)', v: FIG1911.A2, min: 1e-20, max: 1e30, scale: 'log' },
      { k: 'E2', label: 'E₂', unit: 'J/mol', v: FIG1911.E2, min: 0, max: 1e6 },
      { k: 'dH', label: 'ΔH_r used in the energy balance', unit: 'J/mol', v: FIG1911.dH, min: -1e6, max: 1e6 },
      { k: 'Cp', label: 'C_p of A, R and inert', unit: 'J/(mol·K)', v: FIG1911.Cp, min: 1, max: 1e4 },
      { k: 'P', label: 'Pressure π', unit: 'MPa', v: FIG1911.P, min: 1e-4, max: 100 }
    ],
    presets: [
      { name: 'Example 19.1 — one bed, 1 A : 7 inert, 80%', v: { FA0: 100, yA0: 0.125, Xf: 0.8, Tfeed: 300, Tprod: 300, Tmax: 900, sch: 'rec', adMode: 'T1', T1g: 600, recMode: 'given', R: 1 }, book: { W: { v: 3840, tol: 0.1 }, Tout: { v: 800, tol: 0.01 } }, note: 'The book draws the adiabatic from 600 K, reads 1/(−r′_A) ≈ 6 × 8 on it and chooses recycle R = 1: W/F_A0 = 38.4 kg·s/mol, W = 3840 kg; the product leaves at 800 K.' },
      { name: 'Example 19.1 — the best single bed', v: { FA0: 100, yA0: 0.125, Xf: 0.8, Tfeed: 300, Tprod: 300, Tmax: 900, sch: 'best1' }, note: 'Searching all feed temperatures and recycle ratios rather than reading the chart.' },
      { name: 'Example 19.2 — two mixed beds, pure A, 85%', v: { FA0: 100, yA0: 1, Xf: 0.85, Tfeed: 300, Tprod: 300, Tmax: 900, sch: 'mfr', Nst: 2 }, book: { W1: { v: 237.6, tol: 0.05 }, W2: { v: 381.9, tol: 0.05 }, Q1: { v: -3.2, tol: 0.05 }, Q2: { v: -1.8, tol: 0.06 }, Qlast: { v: -1.8, tol: 0.05 } }, note: 'The adiabatic for pure A is so shallow (slope 1/2000) that the rate rises all along it: use mixed beds on the optimum, sized by maximizing rectangles. Book: 0.66 at 820 K, then 0.85 at 750 K; 237.6 and 381.9 kg; −3.2, −1.8 and −1.8 MW.' },
      { name: 'Problem 19.13 — pure A, 85%, optimum progression', v: { FA0: 100, yA0: 1, Xf: 0.85, Tmax: 900, sch: 'otp' } },
      { name: 'Problem 19.14 — 50% A, 70%, one bed', v: { FA0: 100, yA0: 0.5, Xf: 0.7, Tmax: 900, sch: 'best1' } },
      { name: 'Problem 19.15 — 20% A, leaving at 75% and 825 K', v: { FA0: 100, yA0: 0.2, Xf: 0.75, Tmax: 900, sch: 'rec', adMode: 'out', Toutg: 825, recMode: 'best' } },
      { name: 'Problem 19.16 — 5% A, 50%, one bed', v: { FA0: 100, yA0: 0.05, Xf: 0.5, Tmax: 900, sch: 'best1' } },
      { name: 'Problem 19.17 — two beds, pure A, 85%', v: { FA0: 100, yA0: 1, Xf: 0.85, Tmax: 900, sch: 'all', Nst: 2 } },
      { name: 'Problem 19.18 — two beds, T_max = 550 K', v: { FA0: 100, yA0: 1, Xf: 0.85, Tmax: 550, sch: 'all', Nst: 2 } },
      { name: 'Problem 19.19 — two beds, 20% A, 85%', v: { FA0: 100, yA0: 0.2, Xf: 0.85, Tmax: 900, sch: 'all', Nst: 2 } },
      { name: 'Problem 19.20 — two beds, 40% A, 95%', v: { FA0: 100, yA0: 0.4, Xf: 0.95, Tmax: 900, sch: 'all', Nst: 2 } },
      { name: 'Problem 19.21 — 5% A, 95%, cold shot', v: { FA0: 100, yA0: 0.05, Xf: 0.95, Tmax: 900, sch: 'cold', Nst: 2 } }
    ],
    eqs: [
      ['-r_A^\\prime = k_1^\\prime C_A - k_2^\\prime C_R,\\qquad C_A = \\frac{y_{A0}\\pi}{RT}(1-X_A)', 'Fig. 19.11'],
      ['\\frac{W}{F_{A0}} = \\int_0^{X_A}\\frac{dX_A}{-r_A^\\prime}\\ \\text{along}\\ T = T_1 + \\frac{-\\Delta H_r}{C_p^\\prime}X_A', 'Ch. 9, 19', 'Adiabatic line; C′_p per mole of A fed.'],
      ['\\int_{X_{in}}^{X_{out}}\\left(\\frac{\\partial(1/-r_A^\\prime)}{\\partial T}\\right)_{X_A} dX_A = 0', 'Eq. 19.1', 'Where each optimal stage stops (Konoki).'],
      ['(-r_A^\\prime)_{\\text{leaving a stage}} = (-r_A^\\prime)_{\\text{entering the next}}', 'Eq. 19.2', 'How far to cool between stages.'],
      ['\\frac{W}{F_{A0}} = (R+1)\\int_{\\frac{R}{R+1}X_f}^{X_f}\\frac{dX_A}{-r_A^\\prime}', 'Eq. 6.21', 'Bed with recycle.'],
      ['\\text{slope of the adiabatic} = \\frac{C_p^\\prime}{-\\Delta H_r}', 'Ch. 19', 'Small slope (pure gas): recycle or mixed. Large slope (dilute gas, liquid): plug flow.']
    ],
    method: ['The kinetics of Fig. 19.11 are built from k′₁ and k′₂; the energy balance uses the chart’s ΔH_r = −80 kJ/mol and C_p = 40 J/mol·K for every species, so C′_p = 40/y_A0 per mole of A. Every catalyst amount is an exact quadrature of 1/(−r′_A) along the path.', 'Intercooled stages follow Konoki’s two conditions by shooting on the first inlet temperature. Mixed beds sit on the locus of maximum rates with the stage conversions chosen by maximizing rectangles. Beds with recycle and cold-shot beds are optimized by Nelder–Mead over stage conversions, inlet temperatures and feed split, with T_max as a penalty. The best single bed scans feed temperatures and recycle ratios.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const m = pbModel(v);
      const c = { FA0: +v.FA0, Xf: +v.Xf, Tfeed: +v.Tfeed, Tprod: +v.Tprod, Tmin: +v.Tmin, Tmax: +v.Tmax };
      const Xeq = m.Xe(c.Tmin);
      if (!(c.Xf < Xeq)) throw new Error('The conversion wanted is beyond equilibrium even at ' + fmt(c.Tmin, 4) + ' K (X_Ae = ' + fmt(Xeq, 4) + ').');
      const Nst = Math.max(1, Math.round(+v.Nst || 1));
      const slope = m.Cp / -m.dHA(600);
      kpi.push({ k: 'slope', l: 'Slope of the adiabatic dX/dT', v: slope, u: '1/K', s: 'ΔT for full conversion ' + fmt(1 / slope, 4) + ' K — ' + (1 / slope > 400 ? 'shallow: favour recycle or mixed beds' : 1 / slope < 150 ? 'steep: favour plug flow' : 'intermediate') });
      let list = [];
      const T1given = v.adMode === 'T1' ? +v.T1g : v.adMode === 'out' ? m.T1for(c.Xf, +v.Toutg) : NaN;
      if (v.sch === 'otp') list = [sOTP(m, c)];
      else if (v.sch === 'ad') list = [sAD(m, c, T1given)];
      else if (v.sch === 'rec') list = [sREC(m, c, T1given, v.recMode === 'given' ? +v.R : NaN)];
      else if (v.sch === 'mfr') list = [sMFR(m, c, Nst)];
      else if (v.sch === 'stg') list = [sSTG(m, c, Nst)];
      else if (v.sch === 'rst') list = [sRST(m, c, Nst, +v.R)];
      else if (v.sch === 'cold') list = [Nst < 2 ? sAD(m, c) : sCOLD(m, c, Nst)];
      else if (v.sch === 'best1') list = [sAD(m, c), sREC(m, c, NaN, NaN), sMFR(m, c, 1)];
      else list = [sOTP(m, c), sAD(m, c), sREC(m, c, NaN, NaN), sMFR(m, c, 1)].concat(Nst > 1 ? [sMFR(m, c, Nst), sSTG(m, c, Nst), sRST(m, c, Nst, 1), sCOLD(m, c, Nst)] : []);
      const ok = list.filter(s => !s.fail);
      list.filter(s => s.fail).forEach(s => warn.push(s.name + ': ' + s.fail + '.'));
      if (!ok.length) throw new Error('No arrangement meets this duty within the temperature limits.');
      const feasible = ok.filter(s => !(s.Tmaxed > c.Tmax + 0.5) && s.key !== 'otp');
      const pick = v.sch === 'best1' || v.sch === 'all' ? (feasible.length ? feasible : ok).reduce((a, b) => (b.W < a.W ? b : a)) : ok[0];
      kpi.unshift({ k: 'W', l: 'Catalyst, ' + pick.name, v: pick.W, u: 'kg', main: true, s: 'W/F_A0 = ' + fmt(pick.W / c.FA0, 4) + ' kg·s/mol' });
      pick.beds.forEach((b, i) => kpi.push({ k: 'W' + (i + 1), l: 'Bed ' + (i + 1) + (b.type === 'mixed' ? ' (mixed)' : b.type === 'rec' ? ' (recycle)' : ''), v: b.W, u: 'kg', s: 'X ' + fmt(b.Xin, 3) + ' → ' + fmt(b.Xout, 3) + ', T ' + fmt(b.Tin, 4) + ' → ' + fmt(b.Tout, 4) + ' K' + (b.phi !== undefined ? ', ' + fmt(b.phi * 100, 3) + '% of feed added' : '') }));
      if (pick.beds.length === 1) kpi.push({ k: 'Tout', l: 'Leaving the bed', v: pick.beds[0].Tout, u: 'K' });
      pick.duties.forEach((d, i) => kpi.push({ k: i === pick.duties.length - 1 ? 'Qlast' : 'Q' + (i + 1), l: 'Heat duty: ' + d[0], v: d[1] / 1e6, u: 'MW', s: d[1] < 0 ? 'cooling' : 'heating' }));
      if (pick.Tmaxed > c.Tmax + 0.5) warn.push('The hottest point, ' + fmt(pick.Tmaxed, 4) + ' K, is above T_max.');
      if (ok.length > 1) {
        const rows = ok.map(s => [s.name, s.W, s.W / c.FA0, s.beds.length, s.duties.reduce((t, d) => t + Math.abs(d[1]), 0) / 1e6, s.Tmaxed !== undefined ? s.Tmaxed : s.beds.reduce((t, b) => Math.max(t, b.Tout), 0), s.key === 'otp' ? 'ideal' : (s.Tmaxed > c.Tmax + 0.5 ? 'exceeds T_max' : s === pick ? 'chosen' : '')]).sort((a, b) => a[1] - b[1]);
        tables.push({ title: 'Arrangements for this duty, least catalyst first', cols: ['arrangement', 'W, kg', 'W/F_A0', 'beds', 'Σ|Q|, MW', 'hottest, K', ''], rows });
      }
      /* chart */
      const extra = [], marks = [];
      const colors = [3, 5, 6, 7, 8];
      (v.sch === 'all' || v.sch === 'best1' ? ok : [pick]).forEach((s, si) => {
        const col = s === pick ? 1 : colors[si % colors.length];
        (s.paths || []).forEach((p, i) => extra.push({ name: i === 0 ? s.name : undefined, x: p.x, y: p.y, color: col, w: s === pick ? 3 : 1.6, dash: p.dash, legend: i === 0 && !!s.name }));
        (s.marks || []).forEach(mk => marks.push(Object.assign({}, mk, { color: col })));
      });
      marks.push({ t: 'v', x: c.Tmax, color: 'var(--muted)', label: 'T_max' });
      plots.push(XTI.chartSpec(m, { kelvin: true, Ta: c.Tmin, Tb: Math.max(c.Tmax, c.Tmin + 100), ru: 'mol/(kg cat·s)', decades: 6, extra, marks, title: 'Conversion–temperature chart with the design (Fig. 19.11)', cap: 'Contours: −r′_A in mol/(kg cat·s). The chosen arrangement is drawn thick; dotted segments are mixing or cooling at constant conversion.' }));
      /* 1/(-r) along the chosen path */
      if (pick.beds.some(b => b.type !== 'mixed')) {
        const xs = [], ys = [];
        pick.beds.forEach(b => {
          const Tof = lineThrough(m, b.Xout, b.Tout);
          if (b.type === 'plug' || b.type === 'rec') N.linspace(b.Xin, b.Xout, 50).forEach(x => { xs.push(x); ys.push(1 / m.rate(x, pick.key === 'otp' ? m.Topt(x, c.Tmin, c.Tmax) : Tof(x))); });
          xs.push(NaN); ys.push(NaN);
        });
        plots.push({ title: '1/(−r′_A) along the chosen path', x: { label: 'X_A', min: 0, max: c.Xf * 1.02 }, y: { label: '1/(−r′_A)', unit: 'kg·s/mol', min: 0 }, series: [{ name: pick.name, x: xs, y: ys, color: 1, w: 2.2, type: 'area', line: true, fill: 0.14 }], cap: 'Shaded area times F_A0 (and times R + 1 for a recycle bed) is the catalyst.' });
      }
      if (slope < 1 / 400) info.push('Pure or rich gas: the adiabatic is shallow, the rate rises along it, so recycle or mixed beds avoid the slow cold end (Fig. 19.10).');
      if (v.sch === 'cold' || (v.sch === 'all' && Nst > 1)) info.push('Cold shot saves the interstage exchangers; it pays when the feed is much colder than the reaction and the temperature rise per stage is small (Ch. 19).');
      return { kpi, warn, info, plots, tables, diagram: pick.items ? XTI.flowsheet(pick.items, { feed: fmt(c.Tfeed, 4) + ' K', recycle: pick.recycle || null, recLabel: pick.recLabel, aria: pick.name }) : '' };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
