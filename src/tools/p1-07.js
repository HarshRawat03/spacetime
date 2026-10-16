/* Chapter 7 — Design for parallel reactions */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, KIT = ST.KIT, D = ST.D;
  const NAMES = ['R', 'S', 'T'];

  function system(v) {
    const paths = [];
    for (let i = 0; i < 3; i++) {
      if (i === 2 && !v.use3) continue;
      paths.push({ name: NAMES[i], k: +v['k' + i], a: +v['a' + i], b: v.useB ? +v['b' + i] : 0, nu: v.useB ? +v['nu' + i] : 0 });
    }
    const P = (c, e) => (e === 0 ? 1 : Math.pow(Math.max(c, 0), e));
    const rates = (CA, CB) => paths.map(p => p.k * P(CA, p.a) * P(CB, p.b));
    const want = Math.min(+v.want, paths.length - 1);
    const phi = (CA, CB) => { const r = rates(CA, CB), s = N.sum(r); return s > 0 ? r[want] / s : 0; };
    return { paths, rates, phi, want, useB: !!v.useB };
  }

  /* plug flow from (CA0, CB0) until CA = CAf; returns products and tau */
  function pfr(sys, CA0, CB0, CAf, o = {}) {
    const n = sys.paths.length;
    const f = (t, y) => { const r = sys.rates(y[0], y[1]); const s = N.sum(r); return [o.holdA ? 0 : -s, o.holdB ? 0 : -N.sum(r.map((x, i) => x * sys.paths[i].nu)), s].concat(r); };
    const target = o.consume !== undefined ? o.consume : CA0 - CAf;
    let span = 1, sol;
    for (let it = 0; it < 60; it++) {
      sol = N.ode(f, 0, [CA0, CB0, 0].concat(new Array(n).fill(0)), span, { rtol: 1e-9, atol: 1e-13, event: (t, y) => y[2] - target, minPts: 20 });
      if (sol.status === 'event') break;
      span *= 4;
    }
    const y = sol.last();
    return { tau: sol.t[sol.t.length - 1], CA: y[0], CB: y[1], prod: y.slice(3), consumed: y[2], ok: sol.status === 'event', sol };
  }
  /* mixed flow at outlet CA = CAf */
  function mfr(sys, CA0, CB0, CAf) {
    const dA = CA0 - CAf;
    let CB = CB0;
    if (sys.useB && sys.paths.some(p => p.nu > 0)) {
      const g = cb => { const r = sys.rates(CAf, cb), s = N.sum(r); return (CB0 - cb) - dA * N.sum(r.map((x, i) => x * sys.paths[i].nu)) / s; };
      CB = N.root(g, 1e-12 * CB0, CB0, { n: 200 });
      if (!isFinite(CB)) return null;
    }
    const r = sys.rates(CAf, CB), s = N.sum(r);
    const tau = dA / s;
    return { tau, CA: CAf, CB, prod: r.map(x => x * tau) };
  }

  const pathInputs = () => {
    const L = [{ group: 'Paths (A → product)' }];
    const d = [[1, 1, 0.3, 1], [1, 0.5, 1.8, 1], [1, 2, 0, 0]];
    for (let i = 0; i < 3; i++) {
      const cond = i === 2 ? v => !!v.use3 : null;
      L.push({ k: 'k' + i, label: 'k for → ' + NAMES[i], v: d[i][0], min: 1e-6, max: 1e6, scale: 'log', if: cond || undefined });
      L.push({ k: 'a' + i, label: 'order in A, → ' + NAMES[i], v: d[i][1], min: 0, max: 4, step: 0.05, if: cond || undefined });
      L.push({ k: 'b' + i, label: 'order in B, → ' + NAMES[i], v: d[i][2], min: 0, max: 4, step: 0.05, if: i === 2 ? v => !!v.use3 && v.useB : v => v.useB });
      L.push({ k: 'nu' + i, label: 'B used per A, → ' + NAMES[i], v: d[i][3], min: 0, max: 4, step: 0.5, if: i === 2 ? v => !!v.use3 && v.useB : v => v.useB });
    }
    return L;
  };

  ST.tool({
    id: 'parallel', tag: 'H-701', ch: 7, sec: 'Ch. 7', title: 'Parallel reactions: fractional yield and the best contacting',
    lede: 'Up to three competing paths with any orders in A and B. The instantaneous fractional yield curve, overall yields for plug, mixed, staged and mixed-then-plug flow, side feeding of either reactant, and the separation-and-recycle limit.',
    keys: 'parallel reactions fractional yield phi selectivity product distribution contacting pattern side entry semibatch keep concentration high low trambouze denbigh',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'System' },
      { k: 'use3', t: 'check', label: 'Third path (→ T)', v: false },
      { k: 'useB', t: 'check', label: 'Second reactant B takes part', v: true },
      { k: 'want', t: 'select', label: 'Desired product', opts: [['0', 'R'], ['1', 'S'], ['2', 'T']], v: '0' }
    ].concat(pathInputs(), [
      { group: 'Feed and duty' },
      { k: 'CA0', label: '$C_{A0}$ after mixing', unit: '[C]', v: 10, min: 1e-6, max: 1e6, scale: 'log' },
      { k: 'CB0', label: '$C_{B0}$ after mixing', unit: '[C]', v: 10, min: 0, max: 1e6, if: v => v.useB },
      { k: 'CAf', label: 'Outlet $C_{Af}$', unit: '[C]', v: 1, min: 1e-9, max: 1e6, scale: 'log' },
      { k: 'Nt', label: 'Equal mixed tanks for the staged case', v: 3, min: 2, max: 20, step: 1 },
      { k: 'CBs', label: 'Side-fed B held at', unit: '[C]', v: 1, min: 1e-6, max: 1e6, scale: 'log', if: v => v.useB },
      { k: 'CAs', label: 'A entering the side-fed reactor', unit: '[C]', v: 19, min: 1e-6, max: 1e6, scale: 'log', if: v => v.useB, hint: 'The book takes 19 mol/L for Example 7.2(c), ignoring dilution by the B stream.' },
      { k: 'CAh', label: 'Side-fed A held at', unit: '[C]', v: 1, min: 1e-6, max: 1e6, scale: 'log', if: v => v.useB }
    ]),
    presets: [
      { name: 'Example 7.2', v: { use3: false, useB: true, want: '0', k0: 1, a0: 1.5, b0: 0.3, nu0: 1, k1: 1, a1: 0.5, b1: 1.8, nu1: 1, CA0: 10, CB0: 10, CAf: 1, CBs: 1, CAs: 19, CAh: 1, Nt: 3, _cu: 'mol/L', _tu: 'min' },
        book: { PhiP: { v: 0.32, tol: 0.02 }, CRP: { v: 2.86, tol: 0.02 }, PhiM: 0.5, CRM: 4.5, PhiSB: { v: 0.87, tol: 0.01 }, CRSB: { v: 7.85, tol: 0.01 } },
        note: 'Separate A and B streams of 20 mol/L at equal flow, 90% conversion of A: plug flow, mixed flow, and plug flow of A with B fed along the way so that C_B stays at 1.' },
      { name: 'Example 7.3', v: { use3: true, useB: false, want: '1', k0: 1, a0: 0, k1: 2, a1: 1, k2: 1, a2: 2, CA0: 2, CAf: 0.0001, Nt: 3, _cu: 'mol/L' },
        book: { bestM: { v: 0.667, tol: 0.01 }, bestP: { v: 0.864, tol: 0.01 }, phimax: 0.5 },
        note: 'A → R (rate 1), A → S (2C_A), A → T (C_A²), S wanted, C_{A0} = 2. Book: 0.33 and 0.43 mol S per mol A fed for the best mixed and plug reactors, 0.5 with separation and recycle.' },
      { name: 'Example 7.4', v: { use3: true, useB: false, want: '1', k0: 1, a0: 0, k1: 2, a1: 1, k2: 1, a2: 2, CA0: 2, CAf: 0.0001, Nt: 3 }, book: { CMP: { v: 0.886, tol: 0.01 } }, note: 'Mixed flow to C_A = 1, where φ peaks, then plug flow to the end: 0.5 + 0.386 mol/L of S.' },
      { name: 'Example 10.1 (Trambouze)', v: { use3: true, useB: false, want: '1', k0: 0.025, a0: 0, k1: 0.2, a1: 1, k2: 0.4, a2: 2, CA0: 1, CAf: 0.25, Nt: 4 }, book: { phimax: 0.5, Caopt: 0.25 }, note: 'Maximize the fractional yield of S: run every mixed reactor at C_A = 0.25, where φ(S/A) = 0.5. With 100 L/min this takes 750 L in total.' },
      { name: 'Problem 10.7 (max S produced)', v: { use3: true, useB: false, want: '1', k0: 0.025, a0: 0, k1: 0.2, a1: 1, k2: 0.4, a2: 2, CA0: 1, CAf: 0.0001, Nt: 4 }, note: 'Now maximize the production of S rather than its yield.' }
    ],
    eqs: [
      ['\\varphi\\left(\\frac{R}{A}\\right) = \\frac{dC_R}{-dC_A}', 'Eq. 7.7'],
      ['\\Phi_p = \\frac{-1}{C_{A0}-C_{Af}}\\int_{C_{A0}}^{C_{Af}}\\varphi\\,dC_A', 'Eq. 7.9', 'Plug flow: the area under the φ curve.'],
      ['\\Phi_m = \\varphi\\big|_{C_{Af}}', 'Eq. 7.10', 'Mixed flow: a rectangle.'],
      ['\\Phi_{N} = \\frac{\\varphi_1(C_{A0}-C_{A1})+\\varphi_2(C_{A1}-C_{A2})+\\cdots}{C_{A0}-C_{AN}}', 'Eq. 7.12'],
      ['C_{Rf} = \\Phi\\,(C_{A0}-C_{Af})', 'Eq. 7.13']
    ],
    method: ['Every scheme is integrated or solved on the full kinetics (not just φ), so stoichiometric use of B is followed exactly in the premixed cases. Side feeding is idealized as in the book: the side-fed reactant is held at a fixed concentration and dilution by its stream is ignored.', 'The mixed-then-plug optimum and the best single reactors are found by one-dimensional minimization; separation with recycle operates at the maximum of φ.'],
    run(v) {
      const sys = system(v), w = sys.want, nm = NAMES[w];
      const CA0 = +v.CA0, CB0 = sys.useB ? +v.CB0 : 0, CAf = +v.CAf;
      if (!(CAf < CA0)) throw new Error('Outlet C_A must be below the feed.');
      const kpi = [], warn = [], rows = [];
      const dA = CA0 - CAf;
      const P = pfr(sys, CA0, CB0, CAf);
      if (!P.ok) warn.push('Plug flow cannot reach C_Af: B runs out first.');
      const Mx = mfr(sys, CA0, CB0, CAf);
      const PhiP = P.prod[w] / P.consumed, PhiM = Mx ? Mx.prod[w] / dA : NaN;
      kpi.push({ k: 'PhiP', l: 'Φ, plug flow', v: PhiP, s: 'C_' + nm + ' = ' + N.fmt(P.prod[w], 4) }, { k: 'CRP', l: 'C_' + nm + ', plug flow', v: P.prod[w], u: v._cu });
      kpi.push({ k: 'PhiM', l: 'Φ, mixed flow', v: PhiM, s: Mx ? 'C_' + nm + ' = ' + N.fmt(Mx.prod[w], 4) : 'no solution' }, { k: 'CRM', l: 'C_' + nm + ', mixed flow', v: Mx ? Mx.prod[w] : NaN, u: v._cu });
      rows.push(['plug flow, premixed', P.tau, PhiP, P.prod[w]], ['mixed flow, premixed', Mx ? Mx.tau : NaN, PhiM, Mx ? Mx.prod[w] : NaN]);
      /* staged mixed */
      const Nt = Math.round(+v.Nt);
      const stage = tauEach => { let ca = CA0, cb = CB0, prod = new Array(sys.paths.length).fill(0); for (let i = 0; i < Nt; i++) { const g = c => { const m = mfr(sys, ca, cb, c); return m ? m.tau - tauEach : 1; }; const c = N.root(g, 1e-12 * ca, ca * (1 - 1e-12), { n: 160 }); if (!isFinite(c)) return null; const m = mfr(sys, ca, cb, c); prod = prod.map((p, j) => p + m.prod[j]); ca = c; cb = m.CB; } return { ca, prod }; };
      if (Mx) {
        const te = N.root(t => { const s = stage(t); return s ? s.ca - CAf : -1; }, Mx.tau / Nt * 1e-3, Mx.tau, { n: 80, log: true });
        const S = isFinite(te) ? stage(te) : null;
        if (S) { rows.push([Nt + ' equal mixed tanks', te * Nt, S.prod[w] / dA, S.prod[w]]); kpi.push({ k: 'PhiN', l: 'Φ, ' + Nt + ' equal tanks', v: S.prod[w] / dA }); }
      }
      /* mixed then plug with fixed outlet */
      const mpOf = c1 => { const m = mfr(sys, CA0, CB0, c1); if (!m) return null; const p = pfr(sys, c1, m.CB, CAf); return { tau: m.tau + p.tau, prod: m.prod.map((x, j) => x + p.prod[j]) }; };
      const mpo = N.fminScan(c1 => { const r = mpOf(c1); return r ? -r.prod[w] : 0; }, CAf + 1e-9 * CA0, CA0 * (1 - 1e-9), 36);
      const MP = mpOf(mpo.x);
      if (MP) { rows.push(['mixed to ' + N.fmt(mpo.x, 3) + ', then plug', MP.tau, MP.prod[w] / dA, MP.prod[w]]); kpi.push({ k: 'CMP', l: 'C_' + nm + ', mixed then plug', v: MP.prod[w], u: v._cu, s: 'switch at C_A = ' + N.fmt(mpo.x, 3) }); }
      /* best single reactors with free outlet */
      let pm = { x: NaN, f: NaN };
      if (!sys.useB) {
        const bm = N.fminScan(c => { const m = mfr(sys, CA0, CB0, c); return m ? -m.prod[w] : 0; }, CA0 * 1e-6, CA0 * (1 - 1e-6), 120, { log: false });
        const full = pfr(sys, CA0, CB0, CA0 * 1e-7);
        const tt = N.linspace(0, full.tau, 400);
        let bpv = -1, bpc = NaN;
        tt.forEach(t => { const y = full.sol.at(t); if (y[3 + w] > bpv) { bpv = y[3 + w]; bpc = y[0]; } });
        kpi.push({ k: 'bestM', l: 'Most ' + nm + ' from one mixed reactor', v: -bm.f, u: v._cu, s: 'at C_A = ' + N.fmt(bm.x, 3) }, { k: 'bestP', l: 'Most ' + nm + ' from plug flow', v: bpv, u: v._cu, s: 'run to C_A = ' + N.fmt(bpc, 3) });
        pm = N.fminScan(c => -sys.phi(c, 0), CA0 * 1e-6, CA0, 300);
        kpi.push({ k: 'phimax', l: 'Separation and recycle: φ max', v: -pm.f, main: true }, { k: 'Caopt', l: 'Operate mixed flow at C_A', v: pm.x, u: v._cu });
        rows.push(['mixed at the φ maximum, unused A recycled', NaN, -pm.f, NaN]);
      } else {
        let best = { f: -1 };
        N.linspace(CA0 * 0.01, CA0, 40).forEach(ca => N.linspace(Math.max(CB0, 1e-9) * 0.01, Math.max(CB0, 1e-9), 40).forEach(cb => { const f = sys.phi(ca, cb); if (f > best.f) best = { f, ca, cb }; }));
        kpi.push({ k: 'phimax', l: 'Highest φ in the feed range', v: best.f, s: 'at C_A = ' + N.fmt(best.ca, 3) + ', C_B = ' + N.fmt(best.cb, 3) });
      }
      /* side feeding */
      if (sys.useB) {
        const CAs = +v.CAs;
        const SB = pfr(sys, CAs, +v.CBs, CAf, { holdB: true });
        const PhiSB = SB.prod[w] / SB.consumed;
        kpi.push({ k: 'PhiSB', l: 'Φ, plug A with B side-fed', v: PhiSB, main: true, s: 'B held at ' + N.fmt(+v.CBs, 3) }, { k: 'CRSB', l: 'C_' + nm + ' leaving (side-fed B)', v: PhiSB * dA, u: v._cu });
        rows.push(['plug A, B side-fed (C_B = ' + N.fmt(+v.CBs, 3) + ')', SB.tau, PhiSB, PhiSB * dA]);
        const SA = pfr(sys, +v.CAh, CB0, 0, { holdA: true, consume: dA });
        if (SA.ok) { const PhiSA = SA.prod[w] / SA.consumed; rows.push(['plug B, A side-fed (C_A = ' + N.fmt(+v.CAh, 3) + ')', SA.tau, PhiSA, PhiSA * dA]); kpi.push({ k: 'PhiSA', l: 'Φ, plug B with A side-fed', v: PhiSA }); }
      }
      /* phi curve */
      const cs = N.linspace(0, CA0, 240);
      const cbPath = sys.useB ? (P.sol ? null : null) : null;
      const phiP = c => { if (!sys.useB) return sys.phi(c, 0); const dc = CA0 - c; const nuAvg = sys.paths.length ? 1 : 0; return sys.phi(c, Math.max(CB0 - nuAvg * dc, 0)); };
      const series = [{ name: 'φ, premixed plug path', x: cs, y: cs.map(phiP), color: 1 }];
      if (sys.useB) series.push({ name: 'φ at C_B = ' + N.fmt(+v.CBs, 3), x: cs, y: cs.map(c => sys.phi(c, +v.CBs)), color: 6, dash: '5 4' });
      const sx = N.linspace(CAf, CA0, 80);
      series.push({ x: sx, y: sx.map(phiP), type: 'area', fill: 0.18, color: 3, hover: false, legend: false, noRange: true });
      const marks = [{ t: 'rect', x0: CAf, x1: CA0, y0: 0, y1: PhiM, color: 2, fill: 0.12, noRange: true }, { t: 'v', x: CAf, label: 'C_{Af}', color: 8 }];
      if (isFinite(pm.x)) marks.push({ t: 'pt', x: pm.x, y: -pm.f, color: 2, label: 'φ max' });
      const plots = [{ title: 'Instantaneous fractional yield of ' + nm + ' (Fig. 7.3)', x: { label: 'C_A', unit: v._cu, min: 0, max: CA0 }, y: { label: 'φ(' + nm + '/A)', min: 0, max: 1 }, series, marks, cap: 'Shaded area under the curve over the reactor’s range = plug flow product; rectangle = mixed flow product.' }];
      const ts = P.sol ? N.linspace(0, P.tau, 120) : [];
      if (ts.length) plots.push({ title: 'Plug flow profiles', x: { label: 'τ', unit: v._tu, min: 0 }, y: { label: 'C', unit: v._cu, min: 0 }, series: [{ name: 'C_A', x: ts, y: ts.map(t => P.sol.at(t)[0]), color: 1 }].concat(sys.useB ? [{ name: 'C_B', x: ts, y: ts.map(t => P.sol.at(t)[1]), color: 8 }] : [], sys.paths.map((p, j) => ({ name: 'C_' + p.name, x: ts, y: ts.map(t => P.sol.at(t)[3 + j]), color: j + 2 }))) });
      return { kpi, warn, plots, tables: [{ title: 'Schemes compared', cols: ['scheme', 'τ', 'Φ(' + nm + '/A)', 'C_' + nm + ' out'], rows }] };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
