/* Chapter 16 — Earliness of mixing, segregation, and RTD */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, RT = ST.RT, KIT = ST.KIT, K = ST.K, D = ST.D;

  /* seeded random numbers (mulberry32) so Monte Carlo runs repeat exactly */
  const rng = seed => () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };

  /* nth-order batch step, exact */
  const stepN = (c, k, n, dt) => {
    if (c <= 0) return 0;
    if (Math.abs(n - 1) < 1e-12) return c * Math.exp(-k * dt);
    if (Math.abs(n) < 1e-12) return Math.max(0, c - k * dt);
    const v = Math.pow(c, 1 - n) + (n - 1) * k * dt;
    if (n < 1 && v <= 0) return 0;
    return Math.pow(v, 1 / (1 - n));
  };

  /* Curl's coalescence–redispersion model of a mixed reactor (structure of arrays, in-place updates) */
  function curl(o) {
    const R = rng(o.seed || 12345), Np = o.Np || 1000, tau = 1, dt = tau / (o.steps || 80);
    const nsp = o.nsp, C = N.range(nsp).map(() => new Float64Array(Np));
    const f0 = o.feed(0), f1 = o.feed(1);
    for (let i = 0; i < Np; i++) { const f = i % 2 ? f1 : f0; for (let q = 0; q < nsp; q++) C[q][i] = f[q]; }
    const tEnd = o.runs || 10, tAvg = 5;
    const exitN = Np * dt / tau, pairs = o.I / tau * Np / 2 * dt;
    const sum = new Float64Array(nsp), sum2 = new Float64Array(nsp);
    let nOut = 0, carry = 0, carryP = 0, fi = 0;
    for (let t = 0; t < tEnd; t += dt) {
      o.react(C, Np, dt);
      carry += exitN;
      const nx = Math.floor(carry); carry -= nx;
      const rec = t > tEnd - tAvg;
      for (let j = 0; j < nx; j++) {
        const i = Math.floor(R() * Np);
        if (rec) { for (let q = 0; q < nsp; q++) { sum[q] += C[q][i]; sum2[q] += C[q][i] * C[q][i]; } nOut++; }
        const f = fi++ % 2 ? f1 : f0;
        for (let q = 0; q < nsp; q++) C[q][i] = f[q];
      }
      carryP += pairs;
      const np = Math.floor(carryP); carryP -= np;
      for (let j = 0; j < np; j++) {
        const a = Math.floor(R() * Np); let b = Math.floor(R() * Np); if (b === a) b = (a + 1) % Np;
        for (let q = 0; q < nsp; q++) { const mval = 0.5 * (C[q][a] + C[q][b]); C[q][a] = mval; C[q][b] = mval; }
      }
    }
    const mean = Array.from(sum).map(x => x / Math.max(nOut, 1));
    const se = Array.from(sum2).map((x, q) => Math.sqrt(Math.max(0, x / Math.max(nOut, 1) - mean[q] * mean[q]) / Math.max(1, nOut / 40)));
    return { mean, se, n: nOut };
  }

  ST.tool({
    id: 'segregation', tag: 'F-1601', ch: 16, sec: 'Ch. 16', title: 'Segregation and earliness of mixing: micro- and macrofluids, early and late mixing',
    lede: 'When does it matter how fluid mixes, not just how long it stays? Micro- against macrofluid in mixed flow for any order (Fig. 16.2, Table 16.2), the earliest and latest mixing allowed by one RTD (Example 16.1), partial segregation by Curl’s coalescence–redispersion model, the loss of an intermediate when two feeds mix slowly, how long a fluid element keeps its identity, and the exponential-integral table.',
    keys: 'segregation micromixing macrofluid microfluid earliness of mixing early late mixing exponential integral ei Ei table 16.1 16.2 coalescence redispersion curl partial segregation intensity of mixing fast reactions two feeds competitive consecutive element lifetime example 16.1',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'Task' },
      { k: 'what', t: 'select', label: 'Show', opts: [['mfr', 'micro vs macro in mixed flow (Fig. 16.2)'], ['early', 'early vs late mixing in plug + mixed (Ex. 16.1)'], ['partial', 'partial segregation (coalescence model)'], ['twofeed', 'two feeds, fast series–parallel reactions'], ['life', 'how long an element keeps its identity'], ['ei', 'exponential integrals (Table 16.1)']], v: 'early' },
      { k: 'n', label: 'Reaction order n', v: 2, min: 0, max: 4, step: 0.05, if: v => v.what === 'mfr' || v.what === 'partial' },
      { k: 'Rg', label: 'Rate group kC_A0^{n−1}τ', v: 1, min: 1e-4, max: 1e4, scale: 'log', if: v => v.what === 'mfr' || v.what === 'partial' },
      { k: 'tP', label: 'Plug-flow space time', unit: '[t]', v: 1, min: 0, max: 1e9, if: v => v.what === 'early' },
      { k: 'tM', label: 'Mixed-flow space time', unit: '[t]', v: 1, min: 1e-9, max: 1e9, if: v => v.what === 'early' },
      { k: 'Ilist', t: 'text', label: 'Mixing intensities I = (collision frequency)·τ to try', v: '0, 1, 3, 10, 30, 100', mono: true, if: v => v.what === 'partial' || v.what === 'twofeed' },
      { k: 'k1f', label: 'k₁·C_B0·τ for A + B → R', v: 50, min: 1e-4, max: 1e6, scale: 'log', if: v => v.what === 'twofeed' },
      { k: 'k2f', label: 'k₂/k₁ for R + B → S', v: 0.5, min: 0, max: 100, if: v => v.what === 'twofeed' },
      { k: 'MBf', label: 'B fed per A fed (moles)', v: 1, min: 0.01, max: 100, if: v => v.what === 'twofeed' },
      { k: 'Np', label: 'Drops in the simulation', v: 600, min: 200, max: 20000, step: 100, if: v => v.what === 'partial' || v.what === 'twofeed' },
      { k: 'dEl', label: 'Element size', unit: 'm', v: 1e-6, min: 1e-10, max: 1, scale: 'log', if: v => v.what === 'life' },
      { k: 'Dif', label: 'Diffusion coefficient', unit: 'm²/s', v: 1e-9, min: 1e-16, max: 1e-3, scale: 'log', if: v => v.what === 'life' },
      { k: 'tauL', label: 'Compare with a residence time of', unit: 's', v: 60, min: 1e-9, max: 1e9, scale: 'log', if: v => v.what === 'life' },
      { k: 'xs', t: 'text', label: 'x values', v: '0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 1.4, 2', mono: true, if: v => v.what === 'ei' }
    ].concat(KIT.kin({ gas: false, d: { law: 'power', k: 1, a: 2, CA0: 1 } }).map(f => Object.assign({}, f, { if: f.group ? (vv => vv.what === 'early') : f.if ? (vv => vv.what === 'early' && f.if(vv)) : (vv => vv.what === 'early') }))),
    presets: [
      { name: 'Example 16.1', v: { what: 'early', tP: 1, tM: 1, law: 'power', k: 1, a: 2, ob: 0, CA0: 1 }, book: { CA_A: { v: 0.382, tol: 0.005 }, CA_B: { v: 0.366, tol: 0.005 }, CA_C: { v: 0.362, tol: 0.005 } }, note: 'Second order, k = C_A0 = 1, τ = 1 for each unit. Mixed then plug (earliest mixing): 0.382. Plug then mixed: 0.366. Macrofluid, or latest mixing: 0.362.' },
      { name: 'Fig. 16.2, second order', v: { what: 'mfr', n: 2, Rg: 10 }, note: 'Segregation helps a second-order reaction in mixed flow.' },
      { name: 'Fig. 16.2, zero order', v: { what: 'mfr', n: 0, Rg: 0.8 }, note: 'Segregation hurts a zero-order reaction in mixed flow.' },
      { name: 'Table 16.1', v: { what: 'ei', xs: '0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 1.4' }, book: { ei01: { v: 1.8229, tol: 0.001 }, Ei1: { v: 1.8951, tol: 0.001 } } },
      { name: 'Partial segregation (Curl)', v: { what: 'partial', n: 2, Rg: 10, Ilist: '0, 1, 3, 10, 30, 100', Np: 600 }, note: 'A second-order reaction in a mixed reactor whose drops meet I times per residence time: from the macrofluid (I = 0) to the microfluid (I → ∞).' },
      { name: 'Fast reactions, two feeds', v: { what: 'twofeed', k1f: 50, k2f: 0.5, MBf: 1, Ilist: '1, 3, 10, 30, 100, 300', Np: 600 }, note: 'A and B enter as separate streams. When reaction outruns mixing, R forms where B is plentiful and is attacked again: poor mixing lowers the yield of R (Section 16.2).' },
      { name: 'Element lifetime (water, 1 μm)', v: { what: 'life', dEl: 1e-6, Dif: 1e-9, tauL: 60 }, book: { tlife: { v: 1e-3, tol: 0.01 } }, note: 'A 1 μm element of water (𝒟 ≈ 10⁻⁵ cm²/s) loses its identity in about 10⁻³ s.' },
      { name: 'Element lifetime (viscous oil, 1 mm)', v: { what: 'life', dEl: 1e-3, Dif: 1e-11, tauL: 3600 }, book: { tlife: { v: 1e5, tol: 0.01 } }, note: 'A 1 mm element of a liquid 100 times as viscous as water keeps its identity for about 10⁵ s, roughly 30 h.' }
    ],
    eqs: [
      ['1 - \\bar X_A = \\int_0^\\infty \\left(\\frac{C_A}{C_{A0}}\\right)_{batch}\\frac{e^{-t/\\bar t}}{\\bar t}\\,dt', 'Eq. 16.4', 'Macrofluid in mixed flow.'],
      ['\\frac{\\bar C_A}{C_{A0}} = \\frac{1}{R}e^{1/R}\\,\\mathrm{ei}\\!\\left(\\frac{1}{R}\\right),\\quad R = kC_{A0}\\bar t', 'Eq. 16.8', 'Second order, macrofluid.'],
      ['\\frac{C_A}{C_{A0}} = \\frac{-1+\\sqrt{1+4R}}{2R}', 'Eq. 16.9', 'Second order, microfluid.'],
      ['\\frac{\\bar C_A}{C_{A0}} = 1 - R\\left(1 - e^{-1/R}\\right)', 'Eq. 16.5', 'Zero order, macrofluid.'],
      ['\\mathrm{ei}(x) = \\int_x^\\infty \\frac{e^{-u}}{u}du,\\qquad \\mathrm{Ei}(x) = \\int_{-\\infty}^x \\frac{e^u}{u}du', 'Table 16.1'],
      ['t \\approx \\frac{(\\text{element size})^2}{\\mathscr{D}}', 'Eq. 16.13']
    ],
    method: ['Macrofluid results integrate the batch history over the RTD; microfluid results solve the unit balances (mixed then plug, plug then mixed). For a plug–mixed pair, mixed-first is the earliest mixing the RTD allows and segregation the latest.', 'Partial segregation uses Curl’s coalescence–redispersion model as a Monte Carlo simulation (Spielman and Levenspiel, 1965): a population of drops reacts as small batches, drops leave at random (an exponential RTD), and pairs meet and equalize I times per residence time. The random sequence is seeded, so results repeat; statistical scatter is one or two percent with a few hundred drops.', 'Exponential integrals use a power series for small arguments and a continued fraction for large ones.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      if (v.what === 'ei') {
        const xs = String(v.xs).split(/[,;\s]+/).map(parseFloat).filter(x => x > 0).slice(0, 40);
        tables.push({ title: 'Exponential integrals (Table 16.1)', cols: ['x', 'Ei(x)', 'ei(x) = E₁(x)'], rows: xs.map(x => [x, N.Ei(x), N.E1(x)]) });
        const i01 = xs.indexOf(0.1), i1 = xs.indexOf(1);
        if (i01 >= 0) kpi.push({ k: 'ei01', l: 'ei(0.1)', v: N.E1(0.1) });
        if (i1 >= 0) kpi.push({ k: 'Ei1', l: 'Ei(1)', v: N.Ei(1) });
        kpi.push({ k: 'ei2', l: 'ei(2), used in Example 16.1', v: N.E1(2) });
        const xg = N.logspace(0.005, 3, 200);
        plots.push({ title: 'Ei and ei', x: { label: 'x', log: true }, y: { label: 'value', min: -5, max: 5 }, series: [{ name: 'Ei(x)', x: xg, y: xg.map(N.Ei), color: 1 }, { name: 'ei(x)', x: xg, y: xg.map(N.E1), color: 2 }] });
        return { kpi, warn, info, plots, tables };
      }
      if (v.what === 'life') {
        const t = Math.pow(+v.dEl, 2) / +v.Dif;
        kpi.push({ k: 'tlife', l: 'Time to lose identity, d²/𝒟', v: t, u: 's', main: true, s: t > 3600 ? N.fmt(t / 3600, 3) + ' h' : t > 60 ? N.fmt(t / 60, 3) + ' min' : '' });
        const r = t / +v.tauL;
        kpi.push({ k: 'ratio', l: 'Against the residence time', v: r, s: r > 1 ? 'behaves as a macrofluid' : r < 0.01 ? 'behaves as a microfluid' : 'partly segregated' });
        info.push('Large elements are first stretched and folded by the flow (turbulence, baffles, stirrers); diffusion then erases the thin striations. This estimate is the diffusion step alone.');
        return { kpi, warn, info, plots, tables };
      }
      if (v.what === 'mfr') {
        const n = +v.n, Rg = +v.Rg;
        const micro = R => { if (Math.abs(n) < 1e-12) return Math.max(0, 1 - R); if (Math.abs(n - 1) < 1e-12) return 1 / (1 + R); return N.brent(c => 1 - c - R * Math.pow(Math.max(c, 0), n), 0, 1, { tol: 1e-14 }); };
        const macro = R => N.quad(th => stepN(1, R, n, th) * Math.exp(-th), 0, Infinity, { rtol: 1e-10 });
        const ci = micro(Rg), ca = macro(Rg);
        kpi.push({ k: 'Cmicro', l: 'Microfluid C_A/C_A0', v: ci, main: true }, { k: 'Cmacro', l: 'Macrofluid C_A/C_A0', v: ca, main: true }, { k: 'Cplug', l: 'Plug flow C_A/C_A0', v: stepN(1, Rg, n, 1) });
        if (Math.abs(n - 2) < 1e-12) kpi.push({ k: 'Ceq8', l: 'Eq. 16.8 with ei', v: Math.exp(1 / Rg) * N.E1(1 / Rg) / Rg });
        info.push(n > 1 ? 'n > 1: segregation raises conversion.' : n < 1 ? 'n < 1: segregation lowers conversion.' : 'First order: segregation has no effect.');
        const Rs = N.logspace(0.01, 1000, 80);
        plots.push({ title: 'Mixed flow, micro- and macrofluid (Fig. 16.2), n = ' + N.fmt(n, 3), x: { label: 'kC_A0^{n−1}τ', log: true }, y: { label: 'C_A/C_A0', log: true, min: 1e-3, max: 1 }, series: [{ name: 'microfluid', x: Rs, y: Rs.map(micro), color: 1 }, { name: 'macrofluid', x: Rs, y: Rs.map(macro), color: 2 }, { name: 'plug flow', x: Rs, y: Rs.map(R => stepN(1, R, n, 1)), color: 8, dash: '5 4' }], marks: [{ t: 'v', x: Rg, color: 8 }] });
        return { kpi, warn, info, plots, tables };
      }
      if (v.what === 'early') {
        const m = KIT.model(v), tP = +v.tP, tM = +v.tM;
        const XA = K.XPFR(m, tP, K.XMFR(m, tM)).X;
        const XB = K.XMFR(m, tM, K.XPFR(m, tP).X);
        const Xb = RT.batchX(m, tP + 40 * tM);
        const Xseg = N.quad(u => Xb(tP + u) * RT.Emixed(u, tM), 0, Infinity, { rtol: 1e-11 });
        const Lam = l => (l > tM ? 1 / tM : 0);
        /* maximum mixedness for this RTD equals mixed-then-plug */
        kpi.push({ k: 'CA_A', l: 'Mixed then plug (earliest mixing)', v: 1 - XA, main: true }, { k: 'CA_B', l: 'Plug then mixed (later mixing)', v: 1 - XB, main: true }, { k: 'CA_C', l: 'Macrofluid, or latest mixing', v: 1 - Xseg, main: true });
        const n = v.law === 'power' ? +v.a : NaN;
        info.push('All three share one RTD: a delay of ' + N.fmt(tP, 3) + ' followed by an exponential tail. ' + (isFinite(n) ? (n > 1 ? 'For n > 1, later mixing and segregation give more conversion.' : n < 1 ? 'For n < 1, early mixing gives more conversion.' : 'First order: all three agree.') : ''));
        tables.push({ title: 'One RTD, three mixing histories (Example 16.1)', cols: ['scheme', 'C_A/C_A0', 'X_A'], rows: [['A: mixed then plug, microfluid', 1 - XA, XA], ['B: plug then mixed, microfluid', 1 - XB, XB], ['C–E: latest mixing, or macrofluid', 1 - Xseg, Xseg]] });
        const diag = (() => {
          let s = '';
          const row = (y, first, lab) => { s += D.text(10, y + 4, lab, { anchor: 'start', fs: 10 }); s += D.pipe([[96, y], [112, y]], {}); if (first === 'M') { s += D.tank(114, y - 18, 34, 36, { fill: 0.4, stir: false }); s += D.pipe([[150, y], [162, y]], { arrow: false }); s += D.tube(164, y - 8, 60, 16, { prof: [0.4, 0.62] }); s += D.pipe([[226, y], [246, y]], {}); } else { s += D.tube(114, y - 8, 60, 16, { prof: [0, 0.5] }); s += D.pipe([[176, y], [188, y]], { arrow: false }); s += D.tank(190, y - 18, 34, 36, { fill: 0.63, stir: false }); s += D.pipe([[226, y], [246, y]], {}); } };
          row(30, 'M', 'A: early'); row(84, 'P', 'B: later');
          s += D.text(252, 34, N.fmt(1 - XA, 4), { anchor: 'start', fs: 10.5 }) + D.text(252, 88, N.fmt(1 - XB, 4), { anchor: 'start', fs: 10.5 });
          return D.svg(300, 110, s, 'mixed then plug, and plug then mixed');
        })();
        return { kpi, warn, info, plots, tables, diagram: diag };
      }
      const Is = String(v.Ilist).split(/[,;\s]+/).map(parseFloat).filter(x => x >= 0).slice(0, 10);
      let Np = Math.round(+v.Np);
      const budget = (v.what === 'twofeed' ? 1.6e6 : 4e6) / Math.max(1, Is.length);
      if (Np * 800 > budget) { const np2 = Math.max(200, Math.floor(budget / 800)); if (np2 < Np) { info.push('Simulated with ' + np2 + ' drops per run to keep the page responsive; scatter grows as 1/√drops.'); Np = np2; } }
      if (v.what === 'partial') {
        const n = +v.n, Rg = +v.Rg;
        const res = Is.map((I, j) => ({ I, r: curl({ I, Np, nsp: 1, seed: 1000 + j, feed: () => [1], react: (C, n0, dt) => { const a = C[0]; for (let i = 0; i < n0; i++) a[i] = stepN(a[i], Rg, n, dt); } }) }));
        const micro = Math.abs(n - 1) < 1e-12 ? 1 / (1 + Rg) : N.brent(c => 1 - c - Rg * Math.pow(Math.max(c, 0), n), 0, 1, { tol: 1e-14 });
        const macro = N.quad(th => stepN(1, Rg, n, th) * Math.exp(-th), 0, Infinity, { rtol: 1e-10 });
        kpi.push({ k: 'micro', l: 'Microfluid limit', v: micro }, { k: 'macro', l: 'Macrofluid limit', v: macro });
        res.forEach(q => kpi.push({ k: 'C_I' + q.I, l: 'C_A/C_A0 at I = ' + N.fmt(q.I, 3), v: q.r.mean[0] }));
        plots.push({ title: 'Coalescence–redispersion: from macro- to microfluid', x: { label: 'I = collision frequency × τ', log: true, min: Math.max(0.05, Math.min(...Is.filter(x => x > 0)) / 2) }, y: { label: 'C_A/C_A0' }, series: [{ name: 'simulation', x: res.map(q => Math.max(q.I, 0.05)), y: res.map(q => q.r.mean[0]), color: 1, marker: 'o' }, { name: 'microfluid', x: [0.05, 1000], y: [micro, micro], color: 2, dash: '5 4' }, { name: 'macrofluid', x: [0.05, 1000], y: [macro, macro], color: 3, dash: '2 3' }], cap: 'Points at I = 0 are drawn at the left edge.' });
        tables.push({ title: 'Simulation', cols: ['I', 'C_A/C_A0', '± (statistical)', 'drops sampled'], rows: res.map(q => [q.I, q.r.mean[0], q.r.se[0], q.r.n]) });
        return { kpi, warn, info, plots, tables };
      }
      /* two feeds: A in half the drops, B in the other half */
      const Da = +v.k1f, kr = +v.k2f, MB = +v.MBf;
      /* concentrations scaled by the B feed concentration; A drops carry 2·C_A0, B drops 2·C_B0 */
      const CA0 = 2 / MB, CB0 = 2;
      const react = (C, n0, dt) => {
        const A = C[0], B = C[1], Rr = C[2], S = C[3];
        const nsub = Math.max(2, Math.ceil(Da * Math.max(CA0, CB0) * dt * 2)), h = dt / nsub;
        for (let i = 0; i < n0; i++) {
          let a = A[i], b = B[i], r = Rr[i], s2 = S[i];
          if (b <= 0 || (a <= 0 && r <= 0)) continue;
          for (let k = 0; k < nsub; k++) {
            const r1 = Da * a * b, r2 = Da * kr * r * b;
            const am = a - 0.5 * h * r1, bm = b - 0.5 * h * (r1 + r2), rm = r + 0.5 * h * (r1 - r2);
            const q1 = Da * Math.max(am, 0) * Math.max(bm, 0), q2 = Da * kr * Math.max(rm, 0) * Math.max(bm, 0);
            a = Math.max(0, a - h * q1); b = Math.max(0, b - h * (q1 + q2)); r = Math.max(0, r + h * (q1 - q2)); s2 += h * q2;
          }
          A[i] = a; B[i] = b; Rr[i] = r; S[i] = s2;
        }
      };
      const res = Is.map((I, j) => ({ I, r: curl({ I, Np, nsp: 4, seed: 2000 + j, steps: 80, feed: k => (k === 0 ? [CA0, 0, 0, 0] : [0, CB0, 0, 0]), react }) }));
      /* micro limit: premixed feed in a mixed reactor */
      const fm = c => { const [a, b, r] = c; const r1 = Da * a * b, r2 = Da * kr * r * b; return [CA0 / 2 - a - r1, CB0 / 2 - b - r1 - r2, -r - r2 + r1]; };
      const sol = N.newton(fm, [CA0 / 4, CB0 / 4, CA0 / 8], { lo: [0, 0, 0], tol: 1e-13 });
      const Rmicro = sol.x[2] / (CA0 / 2);
      kpi.push({ k: 'Rmicro', l: 'R per A fed, perfectly mixed feed', v: Rmicro, main: true });
      res.forEach(q => kpi.push({ k: 'R_I' + q.I, l: 'R per A fed at I = ' + N.fmt(q.I, 3), v: q.r.mean[2] / (CA0 / 2) }));
      plots.push({ title: 'Yield of R against the intensity of mixing', x: { label: 'I = collision frequency × τ', log: true }, y: { label: 'C_R per A fed', min: 0 }, series: [{ name: 'R (simulation)', x: res.map(q => Math.max(q.I, 0.05)), y: res.map(q => q.r.mean[2] / (CA0 / 2)), color: 1, marker: 'o' }, { name: 'S (simulation)', x: res.map(q => Math.max(q.I, 0.05)), y: res.map(q => q.r.mean[3] / (CA0 / 2)), color: 2, marker: 's' }, { name: 'R, premixed limit', x: [Math.max(0.05, Math.min(...Is)), Math.max(...Is)], y: [Rmicro, Rmicro], color: 3, dash: '5 4' }], cap: 'At low intensity little A meets B at all; at intermediate intensity reaction runs in thin zones where R meets fresh B; vigorous mixing approaches the premixed result.' });
      tables.push({ title: 'Simulation (per mole of A fed)', cols: ['I', 'A left', 'R', 'S'], rows: res.map(q => [q.I, q.r.mean[0] / (CA0 / 2), q.r.mean[2] / (CA0 / 2), q.r.mean[3] / (CA0 / 2)]) });
      info.push('To keep R: mix as hard as possible, disperse B into A rather than A into B, or slow the reaction (Section 16.2).');
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
