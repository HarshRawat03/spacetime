/* Chapter 6 — Design for single reactions */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, K = ST.K, KIT = ST.KIT, D = ST.D;

  /* first/second order N equal tanks, dimensionless rate group K (kτ or kC_A0τ, total over all tanks) */
  function tanksX(order, Ntot, Ktot) {
    const k = Ktot / Ntot;
    let c = 1;
    for (let i = 0; i < Ntot; i++) c = order === 1 ? c / (1 + k) : (-1 + Math.sqrt(1 + 4 * k * c)) / (2 * k);
    return 1 - c;
  }
  function tanksK(order, Ntot, X) {
    if (Ntot === Infinity) return order === 1 ? -Math.log(1 - X) : X / (1 - X);
    if (order === 1) return Ntot * (Math.pow(1 - X, -1 / Ntot) - 1);
    return N.brent(Kt => tanksX(2, Ntot, Kt) - X, 1e-12, 1e9);
  }
  /* nth order with expansion: tau_m/tau_p */
  function ratioMP(n, eps, X) {
    const rate = x => Math.pow((1 - x) / (1 + eps * x), n);
    const tm = X / rate(X), tp = N.quad(x => 1 / rate(x), 0, X);
    return tm / tp;
  }

  ST.tool({
    id: 'mixvsplug', tag: 'H-601', ch: 6, sec: '6.1–6.2', title: 'Mixed against plug flow: size ratios and tanks in series',
    lede: 'The charts of Figs. 6.1, 6.5 and 6.6 as live calculators: how much larger a mixed flow reactor must be, and how quickly N tanks in series close the gap to plug flow.',
    keys: 'size comparison mixed flow plug flow volume ratio nth order expansion tanks in series figure 6.1 6.5 6.6 performance chart k tau',
    inputs: [
      { group: 'Chart' },
      { k: 'mode', t: 'seg', label: 'Show', opts: [['nth', 'single MFR vs PFR (Fig. 6.1)'], ['tanks', 'N tanks vs PFR (Figs. 6.5–6.6)']], v: 'tanks' },
      { k: 'n', label: 'Reaction order n', v: 2, min: 0, max: 3, step: 0.05, if: v => v.mode === 'nth' },
      { k: 'eps', label: 'Expansion factor $\\varepsilon_A$', v: 0, min: -0.9, max: 3, step: 0.05, if: v => v.mode === 'nth' },
      { k: 'order', t: 'seg', label: 'Order', opts: [['1', 'first'], ['2', 'second']], v: '2', if: v => v.mode === 'tanks' },
      { k: 'N', label: 'Number of equal tanks N', v: 2, min: 1, max: 50, step: 1, if: v => v.mode === 'tanks' },
      { group: 'Your operating point' },
      { k: 'X', label: 'Conversion $X_A$', v: 0.9, min: 0.01, max: 0.999, step: 0.001 },
      { k: 'Kt', label: 'Total rate group (kτ or kC_{A0}τ)', v: 180, min: 0.01, max: 1e4, scale: 'log', if: v => v.mode === 'tanks', hint: 'For the whole system: first order kτ, second order kC_{A0}τ.' }
    ],
    presets: [
      { name: 'Example 6.2', v: { mode: 'tanks', order: '2', N: 2, X: 0.9, Kt: 180 }, book: { XN: 0.974, KN: 27.5, gain: { v: 6.6, tol: 0.02 } }, note: '90% in one mixed reactor (kC_{A0}τ = 90). A second identical tank in series: at the same feed rate, or at the same 90% conversion.' },
      { name: 'Problem 6.2', v: { mode: 'tanks', order: '1', N: 2, X: 6 / 7, Kt: 6 }, note: 'Radioactive decay in a holdup tank leaves 1/7 of the activity; a baffle makes two tanks of half the size each.' },
      { name: 'Problem 6.13', v: { mode: 'tanks', order: '1', N: 2, X: 0.5, Kt: 2 * 40000 / 100 * Math.log(2) / 20 }, note: 'Half-life 20 h, 100 L/h through two 40 000 L tanks: total kτ = 800 h × ln2/20 h.' },
      { name: 'Fig. 6.1, n = 2', v: { mode: 'nth', n: 2, eps: 0, X: 0.9 } },
      { name: 'Fig. 6.1, gas, ε = 1', v: { mode: 'nth', n: 1, eps: 1, X: 0.9 } }
    ],
    eqs: [
      ['\\frac{(\\tau C_{A0}^{n-1})_m}{(\\tau C_{A0}^{n-1})_p} = \\frac{\\left[X_A\\frac{(1+\\varepsilon_AX_A)^n}{(1-X_A)^n}\\right]_m}{\\left[\\int_0^{X_A}\\frac{(1+\\varepsilon_AX_A)^n}{(1-X_A)^n}dX_A\\right]_p}', 'Eq. 6.1'],
      ['\\frac{C_N}{C_0} = \\frac{1}{(1+k\\tau_i)^N},\\qquad N\\tau_i = \\tau_{N}', 'Eq. 6.6'],
      ['C_i = \\frac{-1+\\sqrt{1+4k\\tau_iC_{i-1}}}{2k\\tau_i}', 'Eq. 6.8', 'Second order, applied tank by tank.']
    ],
    run(v) {
      const X = +v.X, kpi = [], plots = [];
      if (v.mode === 'nth') {
        const n = +v.n, e = +v.eps;
        const r = ratioMP(n, e, X);
        kpi.push({ k: 'ratio', l: 'V_mixed/V_plug', v: r, main: true, s: 'same feed, same conversion' });
        const xs = N.logspace(0.003, 0.99, 120); /* 1 - X */
        const series = [0, 0.5, 1, 2, 3].map((nn, i) => ({ name: 'n = ' + nn, x: xs, y: xs.map(q => ratioMP(nn, e, 1 - q)), color: i + 1 }));
        if ([0, 0.5, 1, 2, 3].indexOf(n) < 0) series.push({ name: 'n = ' + N.fmt(n, 3), x: xs, y: xs.map(q => ratioMP(n, e, 1 - q)), color: 6, w: 2.6 });
        plots.push({ title: 'Fig. 6.1 for ε_A = ' + N.fmt(e, 3), wide: true, x: { label: '1 − X_A', log: true, min: 0.003, max: 1 }, y: { label: 'τ_m/τ_p', log: true, min: 1, max: 1000 }, series, marks: [{ t: 'pt', x: 1 - X, y: r, color: 2, label: N.fmt(r, 3) }] });
      } else {
        const ord = +v.order, Nn = Math.round(+v.N), Kt = +v.Kt;
        const XN = tanksX(ord, Nn, Kt), KN = tanksK(ord, Nn, X), K1 = tanksK(ord, 1, X), Kp = tanksK(ord, Infinity, X);
        kpi.push({ k: 'XN', l: 'X_A from ' + Nn + ' tank' + (Nn > 1 ? 's' : '') + ' at the given rate group', v: XN, main: true });
        kpi.push({ k: 'KN', l: 'Rate group for X_A with ' + Nn + ' tanks', v: KN }, { k: 'K1', l: 'with one tank', v: K1 }, { k: 'Kp', l: 'with plug flow', v: Kp }, { k: 'ratio', l: 'V_N/V_plug at X_A', v: KN / Kp });
        if (Nn > 1) kpi.push({ k: 'gain', l: 'Throughput gain, N units each the size of one', v: Nn * K1 / KN, s: 'same conversion' });
        const xs = N.logspace(0.002, 1, 120);
        const Ns = [1, 2, 3, 4, 5, 10, 20, Infinity];
        const series = Ns.map((nn, i) => ({ name: nn === Infinity ? 'plug' : 'N = ' + nn, x: xs, y: xs.map(q => tanksK(ord, nn, 1 - q) / tanksK(ord, Infinity, 1 - q)), color: i + 1, w: nn === Nn ? 2.8 : 1.4 }));
        const Ks = ord === 1 ? [0.5, 1, 2, 5, 10, 20, 50] : [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
        const marks = [{ t: 'pt', x: 1 - X, y: KN / Kp, color: 2, label: 'X_A, N = ' + Nn }, { t: 'pt', x: 1 - XN, y: tanksK(ord, Nn, XN) / tanksK(ord, Infinity, XN), color: 6, label: 'given ' + (ord === 1 ? 'kτ' : 'kC_{A0}τ') }];
        Ks.forEach(K => {
          const pts = N.linspace(0, 1, 40).map(s => 1 + s * 39).filter(nn => true);
          const xx = [], yy = [];
          [1, 1.2, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 40, 100, 1000].forEach(nn => { const x1 = tanksX(ord, nn, K); xx.push(1 - x1); yy.push(K / tanksK(ord, Infinity, x1)); });
          series.push({ x: xx, y: yy, color: 8, dash: '3 4', w: 1, legend: false, hover: false, noRange: true });
          if (pts.length) marks.push({ t: 'text', x: xx[0], y: yy[0] * 1.06, s: N.fmt(K, 3), color: 8, anchor: 'middle', noRange: true });
        });
        plots.push({ title: (ord === 1 ? 'Fig. 6.5, first order' : 'Fig. 6.6, second order') + ': dashed lines of constant ' + (ord === 1 ? 'kτ' : 'kC_{A0}τ'), wide: true, h: 420, x: { label: '1 − X_A', log: true, min: 0.002, max: 1 }, y: { label: 'V_N/V_p', log: true, min: 1, max: 100 }, series, marks });
        const rows = [1, 2, 3, 4, 5, 10].map(nn => [nn, tanksK(ord, nn, X), tanksK(ord, nn, X) / Kp, tanksX(ord, nn, Kt)]);
        return { kpi, plots, tables: [{ title: 'Effect of staging at X_A = ' + N.fmt(X, 3), cols: ['N', 'rate group needed', 'V_N/V_p', 'X_A at given rate group'], rows }], diagram: D.train(Array.from({ length: Math.min(Nn, 8) }, (_, i) => ({ type: 'mfr', tag: 'R-60' + (i + 1), X: 1 - Math.pow(1 - XN, (i + 1) / Nn) })), { h: 118 }).svg };
      }
      return { kpi, plots };
    }
  });

  /* ---------------- reactor train with permutations and parallel branches ---------------- */
  function parseUnits(s) {
    const out = [];
    String(s).split(/[,;\n]+/).map(t => t.trim()).filter(Boolean).forEach(t => {
      const m = /^(p|pfr|plug|m|mfr|mixed|cstr)\s*([0-9.eE+-]+)$/i.exec(t.replace(/\s+/g, ' '));
      if (!m) throw new Error('Cannot read “' + t + '”. Write units like “PFR 2, MFR 1.5”.');
      out.push({ type: /^p/i.test(m[1]) ? 'pfr' : 'mfr', tau: parseFloat(m[2]) });
    });
    if (!out.length) throw new Error('No units given.');
    return out;
  }
  function runTrain(m, units, X0 = 0) {
    let X = X0;
    const xs = [X0];
    for (const u of units) {
      X = u.type === 'mfr' ? K.XMFR(m, u.tau, X) : K.XPFR(m, u.tau, X).X;
      xs.push(X);
    }
    return xs;
  }
  function perms(arr) {
    const res = [], seen = new Set();
    const go = (a, rest) => {
      if (!rest.length) { const key = a.map(u => u.type + u.tau).join('|'); if (!seen.has(key)) { seen.add(key); res.push(a); } return; }
      rest.forEach((u, i) => go(a.concat([u]), rest.slice(0, i).concat(rest.slice(i + 1))));
    };
    go([], arr);
    return res;
  }
  const label = units => units.map(u => (u.type === 'pfr' ? 'P' : 'M') + N.fmt(u.tau, 3)).join(' → ');

  ST.tool({
    id: 'train', tag: 'H-602', ch: 6, sec: '6.2', title: 'Reactor trains, every ordering, parallel branches',
    lede: 'Any sequence of plug and mixed flow units for any single rate law: conversion after each unit, every distinct ordering of the same units ranked, and the best split of feed between parallel branches.',
    keys: 'reactors in series different types best arrangement order of reactors permutations parallel branches feed split jones graphical method figure 6.8 6.12 rules',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'Arrangement' },
      { k: 'mode', t: 'seg', label: 'Layout', opts: [['series', 'one train'], ['par', 'parallel branches']], v: 'series' },
      { k: 'units', t: 'text', label: 'Units in order (space-times)', v: 'MFR 1, PFR 1, MFR 3', mono: true, if: v => v.mode === 'series', hint: 'Space-time of each unit based on the feed to the train, e.g. “PFR 2, MFR 0.5”.' },
      { k: 'perm', t: 'check', label: 'Rank every distinct ordering of these units', v: true, if: v => v.mode === 'series' },
      { k: 'br1', t: 'text', label: 'Branch 1 (volumes)', v: 'PFR 50, PFR 30', mono: true, if: v => v.mode === 'par' },
      { k: 'br2', t: 'text', label: 'Branch 2 (volumes)', v: 'PFR 40', mono: true, if: v => v.mode === 'par' },
      { k: 'br3', t: 'text', label: 'Branch 3 (optional)', v: '', mono: true, if: v => v.mode === 'par' },
      { k: 'vt', label: 'Total feed rate v₀', unit: v => KIT.volUnit(v._cu) + '/' + v._tu, v: 40, min: 1e-6, max: 1e9, scale: 'log', if: v => v.mode === 'par' }
    ].concat(KIT.kin({ d: { law: 'power', a: 2, k: 1, CA0: 1 } })),
    presets: [
      { name: 'Example 6.1', v: { mode: 'par', br1: 'PFR 50, PFR 30', br2: 'PFR 40', br3: '', vt: 40, law: 'power', a: 1, ob: 0, k: 0.05, CA0: 1, eps: 0, _cu: 'mol/L', _tu: 'min' }, book: { f1: 0.6667 }, note: 'Branch D (50 L + 30 L) beside branch E (40 L): what fraction of the feed should go to D? Any kinetics give the same answer for plug flow branches.' },
      { name: 'Fig. 6.12, n = 2', v: { mode: 'series', units: 'MFR 1, PFR 1, MFR 3', perm: true, law: 'power', a: 2, ob: 0, k: 1, CA0: 1, eps: 0 }, note: 'For n > 1 the best order is plug, small mixed, large mixed.' },
      { name: 'Fig. 6.12, n = 0.5', v: { mode: 'series', units: 'MFR 1, PFR 1, MFR 3', perm: true, law: 'power', a: 0.5, ob: 0, k: 0.4, CA0: 1, eps: 0 }, note: 'For n < 1 the order reverses.' },
      { name: 'Problem 6.3', v: { mode: 'series', units: 'MFR 3, PFR 9', perm: true, law: 'power', a: 2, ob: 0, k: 1, CA0: 4, eps: 0 }, note: 'Second order from 4 mol/L: the mixed unit brings A to 1 mol/L (τ = 3 with k = 1) and the plug unit is three times larger. The exit comes out at 0.1 mol/L; the ranking shows the reverse order is worse.' },
      { name: 'Problem 6.11', v: { mode: 'series', units: 'MFR 66.67, MFR 66.67, MFR 66.67, MFR 66.67', perm: false, law: 'auto', k: 0.001, CA0: 10, MR: 0.0001, r: 1, eps: 0, _tu: 's' }, note: 'Autocatalytic A → R with four 100 L tanks and 1.5 L/s; try putting tanks in parallel too.' }
    ],
    eqs: [
      ['\\frac{V_i}{F_{A0}} = \\frac{X_i-X_{i-1}}{(-r_A)_i}\\;\\text{(mixed)},\\qquad \\frac{V_i}{F_{A0}} = \\int_{X_{i-1}}^{X_i}\\frac{dX_A}{-r_A}\\;\\text{(plug)}', 'Eqs. 6.9–6.14'],
      ['\\left(\\frac{V}{F}\\right)_{\\text{branch 1}} = \\left(\\frac{V}{F}\\right)_{\\text{branch 2}} = \\cdots', '§6.2', 'Plug flow branches: split the feed so streams that meet have the same composition.']
    ],
    method: ['Each unit is solved in turn on the conversion basis of the original feed, so the same equations hold for gases with expansion. With several units, every distinct ordering is generated (duplicates removed) and solved; this is the general check behind the book’s rules for arranging a set of reactors.', 'For parallel branches the split fractions are optimized numerically to maximize the conversion of the mixed outlet.'],
    run(v) {
      const m = KIT.model(v), warn = [], kpi = [], plots = [], tables = [];
      if (v.mode === 'par') {
        const brs = [v.br1, v.br2, v.br3].filter(s => String(s || '').trim()).map(parseUnits);
        if (brs.length < 2) throw new Error('Give at least two branches.');
        const vt = +v.vt;
        const outX = f => {
          let tot = 0;
          brs.forEach((b, i) => { if (f[i] <= 1e-9) return; const units = b.map(u => ({ type: u.type, tau: u.tau / (f[i] * vt) })); tot += f[i] * runTrain(m, units)[units.length]; });
          return tot;
        };
        const nb = brs.length;
        const toF = z => { const e = z.map(Math.exp), s = 1 + N.sum(e); return [1 / s].concat(e.map(x => x / s)); };
        let best;
        if (nb === 2) { const r = N.fminScan(f => -outX([f, 1 - f]), 1e-4, 1 - 1e-4, 60); best = [r.x, 1 - r.x]; }
        else { const r = N.nelderMead(z => -outX(toF(z)), new Array(nb - 1).fill(0)); best = toF(r.x); }
        const Vs = brs.map(b => N.sum(b.map(u => u.tau)));
        const equal = Vs.map(V => V / N.sum(Vs));
        best.forEach((f, i) => kpi.push({ k: 'f' + (i + 1), l: 'Feed fraction to branch ' + (i + 1), v: f, main: i === 0 }));
        kpi.push({ k: 'X', l: 'Mixed outlet X_A', v: outX(best) }, { k: 'Xeq', l: 'X_A with split ∝ volume', v: outX(equal) });
        const allPlug = brs.every(b => b.every(u => u.type === 'pfr'));
        if (allPlug) warn.push('All branches are plug flow, so the best split puts the same V/F on every branch: fractions ' + equal.map(f => N.fmt(f, 4)).join(', ') + '.');
        tables.push({ title: 'Branches at the best split', cols: ['branch', 'units', 'volume', 'feed fraction', 'V/F', 'X_A leaving'], rows: brs.map((b, i) => { const units = b.map(u => ({ type: u.type, tau: u.tau / (best[i] * vt) })); return [i + 1, label(b), Vs[i], best[i], Vs[i] / (best[i] * vt), runTrain(m, units)[units.length]]; }) });
        const fs = N.linspace(0.02, 0.98, 60);
        if (nb === 2) plots.push({ title: 'Outlet conversion against the split', x: { label: 'fraction to branch 1', min: 0, max: 1 }, y: { label: 'X_A' }, series: [{ name: 'X_A', x: fs, y: fs.map(f => outX([f, 1 - f])), color: 1 }], marks: [{ t: 'v', x: best[0], label: 'best', color: 2 }] });
        let s = D.pipe([[4, 60], [40, 60]], { label: 'feed' }) + D.node(42, 60);
        brs.forEach((b, i) => {
          const y = 30 + i * 62;
          s += D.pipe([[42, 60], [42, y], [62, y]], { flow: true });
          let x = 64;
          const units = b.map(u => ({ type: u.type, tau: u.tau / (best[i] * vt) }));
          const xs = runTrain(m, units);
          b.forEach((u, j) => {
            if (u.type === 'pfr') { s += D.tube(x, y - 9, 80, 18, { prof: [xs[j], xs[j + 1]] }); x += 80; }
            else { s += D.tank(x, y - 22, 40, 44, { fill: xs[j + 1], stir: false }); x += 40; }
            s += D.pipe([[x, y], [x + 22, y]]); x += 24;
          });
          s += D.text(x + 40, y + 4, 'f = ' + N.fmt(best[i], 3), { cls: 'pid-r' });
        });
        return { kpi, warn, plots, tables, diagram: D.svg(470, 40 + nb * 62, s, 'parallel branches') };
      }
      const units = parseUnits(v.units);
      const xs = runTrain(m, units);
      const Xf = xs[xs.length - 1];
      kpi.push({ k: 'X', l: 'Final X_A', v: Xf, main: true }, { k: 'CA', l: 'Final C_A', v: m.conc(Xf).A, u: v._cu }, { k: 'tau', l: 'Total space-time', v: N.sum(units.map(u => u.tau)), u: v._tu });
      /* Levenspiel plot with each unit's area */
      const lev = KIT.levPlot(m, Xf, { tu: v._tu, markX: true, title: 'Each unit on the Levenspiel plot (Fig. 6.12)' });
      units.forEach((u, i) => {
        if (u.type === 'mfr') lev.marks.push({ t: 'rect', x0: xs[i], x1: xs[i + 1], y0: 0, y1: m.CA0 / m.rate(xs[i + 1]), color: 2 + (i % 4), fill: 0.22, noRange: true, label: 'M' + (i + 1) });
        else { const sx = N.linspace(xs[i], xs[i + 1], 60); lev.series.push({ x: sx, y: sx.map(x => m.CA0 / m.rate(x)), type: 'area', fill: 0.22, color: 2 + (i % 4), hover: false, legend: false, noRange: true }); }
      });
      plots.push(lev);
      /* Jones construction for constant density */
      if (!m.eps) {
        const Cs = N.linspace(0, m.CA0, 200);
        const rC = c => m.rate(1 - c / m.CA0);
        const jones = { title: 'Rate against concentration (Fig. 6.8): slope −1/τ_i per mixed unit', x: { label: 'C_A', unit: v._cu, min: 0, max: m.CA0 }, y: { label: '−r_A', unit: v._cu + '/' + v._tu, min: 0 }, series: [{ name: '−r_A', x: Cs, y: Cs.map(rC), color: 1 }], marks: [] };
        units.forEach((u, i) => {
          const c0 = m.CA0 * (1 - xs[i]), c1 = m.CA0 * (1 - xs[i + 1]);
          if (u.type === 'mfr') jones.marks.push({ t: 'seg', x0: c0, y0: 0, x1: c1, y1: rC(c1), color: 2, w: 1.6 }, { t: 'pt', x: c1, y: rC(c1), color: 2, label: 'M' + (i + 1) });
          else jones.marks.push({ t: 'pt', x: c1, y: rC(c1), color: 3, label: 'P' + (i + 1) });
        });
        plots.push(jones);
      }
      tables.push({ title: 'Unit by unit', cols: ['unit', 'type', 'τ', 'X_A in', 'X_A out', 'C_A out'], rows: units.map((u, i) => [i + 1, u.type === 'pfr' ? 'plug' : 'mixed', u.tau, xs[i], xs[i + 1], m.conc(xs[i + 1]).A]) });
      if (v.perm && units.length > 1) {
        if (units.length > 7) warn.push('Ordering comparison is limited to seven units.');
        else {
          const all = perms(units).map(p => ({ p, X: runTrain(m, p)[p.length] })).sort((a, b) => b.X - a.X);
          const yours = all.findIndex(o => label(o.p) === label(units));
          kpi.push({ k: 'bestX', l: 'Best ordering gives', v: all[0].X, s: label(all[0].p) }, { k: 'rank', l: 'Your ordering ranks', v: (yours + 1) + ' of ' + all.length });
          tables.push({ title: 'Every distinct ordering, best first', cols: ['rank', 'ordering', 'final X_A', 'C_A out'], rows: all.slice(0, 40).map((o, i) => [i + 1, label(o.p), o.X, m.conc(o.X).A]) });
          if (all.length > 40) warn.push('Showing the best 40 of ' + all.length + ' orderings.');
        }
      }
      const diag = D.train(units.map((u, i) => ({ type: u.type, X: xs[i + 1], tag: 'R-6' + (i + 1).toString().padStart(2, '0'), prof: u.type === 'pfr' ? N.linspace(0, 1, 12).map(s2 => K.XPFR(m, u.tau * s2, xs[i]).X) : undefined })), { h: 120 }).svg;
      return { kpi, warn, plots, tables, diagram: diag };
    }
  });

  /* ---------------- optimum sizes of mixed reactors in series ---------------- */
  ST.tool({
    id: 'mfropt', tag: 'H-603', ch: 6, sec: '6.2', title: 'Best sizes for mixed flow reactors in series',
    lede: 'The smallest total volume of N mixed reactors in series for a target conversion, found by optimizing the intermediate conversions, compared with equal-size units. For two tanks this is the maximization of rectangles.',
    keys: 'maximization of rectangles optimum intermediate conversion unequal tanks two mixed flow reactors smallest first larger first szepe levenspiel',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'Duty' },
      { k: 'Xf', label: 'Final conversion $X_{Af}$', v: 0.9, min: 0.05, max: 0.999, step: 0.001 },
      { k: 'Nt', label: 'Number of tanks N', v: 2, min: 1, max: 8, step: 1 }
    ].concat(KIT.kin({ d: { law: 'power', a: 2, k: 1, CA0: 1 } }), KIT.flowInputs()),
    presets: [
      { name: 'Problem 6.7', v: { Xf: 0.99, Nt: 2, law: 'sat', sk1: 5, sk2: 5, CA0: 1, basis: 'v0', v0: 10, _cu: 'mol/L', _tu: 'min' }, note: '−r_A = C_A/(0.2 + C_A), 10 L/min to 99% conversion with two mixed reactors.' },
      { name: 'n = 2 (small tank first)', v: { Xf: 0.9, Nt: 2, law: 'power', a: 2, k: 1, CA0: 1 } },
      { name: 'n = 0.5 (large tank first)', v: { Xf: 0.9, Nt: 2, law: 'power', a: 0.5, k: 1, CA0: 1 } }
    ],
    eqs: [
      ['\\frac{\\tau_i}{C_{A0}} = \\frac{X_i-X_{i-1}}{(-r_A)_i}', 'Eq. 5.12'],
      ['\\left(\\frac{dy}{dx}\\right)_M = \\frac{y_M}{x_M}', 'Fig. 6.10', 'Two tanks: the best intermediate point makes the slope of the 1/(−r_A) curve equal to the diagonal of the rectangle.']
    ],
    run(v) {
      const m = KIT.model(v), Xf = Math.min(+v.Xf, Math.min(m.Xe, m.Xmax) * 0.999), Nt = Math.round(+v.Nt);
      const tauOf = xs => { let t = 0, prev = 0; xs.forEach(x => { t += K.tauMFR(m, x, prev); prev = x; }); return t + K.tauMFR(m, Xf, prev); };
      let inter = [];
      if (Nt === 2) inter = [N.fminScan(x => tauOf([x]), 1e-5, Xf * (1 - 1e-6), 80).x];
      else if (Nt > 2) {
        const toX = z => { const e = z.map(Math.exp), s = N.sum(e) + 1; let c = 0; return e.map(x => (c += x / s) * Xf); };
        const r = N.nelderMead(z => tauOf(toX(z)), new Array(Nt - 1).fill(0), { maxIt: 4000 });
        inter = toX(r.x);
      }
      const xs = [0].concat(inter, [Xf]);
      const taus = xs.slice(1).map((x, i) => K.tauMFR(m, x, xs[i]));
      const tauMin = N.sum(taus);
      /* equal tanks */
      const eqX = t => { let X = 0; for (let i = 0; i < Nt; i++) X = K.XMFR(m, t, X); return X; };
      let hi = K.tauMFR(m, Xf);
      const tEq = N.brent(t => eqX(t) - Xf, 1e-12, hi * 1.0001 + 1e-9);
      const kpi = [{ k: 'tau', l: 'Smallest total τ', v: tauMin, u: v._tu, main: true }, { k: 'taueq', l: 'Total τ with equal tanks', v: tEq * Nt, u: v._tu }, { k: 'save', l: 'Saving over equal tanks', v: (1 - tauMin / (tEq * Nt)) * 100, u: '%' }, { k: 'tau1', l: 'One tank', v: K.tauMFR(m, Xf), u: v._tu }, { k: 'taup', l: 'Plug flow', v: K.tauPFR(m, Xf), u: v._tu }];
      const v0 = KIT.v0of(v, m.CA0);
      const rows = taus.map((t, i) => [i + 1, xs[i], xs[i + 1], t].concat(isFinite(v0) ? [t * v0] : []));
      const lev = KIT.levPlot(m, Xf, { tu: v._tu, title: 'Optimal rectangles' });
      taus.forEach((t, i) => lev.marks.push({ t: 'rect', x0: xs[i], x1: xs[i + 1], y0: 0, y1: m.CA0 / m.rate(xs[i + 1]), color: 2 + (i % 4), fill: 0.22, noRange: true, label: String(i + 1) }));
      if (Nt === 2) {
        const x1 = xs[1], y1 = m.CA0 / m.rate(x1), yf = m.CA0 / m.rate(Xf);
        lev.marks.push({ t: 'seg', x0: 0, y0: yf, x1: x1, y1: y1, color: 8, dash: '4 3' });
      }
      const warn = [];
      if (Nt === 2 && v.law === 'power' && +v.ob === 0) warn.push(+v.a > 1 ? 'Order above one: the smaller tank goes first.' : +v.a < 1 ? 'Order below one: the larger tank goes first.' : 'First order: equal tanks are best.');
      return { kpi, warn, plots: [lev], tables: [{ title: 'Optimal train', cols: ['tank', 'X_A in', 'X_A out', 'τ'].concat(isFinite(v0) ? ['V'] : []), rows }], diagram: D.train(xs.slice(1).map((x, i) => ({ type: 'mfr', X: x, tag: 'R-' + (i + 1) })), { h: 118 }).svg };
    }
  });

  /* ---------------- recycle reactor ---------------- */
  ST.tool({
    id: 'recycle', tag: 'H-604', ch: 6, sec: '6.3', title: 'Recycle reactor',
    lede: 'Plug flow with recycle for any rate law and expansion factor: space-time at a given recycle ratio, the optimum recycle ratio, the approach to mixed flow, and the recycle ratio equivalent to N tanks.',
    keys: 'recycle ratio recycle reactor optimum recycle autocatalytic backmixing equivalent tanks in series table',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'Duty' },
      { k: 'Xf', label: 'Conversion $X_{Af}$', v: 0.9, min: 0.01, max: 0.999, step: 0.001 },
      { k: 'R', label: 'Recycle ratio R', v: 2, min: 0, max: 100, scale: 'log' }
    ].concat(KIT.kin({ d: { law: 'power', a: 1, k: 1, CA0: 1 } }), KIT.flowInputs()),
    presets: [
      { name: 'Problem 6.9', v: { Xf: 0.9, R: 2, law: 'power', a: 1, ob: 0, k: 1, CA0: 10, eps: 0, basis: 'none' }, note: 'First order at 90% with R = 2. Shutting off the recycle lets the same reactor treat τ(R=2)/τ(R=0) times the feed.' },
      { name: 'Problem 6.10', v: { Xf: 0.96, R: 1, law: 'auto', k: 1, CA0: 2, MR: 0, r: 1, eps: 0, basis: 'none', _cu: 'mol/L', _tu: 'min' }, note: 'A → R with −r_A = C_AC_R, 2 mol/L of A and no R in the feed, 96% conversion in a 10 L tube: recycle is essential.' },
      { name: 'Second order', v: { Xf: 0.9, R: 1, law: 'power', a: 2, ob: 0, k: 1, CA0: 1, eps: 0 } }
    ],
    eqs: [
      ['\\frac{V}{F_{A0}} = (R+1)\\int_{\\frac{R}{R+1}X_{Af}}^{X_{Af}}\\frac{dX_A}{-r_A}', 'Eq. 6.21'],
      ['\\frac{k\\tau}{R+1} = \\ln\\left[\\frac{C_{A0}+RC_{Af}}{(R+1)C_{Af}}\\right]', 'Eq. 6.22', 'First order, ε_A = 0.'],
      ['\\frac{kC_{A0}\\tau}{R+1} = \\frac{C_{A0}(C_{A0}-C_{Af})}{C_{Af}(C_{A0}+RC_{Af})}', 'Eq. 6.23', 'Second order, ε_A = 0.'],
      ['\\frac{1}{-r_A}\\bigg|_{X_{Ai}} = \\frac{\\int_{X_{Ai}}^{X_{Af}}\\frac{dX_A}{-r_A}}{X_{Af}-X_{Ai}}', 'Eq. 6.28', 'Optimum: the entering 1/(−r_A) equals its average over the reactor.']
    ],
    run(v) {
      const m = KIT.model(v), Xf = Math.min(+v.Xf, Math.min(m.Xe, m.Xmax) * 0.999), R = +v.R;
      const tauR = Rr => { const X1 = Rr / (Rr + 1) * Xf; const t = (Rr + 1) * K.tauPFR(m, Xf, X1); return t; };
      const tR = tauR(R), tP = K.tauPFR(m, Xf), tM = K.tauMFR(m, Xf);
      const opt = N.fminScan(lr => { const t = tauR(Math.pow(10, lr)); return isFinite(t) ? t : 1e300; }, -4, 3, 90);
      let Ropt = Math.pow(10, opt.x), tOpt = opt.f;
      if (tP <= tOpt) { Ropt = 0; tOpt = tP; }
      const kpi = [{ k: 'tau', l: 'τ at R = ' + N.fmt(R, 3), v: tR, u: v._tu, main: true }, { k: 'Ropt', l: 'Optimum R', v: Ropt, s: Ropt === 0 ? 'no recycle: plug flow is best' : '' }, { k: 'tauopt', l: 'τ at optimum', v: tOpt, u: v._tu }, { k: 'taup', l: 'Plug flow (R = 0)', v: tP, u: v._tu }, { k: 'taum', l: 'Mixed flow (R → ∞)', v: tM, u: v._tu }];
      const v0 = KIT.v0of(v, m.CA0);
      if (isFinite(v0)) kpi.push({ k: 'V', l: 'Volume at R', v: tR * v0, u: KIT.volUnit(v._cu) }, { k: 'vr', l: 'Recycle flow', v: R * v0, u: KIT.volUnit(v._cu) + '/' + v._tu });
      if (v.law === 'power' && +v.ob === 0 && !m.eps) {
        const k = +v.k, Cf = m.CA0 * (1 - Xf);
        if (Math.abs(+v.a - 1) < 1e-12) kpi.push({ k: 'chk', l: 'Eq. 6.22 check', v: (R + 1) / k * Math.log((m.CA0 + R * Cf) / ((R + 1) * Cf)), u: v._tu });
        if (Math.abs(+v.a - 2) < 1e-12) kpi.push({ k: 'chk', l: 'Eq. 6.23 check', v: (R + 1) / (k * m.CA0) * m.CA0 * (m.CA0 - Cf) / (Cf * (m.CA0 + R * Cf)), u: v._tu });
      }
      if (m.rate(0) <= 0) kpi.push({ k: 'note', l: 'Feed rate is zero', v: 'needs recycle', state: 'warn' });
      const Rs = N.logspace(1e-3, 1e3, 120);
      const tauPlot = { title: 'Space-time against recycle ratio', x: { label: 'R', log: true }, y: { label: 'τ', unit: v._tu, log: true }, series: [{ name: 'τ(R)', x: Rs, y: Rs.map(tauR), color: 1 }], marks: [{ t: 'h', y: tP, label: 'plug', color: 3 }, { t: 'h', y: tM, label: 'mixed', color: 2 }, { t: 'pt', x: R, y: tR, color: 6, label: 'R' }] };
      if (Ropt > 0) tauPlot.marks.push({ t: 'pt', x: Ropt, y: tOpt, color: 2, label: 'optimum' });
      const X1 = R / (R + 1) * Xf;
      const lev = KIT.levPlot(m, Xf, { tu: v._tu, title: 'Recycle rectangle (Fig. 6.14)' });
      lev.marks.push({ t: 'rect', x0: 0, x1: Xf, y0: 0, y1: tR / Xf, color: 6, fill: 0.2, noRange: true }, { t: 'v', x: X1, label: 'X_{A1}', color: 6, noRange: true });
      const sx = N.linspace(X1, Xf, 60);
      lev.series.push({ x: sx, y: sx.map(x => m.CA0 / m.rate(x)), type: 'area', fill: 0.2, color: 3, hover: false, legend: false, noRange: true });
      /* equivalent R for N tanks, first and second order (book table in §6.3) */
      const eqR = (ord, Nn, X) => {
        const KN = ord === 1 ? Nn * (Math.pow(1 - X, -1 / Nn) - 1) : N.brent(Kt => { const kk = Kt / Nn; let c = 1; for (let i = 0; i < Nn; i++) c = (-1 + Math.sqrt(1 + 4 * kk * c)) / (2 * kk); return 1 - c - X; }, 1e-9, 1e8);
        const Cf = 1 - X;
        const Kr = Rr => (ord === 1 ? (Rr + 1) * Math.log((1 + Rr * Cf) / ((Rr + 1) * Cf)) : (Rr + 1) * (1 - Cf) / (Cf * (1 + Rr * Cf)));
        return N.brent(Rr => Kr(Rr) - KN, 1e-9, 1e6);
      };
      const book = { 1: { 2: [1.0, 2.2, 5.4], 3: [0.5, 1.1, 2.1], 4: [0.33, 0.68, 1.3], 10: [0.11, 0.22, 0.36] }, 2: { 2: [1.0, 2.8, 7.5], 3: [0.5, 1.4, 2.9], 4: [0.33, 0.90, 1.7], 10: [0.11, 0.29, 0.5] } };
      const rows = [];
      [2, 3, 4, 10].forEach(Nn => { const r1 = [0.5, 0.9, 0.99].map(X => eqR(1, Nn, X)), r2 = [0.5, 0.9, 0.99].map(X => eqR(2, Nn, X)); rows.push([Nn].concat(r1.map((r, i) => N.fmt(r, 3) + ' (' + book[1][Nn][i] + ')'), r2.map((r, i) => N.fmt(r, 3) + ' (' + book[2][Nn][i] + ')'))); });
      return {
        kpi, plots: [tauPlot, lev],
        tables: [{ title: 'Recycle ratio matching N equal tanks, computed (book in brackets)', cols: ['N', '1st order X=0.5', '0.90', '0.99', '2nd order X=0.5', '0.90', '0.99'], rows }],
        diagram: (() => {
          let s = D.pipe([[6, 46], [64, 46]], { label: 'feed' }) + D.node(68, 46) + D.pipe([[71, 46], [96, 46]]);
          s += D.tube(98, 35, 180, 22, { prof: N.linspace(0, 1, 14).map(q => X1 + q * (Xf - X1)), tag: 'R-604' });
          s += D.pipe([[280, 46], [330, 46]]) + D.node(334, 46) + D.pipe([[337, 46], [420, 46]]);
          s += D.pipe([[334, 50], [334, 98], [68, 98], [68, 52]]);
          s += D.text(200, 112, 'R = ' + N.fmt(R, 3), { cls: 'pid-r' }) + D.text(380, 38, 'X ' + N.fmt(Xf, 3), { cls: 'pid-r' });
          return D.svg(430, 124, s, 'recycle reactor');
        })()
      };
    }
  });

  /* ---------------- autocatalytic and design from rate data ---------------- */
  ST.tool({
    id: 'autocat', tag: 'H-605', ch: 6, sec: '6.4', title: 'Autocatalytic reactions and design from rate data',
    lede: 'For a U-shaped 1/(−r_A) curve, from A + R → 2R kinetics or straight from mixed flow test data: plug with recycle, one or two mixed tanks, mixed-then-plug, and the separation-and-recycle minimum, side by side.',
    keys: 'autocatalytic reaction U shaped curve best reactor setup mixed then plug recycle optimum example 6.3 enzyme waste water rate data design',
    units: { C: 'mmol/m³', t: 'min' }, slow: true,
    inputs: [
      { group: 'Kinetics' },
      { k: 'src', t: 'seg', label: 'Rate from', opts: [['auto', 'A + R → 2R'], ['data', 'mixed flow test data']], v: 'data' },
      { k: 'k', label: 'k', unit: v => '1/(' + v._cu + '·' + v._tu + ')', v: 0.001, min: 1e-8, max: 1e4, scale: 'log', if: v => v.src === 'auto' },
      { k: 'CR0', label: 'R in the feed, $C_{R0}$', unit: '[C]', v: 0.1, min: 0, max: 1e6, if: v => v.src === 'auto' },
      { k: 'runs', t: 'table', label: 'Mixed flow runs', cols: ['C_{A0}', 'C_A', 'τ'], v: [[2, 0.5, 30], [5, 3, 1], [6, 1, 50], [6, 2, 8], [11, 6, 4], [14, 10, 20], [16, 8, 20], [24, 4, 4]], if: v => v.src === 'data' },
      { group: 'Duty' },
      { k: 'CA0', label: 'Feed $C_{A0}$', unit: '[C]', v: 10, min: 1e-6, max: 1e6, scale: 'log' },
      { k: 'CAf', label: 'Outlet $C_{Af}$', unit: '[C]', v: 1, min: 1e-6, max: 1e6, scale: 'log' },
      { k: 'v0', label: 'Feed rate v₀ (for volumes)', unit: v => KIT.volUnit(v._cu) + '/' + v._tu, v: 0.1, min: 1e-9, max: 1e9, scale: 'log' },
      { k: 'show', t: 'select', label: 'Draw', opts: [['rec', 'plug flow with recycle'], ['m1', 'one mixed tank'], ['m2', 'two mixed tanks'], ['mp', 'mixed then plug'], ['sep', 'mixed at the minimum, with separation']], v: 'mp' }
    ],
    presets: [
      { name: 'Example 6.3', v: { src: 'data', runs: [[2, 0.5, 30], [5, 3, 1], [6, 1, 50], [6, 2, 8], [11, 6, 4], [14, 10, 20], [16, 8, 20], [24, 4, 4]], CA0: 10, CAf: 1, v0: 0.1, show: 'mp', _cu: 'mmol/m³', _tu: 'min' },
        book: { Vrec: { v: 1.08, tol: 0.2 }, Rrec: { v: 0.607, tol: 0.35 }, Vm1: { v: 9.0, tol: 0.02 }, Vm2: { v: 2.19, tol: 0.2 }, Vmp: { v: 0.70, tol: 0.2 } },
        note: 'Enzyme-degraded waste water: 0.1 m³/min from 10 to 1 mmol/m³. The book integrates under a curve drawn by eye through scattered points, so expect agreement within the reading accuracy of that curve.' },
      { name: 'Autocatalytic, Problem 6.11 kinetics', v: { src: 'auto', k: 0.001, CR0: 0, CA0: 10, CAf: 0.5, v0: 1.5, _cu: 'mol/L', _tu: 's', show: 'mp' } }
    ],
    eqs: [
      ['\\text{best without separation: mixed flow to the minimum of }1/(-r_A)\\text{, then plug flow}', 'Fig. 6.21a'],
      ['\\text{with separation and recycle of A: one mixed reactor at the minimum}', 'Fig. 6.21b'],
      ['\\frac{1}{-r_A} = \\frac{\\tau}{C_{A0}-C_A}', 'Eq. 5.13', 'Each mixed flow run gives one point on the curve.']
    ],
    method: ['From test data the curve is built by monotone (shape-preserving) interpolation of ln[1/(−r_A)] against C_A, which follows a U-shaped set of points without overshooting between them. Outside the data range the end values are held.', 'All schemes are then optimized numerically: the recycle ratio, the intermediate concentration of two tanks, and the switch point from mixed to plug flow.'],
    run(v) {
      const CA0 = +v.CA0, CAf = +v.CAf;
      if (!(CAf < CA0)) throw new Error('Outlet concentration must be below the feed.');
      let inv, pts = [];
      if (v.src === 'data') {
        pts = (v.runs || []).filter(r => r.slice(0, 3).every(isFinite) && r[0] > r[1]).map(r => ({ C: r[1], y: r[2] / (r[0] - r[1]) })).sort((a, b) => a.C - b.C);
        const merged = [];
        pts.forEach(p => { const q = merged[merged.length - 1]; if (q && Math.abs(q.C - p.C) < 1e-12) { q.y = Math.sqrt(q.y * p.y); } else merged.push({ C: p.C, y: p.y }); });
        if (merged.length < 3) throw new Error('Need at least three runs.');
        const f = N.pchip(merged.map(p => p.C), merged.map(p => Math.log(p.y)));
        const lo = merged[0].C, hi = merged[merged.length - 1].C;
        inv = C => Math.exp(f(N.clamp(C, lo, hi)));
        if (CAf < lo || CA0 > hi) var outside = true;
        pts = merged;
      } else {
        const k = +v.k, C0 = CA0 + (+v.CR0);
        inv = C => 1 / (k * C * (C0 - C));
      }
      const pfr = (a, b) => N.quad(inv, a, b, { rtol: 1e-9 });
      const mfr = (cin, cout) => (cin - cout) * inv(cout);
      const v0 = +v.v0;
      /* recycle */
      const tauRec = Rr => { const Cin = (CA0 + Rr * CAf) / (Rr + 1); return (Rr + 1) * pfr(CAf, Cin); };
      const ro = N.fminScan(lr => tauRec(Math.pow(10, lr)), -4, 3, 100);
      const Rrec = Math.pow(10, ro.x), tRec = ro.f;
      /* one, two tanks */
      const t1 = mfr(CA0, CAf);
      const two = N.fminScan(c => mfr(CA0, c) + mfr(c, CAf), CAf * 1.0001, CA0 * 0.9999, 120);
      const C1 = two.x, t2 = two.f;
      /* mixed then plug */
      const mp = N.fminScan(c => mfr(CA0, c) + pfr(CAf, c), CAf * 1.0001, CA0 * 0.9999, 120);
      const Cm = mp.x;
      /* minimum of 1/-r over [CAf, CA0] */
      const mn = N.fminScan(inv, CAf, CA0, 200);
      const tPlug = pfr(CAf, CA0);
      const kpi = [
        { k: 'Vrec', l: 'Plug with recycle', v: tRec * v0, u: KIT.volUnit(v._cu), s: 'R = ' + N.fmt(Rrec, 3) + ', τ = ' + N.fmt(tRec, 3) },
        { k: 'Rrec', l: 'Optimum recycle ratio', v: Rrec },
        { k: 'Vm1', l: 'One mixed tank', v: t1 * v0, u: KIT.volUnit(v._cu), s: 'τ = ' + N.fmt(t1, 3) },
        { k: 'Vm2', l: 'Two mixed tanks', v: t2 * v0, u: KIT.volUnit(v._cu), s: 'C_A between = ' + N.fmt(C1, 3) + '; ' + N.fmt(mfr(CA0, C1) * v0, 3) + ' + ' + N.fmt(mfr(C1, CAf) * v0, 3) },
        { k: 'Vmp', l: 'Mixed then plug', v: mp.f * v0, u: KIT.volUnit(v._cu), main: true, s: 'switch at C_A = ' + N.fmt(Cm, 3) + '; ' + N.fmt(mfr(CA0, Cm) * v0, 3) + ' + ' + N.fmt(pfr(CAf, Cm) * v0, 3) },
        { k: 'Vp', l: 'Plug flow alone', v: tPlug * v0, u: KIT.volUnit(v._cu) },
        { k: 'Vsep', l: 'Mixed at minimum, with separation', v: (CA0 - CAf) * inv(mn.x) * v0, u: KIT.volUnit(v._cu), s: 'operate at C_A = ' + N.fmt(mn.x, 3) + '; volume for the same throughput of A' }
      ];
      const warn = [];
      if (typeof outside !== 'undefined' && outside) warn.push('The duty runs outside the range of the test data; the curve is held flat beyond the last points.');
      const Cs = N.linspace(Math.min(CAf, pts.length ? pts[0].C : CAf) * 0.95, Math.max(CA0, pts.length ? pts[pts.length - 1].C : CA0), 260);
      const plot = { title: '1/(−r_A) against C_A', x: { label: 'C_A', unit: v._cu, min: 0 }, y: { label: '1/(−r_A)', unit: v._tu + '/(' + v._cu + ')', min: 0, max: Math.max(...Cs.filter(c => c >= CAf && c <= CA0).map(inv)) * 1.15 }, series: [{ name: '1/(−r_A)', x: Cs, y: Cs.map(inv), color: 1 }], marks: [{ t: 'v', x: CA0, label: 'C_{A0}', color: 8 }, { t: 'v', x: CAf, label: 'C_{Af}', color: 8 }] };
      if (pts.length) plot.series.push({ name: 'runs', type: 'pts', x: pts.map(p => p.C), y: pts.map(p => p.y), color: 2 });
      const area = (a, b, col) => { const sx = N.linspace(a, b, 80); plot.series.push({ x: sx, y: sx.map(inv), type: 'area', fill: 0.22, color: col, hover: false, legend: false, noRange: true }); };
      const rect = (a, b, h, col) => plot.marks.push({ t: 'rect', x0: a, x1: b, y0: 0, y1: h, color: col, fill: 0.2, noRange: true });
      let units = [];
      if (v.show === 'rec') { const Cin = (CA0 + Rrec * CAf) / (Rrec + 1); rect(CAf, CA0, tRec / (CA0 - CAf), 6); area(CAf, Cin, 3); units = [{ type: 'pfr', X: 1 - CAf / CA0, label: 'R = ' + N.fmt(Rrec, 3) }]; }
      else if (v.show === 'm1') { rect(CAf, CA0, inv(CAf), 2); units = [{ type: 'mfr', X: 1 - CAf / CA0 }]; }
      else if (v.show === 'm2') { rect(C1, CA0, inv(C1), 2); rect(CAf, C1, inv(CAf), 4); units = [{ type: 'mfr', X: 1 - C1 / CA0 }, { type: 'mfr', X: 1 - CAf / CA0 }]; }
      else if (v.show === 'mp') { rect(Cm, CA0, inv(Cm), 2); area(CAf, Cm, 3); units = [{ type: 'mfr', X: 1 - Cm / CA0 }, { type: 'pfr', X: 1 - CAf / CA0, prof: N.linspace(1 - Cm / CA0, 1 - CAf / CA0, 10) }]; }
      else { rect(CAf, CA0, inv(mn.x), 5); units = [{ type: 'mfr', X: 1 - mn.x / CA0, label: 'then separate' }]; }
      return { kpi, warn, plots: [plot], diagram: D.train(units.map((u, i) => Object.assign({ tag: 'R-60' + (i + 1) }, u)), { h: 126 }).svg, tables: pts.length ? [{ title: 'Points from the runs', cols: ['C_A', '1/(−r_A)'], rows: pts.map(p => [p.C, p.y]) }] : [] };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
