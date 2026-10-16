/* Chapters 7–10 — reaction network workbench: any scheme in batch, plug, mixed or staged flow; optimum τ; fitting */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, D = ST.D, NR = ST.NR;

  function parseData(text) {
    const lines = String(text || '').trim().split(/\r?\n/).map(l => l.trim()).filter(Boolean);
    if (lines.length < 2) return null;
    const head = lines[0].split(/[\t,; ]+/).filter(Boolean);
    const rows = lines.slice(1).map(l => l.split(/[\t,; ]+/).filter(Boolean).map(Number)).filter(r => r.length === head.length && r.every(isFinite));
    return { head, rows };
  }
  const constList = net => {
    const out = [];
    net.rx.forEach((r, j) => { if (!r.custom || r.loc.k !== undefined) out.push({ j, f: 'k', name: 'k' + (j + 1) }); if (r.rev) out.push({ j, f: 'kr', name: 'kr' + (j + 1) }); });
    return out;
  };
  const getK = (net, c) => (c.f === 'k' ? (net.rx[c.j].custom ? net.rx[c.j].loc.k : net.rx[c.j].k) : net.rx[c.j].kr);
  const setK = (net, c, val) => { if (c.f === 'k') { if (net.rx[c.j].custom) net.rx[c.j].loc.k = val; else net.rx[c.j].k = val; } else net.rx[c.j].kr = val; };

  ST.tool({
    id: 'netlab', tag: 'H-805', ch: 8, sec: 'Ch. 7–10', title: 'Reaction network workbench',
    lede: 'Type any scheme of reactions, elementary or with your own rate expressions, and run it in a batch or plug flow, a mixed reactor or a cascade: histories, maxima of every intermediate, the best space-time, and least-squares rate constants from your own data.',
    keys: 'reaction network multiple reactions scheme simulator series parallel reversible fit rate constants data maxima optimum residence time cascade denbigh trambouze sodium sulfide',
    units: { C: 'mol/L', t: 'min' }, slow: true, sweep: false,
    inputs: [
      { group: 'Scheme' },
      { k: 'net', t: 'area', label: 'Reactions', rows: 6, v: 'A -> S ; k = 0.0070\nA -> R ; k = 0.0108\nA -> T ; k = 0.0027\nS -> R ; k = 0.0099\nT -> R ; k = 0.0163', hint: 'One per line, “A + B -> R ; k = 1”. Reversible “<=>” with kr. Own rate: “; k = 2, r = k*A^2/(1+A)”.' },
      { k: 'feed', t: 'text', label: 'Feed or initial concentrations', v: 'A = 185', mono: true },
      { group: 'Reactor' },
      { k: 'mode', t: 'select', label: 'Run as', opts: [['batch', 'batch or plug flow'], ['mfr', 'mixed flow, range of τ'], ['tanks', 'equal mixed tanks, range of τ'], ['fit', 'fit rate constants to data']], v: 'batch' },
      { k: 'tEnd', label: 'Time or largest τ', unit: '[t]', v: 300, min: 1e-6, max: 1e9, scale: 'log', if: v => v.mode !== 'fit' },
      { k: 'Nt', label: 'Number of tanks', v: 3, min: 2, max: 50, step: 1, if: v => v.mode === 'tanks' },
      { k: 'tRep', label: 'Report compositions at', unit: '[t]', v: 30, min: 0, max: 1e9, if: v => v.mode !== 'fit' },
      { k: 'target', t: 'text', label: 'Species to maximize', v: 'S', mono: true, if: v => v.mode !== 'fit' },
      { k: 'logx', t: 'check', label: 'Log time axis', v: false, if: v => v.mode === 'batch' },
      { group: 'Data to fit', if: v => v.mode === 'fit' },
      { k: 'fitR', t: 'seg', label: 'Data are from', opts: [['batch', 'batch runs (t, C…)'], ['mfr', 'mixed flow runs (τ, C…)']], v: 'mfr', if: v => v.mode === 'fit' },
      { k: 'data', t: 'area', label: 'Data (first row: names)', rows: 6, v: 'tau, A, R, S\n5, 50, 40, 10\n20, 20, 40, 40', if: v => v.mode === 'fit', hint: 'Paste columns from a spreadsheet; the first column is time or space-time, the others are species in the scheme.' }
    ],
    presets: [
      { name: 'Example 8.3', v: { net: 'A -> S ; k = 0.0070\nA -> R ; k = 0.0108\nA -> T ; k = 0.0027\nS -> R ; k = 0.0099\nT -> R ; k = 0.0163', feed: 'A = 185', mode: 'batch', tEnd: 300, tRep: 30, target: 'S', _cu: 'mmol/L', _tu: 'min' }, book: { at_A: { v: 100, tol: 0.01 }, max_S: { v: 32, tol: 0.04 }, max_T: { v: 10, tol: 0.06 } }, note: 'The network the book builds for the oxidation of Na₂S (A) to Na₂S₂O₃ (R) through S and T, with its fitted constants. The checks compare with the curves the constants were fitted to.' },
      { name: 'Problem 8.5 (fit)', v: { net: 'A -> R ; k = 0.1\nR -> S ; k = 0.1', feed: 'A = 100', mode: 'fit', fitR: 'mfr', data: 'tau, A, R, S\n5, 50, 40, 10\n20, 20, 40, 40', _tu: 'min' }, note: 'Two mixed flow runs from pure A at 100: try A → R → S and see that both runs are matched exactly.' },
      { name: 'Reversible series (Fig. 8.9)', v: { net: 'A <=> R ; k = 1, kr = 0.3\nR <=> S ; k = 0.5, kr = 0.2', feed: 'A = 1', mode: 'batch', tEnd: 15, tRep: 2, target: 'R' } },
      { name: 'Trambouze in a cascade', v: { net: 'A -> R ; k = 0.025, r = k*A/(A + 0.0001)\nA -> S ; k = 0.2\nA + A -> T + A ; k = 0.4', feed: 'A = 1', mode: 'tanks', Nt: 4, tEnd: 100, tRep: 7.5, target: 'S' }, note: 'Example 10.1 kinetics: zero-order A → R written as a custom rate; the cascade maximizes S.' },
      { name: 'Problem 10.11 (Denbigh, mixed orders)', v: { net: 'A + A -> R + A ; k = 1\nA -> T ; k = 0.6\nR -> S ; k = 0.6\nR + R -> U + R ; k = 0.1', feed: 'A = 6, R = 0.6', mode: 'batch', tEnd: 5, tRep: 0.5, target: 'R', _tu: 's' }, note: 'Second-order first and last steps written as A + A → R + A so that the rate is k[A]².' }
    ],
    eqs: [
      ['\\frac{dC_i}{dt} = \\sum_j \\nu_{ij}r_j,\\qquad r_j = k_j\\prod C^{\\nu}\\;\\text{(elementary)}', 'Ch. 2, 8', 'Batch or plug flow.'],
      ['C_{i0}-C_i + \\tau\\sum_j\\nu_{ij}r_j(C) = 0', 'Eq. 5.13', 'Mixed flow: solved by Newton’s method with continuation in τ.']
    ],
    method: ['Batch and plug flow use Dormand–Prince with automatic switching to a stiff Rosenbrock method. Mixed flow steady states are found with damped Newton iterations, stepping τ up gradually from very small values so that the solution followed is the one reached on start-up.', 'Fitting adjusts every rate constant (on a log scale, so they stay positive) to minimize the squared deviations of the measured species, each scaled by its largest measured value; 95% intervals come from the Jacobian at the optimum.'],
    run(v) {
      const net = ST.parseNet(v.net);
      const feed = ST.K.parseParams(v.feed);
      Object.keys(feed).forEach(s => { if (net.idx(s) < 0) throw new Error('“' + s + '” is not in the scheme.'); });
      const C0 = net.species.map(s => feed[s] || 0);
      const kpi = [], warn = [], plots = [], tables = [];
      const sp = net.species;
      if (v.mode === 'fit') {
        const d = parseData(v.data);
        if (!d || d.rows.length < 1) throw new Error('Give a header row and at least one data row.');
        const cols = d.head.slice(1).map(h => net.idx(h));
        if (cols.some(c => c < 0)) throw new Error('Column names must be species in the scheme: ' + d.head.slice(1).join(', '));
        const scale = cols.map((c, j) => Math.max(...d.rows.map(r => Math.abs(r[j + 1])), 1e-12));
        const cl = constList(net);
        const p0 = cl.map(c => getK(net, c));
        const model = p => {
          p.forEach((x, i) => setK(net, cl[i], x));
          if (v.fitR === 'mfr') return d.rows.map(r => NR.mfr(net, C0, r[0]).C);
          const tEnd = Math.max(...d.rows.map(r => r[0]));
          const sol = N.ode((t, C) => net.prod(C), 0, C0, tEnd, { rtol: 1e-9, atol: 1e-12 });
          return d.rows.map(r => sol.at(r[0]));
        };
        const res = N.lm(p => { const M = model(p); const out = []; d.rows.forEach((r, i) => cols.forEach((c, j) => out.push((M[i][c] - r[j + 1]) / scale[j]))); return out; }, p0, { pos: p0.map(() => true) });
        const M = model(res.p);
        cl.forEach((c, i) => kpi.push({ k: c.name, l: c.name + ' (' + net.rx[c.j].text + ')', v: res.p[i], main: i === 0, s: isFinite(res.se[i]) ? '± ' + N.fmt(1.96 * res.se[i], 2) : '' }));
        kpi.push({ k: 'rms', l: 'rms deviation (scaled)', v: Math.sqrt(res.ssr / Math.max(1, res.n)) });
        tables.push({ title: 'Data against the fitted scheme', cols: [d.head[0]].concat(d.head.slice(1).map(h => h + ' data'), d.head.slice(1).map(h => h + ' fit')), rows: d.rows.map((r, i) => [r[0]].concat(r.slice(1), cols.map(c => M[i][c]))) });
        if (res.n <= cl.length) warn.push('There are no more data points than constants, so the fit is exact or undetermined and the intervals mean little.');
        if (v.fitR === 'batch') {
          const tEnd = Math.max(...d.rows.map(r => r[0]));
          const sol = N.ode((t, C) => net.prod(C), 0, C0, tEnd, { rtol: 1e-9, atol: 1e-12 });
          const ts = N.linspace(0, tEnd, 160);
          const series = [];
          cols.forEach((c, j) => { series.push({ name: sp[c] + ' fit', x: ts, y: ts.map(t => sol.at(t)[c]), color: j + 1 }, { name: sp[c] + ' data', type: 'pts', x: d.rows.map(r => r[0]), y: d.rows.map(r => r[j + 1]), color: j + 1, legend: false }); });
          plots.push({ title: 'Fit', x: { label: 't', unit: v._tu, min: 0 }, y: { label: 'C', unit: v._cu, min: 0 }, series });
        } else {
          const taus = N.logspace(Math.min(...d.rows.map(r => r[0])) / 10, Math.max(...d.rows.map(r => r[0])) * 10, 80);
          const sts = taus.map(t => NR.mfr(net, C0, t).C);
          const series = [];
          cols.forEach((c, j) => { series.push({ name: sp[c] + ' fit', x: taus, y: sts.map(C => C[c]), color: j + 1 }, { name: sp[c] + ' data', type: 'pts', x: d.rows.map(r => r[0]), y: d.rows.map(r => r[j + 1]), color: j + 1, legend: false }); });
          plots.push({ title: 'Fit', x: { label: 'τ', unit: v._tu, log: true }, y: { label: 'C', unit: v._cu, min: 0 }, series });
        }
        return { kpi, warn, plots, tables };
      }
      const tEnd = +v.tEnd, tgt = net.idx(String(v.target).trim());
      let curve, xlab, xlog = false, at;
      if (v.mode === 'batch') {
        const sol = N.ode((t, C) => net.prod(C), 0, C0, tEnd, { rtol: 1e-9, atol: 1e-12, minPts: 80 });
        const ts = v.logx ? [0].concat(N.logspace(tEnd * 1e-4, tEnd, 220)) : N.linspace(0, tEnd, 240);
        curve = { x: ts, C: ts.map(sol.at) }; xlab = 't'; xlog = !!v.logx;
        at = sol.at(Math.min(+v.tRep, tEnd));
        const fine = N.linspace(0, tEnd, 1200).map(t => [t, sol.at(t)]);
        sp.forEach((s, i) => {
          let best = fine[0];
          fine.forEach(f => { if (f[1][i] > best[1][i]) best = f; });
          if (best[0] > 0 && best[0] < tEnd * 0.999) { const o = N.fmin(t => -sol.at(t)[i], Math.max(0, best[0] - tEnd / 1200), Math.min(tEnd, best[0] + tEnd / 1200)); best = [o.x, sol.at(o.x)]; }
          kpi.push({ k: 'max_' + s, l: 'max C_' + s, v: best[1][i], u: v._cu, s: 'at t = ' + N.fmt(best[0], 4), main: i === tgt });
        });
      } else {
        const taus = N.logspace(tEnd * 1e-4, tEnd, 160);
        const Nn = v.mode === 'tanks' ? Math.round(+v.Nt) : 1;
        let prevs = null;
        const run = (tau, warm) => { let C = C0.slice(); const outs = []; for (let i = 0; i < Nn; i++) { C = NR.mfr(net, C, tau / Nn, warm && prevs ? { guess: prevs[i] } : {}).C; outs.push(C); } if (warm) prevs = outs; return C; };
        curve = { x: taus, C: taus.map(t => run(t, true)) }; xlab = Nn > 1 ? 'total τ' : 'τ'; xlog = true;
        at = run(Math.max(+v.tRep, 1e-12));
        sp.forEach((s, i) => {
          let bi = 0; curve.C.forEach((C, j) => { if (C[i] > curve.C[bi][i]) bi = j; });
          let bt = taus[bi], bv = curve.C[bi][i];
          if (bi > 0 && bi < taus.length - 1) { const o = N.fmin(lt => -run(Math.exp(lt))[i], Math.log(taus[bi - 1]), Math.log(taus[bi + 1])); bt = Math.exp(o.x); bv = -o.f; }
          kpi.push({ k: 'max_' + s, l: 'max C_' + s, v: bv, u: v._cu, s: 'at ' + xlab + ' = ' + N.fmt(bt, 4), main: i === tgt });
        });
      }
      sp.forEach((s, i) => kpi.push({ k: 'at_' + s, l: 'C_' + s + ' at ' + N.fmt(+v.tRep, 4), v: at[i], u: v._cu }));
      if (tgt < 0) warn.push('Species to maximize is not in the scheme.');
      const xs = xlog && v.mode === 'batch' ? curve.x.slice(1) : curve.x, Cs = xlog && v.mode === 'batch' ? curve.C.slice(1) : curve.C;
      plots.push({ title: v.mode === 'batch' ? 'Batch or plug flow' : v.mode === 'mfr' ? 'Mixed flow' : Math.round(+v.Nt) + ' equal mixed tanks', wide: true, x: { label: xlab, unit: v._tu, log: xlog, min: xlog ? undefined : 0 }, y: { label: 'C', unit: v._cu, min: 0 }, series: sp.map((s, i) => ({ name: 'C_' + s, x: xs, y: Cs.map(C => C[i]), color: i + 1, w: i === tgt ? 2.8 : 1.8 })), marks: [{ t: 'v', x: +v.tRep, color: 8 }] });
      tables.push({ title: 'Reactions as read', cols: ['#', 'reaction', 'k', 'kr', 'rate'], rows: net.rx.map((r, j) => [j + 1, r.text, r.custom ? (r.loc.k !== undefined ? r.loc.k : '') : r.k, r.rev ? r.kr : '', r.custom ? r.custom : 'mass action']) });
      return { kpi, warn, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
