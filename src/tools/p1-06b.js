/* Chapter 6 — any series–parallel network of ideal reactors, and every such arrangement of a set of units */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, K = ST.K, KIT = ST.KIT, D = ST.D;

  /* grammar:  net := blk (('>' | ',') blk)*   blk := unit | '(' net ('|' net)+ ')' | '(' net ')'   unit := P|M number */
  function parse(src) {
    const s = String(src).replace(/\s+/g, ' ').trim();
    let i = 0;
    const peek = () => { while (s[i] === ' ') i++; return s[i]; };
    function net() {
      const parts = [blk()];
      while (peek() === '>' || peek() === ',') { i++; parts.push(blk()); }
      return parts.length === 1 ? parts[0] : { t: 's', c: flat('s', parts) };
    }
    function blk() {
      if (peek() === '(') {
        i++;
        const br = [net()];
        while (peek() === '|') { i++; br.push(net()); }
        if (peek() !== ')') throw new Error('Missing “)” near position ' + (i + 1));
        i++;
        return br.length === 1 ? br[0] : { t: 'p', c: flat('p', br) };
      }
      const m = /^(pfr|plug|p|mfr|cstr|mixed|m)\s*([0-9.]+(?:e[-+]?\d+)?)/i.exec(s.slice(i).replace(/^ +/, ''));
      if (!m) throw new Error('Expected a unit like “M 100” or “P 50” near position ' + (i + 1));
      while (s[i] === ' ') i++;
      i += m[0].length;
      return { t: 'u', type: /^p/i.test(m[1]) ? 'pfr' : 'mfr', V: parseFloat(m[2]) };
    }
    const flat = (t, arr) => arr.reduce((a, x) => a.concat(x.t === t ? x.c : [x]), []);
    const tree = net();
    if (peek() !== undefined) throw new Error('Unexpected “' + s[i] + '” at position ' + (i + 1));
    return tree;
  }
  const vol = n => (n.t === 'u' ? n.V : N.sum(n.c.map(vol)));
  const str = n => (n.t === 'u' ? (n.type === 'pfr' ? 'P' : 'M') + N.fmt(n.V, 4) : n.t === 's' ? n.c.map(c => (c.t === 'p' ? str(c) : str(c))).join(' > ') : '(' + n.c.map(str).join(' | ') + ')');
  const canon = n => (n.t === 'u' ? (n.type === 'pfr' ? 'P' : 'M') + n.V : n.t === 's' ? 'S[' + n.c.map(canon).join(',') + ']' : 'P{' + n.c.map(canon).sort().join(',') + '}');
  function parNodes(n, out = []) { if (n.t !== 'u') { if (n.t === 'p') out.push(n); n.c.forEach(c => parNodes(c, out)); } return out; }

  /* evaluate: conversion out for conversion in, fraction of total feed f, total feed v0 */
  let CACHE = new Map();
  function unitX(m, type, tau, Xin) {
    const key = type + '|' + tau.toPrecision(10) + '|' + Xin.toPrecision(10);
    let X = CACHE.get(key);
    if (X === undefined) { X = type === 'mfr' ? K.XMFR(m, tau, Xin) : K.XPFR(m, tau, Xin).X; if (CACHE.size > 20000) CACHE = new Map(); CACHE.set(key, X); }
    return X;
  }
  function evalNet(m, n, Xin, f, v0, splits, trace) {
    if (n.t === 'u') {
      const tau = n.V / (f * v0);
      const X = unitX(m, n.type, tau, Xin);
      if (trace) trace.push({ n, Xin, X, f, tau });
      return X;
    }
    if (n.t === 's') { let X = Xin; n.c.forEach(c => { X = evalNet(m, c, X, f, v0, splits, trace); }); return X; }
    const sp = splits.get(n) || n.c.map(c => vol(c) / vol(n));
    let X = 0;
    n.c.forEach((c, j) => { X += sp[j] > 1e-12 ? sp[j] * evalNet(m, c, Xin, f * sp[j], v0, splits, trace) : sp[j] * Xin; });
    return X;
  }
  function optimize(m, tree, v0, Xin, full, quick) {
    const pn = parNodes(tree);
    const splits = new Map();
    if (!pn.length || !full) return { X: evalNet(m, tree, Xin, 1, v0, splits), splits };
    const dims = pn.map(p => p.c.length - 1), nz = N.sum(dims);
    const apply = z => {
      let k = 0;
      pn.forEach((p, i) => {
        const base = p.c.map(c => Math.log(vol(c)));
        const e = p.c.map((c, j) => Math.exp(base[j] + (j ? z[k + j - 1] : 0)));
        const s = N.sum(e);
        splits.set(p, e.map(x => x / s));
        k += dims[i];
      });
    };
    const obj = z => { apply(z); const X = evalNet(m, tree, Xin, 1, v0, splits); return isFinite(X) ? -X : 0; };
    let best = { x: new Array(nz).fill(0), f: obj(new Array(nz).fill(0)) };
    const starts = [new Array(nz).fill(0)];
    for (let r = 0; r < (quick ? 1 : 3); r++) starts.push(new Array(nz).fill(0).map((_, j) => (((r + j) % 3) - 1) * 1.5));
    starts.forEach(s0 => { const r = N.nelderMead(obj, s0, { maxIt: (quick ? 80 : 300) * nz, step: new Array(nz).fill(0.8), tol: quick ? 1e-7 : 1e-10 }); if (r.f < best.f) best = r; });
    apply(best.x);
    return { X: -best.f, splits };
  }

  /* enumerate series–parallel arrangements of a multiset of units */
  function enumerate(units) {
    const memo = new Map();
    const keyOf = idx => idx.map(i => units[i].type + units[i].V).sort().join(',');
    function gen(idx, notType) {
      const key = keyOf(idx) + '|' + (notType || '');
      if (memo.has(key)) return memo.get(key);
      let out = [];
      if (idx.length === 1) out = [Object.assign({ t: 'u' }, units[idx[0]])];
      else {
        if (notType !== 's') {
          /* ordered splits into first block + rest */
          const n = idx.length, seen = new Set();
          for (let mask = 1; mask < (1 << n) - 1; mask++) {
            const a = [], b = [];
            for (let j = 0; j < n; j++) ((mask >> j) & 1 ? a : b).push(idx[j]);
            const k = keyOf(a) + '/' + keyOf(b);
            if (seen.has(k)) continue;
            seen.add(k);
            gen(a, 's').forEach(x => gen(b, null).forEach(y => {
              const tail = y.t === 's' ? y.c : [y];
              out.push({ t: 's', c: [x].concat(tail) });
            }));
          }
        }
        if (notType !== 'p') {
          const n = idx.length, seen = new Set();
          for (let mask = 1; mask < (1 << n) - 1; mask++) {
            if (!(mask & 1)) continue;
            const a = [], b = [];
            for (let j = 0; j < n; j++) ((mask >> j) & 1 ? a : b).push(idx[j]);
            const k = keyOf(a) + '/' + keyOf(b);
            if (seen.has(k)) continue;
            seen.add(k);
            gen(a, 'p').forEach(x => gen(b, null).forEach(y => {
              const rest = y.t === 'p' ? y.c : [y];
              out.push({ t: 'p', c: [x].concat(rest) });
            }));
          }
        }
      }
      const uniq = new Map();
      out.forEach(o => { const c = canon(o); if (!uniq.has(c)) uniq.set(c, o); });
      out = [...uniq.values()];
      memo.set(key, out);
      return out;
    }
    const all = gen(units.map((_, i) => i), null);
    const uniq = new Map();
    all.forEach(o => { const c = canon(o); if (!uniq.has(c)) uniq.set(c, o); });
    return [...uniq.values()];
  }
  function leaves(n, out = []) { if (n.t === 'u') out.push(n); else n.c.forEach(c => leaves(c, out)); return out; }

  /* schematic of a series-parallel tree */
  function layout(n, splits, map) {
    if (n.t === 'u') { const w = n.type === 'pfr' ? 74 : 44; return { w, h: 58, n }; }
    const kids = n.c.map(c => layout(c, splits, map));
    if (n.t === 's') return { w: N.sum(kids.map(k => k.w)) + 18 * (kids.length - 1), h: Math.max(...kids.map(k => k.h)), n, kids };
    return { w: Math.max(...kids.map(k => k.w)) + 36, h: N.sum(kids.map(k => k.h)) + 6 * (kids.length - 1), n, kids };
  }
  function drawL(L, x, ym, out, xs) {
    const n = L.n;
    if (n.t === 'u') {
      const Xo = xs.get(n);
      if (n.type === 'pfr') out.push(D.tube(x, ym - 10, L.w, 20, { prof: Xo ? [Xo.Xin, Xo.X] : null }));
      else out.push(D.tank(x + 2, ym - 21, L.w - 4, 42, { fill: Xo ? Xo.X : 0, stir: false }));
      out.push(D.text(x + L.w / 2, ym + 34, (n.type === 'pfr' ? 'P' : 'M') + ' ' + N.fmt(n.V, 3), { cls: 'pid-l' }));
      return;
    }
    if (n.t === 's') {
      let cx = x;
      L.kids.forEach((k, i) => { drawL(k, cx, ym, out, xs); cx += k.w; if (i < L.kids.length - 1) { out.push(D.pipe([[cx, ym], [cx + 16, ym]], { arrow: false })); cx += 18; } });
      return;
    }
    let cy = ym - L.h / 2;
    const xr = x + L.w;
    L.kids.forEach(k => {
      const yk = cy + k.h / 2;
      out.push(D.pipe([[x, ym], [x, yk], [x + 18 + (L.w - 36 - k.w) / 2, yk]], { arrow: false }));
      drawL(k, x + 18 + (L.w - 36 - k.w) / 2, yk, out, xs);
      out.push(D.pipe([[x + 18 + (L.w - 36 - k.w) / 2 + k.w, yk], [xr, yk], [xr, ym]], { arrow: false }));
      cy += k.h + 6;
    });
    out.push(D.node(x, ym), D.node(xr, ym));
  }

  ST.tool({
    id: 'network', tag: 'H-606', ch: 6, sec: '6.2', title: 'Any series–parallel network of reactors',
    lede: 'Connect plug and mixed flow units in any combination of series and parallel, with the feed split between branches optimized; or let the tool build every series–parallel arrangement of your units and rank them.',
    keys: 'reactor network series parallel combination arrangement connected as you wish best arrangement four tanks autocatalytic feed split optimization every combination',
    units: { C: 'mol/L', t: 's' }, slow: true, sweep: false,
    inputs: [
      { group: 'Network' },
      { k: 'arr', t: 'text', label: 'Arrangement (volumes)', v: '(M 100 | M 100 | M 100) > M 100', mono: true, hint: 'Units: “M 100” mixed, “P 50” plug. Series: “>”. Parallel: “( … | … )”. Example: (P 50 > P 30 | P 40).' },
      { k: 'v0', label: 'Total feed rate v₀', unit: v => KIT.volUnit(v._cu) + '/' + v._tu, v: 1.5, min: 1e-9, max: 1e9, scale: 'log' },
      { k: 'opt', t: 'check', label: 'Optimize the feed split at every branch point', v: true, hint: 'Off: split in proportion to branch volume (equal V/F).' },
      { k: 'enum', t: 'check', label: 'Rank every series–parallel arrangement of these units', v: true, hint: 'Up to five units.' }
    ].concat(KIT.kin({ d: { law: 'auto', k: 0.001, CA0: 10, MR: 0 } })),
    presets: [
      { name: 'Problem 6.11', v: { arr: 'M 100 > M 100 > M 100 > M 100', v0: 1.5, opt: true, enum: true, law: 'auto', k: 0.001, CA0: 10, MR: 0, r: 1, eps: 0, _cu: 'mol/L', _tu: 's' }, note: 'Four 100 L mixed tanks, 1.5 L/s of 10 mol/L A, −r_A = 0.001 C_AC_R. In series nothing reacts; the ranking finds the arrangement that does.' },
      { name: 'Example 6.1 as a network', v: { arr: '(P 50 > P 30 | P 40)', v0: 40, opt: true, enum: false, law: 'power', a: 1, ob: 0, k: 0.05, CA0: 1, eps: 0, _cu: 'mol/L', _tu: 'min' }, note: 'The optimized split sends two-thirds of the feed to the 80 L branch.' },
      { name: 'Second order, mixed and plug', v: { arr: 'M 2 > P 2 > M 4', v0: 1, opt: true, enum: true, law: 'power', a: 2, ob: 0, k: 1, CA0: 1, eps: 0, _tu: 'min' } }
    ],
    eqs: [
      ['X_{mix} = \\sum_j s_jX_j,\\qquad \\tau_{unit} = \\frac{V}{s\\,v_0}', '§6.2', 'Streams leaving parallel branches mix in proportion to their share s of the feed.'],
      ['\\text{plug branches: } \\frac{V_1}{F_1} = \\frac{V_2}{F_2} \\text{ is optimal}', '§6.2']
    ],
    method: ['The arrangement is parsed into a tree of series and parallel blocks. Every unit is solved on the conversion basis of the original feed, with space-time V/(s v₀) for the share s of the feed it receives. Branch splits are optimized together by Nelder–Mead from several starting points.', 'For the ranking, all distinct series–parallel structures of the units are generated (series blocks ordered, parallel blocks unordered: 10 structures for four identical units, 24 for five), each solved with its own optimized splits.'],
    run(v) {
      CACHE = new Map();
      const m = KIT.model(v), v0 = +v.v0;
      const tree = parse(v.arr);
      const res = optimize(m, tree, v0, 0, !!v.opt);
      const trace = [];
      const Xf = evalNet(m, tree, 0, 1, v0, res.splits, trace);
      const kpi = [{ k: 'X', l: 'Outlet X_A', v: Xf, main: true }, { k: 'CA', l: 'Outlet C_A', v: m.conc(Xf).A, u: v._cu }, { k: 'V', l: 'Total volume', v: vol(tree), u: KIT.volUnit(v._cu) }];
      const warn = [], tables = [];
      const pn = parNodes(tree);
      if (pn.length) tables.push({ title: 'Feed splits', cols: ['branch point', 'branches', 'share of the stream entering it'], rows: pn.map((p, i) => [i + 1, p.c.map(str).join('  |  '), (res.splits.get(p) || p.c.map(c => vol(c) / vol(p))).map(x => N.fmt(x, 3)).join(', ')]) });
      tables.push({ title: 'Units', cols: ['unit', 'share of feed', 'τ', 'X_A in', 'X_A out'], rows: trace.map(t => [str(t.n), t.f, t.tau, t.Xin, t.X]) });
      if (v.enum) {
        const us = leaves(tree);
        if (us.length > 5) warn.push('Ranking is limited to five units.');
        else {
          const all = enumerate(us.map(u => ({ type: u.type, V: u.V })));
          const scored = all.map(a => { const r = optimize(m, a, v0, 0, !!v.opt, true); return { a, X: r.X }; }).sort((p, q) => q.X - p.X);
          kpi.push({ k: 'best', l: 'Best arrangement', v: str(scored[0].a), s: 'X_A = ' + N.fmt(scored[0].X, 4) }, { k: 'count', l: 'Arrangements tried', v: scored.length });
          tables.unshift({ title: 'Every series–parallel arrangement, best first', cols: ['rank', 'arrangement', 'X_A'], rows: scored.map((o, i) => [i + 1, str(o.a), o.X]) });
        }
      }
      const L = layout(tree);
      const out = [];
      const xs = new Map(trace.map(t => [t.n, t]));
      const H = Math.max(80, L.h + 46), ym = H / 2 - 4;
      out.push(D.pipe([[4, ym], [28, ym]], { label: 'feed' }));
      drawL(L, 30, ym, out, xs);
      out.push(D.pipe([[30 + L.w, ym], [30 + L.w + 40, ym]]));
      out.push(D.text(30 + L.w + 44, ym + 4, 'X ' + N.fmt(Xf, 3), { cls: 'pid-r', anchor: 'start' }));
      return { kpi, warn, tables, diagram: D.svg(L.w + 120, H, out.join(''), 'reactor network') };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
