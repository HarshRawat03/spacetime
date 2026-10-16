/* Chapter 10 — Choosing the right kind of reactor: Rules 1–6 applied to any network */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, TH = ST.TH, D = ST.D;
  const C2K = TH.C2K, K2C = TH.K2C, R = TH.R;
  const ADV = ST.ADV = {};

  const nameOf = (net, r) => r.text.replace(/\s+/g, ' ');
  const Eof = r => (r.loc.ER !== undefined ? r.loc.ER * R : r.loc.E !== undefined ? +r.loc.E : null);
  const ordIn = (r, i) => { const t = r.reac.find(q => q[0] === i); return !t ? 0 : r.ord && r.ord[i] !== undefined ? r.ord[i] : t[1]; };

  /* Qualitative analysis of a network for a wanted product (Levenspiel Ch. 7–10). */
  ADV.analyze = function (net, want, key) {
    const iW = net.idx(want), iA = net.idx(key), ns = net.ns;
    if (iW < 0) throw new Error('Wanted product ' + want + ' is not in the network.');
    if (iA < 0) throw new Error('Key reactant ' + key + ' is not in the network.');
    /* species graph: X → Y through reaction j when X is a reactant and Y a product */
    const out = N.range(ns).map(() => []);
    net.rx.forEach((r, j) => r.reac.forEach(([x]) => r.prod.forEach(([y]) => out[x].push([y, j]))));
    const dist = new Array(ns).fill(Infinity);
    dist[iA] = 0;
    const q = [iA];
    while (q.length) { const x = q.shift(); out[x].forEach(([y]) => { if (dist[y] === Infinity) { dist[y] = dist[x] + 1; q.push(y); } }); }
    const reachW = new Array(ns).fill(false);
    reachW[iW] = true;
    let changed = true;
    while (changed) { changed = false; for (let x = 0; x < ns; x++) if (!reachW[x] && out[x].some(([y]) => reachW[y])) { reachW[x] = true; changed = true; } }
    if (dist[iW] === Infinity) throw new Error(want + ' cannot be formed from ' + key + ' in this network.');
    const route = N.range(ns).filter(x => dist[x] < Infinity && reachW[x]);
    const onPath = net.rx.map(r => r.reac.some(([x]) => route.indexOf(x) >= 0 && x !== iW) && r.prod.some(([y]) => reachW[y]));
    const pairs = [];
    route.forEach(x => {
      if (x === iW) return;
      const cons = net.rx.map((r, j) => j).filter(j => net.rx[j].reac.some(([s]) => s === x));
      const P = cons.filter(j => onPath[j]), L = cons.filter(j => !onPath[j]);
      P.forEach(p => L.forEach(l => pairs.push({ kind: 'parallel', at: x, stage: dist[x], good: p, bad: l })));
    });
    const formers = net.rx.map((r, j) => j).filter(j => net.rx[j].prod.some(([s]) => s === iW) && onPath[j]);
    const degr = net.rx.map((r, j) => j).filter(j => net.rx[j].reac.some(([s]) => s === iW));
    formers.forEach(f => degr.forEach(d => pairs.push({ kind: 'series', at: iW, stage: dist[iW], good: f, bad: d })));
    const lines = [], rows = [];
    const reqT = [], reqC = {};
    const sp = i => net.species[i];
    pairs.forEach(pq => {
      const g = net.rx[pq.good], b = net.rx[pq.bad];
      const Eg = Eof(g), Eb = Eof(b);
      let tAdv = 'E not given';
      if (Eg !== null && Eb !== null) {
        if (Math.abs(Eg - Eb) < 1e-9) tAdv = 'equal E: T has no effect';
        else { const hi = Eg > Eb; tAdv = hi ? 'high T' : 'low T'; reqT.push({ stage: pq.stage, dir: hi ? 1 : -1, why: (pq.kind === 'series' ? 'forming ' + sp(iW) + ' (' + nameOf(net, g) + ') against destroying it (' + nameOf(net, b) + ')' : 'at ' + sp(pq.at) + ': ' + nameOf(net, g) + ' against ' + nameOf(net, b)) }); }
      }
      /* concentration levels: compare orders of every species taking part */
      const cAdv = [];
      const seen = new Set();
      g.reac.concat(b.reac).forEach(([i]) => {
        if (seen.has(i)) return; seen.add(i);
        if (pq.kind === 'series' && i === iW) return;
        if (g.fn || b.fn) return;
        const og = pq.kind === 'series' && g.reac.every(([s]) => s !== i) ? 0 : ordIn(g, i), ob = ordIn(b, i);
        if (pq.kind === 'series' && og === 0 && ob === 0) return;
        if (Math.abs(og - ob) < 1e-12) return;
        const hi = og > ob;
        cAdv.push((hi ? 'high C_' : 'low C_') + sp(i));
        (reqC[sp(i)] = reqC[sp(i)] || []).push(hi ? 1 : -1);
      });
      rows.push([pq.kind === 'series' ? 'series' : 'parallel at ' + sp(pq.at), nameOf(net, g), nameOf(net, b), Eg === null ? '–' : N.fmt(Eg / 1000, 4), Eb === null ? '–' : N.fmt(Eb / 1000, 4), tAdv, cAdv.join(', ') || 'none', pq.stage]);
    });
    /* temperature plan */
    let plan = 'any', why = '';
    const noLoss = pairs.length === 0;
    if (noLoss) { plan = 'high'; why = 'Every reaction leads toward ' + sp(iW) + ', so only speed matters: the highest allowable temperature (Ex. 10.2b).'; }
    else if (!reqT.length) { plan = 'any'; why = 'The activation energies do not separate the competing steps, so temperature does not change the product distribution.'; }
    else {
      const stages = Array.from(new Set(reqT.map(q => q.stage))).sort((a, b) => a - b);
      const lev = stages.map(s => { const d = reqT.filter(q => q.stage === s).map(q => q.dir); return d.every(x => x > 0) ? 'high' : d.every(x => x < 0) ? 'low' : 'mid'; });
      const uniq = Array.from(new Set(lev));
      if (uniq.length === 1) { plan = uniq[0] === 'mid' ? 'intermediate' : uniq[0]; why = uniq[0] === 'mid' ? 'Requirements conflict at the same point of the network: an intermediate temperature is best (Eq. 9.32).' : 'All competing steps want a ' + uniq[0] + ' temperature (Rule 6).'; }
      else {
        const first = lev[0], last = lev[lev.length - 1];
        const rank = { low: 0, mid: 1, high: 2 };
        if (rank[first] > rank[last]) { plan = 'falling'; why = 'Early steps want ' + first + ' T, later ones ' + last + ' T: a falling temperature progression (Rule 6, Ex. 10.2c).'; }
        else if (rank[first] < rank[last]) { plan = 'rising'; why = 'Early steps want ' + first + ' T, later ones ' + last + ' T: a rising temperature progression (Ex. 10.2d).'; }
        else { plan = 'varying'; why = 'The requirement changes along the network (' + lev.join(' → ') + ').'; }
      }
    }
    reqT.forEach(q => lines.push('Rule 6, ' + q.why + ': ' + (q.dir > 0 ? 'high' : 'low') + ' T.'));
    /* contacting */
    const series = degr.length > 0 || route.some(x => x !== iA && x !== iW);
    const cLev = {};
    Object.keys(reqC).forEach(s => { const d = reqC[s]; cLev[s] = d.every(x => x > 0) ? 'high' : d.every(x => x < 0) ? 'low' : 'intermediate'; });
    Object.keys(cLev).forEach(s => lines.push('Rule 3: keep C_' + s + ' ' + cLev[s] + (cLev[s] === 'high' ? ' (no dilution, no backmixing; high pressure for gases).' : cLev[s] === 'low' ? ' (backmix, dilute, recycle, or feed it in gradually).' : ' (a particular concentration is best).')));
    if (series) lines.push('Rule 2: ' + (degr.length ? sp(iW) + ' is itself consumed' : 'the route passes through intermediates') + ', so do not mix fluids of different composition: plug flow or batch, no recycle.');
    const keyLev = cLev[key];
    const others = Object.keys(cLev).filter(s => s !== key);
    let flow, flowWhy;
    if (keyLev === 'intermediate') { flow = 'mixed at the best C_' + key; flowWhy = 'The wanted step has an intermediate order in ' + key + ': hold C_' + key + ' where φ peaks in mixed flow (Ex. 10.1), separating and recycling unused ' + key + ' for the best yield; for the most product, follow with plug flow (Ex. 7.4).'; }
    else if (keyLev === 'low' && !series) { flow = 'mixed'; flowWhy = 'Mixed flow keeps C_' + key + ' low (dilution or recycle also help).'; }
    else if (keyLev === 'low' && series) { flow = 'plug, ' + key + ' side-fed'; flowWhy = 'Keep C_' + key + ' low without backmixing the products: feed ' + key + ' gradually along a plug-flow reactor.'; }
    else { flow = 'plug'; flowWhy = 'Plug flow (or batch): ' + (keyLev === 'high' ? 'it keeps C_' + key + ' high' : 'no backmixing') + '.'; }
    others.forEach(s => {
      if (cLev[s] === 'low' && flow.indexOf('plug') === 0 && flow.indexOf('side-fed') < 0) { flow = 'plug, ' + s + ' side-fed'; flowWhy = 'Plug flow with ' + s + ' fed along the reactor keeps C_' + s + ' low while C_' + key + ' stays ' + (keyLev || 'as it is') + ' (Ex. 7.2c).'; }
      else if (cLev[s] === 'high' && (keyLev === 'low' || keyLev === 'intermediate')) { flow = 'plug of ' + s + ', ' + key + ' side-fed'; flowWhy = 'Start with all the ' + s + ' and add ' + key + ' gradually.'; }
      else if (cLev[s] === 'high' && flow === 'plug') flowWhy += ' Feed all the ' + s + ' at the inlet (or in excess).';
    });
    if (!series && keyLev === undefined && !others.length) { flow = noLoss ? 'plug' : 'any'; flowWhy = noLoss ? 'Plug flow for the smallest reactor.' : 'The contacting pattern does not change the distribution: choose for size and cost.'; }
    return { lines, rows, plan, why, flow, flowWhy, route: route.map(sp), pairs };
  };

  const EX102 = 'A -> R ; k = 1, E = 79000\nA -> T ; k = 1, E = 113000\nR -> U ; k = 1, E = 126000\nT -> S ; k = 1, E = 151000\nU -> S ; k = 1, E = 0';
  const ex = (want, plan) => ({ name: 'Example 10.2(' + 'abcd'['RSTU'.indexOf(want)] + '): ' + want + ' wanted', v: { net: EX102, key: 'A', want, feed: 'A = 1', size: false, numOn: true, Tref: 100, Tmin: 50, Tmax: 150, Nseg: 6, tauTot: 6 }, book: { Tplan: plan, flow: 'plug' }, note: 'A → R (79), A → T (113), R → U (126), T → S (151), U → S (0 kJ/mol). Book: R — ' + 'low T, plug flow; S — high T, plug flow; T — falling T, plug flow; U — rising T, plug flow.'.split('; ')['RSTU'.indexOf(want)].replace(/^[RSTU] — /, '') + '.' });

  ST.tool({
    id: 'choose', tag: 'H-1001', ch: 10, sec: 'Ch. 10', title: 'Choosing the right kind of reactor: Rules 1–6 applied to your network',
    lede: 'Write the reactions; name the product you want. Every competing pair is found — parallel at each species on the route, and series loss of the product — and judged by activation energy and by reaction order. Out come the temperature level or progression and the contacting pattern, checked if you like by optimizing a temperature profile numerically.',
    keys: 'choosing reactor rules contacting pattern temperature progression falling rising plug flow mixed flow side feed activation energy order series parallel example 10.2 heuristics advisor',
    units: { t: 'min' }, slow: true,
    inputs: [
      { group: 'Network' },
      { k: 'net', t: 'area', label: 'Reactions with activation energies', v: EX102, rows: 7, hint: 'Orders follow the stoichiometry unless set, e.g. n_A = 2. Give E in J/mol, or ER = E/R in K with k0.' },
      { k: 'key', t: 'text', label: 'Key reactant', v: 'A', mono: true },
      { k: 'want', t: 'text', label: 'Wanted product', v: 'T', mono: true },
      { k: 'size', t: 'check', label: 'Reactor size matters too', v: false },
      { group: 'Numerical check' },
      { k: 'numOn', t: 'check', label: 'Optimize the temperature profile of a plug-flow reactor', v: true },
      { k: 'feed', t: 'text', label: 'Feed concentrations', v: 'A = 1', mono: true, if: v => v.numOn },
      { k: 'Tref', label: 'k values are at', unit: '°C', v: 100, min: -100, max: 2000, if: v => v.numOn },
      { k: 'Tmin', label: 'Lowest allowed', unit: '°C', v: 50, min: -150, max: 2000, if: v => v.numOn },
      { k: 'Tmax', label: 'Highest allowed', unit: '°C', v: 150, min: -150, max: 2500, if: v => v.numOn },
      { k: 'Nseg', label: 'Segments of equal space time', v: 6, min: 2, max: 10, step: 1, if: v => v.numOn },
      { k: 'tauTot', label: 'Total space time', unit: '[t]', v: 6, min: 1e-6, max: 1e8, scale: 'log', if: v => v.numOn }
    ],
    presets: [ex('R', 'low'), ex('S', 'high'), ex('T', 'falling'), ex('U', 'rising'),
      { name: 'Problem 10.1: A + B → R, R + B → S', v: { net: 'A + B -> R ; k = 1, E = 50000\nR + B -> S ; k = 1, E = 50000', key: 'A', want: 'R', feed: 'A = 1, B = 1', size: false, numOn: false }, note: 'R wanted. B is needed by both steps to the same order, so its level does not matter; R is an intermediate, so plug flow of A, with B added any way you like.' },
      { name: 'Problem 10.4: A + B → R, A → S', v: { net: 'A + B -> R ; k = 1, E = 60000\nA -> S ; k = 1, E = 60000', key: 'A', want: 'R', feed: 'A = 1, B = 1', size: false, numOn: false }, note: 'High C_B helps R; A’s order is the same in both. Dump B in first and add A slowly, or plug flow of B with A fed along the way.' },
      { name: 'Trambouze (Ex. 10.1)', v: { net: 'A -> R ; k = 0.025, n_A = 0\nA -> S ; k = 0.2\nA -> T ; k = 0.4, n_A = 2', key: 'A', want: 'S', feed: 'A = 1', size: false, numOn: false }, note: 'Zero-, first- and second-order paths with S of intermediate order: hold C_A at its best value (0.25 mol/L) in mixed flow. Orders that differ from the stoichiometry are written n_A = 0.' }
    ],
    eqs: [
      ['\\text{Rule 1: } n>0 \\Rightarrow C \\text{ high};\\quad n<0 \\Rightarrow C \\text{ low}', 'Rule 1'],
      ['\\text{Rule 2: series} \\Rightarrow \\text{no mixing of fluids of different composition}', 'Rule 2'],
      ['\\text{Rule 3: } C_A \\uparrow \\text{ favours the higher order}', 'Rule 3'],
      ['\\frac{k_1}{k_2}\\propto e^{-(E_1-E_2)/RT} \\Rightarrow T\\uparrow \\text{ favours the larger } E', 'Rule 6'],
      ['\\Phi\\left(\\frac{R}{A}\\right) = \\frac{\\text{moles R formed}}{\\text{moles A consumed}},\\quad (\\text{Prod R}) = \\frac{\\text{moles R formed}}{\\text{moles A fed}}', 'Eqs. 10.1–10.2']
    ],
    method: ['The route to the wanted product is every species reachable from the key reactant from which the product is reachable. At each species on the route, steps that stay on the route compete in parallel with steps that leave it; the product competes in series with the steps that destroy it. Each pair is judged by activation energy (Rule 6) and by reaction order in every species taking part (Rules 1 and 3).', 'The distance from the key reactant along the route sets when a requirement acts: early requirements against late ones decide whether the temperature should fall or rise along a plug-flow reactor.', 'The numerical check splits a plug-flow reactor into equal segments with their own temperatures and maximizes the outlet concentration of the wanted product (exact matrix exponentials for first-order networks).'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const net = ST.parseNet(v.net, { Tref: C2K(+v.Tref || 25) });
      const want = String(v.want).trim(), key = String(v.key).trim();
      const a = ADV.analyze(net, want, key);
      const planTxt = { low: 'low', high: 'high', falling: 'falling', rising: 'rising', intermediate: 'intermediate', varying: 'varying', any: 'any' }[a.plan];
      kpi.push({ k: 'Tplan', l: 'Temperature', v: planTxt, main: true, s: a.plan === 'falling' || a.plan === 'rising' ? 'progression along a plug-flow reactor' : a.plan === 'any' ? 'no effect on the distribution' : 'level' });
      kpi.push({ k: 'flow', l: 'Contacting', v: a.flow, main: true });
      kpi.push({ k: 'route', l: 'Route to ' + want, v: a.route.join(', ') });
      info.push(a.why);
      info.push(a.flowWhy);
      a.lines.forEach(s => info.push(s));
      if (v.size && (a.plan === 'low' || a.plan === 'falling')) info.push('Reactor size matters: low temperatures cost volume. Trade the yield loss at a higher temperature against the extra reactor volume, or let the temperature fall only as the product builds up.');
      const tables = [{ title: 'Competing pairs', cols: ['kind', 'wanted step', 'competing step', 'E wanted, kJ/mol', 'E competing, kJ/mol', 'temperature', 'concentration', 'position on route'], rows: a.rows }];
      const plots = [];
      /* numerical check */
      if (v.numOn) {
        const NT = ST.NETT;
        const lin = NT.isLinear(net);
        const C0 = new Array(net.ns).fill(0);
        String(v.feed).split(/[,;\n]+/).forEach(t => { const mm = /^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*([-+0-9.eE]+)\s*$/.exec(t); if (!t.trim()) return; if (!mm) throw new Error('Feed: write it as “A = 1”.'); const i = net.idx(mm[1]); if (i < 0) throw new Error('Feed: ' + mm[1] + ' is not in the network.'); C0[i] = parseFloat(mm[2]); });
        const iW = net.idx(want);
        const Tmin = C2K(+v.Tmin), Tmax = C2K(+v.Tmax), ns = Math.round(+v.Nseg), dt = +v.tauTot / ns;
        if (!(Tmax > Tmin)) throw new Error('The highest temperature must be above the lowest.');
        const lg = x => 1 / (1 + Math.exp(-x)), Tof = z => Tmin + (Tmax - Tmin) * lg(z);
        const runT = Ts => { let C = C0.slice(); const outs = [C.slice()]; Ts.forEach(T => { C = NT.stage(net, lin, 'P', C, T, dt); outs.push(C.slice()); }); return { C, outs }; };
        const f = z => -(runT(z.map(Tof)).C[iW]);
        let best = null;
        [[-3], [0], [3]].forEach(s0 => { const x0 = new Array(ns).fill(s0[0]); const r = N.nelderMead(f, x0, { tol: 1e-11, maxIt: 300 * ns, step: x0.map(() => 1.5) }); if (!best || r.f < best.f) best = r; });
        [N.linspace(3, -3, ns), N.linspace(-3, 3, ns)].forEach(x0 => { const r = N.nelderMead(f, x0, { tol: 1e-11, maxIt: 300 * ns, step: x0.map(() => 1.0) }); if (r.f < best.f) best = r; });
        best = N.nelderMead(f, best.x, { tol: 1e-13, maxIt: 400 * ns, step: best.x.map(() => 0.3) });
        const Ts = best.x.map(Tof);
        const isoBest = N.fminScan(T => -runT(new Array(ns).fill(T)).C[iW], Tmin, Tmax, 60);
        const d = Ts[ns - 1] - Ts[0], tol = 0.03 * (Tmax - Tmin);
        const gain = (-best.f + isoBest.f) / Math.max(-best.f, 1e-300);
        const lvl = T => (T > Tmax - tol ? 'high' : T < Tmin + tol ? 'low' : 'flat');
        const numPlan = gain < 1e-4 ? lvl(isoBest.x) : Math.abs(d) < tol ? (Ts.every(T => T > Tmax - tol) ? 'high' : Ts.every(T => T < Tmin + tol) ? 'low' : 'flat') : d < 0 ? 'falling' : 'rising';
        kpi.push({ k: 'numPlan', l: 'Numerical optimum profile', v: numPlan, s: Ts.map(T => N.fmt(K2C(T), 3)).join(' → ') + ' °C' });
        kpi.push({ k: 'Cbest', l: 'C_' + want + ' leaving, optimal profile', v: -best.f, u: 'mol/L' }, { k: 'Ciso', l: 'C_' + want + ' leaving, best single T', v: -isoBest.f, u: 'mol/L', s: 'at ' + N.fmt(K2C(isoBest.x), 4) + ' °C' });
        if (numPlan !== a.plan && !(numPlan === 'flat' && (a.plan === 'intermediate' || a.plan === 'any'))) info.push('The numerical optimum (' + numPlan + ') differs from the rule-based plan (' + a.plan + '): with these rate constants and space time, the size of the constants and the reactor length matter as well as the activation energies.');
        const cum = N.range(ns + 1).map(i => i * dt);
        plots.push({ title: 'Optimal temperature profile along the plug-flow reactor', x: { label: 'τ', unit: v._tu || 'min', min: 0, max: +v.tauTot }, y: { label: 'T', unit: '°C', min: +v.Tmin, max: +v.Tmax }, series: [{ name: 'optimal', x: cum, y: Ts.concat([Ts[ns - 1]]).map(K2C), type: 'step', color: 1, w: 2 }, { name: 'best single temperature', x: [0, +v.tauTot], y: [K2C(isoBest.x), K2C(isoBest.x)], color: 8, dash: '5 4' }] });
        const res = runT(Ts);
        plots.push({ title: 'Concentrations along the reactor, optimal profile', x: { label: 'τ', unit: v._tu || 'min', min: 0 }, y: { label: 'C', unit: 'mol/L', min: 0 }, series: net.species.map((s, j) => ({ name: 'C_' + s, x: cum, y: res.outs.map(c => c[j]), color: j + 1, marker: 'o' })) });
      }
      /* network sketch: species laid out by their distance from the key reactant, route highlighted */
      const diag = (() => {
        const n = net.ns, iA = net.idx(key);
        const dist = new Array(n).fill(Infinity); dist[iA] = 0;
        const qq = [iA];
        while (qq.length) { const x = qq.shift(); net.rx.forEach(r => { if (r.reac.some(([s0]) => s0 === x)) r.prod.forEach(([y]) => { if (dist[y] === Infinity) { dist[y] = dist[x] + 1; qq.push(y); } }); }); }
        const maxL = Math.max(0, ...dist.filter(isFinite));
        const lev = dist.map(d => (isFinite(d) ? d : maxL + 1));
        const nL = Math.max(...lev) + 1;
        const W = Math.max(300, 130 * nL + 40), byL = N.range(nL).map(l => N.range(n).filter(i => lev[i] === l));
        const H = Math.max(120, 64 * Math.max(...byL.map(a => a.length)) + 40);
        const pos = new Array(n);
        byL.forEach((ids, l) => ids.forEach((i, k) => { pos[i] = [40 + l * (W - 80) / Math.max(1, nL - 1), H / 2 + (k - (ids.length - 1) / 2) * 64]; }));
        let s = '';
        net.rx.forEach(r => r.reac.forEach(([x]) => r.prod.forEach(([y]) => {
          const [x0, y0] = pos[x], [x1, y1] = pos[y], dx = x1 - x0, dy = y1 - y0, L = Math.hypot(dx, dy) || 1;
          const on = a.route.indexOf(net.species[x]) >= 0 && a.route.indexOf(net.species[y]) >= 0;
          s += D.pipe([[x0 + dx / L * 16, y0 + dy / L * 12], [x1 - dx / L * 18, y1 - dy / L * 13]], { flow: on });
          const E = Eof(r);
          if (E !== null) s += D.text((x0 + x1) / 2 - dy / L * 10, (y0 + y1) / 2 + dx / L * 10 + 3, N.fmt(E / 1000, 3) + (Math.abs(dy) < 1 ? ' kJ' : ''), { fs: 9.5 });
        })));
        pos.forEach((p, i) => { s += D.tag(p[0], p[1], net.species[i]); });
        return D.svg(W, H, s, 'reaction network with the route to the wanted product highlighted');
      })();
      return { kpi, warn, info, plots, tables, diagram: diag };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
