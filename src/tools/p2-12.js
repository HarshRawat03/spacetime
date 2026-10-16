/* Chapter 12 — Compartment models: build any network of flow regions, fit it to a tracer curve, react in it */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, RT = ST.RT, KIT = ST.KIT, K = ST.K, D = ST.D, RTDI = ST.RTDI;

  /* grammar
       series := term ('>' term)*
       term   := unit | '(' branch ('|' branch)* ')' | 'R' val '[' series (';' series)? ']'
       branch := [val ':'] series
       unit   := 'P' val | 'M' val | 'T' val val | 'D' val val | 'L' val | 'B'
       val    := number | name                                                   */
  function parse(src) {
    const s = String(src).replace(/\s+/g, ' ').trim();
    let i = 0;
    const ws = () => { while (s[i] === ' ') i++; };
    const peek = () => { ws(); return s[i]; };
    const err = m => { throw new Error(m + ' (at position ' + (i + 1) + ')'); };
    const val = () => {
      ws();
      const m = /^([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?(?:\/[0-9]*\.?[0-9]+)?|[A-Za-z_][A-Za-z0-9_]*)/.exec(s.slice(i));
      if (!m) err('Expected a number or a parameter name');
      i += m[0].length;
      if (/^[A-Za-z_]/.test(m[0])) return { p: m[0] };
      if (m[0].indexOf('/') > 0) { const [a, b] = m[0].split('/'); return { n: parseFloat(a) / parseFloat(b) }; }
      return { n: parseFloat(m[0]) };
    };
    function series() {
      const c = [term()];
      while (peek() === '>') { i++; c.push(term()); }
      return c.length === 1 ? c[0] : { t: 's', c };
    }
    function branch() {
      const save = i;
      ws();
      const m = /^([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?|[A-Za-z_][A-Za-z0-9_]*)\s*:/.exec(s.slice(i));
      let f = null;
      if (m) { i += m[0].length; f = /^[A-Za-z_]/.test(m[1]) ? { p: m[1] } : { n: parseFloat(m[1]) }; }
      else i = save;
      return { f, n: series() };
    }
    function term() {
      const ch = peek();
      if (ch === '(') {
        i++;
        const br = [branch()];
        while (peek() === '|') { i++; br.push(branch()); }
        if (peek() !== ')') err('Missing “)”');
        i++;
        if (br.length === 1 && !br[0].f) return br[0].n;
        return { t: 'p', c: br.map(b => b.n), f: br.map(b => b.f) };
      }
      ws();
      const m = /^([PMTDLBR])\b/i.exec(s.slice(i));
      if (!m) err('Expected a unit: P, M, T, D, L, B or R');
      i += m[0].length;
      const k = m[1].toUpperCase();
      if (k === 'R') {
        const R = val();
        if (peek() !== '[') err('Recycle needs a block in brackets, e.g. R 1 [ P 1 ]');
        i++;
        const fwd = series();
        let back = null;
        if (peek() === ';') { i++; back = series(); }
        if (peek() !== ']') err('Missing “]”');
        i++;
        return { t: 'r', R, c: fwd, back };
      }
      if (k === 'B') return { t: 'u', k };
      const a = val();
      if (k === 'T' || k === 'D') return { t: 'u', k, a, b: val() };
      return { t: 'u', k, a };
    }
    const tree = series();
    if (peek() !== undefined) err('Unexpected “' + s[i] + '”');
    return tree;
  }
  const valOf = (x, P) => { if (!x) return NaN; if (x.n !== undefined) return x.n; if (!(x.p in P)) throw new Error('Parameter “' + x.p + '” has no value.'); return P[x.p]; };
  function paramsIn(n, out = new Set()) {
    const add = x => { if (x && x.p) out.add(x.p); };
    if (n.t === 'u') { add(n.a); add(n.b); }
    else if (n.t === 'r') { add(n.R); paramsIn(n.c, out); if (n.back) paramsIn(n.back, out); }
    else { n.c.forEach(c => paramsIn(c, out)); (n.f || []).forEach(add); }
    return out;
  }
  function volume(n, P) {
    if (n.t === 'u') return n.k === 'B' ? 0 : valOf(n.a, P);
    if (n.t === 'r') return volume(n.c, P) + (n.back ? volume(n.back, P) : 0);
    return N.sum(n.c.map(c => volume(c, P)));
  }
  function fracs(n, P) {
    const f = n.f.map(x => (x ? valOf(x, P) : NaN));
    const given = f.filter(isFinite);
    if (given.length === 0) { const V = n.c.map(c => volume(c, P)), tot = N.sum(V); return V.map(x => x / tot); }
    if (given.length < f.length) { const rest = Math.max(0, 1 - N.sum(given)), k = f.length - given.length; return f.map(x => (isFinite(x) ? x : rest / k)); }
    const tot = N.sum(f);
    return f.map(x => x / tot);
  }
  /* time needed to capture nearly all the tracer (for the grid) */
  function span(n, f, v, P) {
    if (n.t === 'u') {
      if (n.k === 'B') return 0;
      const tau = valOf(n.a, P) / (f * v);
      if (n.k === 'P') return tau;
      if (n.k === 'M') return 13 * tau;
      if (n.k === 'T') return tau * (1 + 10 / Math.sqrt(valOf(n.b, P)));
      if (n.k === 'D') { const d = valOf(n.b, P); return tau * (1 + 10 * Math.sqrt(2 * d) + (d > 0.2 ? 10 : 0)); }
      return 30 * tau;
    }
    if (n.t === 's') return N.sum(n.c.map(c => span(c, f, v, P)));
    if (n.t === 'p') { const fr = fracs(n, P); return Math.max(...n.c.map((c, j) => (fr[j] > 0 ? span(c, f * fr[j], v, P) : 0))); }
    const R = valOf(n.R, P), g = span(n.c, f * (1 + R), v, P) + (n.back ? span(n.back, f * R, v, P) : 0);
    if (!(R > 0)) return span(n.c, f, v, P);
    const k = Math.min(400, Math.log(1e-4) / Math.log(R / (1 + R)));
    return g * Math.max(1, k);
  }
  function buildE(n, f, v, P, tmax, ng) {
    if (n.t === 'u') {
      if (n.k === 'B') return RT.atom(0, 1, tmax, ng);
      const tau = valOf(n.a, P) / (f * v);
      if (!(tau > 0)) return RT.atom(0, 1, tmax, ng);
      if (n.k === 'P') return RT.atom(tau, 1, tmax, ng);
      if (n.k === 'M') return RT.sample(t => RT.Emixed(t, tau), tmax, ng);
      if (n.k === 'T') { const k = valOf(n.b, P); return RT.sample(t => RT.Etanks(t, tau, k), tmax, ng); }
      if (n.k === 'D') { const d = valOf(n.b, P); return RT.sampleSmooth(t => RT.Eclosed(t, tau, d), tmax, ng); }
      return RT.sample(t => RT.Elam(t, tau), tmax, ng);
    }
    if (n.t === 's') return n.c.map(c => buildE(c, f, v, P, tmax, ng)).reduce((A, B) => RT.conv(A, B));
    if (n.t === 'p') {
      const fr = fracs(n, P);
      let S = RT.grid(tmax, ng);
      n.c.forEach((c, j) => { if (fr[j] > 0) S = RT.add(S, RT.scale(buildE(c, f * fr[j], v, P, tmax, ng), fr[j])); });
      return S;
    }
    const R = valOf(n.R, P);
    const g = buildE(n.c, f * (1 + R), v, P, tmax, ng);
    if (!(R > 0)) return g;
    const h = n.back ? buildE(n.back, f * R, v, P, tmax, ng) : null;
    const a = 1 / (1 + R), b = R / (1 + R);
    let Pk = g, E = RT.scale(g, a), wk = a;
    const loop = h ? RT.conv(h, g) : g;
    for (let k = 2; k < 600; k++) {
      wk *= b;
      Pk = RT.conv(Pk, loop);
      if (wk * RT.mass(Pk) < 1e-7) break;
      E = RT.add(E, RT.scale(Pk, wk));
    }
    return E;
  }
  /* microfluid conversion through the structure */
  function evalX(n, Xin, f, v, P, m, trace) {
    if (n.t === 'u') {
      if (n.k === 'B') return Xin;
      const tau = valOf(n.a, P) / (f * v);
      let X;
      if (n.k === 'P') X = K.XPFR(m, tau, Xin).X;
      else if (n.k === 'M') X = K.XMFR(m, tau, Xin);
      else if (n.k === 'T') { const k = Math.max(1, Math.round(valOf(n.b, P))); X = Xin; for (let i = 0; i < k; i++) X = K.XMFR(m, tau / k, X); }
      else if (n.k === 'D') X = RT.dispersionX(m, tau, valOf(n.b, P), { Xin }).X;
      else { const Xb = RT.batchX(m, 40 * tau, Xin); X = N.quad(t => Xb(t) * RT.Elam(t, tau), tau / 2, Infinity, { rtol: 1e-9 }); }
      if (trace) trace.push({ n, f, tau, Xin, X });
      return X;
    }
    if (n.t === 's') { let X = Xin; n.c.forEach(c => { X = evalX(c, X, f, v, P, m, trace); }); return X; }
    if (n.t === 'p') { const fr = fracs(n, P); return N.sum(n.c.map((c, j) => (fr[j] > 0 ? fr[j] * evalX(c, Xin, f * fr[j], v, P, m, trace) : 0))); }
    const R = valOf(n.R, P);
    if (!(R > 0)) return evalX(n.c, Xin, f, v, P, m, trace);
    const cap = Math.min(m.Xe, m.Xmax) * (1 - 1e-12);
    const g = X => { const Xmix = (Xin + R * (n.back ? evalX(n.back, X, f * R, v, P, m, null) : X)) / (1 + R); return evalX(n.c, Xmix, f * (1 + R), v, P, m, null) - X; };
    const X = N.root(g, Xin, cap, { n: 60 });
    const Xf = isFinite(X) ? X : Xin;
    if (trace) { const Xmix = (Xin + R * (n.back ? evalX(n.back, Xf, f * R, v, P, m, trace) : Xf)) / (1 + R); evalX(n.c, Xmix, f * (1 + R), v, P, m, trace); }
    return Xf;
  }
  const str = (n, P) => {
    const V = x => (x.p ? x.p + '=' + N.fmt(P[x.p], 3) : N.fmt(x.n, 4));
    if (n.t === 'u') return n.k === 'B' ? 'B' : n.k + ' ' + V(n.a) + (n.b ? ' ' + V(n.b) : '');
    if (n.t === 's') return n.c.map(c => str(c, P)).join(' > ');
    if (n.t === 'p') { const fr = fracs(n, P); return '(' + n.c.map((c, j) => N.fmt(fr[j], 3) + ': ' + str(c, P)).join(' | ') + ')'; }
    return 'R ' + V(n.R) + ' [' + str(n.c, P) + (n.back ? ' ; ' + str(n.back, P) : '') + ']';
  };

  /* schematic */
  function lay(n) {
    if (n.t === 'u') { const w = n.k === 'B' ? 40 : n.k === 'M' ? 46 : n.k === 'T' ? 78 : 76; return { w, h: 60, n }; }
    if (n.t === 's') { const k = n.c.map(lay); return { w: N.sum(k.map(q => q.w)) + 18 * (k.length - 1), h: Math.max(...k.map(q => q.h)), n, kids: k }; }
    if (n.t === 'p') { const k = n.c.map(lay); return { w: Math.max(...k.map(q => q.w)) + 40, h: N.sum(k.map(q => q.h)) + 8 * (k.length - 1), n, kids: k }; }
    const a = lay(n.c), b = n.back ? lay(n.back) : null;
    return { w: Math.max(a.w, b ? b.w : 0) + 44, h: a.h + (b ? b.h : 34) + 6, n, kids: [a].concat(b ? [b] : []) };
  }
  function draw(L, x, ym, out, P, f, v, lab) {
    const n = L.n;
    if (n.t === 'u') {
      const tau = n.k === 'B' ? 0 : valOf(n.a, P) / (f * v);
      if (n.k === 'B') out.push(D.pipe([[x, ym], [x + L.w, ym]], { arrow: false }), D.text(x + L.w / 2, ym - 7, 'bypass', { fs: 9.5 }));
      else if (n.k === 'M') out.push(D.tank(x + 2, ym - 22, L.w - 4, 44, { fill: 0.15, stir: true }));
      else if (n.k === 'T') { const k = Math.min(4, Math.round(valOf(n.b, P))); for (let i = 0; i < k; i++) out.push(D.tank(x + i * (L.w / k) + 1, ym - 12, L.w / k - 3, 24, { fill: 0.15, stir: false })); }
      else out.push(D.tube(x, ym - 10, L.w, 20, { prof: n.k === 'L' ? [0.1, 0.5, 0.1] : null, packed: n.k === 'D' }));
      if (n.k !== 'B') out.push(D.text(x + L.w / 2, ym + 36, (n.k === 'T' ? Math.round(valOf(n.b, P)) + ' tanks ' : n.k === 'D' ? 'D/uL ' + N.fmt(valOf(n.b, P), 2) + ' ' : '') + 'τ ' + N.fmt(tau, 3), { fs: 9.5 }));
      return;
    }
    if (n.t === 's') { let cx = x; L.kids.forEach((k, i) => { draw(k, cx, ym, out, P, f, v, lab); cx += k.w; if (i < L.kids.length - 1) { out.push(D.pipe([[cx, ym], [cx + 16, ym]], { arrow: false })); cx += 18; } }); return; }
    if (n.t === 'p') {
      const fr = fracs(n, P);
      let cy = ym - L.h / 2;
      const xr = x + L.w;
      L.kids.forEach((k, j) => {
        const yk = cy + k.h / 2, x0 = x + 20 + (L.w - 40 - k.w) / 2;
        out.push(D.pipe([[x, ym], [x, yk], [x0, yk]], { arrow: false }), D.text(x + 4, yk - 5, N.fmt(fr[j], 2), { anchor: 'start', fs: 9 }));
        draw(k, x0, yk, out, P, f * fr[j], v, lab);
        out.push(D.pipe([[x0 + k.w, yk], [xr, yk], [xr, ym]], { arrow: false }));
        cy += k.h + 8;
      });
      out.push(D.node(x, ym), D.node(xr, ym));
      return;
    }
    const R = valOf(n.R, P);
    const a = L.kids[0], b = L.kids[1];
    const x0 = x + 22, xr = x + L.w - 2;
    out.push(D.pipe([[x, ym], [x0, ym]], { arrow: false }));
    draw(a, x0 + (L.w - 44 - a.w) / 2, ym, out, P, f * (1 + R), v, lab);
    out.push(D.pipe([[x0 + (L.w - 44 - a.w) / 2 + a.w, ym], [xr, ym]], { arrow: false }));
    const yb = ym + a.h / 2 + (b ? b.h / 2 : 18) + 4;
    if (b) { out.push(D.pipe([[xr - 10, ym], [xr - 10, yb], [x0 + (L.w - 44 - b.w) / 2 + b.w, yb]], { arrow: false })); draw(b, x0 + (L.w - 44 - b.w) / 2, yb, out, P, f * R, v, lab); out.push(D.pipe([[x0 + (L.w - 44 - b.w) / 2, yb], [x + 10, yb], [x + 10, ym + 2]])); }
    else out.push(D.pipe([[xr - 10, ym], [xr - 10, yb], [x + 10, yb], [x + 10, ym + 2]]));
    out.push(D.text((x + xr) / 2, yb - 5, 'R = ' + N.fmt(R, 3), { fs: 9.5 }), D.node(x + 10, ym), D.node(xr - 10, ym));
  }

  const FIG121 = [['Plug flow', 'P 1'], ['Mixed flow', 'M 1'], ['Plug then mixed', 'P 0.4 > M 0.6'], ['Mixed with bypass', '(0.2: B | 0.8: M 1)'], ['Parallel plug paths', '(0.5: P 0.3 | 0.5: P 0.7)'], ['Plug with recycle', 'R 2 [ P 1 ]'], ['Two mixed in parallel', '(0.5: M 0.2 | 0.5: M 0.8)'], ['Mixed with plug side branch', '(0.7: M 0.4 | 0.3: P 0.6)']];

  ST.tool({
    id: 'compart', tag: 'F-1201', ch: 12, sec: 'Ch. 12', title: 'Compartment models: build, fit to a tracer curve, and react',
    lede: 'Describe a vessel as plug, mixed, tanks-in-series, dispersed or laminar regions, joined in series, in parallel with any flow split, with bypass, recycle loops (with or without a return line) and dead water. See its E and F curves and moments, fit any of its parameters to a measured tracer curve, and run any reaction through it as a microfluid or a macrofluid — then compare with a proposed fix.',
    keys: 'compartment model dead volume bypass recycle short circuit channeling plug mixed combination fit tracer curve diagnose misbehaving reactor stagnant example 12.2 example 14.3 flow model',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'Flow model' },
      { k: 'model', t: 'text', label: 'Model A (present vessel)', v: 'M 2', mono: true, hint: 'P V plug · M V mixed · T V n tanks · D V d dispersion · L V laminar · B bypass · “>” series · ( f: … | f: … ) parallel with flow fractions · R r [ forward ; return ] recycle. Volumes may be parameter names.' },
      { k: 'Vd', label: 'Dead (stagnant) volume, model A', unit: 'vol', v: 4, min: 0, max: 1e9 },
      { k: 'v', label: 'Flow rate v', unit: 'vol/[t]', v: 1, min: 1e-12, max: 1e12, scale: 'log' },
      { k: 'params', t: 'text', label: 'Parameter values', v: '', mono: true, hint: 'e.g. a = 0.4, b = 0.6, r = 1' },
      { k: 'cmpOn', t: 'check', label: 'Compare with a model B (proposed change)', v: true },
      { k: 'modelB', t: 'text', label: 'Model B', v: 'M 6', mono: true, if: v => v.cmpOn },
      { k: 'VdB', label: 'Dead volume, model B', unit: 'vol', v: 0, min: 0, max: 1e9, if: v => v.cmpOn },
      { group: 'Fit to a tracer curve' },
      { k: 'fitOn', t: 'check', label: 'Fit parameters of model A to pulse data', v: false },
      { k: 'data', t: 'table', label: 'Pulse response (t, C)', cols: ['t', 'C'], v: [], if: v => v.fitOn },
      { k: 'fitp', t: 'text', label: 'Parameters to fit', v: 'a, b', mono: true, if: v => v.fitOn },
      { k: 'fitOnF', t: 'seg', label: 'Match', opts: [['E', 'E curve'], ['F', 'F curve']], v: 'F', if: v => v.fitOn, hint: 'F is safer when the model has plug-flow parts.' },
      { group: 'Reaction' },
      { k: 'rxOn', t: 'check', label: 'React in the model', v: true },
      { k: 'calib', t: 'check', label: 'Choose k so that model A gives a measured conversion', v: true, if: v => v.rxOn },
      { k: 'Xmeas', label: 'Measured conversion in model A', v: 0.75, min: 1e-6, max: 0.999999, if: v => v.rxOn && v.calib }
    ].concat(KIT.kin({ gas: false, only: ['power', 'rev', 'sat', 'auto', 'lh', 'inhib', 'custom'], d: { law: 'power', k: 1, a: 1, CA0: 1 } }).map(f => (f.k || f.group ? Object.assign({}, f, { if: f.if ? (vv => vv.rxOn && f.if(vv)) : (f.group ? (vv => vv.rxOn) : (vv => vv.rxOn)) }) : f))),
    presets: [
      { name: 'Example 12.2', v: { model: 'M 2', Vd: 4, v: 1, params: '', cmpOn: true, modelB: 'M 6', VdB: 0, fitOn: false, rxOn: true, calib: true, Xmeas: 0.75, law: 'power', a: 1, ob: 0, k: 1, CA0: 1 }, book: { XB: 0.9 }, note: 'A 6 m³ tank: only 2 m³ is mixed, 4 m³ is dead water, and conversion is 75% (first order). A proper stirrer makes all 6 m³ active: 90%.' },
      { name: 'Example 14.3 (recycle with 12 + 24 tanks)', v: { model: 'R 1 [ T 1/3 12 ; T 2/3 24 ]', Vd: 0, v: 1, params: '', cmpOn: false, fitOn: false, rxOn: false }, book: { tbar: 1 }, note: 'Peaks at θ = 1/6, 1 and 11/6 with areas halving: forward region 1/3 as 12 tanks, return line 2/3 as 24 tanks, R = 1.' },
      { name: 'Fit: plug + mixed to Example 11.1', v: { model: 'P a > M b', Vd: 0, v: 1, params: 'a = 5, b = 10', cmpOn: false, fitOn: true, data: [[0, 0], [5, 3], [10, 5], [15, 5], [20, 4], [25, 2], [30, 1], [35, 0]], fitp: 'a, b', fitOnF: 'F', rxOn: true, calib: false, law: 'power', a: 1, ob: 0, k: 0.307, CA0: 1, _tu: 'min' }, note: 'A delay plus a mixed region fitted to the tracer curve of Example 11.1, then the first-order reaction of Example 11.4 run through it.' },
      { name: 'Problem 12.11 (second order)', v: { model: '(0.75: M 3 | 0.25: B)', Vd: 3, v: 1, params: '', cmpOn: true, modelB: 'M 6', VdB: 0, rxOn: true, calib: true, Xmeas: 0.6, law: 'power', a: 2, ob: 0, k: 1, CA0: 1 }, note: 'A + B, equimolar, 60% conversion in a 6 m³ tank that behaves as part mixed, part bypassed, part dead (structure as an example; read yours from Fig. P12.11). Model B: what the fully mixed tank would do.' }
    ].concat(FIG121.map(([nm, md]) => ({ name: 'Fig. 12.1: ' + nm, v: { model: md, Vd: 0, v: 1, params: '', cmpOn: false, fitOn: false, rxOn: false }, note: 'One of the elementary combinations of Fig. 12.1, with unit flow and unit total volume.' }))),
    eqs: [
      ['V = V_p + V_m + V_d,\\qquad v = v_a + v_b + v_r', 'Ch. 12'],
      ['\\bar t_{obs} = \\frac{V_{active}}{v} < \\frac{V}{v}\\ \\text{if dead spaces}', 'Ch. 12 (b)'],
      ['E_{series} = E_1 * E_2,\\qquad E_{parallel} = \\sum_j f_jE_j', 'Eqs. 11.10–11.12'],
      ['E_{recycle} = \\sum_{k\\ge 1} \\frac{1}{1+R}\\left(\\frac{R}{1+R}\\right)^{k-1} g*(h*g)^{*(k-1)}', '', 'g: forward region at flow (1+R)v; h: return line at flow Rv.'],
      ['\\frac{C_A}{C_{A0}} = \\int_0^\\infty e^{-kt}E\\,dt', 'Eq. 11.13', 'First order: exact for any structure.']
    ],
    method: ['The model is parsed into a tree and its E curve assembled on a time grid that keeps plug-flow parts as exact Dirac atoms: series blocks are convolved by FFT, parallel blocks weighted by their flow fractions, and recycle loops summed pass by pass until less than 10⁻⁷ of the tracer remains.', 'Fitting adjusts the chosen parameters by Levenberg–Marquardt (log-transformed, so they stay positive) to match the normalized E or F curve of the data.', 'Microfluid conversion is solved through the structure itself, unit by unit, with recycle loops closed by root-finding; the macrofluid result integrates the batch history over the whole E. k can be calibrated so that model A reproduces a measured conversion, then carried to model B.'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const P = v.params && String(v.params).trim() ? K.parseParams(v.params) : {};
      const flow = +v.v;
      const tu = v._tu || 'min';
      const A = parse(v.model);
      const need = Array.from(paramsIn(A)).filter(p => !(p in P));
      if (need.length) throw new Error('Give values for: ' + need.join(', ') + '.');
      /* ---- fitting ---- */
      let fitInfo = null;
      if (v.fitOn) {
        const rows = RTDI.rows(v.data, 'tracer readings');
        const R = RTDI.pulse(rows, { rule: 'pl' });
        const names = String(v.fitp || '').split(/[,\s]+/).filter(Boolean);
        names.forEach(nm => { if (!(nm in P)) throw new Error('Fit parameter “' + nm + '” has no starting value.'); });
        const tdat = R.t, Edat = R.E, Fdat = N.cumtrapz(R.t, R.E);
        const resid = p => {
          const PP = Object.assign({}, P);
          names.forEach((nm, j) => { PP[nm] = p[j]; });
          const tm = Math.max(1.2 * tdat[tdat.length - 1], 1.1 * span(A, 1, flow, PP));
          const E = buildE(A, 1, flow, PP, tm, 900);
          if (v.fitOnF === 'F') { const F = RT.cdf(E); return tdat.map((t, i) => N.interp(N.range(E.n).map(q => q * E.dt), Array.from(F), t) - Fdat[i]); }
          const sc = Math.max(...Edat);
          return tdat.map((t, i) => (RT.at(E, t) - Edat[i]) / sc);
        };
        const fit = N.lm(resid, names.map(nm => P[nm]), { pos: names.map(() => true), maxIt: 80 });
        names.forEach((nm, j) => { P[nm] = fit.p[j]; });
        const r = resid(fit.p), rms = Math.sqrt(N.sum(r.map(q => q * q)) / r.length);
        fitInfo = { names, p: fit.p, rms, R };
        names.forEach((nm, j) => kpi.push({ k: 'fit_' + nm, l: 'Fitted ' + nm, v: fit.p[j], main: true }));
        kpi.push({ k: 'rms', l: 'RMS misfit (' + v.fitOnF + ')', v: rms });
      }
      /* ---- E of model A ---- */
      const VA = volume(A, P), VdA = +v.Vd;
      const tmA = 1.1 * Math.max(span(A, 1, flow, P), fitInfo ? fitInfo.R.t[fitInfo.R.t.length - 1] : 0, 1e-9);
      const EA = buildE(A, 1, flow, P, tmA, 3000);
      const moA = RT.moments(EA);
      kpi.push({ k: 'tbar', l: 't̄ from the E curve, model A', v: moA.mean, u: tu, main: true, s: 'V_active/v = ' + N.fmt(VA / flow, 4) });
      kpi.push({ k: 'tV', l: 'V/v including dead water', v: (VA + VdA) / flow, u: tu, s: VdA > 0 ? 't̄_obs < V/v: ' + N.fmt(100 * VdA / (VA + VdA), 3) + '% dead' : 'no dead water' });
      kpi.push({ k: 's2t', l: 'σ²_θ, model A', v: moA.var / (moA.mean * moA.mean), s: 'N = ' + N.fmt(moA.mean * moA.mean / moA.var, 3) + ' equivalent tanks' });
      if (Math.abs(moA.area - 1) > 0.01) warn.push('Only ' + N.fmt(100 * moA.area, 4) + '% of the tracer leaves within the time window: the tail is very long.');
      /* ---- model B ---- */
      let B = null, EB = null, moB = null;
      if (v.cmpOn && String(v.modelB || '').trim()) {
        B = parse(v.modelB);
        const needB = Array.from(paramsIn(B)).filter(p => !(p in P));
        if (needB.length) throw new Error('Model B needs values for: ' + needB.join(', ') + '.');
        EB = buildE(B, 1, flow, P, 1.1 * span(B, 1, flow, P), 3000);
        moB = RT.moments(EB);
        kpi.push({ k: 'tbarB', l: 't̄, model B', v: moB.mean, u: tu, s: 'σ²_θ = ' + N.fmt(moB.var / (moB.mean * moB.mean), 4) });
        if (Math.abs((volume(B, P) + +v.VdB) - (VA + VdA)) > 1e-6 * (VA + VdA)) info.push('Models A and B have different total volumes (' + N.fmt(VA + VdA, 4) + ' and ' + N.fmt(volume(B, P) + +v.VdB, 4) + ').');
      }
      /* ---- reaction ---- */
      const tables = [];
      if (v.rxOn) {
        let vv = v;
        if (v.calib) {
          if (v.law !== 'power') throw new Error('Calibration scales k of a power-law rate; pick that form or untick calibration.');
          const Xt = +v.Xmeas;
          const g = lk => { const m = KIT.model(Object.assign({}, v, { k: Math.exp(lk) })); return evalX(A, 0, 1, flow, P, m, null) - Xt; };
          const lk = N.root(g, Math.log(+v.k) - 30, Math.log(+v.k) + 30, { n: 120 });
          if (!isFinite(lk)) throw new Error('No rate constant makes model A reach X = ' + N.fmt(Xt, 4) + '.');
          vv = Object.assign({}, v, { k: Math.exp(lk) });
          kpi.push({ k: 'kcal', l: 'k that reproduces X_A = ' + N.fmt(Xt, 4) + ' in model A', v: Math.exp(lk), u: KIT.kUnitFor(v), main: true });
        }
        const m = KIT.model(vv);
        const trace = [];
        const XA = evalX(A, 0, 1, flow, P, m, trace);
        const XbA = RT.batchX(m, tmA);
        const XAseg = RT.integrate(EA, t => XbA(t)) / moA.area;
        kpi.push({ k: 'XA', l: 'X_A in model A, microfluid', v: XA, main: !v.calib }, { k: 'XAseg', l: 'X_A in model A, macrofluid', v: XAseg });
        const rows = [['A: ' + str(A, P), XA, XAseg]];
        if (B) {
          const XB = evalX(B, 0, 1, flow, P, m, null);
          const XbB = RT.batchX(m, EB.tmax);
          const XBseg = RT.integrate(EB, t => XbB(t)) / moB.area;
          kpi.push({ k: 'XB', l: 'X_A in model B, microfluid', v: XB, main: true, s: 'macrofluid ' + N.fmt(XBseg, 4) });
          rows.push(['B: ' + str(B, P), XB, XBseg]);
        }
        tables.push({ title: 'Conversion', cols: ['model', 'microfluid (through the structure)', 'macrofluid (segregated)'], rows });
        tables.push({ title: 'Units of model A', cols: ['unit', 'share of flow', 'τ', 'X_A in', 'X_A out'], rows: trace.map(q => [str(q.n, P), q.f, q.tau, q.Xin, q.X]) });
      }
      /* ---- plots ---- */
      const plots = [];
      const tmP = Math.min(tmA, moA.mean * 4 + Math.sqrt(moA.var) * 3);
      const step = Math.max(1, Math.floor(3000 / 500));
      const tx = [], ex = [];
      for (let i = 0; i < EA.n; i += step) { if (i * EA.dt > tmP) break; tx.push(i * EA.dt); ex.push(EA.c[i]); }
      const ser = [{ name: 'model A', x: tx, y: ex, color: 1 }];
      if (EA.d.length) ser.push({ name: 'model A: plug-flow pulses (area)', x: EA.d.map(a => a[0]), y: EA.d.map(a => a[1] / Math.max(EA.dt * 20, moA.mean * 0.02)), type: 'bar', color: 1, bw: 3, legend: true });
      if (EB) { const tB = [], eB = []; for (let i = 0; i < EB.n; i += step) { if (i * EB.dt > tmP) break; tB.push(i * EB.dt); eB.push(EB.c[i]); } ser.push({ name: 'model B', x: tB, y: eB, color: 2, dash: '5 4' }); }
      if (fitInfo) ser.push({ name: 'data (normalized)', x: fitInfo.R.t, y: fitInfo.R.E, type: 'pts', color: 3 });
      plots.push({ title: 'E curves', x: { label: 't', unit: tu, min: 0, max: tmP }, y: { label: 'E', unit: '1/' + tu, min: 0 }, series: ser, cap: EA.d.length ? 'Bars mark plug-flow pulses; their height is drawn proportional to the tracer they carry.' : undefined });
      const FA = RT.cdf(EA), tf = [], ff = [];
      for (let i = 0; i < EA.n; i += step) { if (i * EA.dt > tmP) break; tf.push(i * EA.dt); ff.push(FA[i]); }
      const fser = [{ name: 'model A', x: tf, y: ff, color: 1 }];
      if (EB) { const FB = RT.cdf(EB), tB = [], fB = []; for (let i = 0; i < EB.n; i += step) { if (i * EB.dt > tmP) break; tB.push(i * EB.dt); fB.push(FB[i]); } fser.push({ name: 'model B', x: tB, y: fB, color: 2, dash: '5 4' }); }
      if (fitInfo) fser.push({ name: 'data', x: fitInfo.R.t, y: N.cumtrapz(fitInfo.R.t, fitInfo.R.E), type: 'pts', color: 3 });
      plots.push({ title: 'F curves', x: { label: 't', unit: tu, min: 0, max: tmP }, y: { label: 'F', min: 0, max: 1 }, series: fser });
      /* diagram */
      const L = lay(A), out = [];
      const H = Math.max(90, L.h + 52), ym = 16 + L.h / 2;
      out.push(D.pipe([[4, ym], [28, ym]], { label: 'v = ' + N.fmt(flow, 3), dy: -8 }));
      draw(L, 30, ym, out, P, 1, flow);
      out.push(D.pipe([[30 + L.w, ym], [30 + L.w + 34, ym]]));
      let W = L.w + 80;
      if (VdA > 0) { out.push('<rect class="pack" x="' + (W - 2) + '" y="' + (ym - 18) + '" width="46" height="36"/>', D.text(W + 21, ym + 4, 'dead', { fs: 9.5 }), D.text(W + 21, ym + 30, 'V ' + N.fmt(VdA, 3), { fs: 9.5 })); W += 56; }
      return { kpi, warn, info, plots, tables, diagram: D.svg(W + 10, H, out.join(''), 'compartment model of the vessel') };
    }
  });

  /* ===================================================================== */
  /* F-1202  phase volumes from tracer tests                               */
  /* ===================================================================== */
  ST.tool({
    id: 'holdup', tag: 'F-1202', ch: 12, sec: 'Ch. 12', title: 'Phase volumes from tracer tests: gas, flowing liquid and stagnant liquid',
    lede: 'One tracer per phase, the flow of each phase, and the vessel’s volume: the mean times give the volume each flowing phase occupies, and whatever is left over is stagnant. Enter each response as a table, or just its mean time.',
    keys: 'holdup gas liquid contactor phase fraction stagnant volume tracer mean residence time example 12.1 bubble column packed tower voidage',
    inputs: [
      { group: 'Vessel' },
      { k: 'vmode', t: 'seg', label: 'Volume available to the fluids', opts: [['geo', 'height × area × voidage'], ['V', 'enter it']], v: 'geo' },
      { k: 'H', label: 'Height', unit: 'm', v: 20, min: 0, max: 1e6, if: v => v.vmode === 'geo' },
      { k: 'A', label: 'Cross-section', unit: 'm²', v: 1, min: 0, max: 1e6, if: v => v.vmode === 'geo' },
      { k: 'eps', label: 'Voidage', v: 0.5, min: 0, max: 1, if: v => v.vmode === 'geo' },
      { k: 'V', label: 'Volume', unit: 'm³', v: 10, min: 0, max: 1e9, if: v => v.vmode === 'V' },
      { group: 'Phase 1' },
      { k: 'n1', t: 'text', label: 'Name', v: 'gas' },
      { k: 'q1', label: 'Flow rate', unit: 'm³/s', v: 0.5, min: 0, max: 1e9 },
      { k: 'm1', t: 'seg', label: 'Response given as', opts: [['tab', 'a table'], ['t', 'its mean time']], v: 'tab' },
      { k: 'd1', t: 'table', label: 'Pulse response (t in s, C)', cols: ['t', 'C'], v: [[6, 0], [9, 1], [15, 0]], if: v => v.m1 === 'tab' },
      { k: 't1', label: 'Mean time', unit: 's', v: 10, min: 0, max: 1e12, if: v => v.m1 === 't' },
      { group: 'Phase 2' },
      { k: 'n2', t: 'text', label: 'Name', v: 'liquid' },
      { k: 'q2', label: 'Flow rate', unit: 'm³/s', v: 0.1, min: 0, max: 1e9 },
      { k: 'm2', t: 'seg', label: 'Response given as', opts: [['tab', 'a table'], ['t', 'its mean time']], v: 't' },
      { k: 'd2', t: 'table', label: 'Pulse response (t in s, C)', cols: ['t', 'C'], v: [[30, 0], [40, 1], [50, 0]], if: v => v.m2 === 'tab' },
      { k: 't2', label: 'Mean time', unit: 's', v: 40, min: 0, max: 1e12, if: v => v.m2 === 't' }
    ],
    presets: [
      { name: 'Example 12.1', v: { vmode: 'geo', H: 20, A: 1, eps: 0.5, n1: 'gas', q1: 0.5, m1: 'tab', d1: [[6, 0], [9, 1], [15, 0]], n2: 'liquid', q2: 0.1, m2: 't', t2: 40 }, book: { t1: 10, V1: 5, V2: 4, fs: 0.1 }, note: 'Packed column 20 m tall, 1 m² cross-section, voidage 0.5. Gas at 0.5 m³/s gives a triangular response from 6 to 15 s peaking at 9 s; liquid at 0.1 m³/s has t̄ = 40 s. Gas 50%, flowing liquid 40%, stagnant 10%.' }
    ],
    eqs: [['V_i = \\bar t_i\\,v_i,\\qquad V_{stag} = V - \\sum_i V_i', 'Ex. 12.1']],
    method: ['Each mean time is the first moment of its piecewise-linear response (or taken as given). The stagnant volume is whatever the flowing phases do not account for; a negative value means a flow or a mean time is wrong.'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const V = v.vmode === 'geo' ? +v.H * +v.A * +v.eps : +v.V;
      const tOf = (mode, tab, tt) => (mode === 't' ? +tt : RTDI.pulse(RTDI.rows(tab, 'readings'), { rule: 'pl' }).mean);
      const t1 = tOf(v.m1, v.d1, v.t1), t2 = tOf(v.m2, v.d2, v.t2);
      const V1 = t1 * +v.q1, V2 = t2 * +v.q2, Vs = V - V1 - V2;
      kpi.push({ k: 'V', l: 'Volume available', v: V, u: 'm³' });
      kpi.push({ k: 't1', l: 't̄, ' + v.n1, v: t1, u: 's' }, { k: 'V1', l: 'Volume of ' + v.n1, v: V1, u: 'm³', main: true, s: N.fmt(100 * V1 / V, 3) + '%' });
      kpi.push({ k: 't2', l: 't̄, ' + v.n2, v: t2, u: 's' }, { k: 'V2', l: 'Volume of flowing ' + v.n2, v: V2, u: 'm³', main: true, s: N.fmt(100 * V2 / V, 3) + '%' });
      kpi.push({ k: 'fs', l: 'Stagnant fraction', v: Vs / V, main: true, s: N.fmt(Vs, 4) + ' m³', state: Vs < -1e-9 * V ? 'bad' : undefined });
      if (Vs < -1e-9 * V) warn.push('The flowing phases occupy more than the vessel: a flow rate, mean time or volume is wrong.');
      const plots = [{ title: 'Where the volume is', x: { label: '' , min: 0, max: 4 }, y: { label: 'share of the volume', min: 0, max: 1 }, series: [{ name: 'share', x: [1, 2, 3], y: [V1 / V, V2 / V, Math.max(0, Vs / V)], type: 'bar', color: 1, bw: 40 }], marks: [{ t: 'text', x: 1, y: V1 / V + 0.04, s: v.n1, anchor: 'middle' }, { t: 'text', x: 2, y: V2 / V + 0.04, s: v.n2, anchor: 'middle' }, { t: 'text', x: 3, y: Math.max(0, Vs / V) + 0.04, s: 'stagnant', anchor: 'middle' }], legend: false }];
      return { kpi, warn, info, plots };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
