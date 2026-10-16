/* Reaction networks: parser, mass-action rates, ideal reactor solvers. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N, X = ST.X;
  const R_GAS = 8.314462618;

  function parseSide(side, line) {
    const terms = side.split('+').map(s => s.trim()).filter(Boolean);
    return terms.map(t => {
      const m = /^(\d*\.?\d*)\s*([A-Za-z][A-Za-z0-9_]*)$/.exec(t);
      if (!m) throw new Error('Line ' + line + ': cannot read "' + t + '"');
      const nu = m[1] === '' ? 1 : parseFloat(m[1]);
      if (!(nu > 0)) throw new Error('Line ' + line + ': bad coefficient in "' + t + '"');
      return [m[2], nu];
    });
  }

  /* Parse text such as
       A + B -> R ; k = 0.5
       A <=> S    ; k = 2, kr = 0.4
       2A -> D    ; k = 1, r = k*A^2/(1+K*A)
       K = 3                      (a global parameter)                          */
  ST.parseNet = function (text, o = {}) {
    const lines = String(text).split(/\n|;;/);
    const species = [], rx = [], params = {};
    const addSp = s => { if (species.indexOf(s) < 0) species.push(s); return species.indexOf(s); };
    (o.species || []).forEach(addSp);
    const pending = [];
    lines.forEach((raw, li) => {
      const line = raw.replace(/(#|\/\/).*$/, '').trim();
      if (!line) return;
      const ln = li + 1;
      const arrow = /<=>|<->|⇌|->|→|=>/.exec(line);
      if (!arrow) {
        const m = /^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$/.exec(line);
        if (!m) throw new Error('Line ' + ln + ': expected a reaction like "A + B -> R ; k = 1"');
        const v = parseFloat(m[2]);
        if (!isFinite(v)) throw new Error('Line ' + ln + ': "' + m[2] + '" is not a number');
        params[m[1]] = v;
        return;
      }
      const semi = line.indexOf(';');
      const eqn = semi >= 0 ? line.slice(0, semi) : line, rest = semi >= 0 ? line.slice(semi + 1) : '';
      const rev = /<=>|<->|⇌/.test(arrow[0]);
      const [lhs, rhs] = eqn.split(arrow[0]);
      const reac = lhs.trim() ? parseSide(lhs, ln) : [], prod = rhs && rhs.trim() ? parseSide(rhs, ln) : [];
      if (!reac.length && !prod.length) throw new Error('Line ' + ln + ': empty reaction');
      const loc = {};
      let custom = null;
      rest.split(',').map(s => s.trim()).filter(Boolean).forEach(kv => {
        const m = /^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$/.exec(kv);
        if (!m) throw new Error('Line ' + ln + ': cannot read "' + kv + '"');
        if (m[1] === 'r') custom = m[2];
        else {
          const v = parseFloat(m[2]);
          if (!isFinite(v)) throw new Error('Line ' + ln + ': "' + m[2] + '" is not a number');
          loc[m[1]] = v;
        }
      });
      reac.forEach(([s]) => addSp(s)); prod.forEach(([s]) => addSp(s));
      pending.push({ ln, reac, prod, rev, loc, custom, text: eqn.trim() });
    });
    if (!pending.length) throw new Error('No reactions found');
    pending.forEach(p => {
      const L = p.loc;
      const k0 = L.k0 !== undefined ? L.k0 : L.A;
      const k = L.k !== undefined ? L.k : (L.kf !== undefined ? L.kf : (L.k1 !== undefined ? L.k1 : undefined));
      if (k === undefined && k0 === undefined && !p.custom) throw new Error('Line ' + p.ln + ': give a rate constant, e.g. "; k = 0.5" or "; k0 = 1e9, ER = 6000"');
      let kr = L.kr !== undefined ? L.kr : (L.kb !== undefined ? L.kb : L.k2);
      const kr0 = L.kr0 !== undefined ? L.kr0 : L.Ar;
      if (p.rev && kr === undefined && L.K !== undefined && k !== undefined) kr = k / L.K;
      if (p.rev && kr === undefined && kr0 === undefined && !p.custom) throw new Error('Line ' + p.ln + ': reversible step needs kr (or K, or kr0 with ERr)');
      const r = {
        text: p.text, rev: p.rev, k: k !== undefined ? k : 0, kr: kr || 0, E: L.E || 0, Er: L.Er || 0, dH: L.dH || 0,
        reac: p.reac.map(([s, nu]) => [species.indexOf(s), nu]),
        prod: p.prod.map(([s, nu]) => [species.indexOf(s), nu]), loc: L
      };
      /* power-law orders that differ from the stoichiometry: n_A = 0, n_B = 1.5 … */
      Object.keys(L).forEach(key => {
        const mo = /^n_([A-Za-z][A-Za-z0-9_]*)$/.exec(key);
        if (!mo) return;
        const i = species.indexOf(mo[1]);
        if (i < 0 || !r.reac.some(([q]) => q === i)) throw new Error('Line ' + p.ln + ': ' + key + ' refers to a species that is not a reactant of this step.');
        (r.ord = r.ord || {})[i] = L[key];
      });
      /* absolute Arrhenius form: k = k0·exp(−ER/T), ER = E/R in kelvin (or give E in J/mol) */
      if (k0 !== undefined) { r.k0 = k0; r.ER = L.ER !== undefined ? L.ER : (L.E || 0) / R_GAS; }
      if (kr0 !== undefined) { r.kr0 = kr0; r.ERr = L.ERr !== undefined ? L.ERr : (L.Er || 0) / R_GAS; }
      if (p.custom) {
        const vars = species.concat(Object.keys(p.loc), Object.keys(params), ['T']);
        r.fn = X.compile(p.custom, vars);
        r.custom = p.custom;
        r.vars = vars;
      }
      rx.push(r);
    });
    const ns = species.length, nr = rx.length;
    const S = N.range(ns).map(() => new Array(nr).fill(0));
    rx.forEach((r, j) => {
      r.reac.forEach(([i, nu]) => { S[i][j] -= nu; });
      r.prod.forEach(([i, nu]) => { S[i][j] += nu; });
    });
    const Tref = o.Tref || params.Tref || 298.15;
    const net = { species, rx, S, params, Tref, ns, nr };
    net.rates = function (C, T) {
      const out = new Array(nr);
      const th = T ? (1 / T - 1 / Tref) : 0;
      for (let j = 0; j < nr; j++) {
        const r = rx[j];
        if (r.fn) {
          const scope = C.map(v => Math.max(v, 0)).concat(Object.values(r.loc), Object.values(params), [T || Tref]);
          let v = r.fn(scope);
          if (r.E && T) v *= Math.exp(-r.E / R_GAS * th);
          out[j] = v;
          continue;
        }
        let f = r.k0 !== undefined ? r.k0 * Math.exp(-r.ER / (T || Tref)) : r.k * (r.E && T ? Math.exp(-r.E / R_GAS * th) : 1);
        for (const [i, nu] of r.reac) { const e = r.ord && r.ord[i] !== undefined ? r.ord[i] : nu; f *= e === 1 ? Math.max(C[i], 0) : e === 0 ? 1 : Math.pow(Math.max(C[i], 0), e); }
        if (r.rev) {
          let b = r.kr0 !== undefined ? r.kr0 * Math.exp(-r.ERr / (T || Tref)) : r.kr * (r.Er && T ? Math.exp(-r.Er / R_GAS * th) : 1);
          for (const [i, nu] of r.prod) b *= nu === 1 ? Math.max(C[i], 0) : Math.pow(Math.max(C[i], 0), nu);
          f -= b;
        }
        out[j] = f;
      }
      return out;
    };
    net.prod = function (C, T) {
      const r = net.rates(C, T), out = new Array(ns).fill(0);
      for (let i = 0; i < ns; i++) for (let j = 0; j < nr; j++) if (S[i][j]) out[i] += S[i][j] * r[j];
      return out;
    };
    net.idx = s => species.indexOf(s);
    return net;
  };

  /* ---------- ideal reactors for a network (constant density) ---------- */
  const R = ST.NR = {};

  /* batch or plug flow: returns {t, C[][]} sampled at n points, plus sol */
  R.batch = function (net, C0, tEnd, o = {}) {
    const f = (t, C) => net.prod(C, o.T);
    const sol = N.ode(f, 0, C0, tEnd, { rtol: 1e-9, atol: 1e-12, event: o.event, minPts: 80 });
    const n = o.n || 160, tEndAct = sol.t[sol.t.length - 1];
    const ts = N.linspace(0, tEndAct, n);
    return { t: ts, C: ts.map(sol.at), sol, end: sol.last(), tEnd: tEndAct, ev: sol.ev };
  };

  /* steady MFR: C0 - C + tau * prod(C) = 0. Continuation in tau from 0. */
  R.mfr = function (net, C0, tau, o = {}) {
    if (tau <= 0) return { C: C0.slice(), ok: true };
    const ns = net.ns, lo = new Array(ns).fill(0);
    const F = t => C => net.prod(C, o.T).map((p, i) => C0[i] - C[i] + t * p);
    if (o.guess) {
      const r0 = N.newton(F(tau), o.guess.slice(), { lo, tol: 1e-12 });
      if (r0.ok && r0.x.every(x => x >= -1e-12)) return { C: r0.x, ok: true };
    }
    let C = C0.slice(), ok = false;
    const steps = o.steps || 24;
    const taus = N.logspace(tau * 1e-4, tau, steps);
    for (const t of taus) {
      const r = N.newton(F(t), C, { lo, tol: 1e-12 });
      if (r.ok) { C = r.x; ok = true; } else { ok = false; }
    }
    if (!ok) {
      const g = (t, c) => net.prod(c, o.T).map((p, i) => (C0[i] - c[i]) / tau + p);
      const s = N.ode(g, 0, C0, 60 * tau, { rtol: 1e-9, atol: 1e-12 });
      C = s.last();
      const r = N.newton(F(tau), C, { lo, tol: 1e-12 });
      if (r.ok) { C = r.x; ok = true; }
    }
    return { C, ok };
  };

  /* chain of reactors: units = [{type:'pfr'|'mfr', tau}] */
  R.train = function (net, C0, units, o = {}) {
    let C = C0.slice();
    const outs = [C.slice()], profiles = [];
    units.forEach(u => {
      if (u.type === 'mfr') {
        C = R.mfr(net, C, u.tau, o).C;
        profiles.push(null);
      } else {
        const b = R.batch(net, C, u.tau, { n: 60, T: o.T });
        profiles.push(b);
        C = b.end.slice();
      }
      outs.push(C.slice());
    });
    return { C, outs, profiles };
  };

  /* plug flow with recycle ratio Rr, total space time tau */
  R.recycle = function (net, C0, tau, Rr, o = {}) {
    if (Rr <= 0) return { C: R.batch(net, C0, tau, { T: o.T }).end };
    const tp = tau / (1 + Rr);
    const map = Cout => {
      const Cin = C0.map((c, i) => (c + Rr * Cout[i]) / (1 + Rr));
      return R.batch(net, Cin, tp, { n: 8, T: o.T }).end;
    };
    const guess = R.mfr(net, C0, tau, o).C;
    const res = N.newton(Co => map(Co).map((v, i) => v - Co[i]), guess, { lo: new Array(net.ns).fill(0), tol: 1e-11 });
    return { C: res.x, ok: res.ok, Cin: C0.map((c, i) => (c + Rr * res.x[i]) / (1 + Rr)) };
  };

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
