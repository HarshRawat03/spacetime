/* Small arithmetic expression compiler. No eval: parses to a tree of closures. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N;
  const X = ST.X = {};

  const FN = {
    exp: Math.exp, ln: Math.log, log: Math.log, log10: Math.log10, sqrt: Math.sqrt, abs: Math.abs,
    sin: Math.sin, cos: Math.cos, tan: Math.tan, sinh: Math.sinh, cosh: Math.cosh, tanh: Math.tanh,
    asin: Math.asin, acos: Math.acos, atan: Math.atan, floor: Math.floor, ceil: Math.ceil, sign: Math.sign,
    min: Math.min, max: Math.max, pow: Math.pow,
    erf: x => N.erf(x), erfc: x => N.erfc(x), step: x => (x >= 0 ? 1 : 0), pos: x => Math.max(0, x)
  };
  const CONST = { pi: Math.PI, Rg: 8.314 };

  function tokenize(src) {
    const s = src.replace(/−/g, '-').replace(/·|×/g, '*').replace(/\*\*/g, '^');
    const out = [];
    let i = 0;
    while (i < s.length) {
      const c = s[i];
      if (/\s/.test(c)) { i++; continue; }
      const num = /^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?/.exec(s.slice(i));
      if (num) { out.push({ t: 'n', v: parseFloat(num[0]), p: i }); i += num[0].length; continue; }
      const id = /^[A-Za-z_Ͱ-Ͽ][A-Za-z0-9_Ͱ-Ͽ]*/.exec(s.slice(i));
      if (id) { out.push({ t: 'id', v: id[0], p: i }); i += id[0].length; continue; }
      if ('+-*/^(),'.indexOf(c) >= 0) { out.push({ t: c, p: i }); i++; continue; }
      throw new Error('Unexpected character "' + c + '" at position ' + (i + 1));
    }
    out.push({ t: 'end', p: s.length });
    return out;
  }

  /* compile(src, vars) -> f(arr) where arr holds values in the order of vars */
  X.compile = function (src, vars) {
    if (!src || !String(src).trim()) throw new Error('Expression is empty');
    const tk = tokenize(String(src));
    let k = 0;
    const used = new Set();
    const peek = () => tk[k], next = () => tk[k++];
    const expect = t => { const x = next(); if (x.t !== t) throw new Error('Expected "' + t + '" at position ' + (x.p + 1)); return x; };
    function expr() {
      let a = term();
      while (peek().t === '+' || peek().t === '-') {
        const op = next().t, b = term(), l = a;
        a = op === '+' ? s => l(s) + b(s) : s => l(s) - b(s);
      }
      return a;
    }
    function term() {
      let a = unary();
      while (peek().t === '*' || peek().t === '/') {
        const op = next().t, b = unary(), l = a;
        a = op === '*' ? s => l(s) * b(s) : s => l(s) / b(s);
      }
      return a;
    }
    function unary() {
      if (peek().t === '-') { next(); const a = unary(); return s => -a(s); }
      if (peek().t === '+') { next(); return unary(); }
      return power();
    }
    function power() {
      const a = primary();
      if (peek().t === '^') {
        next();
        const b = unary();
        return s => Math.pow(a(s), b(s));
      }
      return a;
    }
    function primary() {
      const x = next();
      if (x.t === 'n') { const v = x.v; return () => v; }
      if (x.t === '(') { const a = expr(); expect(')'); return a; }
      if (x.t === 'id') {
        if (peek().t === '(') {
          const fn = FN[x.v];
          if (!fn) throw new Error('Unknown function "' + x.v + '"');
          next();
          const args = [];
          if (peek().t !== ')') { args.push(expr()); while (peek().t === ',') { next(); args.push(expr()); } }
          expect(')');
          if (args.length === 1) { const a0 = args[0]; return s => fn(a0(s)); }
          if (args.length === 2) { const a0 = args[0], a1 = args[1]; return s => fn(a0(s), a1(s)); }
          return s => fn(...args.map(a => a(s)));
        }
        const idx = vars.indexOf(x.v);
        if (idx >= 0) { used.add(x.v); return s => s[idx]; }
        if (x.v in CONST) { const v = CONST[x.v]; return () => v; }
        throw new Error('Unknown name "' + x.v + '". Available: ' + vars.join(', '));
      }
      if (x.t === 'end') throw new Error('Expression ends too early');
      throw new Error('Unexpected "' + x.t + '" at position ' + (x.p + 1));
    }
    const root = expr();
    if (peek().t !== 'end') throw new Error('Unexpected "' + (peek().v || peek().t) + '" at position ' + (peek().p + 1));
    root.used = [...used];
    return root;
  };

  /* names that appear in an expression (identifiers not followed by "(") */
  X.names = function (src) {
    const tk = tokenize(String(src));
    const out = new Set();
    tk.forEach((t, i) => { if (t.t === 'id' && !(tk[i + 1] && tk[i + 1].t === '(') && !(t.v in CONST)) out.add(t.v); });
    return [...out];
  };
  X.FN = FN;

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
