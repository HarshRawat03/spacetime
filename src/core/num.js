/* Numerical routines: ODEs, roots, quadrature, fitting, special functions, Laplace inversion. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N = {};
  const EPS = 2.220446049250313e-16;

  N.linspace = (a, b, n) => {
    if (n < 2) return [a];
    const h = (b - a) / (n - 1), r = new Array(n);
    for (let i = 0; i < n; i++) r[i] = a + h * i;
    r[n - 1] = b;
    return r;
  };
  N.logspace = (a, b, n) => N.linspace(Math.log10(a), Math.log10(b), n).map(v => Math.pow(10, v));
  N.sum = a => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s; };
  N.mean = a => N.sum(a) / a.length;
  N.clamp = (x, a, b) => Math.min(b, Math.max(a, x));
  N.fin = x => typeof x === 'number' && isFinite(x);
  N.range = n => Array.from({ length: n }, (_, i) => i);
  N.maxAbs = a => a.reduce((m, v) => Math.max(m, Math.abs(v)), 0);

  /* ---------- linear algebra (small dense) ---------- */
  N.lu = function (A) {
    const n = A.length, M = A.map(r => r.slice()), p = N.range(n);
    let sing = false;
    for (let k = 0; k < n; k++) {
      let piv = k, mx = Math.abs(M[k][k]);
      for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > mx) { mx = Math.abs(M[i][k]); piv = i; }
      if (mx < 1e-300) { sing = true; continue; }
      if (piv !== k) { [M[k], M[piv]] = [M[piv], M[k]]; [p[k], p[piv]] = [p[piv], p[k]]; }
      for (let i = k + 1; i < n; i++) {
        const f = M[i][k] / M[k][k];
        M[i][k] = f;
        for (let j = k + 1; j < n; j++) M[i][j] -= f * M[k][j];
      }
    }
    return { M, p, sing };
  };
  N.luSolve = function (F, b) {
    const n = b.length, { M, p } = F, x = new Array(n);
    for (let i = 0; i < n; i++) {
      let s = b[p[i]];
      for (let j = 0; j < i; j++) s -= M[i][j] * x[j];
      x[i] = s;
    }
    for (let i = n - 1; i >= 0; i--) {
      let s = x[i];
      for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
      x[i] = s / (M[i][i] || 1e-300);
    }
    return x;
  };
  N.solve = (A, b) => N.luSolve(N.lu(A), b);
  N.inv = function (A) {
    const n = A.length, F = N.lu(A), R = [];
    for (let j = 0; j < n; j++) {
      const e = new Array(n).fill(0); e[j] = 1;
      const c = N.luSolve(F, e);
      for (let i = 0; i < n; i++) (R[i] || (R[i] = []))[j] = c[i];
    }
    return R;
  };
  N.matmul = (A, B) => A.map(r => B[0].map((_, j) => r.reduce((s, v, k) => s + v * B[k][j], 0)));
  N.transpose = A => A[0].map((_, j) => A.map(r => r[j]));
  N.matvec = (A, x) => A.map(r => r.reduce((s, v, k) => s + v * x[k], 0));

  /* Householder QR least squares: minimize |A x - b|. Returns x and (R^T R)^-1 */
  N.lstsq = function (A, b) {
    const m = A.length, n = A[0].length, Q = A.map(r => r.slice()), y = b.slice();
    const d = new Array(n);
    for (let k = 0; k < n; k++) {
      let nrm = 0;
      for (let i = k; i < m; i++) nrm = Math.hypot(nrm, Q[i][k]);
      if (nrm === 0) { d[k] = 0; continue; }
      if (Q[k][k] > 0) nrm = -nrm;
      for (let i = k; i < m; i++) Q[i][k] /= -nrm;
      Q[k][k] += 1;
      for (let j = k + 1; j < n; j++) {
        let s = 0;
        for (let i = k; i < m; i++) s += Q[i][k] * Q[i][j];
        s = -s / Q[k][k];
        for (let i = k; i < m; i++) Q[i][j] += s * Q[i][k];
      }
      let s = 0;
      for (let i = k; i < m; i++) s += Q[i][k] * y[i];
      s = -s / Q[k][k];
      for (let i = k; i < m; i++) y[i] += s * Q[i][k];
      d[k] = nrm;
    }
    const x = new Array(n);
    for (let k = n - 1; k >= 0; k--) {
      let s = y[k];
      for (let j = k + 1; j < n; j++) s -= Q[k][j] * x[j];
      x[k] = d[k] ? s / d[k] : 0;
    }
    const R = N.range(n).map(i => N.range(n).map(j => (j < i ? 0 : j === i ? d[i] : Q[i][j])));
    const Ri = N.inv(R);
    const covU = N.matmul(Ri, N.transpose(Ri));
    return { x, covU };
  };

  /* ---------- statistics ---------- */
  N.lgamma = function (x) {
    if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - N.lgamma(1 - x);
    const g = 7, c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
      -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    x -= 1;
    let a = c[0];
    const t = x + g + 0.5;
    for (let i = 1; i < 9; i++) a += c[i] / (x + i);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  };
  N.gamma = function (x) {
    if (x < 0.5) return Math.PI / (Math.sin(Math.PI * x) * N.gamma(1 - x));
    return Math.exp(N.lgamma(x));
  };
  /* regularized lower incomplete gamma P(a,x) */
  N.gammaP = function (a, x) {
    if (x <= 0) return 0;
    if (x < a + 1) {
      let ap = a, s = 1 / a, d = s;
      for (let i = 0; i < 500; i++) { ap += 1; d *= x / ap; s += d; if (Math.abs(d) < Math.abs(s) * 1e-15) break; }
      return s * Math.exp(-x + a * Math.log(x) - N.lgamma(a));
    }
    return 1 - N.gammaQcf(a, x);
  };
  N.gammaQcf = function (a, x) {
    let b = x + 1 - a, c = 1 / 1e-300, d = 1 / b, h = d;
    for (let i = 1; i < 500; i++) {
      const an = -i * (i - a);
      b += 2;
      d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 1e-15) break;
    }
    return Math.exp(-x + a * Math.log(x) - N.lgamma(a)) * h;
  };
  N.betacf = function (a, b, x) {
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < 1e-300) d = 1e-300;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 300; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = 1 + aa / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = 1 + aa / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 1e-15) break;
    }
    return h;
  };
  N.betaI = function (a, b, x) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const bt = Math.exp(N.lgamma(a + b) - N.lgamma(a) - N.lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
    if (x < (a + 1) / (a + b + 2)) return bt * N.betacf(a, b, x) / a;
    return 1 - bt * N.betacf(b, a, 1 - x) / b;
  };
  N.tcdf = (t, nu) => {
    const x = nu / (nu + t * t), p = 0.5 * N.betaI(nu / 2, 0.5, x);
    return t > 0 ? 1 - p : p;
  };
  N.tinv = function (p, nu) {
    if (nu <= 0) return NaN;
    let lo = -1e3, hi = 1e3;
    for (let i = 0; i < 200; i++) {
      const m = 0.5 * (lo + hi);
      if (N.tcdf(m, nu) < p) lo = m; else hi = m;
      if (hi - lo < 1e-12) break;
    }
    return 0.5 * (lo + hi);
  };

  /* ---------- special functions ---------- */
  N.erf = function (x) {
    const ax = Math.abs(x);
    if (ax >= 3) return Math.sign(x) * (1 - N.erfc(ax));
    let term = ax, s = ax;
    const x2 = ax * ax;
    for (let n = 1; n < 200; n++) { term *= 2 * x2 / (2 * n + 1); s += term; if (term < s * 1e-17) break; }
    return Math.sign(x) * 2 / Math.sqrt(Math.PI) * Math.exp(-x2) * s;
  };
  N.erfc = function (x) {
    if (x < 3) return 1 - N.erf(x);
    let f = 0;
    for (let n = 60; n >= 1; n--) f = (n / 2) / (x + f);
    return Math.exp(-x * x) / Math.sqrt(Math.PI) / (x + f);
  };
  const EULER = 0.5772156649015329;
  /* E1(x) = integral from x to inf of exp(-u)/u du, x > 0 */
  N.E1 = function (x) {
    if (x <= 0) return Infinity;
    if (x <= 1) {
      let s = 0, term = 1;
      for (let k = 1; k < 100; k++) {
        term *= -x / k;
        const d = -term / k;
        s += d;
        if (Math.abs(d) < 1e-17 * Math.abs(s)) break;
      }
      return -EULER - Math.log(x) + s;
    }
    let b = x + 1, c = 1 / 1e-300, d = 1 / b, h = d;
    for (let i = 1; i < 500; i++) {
      const an = -i * i;
      b += 2;
      d = 1 / (an * d + b);
      c = b + an / c;
      const del = c * d;
      h *= del;
      if (Math.abs(del - 1) < 1e-15) break;
    }
    return h * Math.exp(-x);
  };
  /* Ei(x) for x > 0 (principal value); Ei(-x) = -E1(x) */
  N.Ei = function (x) {
    if (x < 0) return -N.E1(-x);
    if (x === 0) return -Infinity;
    if (x < 40) {
      let s = 0, term = 1;
      for (let k = 1; k < 300; k++) {
        term *= x / k;
        const d = term / k;
        s += d;
        if (d < 1e-17 * s) break;
      }
      return EULER + Math.log(x) + s;
    }
    let s = 1, term = 1;
    for (let k = 1; k < 30; k++) { term *= k / x; s += term; if (term < 1e-17) break; }
    return Math.exp(x) / x * s;
  };
  /* Lambert W, principal branch, x >= -1/e */
  N.lambertW = function (x) {
    if (x === 0) return 0;
    if (x < -1 / Math.E) return NaN;
    let w = x < 1 ? (x < -0.3 ? -1 + Math.sqrt(2 * (1 + Math.E * x)) : x * (1 - x)) : Math.log(1 + x) * (1 - Math.log(1 + Math.log(1 + x)) / (2 + Math.log(1 + x)));
    for (let i = 0; i < 60; i++) {
      const ew = Math.exp(w), f = w * ew - x, wp1 = w + 1;
      if (wp1 === 0) break;
      const dw = f / (ew * wp1 - (w + 2) * f / (2 * wp1));
      w -= dw;
      if (Math.abs(dw) < 1e-15 * (1 + Math.abs(w))) break;
    }
    return w;
  };

  /* ---------- complex numbers (for Laplace inversion) ---------- */
  const C = N.C = (re, im = 0) => ({ re, im });
  C.add = (a, b) => ({ re: a.re + b.re, im: a.im + b.im });
  C.sub = (a, b) => ({ re: a.re - b.re, im: a.im - b.im });
  C.mul = (a, b) => ({ re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re });
  C.div = (a, b) => {
    const d = b.re * b.re + b.im * b.im;
    return { re: (a.re * b.re + a.im * b.im) / d, im: (a.im * b.re - a.re * b.im) / d };
  };
  C.scale = (a, s) => ({ re: a.re * s, im: a.im * s });
  C.exp = a => { const e = Math.exp(a.re); return { re: e * Math.cos(a.im), im: e * Math.sin(a.im) }; };
  C.abs = a => Math.hypot(a.re, a.im);
  C.sqrt = a => {
    const r = C.abs(a);
    let re = Math.sqrt(0.5 * (r + a.re)), im = Math.sqrt(0.5 * Math.max(0, r - a.re));
    if (a.im < 0) im = -im;
    return { re, im };
  };
  C.log = a => ({ re: Math.log(C.abs(a)), im: Math.atan2(a.im, a.re) });
  C.pow = (a, p) => (a.re === 0 && a.im === 0 ? C(0) : C.exp(C.scale(C.log(a), p)));
  C.inv = a => C.div(C(1), a);

  /* Fixed Talbot inversion (Abate & Valko). F(s) takes and returns complex. */
  N.ilt = function (F, t, M = 28) {
    if (t <= 0) return 0;
    const r = 2 * M / (5 * t);
    let s = 0.5 * F(C(r)).re * Math.exp(r * t);
    for (let k = 1; k < M; k++) {
      const th = k * Math.PI / M, ct = Math.cos(th) / Math.sin(th);
      const sk = C(r * th * ct, r * th);
      const sig = th + (th * ct - 1) * ct;
      const term = C.mul(C.mul(C.exp(C.scale(sk, t)), F(sk)), C(1, sig));
      s += term.re;
    }
    return r / M * s;
  };

  /* ---------- roots ---------- */
  N.brent = function (f, a, b, o = {}) {
    const tol = o.tol || 1e-12, maxIt = o.maxIt || 200;
    let fa = o.fa !== undefined ? o.fa : f(a), fb = o.fb !== undefined ? o.fb : f(b);
    if (!isFinite(fa) || !isFinite(fb)) return NaN;
    if (fa === 0) return a;
    if (fb === 0) return b;
    if (fa * fb > 0) return NaN;
    let c = a, fc = fa, d = b - a, e = d;
    for (let i = 0; i < maxIt; i++) {
      if (fb * fc > 0) { c = a; fc = fa; d = b - a; e = d; }
      if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
      const tol1 = 2 * EPS * Math.abs(b) + 0.5 * tol, xm = 0.5 * (c - b);
      if (Math.abs(xm) <= tol1 || fb === 0) return b;
      if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
        const s = fb / fa;
        let p, q;
        if (a === c) { p = 2 * xm * s; q = 1 - s; }
        else {
          const qq = fa / fc, r = fb / fc;
          p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1));
          q = (qq - 1) * (r - 1) * (s - 1);
        }
        if (p > 0) q = -q;
        p = Math.abs(p);
        if (2 * p < Math.min(3 * xm * q - Math.abs(tol1 * q), Math.abs(e * q))) { e = d; d = p / q; }
        else { d = xm; e = d; }
      } else { d = xm; e = d; }
      a = b; fa = fb;
      b += Math.abs(d) > tol1 ? d : (xm > 0 ? tol1 : -tol1);
      fb = f(b);
      if (!isFinite(fb)) return NaN;
    }
    return b;
  };
  /* all sign-change roots of f on [a,b] by scanning n panels */
  N.roots = function (f, a, b, n = 400, o = {}) {
    const xs = o.log ? N.logspace(a, b, n + 1) : N.linspace(a, b, n + 1);
    const fs = xs.map(f), out = [];
    for (let i = 0; i < n; i++) {
      const f0 = fs[i], f1 = fs[i + 1];
      if (!isFinite(f0) || !isFinite(f1)) continue;
      if (f0 === 0) { out.push(xs[i]); continue; }
      if (f0 * f1 < 0) {
        const r = N.brent(f, xs[i], xs[i + 1], { fa: f0, fb: f1, tol: o.tol || 1e-13 });
        if (isFinite(r)) out.push(r);
      }
    }
    if (fs[n] === 0) out.push(xs[n]);
    return out;
  };
  N.root = function (f, a, b, o = {}) {
    let fa = f(a), fb = f(b);
    if (isFinite(fa) && isFinite(fb) && fa * fb <= 0) return N.brent(f, a, b, { fa, fb, tol: o.tol });
    const r = N.roots(f, a, b, o.n || 200, o);
    return r.length ? r[0] : NaN;
  };
  N.jac = function (F, x, Fx) {
    const n = x.length, f0 = Fx || F(x), J = f0.map(() => new Array(n));
    for (let j = 0; j < n; j++) {
      const h = 1e-7 * Math.max(Math.abs(x[j]), 1e-6), xj = x.slice();
      xj[j] += h;
      const f1 = F(xj);
      for (let i = 0; i < f0.length; i++) J[i][j] = (f1[i] - f0[i]) / h;
    }
    return J;
  };
  N.newton = function (F, x0, o = {}) {
    const tol = o.tol || 1e-11, maxIt = o.maxIt || 100, lo = o.lo, hi = o.hi;
    let x = x0.slice(), Fx = F(x), nrm = Math.hypot(...Fx), it = 0;
    const clampX = v => v.map((xi, i) => {
      let r = xi;
      if (lo && lo[i] !== undefined) r = Math.max(lo[i], r);
      if (hi && hi[i] !== undefined) r = Math.min(hi[i], r);
      return r;
    });
    for (; it < maxIt && nrm > tol; it++) {
      const J = o.J ? o.J(x) : N.jac(F, x, Fx);
      const dx = N.solve(J, Fx.map(v => -v));
      if (dx.some(v => !isFinite(v))) break;
      let lam = 1, xn, Fn, nn;
      for (let k = 0; k < 30; k++) {
        xn = clampX(x.map((xi, i) => xi + lam * dx[i]));
        Fn = F(xn);
        nn = Math.hypot(...Fn);
        if (isFinite(nn) && nn < (1 - 1e-4 * lam) * nrm) break;
        lam *= 0.5;
      }
      if (!isFinite(nn)) break;
      const step = Math.max(...xn.map((v, i) => Math.abs(v - x[i]) / (Math.abs(v) + 1e-12)));
      x = xn; Fx = Fn; nrm = nn;
      if (step < 1e-14) break;
    }
    return { x, ok: nrm <= Math.max(tol, 1e-8 * (1 + N.maxAbs(x))), it, norm: nrm };
  };

  /* ---------- minimisation ---------- */
  N.fmin = function (f, a, b, o = {}) {
    const tol = o.tol || 1e-10, cg = 0.3819660112501051;
    let x = a + cg * (b - a), w = x, v = x, fx = f(x), fw = fx, fv = fx, d = 0, e = 0;
    for (let i = 0; i < (o.maxIt || 200); i++) {
      const xm = 0.5 * (a + b), tol1 = tol * Math.abs(x) + 1e-12, tol2 = 2 * tol1;
      if (Math.abs(x - xm) <= tol2 - 0.5 * (b - a)) break;
      let useGolden = true;
      if (Math.abs(e) > tol1) {
        let r = (x - w) * (fx - fv), q = (x - v) * (fx - fw), p = (x - v) * q - (x - w) * r;
        q = 2 * (q - r);
        if (q > 0) p = -p;
        q = Math.abs(q);
        const et = e;
        e = d;
        if (!(Math.abs(p) >= Math.abs(0.5 * q * et) || p <= q * (a - x) || p >= q * (b - x))) {
          d = p / q;
          const u = x + d;
          if (u - a < tol2 || b - u < tol2) d = xm >= x ? tol1 : -tol1;
          useGolden = false;
        }
      }
      if (useGolden) { e = (x >= xm ? a : b) - x; d = cg * e; }
      const u = x + (Math.abs(d) >= tol1 ? d : (d > 0 ? tol1 : -tol1)), fu = f(u);
      if (fu <= fx) {
        if (u >= x) a = x; else b = x;
        v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
      } else {
        if (u < x) a = u; else b = u;
        if (fu <= fw || w === x) { v = w; fv = fw; w = u; fw = fu; }
        else if (fu <= fv || v === x || v === w) { v = u; fv = fu; }
      }
    }
    return { x, f: fx };
  };
  /* global-ish 1D minimum: coarse scan then Brent around best */
  N.fminScan = function (f, a, b, n = 60, o = {}) {
    const xs = o.log ? N.logspace(a, b, n) : N.linspace(a, b, n);
    let bi = 0, bf = Infinity;
    xs.forEach((x, i) => { const v = f(x); if (v < bf) { bf = v; bi = i; } });
    const lo = xs[Math.max(0, bi - 1)], hi = xs[Math.min(n - 1, bi + 1)];
    const r = N.fmin(f, lo, hi, o);
    return r.f <= bf ? r : { x: xs[bi], f: bf };
  };
  N.nelderMead = function (f, x0, o = {}) {
    const n = x0.length, maxIt = o.maxIt || 400 * n, tol = o.tol || 1e-10;
    const step = o.step || x0.map(v => (Math.abs(v) > 1e-8 ? 0.1 * v : 0.05));
    let P = [x0.slice()];
    for (let i = 0; i < n; i++) { const p = x0.slice(); p[i] += step[i]; P.push(p); }
    let F = P.map(f);
    for (let it = 0; it < maxIt; it++) {
      const idx = N.range(n + 1).sort((a, b) => F[a] - F[b]);
      P = idx.map(i => P[i]); F = idx.map(i => F[i]);
      if (Math.abs(F[n] - F[0]) <= tol * (Math.abs(F[0]) + 1e-15)) {
        let spread = 0;
        for (let i = 1; i <= n; i++) for (let j = 0; j < n; j++) spread = Math.max(spread, Math.abs(P[i][j] - P[0][j]));
        if (spread < 1e-9 * (1 + N.maxAbs(P[0]))) break;
      }
      const c = new Array(n).fill(0);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) c[j] += P[i][j] / n;
      const xr = c.map((v, j) => v + (v - P[n][j])), fr = f(xr);
      if (fr < F[0]) {
        const xe = c.map((v, j) => v + 2 * (v - P[n][j])), fe = f(xe);
        if (fe < fr) { P[n] = xe; F[n] = fe; } else { P[n] = xr; F[n] = fr; }
      } else if (fr < F[n - 1]) { P[n] = xr; F[n] = fr; }
      else {
        const out = fr < F[n];
        const xc = out ? c.map((v, j) => v + 0.5 * (xr[j] - v)) : c.map((v, j) => v + 0.5 * (P[n][j] - v));
        const fc = f(xc);
        if (fc < (out ? fr : F[n])) { P[n] = xc; F[n] = fc; }
        else {
          for (let i = 1; i <= n; i++) { P[i] = P[i].map((v, j) => P[0][j] + 0.5 * (v - P[0][j])); F[i] = f(P[i]); }
        }
      }
    }
    let bi = 0;
    for (let i = 1; i <= n; i++) if (F[i] < F[bi]) bi = i;
    return { x: P[bi], f: F[bi] };
  };

  /* ---------- quadrature ---------- */
  const XGK = [0.991455371120812639206854697526329, 0.949107912342758524526189684047851, 0.864864423359769072789712788640926,
    0.741531185599394439863864773280788, 0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0];
  const WGK = [0.022935322010529224963732008058970, 0.063092092629978553290700663189204, 0.104790010322250183839876322541518,
    0.140653259715525918745189590510238, 0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714];
  const WG = [0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975,
    0.417959183673469387755102040816327];
  function gk15(f, a, b) {
    const c = 0.5 * (a + b), h = 0.5 * (b - a), fc = f(c);
    let rk = fc * WGK[7], rg = fc * WG[3];
    for (let j = 0; j < 7; j++) {
      const x = h * XGK[j], f1 = f(c - x), f2 = f(c + x);
      rk += WGK[j] * (f1 + f2);
      if (j % 2 === 1) rg += WG[(j - 1) / 2] * (f1 + f2);
    }
    return { I: rk * h, err: Math.abs((rk - rg) * h), a, b };
  }
  N.quad = function (f, a, b, o = {}) {
    if (a === b) return 0;
    if (b === Infinity) {
      const g = u => { const t = a + u / (1 - u); return f(t) / ((1 - u) * (1 - u)); };
      return N.quad(g, 0, 1, o);
    }
    const rtol = o.rtol || 1e-10, atol = o.atol || 1e-13, maxSeg = o.maxSeg || 400;
    let segs = [gk15(f, a, b)];
    for (let k = 0; k < maxSeg; k++) {
      let I = 0, E = 0, wi = 0;
      segs.forEach((s, i) => { I += s.I; E += s.err; if (s.err > segs[wi].err) wi = i; });
      if (E <= Math.max(atol, rtol * Math.abs(I)) || !isFinite(I)) return I;
      const s = segs[wi], m = 0.5 * (s.a + s.b);
      segs.splice(wi, 1, gk15(f, s.a, m), gk15(f, m, s.b));
    }
    return segs.reduce((t, s) => t + s.I, 0);
  };
  N.trapz = function (x, y) {
    let s = 0;
    for (let i = 1; i < x.length; i++) s += 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    return s;
  };
  N.cumtrapz = function (x, y) {
    const r = [0];
    for (let i = 1; i < x.length; i++) r.push(r[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]));
    return r;
  };

  /* ---------- interpolation ---------- */
  N.interp = function (xs, ys, x) {
    const n = xs.length;
    if (n === 0) return NaN;
    if (n === 1) return ys[0];
    const asc = xs[n - 1] >= xs[0];
    if (asc ? x <= xs[0] : x >= xs[0]) return ys[0] + (ys[1] - ys[0]) * (x - xs[0]) / (xs[1] - xs[0]);
    if (asc ? x >= xs[n - 1] : x <= xs[n - 1]) return ys[n - 2] + (ys[n - 1] - ys[n - 2]) * (x - xs[n - 2]) / (xs[n - 1] - xs[n - 2]);
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      const m = (lo + hi) >> 1;
      if (asc ? xs[m] <= x : xs[m] >= x) lo = m; else hi = m;
    }
    const dx = xs[hi] - xs[lo];
    return dx === 0 ? ys[lo] : ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / dx;
  };
  N.interpClamp = function (xs, ys, x) {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
    return N.interp(xs, ys, x);
  };
  /* monotone cubic (Fritsch-Carlson) */
  N.pchip = function (xs, ys) {
    const n = xs.length;
    if (n < 3) return x => N.interp(xs, ys, x);
    const h = [], d = [], m = new Array(n);
    for (let i = 0; i < n - 1; i++) { h[i] = xs[i + 1] - xs[i]; d[i] = (ys[i + 1] - ys[i]) / h[i]; }
    m[0] = d[0]; m[n - 1] = d[n - 2];
    for (let i = 1; i < n - 1; i++) {
      if (d[i - 1] * d[i] <= 0) m[i] = 0;
      else { const w1 = 2 * h[i] + h[i - 1], w2 = h[i] + 2 * h[i - 1]; m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]); }
    }
    return function (x) {
      if (x <= xs[0]) return ys[0] + m[0] * (x - xs[0]);
      if (x >= xs[n - 1]) return ys[n - 1] + m[n - 1] * (x - xs[n - 1]);
      let lo = 0, hi = n - 1;
      while (hi - lo > 1) { const k = (lo + hi) >> 1; if (xs[k] <= x) lo = k; else hi = k; }
      const t = (x - xs[lo]) / h[lo], t2 = t * t, t3 = t2 * t;
      return (2 * t3 - 3 * t2 + 1) * ys[lo] + (t3 - 2 * t2 + t) * h[lo] * m[lo] + (-2 * t3 + 3 * t2) * ys[hi] + (t3 - t2) * h[lo] * m[hi];
    };
  };

  /* ---------- ODE integration ---------- */
  const DPc = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
  const DPa = [[], [1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9], [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656], [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]];
  const DPe = [71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40];

  function hermite(t0, y0, f0, t1, y1, f1, t) {
    const h = t1 - t0, s = (t - t0) / h, s2 = s * s, s3 = s2 * s;
    const h00 = 2 * s3 - 3 * s2 + 1, h10 = s3 - 2 * s2 + s, h01 = -2 * s3 + 3 * s2, h11 = s3 - s2;
    return y0.map((v, i) => h00 * v + h10 * h * f0[i] + h01 * y1[i] + h11 * h * f1[i]);
  }
  function makeSol(T, Y, F, status, ev) {
    const sol = { t: T, y: Y, status, ev };
    sol.at = function (tq) {
      const n = T.length;
      if (n === 1) return Y[0].slice();
      const asc = T[n - 1] >= T[0];
      if (asc ? tq <= T[0] : tq >= T[0]) return Y[0].slice();
      if (asc ? tq >= T[n - 1] : tq <= T[n - 1]) return Y[n - 1].slice();
      let lo = 0, hi = n - 1;
      while (hi - lo > 1) { const m = (lo + hi) >> 1; if (asc ? T[m] <= tq : T[m] >= tq) lo = m; else hi = m; }
      return hermite(T[lo], Y[lo], F[lo], T[hi], Y[hi], F[hi], tq);
    };
    sol.col = i => Y.map(r => r[i]);
    sol.last = () => Y[Y.length - 1];
    sol.sample = function (n) {
      const ts = N.linspace(T[0], T[T.length - 1], n);
      return { t: ts, y: ts.map(sol.at) };
    };
    return sol;
  }
  function errNorm(err, y0, y1, rtol, atol) {
    let s = 0;
    for (let i = 0; i < err.length; i++) {
      const sc = atol + rtol * Math.max(Math.abs(y0[i]), Math.abs(y1[i]));
      s += (err[i] / sc) ** 2;
    }
    return Math.sqrt(s / err.length);
  }
  function locateEvent(o, t0, y0, f0, t1, y1, f1, g0, g1) {
    let a = t0, b = t1, ga = g0;
    for (let i = 0; i < 80; i++) {
      const m = 0.5 * (a + b), ym = hermite(t0, y0, f0, t1, y1, f1, m), gm = o.event(m, ym);
      if (ga * gm <= 0) b = m; else { a = m; ga = gm; }
      if (Math.abs(b - a) < 1e-14 * (1 + Math.abs(m))) break;
    }
    const tE = 0.5 * (a + b);
    return { t: tE, y: hermite(t0, y0, f0, t1, y1, f1, tE) };
  }

  function dpStep(f, t, y, fy, hd) {
    const k = [fy];
    for (let s = 1; s < 7; s++) {
      const ys = y.map((v, i) => { let a = v; for (let j = 0; j < s; j++) a += hd * DPa[s][j] * k[j][i]; return a; });
      k[s] = f(t + DPc[s] * hd, ys);
    }
    return y.map((v, i) => { let a = v; for (let j = 0; j < 6; j++) a += hd * DPa[6][j] * k[j][i]; return a; });
  }
  function refineEvent(f, o, t, y, fy, ev, t1, g1) {
    let ta = ev.t, ga = o.event(ta, dpStep(f, t, y, fy, ta - t)), tb = t1, gb = g1;
    if (ga === 0) return { t: ta, y: dpStep(f, t, y, fy, ta - t) };
    let tc = ta - 1e-6 * (t1 - t), gc = o.event(tc, dpStep(f, t, y, fy, tc - t));
    for (let i = 0; i < 8; i++) {
      if (gc === ga) break;
      const tn = ta - ga * (ta - tc) / (ga - gc);
      if (!isFinite(tn)) break;
      tc = ta; gc = ga; ta = tn; ga = o.event(ta, dpStep(f, t, y, fy, ta - t));
      if (Math.abs(ta - tc) < 1e-15 * (1 + Math.abs(ta))) break;
    }
    if (!(Math.min(t, tb) <= ta && ta <= Math.max(t, tb))) return ev;
    return { t: ta, y: dpStep(f, t, y, fy, ta - t) };
  }

  N.rk45 = function (f, t0, y0, t1, o = {}) {
    const rtol = o.rtol || 1e-8, atol = o.atol || 1e-11, maxSteps = o.maxSteps || 40000;
    const dir = t1 >= t0 ? 1 : -1, span = Math.abs(t1 - t0);
    const hmax = o.hmax || span / (o.minPts || 50);
    let t = t0, y = y0.slice(), fy = f(t, y);
    const T = [t], Y = [y.slice()], F = [fy.slice()];
    let h = o.h0 || Math.min(hmax, span * 1e-3), steps = 0, g0 = o.event ? o.event(t, y) : 0;
    const k = new Array(7);
    while (dir * (t1 - t) > 1e-14 * Math.max(1, Math.abs(t1))) {
      if (++steps > maxSteps) return makeSol(T, Y, F, 'maxsteps', null);
      h = Math.min(h, Math.abs(t1 - t), hmax);
      const hd = dir * h;
      k[0] = fy;
      for (let s = 1; s < 7; s++) {
        const ys = y.map((v, i) => { let a = v; for (let j = 0; j < s; j++) a += hd * DPa[s][j] * k[j][i]; return a; });
        k[s] = f(t + DPc[s] * hd, ys);
      }
      const yn = y.map((v, i) => { let a = v; for (let j = 0; j < 6; j++) a += hd * DPa[6][j] * k[j][i]; return a; });
      const err = y.map((_, i) => { let a = 0; for (let j = 0; j < 7; j++) a += hd * DPe[j] * k[j][i]; return a; });
      const en = errNorm(err, y, yn, rtol, atol);
      if (!isFinite(en)) { h *= 0.25; if (h < 1e-14 * span) return makeSol(T, Y, F, 'fail', null); continue; }
      if (en <= 1) {
        const tn = t + hd, fn = k[6];
        if (o.event) {
          const g1 = o.event(tn, yn);
          if (g0 * g1 < 0 || g1 === 0) {
            const ev = refineEvent(f, o, t, y, fy, locateEvent(o, t, y, fy, tn, yn, fn, g0, g1), tn, g1);
            T.push(ev.t); Y.push(ev.y); F.push(f(ev.t, ev.y));
            return makeSol(T, Y, F, 'event', ev);
          }
          g0 = g1;
        }
        t = tn; y = yn; fy = fn;
        T.push(t); Y.push(y.slice()); F.push(fy.slice());
        if (o.stop && o.stop(t, y)) return makeSol(T, Y, F, 'stopped', null);
      }
      h *= Math.min(5, Math.max(0.2, 0.9 * Math.pow(en || 1e-10, -0.2)));
      if (h < 1e-14 * span) return makeSol(T, Y, F, 'fail', null);
    }
    return makeSol(T, Y, F, 'ok', null);
  };

  /* Rosenbrock ode23s (Shampine & Reichelt) for stiff problems */
  N.ode23s = function (f, t0, y0, t1, o = {}) {
    const rtol = o.rtol || 1e-6, atol = o.atol || 1e-10, maxSteps = o.maxSteps || 40000;
    const dir = t1 >= t0 ? 1 : -1, span = Math.abs(t1 - t0), hmax = o.hmax || span / (o.minPts || 50);
    const d = 1 / (2 + Math.SQRT2), e32 = 6 + Math.SQRT2, n = y0.length;
    let t = t0, y = y0.slice(), F0 = f(t, y), h = o.h0 || span * 1e-4, steps = 0, g0 = o.event ? o.event(t, y) : 0;
    const T = [t], Y = [y.slice()], FF = [F0.slice()];
    while (dir * (t1 - t) > 1e-14 * Math.max(1, Math.abs(t1))) {
      if (++steps > maxSteps) return makeSol(T, Y, FF, 'maxsteps', null);
      h = Math.min(h, Math.abs(t1 - t), hmax);
      const hd = dir * h;
      const J = N.jac(yy => f(t, yy), y, F0);
      const dt = 1e-7 * Math.max(Math.abs(t), 1e-6) * dir;
      const Ft = f(t + dt, y), Tv = Ft.map((v, i) => (v - F0[i]) / dt);
      const W = N.range(n).map(i => N.range(n).map(j => (i === j ? 1 : 0) - hd * d * J[i][j]));
      const L = N.lu(W);
      const k1 = N.luSolve(L, F0.map((v, i) => v + hd * d * Tv[i]));
      const F1 = f(t + 0.5 * hd, y.map((v, i) => v + 0.5 * hd * k1[i]));
      const k2 = N.luSolve(L, F1.map((v, i) => v - k1[i])).map((v, i) => v + k1[i]);
      const yn = y.map((v, i) => v + hd * k2[i]);
      const F2 = f(t + hd, yn);
      const k3 = N.luSolve(L, F2.map((v, i) => v - e32 * (k2[i] - F1[i]) - 2 * (k1[i] - F0[i]) + hd * d * Tv[i]));
      const err = k1.map((v, i) => hd / 6 * (v - 2 * k2[i] + k3[i]));
      const en = errNorm(err, y, yn, rtol, atol);
      if (!isFinite(en)) { h *= 0.25; if (h < 1e-15 * span) return makeSol(T, Y, FF, 'fail', null); continue; }
      if (en <= 1) {
        const tn = t + hd;
        if (o.event) {
          const g1 = o.event(tn, yn);
          if (g0 * g1 < 0 || g1 === 0) {
            const ev = locateEvent(o, t, y, F0, tn, yn, F2, g0, g1);
            T.push(ev.t); Y.push(ev.y); FF.push(f(ev.t, ev.y));
            return makeSol(T, Y, FF, 'event', ev);
          }
          g0 = g1;
        }
        t = tn; y = yn; F0 = F2;
        T.push(t); Y.push(y.slice()); FF.push(F0.slice());
        if (o.stop && o.stop(t, y)) return makeSol(T, Y, FF, 'stopped', null);
      }
      h *= Math.min(5, Math.max(0.2, 0.8 * Math.pow(en || 1e-10, -1 / 3)));
      if (h < 1e-15 * span) return makeSol(T, Y, FF, 'fail', null);
    }
    return makeSol(T, Y, FF, 'ok', null);
  };

  /* default integrator: explicit first, stiff fallback */
  N.ode = function (f, t0, y0, t1, o = {}) {
    if (o.stiff) return N.ode23s(f, t0, y0, t1, o);
    const r = N.rk45(f, t0, y0, t1, Object.assign({}, o, { maxSteps: o.maxSteps || 6000 }));
    if (r.status === 'maxsteps' || r.status === 'fail') {
      const s = N.ode23s(f, t0, y0, t1, Object.assign({}, o, { rtol: Math.max(o.rtol || 0, 2e-7), atol: Math.max(o.atol || 0, 1e-13), maxSteps: 60000 }));
      s.stiff = true;
      return s;
    }
    return r;
  };

  /* ---------- regression ---------- */
  N.linreg = function (x, y, o = {}) {
    const n = x.length;
    if (o.origin) {
      let sxy = 0, sxx = 0;
      for (let i = 0; i < n; i++) { sxy += x[i] * y[i]; sxx += x[i] * x[i]; }
      const b = sxy / sxx, res = y.map((v, i) => v - b * x[i]), ssr = N.sum(res.map(v => v * v));
      const ym = N.mean(y), sst = N.sum(y.map(v => (v - ym) ** 2));
      const s2 = n > 1 ? ssr / (n - 1) : NaN;
      return { a: 0, b, sa: 0, sb: Math.sqrt(s2 / sxx), r2: 1 - ssr / sst, n, res, ssr, dof: n - 1 };
    }
    const xm = N.mean(x), ym = N.mean(y);
    let sxx = 0, sxy = 0, syy = 0;
    for (let i = 0; i < n; i++) { const dx = x[i] - xm, dy = y[i] - ym; sxx += dx * dx; sxy += dx * dy; syy += dy * dy; }
    const b = sxy / sxx, a = ym - b * xm;
    const res = y.map((v, i) => v - a - b * x[i]), ssr = N.sum(res.map(v => v * v));
    const dof = n - 2, s2 = dof > 0 ? ssr / dof : NaN;
    const sb = Math.sqrt(s2 / sxx), sa = Math.sqrt(s2 * (1 / n + xm * xm / sxx));
    return { a, b, sa, sb, r2: syy > 0 ? 1 - ssr / syy : 1, n, res, ssr, dof, s: Math.sqrt(s2), xm, sxx };
  };
  N.polyfit = function (x, y, deg) {
    const A = x.map(v => N.range(deg + 1).map(p => Math.pow(v, p)));
    const r = N.lstsq(A, y);
    const res = y.map((v, i) => v - N.polyval(r.x, x[i]));
    return { c: r.x, res, ssr: N.sum(res.map(v => v * v)) };
  };
  N.polyval = (c, x) => { let s = 0; for (let i = c.length - 1; i >= 0; i--) s = s * x + c[i]; return s; };
  N.polyder = c => c.slice(1).map((v, i) => v * (i + 1));

  /* Levenberg-Marquardt on residual function r(p) -> array.
     o.pos[i] true fits parameter i in log space (keeps it positive). */
  N.lm = function (resid, p0, o = {}) {
    const np = p0.length, pos = o.pos || [], fixed = o.fixed || [];
    const toU = p => p.map((v, i) => (pos[i] ? Math.log(Math.max(v, 1e-300)) : v));
    const toP = u => u.map((v, i) => (pos[i] ? Math.exp(v) : v));
    const free = N.range(np).filter(i => !fixed[i]);
    let u = toU(p0);
    const R = uu => resid(toP(uu));
    let r = R(u), ssr = N.sum(r.map(v => v * v)), lam = 1e-3, it = 0;
    if (!isFinite(ssr)) return { p: p0.slice(), ok: false, ssr, msg: 'model not finite at start' };
    const maxIt = o.maxIt || 300;
    let J = null;
    for (; it < maxIt; it++) {
      J = r.map(() => new Array(free.length));
      free.forEach((pi, jj) => {
        const h = 1e-6 * Math.max(Math.abs(u[pi]), 1e-3), uu = u.slice();
        uu[pi] += h;
        const r1 = R(uu);
        for (let i = 0; i < r.length; i++) J[i][jj] = (r1[i] - r[i]) / h;
      });
      const m = free.length;
      const A = N.range(m).map(a => N.range(m).map(b => { let s = 0; for (let i = 0; i < r.length; i++) s += J[i][a] * J[i][b]; return s; }));
      const g = N.range(m).map(a => { let s = 0; for (let i = 0; i < r.length; i++) s += J[i][a] * r[i]; return s; });
      let improved = false;
      for (let tries = 0; tries < 12; tries++) {
        const Ad = A.map((row, a) => row.map((v, b) => (a === b ? v * (1 + lam) + 1e-12 : v)));
        const dl = N.solve(Ad, g.map(v => -v));
        const un = u.slice();
        free.forEach((pi, jj) => { un[pi] += dl[jj]; });
        const rn = R(un), sn = N.sum(rn.map(v => v * v));
        if (isFinite(sn) && sn < ssr) {
          const rel = (ssr - sn) / (ssr + 1e-300);
          u = un; r = rn; ssr = sn; lam = Math.max(lam / 5, 1e-12); improved = true;
          if (rel < 1e-13) it = maxIt;
          break;
        }
        lam *= 8;
      }
      if (!improved) break;
    }
    const p = toP(u), n = r.length, dof = Math.max(1, n - free.length);
    let se = new Array(np).fill(NaN), cov = null;
    try {
      const m = free.length;
      const A = N.range(m).map(a => N.range(m).map(b => { let s = 0; for (let i = 0; i < n; i++) s += J[i][a] * J[i][b]; return s; }));
      const Ai = N.inv(A), s2 = ssr / dof;
      cov = Ai.map(row => row.map(v => v * s2));
      free.forEach((pi, jj) => {
        const su = Math.sqrt(Math.max(cov[jj][jj], 0));
        se[pi] = pos[pi] ? p[pi] * su : su;
      });
    } catch (e) { /* singular */ }
    const aic = n * Math.log(Math.max(ssr, 1e-300) / n) + 2 * free.length;
    return { p, se, ssr, n, dof, aic, ok: true, it, resid: r };
  };

  /* simple statistics on residuals vs observations */
  N.r2 = function (yObs, yFit) {
    const m = N.mean(yObs);
    let ssr = 0, sst = 0;
    for (let i = 0; i < yObs.length; i++) { ssr += (yObs[i] - yFit[i]) ** 2; sst += (yObs[i] - m) ** 2; }
    return sst > 0 ? 1 - ssr / sst : 1;
  };

  /* ---------- derivatives from data ---------- */
  N.diffData = function (x, y, method = 'central', o = {}) {
    const n = x.length;
    if (method === 'poly') {
      const deg = Math.min(o.deg || 3, n - 1), f = N.polyfit(x, y, deg), dc = N.polyder(f.c);
      return x.map(v => N.polyval(dc, v));
    }
    if (method === 'local') {
      const w = Math.max(1, o.half || 2);
      return x.map((_, i) => {
        const a = Math.max(0, i - w), b = Math.min(n - 1, i + w), xs = x.slice(a, b + 1), ys = y.slice(a, b + 1);
        const deg = Math.min(2, xs.length - 1), x0 = x[i];
        const f = N.polyfit(xs.map(v => v - x0), ys, deg);
        return f.c[1] !== undefined ? f.c[1] : 0;
      });
    }
    return x.map((_, i) => {
      if (i === 0) {
        if (n > 2) {
          const h1 = x[1] - x[0], h2 = x[2] - x[1];
          return (-(2 * h1 + h2) / (h1 * (h1 + h2))) * y[0] + ((h1 + h2) / (h1 * h2)) * y[1] - (h1 / (h2 * (h1 + h2))) * y[2];
        }
        return (y[1] - y[0]) / (x[1] - x[0]);
      }
      if (i === n - 1) {
        if (n > 2) {
          const h1 = x[n - 2] - x[n - 3], h2 = x[n - 1] - x[n - 2];
          return (h2 / (h1 * (h1 + h2))) * y[n - 3] - ((h1 + h2) / (h1 * h2)) * y[n - 2] + ((2 * h2 + h1) / (h2 * (h1 + h2))) * y[n - 1];
        }
        return (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
      }
      const h1 = x[i] - x[i - 1], h2 = x[i + 1] - x[i];
      return (-h2 / (h1 * (h1 + h2))) * y[i - 1] + ((h2 - h1) / (h1 * h2)) * y[i] + (h1 / (h2 * (h1 + h2))) * y[i + 1];
    });
  };

  /* ---------- number formatting ---------- */
  const SUP = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
  N.sup = s => String(s).split('').map(c => SUP[c] || c).join('');
  N.fmt = function (v, sig = 4, o = {}) {
    if (v === null || v === undefined) return '–';
    if (typeof v === 'string') return v;
    if (Number.isNaN(v)) return '–';
    if (v === Infinity) return '∞';
    if (v === -Infinity) return '−∞';
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (Number.isInteger(v) && a < 1e7 && !o.fixed) return (v < 0 ? '−' : '') + String(a);
    if (a >= 1e6 || a < 1e-3) {
      const e = Math.floor(Math.log10(a));
      let m = v / Math.pow(10, e);
      let ms = m.toPrecision(sig);
      if (Math.abs(parseFloat(ms)) >= 10) { m /= 10; ms = m.toPrecision(sig); return fmtSci(ms, e + 1, o); }
      return fmtSci(ms, e, o);
    }
    let s = v.toPrecision(sig);
    if (s.indexOf('e') >= 0) s = String(parseFloat(s));
    if (!o.keepZeros && s.indexOf('.') >= 0) s = s.replace(/(\.\d*?)0+$/, '$1').replace(/\.$/, '');
    return s.replace('-', '−');
  };
  function fmtSci(ms, e, o) {
    ms = ms.replace(/(\.\d*?)0+$/, '$1').replace(/\.$/, '').replace('-', '−');
    if (o.plain) return ms + 'e' + e;
    if (o.html) return ms + ' × 10<sup>' + String(e).replace('-', '−') + '</sup>';
    return ms + '×10' + N.sup(e);
  }
  N.fmtCSV = v => (typeof v === 'number' ? (isFinite(v) ? String(+v.toPrecision(10)) : '') : '"' + String(v).replace(/"/g, '""') + '"');

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
