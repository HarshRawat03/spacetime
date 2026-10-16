/* Residence-time distributions: model E curves, a time-grid representation with Dirac atoms,
   convolution (FFT), recycle, moments, and conversion by segregation or maximum mixedness. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N, C = N.C;
  const RT = ST.RT = {};

  /* ---------- analytic E curves, E(t) with mean residence time tau ---------- */
  RT.Emixed = (t, tau) => (t < 0 ? 0 : Math.exp(-t / tau) / tau);
  RT.Etanks = (t, tau, Nk) => {
    if (t <= 0) return Nk === 1 ? 1 / tau : 0;
    return Math.exp(Nk * Math.log(Nk / tau) + (Nk - 1) * Math.log(t) - Nk * t / tau - N.lgamma(Nk));
  };
  RT.Ftanks = (t, tau, Nk) => (t <= 0 ? 0 : N.gammaP(Nk, Nk * t / tau));
  RT.Egauss = (t, tau, d) => Math.exp(-Math.pow(1 - t / tau, 2) / (4 * d)) / (tau * Math.sqrt(4 * Math.PI * d));
  RT.Eopen = (t, tau, d) => { const th = t / tau; if (th <= 0) return 0; return Math.exp(-Math.pow(1 - th, 2) / (4 * d * th)) / (tau * Math.sqrt(4 * Math.PI * d * th)); };
  RT.Elam = (t, tau) => (t >= tau / 2 ? tau * tau / (2 * t * t * t) : 0);
  RT.Flam = (t, tau) => (t >= tau / 2 ? 1 - tau * tau / (4 * t * t) : 0);
  /* variances in θ units */
  RT.varClosed = d => 2 * d - 2 * d * d * (1 - Math.exp(-1 / d));
  RT.varOpen = d => 2 * d + 8 * d * d;
  /* closed–closed dispersion: transfer function written to avoid overflow, inverted by fixed Talbot */
  const Gclosed = Pe => s => {
    const a = C.sqrt(C.add(C(1), C.scale(s, 4 / Pe)));
    const one = C(1);
    const ap = C.add(one, a), am = C.sub(one, a);
    const num = C.scale(C.mul(a, C.exp(C.scale(am, Pe / 2))), 4);
    const den = C.sub(C.mul(ap, ap), C.mul(C.mul(am, am), C.exp(C.scale(a, -Pe))));
    return C.div(num, den);
  };
  RT.EthClosed = function (th, d, M = 32) {
    if (th <= 0) return 0;
    if (d < 0.008) { const v = RT.varClosed(d); return Math.exp(-Math.pow(1 - th, 2) / (2 * v)) / Math.sqrt(2 * Math.PI * v); }
    const e = N.ilt(Gclosed(1 / d), th, M);
    return e > 0 ? e : 0;
  };
  RT.Eclosed = (t, tau, d) => RT.EthClosed(t / tau, d) / tau;
  /* first-order conversion with dispersion, closed vessel (Wehner–Wilhelm), stable form */
  RT.wehner = function (kt, d) {
    if (d <= 0) return Math.exp(-kt);
    const Pe = 1 / d, a = Math.sqrt(1 + 4 * kt * d);
    return 4 * a * Math.exp((1 - a) * Pe / 2) / ((1 + a) * (1 + a) - (1 - a) * (1 - a) * Math.exp(-a * Pe));
  };

  /* ---------- time grid with atoms ---------- */
  RT.grid = (tmax, n) => ({ tmax, n, dt: tmax / (n - 1), c: new Float64Array(n), d: [] });
  RT.sample = function (f, tmax, n) { const g = RT.grid(tmax, n); for (let i = 0; i < n; i++) g.c[i] = f(i * g.dt); return g; };
  /* smooth curves that are costly to evaluate: m knots, monotone cubic in between */
  RT.sampleSmooth = function (f, tmax, n, m = 360) {
    const xs = N.linspace(0, tmax, m), ys = xs.map((x, i) => f(i ? x : 1e-7 * tmax)), p = N.pchip(xs, ys), g = RT.grid(tmax, n);
    for (let i = 0; i < n; i++) g.c[i] = Math.max(0, p(i * g.dt));
    return g;
  };
  RT.atom = function (t, w, tmax, n) { const g = RT.grid(tmax, n); if (t <= tmax) g.d.push([t, w]); return g; };
  const like = A => RT.grid(A.tmax, A.n);
  RT.scale = function (A, s) { const B = like(A); for (let i = 0; i < A.n; i++) B.c[i] = A.c[i] * s; B.d = A.d.map(([t, w]) => [t, w * s]); return B; };
  RT.add = function (A, B) { const S = like(A); for (let i = 0; i < A.n; i++) S.c[i] = A.c[i] + B.c[i]; S.d = mergeAtoms(A.d.concat(B.d), A.dt); return S; };
  function mergeAtoms(list, dt) {
    const out = [];
    list.slice().sort((a, b) => a[0] - b[0]).forEach(([t, w]) => {
      if (Math.abs(w) < 1e-14) return;
      const last = out[out.length - 1];
      if (last && Math.abs(last[0] - t) < 0.25 * dt) { const W = last[1] + w; last[0] = (last[0] * last[1] + t * w) / W; last[1] = W; }
      else out.push([t, w]);
    });
    return out;
  }
  const interp = (c, dt, t) => { if (t < 0) return 0; const x = t / dt, i = Math.floor(x); if (i >= c.length - 1) return i === c.length - 1 ? c[i] : 0; const f = x - i; return c[i] * (1 - f) + c[i + 1] * f; };
  RT.at = (A, t) => interp(A.c, A.dt, t);

  /* radix-2 FFT, in place */
  function fft(re, im, inv) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = 2 * Math.PI / len * (inv ? -1 : 1), wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let j = 0; j < len / 2; j++) {
          const ur = re[i + j], ui = im[i + j], k = i + j + len / 2;
          const vr = re[k] * cr - im[k] * ci, vi = re[k] * ci + im[k] * cr;
          re[i + j] = ur + vr; im[i + j] = ui + vi; re[k] = ur - vr; im[k] = ui - vi;
          const nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
        }
      }
    }
    if (inv) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  }
  /* trapezoidal convolution of two sampled densities */
  function convCC(a, b, dt) {
    const n = a.length;
    let m = 1;
    while (m < 2 * n) m <<= 1;
    const ar = new Float64Array(m), ai = new Float64Array(m), br = new Float64Array(m), bi = new Float64Array(m);
    ar.set(a); br.set(b);
    fft(ar, ai, false); fft(br, bi, false);
    for (let i = 0; i < m; i++) { const r = ar[i] * br[i] - ai[i] * bi[i], im = ar[i] * bi[i] + ai[i] * br[i]; ar[i] = r; ai[i] = im; }
    fft(ar, ai, true);
    const out = new Float64Array(n);
    for (let i = 1; i < n; i++) out[i] = (ar[i] - 0.5 * (a[0] * b[i] + a[i] * b[0])) * dt;
    return out;
  }
  RT.conv = function (A, B) {
    const S = like(A), n = A.n, dt = A.dt;
    const anyA = A.c.some(x => x !== 0), anyB = B.c.some(x => x !== 0);
    if (anyA && anyB) S.c = convCC(A.c, B.c, dt);
    const extra = [];
    const shiftAdd = (src, t0, w) => {
      for (let i = 0; i < n; i++) { const t = i * dt - t0; if (t >= 0) S.c[i] += w * interp(src, dt, t); }
      /* keep the area exact when the shifted density starts with a jump */
      if (src[0] !== 0 && t0 > 0) { const i0 = Math.ceil(t0 / dt - 1e-12), del = i0 * dt - t0; extra.push([t0, -w * src[0] * (0.5 * dt - del)]); }
    };
    A.d.forEach(([t, w]) => { if (anyB) shiftAdd(B.c, t, w); });
    B.d.forEach(([t, w]) => { if (anyA) shiftAdd(A.c, t, w); });
    const atoms = [];
    A.d.forEach(([t1, w1]) => B.d.forEach(([t2, w2]) => { if (t1 + t2 <= A.tmax) atoms.push([t1 + t2, w1 * w2]); }));
    S.d = mergeAtoms(atoms.concat(extra), dt);
    return S;
  };
  RT.mass = A => { let s = 0; for (let i = 1; i < A.n; i++) s += 0.5 * (A.c[i] + A.c[i - 1]) * A.dt; return s + N.sum(A.d.map(a => a[1])); };
  /* recycle ratio R around block g: E = Σ (1/(1+R)) (R/(1+R))^(k−1) g^{*k} */
  RT.recycle = function (g, R) {
    if (!(R > 0)) return g;
    const a = 1 / (1 + R), b = R / (1 + R);
    let P = g, E = RT.scale(g, a), wk = a;
    for (let k = 2; k < 400; k++) {
      wk *= b;
      P = RT.conv(P, g);
      const m = RT.mass(P);
      if (wk * m < 1e-7) break;
      E = RT.add(E, RT.scale(P, wk));
    }
    return E;
  };
  RT.moments = function (A) {
    let m0 = 0, m1 = 0, m2 = 0, m3 = 0;
    for (let i = 1; i < A.n; i++) {
      const t0 = (i - 1) * A.dt, t1 = i * A.dt, e0 = A.c[i - 1], e1 = A.c[i];
      const h = A.dt;
      m0 += 0.5 * h * (e0 + e1); m1 += 0.5 * h * (t0 * e0 + t1 * e1); m2 += 0.5 * h * (t0 * t0 * e0 + t1 * t1 * e1); m3 += 0.5 * h * (t0 ** 3 * e0 + t1 ** 3 * e1);
    }
    A.d.forEach(([t, w]) => { m0 += w; m1 += w * t; m2 += w * t * t; m3 += w * t ** 3; });
    const mean = m1 / m0, vr = m2 / m0 - mean * mean;
    const sk = (m3 / m0 - 3 * mean * vr - mean ** 3) / Math.pow(vr, 1.5);
    return { area: m0, mean, var: vr, skew: sk };
  };
  RT.cdf = function (A) {
    const F = new Float64Array(A.n);
    const atoms = A.d.slice().sort((a, b) => a[0] - b[0]);
    let s = 0, k = 0;
    for (let i = 0; i < A.n; i++) {
      if (i) s += 0.5 * (A.c[i] + A.c[i - 1]) * A.dt;
      let extra = 0;
      while (k < atoms.length && atoms[k][0] <= i * A.dt + 1e-12) { extra += atoms[k][1]; k++; }
      s += extra;
      F[i] = s;
    }
    return F;
  };
  /* ∫ f(t) E(t) dt over the grid, plus atoms */
  RT.integrate = function (A, f) {
    let s = 0, f0 = f(0);
    for (let i = 1; i < A.n; i++) { const f1 = f(i * A.dt); s += 0.5 * A.dt * (A.c[i - 1] * f0 + A.c[i] * f1); f0 = f1; }
    A.d.forEach(([t, w]) => { s += w * f(t); });
    return s;
  };

  /* ---------- piecewise-linear E from a table (allows vertical steps) ---------- */
  RT.pl = function (ts, es) {
    const segs = [];
    for (let i = 1; i < ts.length; i++) if (ts[i] > ts[i - 1]) segs.push([ts[i - 1], ts[i], es[i - 1], es[i]]);
    const E = t => { for (const s of segs) if (t >= s[0] && t <= s[1]) return s[2] + (s[3] - s[2]) * (t - s[0]) / (s[1] - s[0]); return 0; };
    const integ = (f, rtol) => { let s = 0; segs.forEach(sg => { if (sg[2] === 0 && sg[3] === 0) return; s += N.quad(t => f(t) * (sg[2] + (sg[3] - sg[2]) * (t - sg[0]) / (sg[1] - sg[0])), sg[0], sg[1], { rtol: rtol || 1e-9 }); }); return s; };
    const area = N.sum(segs.map(s => 0.5 * (s[2] + s[3]) * (s[1] - s[0])));
    return { segs, E, integ, area, tmax: ts[ts.length - 1] };
  };

  /* ---------- conversion with an RTD ---------- */
  /* batch history X(t) for a constant-density single-reaction model (ST.K) */
  RT.batchX = function (m, tEnd, X0 = 0) {
    const cap = Math.min(m.Xe, m.Xmax) * (1 - 1e-12);
    if (!(tEnd > 0) || X0 >= cap) return () => Math.min(X0, cap);
    const sol = N.ode((t, y) => [Math.max(0, m.rate(Math.min(y[0], cap))) / m.CA0], 0, [X0], tEnd, { rtol: 1e-10, atol: 1e-13, minPts: 200, stop: (t, y) => y[0] >= cap });
    const tl = sol.t[sol.t.length - 1], Xl = Math.min(sol.last()[0], cap);
    return t => (t <= 0 ? X0 : t >= tl ? Xl : Math.min(sol.at(t)[0], cap));
  };
  /* Zwietering's maximum mixedness: dX/dλ = Λ(λ)X − r(X)/C_A0 integrated from large λ down to 0 */
  RT.maxMix = function (m, Lam, lmax) {
    const cap = Math.min(m.Xe, m.Xmax) * (1 - 1e-12);
    const L0 = Lam(lmax);
    let X0 = 0;
    if (L0 > 0) {
      const g = x => m.CA0 * L0 * x - m.rate(x);
      const r = N.roots(g, 0, cap, 200, { tol: 1e-13 });
      X0 = r.length ? r[r.length - 1] : 0;
    }
    const sol = N.ode((u, y) => { const x = Math.min(Math.max(y[0], 0), cap); return [m.rate(x) / m.CA0 - Lam(lmax - u) * x]; }, 0, [X0], lmax, { rtol: 1e-8, atol: 1e-11, minPts: 200 });
    return Math.min(Math.max(sol.last()[0], 0), cap);
  };

  /* Axial dispersion with reaction, closed vessel (Danckwerts conditions), any single-reaction kinetics:
       d·X'' − X' + (τ/C_A0)(−r_A) = 0,   X(0) = d·X'(0),   X'(1) = 0.
     Finite differences (central, upwinded where the cell Péclet number is large) and damped Newton
     on the tridiagonal system. Returns {X (exit), z, x (profile)}. */
  RT.dispersionX = function (m, tau, d, o = {}) {
    const cap = Math.min(m.Xe, m.Xmax) * (1 - 1e-12);
    const q = tau / m.CA0;
    const M = Math.min(20000, Math.max(o.M || 400, Math.ceil(4 / Math.max(d, 1e-9))));
    const h = 1 / M;
    const up = h / Math.max(d, 1e-300) > 2;
    const r = x => m.rate(Math.min(Math.max(x, 0), cap));
    /* initial guess: plug flow profile */
    const Xin = o.Xin || 0;
    let X = new Float64Array(M + 1);
    if (o.guess && o.guess.length === M + 1) X.set(o.guess);
    else { const pf = RT.batchX(m, tau, Xin); for (let i = 0; i <= M; i++) X[i] = Math.min(pf(i * h * tau), cap); }
    const res = new Float64Array(M + 1), a = new Float64Array(M + 1), b = new Float64Array(M + 1), c = new Float64Array(M + 1);
    const build = () => {
      for (let i = 0; i <= M; i++) {
        let xm = i > 0 ? X[i - 1] : X[1] - 2 * h * (X[0] - Xin) / d;
        const xp = i < M ? X[i + 1] : X[M - 1];
        const ri = r(X[i]), dr = (r(X[i] + 1e-7) - r(X[i] - 1e-7)) / 2e-7;
        let lo, di, hi;
        if (up) {
          res[i] = d * (xp - 2 * X[i] + xm) / (h * h) - (X[i] - xm) / h + q * ri;
          lo = d / (h * h) + 1 / h; di = -2 * d / (h * h) - 1 / h + q * dr; hi = d / (h * h);
        } else {
          res[i] = d * (xp - 2 * X[i] + xm) / (h * h) - (xp - xm) / (2 * h) + q * ri;
          lo = d / (h * h) + 1 / (2 * h); di = -2 * d / (h * h) + q * dr; hi = d / (h * h) - 1 / (2 * h);
        }
        if (i === 0) { /* ghost: X_-1 = X_1 − 2hX_0/d */ di += lo * (-2 * h / d); hi += lo; lo = 0; }
        if (i === M) { lo += hi; hi = 0; }
        a[i] = lo; b[i] = di; c[i] = hi;
      }
    };
    let ok = false;
    for (let it = 0; it < 60; it++) {
      build();
      let nr = 0;
      for (let i = 0; i <= M; i++) nr = Math.max(nr, Math.abs(res[i]));
      if (nr < 1e-10 * (1 + q * Math.abs(r(0)))) { ok = true; break; }
      /* Thomas algorithm for J·dx = −res */
      const cp = new Float64Array(M + 1), dp = new Float64Array(M + 1);
      cp[0] = c[0] / b[0]; dp[0] = -res[0] / b[0];
      for (let i = 1; i <= M; i++) { const den = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / den; dp[i] = (-res[i] - a[i] * dp[i - 1]) / den; }
      const dx = new Float64Array(M + 1);
      dx[M] = dp[M];
      for (let i = M - 1; i >= 0; i--) dx[i] = dp[i] - cp[i] * dx[i + 1];
      let lam = 1, mx = 0;
      for (let i = 0; i <= M; i++) mx = Math.max(mx, Math.abs(dx[i]));
      if (mx > 0.25) lam = 0.25 / mx;
      for (let i = 0; i <= M; i++) X[i] = Math.min(Math.max(X[i] + lam * dx[i], -0.5), cap);
      if (mx * lam < 1e-13) { ok = true; break; }
    }
    const step = Math.max(1, Math.floor(M / 200));
    const z = [], x = [];
    for (let i = 0; i <= M; i += step) { z.push(i * h); x.push(X[i]); }
    if (z[z.length - 1] !== 1) { z.push(1); x.push(X[M]); }
    return { X: X[M], X0: X[0], z, x, ok, prof: X };
  };

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
