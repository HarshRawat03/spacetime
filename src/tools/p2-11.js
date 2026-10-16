/* Chapter 11 — Basics of non-ideal flow: the RTD from tracer tests, convolution, conversion from the RTD */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, RT = ST.RT, KIT = ST.KIT, K = ST.K, D = ST.D;
  const RTDI = ST.RTDI = {};

  /* ---------- tracer data ---------- */
  RTDI.rows = function (rows, what) {
    const r = (rows || []).filter(q => q && isFinite(q[0]) && isFinite(q[1])).map(q => [+q[0], +q[1]]);
    if (r.length < 2) throw new Error('Enter at least two rows of ' + (what || 'data') + '.');
    for (let i = 1; i < r.length; i++) if (r[i][0] < r[i - 1][0]) throw new Error('Times must not decrease (row ' + (i + 1) + ').');
    return r;
  };
  /* exact integrals of t^k·C(t) for a piecewise-linear C */
  function plMoments(t, c) {
    const m = [0, 0, 0, 0];
    for (let i = 1; i < t.length; i++) {
      const a = t[i - 1], b = t[i];
      if (b <= a) continue;
      const s = (c[i] - c[i - 1]) / (b - a), c0 = c[i - 1] - s * a;
      for (let k = 0; k < 4; k++) m[k] += c0 * (b ** (k + 1) - a ** (k + 1)) / (k + 1) + s * (b ** (k + 2) - a ** (k + 2)) / (k + 2);
    }
    return m;
  }
  function sumMoments(t, c) {
    const m = [0, 0, 0, 0], n = t.length;
    for (let i = 0; i < n; i++) {
      const dt = n === 1 ? 1 : i === 0 ? (t[1] - t[0]) / 2 : i === n - 1 ? (t[n - 1] - t[n - 2]) / 2 : (t[i + 1] - t[i - 1]) / 2;
      const w = (i === 0 || i === n - 1) && c[i] !== 0 ? dt : dt;
      for (let k = 0; k < 4; k++) m[k] += c[i] * Math.pow(t[i], k) * w;
    }
    return m;
  }
  /* exponential tail fitted to the last points */
  function tailFit(t, c) {
    const n = t.length, pts = [];
    for (let i = n - 1; i >= 0 && pts.length < 4; i--) { if (c[i] > 0) pts.unshift([t[i], Math.log(c[i])]); else if (pts.length) break; }
    if (pts.length < 3) return null;
    const lr = N.linreg(pts.map(p => p[0]), pts.map(p => p[1]));
    const slope = lr.slope !== undefined ? lr.slope : lr.b, icpt = lr.intercept !== undefined ? lr.intercept : lr.a;
    if (!(slope < 0)) return null;
    const tau = -1 / slope, tl = t[n - 1], cl = c[n - 1];
    if (!(cl > 0)) return null;
    /* ∫_{tl}^{∞} t^k cl e^{−(t−tl)/τ} dt */
    const m = [cl * tau, cl * tau * (tl + tau), cl * tau * (tl * tl + 2 * tl * tau + 2 * tau * tau), cl * tau * (tl ** 3 + 3 * tl * tl * tau + 6 * tl * tau * tau + 6 * tau ** 3)];
    return { tau, tl, cl, m, icpt };
  }
  /* analyse a pulse curve (t, C): returns E as piecewise-linear plus moments */
  RTDI.pulse = function (rows, o = {}) {
    const t = rows.map(r => r[0]), c = rows.map(r => r[1]);
    let m = o.rule === 'sum' ? sumMoments(t, c) : plMoments(t, c);
    let tail = null;
    if (o.tail) { tail = tailFit(t, c); if (tail) m = m.map((x, k) => x + tail.m[k]); }
    const area = m[0], mean = m[1] / area, vr = m[2] / area - mean * mean;
    const skew = (m[3] / area - 3 * mean * vr - mean ** 3) / Math.pow(vr, 1.5);
    const E = c.map(x => x / area);
    const pl = RT.pl(t, E);
    return { t, c, E, area, mean, var: vr, skew, tail, pl };
  };
  /* analyse a step curve (t, C): F from C/Cmax; moments from 1 − F */
  RTDI.step = function (rows, o = {}) {
    const t = rows.map(r => r[0]), c = rows.map(r => r[1]);
    const Cmax = o.Cmax > 0 ? o.Cmax : Math.max(...c);
    const F = c.map(x => x / Cmax);
    let I0 = 0, I1 = 0, I2 = 0;
    for (let i = 1; i < t.length; i++) {
      const a = t[i - 1], b = t[i], ga = 1 - F[i - 1], gb = 1 - F[i];
      I0 += 0.5 * (ga + gb) * (b - a); I1 += 0.5 * (a * ga + b * gb) * (b - a); I2 += 0.5 * (a * a * ga + b * b * gb) * (b - a);
    }
    if (t[0] > 0) { I0 += t[0]; I1 += t[0] * t[0] / 2; I2 += t[0] ** 3 / 3; }
    const mean = I0, vr = 2 * I1 - mean * mean;
    const E = t.map((x, i) => { const a = Math.max(0, i - 1), b = Math.min(t.length - 1, i + 1); return t[b] > t[a] ? Math.max(0, (F[b] - F[a]) / (t[b] - t[a])) : 0; });
    const third = 3 * I2;
    const skew = (third - 3 * mean * vr - mean ** 3) / Math.pow(vr, 1.5);
    return { t, c, F, E, Cmax, area: NaN, mean, var: vr, skew, pl: RT.pl(t, E) };
  };
  /* closed-vessel dispersion number from a dimensionless variance */
  RTDI.dFromVar = function (s2, bc) {
    if (bc === 'open') return s2 > 0 ? (-2 + Math.sqrt(4 + 32 * s2)) / 16 : NaN;
    if (bc === 'small') return s2 / 2;
    if (!(s2 > 0) || s2 >= 1) return NaN;
    return N.brent(d => RT.varClosed(d) - s2, 1e-9, 1e6, { tol: 1e-14 });
  };
  const VOLU = { L: 1e-3, 'm³': 1 };
  const FLOWU = { 'L/s': [1e-3, 1], 'L/min': [1e-3, 60], 'L/h': [1e-3, 3600], 'm³/s': [1, 1], 'm³/min': [1, 60], 'm³/h': [1, 3600] };
  const TSEC = { s: 1, min: 60, h: 3600, day: 86400 };
  const CVOL = { 'g/L': 'L', 'mol/L': 'L', 'mg/L': 'L', 'kg/m³': 'm³', 'mol/m³': 'm³', 'g/m³': 'm³' };
  RTDI.flowIn = (v, vu, volu, tu) => { const [vf, ts] = FLOWU[vu]; return v * vf / ts * TSEC[tu] / VOLU[volu]; };

  /* ===================================================================== */
  /* F-1101  the RTD from a tracer test                                     */
  /* ===================================================================== */
  const EX111 = [[0, 0], [5, 3], [10, 5], [15, 5], [20, 4], [25, 2], [30, 1], [35, 0]];
  const EX112 = N.range(10).map(k => [2 * (k + 1), 0.375 * Math.pow(0.25, k)]);
  ST.tool({
    id: 'rtddata', tag: 'F-1101', ch: 11, sec: 'Ch. 11.1', title: 'The RTD from a tracer test: E, F, moments and a diagnosis',
    lede: 'Paste a pulse or step response, or a train of ideal peaks. Out come the material-balance check, mean, variance and skewness, the E, F and E_θ curves, the intensity function, the flowing volume and its share of the vessel, equivalent tanks and dispersion number, and what the curve says about the flow inside.',
    keys: 'residence time distribution rtd tracer pulse step e curve f curve mean variance skewness material balance dead volume bypass channeling recirculation intensity function internal age e theta example 11.1 11.2',
    inputs: [
      { group: 'Tracer response' },
      { k: 'kind', t: 'seg', label: 'Experiment', opts: [['pulse', 'pulse'], ['step', 'step'], ['peaks', 'ideal peaks']], v: 'pulse' },
      { k: 'data', t: 'table', label: 'Outlet readings', cols: ['t', 'C'], v: EX111, hint: 'For ideal peaks, give each peak’s time and area.' },
      { k: 'tu', t: 'select', label: 'Time unit of the data', opts: [['s', 's'], ['min', 'min'], ['h', 'h'], ['day', 'day']], v: 'min' },
      { k: 'cu', t: 'select', label: 'Concentration unit', opts: [['g/L', 'g/L'], ['mg/L', 'mg/L'], ['mol/L', 'mol/L'], ['kg/m³', 'kg/m³'], ['mol/m³', 'mol/m³'], ['g/m³', 'g/m³'], ['arb', 'arbitrary']], v: 'g/L' },
      { k: 'rule', t: 'seg', label: 'Integrate', opts: [['pl', 'piecewise linear'], ['sum', 'Σ C Δt (book)']], v: 'pl', if: v => v.kind === 'pulse' },
      { k: 'tail', t: 'check', label: 'Extrapolate an exponential tail', v: false, if: v => v.kind === 'pulse' },
      { k: 'Cmax', label: 'Step height (blank or 0 = largest reading)', v: 0, min: 0, max: 1e9, if: v => v.kind === 'step' },
      { group: 'Flow and tracer (optional)' },
      { k: 'flowOn', t: 'check', label: 'I know the flow rate', v: true },
      { k: 'vf', label: 'Flow rate v', v: 1, min: 1e-12, max: 1e9, scale: 'log', if: v => v.flowOn },
      { k: 'vu', t: 'select', label: 'Flow unit', opts: Object.keys(FLOWU).map(u => [u, u]), v: 'L/s', if: v => v.flowOn },
      { k: 'M', label: 'Tracer injected, M (same mass unit as C)', v: 0, min: 0, max: 1e12, if: v => v.flowOn && v.kind !== 'step', hint: 'Zero to skip the material-balance check.' },
      { k: 'Vv', label: 'Vessel volume (0 if unknown)', v: 0, min: 0, max: 1e12, if: v => v.flowOn },
      { k: 'Vu', t: 'select', label: 'Volume unit', opts: [['L', 'L'], ['m³', 'm³']], v: 'L', if: v => v.flowOn }
    ],
    presets: [
      { name: 'Example 11.1', v: { kind: 'pulse', data: EX111, tu: 'min', cu: 'g/L', rule: 'sum', tail: false, flowOn: false }, book: { tbar: 15, area: 100, E10: 0.05 }, note: 'Pulse response of a closed vessel: t̄ = 15 min and E = C/100 per minute.' },
      { name: 'Example 11.2', v: { kind: 'peaks', data: EX112, tu: 'min', cu: 'g/L', flowOn: true, vf: 5, vu: 'L/s', M: 150, Vv: 860, Vu: 'L' }, book: { area: 0.5, tbar: { v: 2.67, tol: 0.01 }, Vflow: { v: 800, tol: 0.01 }, frac: { v: 0.93, tol: 0.01 } }, note: 'Gas–liquid contactor, 860 L, liquid at 5 L/s, 150 g of tracer: peaks every 2 min, each a quarter of the last. Liquid fills 800 L (93%) and recirculates strongly.' },
      { name: 'Problem 14.10 (baffled reactor)', v: { kind: 'pulse', data: [[0, 0], [10, 35], [20, 38], [30, 40], [40, 40], [50, 39], [60, 37], [70, 36], [80, 35], [90, 0]], tu: 'min', cu: 'arb', rule: 'pl', tail: false, flowOn: false }, note: 'A long flat response: much more spread than the mean suggests.' },
      { name: 'Problem 14.11 (non-normalized)', v: { kind: 'pulse', data: [[0, 0], [1, 9], [2, 57], [3, 81], [4, 90], [5, 90], [6, 86], [8, 77], [10, 67], [15, 47], [20, 32], [30, 15], [41, 7], [52, 3], [67, 1], [70, 0]], tu: 'min', cu: 'arb', rule: 'pl', tail: false, flowOn: false }, note: 'Fit tanks-in-series or dispersion to this reactor, then use F-1103 for conversion.' },
      { name: 'Step response', v: { kind: 'step', data: [[0, 0], [1, 0.02], [2, 0.12], [3, 0.3], [4, 0.5], [5, 0.66], [6, 0.78], [8, 0.91], [10, 0.966], [12, 0.988], [15, 0.998], [20, 1]], tu: 'min', cu: 'arb', Cmax: 0, flowOn: true, vf: 1, vu: 'm³/min', Vv: 5, Vu: 'm³' }, note: 'A step switched on at t = 0; F = C/C_max, and the mean comes from the area above the F curve.' }
    ],
    eqs: [
      ['\\int_0^\\infty E\\,dt = 1,\\qquad E = \\frac{C_{pulse}}{M/v}', 'Eqs. 11.1, 11.5'],
      ['\\bar t = \\frac{\\int tC\\,dt}{\\int C\\,dt} \\approx \\frac{\\sum t_iC_i\\Delta t_i}{\\sum C_i\\Delta t_i} = \\frac{V}{v}', 'Eq. 11.4'],
      ['\\sigma^2 = \\frac{\\int (t-\\bar t)^2 C\\,dt}{\\int C\\,dt}', 'Eq. 13.2'],
      ['F = \\frac{C_{step}}{C_{max}} = \\int_0^t E\\,dt,\\qquad \\bar t = \\int_0^\\infty (1-F)\\,dt', 'Eqs. 11.6–11.8'],
      ['E_\\theta = \\bar t E,\\quad \\theta = t/\\bar t', 'Eq. 11.6'],
      ['\\Lambda(t) = \\frac{E}{1-F},\\qquad I(t) = \\frac{1-F}{\\bar t}', '', 'Intensity function and internal age distribution (diagnostics).']
    ],
    method: ['Areas and moments are exact integrals of the piecewise-linear curve through the readings, or the book’s Σ C Δt sums. An optional exponential tail is fitted to the last three or four readings and integrated analytically.', 'Step data are read through F; the mean and variance come from areas above the F curve, which is more robust than differentiating. Dispersion numbers come from the closed-vessel variance relation solved exactly.'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const rows = RTDI.rows(v.data, 'readings');
      const tu = v.tu;
      let R;
      if (v.kind === 'peaks') {
        const area = N.sum(rows.map(r => r[1]));
        const mean = N.sum(rows.map(r => r[0] * r[1])) / area;
        const vr = N.sum(rows.map(r => r[1] * Math.pow(r[0] - mean, 2))) / area;
        const sk = N.sum(rows.map(r => r[1] * Math.pow(r[0] - mean, 3))) / area / Math.pow(vr, 1.5);
        R = { t: rows.map(r => r[0]), c: rows.map(r => r[1]), area, mean, var: vr, skew: sk, peaks: true };
      } else if (v.kind === 'step') R = RTDI.step(rows, { Cmax: +v.Cmax });
      else R = RTDI.pulse(rows, { rule: v.rule, tail: v.tail });
      const tb = R.mean, s2 = R.var, s2t = s2 / (tb * tb);
      if (!(tb > 0)) throw new Error('The mean time came out zero or negative: check the readings.');
      const cu = v.cu === 'arb' ? '' : v.cu;
      if (v.kind !== 'step') kpi.push({ k: 'area', l: 'Area under the curve', v: R.area, u: cu ? cu.split('/')[0] + '·' + tu + '/' + cu.split('/')[1] : '' });
      kpi.push({ k: 'tbar', l: 'Mean residence time t̄', v: tb, u: tu, main: true }, { k: 'var', l: 'Variance σ²', v: s2, u: tu + '²', s: 'σ = ' + N.fmt(Math.sqrt(s2), 4) + ' ' + tu }, { k: 's2t', l: 'σ²_θ = σ²/t̄²', v: s2t });
      if (isFinite(R.skew)) kpi.push({ k: 'skew', l: 'Skewness', v: R.skew, s: 'mixed flow: 2; gaussian: 0' });
      kpi.push({ k: 'Neq', l: 'Equivalent tanks in series', v: 1 / s2t, s: 'N = t̄²/σ²' });
      const dC = RTDI.dFromVar(s2t, 'closed');
      kpi.push({ k: 'dC', l: 'D/uL, closed vessel', v: dC, s: isFinite(dC) ? (dC < 0.01 ? 'small deviation from plug' : dC > 1 ? 'large: question the model' : 'intermediate') : 'σ²_θ ≥ 1: no dispersion model fits' });
      if (R.tail) info.push('Tail fitted from the last readings: time constant ' + N.fmt(R.tail.tau, 4) + ' ' + tu + '; it adds ' + N.fmt(100 * R.tail.m[0] / R.area, 3) + '% to the area.');
      let Ev = null;
      if (v.flowOn) {
        const volu = CVOL[v.cu] || v.Vu;
        if (v.kind !== 'step' && +v.M > 0 && v.cu !== 'arb') {
          const vIn = RTDI.flowIn(+v.vf, v.vu, volu, tu);
          const Mv = +v.M / vIn;
          const err = (R.area - Mv) / Mv;
          kpi.push({ k: 'Mv', l: 'M/v from the injection', v: Mv, u: cu.split('/')[0] + '·' + tu + '/' + cu.split('/')[1], s: Math.abs(err) < 0.05 ? 'agrees with the area (' + N.fmt(100 * err, 2) + '%)' : 'disagrees with the area by ' + N.fmt(100 * err, 3) + '%', state: Math.abs(err) < 0.05 ? undefined : 'bad' });
          if (Math.abs(err) >= 0.05) warn.push('Material balance does not close: tracer lost (adsorbed, reacted, short sampling) or flow/M mismeasured. Moments may still be useful, but treat V and E with care.');
        }
        const vV = RTDI.flowIn(+v.vf, v.vu, v.Vu, tu);
        const Vf = vV * tb;
        kpi.push({ k: 'Vflow', l: 'Volume swept by the flow, v·t̄', v: Vf, u: v.Vu, main: true });
        if (+v.Vv > 0) {
          const fr = Vf / +v.Vv;
          Ev = fr;
          kpi.push({ k: 'frac', l: 'Share of the vessel volume', v: fr, s: fr < 0.97 ? N.fmt(100 * (1 - fr), 3) + '% not swept (dead water or another phase)' : fr > 1.03 ? 'more than the vessel: check flow meter or tracer' : 'all of it', state: fr > 1.03 ? 'bad' : undefined });
          if (fr > 1.03) warn.push('Late tracer: t̄ exceeds V/v. The material balance says this cannot happen, so the flow rate or volume is wrong, the tracer is not inert, or the vessel is far from closed (Fig. 12.3).');
          if (fr < 0.9) info.push('Early mean: about ' + N.fmt(100 * (1 - fr), 3) + '% of the vessel is not swept — stagnant zones, or space taken by a second phase.');
        }
      }
      /* shape diagnosis */
      let peaksN = 0, early = false;
      if (!R.peaks) {
        const E = R.E;
        const mx = Math.max(...E);
        for (let i = 1; i < E.length - 1; i++) if (E[i] > E[i - 1] && E[i] >= E[i + 1] && E[i] > 0.15 * mx) peaksN++;
        const ip = E.indexOf(mx), thp = R.t[ip] / tb;
        early = thp < 0.35 && s2t < 0.9 && E[0] < 0.5 * mx;
        if (peaksN > 1) info.push(peaksN + ' peaks: regularly spaced peaks mean internal recirculation, two peaks mean parallel paths or channeling (Fig. 12.3).');
        if (early) info.push('The peak arrives at θ = ' + N.fmt(thp, 3) + ', well before the mean: a fast channel or bypass plus a slow, long-tailed region.');
        if (R.skew > 2.5) info.push('Strongly skewed (skewness ' + N.fmt(R.skew, 3) + '): a long tail from slow exchange with stagnant pockets.');
        if (s2t > 1.05) info.push('σ²_θ > 1: wider than mixed flow — think bypass with dead water, or recirculation through a slow loop.');
        if (s2t < 0.05 && peaksN <= 1) info.push('A narrow, single peak: close to plug flow (D/uL ' + N.fmt(dC, 3) + ').');
      } else info.push('Ideal peaks at regular intervals, each a fraction of the last: liquid recirculating through the vessel.');
      /* plots */
      const plots = [];
      if (R.peaks) {
        const area = R.area;
        plots.push({ title: 'Ideal peaks (areas)', x: { label: 't', unit: tu, min: 0 }, y: { label: 'peak area', min: 0 }, series: [{ name: 'peaks', x: R.t, y: R.c, type: 'bar', color: 1 }] });
        const F = []; let s = 0; R.c.forEach(a => { s += a / area; F.push(s); });
        plots.push({ title: 'F curve', x: { label: 't', unit: tu, min: 0 }, y: { label: 'F', min: 0, max: 1 }, series: [{ name: 'F', x: [0].concat(R.t), y: [0].concat(F), type: 'step', color: 1 }] });
      } else {
        const ts = R.t, tmax = ts[ts.length - 1];
        const fine = N.linspace(0, tmax * (R.tail ? 1.5 : 1), 300);
        const Efun = t => (t <= tmax ? R.pl.E(t) : R.tail ? R.tail.cl * Math.exp(-(t - R.tail.tl) / R.tail.tau) / R.area : 0);
        const Ffun = (() => { const ys = fine.map(Efun), cum = N.cumtrapz(fine, ys); return t => N.interp(fine, cum, t); })();
        if (v.kind === 'step') {
          plots.push({ title: 'Step response, F = C/C_max', x: { label: 't', unit: tu, min: 0 }, y: { label: 'F', min: 0, max: 1.05 }, series: [{ name: 'F', x: ts, y: R.F, color: 1, marker: 'o' }], marks: [{ t: 'v', x: tb, color: 8, label: 't̄' }] });
          plots.push({ title: 'E = dF/dt', x: { label: 't', unit: tu, min: 0 }, y: { label: 'E', unit: '1/' + tu, min: 0 }, series: [{ name: 'E', x: ts, y: R.E, color: 2, marker: 'o' }] });
        } else {
          plots.push({ title: 'Exit age distribution E (Fig. 11.6)', x: { label: 't', unit: tu, min: 0 }, y: { label: 'E', unit: '1/' + tu, min: 0 }, series: [{ name: 'E from the data', x: fine, y: fine.map(Efun), color: 1 }, { x: fine, y: fine.map(Efun), type: 'area', fill: 0.15, color: 1, hover: false, legend: false }, { name: 'readings / area', x: ts, y: R.E, type: 'pts', color: 1 }], marks: [{ t: 'v', x: tb, color: 2, label: 't̄ = ' + N.fmt(tb, 4) }] });
          plots.push({ title: 'F curve', x: { label: 't', unit: tu, min: 0 }, y: { label: 'F', min: 0, max: 1 }, series: [{ name: 'F', x: fine, y: fine.map(Ffun), color: 1 }] });
        }
        /* E_theta with models of equal variance */
        const th = N.linspace(0, 3, 300);
        const ser = [];
        if (v.kind !== 'step') ser.push({ name: 'data', x: fine.map(t => t / tb), y: fine.map(t => tb * Efun(t)), color: 1, w: 2.2 });
        else ser.push({ name: 'data', x: ts.map(t => t / tb), y: R.E.map(e => e * tb), color: 1, marker: 'o' });
        const Nq = 1 / s2t;
        if (Nq >= 0.5 && Nq < 400) ser.push({ name: 'N = ' + N.fmt(Nq, 3) + ' tanks', x: th, y: th.map(q => RT.Etanks(q, 1, Nq)), color: 3, dash: '5 4' });
        if (isFinite(dC)) ser.push({ name: 'dispersion, D/uL = ' + N.fmt(dC, 3), x: th, y: th.map(q => RT.EthClosed(q, dC)), color: 4, dash: '2 3' });
        ser.push({ name: 'mixed flow', x: th, y: th.map(q => Math.exp(-q)), color: 8, dash: '1 3' });
        plots.push({ title: 'E_θ against θ, with models of the same variance', x: { label: 'θ = t/t̄', min: 0, max: 3 }, y: { label: 'E_θ', min: 0 }, series: ser, cap: 'Equal variance does not mean equal shape: if a model curve sits far from the data, that model should not be used (Ch. 13, comment d).' });
        if (v.kind !== 'step') {
          const tt = fine.filter(t => Ffun(t) < 0.995);
          plots.push({ title: 'Intensity function Λ = E/(1−F)', x: { label: 'θ', min: 0 }, y: { label: 'Λ·t̄' , min: 0 }, series: [{ name: 'data', x: tt.map(t => t / tb), y: tt.map(t => tb * Efun(t) / Math.max(1 - Ffun(t), 1e-9)), color: 1 }, { name: 'mixed flow', x: [0, 3], y: [1, 1], color: 8, dash: '4 4' }], cap: 'Λ rising to a peak and falling again is the signature of stagnancy; a monotone rise is typical of near-plug flow.' });
        }
      }
      const tables = v.kind === 'peaks'
        ? [{ title: 'Peaks', cols: ['t, ' + tu, 'area', 'share of the tracer', 'θ'], rows: R.t.map((t, i) => [t, R.c[i], R.c[i] / R.area, t / tb]) }]
        : [{ title: 'Readings and E', cols: ['t, ' + tu, 'C', v.kind === 'step' ? 'F' : 'E', 'θ', 'E_θ'], rows: R.t.map((t, i) => [t, R.c[i], v.kind === 'step' ? R.F[i] : R.E[i], t / tb, R.E[i] * tb]) }];
      if (v.kind === 'pulse') { const at10 = R.t.indexOf(10); if (at10 >= 0) kpi.push({ k: 'E10', l: 'E at t = 10 ' + tu, v: R.E[at10], u: '1/' + tu }); }
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* F-1102  convolution and deconvolution                                  */
  /* ===================================================================== */
  const CIN113 = [[0, 0], [1, 0], [2, 8], [3, 4], [4, 6], [5, 0]];
  const E113 = [[5, 0], [6, 0.05], [7, 0.5], [8, 0.35], [9, 0.1], [10, 0]];
  ST.tool({
    id: 'convolve', tag: 'F-1102', ch: 11, sec: 'Ch. 11.1', title: 'Convolution and deconvolution of tracer signals',
    lede: 'Push any input signal through a vessel of known E, or chain several vessels, with the book’s discrete sums or a continuous interpretation. Going the other way, recover a vessel’s E from sloppy input and output signals: by moments for a one-parameter model (Example 14.4), by least squares on the whole curve, or non-parametrically.',
    keys: 'convolution integral deconvolution tracer input output signal e curve sloppy pulse example 11.3 example 14.4 additivity of variances',
    inputs: [
      { group: 'What to find' },
      { k: 'mode', t: 'seg', label: 'Find', opts: [['out', 'output from input and E'], ['E', 'E from input and output']], v: 'out' },
      { k: 'cin', t: 'table', label: 'Input signal C_in', cols: ['t', 'C_in'], v: CIN113 },
      { k: 'Esrc', t: 'seg', label: 'Vessel E given as', opts: [['table', 'a table'], ['model', 'a model']], v: 'table', if: v => v.mode === 'out' },
      { k: 'etab', t: 'table', label: 'E of the vessel', cols: ['t', 'E'], v: E113, if: v => v.mode === 'out' && v.Esrc === 'table' },
      { k: 'model', t: 'select', label: 'Vessel model', opts: [['mixed', 'mixed flow'], ['tanks', 'N tanks in series'], ['disp', 'dispersion, closed'], ['plugmix', 'plug + mixed in series'], ['lam', 'laminar flow']], v: 'tanks', if: v => v.mode === 'out' && v.Esrc === 'model' },
      { k: 'mtau', label: 'Mean time of the vessel', unit: '[t]', v: 4, min: 1e-9, max: 1e9, scale: 'log', if: v => v.mode === 'out' && v.Esrc === 'model' },
      { k: 'mN', label: 'N (tanks)', v: 4, min: 0.5, max: 500, scale: 'log', if: v => v.mode === 'out' && v.Esrc === 'model' && v.model === 'tanks' },
      { k: 'md', label: 'D/uL', v: 0.05, min: 1e-5, max: 50, scale: 'log', if: v => v.mode === 'out' && v.Esrc === 'model' && v.model === 'disp' },
      { k: 'mfp', label: 'Plug-flow share of the mean time', v: 0.5, min: 0, max: 0.999, if: v => v.mode === 'out' && v.Esrc === 'model' && v.model === 'plugmix' },
      { k: 'cout', t: 'table', label: 'Output signal C_out', cols: ['t', 'C_out'], v: [], if: v => v.mode === 'E' },
      { k: 'fitm', t: 'select', label: 'Model for E', opts: [['tanks', 'tanks in series'], ['disp', 'dispersion (closed)'], ['free', 'no model: smooth E (regularized)']], v: 'tanks', if: v => v.mode === 'E' },
      { k: 'how', t: 'seg', label: 'Fit by', opts: [['mom', 'moments (Δt̄, Δσ²)'], ['ls', 'least squares on the output']], v: 'mom', if: v => v.mode === 'E' && v.fitm !== 'free' },
      { k: 'reg', label: 'Smoothing weight', v: 1, min: 1e-4, max: 1e4, scale: 'log', if: v => v.mode === 'E' && v.fitm === 'free' },
      { group: 'Arithmetic' },
      { k: 'disc', t: 'seg', label: 'Treat signals as', opts: [['disc', 'discrete slices (book)'], ['cont', 'continuous curves']], v: 'disc', if: v => v.mode === 'out' },
      { k: 'tu', t: 'select', label: 'Time unit', opts: [['s', 's'], ['min', 'min'], ['h', 'h']], v: 'min' }
    ],
    presets: [
      { name: 'Example 11.3', v: { mode: 'out', cin: CIN113, Esrc: 'table', etab: E113, disc: 'disc', tu: 'min' }, book: { Cout9: 4.2, Cout10: 5.1, Cout11: 5.2 }, note: 'Input of 8, 4, 6 at 2, 3, 4 min through a vessel with E = 0.05, 0.5, 0.35, 0.1 at 6–9 min. The tracer leaves between 8 and 13 min, total area 18.' },
      { name: 'Example 14.4 (sloppy input)', v: { mode: 'E', cin: [[200, 0], [203, 1], [237, 1], [240, 0]], cout: [], fitm: 'tanks', how: 'mom', tu: 's', _ex144: true }, book: { N: { v: 4, tol: 0.01 } }, note: 'Input with t̄ = 220 s, σ² = 100 s²; output with t̄ = 280 s, σ² = 1000 s². For tanks in series N = (Δt̄)²/Δσ² = 60²/900 = 4.' },
      { name: 'Chain: input through 4 tanks', v: { mode: 'out', cin: CIN113, Esrc: 'model', model: 'tanks', mtau: 4, mN: 4, disc: 'cont', tu: 'min' }, note: 'The same input through four tanks with a 4-minute mean: means add, variances add.' }
    ],
    eqs: [
      ['C_{out}(t) = \\int_0^t C_{in}(t-t\')\\,E(t\')\\,dt\' = C_{in} * E', 'Eq. 11.10'],
      ['\\bar t_{out} - \\bar t_{in} = \\bar t_{vessel},\\qquad \\sigma^2_{out} - \\sigma^2_{in} = \\sigma^2_{vessel}', 'Eq. 13.9, 14.4'],
      ['N = \\frac{(\\Delta \\bar t)^2}{\\Delta\\sigma^2}', 'Ex. 14.4'],
      ['\\min_E\\ \\|C_{in} * E - C_{out}\\|^2 + \\lambda\\|D^2E\\|^2,\\quad E\\ge 0', '', 'Regularized deconvolution when no model is assumed.']
    ],
    method: ['Discrete mode reproduces the book’s slice-by-slice sums. Continuous mode treats both signals as piecewise-linear curves and convolves them on a fine grid by FFT.', 'Deconvolution by moments uses the additivity of means and variances; least squares fits the model parameter and mean time to the whole output curve with Levenberg–Marquardt; the model-free route solves a non-negative, curvature-penalized least-squares problem on the grid.'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const tu = v.tu;
      let cin = RTDI.rows(v.cin, 'input');
      const mom = rows => { const R = RTDI.pulse(rows, { rule: 'pl' }); return R; };
      if (v.mode === 'out') {
        /* E as a function */
        let Ef, tE = 0, Elist = null;
        if (v.Esrc === 'table') { Elist = RTDI.rows(v.etab, 'E'); const pl = RT.pl(Elist.map(r => r[0]), Elist.map(r => r[1])); Ef = pl.E; tE = Elist[Elist.length - 1][0]; }
        else {
          const tau = +v.mtau;
          if (v.model === 'mixed') Ef = t => RT.Emixed(t, tau);
          else if (v.model === 'tanks') Ef = t => RT.Etanks(t, tau, +v.mN);
          else if (v.model === 'disp') { const g = RT.sampleSmooth(t => RT.Eclosed(t, tau, +v.md), tau * 8, 2000); Ef = t => RT.at(g, t); }
          else if (v.model === 'plugmix') { const tp = +v.mfp * tau, tm = tau - tp; Ef = t => (t >= tp ? RT.Emixed(t - tp, tm) : 0); }
          else Ef = t => RT.Elam(t, tau);
          tE = tau * 8;
        }
        let tOut = [], cOut = [];
        if (v.disc === 'disc') {
          /* discrete slices: both signals at the input's spacing */
          const dt = cin.length > 1 ? cin[1][0] - cin[0][0] : 1;
          const Ein = Elist ? Elist : N.range(Math.ceil(tE / dt) + 1).map(i => [i * dt, Ef(i * dt)]);
          const tl = cin[cin.length - 1][0] + Ein[Ein.length - 1][0];
          for (let t = cin[0][0]; t <= tl + 1e-9; t += dt) {
            let s = 0;
            cin.forEach(([ti, ci]) => { const tp = t - ti; const e = Ein.find(r => Math.abs(r[0] - tp) < 1e-9 * Math.max(1, dt)); if (e) s += ci * e[1] * dt; });
            tOut.push(+t.toFixed(10)); cOut.push(s);
          }
        } else {
          const tmax = cin[cin.length - 1][0] + tE * 1.05, n = 2400;
          const plIn = RT.pl(cin.map(r => r[0]), cin.map(r => r[1]));
          const A = RT.sample(plIn.E, tmax, n), B = RT.sample(Ef, tmax, n);
          const Cc = RT.conv(A, B);
          for (let i = 0; i < n; i += 8) { tOut.push(i * Cc.dt); cOut.push(Cc.c[i]); }
        }
        const areaIn = v.disc === 'disc' ? N.sum(cin.map(r => r[1])) * (cin.length > 1 ? cin[1][0] - cin[0][0] : 1) : RT.pl(cin.map(r => r[0]), cin.map(r => r[1])).area;
        const areaOut = v.disc === 'disc' ? N.sum(cOut) * (tOut.length > 1 ? tOut[1] - tOut[0] : 1) : N.trapz(tOut, cOut);
        kpi.push({ k: 'areaIn', l: 'Area under C_in', v: areaIn }, { k: 'areaOut', l: 'Area under C_out', v: areaOut, s: Math.abs(areaOut - areaIn) < 0.01 * areaIn ? 'equal, as it must be' : 'differs: E may not have unit area' });
        [9, 10, 11].forEach(tq => { const i = tOut.findIndex(t => Math.abs(t - tq) < 1e-9); if (i >= 0) kpi.push({ k: 'Cout' + tq, l: 'C_out at t = ' + tq + ' ' + tu, v: cOut[i] }); });
        const mi = mom(cin), mo = mom(tOut.map((t, i) => [t, cOut[i]]));
        kpi.push({ k: 'dt', l: 'Mean shift t̄_out − t̄_in', v: mo.mean - mi.mean, u: tu, main: true }, { k: 'dvar', l: 'Variance added', v: mo.var - mi.var, u: tu + '²' });
        const plots = [{ title: 'Input, vessel and output (Fig. E11.3b)', x: { label: 't', unit: tu, min: 0 }, y: { label: 'C' , min: 0 }, series: [{ name: 'C_in', x: cin.map(r => r[0]), y: cin.map(r => r[1]), color: 1, type: v.disc === 'disc' ? 'bar' : undefined, marker: v.disc === 'disc' ? undefined : 'o' }, { name: 'C_out', x: tOut, y: cOut, color: 2, type: v.disc === 'disc' ? 'bar' : undefined }] }];
        const tg = N.linspace(0, tE, 200);
        plots.push({ title: 'E of the vessel', x: { label: 't', unit: tu, min: 0 }, y: { label: 'E', unit: '1/' + tu, min: 0 }, series: [{ name: 'E', x: Elist ? Elist.map(r => r[0]) : tg, y: Elist ? Elist.map(r => r[1]) : tg.map(Ef), color: 3, type: Elist && v.disc === 'disc' ? 'bar' : undefined }] });
        return { kpi, warn, info, plots, tables: [{ title: 'Output signal', cols: ['t, ' + tu, 'C_out'], rows: tOut.map((t, i) => [t, cOut[i]]).filter((r, i) => v.disc === 'disc' || i % 6 === 0) }] };
      }
      /* ---------- deconvolution ---------- */
      let cout = (v.cout || []).filter(q => q && isFinite(q[0]) && isFinite(q[1]));
      if (v._ex144 && cout.length < 2) {
        /* the book gives only moments; synthesize an output with t̄ = 280 s, σ² = 1000 s² by passing the input through N = 4 tanks */
        const tmax = 520, n = 2600, plIn = RT.pl(cin.map(r => r[0]), cin.map(r => r[1]));
        const A = RT.sample(plIn.E, tmax, n), B = RT.sample(t => RT.Etanks(t, 60, 4), tmax, n), Cc = RT.conv(A, B);
        cout = []; for (let i = 0; i < n; i += 4) cout.push([i * Cc.dt, Cc.c[i]]);
        info.push('The book gives only the moments of the output (t̄ = 280 s, σ² = 1000 s²); an output with exactly those moments has been generated from the input for this example.');
      }
      cout = RTDI.rows(cout, 'output');
      const mi = mom(cin), mo = mom(cout);
      const dT = mo.mean - mi.mean, dV = mo.var - mi.var;
      kpi.push({ k: 'dT', l: 'Vessel mean time Δt̄', v: dT, u: tu, main: true }, { k: 'dV', l: 'Vessel variance Δσ²', v: dV, u: tu + '²' });
      if (!(dT > 0) || !(dV > 0)) throw new Error('The output must come later and be wider than the input (Δt̄ > 0, Δσ² > 0).');
      const tmax = cout[cout.length - 1][0] * 1.05, n = 1600;
      const plIn = RT.pl(cin.map(r => r[0]), cin.map(r => r[1])), A = RT.sample(plIn.E, tmax, n);
      const plOut = RT.pl(cout.map(r => r[0]), cout.map(r => r[1]));
      const s2 = dV / (dT * dT);
      let Ef, label, par = '';
      const predict = Efn => { const B = RT.sample(Efn, tmax, n); return RT.conv(A, B); };
      const obsAt = t => plOut.E(t);
      if (v.fitm === 'free') {
        /* non-negative, curvature-penalized least squares on a coarse grid */
        const m = 120, dt = tmax / (m - 1);
        const cA = N.range(m).map(i => plIn.E(i * dt));
        const Am = N.range(m).map(i => N.range(m).map(j => (j <= i ? cA[i - j] * dt * (j === 0 || j === i ? 0.5 : 1) : 0)));
        const y = N.range(m).map(i => obsAt(i * dt));
        const lam = +v.reg * N.sum(y.map(q => q * q)) / m * 1e-3;
        let e = new Array(m).fill(1 / (dT * 4));
        const loss = ee => { let s = 0; for (let i = 0; i < m; i++) { let p = 0; for (let j = 0; j < m; j++) p += Am[i][j] * ee[j]; s += (p - y[i]) ** 2; } for (let i = 1; i < m - 1; i++) s += lam * (ee[i - 1] - 2 * ee[i] + ee[i + 1]) ** 2; return s; };
        /* projected gradient with Barzilai–Borwein steps */
        const grad = ee => { const g = new Array(m).fill(0), r = N.range(m).map(i => { let p = 0; for (let j = 0; j < m; j++) p += Am[i][j] * ee[j]; return p - y[i]; }); for (let j = 0; j < m; j++) { let s = 0; for (let i = 0; i < m; i++) s += Am[i][j] * r[i]; g[j] = 2 * s; } for (let i = 1; i < m - 1; i++) { const c = 2 * lam * (ee[i - 1] - 2 * ee[i] + ee[i + 1]); g[i - 1] += c; g[i] -= 2 * c; g[i + 1] += c; } return g; };
        let g = grad(e), step = 1e-3;
        for (let it = 0; it < 1500; it++) {
          const en = e.map((q, i) => Math.max(0, q - step * g[i]));
          const gn = grad(en);
          const sv = en.map((q, i) => q - e[i]), yv = gn.map((q, i) => q - g[i]);
          const sy = N.sum(sv.map((q, i) => q * yv[i])), ss = N.sum(sv.map(q => q * q));
          step = sy > 0 ? ss / sy : step * 0.5;
          e = en; g = gn;
          if (ss < 1e-24) break;
        }
        const area = N.trapz(N.range(m).map(i => i * dt), e);
        const xs = N.range(m).map(i => i * dt);
        Ef = t => N.interp(xs, e, t) / Math.max(area, 1e-300);
        label = 'regularized deconvolution';
        kpi.push({ k: 'Earea', l: 'Area of the recovered E', v: area, s: 'should be near 1' });
        info.push('Model-free deconvolution is ill-posed: raise the smoothing weight if E looks noisy, lower it if E looks too smooth. Final loss ' + N.fmt(loss(e), 3) + '.');
      } else if (v.how === 'mom') {
        if (v.fitm === 'tanks') { const Nq = 1 / s2; Ef = t => RT.Etanks(t, dT, Nq); label = 'N = ' + N.fmt(Nq, 4) + ' tanks'; kpi.push({ k: 'N', l: 'Tanks in series, N = Δt̄²/Δσ²', v: Nq, main: true }); par = 'E = ' + N.fmt(Math.pow(Nq / dT, Nq) / Math.exp(N.lgamma(Nq)), 5) + ' t^{' + N.fmt(Nq - 1, 3) + '} e^{−' + N.fmt(Nq / dT, 4) + 't}'; }
        else { const d = RTDI.dFromVar(s2, 'closed'); if (!isFinite(d)) throw new Error('Δσ²/Δt̄² ≥ 1: no dispersion number fits.'); const g = RT.sampleSmooth(t => RT.Eclosed(t, dT, d), dT * 8, 2000); Ef = t => RT.at(g, t); label = 'D/uL = ' + N.fmt(d, 4); kpi.push({ k: 'd', l: 'D/uL (closed vessel)', v: d, main: true }); }
      } else {
        /* least squares on the output curve */
        const ts = cout.map(r => r[0]), ys = cout.map(r => r[1]);
        const scaleOut = N.sum(ys.map(q => q * q));
        const mk = p => (v.fitm === 'tanks' ? (t => RT.Etanks(t, p[0], p[1])) : (() => { const g = RT.sampleSmooth(t => RT.Eclosed(t, p[0], p[1]), p[0] * 8, 1200); return t => RT.at(g, t); })());
        const resid = p => { const Cc = predict(mk(p)); return ts.map((t, i) => (RT.at(Cc, t) - ys[i]) / Math.sqrt(scaleOut)); };
        const p0 = [dT, v.fitm === 'tanks' ? 1 / s2 : Math.max(1e-4, RTDI.dFromVar(s2, 'closed') || 0.05)];
        const fit = N.lm(resid, p0, { pos: [true, true], maxIt: 60 });
        const p = fit.p || fit.x || p0;
        Ef = mk(p);
        label = v.fitm === 'tanks' ? 'N = ' + N.fmt(p[1], 4) + ' tanks (least squares)' : 'D/uL = ' + N.fmt(p[1], 4) + ' (least squares)';
        kpi.push({ k: v.fitm === 'tanks' ? 'N' : 'd', l: v.fitm === 'tanks' ? 'Tanks in series (least squares)' : 'D/uL (least squares)', v: p[1], main: true }, { k: 'tfit', l: 'Vessel mean time (least squares)', v: p[0], u: tu });
      }
      if (par) info.push('Recovered E: ' + par + ' (t in ' + tu + ').');
      const Cp = predict(Ef);
      const tg = N.linspace(0, tmax, 300);
      const plots = [{ title: 'Signals: measured output against input ∗ E', x: { label: 't', unit: tu, min: 0 }, y: { label: 'C', min: 0 }, series: [{ name: 'C_in', x: cin.map(r => r[0]), y: cin.map(r => r[1]), color: 1 }, { name: 'C_out measured', x: cout.map(r => r[0]), y: cout.map(r => r[1]), color: 2, type: 'pts' }, { name: 'C_in ∗ E', x: tg, y: tg.map(t => RT.at(Cp, t)), color: 3 }] }];
      const te = N.linspace(0, dT * 3.5, 300);
      plots.push({ title: 'Recovered E of the vessel (Fig. E14.4b)', x: { label: 't', unit: tu, min: 0 }, y: { label: 'E', unit: '1/' + tu, min: 0 }, series: [{ name: label, x: te, y: te.map(Ef), color: 3 }, { x: te, y: te.map(Ef), type: 'area', fill: 0.15, color: 3, hover: false, legend: false }] });
      return { kpi, warn, info, plots, tables: [{ title: 'Moments', cols: ['signal', 't̄', 'σ²'], rows: [['input', mi.mean, mi.var], ['output', mo.mean, mo.var], ['vessel (difference)', dT, dV]] }] };
    }
  });

  /* ===================================================================== */
  /* F-1103  conversion from the RTD: segregation and maximum mixedness     */
  /* ===================================================================== */
  ST.tool({
    id: 'rtdconv', tag: 'F-1103', ch: 11, sec: 'Ch. 11.2, 16', title: 'Conversion from the RTD: macrofluid, latest and earliest mixing',
    lede: 'Any single reaction through any RTD — measured, or from a model. The segregated (macrofluid, or latest-mixing microfluid) conversion and Zwietering’s maximum-mixedness bound for the earliest-mixing microfluid, beside plug flow, mixed flow, and the dispersion and tanks-in-series models of the same variance.',
    keys: 'conversion non ideal reactor segregation model macrofluid microfluid maximum mixedness zwietering rtd e curve example 11.4 11.5 13.4 bounds earliest latest mixing',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'RTD' },
      { k: 'src', t: 'select', label: 'E from', opts: [['data', 'pulse data (table)'], ['peaks', 'ideal peaks (table)'], ['pl', 'a drawn E curve (t, E corners)'], ['mixed', 'mixed flow'], ['tanks', 'N tanks in series'], ['disp', 'dispersion, closed'], ['lam', 'laminar flow'], ['pm', 'plug + mixed in series']], v: 'data' },
      { k: 'data', t: 'table', label: 'Readings (t, C) or peaks (t, area)', cols: ['t', 'C'], v: EX111, if: v => v.src === 'data' || v.src === 'peaks' },
      { k: 'rule', t: 'seg', label: 'Integrate', opts: [['pl', 'piecewise linear'], ['sum', 'Σ f E Δt (book)']], v: 'sum', if: v => v.src === 'data' },
      { k: 'pl', t: 'table', label: 'E curve corners (t, E)', cols: ['t', 'E'], v: [[1, 0], [1, 0.5], [3, 0.5], [3, 0]], if: v => v.src === 'pl', hint: 'Repeat a time to draw a vertical edge.' },
      { k: 'tau', label: 'Mean residence time t̄', unit: '[t]', v: 15, min: 1e-9, max: 1e9, scale: 'log', if: v => ['mixed', 'tanks', 'disp', 'lam', 'pm'].indexOf(v.src) >= 0 },
      { k: 'Nt', label: 'N', v: 4, min: 0.5, max: 500, scale: 'log', if: v => v.src === 'tanks' },
      { k: 'd', label: 'D/uL', v: 0.12, min: 1e-5, max: 50, scale: 'log', if: v => v.src === 'disp' },
      { k: 'fp', label: 'Plug-flow share of t̄', v: 0.5, min: 0, max: 0.999, if: v => v.src === 'pm' }
    ].concat(KIT.kin({ gas: false, d: { law: 'power', k: 0.307, a: 1, CA0: 1 } })),
    presets: [
      { name: 'Example 11.4', v: { src: 'data', data: EX111, rule: 'sum', law: 'power', k: 0.307, a: 1, ob: 0, CA0: 1, _cu: 'mol/L', _tu: 'min' }, book: { unSeg: { v: 0.0469, tol: 0.02 }, unPlug: { v: 0.01, tol: 0.03 }, unDisp: { v: 0.035, tol: 0.04 } }, note: 'The vessel of Example 11.1 with a first-order decomposition, k = 0.307 min⁻¹: 4.7% unconverted against 1.0% in plug flow; the dispersion model of equal variance (Example 13.4) leaves about 3.5%.' },
      { name: 'Example 11.5', v: { src: 'pl', pl: [[1, 0], [1, 0.5], [3, 0.5], [3, 0]], law: 'power', k: 0.5, a: 2, ob: 0, CA0: 2, _cu: 'mol/L', _tu: 'min' }, book: { CSeg: { v: 0.347, tol: 0.01 } }, note: 'Noncoalescing droplets (a macrofluid), C_A0 = 2 mol/L, −r_A = 0.5 C_A², E = 0.5 between 1 and 3 min: mean C_A/C_A0 = 0.5 ln 2 = 0.347.' },
      { name: 'Example 16.1 (plug + mixed)', v: { src: 'pm', tau: 2, fp: 0.5, law: 'power', k: 1, a: 2, ob: 0, CA0: 1, _cu: 'mol/L', _tu: 'min' }, book: { CSeg: { v: 0.362, tol: 0.01 } }, note: 'Equal plug and mixed reactors, τ = 1 min each, second order with k = C_A0 = 1. Late mixing or a macrofluid gives C = 0.362; Ex. 16.1 also finds 0.382 (mixed first) and 0.366 (plug first) — see F-1601.' },
      { name: 'Problem 11.7', v: { src: 'pl', pl: [[0, 0], [0, 0.5], [2, 0.5], [2, 0]], law: 'power', k: 2, a: 0.5, ob: 0, CA0: 1, _cu: 'mol/L', _tu: 'min' }, note: 'A half-order macrofluid through a box-shaped RTD (sketch read from Fig. P11.7).' }
    ],
    eqs: [
      ['\\frac{\\bar C_A}{C_{A0}} = \\int_0^\\infty \\left(\\frac{C_A}{C_{A0}}\\right)_{batch} E\\,dt', 'Eq. 11.13', 'Macrofluid, or a microfluid mixing as late as possible.'],
      ['\\frac{dC_A}{d\\lambda} = (-r_A) + \\frac{E(\\lambda)}{1-F(\\lambda)}\\,(C_A - C_{A0})', 'Zwietering (1959)', 'Maximum mixedness: earliest mixing; integrated from λ = ∞ to 0.'],
      ['\\left(\\frac{C_A}{C_{A0}}\\right)_{batch} = e^{-kt},\\ \\frac{1}{1+kC_{A0}t},\\ \\left[1+(n-1)C_{A0}^{n-1}kt\\right]^{1/(1-n)}', 'Eqs. 11.14–11.16'],
      ['\\frac{C}{C_0} = \\frac{4a\\,e^{Pe/2}}{(1+a)^2e^{aPe/2}-(1-a)^2e^{-aPe/2}},\\quad a = \\sqrt{1+4k\\tau\\,D/uL}', 'Eq. 13.19', 'Dispersion model of equal variance, first order.']
    ],
    method: ['The batch history X(t) is integrated once with dense output, then weighted by E: exactly on each linear piece of a drawn E, or by the book’s sums for tabulated data, plus discrete atoms for ideal peaks.', 'Maximum mixedness integrates Zwietering’s equation in the life-expectancy variable from the end of the curve, where it starts from the algebraic balance with the local intensity Λ = E/(1−F); for plug-flow parts (atoms) the bound is not defined and only the segregated value is given.', 'Comparison models share the mean and variance of the RTD: N = t̄²/σ² mixed tanks solved in turn, and the closed dispersion model solved as a boundary-value problem.'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const m = KIT.model(v);
      const tu = v._tu || 'min';
      if (m.eps) warn.push('Segregation analysis assumes constant density.');
      /* build E */
      let E, atoms = [], tb, s2, tmaxE, plE = null, sumRule = false, sumPts = null;
      if (v.src === 'data') {
        const R = RTDI.pulse(RTDI.rows(v.data, 'readings'), { rule: v.rule });
        tb = R.mean; s2 = R.var; plE = R.pl; tmaxE = R.t[R.t.length - 1];
        if (v.rule === 'sum') { sumRule = true; sumPts = R.t.map((t, i) => { const n = R.t.length, dt = i === 0 ? (R.t[1] - R.t[0]) / 2 : i === n - 1 ? (R.t[n - 1] - R.t[n - 2]) / 2 : (R.t[i + 1] - R.t[i - 1]) / 2; return [t, R.E[i], dt]; }); }
        E = plE.E;
      } else if (v.src === 'peaks') {
        const rows = RTDI.rows(v.data, 'peaks'), area = N.sum(rows.map(r => r[1]));
        atoms = rows.map(r => [r[0], r[1] / area]);
        tb = N.sum(atoms.map(a => a[0] * a[1])); s2 = N.sum(atoms.map(a => a[1] * (a[0] - tb) ** 2)); tmaxE = atoms[atoms.length - 1][0];
        E = () => 0;
      } else if (v.src === 'pl') {
        const rows = (v.pl || []).filter(q => q && isFinite(q[0]) && isFinite(q[1]));
        if (rows.length < 2) throw new Error('Give the corners of the E curve.');
        plE = RT.pl(rows.map(r => r[0]), rows.map(r => r[1]));
        if (Math.abs(plE.area - 1) > 0.02) warn.push('The drawn E curve has area ' + N.fmt(plE.area, 4) + ', not 1: it has been normalized.');
        const A = plE.area;
        const pl2 = RT.pl(rows.map(r => r[0]), rows.map(r => r[1] / A));
        plE = pl2; E = pl2.E; tmaxE = rows[rows.length - 1][0];
        tb = pl2.integ(t => t); s2 = pl2.integ(t => (t - tb) * (t - tb));
      } else {
        const tau = +v.tau;
        tb = tau; tmaxE = tau * 12;
        if (v.src === 'mixed') { E = t => RT.Emixed(t, tau); s2 = tau * tau; }
        else if (v.src === 'tanks') { E = t => RT.Etanks(t, tau, +v.Nt); s2 = tau * tau / +v.Nt; tmaxE = tau * (1 + 12 / Math.sqrt(+v.Nt)); }
        else if (v.src === 'disp') { const g = RT.sampleSmooth(t => RT.Eclosed(t, tau, +v.d), tau * (4 + 20 * +v.d), 3000); E = t => RT.at(g, t); s2 = tau * tau * RT.varClosed(+v.d); tmaxE = g.tmax; }
        else if (v.src === 'lam') { E = t => RT.Elam(t, tau); s2 = NaN; tmaxE = tau * 60; }
        else { const tp = +v.fp * tau, tm = tau - tp; E = t => (t >= tp ? RT.Emixed(t - tp, tm) : 0); s2 = tm * tm; tmaxE = tp + 14 * tm; }
      }
      const Xb = RT.batchX(m, tmaxE * 1.01 + 1e-9);
      const cOf = X => 1 - X;
      /* segregated */
      let unSeg;
      if (atoms.length) unSeg = N.sum(atoms.map(([t, w]) => w * cOf(Xb(t))));
      else if (sumRule) unSeg = N.sum(sumPts.map(([t, e, dt]) => cOf(Xb(t)) * e * dt));
      else if (plE) unSeg = plE.integ(t => cOf(Xb(t)), 1e-10);
      else if (v.src === 'lam') { const tau = +v.tau; unSeg = N.quad(t => cOf(Xb(t)) * RT.Elam(t, tau), tau / 2, Infinity, { rtol: 1e-9 }); }
      else if (v.src === 'pm') { const tau = +v.tau, tp = +v.fp * tau, tm = tau - tp; unSeg = N.quad(u => cOf(Xb(tp + u)) * RT.Emixed(u, tm), 0, Infinity, { rtol: 1e-10 }); }
      else unSeg = N.quad(t => cOf(Xb(t)) * E(t), 0, tmaxE, { rtol: 1e-9 });
      let unSegC = unSeg;
      if (sumRule) {
        unSegC = plE.integ(t => cOf(Xb(t)), 1e-10);
        info.push('With the book’s sums the segregated result is ' + N.fmt(unSeg, 4) + '; reading the same data as a continuous piecewise-linear curve gives ' + N.fmt(unSegC, 4) + '. The gap comes from the coarse spacing of the readings, not from the physics.');
      }
      kpi.push({ k: 'unSeg', l: 'Unconverted, macrofluid (segregated)', v: unSeg, main: true, s: 'X_A = ' + N.fmt(1 - unSeg, 4) });
      kpi.push({ k: 'CSeg', l: 'Mean C_A/C_A0 leaving, segregated', v: unSeg });
      /* maximum mixedness (smooth E only) */
      let unMM = NaN;
      if (!atoms.length && v.src !== 'lam') {
        const ts = N.linspace(0, tmaxE, 1500), es = ts.map(t => E(t)), Fs = N.cumtrapz(ts, es);
        const Ftot = Fs[Fs.length - 1];
        let li = ts.length - 1;
        while (li > 10 && 1 - Fs[li] / Ftot < 1e-6) li--;
        const Lam = l => { const e = N.interp(ts, es, l) / Ftot, F = N.interp(ts, Fs, l) / Ftot; return e / Math.max(1 - F, 1e-9); };
        try { unMM = 1 - RT.maxMix(m, Lam, ts[li]); } catch (e) { unMM = NaN; }
        if (isFinite(unMM)) kpi.push({ k: 'unMM', l: 'Unconverted, maximum mixedness (earliest-mixing microfluid)', v: unMM, s: 'X_A = ' + N.fmt(1 - unMM, 4) });
      } else info.push('Maximum mixedness needs a smooth E; with ideal peaks or laminar flow only the segregated value is shown.');
      if (isFinite(unMM)) {
        const n = v.law === 'power' ? +v.a : NaN;
        info.push(Math.abs(unSegC - unMM) < 1e-4 * Math.max(unSegC, 1e-9) + 1e-9 ? 'Segregated and maximum-mixedness results agree: for first-order kinetics mixing does not matter, only the RTD.' : 'The real reactor lies between ' + N.fmt(Math.min(unSegC, unMM), 4) + ' and ' + N.fmt(Math.max(unSegC, unMM), 4) + ' unconverted. ' + (isFinite(n) ? (n > 1 ? 'Order above one: segregation and late mixing help.' : n < 1 ? 'Order below one: early mixing helps.' : '') : ''));
      }
      /* references */
      const unPlug = cOf(Xb(tb));
      const tauC = tb;
      const unMix = 1 - K.XMFR(m, tauC);
      kpi.push({ k: 'unPlug', l: 'Unconverted, plug flow of the same t̄', v: unPlug }, { k: 'unMix', l: 'Unconverted, mixed flow of the same t̄', v: unMix });
      const rows = [['macrofluid / latest mixing' + (sumRule ? ' (book’s sums)' : ''), unSeg, 1 - unSeg]].concat(sumRule ? [['macrofluid, continuous curve', unSegC, 1 - unSegC]] : [], [['maximum mixedness', unMM, 1 - unMM]]).concat([ ['plug flow, same t̄', unPlug, 1 - unPlug], ['mixed flow, same t̄', unMix, 1 - unMix]]);
      if (isFinite(s2) && s2 > 0) {
        const s2t = s2 / (tb * tb), Nq = 1 / s2t;
        const d = RTDI.dFromVar(s2t, 'closed');
        if (isFinite(d)) { const bv = RT.dispersionX(m, tb, d); kpi.push({ k: 'unDisp', l: 'Unconverted, dispersion model of equal variance', v: 1 - bv.X, s: 'D/uL = ' + N.fmt(d, 4) }); rows.push(['dispersion model, D/uL = ' + N.fmt(d, 3), 1 - bv.X, bv.X]); }
        if (Nq >= 1 && Nq <= 200) {
          const nInt = Math.round(Nq);
          let X = 0; for (let i = 0; i < nInt; i++) X = K.XMFR(m, tb / nInt, X);
          kpi.push({ k: 'unTanks', l: 'Unconverted, ' + nInt + ' tanks in series (microfluid)', v: 1 - X, s: 'N = t̄²/σ² = ' + N.fmt(Nq, 3) });
          rows.push([nInt + ' tanks in series', 1 - X, X]);
        }
        kpi.push({ k: 'tbar', l: 't̄ and σ²_θ of the RTD', v: tb, u: tu, s: 'σ²_θ = ' + N.fmt(s2t, 4) });
      } else kpi.push({ k: 'tbar', l: 'Mean residence time', v: tb, u: tu });
      /* plots */
      const ts = N.linspace(0, Math.min(tmaxE, tb * 6), 300);
      const plots = [];
      const ser = [];
      if (atoms.length) ser.push({ name: 'E (peaks, areas)', x: atoms.map(a => a[0]), y: atoms.map(a => a[1]), type: 'bar', color: 1 });
      else ser.push({ name: 'E', x: ts, y: ts.map(E), color: 1 }, { x: ts, y: ts.map(E), type: 'area', fill: 0.12, color: 1, hover: false, legend: false });
      plots.push({ title: 'RTD and the batch history it weights', x: { label: 't', unit: tu, min: 0 }, y: { label: 'E', unit: '1/' + tu, min: 0 }, series: ser });
      plots.push({ title: 'Batch concentration of an element against its age', x: { label: 't', unit: tu, min: 0 }, y: { label: 'C_A/C_A0', min: 0, max: 1 }, series: [{ name: '(C_A/C_A0)_batch', x: ts, y: ts.map(t => cOf(Xb(t))), color: 2 }, { name: 'E·(C_A/C_A0)_batch (scaled)', x: ts, y: ts.map(t => { const e = atoms.length ? 0 : E(t); return e * cOf(Xb(t)) * tb; }), color: 3, dash: '4 3' }], cap: 'The unconverted fraction is the area under E·(C_A/C_A0)_batch: early leavers carry most of the unreacted A (Ex. 11.4).' });
      const tables = [{ title: 'Unconverted fraction by flow assumption', cols: ['model', 'C_A/C_A0', 'X_A'], rows }];
      if (sumRule) tables.push({ title: 'Book’s sum (Table E11.4)', cols: ['t, ' + tu, 'E', 'C/C₀ batch', 'C/C₀·E·Δt'], rows: sumPts.map(([t, e, dt]) => [t, e, cOf(Xb(t)), cOf(Xb(t)) * e * dt]) });
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
