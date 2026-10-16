/* Chapter 13 — The dispersion model: curves and parameter estimation, reaction with dispersion, correlations */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, RT = ST.RT, KIT = ST.KIT, K = ST.K, D = ST.D, RTDI = ST.RTDI;

  /* ===================================================================== */
  /* F-1301  dispersion model: E and F curves, D/uL from measurements       */
  /* ===================================================================== */
  ST.tool({
    id: 'dispersion', tag: 'F-1301', ch: 13, sec: 'Ch. 13.1', title: 'The dispersion model: E and F curves and D/uL from any measurement',
    lede: 'Closed, open and small-deviation tracer curves for any vessel dispersion number, with their exact means and variances; and every way the book extracts D/uL — from a variance, from the change in variance between two probes, from the 16% and 84% points of a step response, or by fitting a whole pulse curve. Plus the √L and L spreading laws that tell dispersion from convection and locate a source upstream.',
    keys: 'dispersion model vessel dispersion number D/uL peclet axial dispersion closed open boundary condition variance gaussian probability paper step response 16 84 percentile spreading river pollution source winery example 13.1 13.2 13.3 14.1 14.2',
    inputs: [
      { group: 'Task' },
      { k: 'what', t: 'select', label: 'Find', opts: [['curves', 'E and F curves for given D/uL'], ['var', 'D/uL from a variance'], ['two', 'D/uL from two probes (Δσ²)'], ['step', 'D/uL from the 16% and 84% points of a step'], ['data', 'D/uL by fitting a pulse curve'], ['spread', 'spreading with distance: predict, or find the source']], v: 'curves' },
      { k: 'dlist', t: 'text', label: 'D/uL values to draw', v: '0.002, 0.025, 0.2, 1', mono: true, if: v => v.what === 'curves' },
      { k: 'bc', t: 'seg', label: 'Vessel boundaries', opts: [['closed', 'closed'], ['open', 'open–open'], ['small', 'small deviation (gaussian)']], v: 'closed', if: v => v.what === 'curves' || v.what === 'var' || v.what === 'data' },
      { k: 's2', label: 'Variance σ²', unit: '[t]²', v: 47.5, min: 0, max: 1e12, if: v => v.what === 'var' },
      { k: 'tb', label: 'Mean time t̄', unit: '[t]', v: 15, min: 1e-9, max: 1e12, if: v => v.what === 'var' },
      { k: 'sa', label: 'Variance at the first probe', unit: '[t]²', v: 39, min: 0, max: 1e12, if: v => v.what === 'two' },
      { k: 'sb', label: 'Variance at the second probe', unit: '[t]²', v: 64, min: 0, max: 1e12, if: v => v.what === 'two' },
      { k: 'tsrc', t: 'seg', label: 'Time between probes from', opts: [['t', 'Δt̄ directly'], ['geo', 'distance, voidage and superficial velocity']], v: 'geo', if: v => v.what === 'two' || v.what === 'step' },
      { k: 'dtb', label: 'Δt̄ between probes (or t̄ of the vessel)', unit: '[t]', v: 30, min: 1e-9, max: 1e12, if: v => (v.what === 'two' || v.what === 'step') && v.tsrc === 't' },
      { k: 'Lg', label: 'Distance (or vessel length)', unit: 'length', v: 90, min: 1e-9, max: 1e12, if: v => (v.what === 'two' || v.what === 'step') && v.tsrc === 'geo' },
      { k: 'eps', label: 'Voidage (1 for an empty pipe)', v: 0.4, min: 0.01, max: 1, if: v => (v.what === 'two' || v.what === 'step') && v.tsrc === 'geo' },
      { k: 'us', label: 'Superficial velocity', unit: 'length/[t]', v: 1.2, min: 1e-12, max: 1e9, scale: 'log', if: v => (v.what === 'two' || v.what === 'step') && v.tsrc === 'geo' },
      { k: 't16', label: 'Time at 16% of the step', unit: '[t]', v: 178550, min: 0, max: 1e15, if: v => v.what === 'step' },
      { k: 't84', label: 'Time at 84% of the step', unit: '[t]', v: 187750, min: 0, max: 1e15, if: v => v.what === 'step' },
      { k: 'data', t: 'table', label: 'Pulse response (t, C)', cols: ['t', 'C'], v: [[0, 0], [5, 3], [10, 5], [15, 5], [20, 4], [25, 2], [30, 1], [35, 0]], if: v => v.what === 'data' },
      { k: 'spMode', t: 'seg', label: 'Spreading', opts: [['predict', 'predict at a new distance'], ['source', 'locate the source']], v: 'source', if: v => v.what === 'spread' },
      { k: 'w1', label: 'Spread (width or σ) at station 1', v: 10.5, min: 1e-12, max: 1e12, if: v => v.what === 'spread' },
      { k: 'L1', label: 'Distance of station 1 from the source', v: 32, min: 1e-12, max: 1e12, if: v => v.what === 'spread' && v.spMode === 'predict' },
      { k: 'L2', label: 'Distance at which to predict', v: 50, min: 1e-12, max: 1e12, if: v => v.what === 'spread' && v.spMode === 'predict' },
      { k: 'w2', label: 'Spread at station 2 (downstream)', v: 14, min: 1e-12, max: 1e12, if: v => v.what === 'spread' && v.spMode === 'source' },
      { k: 'dL', label: 'Distance between the stations', v: 119, min: 1e-12, max: 1e12, if: v => v.what === 'spread' && v.spMode === 'source' },
      { k: 'tu', t: 'select', label: 'Time unit', opts: [['s', 's'], ['min', 'min'], ['h', 'h'], ['day', 'day']], v: 'min' }
    ],
    presets: [
      { name: 'Example 13.1', v: { what: 'var', bc: 'closed', s2: 47.5, tb: 15, tu: 'min' }, book: { d: { v: 0.12, tol: 0.01 } }, note: 'The vessel of Example 11.1: σ² = 47.5 min², t̄ = 15 min, closed vessel. First guess σ²_θ/2 = 0.106; corrected, D/uL = 0.12.' },
      { name: 'Example 13.2', v: { what: 'step', tsrc: 'geo', Lg: 1219, eps: 1, us: 0.0067, t16: 178550, t84: 187750, tu: 's' }, book: { d: { v: 0.00032, tol: 0.02 } }, note: 'von Rosenberg’s packed column, 1219 mm long, u = 0.0067 mm/s: the 16% and 84% points of the F curve are 4600 s either side of the mean. D/uL = 0.00032.' },
      { name: 'Example 13.3', v: { what: 'two', tsrc: 'geo', sa: 39, sb: 64, Lg: 90, eps: 0.4, us: 1.2, tu: 's' }, book: { d: { v: 0.0139, tol: 0.01 } }, note: 'Two Geiger counters 90 cm apart in a bed of voidage 0.4, superficial velocity 1.2 cm/s; signal variances 39 and 64 s². D/uL = 1/72.' },
      { name: 'Example 14.1 (winery)', v: { what: 'spread', spMode: 'predict', w1: 8, L1: 32, L2: 50 }, book: { wPred: { v: 10, tol: 0.01 } }, note: 'Eight bottles of rosé from a 32 m line; the line becomes 50 m. Spread grows as √L: 10 bottles.' },
      { name: 'Example 14.2 (river)', v: { what: 'spread', spMode: 'source', w1: 10.5, w2: 14, dL: 119 }, book: { Lsrc: { v: 272, tol: 0.01 } }, note: 'The phenol slug took 10.5 h to pass Portsmouth and 14 h to pass Cincinnati, 119 miles further: the source is about 272 miles above Cincinnati.' },
      { name: 'Example 13.4 vessel by curve fit', v: { what: 'data', bc: 'closed', data: [[0, 0], [5, 3], [10, 5], [15, 5], [20, 4], [25, 2], [30, 1], [35, 0]], tu: 'min' }, note: 'The data of Example 11.1 fitted by the closed dispersion model, by moments and by least squares on the whole curve.' },
      { name: 'E curves, closed vessel (Fig. 13.8)', v: { what: 'curves', bc: 'closed', dlist: '0.002, 0.025, 0.2, 1' } },
      { name: 'E curves, open vessel (Fig. 13.10)', v: { what: 'curves', bc: 'open', dlist: '0.002, 0.025, 0.2, 1' } }
    ],
    eqs: [
      ['E_\\theta = \\frac{1}{\\sqrt{4\\pi(D/uL)}}\\exp\\left[-\\frac{(1-\\theta)^2}{4(D/uL)}\\right],\\quad \\sigma^2_\\theta = 2\\frac{D}{uL}', 'Eqs. 13.7–13.8', 'Small deviation from plug flow.'],
      ['\\sigma^2_\\theta = 2\\frac{D}{uL} - 2\\left(\\frac{D}{uL}\\right)^2\\left(1-e^{-uL/D}\\right),\\quad \\bar\\theta = 1', 'Eq. 13.13', 'Closed vessel.'],
      ['E_{\\theta,oo} = \\frac{1}{\\sqrt{4\\pi\\theta(D/uL)}}\\exp\\left[-\\frac{(1-\\theta)^2}{4\\theta(D/uL)}\\right],\\ \\bar\\theta = 1+2\\frac{D}{uL},\\ \\sigma^2_\\theta = 2\\frac{D}{uL}+8\\left(\\frac{D}{uL}\\right)^2', 'Eqs. 13.14–13.15', 'Open vessel.'],
      ['\\frac{\\Delta\\sigma^2}{(\\Delta\\bar t)^2} = 2\\frac{D}{uL}', 'Eq. 13.11', 'Any one-shot input, two measuring points (Aris).'],
      ['\\sigma = \\frac{t_{84\\%}-t_{16\\%}}{2}', 'Fig. 13.12'],
      ['\\sigma \\propto \\sqrt{L}\\ \\text{(dispersion, tanks)},\\qquad \\sigma \\propto L\\ \\text{(convection)}', 'Eqs. 15.1–15.2']
    ],
    method: ['The closed-vessel E curve has no closed form; it is computed by numerical inversion (fixed Talbot contour) of its Laplace transform, written so that it does not overflow, and checked against the exact mean and variance. For D/uL below 0.008 the gaussian with the exact closed-vessel variance is used.', 'Variance relations are solved exactly for D/uL. Curve fitting adjusts t̄ and D/uL by Levenberg–Marquardt on the normalized E curve.'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const tu = v.tu, plots = [], tables = [];
      const Eof = (bc, d) => (bc === 'open' ? th => RT.Eopen(th, 1, d) : bc === 'small' ? th => RT.Egauss(th, 1, d) : th => RT.EthClosed(th, d));
      const tbFromGeo = () => (v.tsrc === 't' ? +v.dtb : +v.Lg * +v.eps / +v.us);
      if (v.what === 'curves') {
        const ds = String(v.dlist).split(/[,;\s]+/).map(parseFloat).filter(x => x > 0).slice(0, 8);
        if (!ds.length) throw new Error('List one or more D/uL values.');
        const th = N.linspace(0, 3, 400);
        const ser = [], fser = [], rows = [];
        ds.forEach((d, i) => {
          const f = Eof(v.bc, d), ys = th.map(f), F = N.cumtrapz(th, ys);
          ser.push({ name: 'D/uL = ' + N.fmt(d, 3), x: th, y: ys, color: i + 1 });
          fser.push({ name: 'D/uL = ' + N.fmt(d, 3), x: th, y: F.map(q => Math.min(q, 1)), color: i + 1 });
          const mean = v.bc === 'open' ? 1 + 2 * d : 1, vr = v.bc === 'open' ? RT.varOpen(d) : v.bc === 'small' ? 2 * d : RT.varClosed(d);
          rows.push([d, 1 / d, mean, vr, 1 / vr, d < 0.01 ? 'small deviation' : d > 1 ? 'question the model' : '']);
        });
        plots.push({ title: 'E_θ for the dispersion model (' + (v.bc === 'open' ? 'open vessel, Fig. 13.10' : v.bc === 'small' ? 'small deviation, Fig. 13.4' : 'closed vessel, Fig. 13.8') + ')', x: { label: 'θ = t/t̄', min: 0, max: 3 }, y: { label: 'E_θ', min: 0 }, series: ser, hmax: 480 });
        plots.push({ title: 'F curves (Figs. 13.11, 13.13)', x: { label: 'θ', min: 0, max: 3 }, y: { label: 'F', min: 0, max: 1 }, series: fser });
        tables.push({ title: 'Moments', cols: ['D/uL', 'Pe = uL/D', 'θ̄', 'σ²_θ', 'equivalent N = 1/σ²_θ', ''], rows });
        kpi.push({ k: 'n', l: 'Curves drawn', v: ds.length, s: v.bc + ' boundaries' });
        return { kpi, warn, info, plots, tables };
      }
      let d = NaN, s2t = NaN;
      if (v.what === 'var') {
        const tb = +v.tb; s2t = +v.s2 / (tb * tb);
        d = RTDI.dFromVar(s2t, v.bc);
        kpi.push({ k: 's2t', l: 'σ²_θ = σ²/t̄²', v: s2t }, { k: 'dfirst', l: 'First guess σ²_θ/2', v: s2t / 2 });
        if (!isFinite(d)) throw new Error('σ²_θ ≥ 1 cannot come from a closed dispersion vessel.');
      } else if (v.what === 'two') {
        const dt = tbFromGeo(); s2t = (+v.sb - +v.sa) / (dt * dt);
        if (!(s2t > 0)) throw new Error('The second probe must see the wider signal.');
        d = s2t / 2;
        kpi.push({ k: 'dtb', l: 'Mean time between the probes', v: dt, u: tu }, { k: 's2t', l: 'Δσ²_θ', v: s2t });
        info.push('Between two probes inside a bed or pipe the boundary conditions are effectively open, and Δσ²/Δt̄² = 2 D/uL holds whatever the input shape.');
      } else if (v.what === 'step') {
        const tb = tbFromGeo(), sg = (+v.t84 - +v.t16) / 2;
        s2t = Math.pow(sg / tb, 2); d = s2t / 2;
        kpi.push({ k: 'sigma', l: 'σ = (t₈₄ − t₁₆)/2', v: sg, u: tu }, { k: 'tbs', l: 't̄ = L/u', v: tb, u: tu }, { k: 's2t', l: 'σ²_θ', v: s2t });
        if (d > 0.01) warn.push('D/uL above 0.01: the probability-paper method assumes a gaussian F curve and is no longer accurate.');
      } else if (v.what === 'data') {
        const R = RTDI.pulse(RTDI.rows(v.data, 'readings'), { rule: 'pl' });
        s2t = R.var / (R.mean * R.mean);
        d = RTDI.dFromVar(s2t, v.bc);
        kpi.push({ k: 'tbar', l: 't̄ of the data', v: R.mean, u: tu }, { k: 's2t', l: 'σ²_θ of the data', v: s2t }, { k: 'dmom', l: 'D/uL by moments', v: d });
        const ts = R.t, Es = R.E, sc = Math.max(...Es);
        const model = (p, t) => (v.bc === 'open' ? RT.Eopen(t, p[0], p[1]) : v.bc === 'small' ? RT.Egauss(t, p[0], p[1]) : RT.Eclosed(t, p[0], p[1]));
        const fit = N.lm(p => ts.map((t, i) => (model(p, t) - Es[i]) / sc), [R.mean, isFinite(d) ? Math.max(d, 1e-4) : 0.1], { pos: [true, true], maxIt: 60 });
        kpi.push({ k: 'dls', l: 'D/uL by least squares on the curve', v: fit.p[1], main: true, s: 't̄ = ' + N.fmt(fit.p[0], 4) + ' ' + tu }, { k: 'd', l: 'D/uL (moments)', v: d });
        const tg = N.linspace(0, ts[ts.length - 1] * 1.2, 300);
        plots.push({ title: 'Data against the dispersion model (Fig. E13.4)', x: { label: 't', unit: tu, min: 0 }, y: { label: 'E', unit: '1/' + tu, min: 0 }, series: [{ name: 'data', x: ts, y: Es, type: 'pts', color: 1 }, { name: 'moments: D/uL = ' + N.fmt(d, 3), x: tg, y: tg.map(t => (isFinite(d) ? model([R.mean, d], t) : NaN)), color: 2 }, { name: 'least squares: D/uL = ' + N.fmt(fit.p[1], 3), x: tg, y: tg.map(t => model(fit.p, t)), color: 3, dash: '5 4' }], cap: 'Matching the variance alone can misplace the bulk of the curve when the data have a long tail or early leakage.' });
        if (d > 1) warn.push('D/uL > 1: the dispersion model is questionable here (Ch. 13, comment c).');
        return { kpi, warn, info, plots, tables };
      } else if (v.what === 'spread') {
        const w1 = +v.w1;
        if (v.spMode === 'predict') {
          const L1 = +v.L1, L2 = +v.L2;
          kpi.push({ k: 'wPred', l: 'Spread at the new distance, dispersion (σ ∝ √L)', v: w1 * Math.sqrt(L2 / L1), main: true }, { k: 'wConv', l: 'Spread if convection-dominated (σ ∝ L)', v: w1 * L2 / L1 });
          const Ls = N.linspace(0, Math.max(L1, L2) * 1.3, 120);
          plots.push({ title: 'Spread against distance', x: { label: 'distance from the source', min: 0 }, y: { label: 'spread', min: 0 }, series: [{ name: 'dispersion or tanks: √L', x: Ls, y: Ls.map(L => w1 * Math.sqrt(L / L1)), color: 1 }, { name: 'pure convection: L', x: Ls, y: Ls.map(L => w1 * L / L1), color: 2, dash: '5 4' }], marks: [{ t: 'pt', x: L1, y: w1, color: 8, label: 'known' }] });
        } else {
          const w2 = +v.w2, dL = +v.dL;
          if (!(w2 > w1)) throw new Error('The downstream spread must be the larger one.');
          const Ld = dL * w2 * w2 / (w2 * w2 - w1 * w1), Lc = dL * w2 / (w2 - w1);
          kpi.push({ k: 'Lsrc', l: 'Source upstream of station 2, dispersion (σ ∝ √L)', v: Ld, main: true, s: N.fmt(Ld - dL, 4) + ' above station 1' }, { k: 'Lconv', l: 'If convection-dominated (σ ∝ L)', v: Lc });
          info.push('If the release was not instantaneous, the source can only be closer than this: anything within ' + N.fmt(Ld, 4) + ' upstream is suspect.');
          const Ls = N.linspace(0, Ld * 1.2, 160);
          plots.push({ title: 'Spread against distance from the source', x: { label: 'distance upstream of station 2', min: 0 }, y: { label: 'spread', min: 0 }, series: [{ name: 'σ ∝ √L through both stations', x: Ls, y: Ls.map(L => w2 * Math.sqrt(L / Ld)), color: 1 }], marks: [{ t: 'pt', x: Ld - dL, y: w1, color: 2, label: 'station 1' }, { t: 'pt', x: Ld, y: w2, color: 3, label: 'station 2' }] });
        }
        return { kpi, warn, info, plots, tables };
      }
      kpi.push({ k: 'd', l: 'Vessel dispersion number D/uL', v: d, main: true, s: d < 0.01 ? 'small: gaussian shape, boundaries do not matter' : d > 1 ? 'large: question the model' : 'intermediate' });
      kpi.push({ k: 'Pe', l: 'Pe = uL/D', v: 1 / d }, { k: 'Neq', l: 'Equivalent tanks in series', v: 1 / (v.what === 'var' && v.bc === 'closed' ? RT.varClosed(d) : 2 * d) });
      const th = N.linspace(0, 2.5, 300);
      plots.push({ title: 'E_θ for D/uL = ' + N.fmt(d, 3), x: { label: 'θ', min: 0 }, y: { label: 'E_θ', min: 0 }, series: [{ name: 'closed', x: th, y: th.map(q => RT.EthClosed(q, d)), color: 1 }, { name: 'open–open', x: th, y: th.map(q => RT.Eopen(q, 1, d)), color: 2, dash: '5 4' }, { name: 'gaussian', x: th, y: th.map(q => RT.Egauss(q, 1, d)), color: 3, dash: '2 3' }] });
      return { kpi, warn, info, plots, tables };
    }
  });

  /* ===================================================================== */
  /* F-1302  reaction with axial dispersion                                 */
  /* ===================================================================== */
  ST.tool({
    id: 'disprxn', tag: 'F-1302', ch: 13, sec: 'Ch. 13.3', title: 'Reactors with axial dispersion: any kinetics, profiles and the design charts',
    lede: 'Solves the dispersion-reaction boundary-value problem with Danckwerts conditions for any single reaction — first order matched to Wehner and Wilhelm’s closed form — giving outlet conversion, the inlet jump and the axial profile, the volume needed against plug flow, and the book’s comparison charts (Figs. 13.19 and 13.20) redrawn for your kinetics.',
    keys: 'dispersion reaction wehner wilhelm danckwerts boundary conditions axial dispersion reactor conversion design chart second order plug flow comparison volume ratio example 13.4 problem 13.11 13.12',
    units: { C: 'mol/L', t: 'min' }, slow: true,
    inputs: [
      { group: 'Reactor' },
      { k: 'd', label: 'Vessel dispersion number D/uL', v: 0.12, min: 1e-6, max: 100, scale: 'log' },
      { k: 'find', t: 'seg', label: 'Find', opts: [['X', 'conversion for a space time'], ['tau', 'space time for a conversion']], v: 'X' },
      { k: 'tau', label: 'Space time τ', unit: '[t]', v: 15, min: 1e-9, max: 1e9, scale: 'log', if: v => v.find === 'X' },
      { k: 'Xt', label: 'Conversion wanted', v: 0.99, min: 1e-6, max: 0.999999, if: v => v.find === 'tau' },
      { k: 'chart', t: 'check', label: 'Draw the design chart for this reaction', v: true }
    ].concat(KIT.kin({ gas: false, d: { law: 'power', k: 0.307, a: 1, CA0: 1 } })),
    presets: [
      { name: 'Example 13.4', v: { d: 0.12, find: 'X', tau: 15, law: 'power', k: 0.307, a: 1, ob: 0, CA0: 1, chart: true, _tu: 'min' }, book: { unconv: { v: 0.035, tol: 0.04 } }, note: 'kτ = 4.6 and D/uL = 0.12 (from Example 13.1): about 3.5% unconverted, against 1% in plug flow.' },
      { name: 'Second order (Fig. 13.20)', v: { d: 0.05, find: 'X', tau: 10, law: 'power', k: 1, a: 2, ob: 0, CA0: 1, chart: true }, note: 'kC_A0τ = 10: compare the closed-vessel dispersion result with plug and mixed flow.' },
      { name: 'Problem 13.11 (packed laboratory reactor)', v: { d: 0.0104, find: 'X', tau: 1, law: 'power', k: 4.6, a: 1, ob: 0, CA0: 1, chart: false, _tu: 's' }, note: '48 cm of 5 mm pellets, Re_p = 200: D/ud_p ≈ 1 from Fig. 13.17 gives D/uL ≈ 0.5/48. With 99% conversion in 1 s, compare k from plug flow (4.6 s⁻¹) with k that accounts for dispersion — use “space time for a conversion”.' }
    ],
    eqs: [
      ['\\frac{D}{uL}\\frac{d^2C_A}{dz^2} - \\frac{dC_A}{dz} - k\\tau C_{A0}^{n-1}C_A^n = 0', 'Eq. 13.18'],
      ['C_{A0} = C_A - \\frac{D}{uL}\\frac{dC_A}{dz}\\Big|_{0^+},\\qquad \\frac{dC_A}{dz}\\Big|_{1} = 0', '', 'Danckwerts conditions, closed vessel.'],
      ['\\frac{C_A}{C_{A0}} = \\frac{4a\\exp\\left(\\frac{1}{2}\\frac{uL}{D}\\right)}{(1+a)^2\\exp\\left(\\frac{a}{2}\\frac{uL}{D}\\right)-(1-a)^2\\exp\\left(-\\frac{a}{2}\\frac{uL}{D}\\right)},\\ a = \\sqrt{1+4k\\tau\\frac{D}{uL}}', 'Eq. 13.19', 'First order (Wehner and Wilhelm).'],
      ['\\frac{C_A}{C_{A0}} = \\exp\\left[-k\\tau + (k\\tau)^2\\frac{D}{uL}\\right],\\qquad \\frac{L}{L_p} = 1 + (k\\tau)\\frac{D}{uL}', 'Eqs. 13.20–13.21', 'Small deviation from plug flow.']
    ],
    method: ['Second-order finite differences on up to 20 000 nodes (upwinded where the cell Péclet number exceeds 2), with the inlet condition imposed through a ghost node, and damped Newton iterations on the tridiagonal system started from the plug-flow profile. The first-order result is checked against Eq. 13.19 to better than 0.1% in the tests.', 'Space time for a target conversion is found by bracketing and Brent’s method on log τ.'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const m = KIT.model(v), d = +v.d, tu = v._tu || 'min';
      let tau = +v.tau, sol;
      if (v.find === 'tau') {
        const Xt = +v.Xt;
        if (!(Xt < Math.min(m.Xe, m.Xmax))) throw new Error('That conversion is beyond equilibrium or the stoichiometric limit.');
        const tp = K.tauPFR(m, Xt);
        let gs = null;
        const lt = N.root(lt0 => { const r = RT.dispersionX(m, Math.exp(lt0), d, { M: 300, guess: gs }); gs = r.prof; return r.X - Xt; }, Math.log(tp * 0.5), Math.log(tp * 1e4), { n: 40 });
        if (!isFinite(lt)) throw new Error('No space time reaches that conversion.');
        tau = Math.exp(lt);
        kpi.push({ k: 'tau', l: 'Space time needed', v: tau, u: tu, main: true }, { k: 'tauP', l: 'Plug flow would need', v: tp, u: tu, s: 'ratio ' + N.fmt(tau / tp, 4) });
        if (v.law === 'power' && Math.abs(+v.a - 1) < 1e-12) {
          kpi.push({ k: 'ratioS', l: 'Small-deviation estimate L/L_p = 1 + kτ·D/uL', v: 1 + (+v.k) * tp * d });
        }
      }
      sol = RT.dispersionX(m, tau, d);
      const X = sol.X;
      kpi.push({ k: 'X', l: 'Conversion leaving', v: X, main: v.find === 'X' }, { k: 'unconv', l: 'C_A/C_A0 leaving', v: 1 - X, main: v.find === 'X' });
      const Xp = K.XPFR(m, tau).X, Xm = K.XMFR(m, tau);
      kpi.push({ k: 'Xp', l: 'Plug flow, same τ', v: Xp }, { k: 'Xm', l: 'Mixed flow, same τ', v: Xm }, { k: 'jump', l: 'Conversion just inside the inlet', v: sol.X0, s: 'Danckwerts jump' });
      if (v.law === 'power' && Math.abs(+v.a - 1) < 1e-12) {
        const ww = RT.wehner(+v.k * tau, d);
        kpi.push({ k: 'ww', l: 'Wehner–Wilhelm (Eq. 13.19)', v: ww, s: 'numerical: ' + N.fmt(1 - X, 5) });
        if (d < 0.05) kpi.push({ k: 'small', l: 'Small-deviation estimate (Eq. 13.20)', v: Math.exp(-(+v.k) * tau + Math.pow(+v.k * tau, 2) * d) });
      }
      const Nq = 1 / RT.varClosed(d);
      if (Nq >= 1 && Nq < 300) { const nI = Math.round(Nq); let Xt = 0; for (let i = 0; i < nI; i++) Xt = K.XMFR(m, tau / nI, Xt); kpi.push({ k: 'Xtanks', l: nI + ' tanks in series (equal variance)', v: Xt }); }
      if (!sol.ok) warn.push('The boundary-value solver did not fully converge; treat the last digits with caution.');
      const plots = [{ title: 'Concentration along the reactor', x: { label: 'z = x/L', min: 0, max: 1 }, y: { label: 'C_A/C_A0', min: 0, max: 1 }, series: [{ name: 'dispersion, D/uL = ' + N.fmt(d, 3), x: [0].concat(sol.z), y: [1].concat(sol.x.map(x => 1 - x)), color: 1 }, { name: 'plug flow', x: N.linspace(0, 1, 60), y: N.linspace(0, 1, 60).map(z => 1 - K.XPFR(m, z * tau).X), color: 8, dash: '5 4' }, { name: 'mixed flow (outlet)', x: [0, 1], y: [1 - Xm, 1 - Xm], color: 2, dash: '2 3' }], cap: 'The step at z = 0 is the Danckwerts inlet jump: dispersion carries reaction back into the feed.' }];
      if (v.chart) {
        /* Fig. 13.19/13.20 style: V/Vp against C/C0 for several D/uL. Each dispersion line is traced by
           sweeping τ with warm-started solutions, then read at the target exit concentrations. */
        const first = v.law === 'power' && Math.abs(+v.a - 1) < 1e-12 && !m.useB;
        const dl = [0, 0.0625, 0.25, 1, Infinity].concat(d > 0 && [0.0625, 0.25, 1].indexOf(d) < 0 ? [d] : []);
        const fr = N.logspace(1e-3, 0.95, 30);
        const top = Math.min(m.Xe, m.Xmax);
        const tpOf = cf => (1 - cf < top ? K.tauPFR(m, 1 - cf) : NaN);
        const ser = dl.map((dd, i) => {
          let pts;
          if (dd === 0) pts = fr.map(cf => [cf, 1]);
          else if (dd === Infinity) pts = fr.map(cf => [cf, (1 - cf < top ? K.tauMFR(m, 1 - cf) : NaN) / tpOf(cf)]);
          else if (first) { const kk = +v.k; pts = fr.map(cf => { const lt = N.root(q => RT.wehner(kk * Math.exp(q), dd) - cf, Math.log(tpOf(cf) * 0.5), Math.log(tpOf(cf) * 1e5), { n: 60 }); return [cf, Math.exp(lt) / tpOf(cf)]; }); }
          else {
            const tHi = K.tauMFR(m, Math.min(1 - fr[0], top * 0.999)) * 2, tLo = tpOf(fr[fr.length - 1]) * 0.5;
            const taus = N.logspace(tLo, tHi, 46), Xs = [];
            let g = null;
            taus.forEach(t => { const r = RT.dispersionX(m, t, dd, { M: 240, guess: g }); g = r.prof; Xs.push(r.X); });
            pts = fr.map(cf => { const Xg = 1 - cf; if (!(Xg < top) || Xg > Xs[Xs.length - 1] || Xg < Xs[0]) return [cf, NaN]; const t2 = Math.exp(N.interp(Xs, taus.map(Math.log), Xg)); return [cf, t2 / tpOf(cf)]; });
          }
          return { name: dd === 0 ? 'plug flow' : dd === Infinity ? 'mixed flow' : 'D/uL = ' + N.fmt(dd, 3), x: pts.map(p => p[0]), y: pts.map(p => p[1]), color: dd === d ? 2 : i + 1, w: dd === d ? 2.4 : 1.4, dash: dd === Infinity ? '5 4' : undefined };
        });
        plots.push({ title: 'Volume relative to plug flow (Figs. 13.19, 13.20)', x: { label: 'C_A/C_A0 leaving', log: true, min: 1e-3, max: 1 }, y: { label: 'V/V_p', log: true, min: 1 }, series: ser, cap: 'Each line is a fixed dispersion number, computed for this rate law.' });
      }
      return { kpi, warn, info, plots, tables: [{ title: 'Profile', cols: ['z', 'X_A', 'C_A/C_A0'], rows: sol.z.filter((z, i) => i % 10 === 0 || i === sol.z.length - 1).map(z => { const i = sol.z.indexOf(z); return [z, sol.x[i], 1 - sol.x[i]]; }) }] };
    }
  });

  /* ===================================================================== */
  /* F-1303  correlations for axial dispersion                              */
  /* ===================================================================== */
  const wenFan = Re => 3e7 / Math.pow(Re, 2.1) + 1.35 / Math.pow(Re, 0.125);
  const taylorTurb = Re => 3.57 * Math.sqrt(0.0791 * Math.pow(Re, -0.25));
  const taylorAris = Bo => 1 / Bo + Bo / 192;
  const bedGas = (Re, Sc) => { const pm = Re * Sc; return 0.73 / pm + 0.5 / (1 + 9.7 / pm); };
  const bedLiq = (Re, eps) => eps / (0.2 + 0.011 * Math.pow(Re, 0.48));
  ST.tool({
    id: 'dispcorr', tag: 'F-1303', ch: 13, sec: 'Ch. 13.2', title: 'Axial dispersion from correlations: pipes, packed beds, and pipeline interfaces',
    lede: 'D/ud for laminar and turbulent pipe flow and for gases and liquids in packed beds, as in Figs. 13.15–13.17, turned into D/uL for your vessel; the regime check that says whether the dispersion model applies to laminar flow at all; and the contaminated width between two products pumped one after the other.',
    keys: 'dispersion correlation pipe turbulent laminar taylor aris packed bed gas liquid bodenstein reynolds schmidt pipeline interface contaminated width problem 13.4 13.6 13.7 13.12',
    inputs: [
      { group: 'Vessel' },
      { k: 'geo', t: 'seg', label: 'Flow in', opts: [['pipe', 'an empty pipe'], ['bed', 'a packed bed']], v: 'pipe' },
      { k: 'dt', label: 'Pipe diameter d_t', unit: 'm', v: 0.255, min: 1e-6, max: 100, scale: 'log', if: v => v.geo === 'pipe' },
      { k: 'dp', label: 'Particle diameter d_p', unit: 'm', v: 0.005, min: 1e-7, max: 1, scale: 'log', if: v => v.geo === 'bed' },
      { k: 'eps', label: 'Bed voidage', v: 0.4, min: 0.2, max: 0.95, if: v => v.geo === 'bed' },
      { k: 'L', label: 'Length', unit: 'm', v: 293000, min: 1e-6, max: 1e8, scale: 'log' },
      { k: 'u', label: 'Velocity (mean in a pipe, superficial in a bed)', unit: 'm/s', v: 0.817, min: 1e-9, max: 100, scale: 'log' },
      { group: 'Fluid' },
      { k: 'phase', t: 'seg', label: 'Fluid', opts: [['liq', 'liquid'], ['gas', 'gas']], v: 'liq' },
      { k: 'rho', label: 'Density', unit: 'kg/m³', v: 850, min: 1e-3, max: 2e4, scale: 'log' },
      { k: 'mu', label: 'Viscosity', unit: 'Pa·s', v: 7.4e-3, min: 1e-7, max: 1e3, scale: 'log' },
      { k: 'Dm', label: 'Molecular diffusivity 𝒟', unit: 'm²/s', v: 1e-9, min: 1e-14, max: 1e-2, scale: 'log' },
      { k: 'turb', t: 'seg', label: 'Turbulent pipe correlation', opts: [['wf', 'experimental band (Wen & Fan fit)'], ['taylor', 'Taylor theory']], v: 'wf', if: v => v.geo === 'pipe' }
    ],
    presets: [
      { name: 'Problem 13.4 (Rangely–Salt Lake City line)', v: { geo: 'pipe', dt: 0.255, L: 293000, u: 0.817, phase: 'liq', rho: 850, mu: 0.00737, Dm: 1e-9, turb: 'wf' }, note: 'A 10 in. line, 293 km, u = 81.7 cm/s, Re = 24 000. Fluid with more than half the peak tracer passes in about ±1.18σ; the reported passage time was 895 s.' },
      { name: 'Problem 13.6 (refinery products)', v: { geo: 'pipe', dt: 0.1, L: 100000, u: 0.2, phase: 'liq', rho: 850, mu: 1.7e-3, Dm: 1e-9, turb: 'wf' }, note: 'Products A and B pumped in turn through 100 km of 10 cm line at 20 cm/s: the 16%–84% contaminated width.' },
      { name: 'Problem 13.11 (packed bed, liquid)', v: { geo: 'bed', dp: 0.005, eps: 0.4, L: 0.48, u: 0.04, phase: 'liq', rho: 1000, mu: 1e-3, Dm: 1e-9 }, note: '48 cm of 5 mm nonporous pellets with liquid at Re_p = 200 (u_s = 0.04 m/s with water).' },
      { name: 'Laminar gas in a capillary', v: { geo: 'pipe', dt: 0.002, L: 2, u: 0.5, phase: 'gas', rho: 1.2, mu: 1.8e-5, Dm: 2e-5, turb: 'wf' }, note: 'Re ≈ 67, Bo = Re·Sc ≈ 50: the Taylor–Aris line of Fig. 13.16 applies if the tube is long enough.' }
    ],
    eqs: [
      ['\\frac{D}{uL} = \\left(\\frac{D}{ud}\\right)\\left(\\frac{d}{L}\\right)', 'Eq. 13.22'],
      ['D = \\mathscr{D} + \\frac{u^2d_t^2}{192\\,\\mathscr{D}},\\quad \\frac{D}{ud_t} = \\frac{1}{Re\\,Sc} + \\frac{Re\\,Sc}{192}', 'Fig. 13.16', 'Streamline flow (Taylor, Aris).'],
      ['\\frac{D}{ud_t} = \\frac{3\\times10^7}{Re^{2.1}} + \\frac{1.35}{Re^{1/8}}', 'Fig. 13.15', 'Fit to the experimental band for turbulent pipe flow (Wen and Fan, 1975).'],
      ['\\frac{D}{ud_t} = 3.57\\sqrt{f},\\quad f = 0.0791\\,Re^{-1/4}', 'Taylor (1954)'],
      ['\\frac{D}{ud_p} = \\frac{0.73}{Re\\,Sc} + \\frac{0.5}{1+9.7/(Re\\,Sc)}', 'Fig. 13.17', 'Gases in packed beds (Edwards and Richardson, 1968); u interstitial.'],
      ['\\frac{D}{ud_p} = \\frac{\\varepsilon}{0.2 + 0.011\\,Re_p^{0.48}}', 'Fig. 13.17', 'Liquids in packed beds (Chung and Wen, 1968).'],
      ['\\text{width}_{16\\text{–}84\\%} = 2\\sigma_L = 2\\sqrt{2DL/u}', 'Fig. 13.12']
    ],
    method: ['Dimensionless groups: Re = d u ρ/μ (pipe: mean velocity; bed: superficial velocity for Re_p), Sc = μ/(ρ𝒟). In packed beds the dispersion number is referred to the interstitial velocity u/ε, as the book does.', 'The correlations are smooth fits that reproduce the curves and bands of Figs. 13.15–13.17; expect the scatter of those figures (a factor of about two) in any single prediction.', 'For laminar pipe flow the Taylor–Aris result holds only when molecular diffusion has time to even out the radial profile: the tool checks 𝒟t̄/R² against the limits of Fig. 15.2.'],
    run(v) {
      const kpi = [], warn = [], info = [];
      const rho = +v.rho, mu = +v.mu, Dm = +v.Dm, L = +v.L, us = +v.u;
      const Sc = mu / (rho * Dm);
      let Dud, dchar, uref, regime = '';
      if (v.geo === 'pipe') {
        dchar = +v.dt; uref = us;
        const Re = dchar * us * rho / mu;
        kpi.push({ k: 'Re', l: 'Reynolds number', v: Re }, { k: 'Sc', l: 'Schmidt number', v: Sc });
        if (Re < 2100) {
          const Bo = Re * Sc;
          Dud = taylorAris(Bo);
          regime = 'laminar (Taylor–Aris)';
          const tD = Dm * (L / us) / Math.pow(dchar / 2, 2);
          kpi.push({ k: 'Bo', l: 'Bodenstein number Re·Sc', v: Bo }, { k: 'tD', l: '𝒟t̄/R² (radial diffusion time ratio)', v: tD, s: tD > 0.8 ? 'dispersion model applies' : tD < 0.003 ? 'pure convection: use Ch. 15 instead' : 'in between: try both models' });
          if (tD < 0.8) warn.push('Laminar flow with too little time for radial diffusion (𝒟t̄/R² = ' + N.fmt(tD, 3) + '): the dispersion model does not describe this tube; the convection model of Chapter 15 does' + (tD > 0.003 ? ', or average the two' : '') + '.');
        } else {
          Dud = v.turb === 'taylor' ? taylorTurb(Re) : wenFan(Re);
          regime = 'turbulent (' + (v.turb === 'taylor' ? 'Taylor' : 'experimental band') + ')';
        }
      } else {
        dchar = +v.dp; const eps = +v.eps; uref = us / eps;
        const Re = dchar * us * rho / mu;
        kpi.push({ k: 'Re', l: 'Particle Reynolds number (superficial)', v: Re }, { k: 'Sc', l: 'Schmidt number', v: Sc });
        Dud = v.phase === 'gas' ? bedGas(Re / eps, Sc) : bedLiq(Re, eps);
        regime = v.phase === 'gas' ? 'packed bed, gas' : 'packed bed, liquid';
        if (L / dchar < 30) warn.push('Fewer than 30 particle diameters of bed: wall and end effects compete with the dispersion model.');
      }
      const Dax = Dud * uref * dchar, dL = Dud * dchar / L;
      const tb = L / uref;
      kpi.push({ k: 'Dud', l: 'D/ud (' + regime + ')', v: Dud, main: true }, { k: 'D', l: 'Axial dispersion coefficient D', v: Dax, u: 'm²/s' }, { k: 'dL', l: 'Vessel dispersion number D/uL', v: dL, main: true, s: dL < 0.01 ? 'small deviation from plug flow' : 'large: use the closed/open-vessel forms' });
      const sT = tb * Math.sqrt(2 * dL);
      kpi.push({ k: 'tb', l: 'Mean passage time L/u', v: tb, u: 's' }, { k: 'sig', l: 'Spread σ of a pulse at the exit', v: sT, u: 's', s: '±1.18σ above half height: ' + N.fmt(2.36 * sT, 4) + ' s' });
      const width = 2 * Math.sqrt(2 * Dax * L / uref);
      kpi.push({ k: 'width', l: 'Contaminated width 16%–84%', v: width, u: 'm', s: '5%–95%: ' + N.fmt(width * 1.645, 4) + ' m; volume ' + N.fmt(width * (v.geo === 'pipe' ? Math.PI * dchar * dchar / 4 : 1), 4) + (v.geo === 'pipe' ? ' m³' : ' m³ per m² of bed') });
      /* chart of the correlation with the operating point */
      const plots = [];
      if (v.geo === 'pipe') {
        const Res = N.logspace(1, 1e6, 300);
        const lam = Res.filter(r => r < 2100), tur = Res.filter(r => r >= 2100);
        plots.push({ title: 'Dispersion in pipes (Fig. 13.15)', x: { label: 'Re', log: true }, y: { label: 'D/ud_t', log: true, min: 0.05, max: 1e4 }, series: [{ name: 'laminar, Sc = ' + N.fmt(Sc, 3), x: lam, y: lam.map(r => taylorAris(r * Sc)), color: 1 }, { name: 'turbulent, experimental band', x: tur, y: tur.map(wenFan), color: 2 }, { name: 'turbulent, Taylor theory', x: tur, y: tur.map(taylorTurb), color: 3, dash: '5 4' }], marks: [{ t: 'pt', x: dchar * us * rho / mu, y: Dud, color: 8, label: 'this flow' }] });
      } else {
        const Res = N.logspace(0.1, 2000, 200);
        plots.push({ title: 'Dispersion in packed beds (Fig. 13.17)', x: { label: 'Re_p', log: true }, y: { label: 'D/ud_p', log: true, min: 0.1, max: 20 }, series: [{ name: 'liquids', x: Res, y: Res.map(r => bedLiq(r, +v.eps)), color: 1 }, { name: 'gases, Sc = 0.2', x: Res, y: Res.map(r => bedGas(r / +v.eps, 0.2)), color: 2, dash: '5 4' }, { name: 'gases, Sc = 1', x: Res, y: Res.map(r => bedGas(r / +v.eps, 1)), color: 3 }], marks: [{ t: 'pt', x: dchar * us * rho / mu, y: Dud, color: 8, label: 'this flow' }] });
      }
      return { kpi, warn, info, plots };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
