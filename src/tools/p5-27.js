/* Chapter 27 — enzyme fermentation: Michaelis–Menten kinetics in batch, plug, mixed and staged
   reactors, competitive, noncompetitive and mixed inhibition (Eq. 27.15), and the rate constants
   and type of inhibition found from batch or mixed-flow runs by the book's straight-line plots. */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N;
  const fmt = (x, s = 4) => N.fmt(x, s);
  const ENZ = ST.ENZ = {};

  /* apparent constants under inhibition (Eq. 27.15): k3' = k3/(1 + L C_B), C_M' = C_M (1 + N C_B)/(1 + L C_B) */
  ENZ.app = (k3, CM, Nn, L, CB) => ({ k: k3 / (1 + L * CB), CM: CM * (1 + Nn * CB) / (1 + L * CB) });
  /* batch or plug: C_M ln(C_A0/C_A) + (C_A0 − C_A) = k C_E0 t (Eq. 27.5) */
  ENZ.tPlug = (k, CM, CE0, CA0, CA) => (CM * Math.log(CA0 / CA) + (CA0 - CA)) / (k * CE0);
  ENZ.CAplug = (k, CM, CE0, CA0, t) => {
    if (!(t > 0)) return CA0;
    const g = c => ENZ.tPlug(k, CM, CE0, CA0, c) - t;
    return N.brent(g, 1e-300, CA0, { tol: 1e-15 * CA0 });
  };
  /* mixed: k C_E0 τ C_A = (C_A0 − C_A)(C_M + C_A) (Eq. 27.7), positive root */
  ENZ.CAmix = (k, CM, CE0, CA0, tau) => {
    const a = k * CE0 * tau;
    /* C_A² + (C_M − C_A0 + a) C_A − C_A0 C_M = 0 */
    const b = CM - CA0 + a, c = -CA0 * CM;
    return b >= 0 ? -2 * c / (b + Math.sqrt(b * b - 4 * c)) : (-b + Math.sqrt(b * b - 4 * c)) / 2;
  };
  ENZ.tauMix = (k, CM, CE0, CA0, CA) => (CA0 - CA) * (CM + CA) / (k * CE0 * CA);

  ST.tool({
    id: 'enzyme', tag: 'B-2701', ch: 27, sec: 'Ch. 27', title: 'Enzyme fermentation: Michaelis–Menten kinetics, inhibition and reactors',
    lede: 'Michaelis–Menten kinetics in batch, plug, mixed and staged reactors, with an inhibitor acting competitively, noncompetitively or in the general mixed way of Eq. 27.15. Given data from batch or mixed-flow runs the tool draws the book’s straight lines (Eqs. 27.6 and 27.8, which use every point directly), and the Lineweaver–Burk and Eadie plots biologists use, extracts k₃ and C_M, and — with runs at several inhibitor levels — tells competitive from noncompetitive inhibition and gives its constants.',
    keys: 'enzyme fermentation michaelis menten kinetics michaelis constant k3 CE0 batch plug mixed flow fermentor inhibition competitive noncompetitive uncompetitive mixed inhibitor constant lineweaver burk eadie hofstee straight line plot rate constants sucrose sucrase cellulose cellulase cellobiose glucose problem 27.1 27.2 27.3 27.4 27.5 27.6 27.7 27.8',
    units: { t: 'h' },
    inputs: [
      { group: 'Task' },
      { k: 'mode', t: 'seg', label: 'Do', opts: [['design', 'reactor performance'], ['fit', 'rate equation from runs']], v: 'design' },
      { group: 'Kinetics, −r_A = k₃C_E0C_A/(C_M + C_A)', if: v => v.mode === 'design' },
      { k: 'k3', label: 'k₃', unit: '1/[t]', v: 20, min: 1e-12, max: 1e12, scale: 'log', if: v => v.mode === 'design' },
      { k: 'CM', label: 'Michaelis constant C_M', unit: 'mol/m³', v: 0.2, min: 0, max: 1e9, if: v => v.mode === 'design' },
      { k: 'CE0', label: 'Total enzyme C_E0', unit: 'mol/m³', v: 0.01, min: 1e-15, max: 1e9, scale: 'log', if: v => v.mode === 'design' },
      { k: 'CA0', label: 'C_A0', unit: 'mol/m³', v: 1, min: 1e-12, max: 1e9, if: v => v.mode === 'design' },
      { k: 'inh', t: 'select', label: 'Inhibitor B', opts: [['none', 'none'], ['comp', 'competitive: C_M → C_M(1 + N C_B0)'], ['non', 'noncompetitive (Eq. 27.15, N = L)'], ['mixed', 'general (Eq. 27.15, any N and L)'], ['unc', 'acts on the complex only (N = 0)']], v: 'none', if: v => v.mode === 'design' },
      { k: 'CB0', label: 'C_B0', unit: 'mol/m³', v: 10, min: 0, max: 1e9, if: v => v.mode === 'design' && v.inh !== 'none' },
      { k: 'Ni', label: 'N = k₄/k₅ (B on free enzyme)', unit: 'm³/mol', v: 0.05, min: 0, max: 1e9, if: v => v.mode === 'design' && (v.inh === 'comp' || v.inh === 'non' || v.inh === 'mixed') },
      { k: 'Li', label: 'L = k₆/k₇ (B on the complex)', unit: 'm³/mol', v: 0.05, min: 0, max: 1e9, if: v => v.mode === 'design' && (v.inh === 'mixed' || v.inh === 'unc') },
      { group: 'Reactor', if: v => v.mode === 'design' },
      { k: 'find', t: 'seg', label: 'Find', opts: [['X', 'conversion for a time'], ['t', 'time for a conversion']], v: 'X', if: v => v.mode === 'design' },
      { k: 'tt', label: 'Batch time or space time τ', unit: '[t]', v: 6, min: 0, max: 1e12, if: v => v.mode === 'design' && v.find === 'X' },
      { k: 'Xt', label: 'Conversion of A', v: 0.95, min: 0.0001, max: 0.999999, if: v => v.mode === 'design' && v.find === 't' },
      { k: 'Nst', label: 'Mixed stages in series (for the comparison)', v: 3, min: 1, max: 50, step: 1, if: v => v.mode === 'design' },
      { group: 'Runs', if: v => v.mode === 'fit' },
      { k: 'rtype', t: 'seg', label: 'Runs made in', opts: [['batch', 'batch or plug flow (t)'], ['mixed', 'mixed flow (τ)']], v: 'batch', if: v => v.mode === 'fit' },
      { k: 'runs', t: 'table', label: 'Runs', cols: ['C_B0', 'C_E0', 'C_A0', 'C_A', 't or τ'], v: [[0, 0.01, 1, 0.68, 2], [0, 0.01, 1, 0.16, 6], [0, 0.01, 1, 0.006, 10]], if: v => v.mode === 'fit' },
      { k: 'flowIn', t: 'check', label: 'Last column is a flow rate v (τ = V/v)', v: false, if: v => v.mode === 'fit' && v.rtype === 'mixed' },
      { k: 'Vr', label: 'Reactor volume V', unit: 'same volume as v', v: 6, min: 1e-12, max: 1e12, if: v => v.mode === 'fit' && v.rtype === 'mixed' && v.flowIn }
    ],
    presets: [
      { name: 'Problem 27.2 — sucrose by sucrase, batch', v: { mode: 'fit', rtype: 'batch', runs: [[0, 0.01, 1, 0.68, 2], [0, 0.01, 1, 0.16, 6], [0, 0.01, 1, 0.006, 10]], _tu: 'h' }, check: { k3: { v: 20, tol: 0.03 }, CM: { v: 0.2, tol: 0.08 } }, note: 'C_A0 = 1 mol/m³, C_E0 = 0.01 mol/m³; the three readings fall on one straight line of Eq. 27.6: k₃ ≈ 20 hr⁻¹ and C_M ≈ 0.2 mol/m³ (the rounded readings give 19.6 and 0.19).' },
      { name: 'Problem 27.1 — mixed-flow runs', v: { mode: 'fit', rtype: 'mixed', flowIn: true, Vr: 6, runs: [[0, 0.02, 0.2, 0.04, 3], [0, 0.01, 0.3, 0.15, 4], [0, 0.001, 0.69, 0.6, 1.2]], _tu: 'h' }, note: 'V = 6 L, flows in L/hr, concentrations in mol/L.' },
      { name: 'Problem 27.3 — batch runs, different enzyme levels', v: { mode: 'fit', rtype: 'batch', runs: [[0, 3, 400, 10, 1], [0, 2, 200, 5, 1], [0, 1, 20, 1, 1]], _tu: 'h' } },
      { name: 'Problem 27.4 — mixed flow, with and without B', v: { mode: 'fit', rtype: 'mixed', flowIn: true, Vr: 240, runs: [[0, 12.5, 200, 50, 80], [0, 5, 900, 300, 24], [0, 5, 1200, 800, 48], [33.3, 33.3, 700, 33.3, 24], [33.3, 10, 200, 80, 80], [33.3, 20, 900, 500, 120]], _tu: 'min' }, check: { k3: { v: 20, tol: 0.02 }, CM: { v: 200, tol: 0.03 }, type: 'noncompetitive' }, note: 'V = 240 cm³, flows in cm³/min. Without B the three runs give k₃ = 20 min⁻¹ and C_M = 200 mol/m³; with 33.3 mol/m³ of B, C_M is unchanged and k₃ falls to 14.' },
      { name: 'Problem 27.5 — batch runs with and without B', v: { mode: 'fit', rtype: 'batch', runs: [[0, 8, 600, 350, 1], [0, 8, 600, 160, 2], [0, 8, 600, 40, 3], [0, 8, 600, 10, 4], [100, 8, 800, 560, 1], [100, 8, 800, 340, 2], [100, 8, 800, 180, 3], [100, 8, 800, 80, 4], [100, 8, 800, 30, 5]], _tu: 'h' }, check: { type: 'competitive' }, note: 'C_E0 = 8 g/m³ in both runs; 100 mol/m³ of B in the second leaves the slope k₃ unchanged and raises C_M.' },
      { name: 'Problems 27.6–27.8 — cellulose, two inhibitors', v: { mode: 'fit', rtype: 'mixed', flowIn: false, runs: [[0, 0.01, 25, 1.5, 587], [0, 0.01, 25, 4.5, 279], [0, 0.01, 25, 9, 171], [0, 0.01, 25, 21, 36], [5, 0.01, 25, 1.5, 940], [5, 0.01, 25, 4.5, 387], [5, 0.01, 25, 9, 213], [5, 0.01, 25, 21, 40], [10, 0.01, 25, 1.5, 1020], [10, 0.01, 25, 4.5, 433], [10, 0.01, 25, 9, 250], [10, 0.01, 25, 21, 50]], _tu: 'min' }, note: 'Cellulose at 25 kg/m³ and 50 °C, C_E0 = 0.01 kg/m³. Series 2 (C_B0 = 5) is cellobiose and series 3 (entered as 10) is glucose; each is compared with the uninhibited series.' },
      { name: 'Sucrose hydrolysis: batch, mixed and staged', v: { mode: 'design', k3: 20, CM: 0.2, CE0: 0.01, CA0: 1, inh: 'none', find: 't', Xt: 0.99, Nst: 3, _tu: 'h' } },
      { name: 'Competitive inhibitor', v: { mode: 'design', k3: 20, CM: 200, CE0: 5, CA0: 900, inh: 'comp', CB0: 33.3, Ni: 0.03, find: 'X', tt: 10, Nst: 3, _tu: 'min' } },
      { name: 'Noncompetitive inhibitor (Problem 27.4 constants)', v: { mode: 'design', k3: 20, CM: 200, CE0: 5, CA0: 900, inh: 'non', CB0: 33.3, Ni: 0.0129, find: 'X', tt: 10, Nst: 3, _tu: 'min' } }
    ],
    eqs: [
      ['-r_A = r_R = \\frac{k_3C_{E0}C_A}{C_M + C_A},\\qquad C_M = \\frac{k_2+k_3}{k_1}', 'Eqs. 27.3–27.4'],
      ['C_M\\ln\\frac{C_{A0}}{C_A} + (C_{A0}-C_A) = k_3C_{E0}t', 'Eq. 27.5', 'Batch or plug flow.'],
      ['\\frac{C_{A0}-C_A}{\\ln(C_{A0}/C_A)} = -C_M + k_3\\frac{C_{E0}t}{\\ln(C_{A0}/C_A)}', 'Eq. 27.6', 'Straight line from batch data.'],
      ['k_3C_{E0}\\tau = \\frac{(C_{A0}-C_A)(C_M+C_A)}{C_A},\\qquad C_A = -C_M + k_3\\frac{C_{E0}C_A\\tau}{C_{A0}-C_A}', 'Eqs. 27.7–27.8', 'Mixed flow, and its straight line.'],
      ['r_R = \\frac{k_3C_{E0}C_A}{C_M(1+NC_{B0}) + C_A},\\qquad N = \\frac{k_4}{k_5}', 'Eq. 27.11', 'Competitive inhibition.'],
      ['r_R = \\frac{\\frac{k_3}{1+LC_{B0}}C_{E0}C_A}{C_M\\frac{1+NC_{B0}}{1+LC_{B0}} + C_A},\\qquad L = \\frac{k_6}{k_7}', 'Eq. 27.15', 'Noncompetitive inhibition.'],
      ['\\frac{1}{-r_A} = \\frac{1}{k_3C_{E0}} + \\frac{C_M}{k_3C_{E0}}\\frac{1}{C_A},\\qquad -r_A = k_3C_{E0} - C_M\\frac{-r_A}{C_A}', 'Ch. 27', 'Lineweaver–Burk and Eadie plots.']
    ],
    method: ['Batch and plug flow use Eq. 27.5, inverted for C_A by Brent’s method; mixed flow solves the quadratic of Eq. 27.7 for its positive root; mixed stages in series apply it stage by stage. An inhibitor enters only through the apparent constants of Eq. 27.15 — competitive inhibition is L = 0, the classical noncompetitive case is N = L, and N = 0 is the case where B attacks only the enzyme–substrate complex.', 'Fitting: within each inhibitor level every reading is turned into a point of Eq. 27.6 (batch) or 27.8 (mixed) and a least-squares line gives k₃′ (slope) and C_M′ (minus the intercept). Comparing each inhibited level with the uninhibited one gives 1 + LC_B0 = k₃/k₃′ and 1 + NC_B0 = (C_M′/C_M)(1 + LC_B0). Equal slopes (within 8%) signal competitive inhibition.'],
    run(v) {
      const kpi = [], warn = [], info = [], plots = [], tables = [];
      const tu = v._tu || 'h';
      if (v.mode === 'fit') {
        const raw = (v.runs || []).filter(r => r && r.slice(0, 5).every(z => z !== '' && isFinite(+z))).map(r => ({ CB: +r[0], CE: +r[1], CA0: +r[2], CA: +r[3], t: +r[4] }));
        if (raw.length < 2) throw new Error('Enter at least two runs.');
        raw.forEach(r => { if (v.rtype === 'mixed' && v.flowIn) r.t = +v.Vr / r.t; });
        const bad = raw.filter(r => !(r.CA > 0 && r.CA < r.CA0));
        if (bad.length) throw new Error('Every run needs 0 < C_A < C_A0.');
        const groups = {};
        raw.forEach(r => { (groups[r.CB] = groups[r.CB] || []).push(r); });
        const lv = Object.keys(groups).map(Number).sort((a, b) => a - b);
        const fits = lv.map(CB => {
          const g = groups[CB];
          const pts = g.map(r => v.rtype === 'batch' ? { x: r.CE * r.t / Math.log(r.CA0 / r.CA), y: (r.CA0 - r.CA) / Math.log(r.CA0 / r.CA) } : { x: r.CE * r.CA * r.t / (r.CA0 - r.CA), y: r.CA });
          const L = g.length >= 2 ? N.linreg(pts.map(p => p.x), pts.map(p => p.y)) : null;
          return { CB, g, pts, L, k: L ? L.b : NaN, CM: L ? -L.a : NaN };
        });
        const ref = fits.find(f => f.CB === 0) || fits[0];
        if (!ref.L) throw new Error('The reference level needs at least two runs.');
        kpi.push({ k: 'k3', l: 'k₃' + (ref.CB > 0 ? ' at C_B0 = ' + fmt(ref.CB, 3) : ''), v: ref.k, u: '1/' + tu, main: true, s: '± ' + fmt(ref.L.sb, 2) + (ref.L.r2 !== undefined ? ', R² = ' + fmt(ref.L.r2, 4) : '') });
        kpi.push({ k: 'CM', l: 'Michaelis constant C_M', v: ref.CM, u: 'same units as C_A', s: '± ' + fmt(ref.L.sa, 2) });
        if (ref.CM < 0) warn.push('The intercept gives a negative C_M: the Michaelis–Menten form does not fit these runs well (or C_M is too small to resolve).');
        /* inhibition */
        const rows = [];
        let type = '';
        fits.filter(f => f !== ref && f.L).forEach(f => {
          const rk = f.k / ref.k, rm = f.CM / ref.CM;
          let tp, Nn = NaN, Li = NaN;
          if (Math.abs(rk - 1) < 0.08 && Math.abs(rm - 1) < 0.08) tp = 'no effect';
          else if (Math.abs(rk - 1) < 0.08) { tp = 'competitive'; Nn = (rm - 1) / f.CB; Li = 0; }
          else {
            Li = (1 / rk - 1) / f.CB; Nn = (rm * (1 + Li * f.CB) - 1) / f.CB;
            tp = Math.abs(rm - 1) < 0.08 ? 'noncompetitive' : Math.abs(Nn) < 0.1 * Math.abs(Li) ? 'on the complex only' : 'mixed (Eq. 27.15)';
          }
          type = type || tp;
          rows.push([f.CB, f.k, f.CM, rk, rm, tp, Nn, Li]);
        });
        if (rows.length) {
          kpi.push({ k: 'type', l: 'Role of B', v: type, s: rows.length > 1 ? 'first inhibitor level; see the table for all' : '' });
          const r0 = rows[0];
          if (isFinite(r0[6])) kpi.push({ k: 'Ninh', l: 'N = k₄/k₅', v: r0[6], u: '1/(conc.)' });
          if (isFinite(r0[7]) && r0[7] !== 0) kpi.push({ k: 'Linh', l: 'L = k₆/k₇', v: r0[7], u: '1/(conc.)' });
          tables.push({ title: 'Each inhibitor level against the uninhibited runs (Figs. 27.8–27.9)', cols: ['C_B0', 'k₃′', 'C_M′', 'k₃′/k₃', 'C_M′/C_M', 'inhibition', 'N', 'L'], rows });
          info.push('Competitive: the lines stay parallel (same k₃) and move right. Noncompetitive: the slope k₃ falls by 1 + LC_B0 and C_M becomes C_M(1 + NC_B0)/(1 + LC_B0) (Eq. 27.15).');
        }
        tables.push({ title: 'Lines of Eq. ' + (v.rtype === 'batch' ? '27.6' : '27.8'), cols: ['C_B0', 'runs', 'slope k₃′', 'C_M′ = −intercept', 'R²'], rows: fits.map(f => [f.CB, f.g.length, f.k, f.CM, f.L ? f.L.r2 : NaN]) });
        /* plot of the straight lines */
        const ser = [], xmax = Math.max(...fits.flatMap(f => f.pts.map(p => p.x))) * 1.1;
        fits.forEach((f, i) => {
          ser.push({ name: 'C_B0 = ' + fmt(f.CB, 3), x: f.pts.map(p => p.x), y: f.pts.map(p => p.y), type: 'pts', color: i + 1 });
          if (f.L) ser.push({ x: [0, xmax], y: [f.L.a, f.L.a + f.L.b * xmax], color: i + 1, legend: false, dash: '5 4' });
        });
        plots.push({ title: v.rtype === 'batch' ? 'Eq. 27.6: (C_A0 − C_A)/ln(C_A0/C_A) against C_E0t/ln(C_A0/C_A) (Fig. 27.5)' : 'Eq. 27.8: C_A against C_E0C_Aτ/(C_A0 − C_A) (Fig. 27.6)', x: { label: v.rtype === 'batch' ? 'C_E0t / ln(C_A0/C_A)' : 'C_E0C_Aτ / (C_A0 − C_A)', min: 0 }, y: { label: v.rtype === 'batch' ? '(C_A0 − C_A) / ln(C_A0/C_A)' : 'C_A' }, series: ser, marks: [{ t: 'h', y: 0, color: 'var(--muted)' }], cap: 'Slope k₃; the line cuts the vertical axis at −C_M.' });
        /* Lineweaver–Burk and Eadie for mixed data (rates are measured directly) */
        if (v.rtype === 'mixed') {
          const lb = [], ea = [];
          fits.forEach((f, i) => {
            const rr = f.g.map(r => ({ r: (r.CA0 - r.CA) / (r.t * r.CE), C: r.CA }));
            lb.push({ name: 'C_B0 = ' + fmt(f.CB, 3), x: rr.map(q => 1 / q.C), y: rr.map(q => 1 / q.r), type: 'pts', color: i + 1 });
            ea.push({ name: 'C_B0 = ' + fmt(f.CB, 3), x: rr.map(q => q.r / q.C), y: rr.map(q => q.r), type: 'pts', color: i + 1 });
            if (f.L) {
              const xs = [0, Math.max(...rr.map(q => 1 / q.C)) * 1.1];
              lb.push({ x: xs, y: xs.map(x => 1 / f.k + f.CM / f.k * x), color: i + 1, legend: false, dash: '5 4' });
              const xe = [0, Math.max(...rr.map(q => q.r / q.C)) * 1.1];
              ea.push({ x: xe, y: xe.map(x => f.k - f.CM * x), color: i + 1, legend: false, dash: '5 4' });
            }
          });
          plots.push({ title: 'Lineweaver–Burk plot (rates per unit enzyme)', x: { label: '1/C_A', min: 0 }, y: { label: 'C_E0/(−r_A)', min: 0 }, series: lb, cap: 'Dashed lines are the Eq. 27.8 fits, not separate fits: the book’s method uses the data without forming rates first.' });
          plots.push({ title: 'Eadie plot', x: { label: '(−r_A)/(C_E0C_A)', min: 0 }, y: { label: '(−r_A)/C_E0', min: 0 }, series: ea });
        }
        return { kpi, warn, info, plots, tables };
      }
      /* ---------- design ---------- */
      const k3 = +v.k3, CM = +v.CM, CE0 = +v.CE0, CA0 = +v.CA0, CB = v.inh === 'none' ? 0 : +v.CB0;
      const Nn = v.inh === 'unc' || v.inh === 'none' ? 0 : +v.Ni, Li = v.inh === 'mixed' || v.inh === 'unc' ? +v.Li : v.inh === 'non' ? +v.Ni : 0;
      const A = ENZ.app(k3, CM, Nn, Li, CB), k = A.k, cm = A.CM;
      const Nst = Math.max(1, Math.round(+v.Nst));
      const stagesOut = tau => { let c = CA0; for (let i = 0; i < Nst; i++) c = ENZ.CAmix(k, cm, CE0, c, tau / Nst); return c; };
      if (v.inh !== 'none') kpi.push({ k: 'app', l: 'Apparent constants with B', v: k, u: '1/' + tu, s: 'k₃′ (was ' + fmt(k3, 4) + '); C_M′ = ' + fmt(cm, 4) + ' (was ' + fmt(CM, 4) + ')' });
      if (v.find === 'X') {
        const t = +v.tt;
        const Cp = ENZ.CAplug(k, cm, CE0, CA0, t), Cm = ENZ.CAmix(k, cm, CE0, CA0, t), Cn = stagesOut(t);
        kpi.push({ k: 'Xp', l: 'Conversion, batch or plug flow', v: 1 - Cp / CA0, main: true, s: 'C_A = ' + fmt(Cp, 4) });
        kpi.push({ k: 'Xm', l: 'Conversion, mixed flow', v: 1 - Cm / CA0, s: 'C_A = ' + fmt(Cm, 4) });
        kpi.push({ k: 'Xn', l: 'Conversion, ' + Nst + ' mixed stages', v: 1 - Cn / CA0 });
        if (v.inh !== 'none') { const C0 = ENZ.CAplug(k3, CM, CE0, CA0, t); kpi.push({ k: 'X0', l: 'Batch conversion without the inhibitor', v: 1 - C0 / CA0 }); }
      } else {
        const Xt = +v.Xt, CA = CA0 * (1 - Xt);
        const tp = ENZ.tPlug(k, cm, CE0, CA0, CA), tm = ENZ.tauMix(k, cm, CE0, CA0, CA);
        const tn = N.brent(tt => stagesOut(tt) - CA, tp * 0.999, tm * 1.001 + tp, { tol: 1e-12 * tm });
        kpi.push({ k: 'tp', l: 'Batch time or plug-flow τ for X_A = ' + fmt(Xt, 4), v: tp, u: tu, main: true });
        kpi.push({ k: 'tm', l: 'Mixed-flow τ', v: tm, u: tu, s: fmt(tm / tp, 3) + ' × plug flow' });
        kpi.push({ k: 'tn', l: 'Total τ, ' + Nst + ' mixed stages', v: tn, u: tu, s: fmt(tn / tp, 3) + ' × plug flow' });
      }
      kpi.push({ k: 'rmax', l: 'Maximum rate r_max = k₃′C_E0', v: k * CE0, u: 'conc./' + tu, s: 'half of it at C_A = C_M′ = ' + fmt(cm, 4) });
      /* C_A against time, Fig. 27.4 */
      const tEnd = v.find === 'X' ? Math.max(+v.tt, 1e-9) * 1.3 : ENZ.tPlug(k, cm, CE0, CA0, CA0 * (1 - +v.Xt)) * 1.3;
      const ts = N.linspace(0, tEnd, 160);
      const ser = [{ name: 'batch or plug (Eq. 27.5)', x: ts, y: ts.map(t => ENZ.CAplug(k, cm, CE0, CA0, t)), color: 1, w: 2.2 }, { name: 'mixed flow, τ = t', x: ts, y: ts.map(t => ENZ.CAmix(k, cm, CE0, CA0, t)), color: 2, dash: '5 4' }];
      ser.push({ name: 'zero order (C_A ≫ C_M)', x: ts, y: ts.map(t => Math.max(CA0 - k * CE0 * t, NaN)), color: 'var(--muted)', dash: '2 3', noRange: true });
      ser.push({ name: 'first order (C_A ≪ C_M)', x: ts, y: ts.map(t => CA0 * Math.exp(-k * CE0 * t / Math.max(cm, 1e-300))), color: 'var(--muted)', dash: '6 3', noRange: true });
      if (v.inh !== 'none') ser.push({ name: 'batch without the inhibitor', x: ts, y: ts.map(t => ENZ.CAplug(k3, CM, CE0, CA0, t)), color: 3, dash: '3 3' });
      plots.push({ title: 'Concentration against time (Fig. 27.4)', x: { label: 't or τ', unit: tu, min: 0, max: tEnd }, y: { label: 'C_A', min: 0, max: CA0 * 1.02 }, series: ser, cap: 'Zero order while C_A ≫ C_M, first order once C_A ≪ C_M; the shift comes near C_A = C_M.' });
      const cs = N.linspace(0, Math.max(CA0, 4 * cm) * 1.05, 160);
      const rser = [{ name: 'rate', x: cs, y: cs.map(c => k * CE0 * c / (cm + c)), color: 1, w: 2.2 }];
      if (v.inh !== 'none') rser.push({ name: 'without the inhibitor', x: cs, y: cs.map(c => k3 * CE0 * c / (CM + c)), color: 3, dash: '5 4' });
      plots.push({ title: 'Rate against concentration (Fig. 27.3)', x: { label: 'C_A', min: 0 }, y: { label: '−r_A', unit: 'conc./' + tu, min: 0 }, series: rser, marks: [{ t: 'h', y: k * CE0, color: 'var(--muted)', label: 'r_max' }, { t: 'v', x: cm, color: 'var(--muted)', label: 'C_M' }] });
      const Xs = N.linspace(0.001, 0.99, 140);
      plots.push({ title: 'Levenspiel plot: C_A0/(−r_A) against X_A', x: { label: 'X_A', min: 0, max: 1 }, y: { label: 'C_A0/(−r_A)', unit: tu, min: 0 }, series: [{ name: 'C_A0/(−r_A)', x: Xs, y: Xs.map(X => { const c = CA0 * (1 - X); return CA0 * (cm + c) / (k * CE0 * c); }), color: 1, w: 2.2 }], cap: 'The rate only falls with conversion, so plug flow always needs less space time than mixed flow.' });
      return { kpi, warn, info, plots, tables };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
