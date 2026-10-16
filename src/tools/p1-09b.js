/* Chapter 9 — nonisothermal reactor design on the X–T chart (Examples 9.4–9.7 and beyond) */
(function (G) {
  'use strict';
  const ST = G.ST, N = ST.N, TH = ST.TH, D = ST.D, XTI = ST.XTI;
  const C2K = TH.C2K, K2C = TH.K2C;

  /* left-to-right flowsheet: items [{hx:'label'} | {tank:{X,tag,label,jacket}} | {tube:{prof,tag,label,jacket}} | {text}] */
  function flowsheet(items, o = {}) {
    const y = 58;
    let x = 34, s = '';
    const xs = [];
    s += D.pipe([[4, y], [x - 4, y]], { label: o.feed || 'feed', dy: -8 });
    items.forEach((it, i) => {
      let w;
      if (it.hx !== undefined) { w = 26; s += D.hx(x + 13, y, 12, { label: it.hx }); }
      else if (it.tank) { w = 56; s += D.tank(x, y - 30, w, 64, { fill: it.tank.X, tag: it.tank.tag, label: it.tank.label, jacket: it.tank.jacket }); }
      else if (it.tube) { w = it.tube.w || 110; s += D.tube(x, y - 11, w, 22, { prof: it.tube.prof, tag: it.tube.tag, label: it.tube.label }); if (it.tube.jacket) s += '<rect class="jacket" x="' + (x + 4) + '" y="' + (y - 17) + '" width="' + (w - 8) + '" height="34" rx="12" style="fill:none"/>'; }
      xs.push([x, x + w]);
      x += w;
      const gap = it.gap || 30;
      s += D.pipe([[x + 2, y], [x + gap - 4, y]], { label: it.after, dy: -8, at: 0, dx: 2 });
      x += gap;
    });
    if (o.recycle) {
      const [a, b] = o.recycle, x0 = xs[a][0] - 12, x1 = xs[b][1] + 14;
      s += D.pipe([[x1, y], [x1, y + 46], [x0, y + 46], [x0, y + 3]], { label: o.recLabel || 'recycle', at: 2, dx: 8, dy: -6, anchor: 'start' });
      s += D.node(x1, y) + D.node(x0, y);
    }
    return D.svg(x + 30, o.recycle ? 150 : 128, s, o.aria || 'flowsheet');
  }

  XTI.flowsheet = flowsheet;

  const EX94 = { form: 'A=R', rev: true, phase: 'liq', CA0: 4, MR: 0, ksrc: 'arr', A1: Math.exp(17.34), E1: 48900, Ksrc: 'arr2', A2: Math.exp(42.04), E2: 124200, dH: -75300, dCp: 0, Cp: 1046, Tlo: 5, Thi: 95 };

  ST.tool({
    id: 'nonisodesign', tag: 'H-903', ch: 9, sec: 'Ch. 9.1', title: 'Nonisothermal design: optimum progression, mixed flow, adiabatic, recycle, cooled and staged',
    lede: 'Every way of meeting one duty with one feed, sized on the same X–T chart: the optimal temperature progression, mixed flow on the locus of maximum rates with its pre- and post-coolers, the best adiabatic plug flow, adiabatic recycle, isothermal reactors, a wall-cooled tube, N mixed reactors at their best temperatures, and N adiabatic beds with interstage cooling.',
    keys: 'optimum temperature progression otp adiabatic plug flow recycle mixed flow locus maximum rates heat duty precooler postcooler interstage cooling multistage staged adiabatic konoki cooled tube hot spot examples 9.4 9.5 9.6 9.7 nonisothermal reactor volume',
    units: { t: 'min' }, slow: true,
    inputs: XTI.inputs(EX94).concat([
      { group: 'Duty' },
      { k: 'FA0', label: '$F_{A0}$', unit: 'mol/[t]', v: 1000, min: 1e-6, max: 1e9, scale: 'log' },
      { k: 'Xf', label: 'Conversion wanted, $X_{Af}$', v: 0.8, min: 0.001, max: 0.9999 },
      { k: 'Tfeed', label: 'Feed available at', unit: '°C', v: 25, min: -150, max: 2000 },
      { k: 'Tprod', label: 'Product delivered at', unit: '°C', v: 25, min: -150, max: 2000 },
      { group: 'Schemes to size' },
      { k: 'sOTP', t: 'check', label: 'Plug flow, optimum temperature progression', v: true },
      { k: 'sMFR', t: 'check', label: 'Mixed flow on the locus of maximum rates', v: true },
      { k: 'sAD', t: 'check', label: 'Adiabatic plug flow', v: true },
      { k: 'adMode', t: 'seg', label: 'Adiabatic feed temperature', opts: [['best', 'best'], ['given', 'given']], v: 'best', if: v => v.sAD },
      { k: 'T1ad', label: 'Adiabatic feed at', unit: '°C', v: 16, min: -150, max: 2000, if: v => v.sAD && v.adMode === 'given' },
      { k: 'sREC', t: 'check', label: 'Adiabatic plug flow with recycle', v: true },
      { k: 'recMode', t: 'seg', label: 'Recycle', opts: [['line', 'best R on the adiabatic line'], ['both', 'best R and feed temperature'], ['given', 'given R']], v: 'line', if: v => v.sREC },
      { k: 'Rgiven', label: 'Recycle ratio R', v: 1, min: 0, max: 1000, scale: 'log', if: v => v.sREC && v.recMode === 'given' },
      { k: 'sISO', t: 'check', label: 'Isothermal plug and mixed flow', v: false },
      { k: 'Tiso', label: 'Isothermal at', unit: '°C', v: 40, min: -150, max: 2000, if: v => v.sISO },
      { k: 'sCOOL', t: 'check', label: 'Wall-cooled plug flow', v: false },
      { k: 'T1c', label: 'Cooled tube: feed at', unit: '°C', v: 60, min: -150, max: 2000, if: v => v.sCOOL },
      { k: 'Ua', label: 'Cooled tube: Ua per reactor volume', unit: 'W/(L·K)', v: 2, min: 1e-4, max: 1e4, scale: 'log', if: v => v.sCOOL },
      { k: 'Tc', label: 'Coolant at', unit: '°C', v: 20, min: -150, max: 2000, if: v => v.sCOOL },
      { k: 'sNM', t: 'check', label: 'Mixed reactors in series, each at its best temperature', v: false },
      { k: 'Nm', label: 'Number of mixed reactors', v: 3, min: 2, max: 8, step: 1, if: v => v.sNM },
      { k: 'sNA', t: 'check', label: 'Adiabatic stages with interstage cooling', v: false },
      { k: 'Na', label: 'Number of adiabatic stages', v: 3, min: 1, max: 6, step: 1, if: v => v.sNA },
      { k: 'focus', t: 'select', label: 'Flowsheet and profiles for', opts: [['otp', 'optimum progression'], ['mfr', 'mixed flow'], ['ad', 'adiabatic'], ['rec', 'recycle'], ['iso', 'isothermal'], ['cool', 'cooled tube'], ['nm', 'mixed reactors in series'], ['na', 'staged adiabatic']], v: 'mfr' }
    ]),
    presets: [
      { name: 'Examples 9.4–9.7', v: Object.assign({}, EX94, { FA0: 1000, Xf: 0.8, Tfeed: 25, Tprod: 25, sOTP: true, sMFR: true, sAD: true, adMode: 'best', sREC: true, recMode: 'line', sISO: false, sCOOL: false, sNM: true, Nm: 3, sNA: true, Na: 3, focus: 'mfr', _tu: 'min' }),
        book: { VOTP: { v: 405, tol: 0.03 }, XB: { v: 0.27, tol: 0.04 }, VMFR: { v: 2000, tol: 0.03 }, rM: { v: 0.4, tol: 0.03 }, Qpre: { v: 348.7, tol: 0.04 }, Qpost: { v: 645, tol: 0.03 }, T1ad: { v: 16.0, tol: 0.03 }, Toutad: { v: 73.6, tol: 0.01 }, VAD: { v: 1720, tol: 0.03 }, VREC: { v: 1200, tol: 0.03 } },
        note: 'C_A0 = 4 mol/L, F_A0 = 1000 mol/min, 80% conversion, 5–95 °C allowed, C′_p = 250 cal/mol A·K. Book: 405 L (optimum progression), 2000 L (mixed flow), 1720 L (adiabatic plug flow fed at 16 °C, leaving at 73.6 °C) and 1200 L (adiabatic with recycle).' },
      { name: 'Problem 9.1 (60%)', v: Object.assign({}, EX94, { FA0: 1000, Xf: 0.6, sOTP: true, sMFR: false, sAD: false, sREC: false, focus: 'otp' }), note: 'Optimum progression in plug flow to 60% conversion, and the exit temperature.' },
      { name: 'Problem 9.2 (70%, mixed)', v: Object.assign({}, EX94, { FA0: 1000, Xf: 0.7, sOTP: false, sMFR: true, sAD: false, sREC: false, focus: 'mfr' }), note: 'Smallest mixed flow reactor for 70% conversion: inlet and outlet temperatures and τ.' },
      { name: 'Problem 9.4 (40 °C, 90%)', v: Object.assign({}, EX94, { FA0: 1000, Xf: 0.9, sOTP: true, sMFR: false, sAD: false, sREC: false, sISO: true, Tiso: 40, focus: 'iso' }), note: 'Plug flow held at 40 °C throughout to 90% conversion.' },
      { name: 'Problem 9.6 (C_A0 = 1)', v: Object.assign({}, EX94, { CA0: 1, FA0: 1000, Xf: 0.8, sOTP: true, sMFR: true, sAD: true, sREC: true, focus: 'mfr' }), note: 'Examples 9.4–9.7 redone with a feed four times more dilute at the same molar feed rate.' },
      { name: 'Problem 9.9 (mixed, 95%, C_A0 = 10)', v: Object.assign({}, EX94, { CA0: 10, FA0: 1000, Xf: 0.95, sOTP: false, sMFR: true, sAD: false, sREC: false, focus: 'mfr' }), note: 'v = 100 L/min of a 10 mol/L feed, 95% conversion in mixed flow.' },
      { name: 'Cooled tube and hot spot', v: Object.assign({}, EX94, { FA0: 1000, Xf: 0.8, sOTP: true, sMFR: false, sAD: true, sREC: false, sCOOL: true, T1c: 60, Ua: 2, Tc: 20, focus: 'cool' }), note: 'Wall cooling bends the adiabatic line toward the locus of maximum rates. Raise Ua or lower the coolant temperature and watch the hot spot and volume.' }
    ],
    eqs: [
      ['\\frac{V}{F_{A0}} = \\int_0^{X_{Af}} \\frac{dX_A}{(-r_A)_{\\max}}', 'Ex. 9.4', 'Optimum temperature progression.'],
      ['\\frac{V}{F_{A0}} = \\frac{X_{Af}}{(-r_A)_{f}}', 'Eq. 5.11', 'Mixed flow, operated on the locus of maximum rates.'],
      ['X_A = \\frac{C\'_p\\,(T-T_1)}{-\\Delta H_r}', 'Eq. 24', 'Adiabatic operating line; slope C′_p/(−ΔH_r).'],
      ['\\frac{V}{F_{A0}} = (R+1)\\int_{\\frac{R}{R+1}X_{Af}}^{X_{Af}} \\frac{dX_A}{-r_A}', 'Eq. 6.21', 'Recycle on the adiabatic line.'],
      ['F_{A0}C\'_p\\frac{dT}{dV} = (-\\Delta H_r)(-r_A) - Ua\\,(T-T_c)', '', 'Wall-cooled tube.'],
      ['\\int_{X_{i-1}}^{X_i} \\frac{\\partial}{\\partial T}\\left(\\frac{1}{-r_A}\\right) dX_A = 0,\\qquad (-r_A)_{\\text{out},i} = (-r_A)_{\\text{in},i+1}', 'Ch. 19', 'Optimal adiabatic stages with interstage cooling.']
    ],
    method: ['The locus of maximum rates is found at every conversion by a bracketed scan and Brent minimization in T within the allowable range; the optimum progression integrates 1/(−r_A)max by adaptive Gauss–Kronrod quadrature.', 'The best adiabatic feed temperature minimizes V/F_A0 over the feasible range (the line must reach X_Af short of equilibrium and below the upper limit). Recycle is optimized on that line, or jointly with the feed temperature.', 'N mixed reactors: intermediate conversions by Nelder–Mead, each reactor at the temperature of maximum rate for its outlet. Staged adiabatic beds: the two classical optimality conditions, solved by shooting on the first inlet temperature.', 'Heat duties are per mole of A fed, converted with F_A0. Coolers are sized as sensible heat on C′_p; the reactor itself is adiabatic unless stated.'],
    run(v) {
      const m = XTI.model(v), kpi = [], warn = [], info = [];
      const Xf = +v.Xf, FA0 = +v.FA0, tsec = m.tsec;
      const Tlo = m.Tlo, Thi = m.Thi, Tfeed = C2K(+v.Tfeed), Tprod = C2K(+v.Tprod);
      const kW = qJ => qJ * FA0 / tsec / 1000;
      const ru = 'mol/(L·' + (v._tu || 'min') + ')';
      const v0 = FA0 / m.CA0f(Tfeed);
      const Xcap = m.Xe(Tlo);
      if (!(Xf < Xcap)) throw new Error('X_Af = ' + N.fmt(Xf, 4) + ' cannot be reached: even at the lowest allowable temperature equilibrium stops at ' + N.fmt(Xcap, 4) + '.');
      const rows = [], extra = [], marks = [], lev = [], levMarks = [];
      const S = {};
      const xsF = N.linspace(0, Xf, 160);
      /* ---- optimum progression ---- */
      if (v.sOTP) {
        const o = TH.otp(m, Xf);
        S.otp = o;
        kpi.push({ k: 'VOTP', l: 'Optimum progression, plug flow', v: o.vf * FA0, u: 'L', main: true, s: 'τ = ' + N.fmt(o.vf * FA0 / v0, 4) + ' ' + (v._tu || 'min') });
        if (isFinite(o.XB)) kpi.push({ k: 'XB', l: 'Leaves T_max at X_A', v: o.XB, s: 'after ' + N.fmt(100 * o.fB, 3) + '% of the volume' });
        kpi.push({ k: 'ToutOTP', l: 'Exit temperature, optimum progression', v: K2C(o.T[o.T.length - 1]), u: '°C' });
        rows.push(['plug flow, optimum progression', o.vf * FA0, o.vf * FA0 / v0, K2C(o.T[0]), K2C(o.T[o.T.length - 1]), 'heat removed along the tube']);
        extra.push({ name: 'optimum progression', x: o.T.map(K2C), y: o.X, color: 3, w: 3 });
        lev.push({ name: 'optimum progression', x: o.X, y: o.inv, color: 3 }, { x: o.X, y: o.inv, type: 'area', fill: 0.15, color: 3, hover: false, legend: false, noRange: true });
        const prof = []; for (let f = 0; f <= 1.0001; f += 0.1) { const i = o.frac.findIndex(q => q >= f - 1e-9); prof.push([Math.round(f * 100) + '%', o.X[i], K2C(o.T[i])]); }
        S.otpTable = { title: 'Optimum progression through the reactor (Fig. E9.4c)', cols: ['fraction of volume', 'X_A', 'T, °C'], rows: prof };
      }
      /* ---- mixed flow on the locus ---- */
      if (v.sMFR) {
        const Tm = m.Topt(Xf), rM = m.rate(Xf, Tm), vf = Xf / rM;
        const T1 = m.T1for(Xf, Tm);
        const Cpp = m.Cp + (m.dHA(Tm) - m.dHA(T1)) / Math.max(Tm - T1, 1e-9);
        const qpre = m.Cp * (Tfeed - T1), qpost = Cpp * (Tm - Tprod);
        S.mfr = { Tm, rM, vf, T1, qpre, qpost };
        kpi.push({ k: 'VMFR', l: 'Mixed flow on the locus', v: vf * FA0, u: 'L', main: true, s: 'at ' + N.fmt(K2C(Tm), 4) + ' °C' });
        kpi.push({ k: 'rM', l: '−r_A in the mixed reactor', v: rM, u: ru });
        kpi.push({ k: 'TM', l: 'Mixed reactor temperature', v: K2C(Tm), u: '°C' }, { k: 'T1M', l: 'Adiabatic feed temperature for it', v: K2C(T1), u: '°C' });
        kpi.push({ k: 'Qpre', l: qpre >= 0 ? 'Precooler duty' : 'Preheater duty', v: Math.abs(kW(qpre)), u: 'kW', s: N.fmt(Math.abs(qpre) / TH.CAL, 4) + ' cal per mol A; ' + N.fmt(Math.abs(Tfeed - T1), 4) + ' K' });
        kpi.push({ k: 'Qpost', l: qpost >= 0 ? 'Product cooler duty' : 'Product heater duty', v: Math.abs(kW(qpost)), u: 'kW', s: N.fmt(Math.abs(qpost) / TH.CAL, 4) + ' cal per mol A; ' + N.fmt(Math.abs(Tm - Tprod), 4) + ' K' });
        rows.push(['mixed flow on the locus', vf * FA0, vf * FA0 / v0, K2C(T1), K2C(Tm), 'precool ' + N.fmt(Math.abs(kW(qpre)), 4) + ' kW, postcool ' + N.fmt(Math.abs(kW(qpost)), 4) + ' kW']);
        const TT = N.linspace(T1, Tm, 30);
        extra.push({ name: 'mixed flow: energy line', x: TT.map(K2C), y: TT.map(T => m.Xline(T, T1)), color: 2, dash: '5 4' });
        marks.push({ t: 'pt', x: K2C(Tm), y: Xf, color: 2, label: 'mixed flow' });
        levMarks.push({ t: 'rect', x0: 0, x1: Xf, y0: 0, y1: 1 / rM, color: 2, fill: 0.08, noRange: true, label: 'mixed' });
        if (T1 < Tlo - 1e-6) info.push('Mixed flow: the adiabatic feed temperature ' + N.fmt(K2C(T1), 4) + ' °C is below the lowest allowable temperature; feed warmer and remove heat inside the reactor instead (' + N.fmt(Math.abs(kW(qpre)), 4) + ' kW).');
      }
      /* ---- adiabatic plug flow ---- */
      let adT1 = NaN;
      if (v.sAD || v.sREC) {
        if (v.sAD && v.adMode === 'given') adT1 = C2K(+v.T1ad);
        else { const b = TH.adiabBest(m, Xf, 0, { Tcap: Thi }); if (b) adT1 = b.T1; }
      }
      if (v.sAD) {
        if (!isFinite(adT1)) warn.push('No adiabatic plug-flow line reaches X_Af within the allowable temperatures.');
        else {
          const Tof = X => m.Tline(X, adT1);
          const ok = m.rate(Xf, Tof(Xf)) > 0;
          const vf = ok ? TH.vf(m, Tof, 0, Xf) : NaN;
          const Tout = Tof(Xf);
          S.ad = { T1: adT1, vf, Tout, Tof };
          kpi.push({ k: 'VAD', l: 'Adiabatic plug flow', v: vf * FA0, u: 'L', main: true, s: v.adMode === 'given' ? 'feed at ' + N.fmt(+v.T1ad, 4) + ' °C' : 'best feed temperature' });
          kpi.push({ k: 'T1ad', l: 'Adiabatic feed temperature', v: K2C(adT1), u: '°C' }, { k: 'Toutad', l: 'Adiabatic exit temperature', v: K2C(Tout), u: '°C' });
          if (!ok) warn.push('The adiabatic line from ' + N.fmt(K2C(adT1), 4) + ' °C meets equilibrium before X_Af.');
          if (m.exo && Tout > Thi + 1e-6) warn.push('The adiabatic exit temperature exceeds the upper limit.');
          rows.push(['adiabatic plug flow', vf * FA0, vf * FA0 / v0, K2C(adT1), K2C(Tout), 'feed heat exchange ' + N.fmt(kW(m.Cp * (Tfeed - adT1)), 4) + ' kW removed']);
          const xs = xsF.filter(X => X <= Xf);
          extra.push({ name: 'adiabatic plug flow', x: xs.map(X => K2C(Tof(X))), y: xs, color: 6, w: 2 });
          lev.push({ name: 'adiabatic', x: xs, y: xs.map(X => 1 / m.rate(X, Tof(X))), color: 6 });
        }
      }
      /* ---- recycle ---- */
      if (v.sREC) {
        let T1 = adT1, rc = null;
        if (v.recMode === 'both') {
          const rg = TH.adiabRange(m, Xf, 0, { Tcap: Thi });
          if (rg) {
            const f = t => Math.log(TH.recycleBest(m, Xf, t).vf);
            const sc = N.fminScan(f, rg.lo, rg.hi * (1 - 1e-9), 18);
            T1 = sc.x;
          }
        }
        if (isFinite(T1)) {
          if (v.recMode === 'given') {
            const R = +v.Rgiven, Tof = X => m.Tline(X, T1), Xin = R / (1 + R) * Xf;
            rc = { R, vf: (1 + R) * TH.vf(m, Tof, Xin, Xf), Xin, Tin: Tof(Xin) };
          } else rc = TH.recycleBest(m, Xf, T1);
          S.rec = Object.assign({ T1 }, rc);
          kpi.push({ k: 'VREC', l: 'Adiabatic plug flow with recycle', v: rc.vf * FA0, u: 'L', main: true, s: 'R = ' + N.fmt(rc.R, 3) + ', fresh feed at ' + N.fmt(K2C(T1), 4) + ' °C' });
          kpi.push({ k: 'Rrec', l: 'Recycle ratio', v: rc.R });
          rows.push(['adiabatic plug flow, recycle R = ' + N.fmt(rc.R, 3), rc.vf * FA0, rc.vf * FA0 / v0, K2C(rc.Tin), K2C(m.Tline(Xf, T1)), 'reactor inlet after mixing at X_A = ' + N.fmt(rc.Xin, 3)]);
          const Tof = X => m.Tline(X, T1);
          const xs = N.linspace(rc.Xin, Xf, 60);
          extra.push({ name: 'recycle reactor', x: xs.map(X => K2C(Tof(X))), y: xs, color: 5, w: 4 });
          levMarks.push({ t: 'rect', x0: 0, x1: Xf, y0: 0, y1: rc.vf / Xf, color: 5, fill: 0.1, noRange: true, label: 'recycle' });
          if (!v.sAD) { const xa = xsF; lev.push({ name: 'adiabatic line used for recycle', x: xa, y: xa.map(X => 1 / m.rate(X, Tof(X))), color: 5, dash: '4 3' }); }
        }
      }
      /* ---- isothermal ---- */
      if (v.sISO) {
        const T = C2K(+v.Tiso);
        if (!(Xf < m.Xe(T))) warn.push('At ' + N.fmt(+v.Tiso, 4) + ' °C equilibrium stops at X = ' + N.fmt(m.Xe(T), 4) + ': no isothermal reactor reaches X_Af.');
        else {
          const vp = TH.vf(m, () => T, 0, Xf), vm = Xf / m.rate(Xf, T);
          S.iso = { T, vp, vm };
          kpi.push({ k: 'Viso', l: 'Isothermal plug flow at ' + N.fmt(+v.Tiso, 4) + ' °C', v: vp * FA0, u: 'L', main: true }, { k: 'VisoM', l: 'Isothermal mixed flow at ' + N.fmt(+v.Tiso, 4) + ' °C', v: vm * FA0, u: 'L' });
          kpi.push({ k: 'Qiso', l: 'Heat to remove, isothermal', v: kW(-m.dHA(T) * Xf + m.Cp * (T - Tfeed)), u: 'kW', s: 'feed brought to ' + N.fmt(+v.Tiso, 4) + ' °C inside' });
          rows.push(['isothermal plug flow', vp * FA0, vp * FA0 / v0, +v.Tiso, +v.Tiso, ''], ['isothermal mixed flow', vm * FA0, vm * FA0 / v0, +v.Tiso, +v.Tiso, '']);
          extra.push({ name: 'isothermal', x: [+v.Tiso, +v.Tiso], y: [0, Xf], color: 7, w: 2 });
          lev.push({ name: 'isothermal', x: xsF, y: xsF.map(X => 1 / m.rate(X, T)), color: 7 });
        }
      }
      /* ---- cooled tube ---- */
      if (v.sCOOL) {
        const cr = TH.cooled(m, Xf, C2K(+v.T1c), +v.Ua, C2K(+v.Tc), { vmax: 1e4 * (S.otp ? S.otp.vf : Xf / m.rate(0, C2K(+v.T1c))), Trun: Thi + 400 });
        S.cool = cr;
        if (cr.ok) {
          kpi.push({ k: 'Vcool', l: 'Wall-cooled plug flow', v: cr.vf * FA0, u: 'L', main: true, s: 'Ua = ' + N.fmt(+v.Ua, 3) + ' W/(L·K), coolant ' + N.fmt(+v.Tc, 4) + ' °C' });
          rows.push(['wall-cooled plug flow', cr.vf * FA0, cr.vf * FA0 / v0, +v.T1c, K2C(cr.Tend), 'hot spot ' + N.fmt(K2C(cr.Thot), 4) + ' °C']);
        } else warn.push('The cooled tube stalls at X = ' + N.fmt(cr.Xend, 4) + ' (' + N.fmt(K2C(cr.Tend), 4) + ' °C): it is quenched or pinned at equilibrium. Feed hotter, cool less, or move the coolant temperature.');
        kpi.push({ k: 'Thot', l: 'Hot spot, cooled tube', v: K2C(cr.Thot), u: '°C', state: cr.Thot > Thi ? 'bad' : undefined, s: 'at ' + N.fmt(100 * cr.vHot / (cr.ok ? cr.vf : cr.sol.t[cr.sol.t.length - 1]), 3) + '% of the length' });
        if (cr.Thot > Thi + 1e-6) warn.push('The cooled tube overshoots the upper temperature limit.');
        const n = cr.sol.t.length;
        extra.push({ name: 'cooled tube', x: cr.sol.y.map(y => K2C(y[1])), y: cr.sol.y.map(y => y[0]), color: 8, w: 2 });
        if (n > 2) lev.push({ name: 'cooled tube', x: cr.sol.y.map(y => y[0]), y: cr.sol.y.map(y => 1 / m.rate(y[0], y[1])), color: 8 });
      }
      /* ---- mixed reactors in series ---- */
      if (v.sNM) {
        const Nn = Math.round(+v.Nm), sr = TH.mfrSeries(m, Xf, Nn);
        S.nm = sr;
        kpi.push({ k: 'VNM', l: Nn + ' mixed reactors at their best temperatures', v: sr.vf * FA0, u: 'L', main: true, s: 'sizes ' + sr.each.map(e => N.fmt(e * FA0, 3)).join(' + ') + ' L' });
        rows.push([Nn + ' mixed reactors, each on the locus', sr.vf * FA0, sr.vf * FA0 / v0, K2C(sr.T[0]), K2C(sr.T[Nn - 1]), 'outlets X = ' + sr.X.slice(1).map(x => N.fmt(x, 3)).join(', ')]);
        extra.push({ name: Nn + ' mixed reactors', x: sr.T.map(K2C), y: sr.X.slice(1), type: 'pts', color: 1, r: 4.5 });
        sr.X.slice(1).forEach((x, i) => levMarks.push({ t: 'rect', x0: sr.X[i], x1: x, y0: 0, y1: 1 / m.rate(x, sr.T[i]), color: 1, fill: 0.1, noRange: true }));
      }
      /* ---- staged adiabatic with interstage cooling ---- */
      if (v.sNA) {
        if (!(m.rev && m.exo)) warn.push('Optimal adiabatic staging with interstage cooling applies to reversible exothermic reactions.');
        else {
          const Nn = Math.round(+v.Na), sg = TH.stagesAd(m, Xf, Nn);
          if (!sg) warn.push('No set of ' + Nn + ' adiabatic stages satisfies the optimality conditions for this duty.');
          else {
            S.na = sg;
            kpi.push({ k: 'VNA', l: Nn + ' adiabatic stages, interstage cooling', v: sg.vf * FA0, u: 'L', main: true, s: 'inlets ' + sg.stages.map(q => N.fmt(K2C(q.Tin), 3)).join(', ') + ' °C' });
            const Tmx = Math.max(...sg.stages.map(q => q.Tout));
            if (Tmx > Thi + 1e-6) warn.push('Staged design: a bed outlet reaches ' + N.fmt(K2C(Tmx), 4) + ' °C, above the upper limit (the optimum ignores that limit).');
            rows.push([Nn + ' adiabatic stages', sg.vf * FA0, sg.vf * FA0 / v0, K2C(sg.stages[0].Tin), K2C(sg.stages[Nn - 1].Tout), 'stage outlets X = ' + sg.stages.map(q => N.fmt(q.Xout, 3)).join(', ')]);
            const xs = [], ys = [];
            sg.stages.forEach((q, i) => {
              const xx = N.linspace(q.Xin, q.Xout, 30);
              xx.forEach(x => { xs.push(K2C(q.Tof(x))); ys.push(x); });
              if (i < Nn - 1) { xs.push(K2C(sg.stages[i + 1].Tin)); ys.push(q.Xout); }
              const lx = N.linspace(q.Xin, q.Xout, 40);
              lev.push({ name: i ? undefined : 'adiabatic stages', legend: i === 0, x: lx, y: lx.map(x => 1 / m.rate(x, q.Tof(x))), color: 4 });
            });
            extra.push({ name: 'staged adiabatic', x: xs, y: ys, color: 4, w: 2 });
          }
        }
      }
      /* ---- plots ---- */
      const Ta = Math.min(Tlo, isFinite(adT1) ? adT1 - 5 : Tlo, S.mfr ? S.mfr.T1 - 5 : Tlo), Tb = Thi + 5;
      const chart = XTI.chartSpec(m, { Ta, Tb, extra, marks: marks.concat([{ t: 'h', y: Xf, color: 8, label: 'X_{Af}' }, { t: 'v', x: K2C(Thi), color: 8, label: 'T_max' }]), ru, title: 'Operating lines on the X–T chart' });
      chart.x.min = K2C(Ta); chart.x.max = K2C(Tb);
      const plots = [chart];
      const levY = lev.filter(s => s.type !== 'area').map(s => Math.max(...s.y.filter(isFinite)));
      plots.push({ title: 'Levenspiel plots along each path (Figs. E9.4b, E9.6, E9.7)', x: { label: 'X_A', min: 0, max: Xf * 1.02 }, y: { label: '1/(−r_A)', unit: 'L·' + (v._tu || 'min') + '/mol', min: 0, max: Math.min(Math.max(...levY, ...levMarks.map(q => q.y1)) * 1.08, 12 * Math.min(...levY.filter(isFinite).concat([1e9]))) }, series: lev, marks: levMarks, cap: 'Area under a curve × F_A0 is a plug-flow volume; rectangles are mixed reactors or the recycle reactor.' });
      /* profiles */
      const prof = [];
      const profX = [];
      if (S.otp) { prof.push({ name: 'optimum progression', x: S.otp.frac, y: S.otp.T.map(K2C), color: 3 }); profX.push({ name: 'optimum progression', x: S.otp.frac, y: S.otp.X, color: 3 }); }
      if (S.ad) { const xs = xsF; const inv = xs.map(X => 1 / m.rate(X, S.ad.Tof(X))); const cum = N.cumtrapz(xs, inv), tot = cum[cum.length - 1]; prof.push({ name: 'adiabatic', x: cum.map(c => c / tot), y: xs.map(X => K2C(S.ad.Tof(X))), color: 6 }); profX.push({ name: 'adiabatic', x: cum.map(c => c / tot), y: xs, color: 6 }); }
      if (S.cool && S.cool.ok) { const t = S.cool.sol.t, tot = t[t.length - 1]; prof.push({ name: 'cooled tube', x: t.map(q => q / tot), y: S.cool.sol.y.map(y => K2C(y[1])), color: 8 }); profX.push({ name: 'cooled tube', x: t.map(q => q / tot), y: S.cool.sol.y.map(y => y[0]), color: 8 }); }
      if (S.na) {
        let acc = 0; const tx = [], ty = [], tX = [];
        S.na.stages.forEach(q => { const xx = N.linspace(q.Xin, q.Xout, 40), inv = xx.map(x => 1 / m.rate(x, q.Tof(x))), cum = N.cumtrapz(xx, inv); xx.forEach((x, i) => { tx.push((acc + cum[i]) / S.na.vf); ty.push(K2C(q.Tof(x))); tX.push(x); }); acc += cum[cum.length - 1]; });
        prof.push({ name: 'staged adiabatic', x: tx, y: ty, color: 4 }); profX.push({ name: 'staged adiabatic', x: tx, y: tX, color: 4 });
      }
      if (prof.length) {
        plots.push({ title: 'Temperature through the reactor', x: { label: 'fraction of reactor volume', min: 0, max: 1 }, y: { label: 'T', unit: '°C' }, series: prof, marks: [{ t: 'h', y: K2C(Thi), color: 8, label: 'T_max' }] });
        plots.push({ title: 'Conversion through the reactor', x: { label: 'fraction of reactor volume', min: 0, max: 1 }, y: { label: 'X_A', min: 0, max: 1 }, series: profX });
      }
      const tables = [{ title: 'Schemes compared for the same duty', cols: ['scheme', 'V, L', 'τ, ' + (v._tu || 'min'), 'T in, °C', 'T out, °C', 'notes'], rows }];
      if (S.otpTable) tables.push(S.otpTable);
      if (S.nm) tables.push({ title: 'Mixed reactors in series', cols: ['reactor', 'X_A out', 'T, °C', 'V, L', 'heat to remove, kW'], rows: S.nm.T.map((T, i) => { const Tin = i ? S.nm.T[i - 1] : Tfeed; const q = m.Cp * (T - Tin) + m.dHA(T) * (S.nm.X[i + 1] - S.nm.X[i]); return ['R-' + (i + 1), S.nm.X[i + 1], K2C(T), S.nm.each[i] * FA0, -kW(q)]; }) });
      if (S.na) tables.push({ title: 'Adiabatic stages', cols: ['stage', 'X in', 'X out', 'T in, °C', 'T out, °C', 'V, L', 'cooling after, kW'], rows: S.na.stages.map((q, i) => [i + 1, q.Xin, q.Xout, K2C(q.Tin), K2C(q.Tout), q.vf * FA0, i < S.na.stages.length - 1 ? kW(m.Cp * (q.Tout - S.na.stages[i + 1].Tin)) : NaN]) });
      /* ---- flowsheet ---- */
      let diagram = '';
      const f = v.focus, fmtT = T => N.fmt(K2C(T), 3) + ' °C';
      if (f === 'mfr' && S.mfr) diagram = flowsheet([{ hx: N.fmt(Math.abs(kW(S.mfr.qpre)), 3) + ' kW', after: fmtT(S.mfr.T1) }, { tank: { X: Xf, tag: 'R-901', label: fmtT(S.mfr.Tm) }, after: 'X ' + N.fmt(Xf, 3) }, { hx: N.fmt(Math.abs(kW(S.mfr.qpost)), 3) + ' kW', after: fmtT(Tprod) }], { feed: fmtT(Tfeed), aria: 'precooler, mixed reactor, product cooler' });
      else if ((f === 'ad' && S.ad) || (f === 'rec' && S.rec)) {
        const T1 = f === 'ad' ? S.ad.T1 : S.rec.T1, Tof = X => m.Tline(X, T1), Xa = f === 'rec' ? S.rec.Xin : 0;
        const pr = N.linspace(Xa, Xf, 16);
        diagram = flowsheet([{ hx: N.fmt(Math.abs(kW(m.Cp * (Tfeed - T1))), 3) + ' kW', after: fmtT(T1) }, { tube: { prof: pr, tag: 'R-902', label: 'adiabatic' }, after: fmtT(Tof(Xf)) }, { hx: N.fmt(Math.abs(kW(m.Cp * (Tof(Xf) - Tprod))), 3) + ' kW', after: fmtT(Tprod) }], { feed: fmtT(Tfeed), recycle: f === 'rec' ? [1, 1] : null, recLabel: 'R = ' + (S.rec ? N.fmt(S.rec.R, 3) : ''), aria: 'adiabatic tube' });
      } else if (f === 'otp' && S.otp) diagram = flowsheet([{ hx: '', after: fmtT(S.otp.T[0]) }, { tube: { prof: S.otp.X.filter((x, i) => i % 10 === 0), tag: 'R-903', label: 'programmed cooling', jacket: true, w: 160 }, after: fmtT(S.otp.T[S.otp.T.length - 1]) }], { feed: fmtT(Tfeed), aria: 'tube with programmed cooling' });
      else if (f === 'cool' && S.cool) diagram = flowsheet([{ hx: '', after: N.fmt(+v.T1c, 3) + ' °C' }, { tube: { prof: S.cool.sol.y.filter((y, i) => i % 4 === 0).map(y => y[0]), tag: 'R-904', label: 'coolant ' + N.fmt(+v.Tc, 3) + ' °C', jacket: true, w: 160 }, after: fmtT(S.cool.Tend) }], { feed: fmtT(Tfeed), aria: 'wall-cooled tube' });
      else if (f === 'nm' && S.nm) diagram = flowsheet(S.nm.T.map((T, i) => ({ tank: { X: S.nm.X[i + 1], tag: 'R-90' + (i + 1), label: fmtT(T), jacket: true }, after: 'X ' + N.fmt(S.nm.X[i + 1], 3), gap: 44 })), { feed: fmtT(Tfeed), aria: 'mixed reactors in series' });
      else if (f === 'na' && S.na) { const it = []; S.na.stages.forEach((q, i) => { it.push({ hx: '', after: fmtT(q.Tin) }); it.push({ tube: { prof: N.linspace(q.Xin, q.Xout, 8), tag: 'B-' + (i + 1), w: 70 }, after: fmtT(q.Tout), gap: 52 }); }); diagram = flowsheet(it, { feed: fmtT(Tfeed), aria: 'adiabatic stages with interstage coolers' }); }
      else if (f === 'iso' && S.iso) diagram = flowsheet([{ hx: '', after: N.fmt(+v.Tiso, 3) + ' °C' }, { tube: { prof: N.linspace(0, Xf, 12), tag: 'R-905', label: 'held at ' + N.fmt(+v.Tiso, 3) + ' °C', jacket: true, w: 150 }, after: 'X ' + N.fmt(Xf, 3) }], { feed: fmtT(Tfeed), aria: 'isothermal tube' });
      if (!diagram) info.push('Turn on the scheme chosen under “Flowsheet and profiles for” to see its flowsheet.');
      return { kpi, warn, info, plots, tables, diagram };
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
