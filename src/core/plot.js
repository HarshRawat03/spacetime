/* SVG charts: line/points/area/step/bars, log axes, hover readout, legend toggles, drag handle. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N;
  const P = ST.Plot = {};
  let uid = 0;

  const fin = v => typeof v === 'number' && isFinite(v);
  const col = c => (c === undefined ? 'var(--s1)' : typeof c === 'number' ? 'var(--s' + (((c - 1) % 8) + 1) + ')' : c);

  function niceStep(span, target) {
    const raw = span / Math.max(1, target), p = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / p;
    return (m < 1.5 ? 1 : m < 2.25 ? 2 : m < 3.5 ? 2.5 : m < 7.5 ? 5 : 10) * p;
  }
  function linTicks(lo, hi, target) {
    const st = niceStep(hi - lo, target), out = [];
    const a = Math.ceil(lo / st - 1e-9) * st;
    for (let v = a; v <= hi + st * 1e-9; v += st) out.push(Math.abs(v) < st * 1e-9 ? 0 : v);
    return { ticks: out, step: st };
  }
  function logTicks(lo, hi) {
    const a = Math.floor(Math.log10(lo)), b = Math.ceil(Math.log10(hi)), major = [], minor = [];
    const dec = b - a;
    for (let e = a; e <= b; e++) {
      const base = Math.pow(10, e);
      if (base >= lo * 0.999 && base <= hi * 1.001) major.push(base);
      for (let m = 2; m < 10; m++) {
        const v = m * base;
        if (v < lo * 0.999 || v > hi * 1.001) continue;
        if (dec <= 1 && (m === 2 || m === 5)) major.push(v);
        else if (dec <= 1) minor.push(v);
        else if (dec <= 2 && (m === 2 || m === 5)) major.push(v);
        else minor.push(v);
      }
    }
    major.sort((x, y) => x - y);
    if (major.length < 2) { const l = linTicks(lo, hi, 4); return { ticks: l.ticks.filter(v => v > 0), minor: [] }; }
    if (dec > 8) return { ticks: major.filter(v => Math.abs(Math.log10(v) - Math.round(Math.log10(v))) < 1e-9).filter((v, i) => i % Math.ceil(dec / 6) === 0), minor: [] };
    return { ticks: major, minor };
  }
  function tickLabel(v, step, log, mult) {
    if (log) {
      const e = Math.log10(v);
      if (Math.abs(e - Math.round(e)) < 1e-9 && (Math.abs(e) >= 4)) return '10' + N.sup(Math.round(e));
      return N.fmt(v, 3);
    }
    const x = v / mult;
    if (x === 0) return '0';
    const st = step / mult;
    const dec = Math.max(0, Math.min(8, -Math.floor(Math.log10(st) + 1e-9) + (Math.abs(st / Math.pow(10, Math.floor(Math.log10(st))) - 2.5) < 1e-6 ? 1 : 0)));
    return x.toFixed(dec).replace('-', '−');
  }

  function range(vals, ax, log) {
    let lo = Infinity, hi = -Infinity;
    vals.forEach(v => { if (fin(v) && (!log || v > 0)) { if (v < lo) lo = v; if (v > hi) hi = v; } });
    if (fin(ax.min) && (!log || ax.min > 0)) lo = ax.min;
    if (fin(ax.max)) hi = ax.max;
    if (!fin(lo) || !fin(hi)) { lo = log ? 0.1 : 0; hi = log ? 10 : 1; }
    if (hi === lo) { if (log) { lo /= 2; hi *= 2; } else { const d = Math.abs(lo) || 1; lo -= 0.5 * d; hi += 0.5 * d; } }
    if (log) {
      if (!fin(ax.min)) lo = Math.pow(10, Math.floor(Math.log10(lo) * 4 + 1e-9) / 4);
      if (!fin(ax.max)) hi = Math.pow(10, Math.ceil(Math.log10(hi) * 4 - 1e-9) / 4);
      return [lo, hi];
    }
    const pad = ax.tight ? 0 : 0.04 * (hi - lo);
    if (!fin(ax.min)) lo = lo >= 0 && lo - pad < 0 && ax.zero !== false ? 0 : lo - pad;
    if (!fin(ax.max)) hi = hi + pad;
    if (ax.zero && lo > 0) lo = 0;
    return [lo, hi];
  }

  /* layout + draw; returns {svg, geom} */
  P.build = function (spec, W, holds) {
    const id = 'p' + (++uid);
    const fs = 12;
    const H = spec.h || Math.round(N.clamp(W * (spec.aspect || 0.62), 220, spec.hmax || 440));
    let xa = spec.x || {};
    const ya = spec.y || {};
    if (xa.cats) xa = Object.assign({ min: -0.6, max: xa.cats.length - 0.4 }, xa);
    const ser = (spec.series || []).filter(s => s && s.x && s.y);
    const hold = (holds || []).filter(s => s && s.x);
    const vis = ser.filter(s => !s.hidden && !s.noRange);
    const xs = [], ys = [];
    vis.concat(hold).forEach(s => {
      for (let i = 0; i < s.x.length; i++) {
        if (fin(s.x[i]) && fin(s.y[i])) { xs.push(s.x[i]); ys.push(s.y[i]); if (s.base && fin(s.base[i])) ys.push(s.base[i]); }
      }
    });
    (spec.marks || []).forEach(m => {
      if (m.noRange) return;
      if (m.t === 'pt' || m.t === 'text') { xs.push(m.x); ys.push(m.y); }
      if (m.t === 'rect') { xs.push(m.x0, m.x1); ys.push(m.y0, m.y1); }
      if (m.t === 'v') xs.push(m.x);
      if (m.t === 'h') ys.push(m.y);
      if (m.t === 'seg' || m.t === 'arrow') { xs.push(m.x0, m.x1); ys.push(m.y0, m.y1); }
    });
    if (fin(ya.base0) || ya.zero) ys.push(0);
    const [x0, x1] = range(xs, xa, xa.log), [y0, y1] = range(ys, ya, ya.log);
    /* axis multipliers for big/small linear numbers */
    const mulOf = (lo, hi, log) => {
      if (log) return 1;
      const m = Math.max(Math.abs(lo), Math.abs(hi));
      if (m >= 1e5 || (m > 0 && m < 1e-3)) return Math.pow(10, Math.floor(Math.log10(m)));
      return 1;
    };
    const xm = mulOf(x0, x1, xa.log), ym = mulOf(y0, y1, ya.log);
    const yT = ya.log ? logTicks(y0, y1) : linTicks(y0, y1, Math.max(3, Math.round(H / 55)));
    const yl = yT.ticks.map(v => tickLabel(v, yT.step, ya.log, ym));
    const maxYL = Math.max(1, ...yl.map(s => s.length));
    const legend = spec.legend !== false && ser.filter(s => s.name && s.legend !== false).length + hold.length > 1;
    /* category labels on x: wrap long ones onto two lines */
    const catLines = xa.cats ? xa.cats.map(c => {
      const t = String(c);
      if (t.length <= 14 || t.indexOf(' ') < 0) return [t];
      const mid = t.length / 2;
      let best = -1;
      for (let i = 0; i < t.length; i++) if (t[i] === ' ' && (best < 0 || Math.abs(i - mid) < Math.abs(best - mid))) best = i;
      return [t.slice(0, best), t.slice(best + 1)];
    }) : null;
    const catExtra = catLines && catLines.some(l => l.length > 1) ? 13 : 0;
    const L = 16 + (ya.label ? 18 : 0) + maxYL * 6.6, Rr = spec.y2 ? 54 : 14, T = 10 + (spec.title ? 16 : 0), B = 34 + (xa.label ? 14 : 0) + catExtra;
    const pw = Math.max(40, W - L - Rr), ph = Math.max(40, H - T - B);
    const xT = xa.cats ? { ticks: xa.cats.map((c, i) => i), minor: [] } : xa.log ? logTicks(x0, x1) : linTicks(x0, x1, Math.max(3, Math.round(pw / 85)));
    if (xa.ticks === false) { xT.ticks = []; xT.minor = []; }
    const sx = xa.log ? v => L + (Math.log10(v) - Math.log10(x0)) / (Math.log10(x1) - Math.log10(x0)) * pw : v => L + (v - x0) / (x1 - x0) * pw;
    const sy = ya.log ? v => T + ph - (Math.log10(v) - Math.log10(y0)) / (Math.log10(y1) - Math.log10(y0)) * ph : v => T + ph - (v - y0) / (y1 - y0) * ph;
    const ix = xa.log ? px => Math.pow(10, Math.log10(x0) + (px - L) / pw * (Math.log10(x1) - Math.log10(x0))) : px => x0 + (px - L) / pw * (x1 - x0);
    const clipY = v => N.clamp(v, T - 2000, T + ph + 2000);
    let o = '<svg class="plot" xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' + (spec.title ? ST.tlPlain(spec.title) : 'chart') + '">';
    o += '<defs><clipPath id="' + id + 'c"><rect x="' + L + '" y="' + (T - 1) + '" width="' + pw + '" height="' + (ph + 2) + '"/></clipPath></defs>';
    if (spec.title) o += '<text class="p-title" x="' + L + '" y="' + (T - 6) + '">' + ST.tlSvg(spec.title, 12) + '</text>';
    o += '<rect class="p-bg" x="' + L + '" y="' + T + '" width="' + pw + '" height="' + ph + '"/>';
    /* grid */
    (yT.minor || []).forEach(v => { o += '<line class="p-gm" x1="' + L + '" x2="' + (L + pw) + '" y1="' + sy(v).toFixed(1) + '" y2="' + sy(v).toFixed(1) + '"/>'; });
    (xT.minor || []).forEach(v => { o += '<line class="p-gm" y1="' + T + '" y2="' + (T + ph) + '" x1="' + sx(v).toFixed(1) + '" x2="' + sx(v).toFixed(1) + '"/>'; });
    if (ya.ticks === false) { yT.ticks = []; yT.minor = []; }
    yT.ticks.forEach((v, i) => {
      const y = sy(v).toFixed(1);
      o += '<line class="p-g" x1="' + L + '" x2="' + (L + pw) + '" y1="' + y + '" y2="' + y + '"/>';
      o += '<text class="p-tk" x="' + (L - 6) + '" y="' + (+y + 4) + '" text-anchor="end">' + yl[i] + '</text>';
    });
    xT.ticks.forEach(v => {
      const x = sx(v).toFixed(1);
      if (catLines) {
        catLines[v].forEach((ln, j) => { o += '<text class="p-tk p-cat" x="' + x + '" y="' + (T + ph + 16 + 13 * j) + '" text-anchor="middle">' + ST.tlSvg(ln, 11) + '</text>'; });
        return;
      }
      o += '<line class="p-g" y1="' + T + '" y2="' + (T + ph) + '" x1="' + x + '" x2="' + x + '"/>';
      o += '<text class="p-tk" x="' + x + '" y="' + (T + ph + 16) + '" text-anchor="middle">' + tickLabel(v, xT.step, xa.log, xm) + '</text>';
    });
    o += '<rect class="p-frame" x="' + L + '" y="' + T + '" width="' + pw + '" height="' + ph + '"/>';
    const axLab = (a, m) => ST.tlSvg(a.label || '', fs) + (a.unit ? '<tspan class="p-unit">' + ', ' + esc(a.unit) + '</tspan>' : '') + (m !== 1 ? '<tspan class="p-unit"> (×10' + N.sup(Math.round(Math.log10(m))) + ')</tspan>' : '');
    if (xa.label) o += '<text class="p-ax" x="' + (L + pw / 2) + '" y="' + (H - 8) + '" text-anchor="middle">' + axLab(xa, xm) + '</text>';
    if (ya.label) o += '<text class="p-ax" transform="translate(13 ' + (T + ph / 2) + ') rotate(-90)" text-anchor="middle">' + axLab(ya, ym) + '</text>';

    /* data */
    o += '<g clip-path="url(#' + id + 'c)">';
    const path = (s, xsA, ysA) => {
      let d = '', pen = false;
      for (let i = 0; i < xsA.length; i++) {
        const xv = xsA[i], yv = ysA[i];
        if (!fin(xv) || !fin(yv) || (xa.log && xv <= 0) || (ya.log && yv <= 0)) { pen = false; continue; }
        const X = sx(xv).toFixed(2), Y = clipY(sy(yv)).toFixed(2);
        if (s.type === 'step' && pen) d += 'H' + X + 'V' + Y;
        else d += (pen ? 'L' : 'M') + X + ' ' + Y;
        pen = true;
      }
      return d;
    };
    (spec.marks || []).filter(m => m.t === 'rect' || m.t === 'band').forEach(m => {
      if (m.t === 'band') {
        const a = sx(m.x0), b = sx(m.x1);
        o += '<rect class="p-band" x="' + Math.min(a, b).toFixed(1) + '" y="' + T + '" width="' + Math.abs(b - a).toFixed(1) + '" height="' + ph + '" style="fill:' + col(m.color || 8) + '"/>';
        return;
      }
      const a = sx(m.x0), b = sx(m.x1), c = sy(m.y0), d = sy(m.y1);
      o += '<rect class="p-rect" x="' + Math.min(a, b).toFixed(1) + '" y="' + Math.min(c, d).toFixed(1) + '" width="' + Math.abs(b - a).toFixed(1) + '" height="' + Math.abs(d - c).toFixed(1) + '" style="fill:' + col(m.color) + ';fill-opacity:' + (m.fill !== undefined ? m.fill : 0.16) + ';stroke:' + col(m.color) + '"/>';
    });
    hold.forEach(s => { o += '<path class="p-hold" d="' + path(s, s.x, s.y) + '"' + (s.color ? ' style="stroke:' + col(s.color) + '"' : '') + '/>'; });
    ser.forEach((s, si) => {
      if (s.hidden) return;
      const c = col(s.color !== undefined ? s.color : si + 1);
      if (s.type === 'area' || s.fill) {
        const base = s.base || null;
        let d = '';
        const n = s.x.length;
        const pts = [];
        for (let i = 0; i < n; i++) if (fin(s.x[i]) && fin(s.y[i])) pts.push(i);
        if (pts.length > 1) {
          d = 'M' + sx(s.x[pts[0]]).toFixed(2) + ' ' + clipY(sy(base ? base[pts[0]] : (ya.log ? y0 : Math.max(y0, Math.min(0, y1))))).toFixed(2);
          pts.forEach(i => { d += 'L' + sx(s.x[i]).toFixed(2) + ' ' + clipY(sy(s.y[i])).toFixed(2); });
          for (let k = pts.length - 1; k >= 0; k--) { const i = pts[k]; d += 'L' + sx(s.x[i]).toFixed(2) + ' ' + clipY(sy(base ? base[i] : (ya.log ? y0 : Math.max(y0, Math.min(0, y1))))).toFixed(2); }
          o += '<path class="p-area" d="' + d + 'Z" style="fill:' + c + ';fill-opacity:' + (typeof s.fill === 'number' ? s.fill : 0.16) + '"/>';
        }
        if (s.type === 'area' && !s.line) return;
      }
      if (s.type === 'bar') {
        const w = s.bw || (pw / Math.max(8, s.x.length * 1.6));
        s.x.forEach((xv, i) => {
          if (!fin(xv) || !fin(s.y[i])) return;
          const yb = sy(ya.log ? y0 : Math.max(y0, 0)), yt = sy(s.y[i]);
          o += '<rect class="p-bar" x="' + (sx(xv) - w / 2).toFixed(1) + '" y="' + Math.min(yb, yt).toFixed(1) + '" width="' + w.toFixed(1) + '" height="' + Math.abs(yb - yt).toFixed(1) + '" style="fill:' + c + '"/>';
        });
        return;
      }
      if (s.type !== 'pts') {
        o += '<path class="p-line" d="' + path(s, s.x, s.y) + '" style="stroke:' + c + (s.w ? ';stroke-width:' + s.w : '') + (s.dash ? ';stroke-dasharray:' + s.dash : '') + '"/>';
      }
      if (s.type === 'pts' || s.marker) {
        s.x.forEach((xv, i) => {
          const yv = s.y[i];
          if (!fin(xv) || !fin(yv) || (xa.log && xv <= 0) || (ya.log && yv <= 0)) return;
          const X = sx(xv).toFixed(1), Y = sy(yv).toFixed(1), r = s.r || 3.6;
          const mk = s.marker || 'o';
          const st = 'stroke:' + c + ';fill:' + (s.open ? 'var(--surface)' : c);
          if (mk === 's') o += '<rect class="p-pt" x="' + (X - r) + '" y="' + (Y - r) + '" width="' + 2 * r + '" height="' + 2 * r + '" style="' + st + '"/>';
          else if (mk === '^') o += '<path class="p-pt" d="M' + X + ' ' + (Y - r * 1.2) + 'L' + (+X + r * 1.1) + ' ' + (+Y + r * 0.8) + 'L' + (X - r * 1.1) + ' ' + (+Y + r * 0.8) + 'Z" style="' + st + '"/>';
          else o += '<circle class="p-pt" cx="' + X + '" cy="' + Y + '" r="' + r + '" style="' + st + '"/>';
        });
      }
    });
    (spec.marks || []).forEach(m => {
      const c = col(m.color || 8);
      if (m.t === 'v' && fin(m.x)) {
        const X = sx(m.x).toFixed(1);
        o += '<line class="p-ref" x1="' + X + '" x2="' + X + '" y1="' + T + '" y2="' + (T + ph) + '" style="stroke:' + c + (m.dash === false ? ';stroke-dasharray:none' : '') + '"/>';
        if (m.label) { const right = +X > L + pw - 46; o += '<text class="p-lab" x="' + (right ? +X - 4 : +X + 4) + '" y="' + (T + 12 + (m.dy || 0)) + '"' + (right ? ' text-anchor="end"' : '') + ' style="fill:' + c + '">' + ST.tlSvg(m.label, 11) + '</text>'; }
      }
      if (m.t === 'h' && fin(m.y)) {
        const Y = sy(m.y).toFixed(1);
        o += '<line class="p-ref" x1="' + L + '" x2="' + (L + pw) + '" y1="' + Y + '" y2="' + Y + '" style="stroke:' + c + (m.dash === false ? ';stroke-dasharray:none' : '') + '"/>';
        if (m.label) o += '<text class="p-lab" x="' + (L + pw - 4) + '" y="' + (Y - 4) + '" text-anchor="end" style="fill:' + c + '">' + ST.tlSvg(m.label, 11) + '</text>';
      }
      if ((m.t === 'seg' || m.t === 'arrow') && fin(m.x0) && fin(m.y0) && fin(m.x1) && fin(m.y1)) {
        o += '<line class="p-seg" x1="' + sx(m.x0).toFixed(1) + '" y1="' + sy(m.y0).toFixed(1) + '" x2="' + sx(m.x1).toFixed(1) + '" y2="' + sy(m.y1).toFixed(1) + '" style="stroke:' + c + (m.w ? ';stroke-width:' + m.w : '') + (m.dash ? ';stroke-dasharray:' + m.dash : '') + '"/>';
      }
      if (m.t === 'pt' && fin(m.x) && fin(m.y)) {
        o += '<circle class="p-mk" cx="' + sx(m.x).toFixed(1) + '" cy="' + sy(m.y).toFixed(1) + '" r="' + (m.r || 4.5) + '" style="stroke:' + c + ';fill:' + (m.open ? 'var(--surface)' : c) + '"/>';
        if (m.label) o += '<text class="p-lab" x="' + (sx(m.x) + (m.dx !== undefined ? m.dx : 7)).toFixed(1) + '" y="' + (sy(m.y) + (m.dy !== undefined ? m.dy : -7)).toFixed(1) + '"' + (m.anchor ? ' text-anchor="' + m.anchor + '"' : '') + ' style="fill:' + c + '">' + ST.tlSvg(m.label, 11) + '</text>';
      }
      if (m.t === 'text' && fin(m.x) && fin(m.y)) {
        o += '<text class="p-lab" x="' + sx(m.x).toFixed(1) + '" y="' + sy(m.y).toFixed(1) + '" text-anchor="' + (m.anchor || 'start') + '" style="fill:' + c + '">' + ST.tlSvg(m.s, 11) + '</text>';
      }
      if (m.t === 'rect' && m.label) {
        o += '<text class="p-lab" x="' + ((sx(m.x0) + sx(m.x1)) / 2).toFixed(1) + '" y="' + (Math.min(sy(m.y0), sy(m.y1)) + 14).toFixed(1) + '" text-anchor="middle" style="fill:' + col(m.color) + '">' + ST.tlSvg(m.label, 11) + '</text>';
      }
    });
    o += '</g>';
    /* legend */
    let leg = '';
    if (legend) {
      leg = '<div class="p-leg">' + ser.map((s, si) => (s.name && s.legend !== false ? '<button type="button" class="p-key' + (s.hidden ? ' off' : '') + '" data-si="' + si + '"><i style="background:' + col(s.color !== undefined ? s.color : si + 1) + (s.type === 'pts' ? ';border-radius:50%;width:8px;height:8px' : '') + '"></i><span class="pk-t">' + ST.tl(s.name) + '</span></button>' : '')).join('') +
        (hold.length ? '<span class="p-key held"><i></i>held</span>' : '') + '</div>';
    }
    if (spec.drag && fin(spec.drag.x)) {
      const X = sx(spec.drag.x).toFixed(1);
      o += '<g class="p-drag"><line x1="' + X + '" x2="' + X + '" y1="' + T + '" y2="' + (T + ph) + '"/><rect x="' + (X - 7) + '" y="' + (T + ph - 18) + '" width="14" height="18" rx="2"/><path d="M' + (X - 2) + ' ' + (T + ph - 13) + 'v8M' + (+X + 2) + ' ' + (T + ph - 13) + 'v8"/></g>';
    }
    o += '<line class="p-cross" x1="0" x2="0" y1="' + T + '" y2="' + (T + ph) + '" visibility="hidden"/><g class="p-dots"></g>';
    o += '</svg>';
    return { svg: o, leg, geom: { L, T, pw, ph, W, H, sx, sy, ix, x0, x1, y0, y1 } };
  };
  const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;');

  function interpAt(s, xq) {
    const n = s.x.length;
    if (n < 2) return NaN;
    let asc = true;
    for (let i = 1; i < n; i++) if (fin(s.x[i]) && fin(s.x[i - 1]) && s.x[i] < s.x[i - 1]) { asc = false; break; }
    if (!asc) return NaN;
    if (xq < s.x[0] || xq > s.x[n - 1]) return NaN;
    if (s.type === 'step') { let i = 0; while (i < n - 1 && s.x[i + 1] <= xq) i++; return s.y[i]; }
    return N.interp(s.x, s.y, xq);
  }

  /* render into host element; keeps state for resize */
  P.render = function (host, spec, o = {}) {
    host.classList.add('plot-host');
    const draw = () => {
      const W = Math.max(260, Math.floor(host.clientWidth || 600));
      const b = P.build(spec, W, o.holds);
      host.innerHTML = '<div class="p-wrap">' + b.svg + '<div class="p-tip" hidden></div></div>' + b.leg;
      const svg = host.querySelector('svg'), tip = host.querySelector('.p-tip'), cross = svg.querySelector('.p-cross'), dots = svg.querySelector('.p-dots');
      const g = b.geom;
      host.querySelectorAll('.p-key[data-si]').forEach(k => k.addEventListener('click', () => {
        const s = spec.series[+k.dataset.si];
        s.hidden = !s.hidden;
        draw();
      }));
      const ser = (spec.series || []).filter(s => s && s.x && !s.hidden && s.hover !== false);
      const hide = () => { cross.setAttribute('visibility', 'hidden'); dots.innerHTML = ''; tip.hidden = true; };
      let dragging = false;
      const toLocal = ev => {
        const r = svg.getBoundingClientRect();
        return { px: (ev.clientX - r.left) * (g.W / r.width), py: (ev.clientY - r.top) * (g.H / r.height), r };
      };
      svg.addEventListener('pointermove', ev => {
        const { px, py, r } = toLocal(ev);
        if (dragging && spec.drag) {
          const xv = N.clamp(g.ix(px), spec.drag.min !== undefined ? spec.drag.min : g.x0, spec.drag.max !== undefined ? spec.drag.max : g.x1);
          spec.drag.onDrag(xv);
          return;
        }
        if (px < g.L || px > g.L + g.pw || py < g.T || py > g.T + g.ph) { hide(); return; }
        const xq = g.ix(px);
        let rows = [], dh = '';
        ser.forEach((s, si) => {
          const c = col(s.color !== undefined ? s.color : spec.series.indexOf(s) + 1);
          if (s.type === 'pts' || s.type === 'bar') {
            let best = -1, bd = 14;
            s.x.forEach((xv, i) => {
              if (!fin(xv) || !fin(s.y[i])) return;
              const d = Math.hypot(g.sx(xv) - px, g.sy(s.y[i]) - py);
              if (d < bd) { bd = d; best = i; }
            });
            if (best >= 0) {
              rows.push('<div><i style="background:' + c + '"></i>' + ST.tl(s.name || 'data') + ' <b>(' + N.fmt(s.x[best], 4) + ', ' + N.fmt(s.y[best], 4) + ')</b></div>');
            }
            return;
          }
          const yv = interpAt(s, xq);
          if (!fin(yv)) return;
          const Y = g.sy(yv);
          if (Y < g.T - 1 || Y > g.T + g.ph + 1) return;
          dh += '<circle cx="' + px.toFixed(1) + '" cy="' + Y.toFixed(1) + '" r="3.2" style="fill:' + c + '"/>';
          rows.push('<div><i style="background:' + c + '"></i>' + ST.tl(s.name || spec.y && spec.y.label || 'y') + ' <b>' + N.fmt(yv, 4) + '</b></div>');
        });
        cross.setAttribute('x1', px); cross.setAttribute('x2', px); cross.setAttribute('visibility', 'visible');
        dots.innerHTML = dh;
        if (!rows.length) { tip.hidden = true; return; }
        tip.innerHTML = '<div class="p-tx">' + ST.tl((spec.x && spec.x.label) || 'x') + ' = <b>' + N.fmt(xq, 4) + '</b></div>' + rows.join('');
        tip.hidden = false;
        const tw = tip.offsetWidth, scale = r.width / g.W;
        let left = px * scale + 14;
        if (left + tw > r.width - 4) left = px * scale - tw - 14;
        tip.style.left = Math.max(0, left) + 'px';
        tip.style.top = Math.max(0, g.T * scale + 4) + 'px';
      });
      svg.addEventListener('pointerleave', () => { if (!dragging) hide(); });
      if (spec.drag) {
        svg.style.touchAction = 'none';
        svg.addEventListener('pointerdown', ev => {
          const { px, py } = toLocal(ev);
          if (px < g.L - 10 || px > g.L + g.pw + 10 || py < g.T || py > g.T + g.ph + 4) return;
          dragging = true;
          svg.setPointerCapture(ev.pointerId);
          hide();
          const xv = N.clamp(g.ix(px), spec.drag.min !== undefined ? spec.drag.min : g.x0, spec.drag.max !== undefined ? spec.drag.max : g.x1);
          spec.drag.onDrag(xv);
        });
        const end = ev => { dragging = false; try { svg.releasePointerCapture(ev.pointerId); } catch (e) { } };
        svg.addEventListener('pointerup', end);
        svg.addEventListener('pointercancel', end);
      }
    };
    draw();
    if (!host._ro && typeof ResizeObserver !== 'undefined') {
      let lastW = host.clientWidth, raf = 0;
      host._ro = new ResizeObserver(() => {
        if (Math.abs(host.clientWidth - lastW) < 2) return;
        lastW = host.clientWidth;
        cancelAnimationFrame(raf);
        raf = requestAnimationFrame(() => host._redraw && host._redraw());
      });
      host._ro.observe(host);
    }
    host._redraw = draw;
    return { redraw: draw };
  };

  /* marching squares: Z[j][i] at (xs[i], ys[j]); returns [{level, x:[], y:[]}] with NaN breaks between segments */
  P.contour = function (xs, ys, Z, levels) {
    const out = [];
    levels.forEach(L => {
      const X = [], Y = [];
      for (let j = 0; j < ys.length - 1; j++) {
        for (let i = 0; i < xs.length - 1; i++) {
          const z = [Z[j][i], Z[j][i + 1], Z[j + 1][i + 1], Z[j + 1][i]];
          if (z.some(q => !isFinite(q))) continue;
          const px = [xs[i], xs[i + 1], xs[i + 1], xs[i]], py = [ys[j], ys[j], ys[j + 1], ys[j + 1]];
          const pts = [];
          for (let e = 0; e < 4; e++) {
            const a = e, b = (e + 1) % 4, za = z[a] - L, zb = z[b] - L;
            if ((za < 0) !== (zb < 0)) { const t = za / (za - zb); pts.push([px[a] + t * (px[b] - px[a]), py[a] + t * (py[b] - py[a])]); }
          }
          if (pts.length === 2) { X.push(pts[0][0], pts[1][0], NaN); Y.push(pts[0][1], pts[1][1], NaN); }
          else if (pts.length === 4) { X.push(pts[0][0], pts[1][0], NaN, pts[2][0], pts[3][0], NaN); Y.push(pts[0][1], pts[1][1], NaN, pts[2][1], pts[3][1], NaN); }
        }
      }
      out.push({ level: L, x: X, y: Y });
    });
    return out;
  };

  /* CSV of all series */
  P.csv = function (spec) {
    const out = [];
    (spec.series || []).forEach(s => {
      if (!s || !s.x) return;
      out.push('# ' + ST.tlPlain(s.name || 'series'));
      out.push(N.fmtCSV(ST.tlPlain((spec.x && spec.x.label) || 'x')) + ',' + N.fmtCSV(ST.tlPlain(s.name || (spec.y && spec.y.label) || 'y')));
      s.x.forEach((x, i) => { if (fin(x) && fin(s.y[i])) out.push(N.fmtCSV(x) + ',' + N.fmtCSV(s.y[i])); });
      out.push('');
    });
    return out.join('\n');
  };

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
