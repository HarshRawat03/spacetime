/* Process schematics in the spirit of a P&ID: vessels, tubes, beds, towers, pipes, tags. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N;
  const D = ST.D = {};
  let uid = 0;
  const f1 = v => (Math.round(v * 10) / 10).toString();
  const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;');

  /* fluid colour for a conversion 0..1 */
  D.mix = p => 'color-mix(in oklab, var(--prod) ' + Math.round(N.clamp(p, 0, 1) * 100) + '%, var(--feed))';

  D.svg = function (w, h, body, label, scale = 1.25) {
    const id = 'd' + (++uid);
    return '<svg class="pid" width="' + f1(w * scale) + '" height="' + f1(h * scale) + '" viewBox="0 0 ' + f1(w) + ' ' + f1(h) + '" role="img" aria-label="' + esc(label || 'process schematic') + '" data-id="' + id + '">' +
      '<defs><marker id="' + id + 'a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 1L10 5L0 9z" class="ah"/></marker></defs>' +
      body.replace(/MARK/g, id + 'a') + '</svg>';
  };
  D.pipe = function (pts, o = {}) {
    const d = pts.map((p, i) => (i ? 'L' : 'M') + f1(p[0]) + ' ' + f1(p[1])).join('');
    let s = '<path class="pipe' + (o.cls ? ' ' + o.cls : '') + '" d="' + d + '"' + (o.arrow === false ? '' : ' marker-end="url(#MARK)"') + '/>';
    if (o.flow !== false) s += '<path class="flow' + (o.cls ? ' ' + o.cls : '') + '" d="' + d + '"/>';
    if (o.label) {
      const p = o.at !== undefined ? pts[o.at] : pts[0];
      s += '<text class="pid-l" x="' + f1(p[0] + (o.dx || 0)) + '" y="' + f1(p[1] + (o.dy !== undefined ? o.dy : -6)) + '"' + (o.anchor ? ' text-anchor="' + o.anchor + '"' : '') + '>' + ST.tlSvg(o.label, 10.5) + '</text>';
    }
    return s;
  };
  D.text = (x, y, s, o = {}) => '<text class="' + (o.cls || 'pid-l') + '" x="' + f1(x) + '" y="' + f1(y) + '" text-anchor="' + (o.anchor || 'middle') + '">' + ST.tlSvg(s, o.fs || 10.5) + '</text>';
  D.tag = function (x, y, s) {
    const w = 7 + String(s).length * 6.4;
    return '<g class="pid-tag"><rect x="' + f1(x - w / 2) + '" y="' + f1(y - 8) + '" width="' + f1(w) + '" height="15"/><text x="' + f1(x) + '" y="' + f1(y + 3.5) + '" text-anchor="middle">' + esc(s) + '</text></g>';
  };
  /* instrument bubble with a readout box */
  D.instr = function (x, y, code, val) {
    let s = '<g class="pid-ins"><circle cx="' + f1(x) + '" cy="' + f1(y) + '" r="11"/><text x="' + f1(x) + '" y="' + f1(y + 3.5) + '" text-anchor="middle">' + esc(code) + '</text>';
    if (val !== undefined && val !== null) {
      const w = 8 + String(val).length * 6.2;
      s += '<rect x="' + f1(x + 13) + '" y="' + f1(y - 8) + '" width="' + f1(w) + '" height="16"/><text class="v" x="' + f1(x + 17) + '" y="' + f1(y + 4) + '">' + esc(val) + '</text>';
    }
    return s + '</g>';
  };
  D.node = (x, y) => '<circle class="pid-node" cx="' + f1(x) + '" cy="' + f1(y) + '" r="3"/>';

  /* stirred vessel; x,y top-left; fill p = conversion (or null for empty) */
  D.tank = function (x, y, w, h, o = {}) {
    const r = Math.min(w * 0.22, 12), lv = y + h * (o.level !== undefined ? 1 - o.level : 0.22);
    let s = '<g class="pid-v">';
    const shell = 'M' + f1(x) + ' ' + f1(y + r) + 'Q' + f1(x) + ' ' + f1(y) + ' ' + f1(x + r) + ' ' + f1(y) + 'H' + f1(x + w - r) + 'Q' + f1(x + w) + ' ' + f1(y) + ' ' + f1(x + w) + ' ' + f1(y + r) +
      'V' + f1(y + h - r) + 'Q' + f1(x + w) + ' ' + f1(y + h) + ' ' + f1(x + w - r) + ' ' + f1(y + h) + 'H' + f1(x + r) + 'Q' + f1(x) + ' ' + f1(y + h) + ' ' + f1(x) + ' ' + f1(y + h - r) + 'Z';
    if (o.fill !== null && o.fill !== undefined) {
      const id = 'tk' + (++uid);
      s += '<clipPath id="' + id + '"><path d="' + shell + '"/></clipPath><rect class="liq" clip-path="url(#' + id + ')" x="' + f1(x) + '" y="' + f1(lv) + '" width="' + f1(w) + '" height="' + f1(y + h - lv) + '" style="fill:' + D.mix(o.fill) + '"/>';
      s += '<line class="lvl" x1="' + f1(x + 3) + '" x2="' + f1(x + w - 3) + '" y1="' + f1(lv) + '" y2="' + f1(lv) + '"/>';
    }
    s += '<path class="shell" d="' + shell + '"/>';
    if (o.stir !== false) {
      const cx = x + w / 2;
      s += '<line class="shaft" x1="' + f1(cx) + '" x2="' + f1(cx) + '" y1="' + f1(y - 10) + '" y2="' + f1(y + h * 0.72) + '"/>' +
        '<g class="imp"><path d="M' + f1(cx - w * 0.26) + ' ' + f1(y + h * 0.72) + 'h' + f1(w * 0.52) + '"/><path d="M' + f1(cx - w * 0.26) + ' ' + f1(y + h * 0.66) + 'v' + f1(h * 0.12) + 'M' + f1(cx + w * 0.26) + ' ' + f1(y + h * 0.66) + 'v' + f1(h * 0.12) + '"/></g>' +
        '<circle class="motor" cx="' + f1(cx) + '" cy="' + f1(y - 16) + '" r="7"/><text class="mt" x="' + f1(cx) + '" y="' + f1(y - 12.6) + '" text-anchor="middle">M</text>';
    }
    if (o.jacket) s += '<path class="jacket" d="M' + f1(x - 5) + ' ' + f1(y + h * 0.3) + 'V' + f1(y + h + 5) + 'H' + f1(x + w + 5) + 'V' + f1(y + h * 0.3) + '"/>';
    if (o.tag) s += D.tag(x + w / 2, y + h + 16, o.tag);
    if (o.label) s += D.text(x + w / 2, y + h + (o.tag ? 34 : 16), o.label);
    return s + '</g>';
  };

  /* tubular reactor; prof = array of conversions along the tube */
  D.tube = function (x, y, w, h, o = {}) {
    const id = 'g' + (++uid);
    let s = '<g class="pid-v">';
    const prof = o.prof && o.prof.length ? o.prof : null;
    if (prof) {
      s += '<linearGradient id="' + id + '" x1="0" x2="1" y1="0" y2="0">';
      const n = Math.min(prof.length, 24);
      for (let i = 0; i < n; i++) {
        const k = Math.round(i * (prof.length - 1) / Math.max(1, n - 1));
        s += '<stop offset="' + (i / Math.max(1, n - 1)).toFixed(3) + '" style="stop-color:' + D.mix(prof[k]) + '"/>';
      }
      s += '</linearGradient><rect class="liq" x="' + f1(x) + '" y="' + f1(y) + '" width="' + f1(w) + '" height="' + f1(h) + '" rx="' + f1(h / 2) + '" style="fill:url(#' + id + ')"/>';
    }
    s += '<rect class="shell" x="' + f1(x) + '" y="' + f1(y) + '" width="' + f1(w) + '" height="' + f1(h) + '" rx="' + f1(h / 2) + '"/>';
    const k = Math.max(2, Math.round(w / 26));
    for (let i = 1; i < k; i++) s += '<line class="slice" x1="' + f1(x + i * w / k) + '" x2="' + f1(x + i * w / k) + '" y1="' + f1(y + 3) + '" y2="' + f1(y + h - 3) + '"/>';
    if (o.packed) s += '<rect class="pack" x="' + f1(x + 6) + '" y="' + f1(y + 3) + '" width="' + f1(w - 12) + '" height="' + f1(h - 6) + '"/>';
    if (o.tag) s += D.tag(x + w / 2, y + h + 16, o.tag);
    if (o.label) s += D.text(x + w / 2, y + h + (o.tag ? 34 : 16), o.label);
    return s + '</g>';
  };

  /* vertical column: kind 'packed' | 'fluid' | 'tower' | 'bubble' | 'trickle' */
  D.column = function (x, y, w, h, o = {}) {
    const r = w / 2.4;
    let s = '<g class="pid-v">';
    const shell = 'M' + f1(x) + ' ' + f1(y + r) + 'A' + f1(w / 2) + ' ' + f1(r) + ' 0 0 1 ' + f1(x + w) + ' ' + f1(y + r) + 'V' + f1(y + h - r) + 'A' + f1(w / 2) + ' ' + f1(r) + ' 0 0 1 ' + f1(x) + ' ' + f1(y + h - r) + 'Z';
    const id = 'c' + (++uid);
    s += '<clipPath id="' + id + '"><path d="' + shell + '"/></clipPath>';
    const kind = o.kind || 'packed';
    const top = y + r, bot = y + h - r;
    if (o.prof && o.prof.length) {
      const gid = 'cg' + (++uid);
      s += '<linearGradient id="' + gid + '" x1="0" x2="0" y1="1" y2="0">';
      const n = Math.min(o.prof.length, 20);
      for (let i = 0; i < n; i++) {
        const k = Math.round(i * (o.prof.length - 1) / Math.max(1, n - 1));
        s += '<stop offset="' + (i / Math.max(1, n - 1)).toFixed(3) + '" style="stop-color:' + D.mix(o.prof[k]) + '"/>';
      }
      s += '</linearGradient><rect class="liq" clip-path="url(#' + id + ')" x="' + f1(x) + '" y="' + f1(y) + '" width="' + f1(w) + '" height="' + f1(h) + '" style="fill:url(#' + gid + ')"/>';
    }
    if (kind === 'packed' || kind === 'trickle' || kind === 'tower') {
      s += '<g clip-path="url(#' + id + ')"><path class="hatch" d="';
      const bh = (bot - top) * (o.bed || 0.78), b0 = bot - bh - (kind === 'packed' ? 0 : 6);
      for (let i = -h; i < w + h; i += 9) s += 'M' + f1(x + i) + ' ' + f1(b0 + bh) + 'L' + f1(x + i + bh) + ' ' + f1(b0);
      s += '"/></g><line class="support" x1="' + f1(x) + '" x2="' + f1(x + w) + '" y1="' + f1(b0 + bh) + '" y2="' + f1(b0 + bh) + '"/><line class="support" x1="' + f1(x) + '" x2="' + f1(x + w) + '" y1="' + f1(b0) + '" y2="' + f1(b0) + '"/>';
    }
    if (kind === 'fluid' || kind === 'bubble') {
      const bedTop = kind === 'fluid' ? bot - (bot - top) * (o.bed || 0.62) : top + 8;
      s += '<g clip-path="url(#' + id + ')">';
      if (kind === 'fluid') s += '<rect class="dense" x="' + f1(x) + '" y="' + f1(bedTop) + '" width="' + f1(w) + '" height="' + f1(bot - bedTop) + '"/>';
      let seed = 7;
      const rnd = () => { seed = (seed * 9301 + 49297) % 233280; return seed / 233280; };
      const nb = Math.round((bot - bedTop) * w / 260);
      for (let i = 0; i < nb; i++) {
        const by = bot - rnd() * (bot - bedTop - 6), br = 2.5 + 4.5 * (1 - (by - bedTop) / (bot - bedTop + 1));
        s += '<circle class="bub" cx="' + f1(x + 5 + rnd() * (w - 10)) + '" cy="' + f1(by) + '" r="' + f1(Math.max(2, br)) + '"/>';
      }
      s += '</g><line class="support" x1="' + f1(x) + '" x2="' + f1(x + w) + '" y1="' + f1(bot) + '" y2="' + f1(bot) + '" stroke-dasharray="3 2"/>';
    }
    s += '<path class="shell" d="' + shell + '"/>';
    if (o.tag) s += D.tag(x + w / 2, y + h + 16, o.tag);
    if (o.label) s += D.text(x + w / 2, y + h + (o.tag ? 34 : 16), o.label);
    return s + '</g>';
  };

  D.hx = function (x, y, r = 12, o = {}) {
    return '<g class="pid-v"><circle class="shell" cx="' + f1(x) + '" cy="' + f1(y) + '" r="' + f1(r) + '"/><path class="zig" d="M' + f1(x - r) + ' ' + f1(y) + 'l' + f1(r * 0.4) + ' ' + f1(-r * 0.5) + 'l' + f1(r * 0.4) + ' ' + f1(r) + 'l' + f1(r * 0.4) + ' ' + f1(-r) + 'l' + f1(r * 0.4) + ' ' + f1(r) + 'l' + f1(r * 0.4) + ' ' + f1(-r * 0.5) + '"/>' +
      (o.label ? D.text(x, y + r + 13, o.label) : '') + '</g>';
  };

  /* left-to-right train of reactors with stream readouts.
     units: [{type:'mfr'|'pfr', tag, X (outlet conversion), prof, label}] */
  D.train = function (units, o = {}) {
    const n = units.length;
    const tw = o.tw || 54, th = o.th || 62, pl = o.pl || 92, ph = o.ph || 22, gap = o.gap || 46;
    let x = 30, s = '';
    const yMid = 52;
    const w0 = 30 + units.reduce((t, u) => t + (u.type === 'pfr' ? pl : tw) + gap, 0) + 20;
    s += D.pipe([[4, yMid], [x - 2, yMid]], { label: o.feedLabel || 'feed', dy: -7 });
    units.forEach((u, i) => {
      const w = u.type === 'pfr' ? pl : tw;
      if (u.type === 'pfr') s += D.tube(x, yMid - ph / 2, w, ph, { prof: u.prof, tag: u.tag, label: u.label });
      else s += D.tank(x, yMid - th * 0.46, w, th, { fill: u.X, tag: u.tag, label: u.label });
      const xo = x + w;
      x = xo + gap;
      s += D.pipe([[xo + 1, yMid], [x - 2, yMid]], { arrow: i === n - 1 ? true : true });
      if (u.X !== undefined && o.readout !== false) s += D.text(xo + gap / 2, yMid - 9, (o.rlab || 'X') + ' ' + N.fmt(u.X, 3), { cls: 'pid-r', fs: 10 });
    });
    return { svg: D.svg(Math.max(w0, x + 4), o.h || 132, s, o.aria || 'reactor train'), w: x };
  };

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
