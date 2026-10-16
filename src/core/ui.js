/* Tool framework: registry, input panel, outputs, presets with book checks, sweeps, held curves. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const N = ST.N;
  ST.TOOLS = ST.TOOLS || [];
  ST.tool = spec => { ST.TOOLS.push(spec); return spec; };
  if (typeof document === 'undefined') { if (typeof module !== 'undefined') module.exports = ST; return; }

  const UI = ST.UI = {};
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));
  const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  UI.esc = esc;
  const store = {
    get(k) { try { const s = localStorage.getItem(k); return s ? JSON.parse(s) : null; } catch (e) { return null; } },
    set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { /* storage unavailable */ } },
    del(k) { try { localStorage.removeItem(k); } catch (e) { /* ignore */ } }
  };
  UI.store = store;

  UI.parseNum = function (s) {
    if (typeof s === 'number') return s;
    let t = String(s).trim().replace(/−/g, '-').replace(/,/g, '');
    t = t.replace(/\s*[×x*]\s*10\^?\s*\(?([-+]?\d+)\)?/, 'e$1');
    if (!/^[-+]?(\d+\.?\d*|\.\d+)(e[-+]?\d+)?$/i.test(t)) return NaN;
    return parseFloat(t);
  };
  const showNum = v => {
    if (!isFinite(v)) return String(v);
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-4 || a >= 1e7)) return v.toExponential(4).replace(/\.?0+e/, 'e');
    return String(+v.toPrecision(8));
  };

  const texCache = new Map();
  UI.texHTML = function (tex, display) {
    const key = (display ? 'D:' : 'I:') + tex;
    if (ST.TEX && ST.TEX[key]) return ST.TEX[key];
    if (texCache.has(key)) return texCache.get(key);
    if (G.katex) { try { const h = G.katex.renderToString(tex, { displayMode: !!display, throwOnError: false, strict: false }); texCache.set(key, h); return h; } catch (e) { /* fall through */ } }
    return '<code class="tex-src">' + esc(tex) + '</code>';
  };
  /* label with $..$ segments → HTML using light math markup */
  UI.lab = s => ST.tl(s || '');
  /* '§9.1' from 'Ch. 9.1'; nothing when the reference is just the whole chapter */
  UI.secRef = sec => { const r = String(sec || '').replace(/^Ch\.\s*/, '').trim(); return !r || /^\d+$/.test(r) ? '' : '§' + r; };

  UI.copy = function (text, btn) {
    const done = ok => {
      if (!btn) return;
      const old = btn.textContent;
      btn.textContent = ok ? 'Copied' : 'Select and copy';
      setTimeout(() => { btn.textContent = old; }, 1400);
    };
    try {
      navigator.clipboard.writeText(text).then(() => done(true), () => { fallback(); done(false); });
    } catch (e) { fallback(); done(false); }
    function fallback() {
      const ta = document.createElement('textarea');
      ta.value = text; ta.className = 'copy-buf';
      document.body.appendChild(ta); ta.select();
      try { document.execCommand('copy'); } catch (e) { /* ignore */ }
      setTimeout(() => ta.remove(), 4000);
    }
  };

  const unitSub = (u, v) => String((typeof u === 'function' ? u(v) : u) || '').replace(/\[C\]/g, v._cu || 'mol/L').replace(/\[t\]/g, v._tu || 'min').replace(/\[V\]/g, v._vu || 'L');

  /* ---------- inputs ---------- */
  function fieldHTML(f, v, tid) {
    const id = 'in-' + tid + '-' + f.k;
    const hint = f.hint ? '<small class="fh">' + UI.lab(f.hint) + '</small>' : '';
    const unit = f.unit ? '<span class="fu">' + ST.tu(unitSub(f.unit, v)) + '</span>' : '';
    const lab = '<span class="fl">' + UI.lab((typeof f.label === 'function' ? f.label(v) : f.label) || f.k) + '</span>';
    switch (f.t) {
      case 'select':
        return '<div class="f fsel" data-k="' + f.k + '"><label for="' + id + '">' + lab + unit + '</label><select id="' + id + '">' +
          f.opts.map(o => '<option value="' + esc(o[0]) + '"' + (String(v[f.k]) === String(o[0]) ? ' selected' : '') + '>' + esc(ST.tlPlain(o[1])) + '</option>').join('') + '</select>' + hint + '</div>';
      case 'seg':
        return '<div class="f fseg" data-k="' + f.k + '" role="radiogroup" aria-label="' + esc(ST.tlPlain(f.label || f.k)) + '"><div class="fl-row">' + lab + '</div><div class="seg">' +
          f.opts.map((o, i) => '<label><input type="radio" name="' + id + '" id="' + id + '-' + i + '" value="' + esc(o[0]) + '"' + (String(v[f.k]) === String(o[0]) ? ' checked' : '') + '><span>' + UI.lab(o[1]) + '</span></label>').join('') + '</div>' + hint + '</div>';
      case 'check':
        return '<div class="f fchk" data-k="' + f.k + '"><label class="chk" for="' + id + '"><input type="checkbox" id="' + id + '"' + (v[f.k] ? ' checked' : '') + '><span>' + UI.lab(f.label) + '</span></label>' + hint + '</div>';
      case 'text':
        return '<div class="f ftxt" data-k="' + f.k + '"><label for="' + id + '">' + lab + unit + '</label><input id="' + id + '" type="text" spellcheck="false" autocomplete="off" class="' + (f.mono ? 'mono' : '') + '" value="' + esc(v[f.k]) + '">' + hint + '<small class="ferr" hidden></small></div>';
      case 'area':
        return '<div class="f farea" data-k="' + f.k + '"><label for="' + id + '">' + lab + '</label><textarea id="' + id + '" rows="' + (f.rows || 5) + '" spellcheck="false" class="mono">' + esc(v[f.k]) + '</textarea>' + hint + '<small class="ferr" hidden></small></div>';
      case 'table': return tableHTML(f, v, id);
      default: {
        const val = v[f.k];
        let rng = '';
        if (f.min !== undefined && f.max !== undefined && f.slider !== false) {
          const lg = f.scale === 'log';
          const lo = lg ? Math.log10(f.min > 0 ? f.min : f.max * 1e-5) : f.min, hi = lg ? Math.log10(f.max) : f.max;
          const st = lg ? (hi - lo) / 400 : (f.step || (hi - lo) / 200);
          const sv = lg ? Math.log10(Math.max(val, Math.pow(10, lo))) : val;
          rng = '<input type="range" class="rng" tabindex="-1" aria-hidden="true" min="' + lo + '" max="' + hi + '" step="' + st + '" value="' + sv + '">';
        }
        return '<div class="f fnum" data-k="' + f.k + '"><label for="' + id + '">' + lab + unit + '</label><input id="' + id + '" class="num" type="text" inputmode="decimal" autocomplete="off" spellcheck="false" value="' + esc(showNum(val)) + '">' + rng + hint + '</div>';
      }
    }
  }
  function tableHTML(f, v, id) {
    const rows = v[f.k] || [];
    const cols = f.cols;
    let h = '<div class="f ftab" data-k="' + f.k + '"><div class="fl-row"><span class="fl">' + UI.lab((typeof f.label === 'function' ? f.label(v) : f.label) || 'Data') + '</span><span class="fu">' + rows.length + ' rows</span></div>' +
      '<div class="tab-wrap"><table class="dt"><thead><tr><th class="ri">#</th>' + cols.map((c, j) => '<th>' + UI.lab(c) + (f.units && f.units[j] ? '<span class="tu">' + ST.tu(unitSub(f.units[j], v)) + '</span>' : '') + '</th>').join('') + '<th></th></tr></thead><tbody>';
    rows.forEach((r, i) => {
      h += '<tr><td class="ri">' + (i + 1) + '</td>' + cols.map((c, j) => '<td><input id="' + id + '-' + i + '-' + j + '" data-r="' + i + '" data-c="' + j + '" type="text" inputmode="decimal" autocomplete="off" value="' + (r[j] === null || r[j] === undefined || Number.isNaN(r[j]) ? '' : esc(showNum(r[j]))) + '" aria-label="' + esc(ST.tlPlain(c)) + ' row ' + (i + 1) + '"></td>').join('') +
        '<td><button type="button" class="rm" data-r="' + i + '" aria-label="Remove row ' + (i + 1) + '">×</button></td></tr>';
    });
    h += '</tbody></table></div><div class="tab-act"><button type="button" class="btn-s add">Add row</button><button type="button" class="btn-s clr">Clear</button><span class="fh">Paste columns from a spreadsheet into any cell.</span></div></div>';
    return h;
  }

  UI.mountTool = function (host, spec) {
    const defaults = {};
    spec.inputs.forEach(f => { if (f.k) defaults[f.k] = JSON.parse(JSON.stringify(f.v)); });
    if (spec.units) { defaults._cu = spec.units.C || 'mol/L'; defaults._tu = spec.units.t || 'min'; }
    const saved = store.get('st:v:' + spec.id);
    const v = Object.assign({}, defaults, saved || {});
    Object.keys(v).forEach(k => { if (!(k in defaults)) delete v[k]; });
    const state = { v, holds: {}, preset: null, last: null, plots: [] };
    host.innerHTML = headHTML(spec) + '<div class="tool-body"><form class="inputs" autocomplete="off" novalidate></form><section class="outputs" aria-live="polite"></section></div>' + methodHTML(spec) + sweepHTML(spec);
    const form = $('.inputs', host), out = $('.outputs', host);

    const unitsBar = () => {
      if (!spec.units) return '';
      const cu = ['mol/L', 'mol/m³', 'kmol/m³', 'mol/cm³', 'mmol/L', 'g/L', 'kg/m³'];
      const tu = ['s', 'min', 'h', 'day'];
      return '<div class="units">' + (spec.units.C ? '<label>Concentration <select data-u="_cu">' + cu.map(u => '<option' + (v._cu === u ? ' selected' : '') + '>' + u + '</option>').join('') + '</select></label>' : '') +
        (spec.units.t ? '<label>Time <select data-u="_tu">' + tu.map(u => '<option' + (v._tu === u ? ' selected' : '') + '>' + u + '</option>').join('') + '</select></label>' : '') + '</div>';
    };
    const buildForm = () => {
      let h = unitsBar(), open = false;
      spec.inputs.forEach(f => {
        if (f.group) {
          if (open) h += '</fieldset>';
          h += '<fieldset class="grp"' + (f.if ? ' data-gif="1"' : '') + '><legend><span class="gl">' + UI.lab(f.group) + '</span></legend>';
          open = true;
          return;
        }
        if (!open) { h += '<fieldset class="grp">'; open = true; }
        h += fieldHTML(f, v, spec.id);
      });
      if (open) h += '</fieldset>';
      h += '<div class="form-act"><button type="button" class="btn-s reset">Reset to defaults</button></div>';
      form.innerHTML = h;
      vis();
    };
    const vis = () => {
      spec.inputs.forEach((f, i) => {
        if (f.k && f.if) {
          const el = $('.f[data-k="' + f.k + '"]', form);
          if (el) el.hidden = !f.if(v);
        }
        if (f.group && f.if) {
          const gs = $$('fieldset.grp', form);
          let gi = -1;
          for (let j = 0; j <= i; j++) if (spec.inputs[j].group) gi++;
          if (gs[gi + (spec.units ? 0 : 0)]) gs[gi].hidden = !f.if(v);
        }
      });
    };
    let raf = 0, timer = 0;
    const schedule = () => {
      if (spec.slow) { clearTimeout(timer); timer = setTimeout(compute, 180); return; }
      cancelAnimationFrame(raf);
      raf = requestAnimationFrame(compute);
    };
    const save = () => store.set('st:v:' + spec.id, v);

    form.addEventListener('input', e => {
      const t = e.target, fEl = t.closest('.f');
      if (t.dataset.u) { v[t.dataset.u] = t.value; buildForm(); schedule(); save(); return; }
      if (!fEl) return;
      const k = fEl.dataset.k, f = spec.inputs.find(x => x.k === k);
      if (!f) return;
      if (f.t === 'table') {
        if (t.dataset.r !== undefined) {
          const r = +t.dataset.r, c = +t.dataset.c, num = UI.parseNum(t.value);
          if (!v[k][r]) v[k][r] = [];
          v[k][r][c] = t.value.trim() === '' ? NaN : num;
          t.classList.toggle('bad', t.value.trim() !== '' && !isFinite(num));
        }
      } else if (f.t === 'select') v[k] = t.value;
      else if (f.t === 'seg') { if (t.checked) v[k] = t.value; }
      else if (f.t === 'check') v[k] = t.checked;
      else if (f.t === 'text' || f.t === 'area') v[k] = t.value;
      else if (t.classList.contains('rng')) {
        const val = f.scale === 'log' ? Math.pow(10, +t.value) : +t.value;
        const rounded = +val.toPrecision(f.scale === 'log' ? 3 : 6);
        v[k] = rounded;
        const box = $('input.num', fEl);
        box.value = showNum(rounded);
        box.classList.remove('bad');
      } else {
        const num = UI.parseNum(t.value);
        const okv = isFinite(num) && (f.lo === undefined || num >= f.lo) && (f.hi === undefined || num <= f.hi);
        t.classList.toggle('bad', !okv);
        if (!okv) return;
        v[k] = num;
        const rng = $('.rng', fEl);
        if (rng) rng.value = f.scale === 'log' ? Math.log10(Math.max(num, 1e-300)) : num;
      }
      vis();
      refreshUnits();
      save();
      schedule();
    });
    const refreshUnits = () => {
      spec.inputs.forEach(f => {
        if (f.k && typeof f.label === 'function') { const lb = $('.f[data-k="' + f.k + '"] .fl', form); if (lb) lb.innerHTML = UI.lab(f.label(v)); }
        if (!f.k || typeof f.unit !== 'function') return;
        const el = $('.f[data-k="' + f.k + '"] .fu', form);
        if (el) el.innerHTML = ST.tu(unitSub(f.unit, v));
      });
    };
    form.addEventListener('change', e => {
      const t = e.target;
      if (t.tagName === 'SELECT' && !t.dataset.u) { vis(); refreshUnits(); }
    });
    form.addEventListener('click', e => {
      const b = e.target.closest('button');
      if (!b) return;
      const fEl = b.closest('.f');
      if (b.classList.contains('reset')) {
        Object.keys(v).forEach(k => delete v[k]);
        Object.assign(v, JSON.parse(JSON.stringify(defaults)));
        state.preset = null;
        store.del('st:v:' + spec.id);
        markPreset();
        buildForm(); compute();
        return;
      }
      if (!fEl) return;
      const k = fEl.dataset.k;
      if (b.classList.contains('add')) { v[k].push(new Array(spec.inputs.find(x => x.k === k).cols.length).fill(NaN)); }
      else if (b.classList.contains('clr')) { v[k] = [[NaN, NaN].slice(0, spec.inputs.find(x => x.k === k).cols.length)]; }
      else if (b.classList.contains('rm')) { v[k].splice(+b.dataset.r, 1); if (!v[k].length) v[k].push([NaN, NaN]); }
      else return;
      const html = tableHTML(spec.inputs.find(x => x.k === k), v, 'in-' + spec.id + '-' + k);
      fEl.outerHTML = html;
      save(); schedule();
    });
    form.addEventListener('paste', e => {
      const t = e.target;
      if (!t.dataset || t.dataset.r === undefined) return;
      const text = (e.clipboardData || G.clipboardData).getData('text');
      if (!/[\t\n,;]/.test(text.trim())) return;
      e.preventDefault();
      const fEl = t.closest('.f'), k = fEl.dataset.k, f = spec.inputs.find(x => x.k === k);
      const lines = text.trim().split(/\r?\n/).map(l => l.split(/\t|;|,(?=\s*[-+.\d])|\s{2,}/).map(c => UI.parseNum(c)));
      const r0 = +t.dataset.r, c0 = +t.dataset.c;
      lines.forEach((cells, i) => {
        if (cells.every(c => !isFinite(c))) return;
        const r = r0 + i;
        while (v[k].length <= r) v[k].push(new Array(f.cols.length).fill(NaN));
        cells.forEach((c, j) => { if (c0 + j < f.cols.length) v[k][r][c0 + j] = c; });
      });
      v[k] = v[k].filter((row, i) => i < r0 || row.some(x => isFinite(x)));
      fEl.outerHTML = tableHTML(f, v, 'in-' + spec.id + '-' + k);
      save(); schedule();
    });

    /* presets */
    const markPreset = () => {
      $$('.chip[data-p]', host).forEach(c => c.setAttribute('aria-pressed', String(state.preset !== null && +c.dataset.p === state.preset)));
    };
    $$('.chip[data-p]', host).forEach(c => c.addEventListener('click', () => {
      const p = spec.presets[+c.dataset.p];
      Object.keys(v).forEach(k => delete v[k]);
      Object.assign(v, JSON.parse(JSON.stringify(defaults)), JSON.parse(JSON.stringify(p.v || {})));
      state.preset = +c.dataset.p;
      markPreset();
      save();
      buildForm(); compute();
    }));

    /* outputs */
    function compute() {
      let r;
      const t0 = performance.now();
      try { r = spec.run(v); } catch (err) {
        out.innerHTML = '<div class="note bad"><b>Can’t compute.</b> ' + esc(err.message || String(err)) + '</div>' + (state.lastHTML ? '<div class="stale">' + state.lastHTML + '</div>' : '');
        return;
      }
      state.last = r;
      state.ms = performance.now() - t0;
      renderOut(r);
    }
    function kpiHTML(k) {
      const val = typeof k.v === 'number' ? N.fmt(k.v, k.d || 4, { html: true }) : esc(k.v === undefined ? '–' : k.v);
      return '<div class="kpi' + (k.state ? ' ' + k.state : '') + (k.main ? ' main' : '') + '"><div class="kl">' + UI.lab(k.l) + '</div><div class="kv">' + val + (k.u ? '<span class="ku">' + ST.tu(unitSub(k.u, v)) + '</span>' : '') + '</div>' + (k.s ? '<div class="ks">' + UI.lab(k.s) + '</div>' : '') + '</div>';
    }
    function checkHTML(r) {
      if (state.preset === null) return '';
      const p = spec.presets[state.preset];
      if (!p || (!p.book && !p.expect && !p.check)) return p && p.note ? '<div class="bookcheck"><b>' + UI.lab(p.name) + '.</b> ' + UI.lab(p.note) + '</div>' : '';
      const rowsOf = (src, word) => Object.keys(src || {}).map(key => {
        const b = src[key], kpi = (r.kpi || []).find(k => k.k === key);
        const bv = typeof b === 'object' ? b.v : b, tol = (typeof b === 'object' && b.tol) || 0.03;
        if (!kpi) return '';
        if (typeof bv === 'string') { const okS = String(kpi.v) === bv; return '<span class="bc ' + (okS ? 'ok' : 'off') + '">' + UI.lab(kpi.l) + ': ' + word + ' “' + esc(bv) + '”, here “' + esc(String(kpi.v)) + '”' + (okS ? ' <b>✓</b>' : '') + '</span>'; }
        if (typeof kpi.v !== 'number') return '';
        const rel = Math.abs(kpi.v - bv) / Math.max(Math.abs(bv), 1e-12);
        const ok = rel <= tol;
        return '<span class="bc ' + (ok ? 'ok' : 'off') + '">' + UI.lab(kpi.l) + ': ' + word + ' ' + N.fmt(bv, 4, { html: true }) + ', here ' + N.fmt(kpi.v, 4, { html: true }) + (ok ? ' <b>✓</b>' : ' <b>(' + (rel * 100).toFixed(1) + '% off)</b>') + '</span>';
      }).join('');
      const rows = rowsOf(p.book, 'book') + rowsOf(p.expect, 'book’s formula') + rowsOf(p.check, 'expected');
      return '<div class="bookcheck"><b>' + UI.lab(p.name) + '.</b> ' + (p.note ? UI.lab(p.note) + ' ' : '') + rows + '</div>';
    }
    function tableOut(t, i) {
      const cols = t.cols.map(c => '<th>' + UI.lab(c) + '</th>').join('');
      const body = t.rows.map(r => '<tr>' + r.map(c => '<td>' + (typeof c === 'number' ? N.fmt(c, t.d || 4, { html: true }) : UI.lab(String(c))) + '</td>').join('') + '</tr>').join('');
      return '<figure class="tb"><figcaption><span>' + UI.lab(t.title || 'Table') + '</span><span class="acts"><button type="button" class="btn-s" data-csv="' + i + '">Copy CSV</button><button type="button" class="btn-s" data-tex="' + i + '">Copy LaTeX</button></span></figcaption><div class="tscroll"><table class="rt"><thead><tr>' + cols + '</tr></thead><tbody>' + body + '</tbody></table></div></figure>';
    }
    function renderOut(r) {
      let h = '';
      h += checkHTML(r);
      if (r.kpi && r.kpi.length) h += '<div class="kpis">' + r.kpi.map(kpiHTML).join('') + '</div>';
      (r.warn || []).forEach(w => { h += '<div class="note warn"><span class="hz" aria-hidden="true"></span><span>' + UI.lab(w) + '</span></div>'; });
      (r.info || []).forEach(w => { h += '<div class="note info"><span>' + UI.lab(w) + '</span></div>'; });
      if (r.diagram) h += '<figure class="diag">' + r.diagram + '</figure>';
      if (r.plots && r.plots.length) {
        h += '<div class="plots n' + Math.min(r.plots.length, 2) + '">' + r.plots.map((p, i) => '<figure class="pl' + (p.wide ? ' wide' : '') + '"><figcaption><span>' + UI.lab(p.title || '') + '</span><span class="acts"><button type="button" class="btn-s" data-hold="' + i + '">Hold</button>' + (state.holds[i] && state.holds[i].length ? '<button type="button" class="btn-s" data-unhold="' + i + '">Clear held</button>' : '') + '<button type="button" class="btn-s" data-pcsv="' + i + '">Copy data</button></span></figcaption><div class="plot-host" data-pi="' + i + '"></div>' + (p.cap ? '<p class="cap">' + UI.lab(p.cap) + '</p>' : '') + '</figure>').join('') + '</div>';
      }
      (r.tables || (r.table ? [r.table] : [])).forEach((t, i) => { h += tableOut(t, i); });
      if (r.html) h += r.html;
      out.innerHTML = h;
      state.lastHTML = '';
      (r.plots || []).forEach((p, i) => {
        const hostEl = $('.plot-host[data-pi="' + i + '"]', out);
        const ps = Object.assign({}, p);
        delete ps.title;
        ST.Plot.render(hostEl, ps, { holds: state.holds[i] });
      });
      if (r.after) r.after(out);
    }
    out.addEventListener('click', e => {
      const b = e.target.closest('button');
      if (!b || !state.last) return;
      const r = state.last;
      if (b.dataset.hold !== undefined) {
        const i = +b.dataset.hold, p = r.plots[i];
        state.holds[i] = (state.holds[i] || []).concat((p.series || []).filter(s => !s.hidden && s.type !== 'pts' && s.hold !== false).map(s => ({ x: s.x.slice(), y: s.y.slice(), color: s.color })));
        renderOut(r);
      } else if (b.dataset.unhold !== undefined) {
        state.holds[+b.dataset.unhold] = [];
        renderOut(r);
      } else if (b.dataset.pcsv !== undefined) {
        UI.copy(ST.Plot.csv(r.plots[+b.dataset.pcsv]), b);
      } else if (b.dataset.csv !== undefined) {
        const t = (r.tables || [r.table])[+b.dataset.csv];
        UI.copy([t.cols.map(c => N.fmtCSV(ST.tlPlain(c))).join(',')].concat(t.rows.map(row => row.map(c => N.fmtCSV(typeof c === 'string' ? ST.tlPlain(c) : c)).join(','))).join('\n'), b);
      } else if (b.dataset.tex !== undefined) {
        const t = (r.tables || [r.table])[+b.dataset.tex];
        const texCell = c => (typeof c === 'number' ? (isFinite(c) ? N.fmt(c, t.d || 4, { plain: true }).replace('−', '-') : '--') : ST.tlPlain(String(c)));
        const lines = ['\\begin{tabular}{' + 'l'.repeat(t.cols.length) + '}', '\\hline', t.cols.map(c => '$' + c + '$').join(' & ') + ' \\\\', '\\hline']
          .concat(t.rows.map(row => row.map(texCell).join(' & ') + ' \\\\'), ['\\hline', '\\end{tabular}']);
        UI.copy(lines.join('\n'), b);
      }
    });

    /* sweep */
    const sw = $('.sweep', host);
    if (sw) setupSweep(sw, spec, v, () => state.last);

    buildForm();
    markPreset();
    compute();
    return { recompute: compute, v };
  };

  function headHTML(spec) {
    const B = ST.BOOK, ch = B && B.ch[spec.ch];
    const chips = (spec.presets || []).map((p, i) => '<button type="button" class="chip" data-p="' + i + '" aria-pressed="false">' + UI.lab(p.name) + '</button>').join('');
    return '<header class="tool-head"><div class="th-row"><span class="tag tag-l">' + esc(spec.tag) + '</span><span class="th-ref">' + (spec.ch === 31 ? 'Appendix' : 'Chapter ' + spec.ch) + (ch ? ' · ' + esc(ch.t) : '') + (UI.secRef(spec.sec) ? ' · ' + esc(UI.secRef(spec.sec)) : '') + '</span></div>' +
      '<h1>' + esc(spec.title) + '</h1><p class="lede">' + UI.lab(spec.lede || '') + '</p>' +
      (chips ? '<div class="presets"><span class="pl-l">Load</span>' + chips + '</div>' : '') + '</header>';
  }
  function methodHTML(spec) {
    if (!spec.eqs && !spec.method) return '';
    let h = '<section class="method"><h2>Method</h2><div class="m-grid">';
    if (spec.eqs) h += '<div class="eqs">' + spec.eqs.map(e => '<div class="eq"><div class="eq-m">' + UI.texHTML(e.tex || e[0], true) + '</div>' + ((e.ref || e[1]) ? '<div class="eq-r">' + esc(e.ref || e[1]) + '</div>' : '') + ((e.note || e[2]) ? '<p class="eq-n">' + UI.lab(e.note || e[2]) + '</p>' : '') + '</div>').join('') + '</div>';
    if (spec.method) h += '<div class="m-txt">' + spec.method.map(p => '<p>' + UI.lab(p) + '</p>').join('') + '</div>';
    return h + '</div></section>';
  }
  function sweepHTML(spec) {
    if (spec.sweep === false) return '';
    return '<details class="sweep"><summary>Parameter sweep</summary><div class="sw-body"></div></details>';
  }
  function setupSweep(el, spec, v, getLast) {
    const body = $('.sw-body', el);
    const nums = spec.inputs.filter(f => f.k && !f.t && f.min !== undefined);
    if (!nums.length) { el.remove(); return; }
    let built = false;
    el.addEventListener('toggle', () => {
      if (!el.open || built) return;
      built = true;
      const last = getLast();
      const kp = (last && last.kpi || []).filter(k => typeof k.v === 'number' && k.k);
      body.innerHTML = '<p class="fh">Vary one input across a range and plot any result. Add a second input to draw a family of curves.</p><div class="sw-grid">' +
        '<label>Vary<select class="sw-x">' + nums.map(f => '<option value="' + f.k + '">' + esc(ST.tlPlain(f.label)) + '</option>').join('') + '</select></label>' +
        '<label>From<input class="sw-a" type="text" inputmode="decimal"></label><label>To<input class="sw-b" type="text" inputmode="decimal"></label>' +
        '<label>Points<input class="sw-n" type="text" inputmode="numeric" value="30"></label>' +
        '<label class="chk"><input type="checkbox" class="sw-log"><span>Log spacing</span></label>' +
        '<label>Result<select class="sw-y">' + kp.map(k => '<option value="' + k.k + '">' + esc(ST.tlPlain(k.l)) + '</option>').join('') + '</select></label>' +
        '<label>Family of<select class="sw-f"><option value="">(none)</option>' + nums.map(f => '<option value="' + f.k + '">' + esc(ST.tlPlain(f.label)) + '</option>').join('') + '</select></label>' +
        '<label>Family values<input class="sw-fv" type="text" placeholder="e.g. 0.5, 1, 2"></label>' +
        '</div><div class="sw-run"><button type="button" class="btn">Run sweep</button><span class="sw-msg fh"></span></div><div class="sw-out"></div>';
      const setRange = () => {
        const f = nums.find(x => x.k === $('.sw-x', body).value);
        $('.sw-a', body).value = showNum(f.min > 0 || f.scale !== 'log' ? f.min : f.max * 1e-3);
        $('.sw-b', body).value = showNum(f.max);
        $('.sw-log', body).checked = f.scale === 'log';
      };
      $('.sw-x', body).addEventListener('change', setRange);
      setRange();
      $('.sw-run .btn', body).addEventListener('click', () => {
        const xk = $('.sw-x', body).value, yk = $('.sw-y', body).value, fk = $('.sw-f', body).value;
        const a = UI.parseNum($('.sw-a', body).value), b = UI.parseNum($('.sw-b', body).value), n = Math.min(120, Math.max(3, parseInt($('.sw-n', body).value, 10) || 30));
        const lg = $('.sw-log', body).checked;
        const fvals = fk ? $('.sw-fv', body).value.split(/[,\s]+/).map(UI.parseNum).filter(isFinite).slice(0, 8) : [null];
        if (!isFinite(a) || !isFinite(b) || (lg && (a <= 0 || b <= 0))) { $('.sw-msg', body).textContent = 'Check the range (log spacing needs positive limits).'; return; }
        const xs = lg ? N.logspace(a, b, n) : N.linspace(a, b, n);
        const t0 = performance.now();
        const series = [];
        let fails = 0;
        (fvals.length ? fvals : [null]).forEach((fv, fi) => {
          const ys = xs.map(x => {
            const vv = Object.assign({}, v, { [xk]: x });
            if (fk && fv !== null) vv[fk] = fv;
            try { const r = spec.run(vv), k = (r.kpi || []).find(q => q.k === yk); return k && typeof k.v === 'number' ? k.v : NaN; } catch (e) { fails++; return NaN; }
          });
          const ff = nums.find(q => q.k === fk);
          series.push({ name: fk && fv !== null ? ST.tlPlain(ff.label).replace(/,.*$/, '') + ' = ' + N.fmt(fv, 3) : (ST.tlPlain((getLast().kpi.find(q => q.k === yk) || {}).l || yk)), x: xs, y: ys, color: fi + 1 });
        });
        const fx = nums.find(q => q.k === xk), ky = (getLast().kpi || []).find(q => q.k === yk) || {};
        const host = $('.sw-out', body);
        host.innerHTML = '<div class="plot-host"></div><div class="sw-acts"><button type="button" class="btn-s">Copy data</button></div>';
        const pspec = { x: { label: fx.label.replace(/\$/g, ''), log: lg }, y: { label: (ky.l || yk).replace(/\$/g, ''), unit: ky.u }, series };
        ST.Plot.render($('.plot-host', host), pspec);
        $('.sw-acts .btn-s', host).addEventListener('click', e => UI.copy(ST.Plot.csv(pspec), e.target));
        $('.sw-msg', body).textContent = (n * fvals.length) + ' runs in ' + Math.round(performance.now() - t0) + ' ms' + (fails ? ' · ' + fails + ' runs failed' : '');
      });
    });
  }

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
