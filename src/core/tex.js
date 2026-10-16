/* Light math markup for labels: C_{A0}, X_A, \tau, k^{-1}. HTML or SVG tspans. */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  const GREEK = {
    alpha: 'α', beta: 'β', gamma: 'γ', delta: 'δ', epsilon: 'ϵ', varepsilon: 'ε', zeta: 'ζ', eta: 'η', theta: 'θ',
    vartheta: 'ϑ', kappa: 'κ', lambda: 'λ', mu: 'μ', nu: 'ν', xi: 'ξ', pi: 'π', rho: 'ρ', sigma: 'σ', tau: 'τ',
    upsilon: 'υ', phi: 'ϕ', varphi: 'φ', chi: 'χ', psi: 'ψ', omega: 'ω',
    Gamma: 'Γ', Delta: 'Δ', Theta: 'Θ', Lambda: 'Λ', Xi: 'Ξ', Pi: 'Π', Sigma: 'Σ', Phi: 'Φ', Psi: 'Ψ', Omega: 'Ω'
  };
  const SYM = {
    cdot: '·', times: '×', to: '→', infty: '∞', le: '≤', ge: '≥', leq: '≤', geq: '≥', approx: '≈', pm: '±',
    partial: '∂', sqrt: '√', int: '∫', sum: 'Σ', rightleftharpoons: '⇌', prime: '′', circ: '°', ',': ' ', ' ': ' ',
    ell: 'ℓ', propto: '∝', deg: '°', bar: '', ldots: '…', neq: '≠', sim: '∼', degree: '°'
  };
  const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

  /* tokens: {k:'ital'|'up'|'sub'|'sup', s} */
  function lex(src) {
    const s = String(src);
    const explicit = s.indexOf('$') >= 0;
    const prose = !explicit && (s.length > 30 || (s.match(/ /g) || []).length > 3);
    let math = !explicit && !prose;
    const out = [];
    let i = 0;
    const isL = c => /[A-Za-z]/.test(c);
    function group() {
      if (s[i] === '{') {
        let depth = 1, j = i + 1;
        while (j < s.length && depth) { if (s[j] === '{') depth++; else if (s[j] === '}') depth--; j++; }
        const inner = s.slice(i + 1, j - 1);
        i = j;
        return inner;
      }
      if (s[i] === '\\') {
        const m = /^\\([A-Za-z]+)/.exec(s.slice(i));
        if (m) { i += m[0].length; return m[0]; }
      }
      /* unbraced: a run of letters and digits, stopping where a new symbol plainly starts
         (digit → letter, lower → upper case, or a letter that carries its own subscript) */
      const c = s[i];
      if (!/[A-Za-z0-9]/.test(c || '')) { i++; return c || ''; }
      let j = i + 1;
      while (j < s.length && /[A-Za-z0-9]/.test(s[j])) {
        const p = s[j - 1], q = s[j];
        if (/[0-9]/.test(p) && /[A-Za-z]/.test(q)) break;
        if (/[a-z]/.test(p) && /[A-Z]/.test(q)) break;
        j++;
      }
      if (j - i >= 2 && s[j] === '_') j--;
      const g = s.slice(i, j);
      i = j;
      return g;
    }
    while (i < s.length) {
      const c = s[i];
      if (c === '$') { math = !math; i++; continue; }
      if (c === '\\') {
        const m = /^\\([A-Za-z]+|[, ])/.exec(s.slice(i));
        if (m) {
          i += m[0].length;
          const w = m[1];
          if (GREEK[w]) out.push({ k: /[A-Z]/.test(w[0]) ? 'up' : 'ital', s: GREEK[w] });
          else if (w in SYM) out.push({ k: 'up', s: SYM[w] });
          else if (w === 'mathrm' || w === 'text' || w === 'rm') { const g = group(); out.push({ k: 'up', s: g }); }
          else out.push({ k: 'up', s: w });
          continue;
        }
      }
      if (c === '_' || c === '^') {
        i++;
        const g = group();
        out.push({ k: c === '_' ? 'sub' : 'sup', s: g });
        continue;
      }
      if (c === '{' || c === '}') { i++; continue; }
      if (isL(c)) {
        let j = i;
        while (j < s.length && isL(s[j])) j++;
        const word = s.slice(i, j);
        if (math && word.length === 1) out.push({ k: 'ital', s: word });
        else if (prose && word.length === 1 && (s[j] === '_' || s[j] === '^' || s[j] === '*')) out.push({ k: 'ital', s: word });
        else if (math && explicit) word.split('').forEach(ch => out.push({ k: 'ital', s: ch }));
        else out.push({ k: 'up', s: word });
        i = j;
        continue;
      }
      if (c === '-' && prose && /[rR]/.test(s[i + 1] || '') && (s[i + 2] === '_')) { out.push({ k: 'up', s: '\u2212' }); i++; continue; }
      if (c === '-' && math) {
        const prev = s[i - 1], next = s[i + 1];
        out.push({ k: 'up', s: prev && next && isL(prev) && isL(next) && !explicit ? '-' : '−' });
        i++;
        continue;
      }
      if (c === '*') {
        /* C_R*, k*: a starred symbol when it closes a symbol; otherwise a product in math */
        let j = i;
        while (s[j] === '*') j++;
        const nx = s[j];
        if (i > 0 && /[A-Za-z0-9}Ͱ-Ͽ]/.test(s[i - 1]) && (nx === undefined || /[\s)\/,;:\]}.^_\-−+=]/.test(nx))) { out.push({ k: 'sup', s: '∗'.repeat(j - i) }); i = j; continue; }
        if (math) { out.push({ k: 'up', s: '·' }); i++; continue; }
      }
      out.push({ k: 'up', s: c });
      i++;
    }
    return out;
  }

  /* HTML */
  ST.tl = function (src) {
    if (src === undefined || src === null) return '';
    return lex(src).map(t => {
      if (t.k === 'ital') return '<i class="mi">' + esc(t.s) + '</i>';
      if (t.k === 'sub') return '<sub>' + ST.tlSub(t.s) + '</sub>';
      if (t.k === 'sup') return '<sup>' + ST.tlSub(t.s, true) + '</sup>';
      return esc(t.s).replace(/−/g, '<span class="minus">−</span>').replace(/([\u0370-\u03FF])/g, '<span class="gk">$1</span>');
    }).join('');
  };
  /* units: upright letters, only sub/superscripts processed */
  ST.tu = function (src) {
    if (src === undefined || src === null) return '';
    return String(src).replace(/&/g, '&amp;').replace(/</g, '&lt;')
      .replace(/\^\{([^}]*)\}/g, '<sup>$1</sup>').replace(/\^(-?[0-9.]+)/g, '<sup>$1</sup>')
      .replace(/_\{([^}]*)\}/g, '<sub>$1</sub>').replace(/-/g, '\u2212');
  };
  /* subscripts upright except greek; superscripts keep italic letters */
  ST.tlSub = function (src, sup) {
    return lex(src).map(t => {
      if (t.k === 'ital') return sup || /[^A-Za-z]/.test(t.s) ? '<i class="mi">' + esc(t.s) + '</i>' : esc(t.s);
      if (t.k === 'sub') return '<sub>' + ST.tlSub(t.s) + '</sub>';
      if (t.k === 'sup') return '<sup>' + ST.tlSub(t.s, true) + '</sup>';
      return esc(t.s);
    }).join('');
  };
  /* plain text (for CSV headers, aria labels) */
  ST.tlPlain = function (src) {
    return lex(src).map(t => (t.k === 'sub' ? '_' + ST.tlPlain(t.s) : t.k === 'sup' ? '^' + ST.tlPlain(t.s) : t.s)).join('');
  };
  /* SVG tspans; fs = base font size in px */
  ST.tlSvg = function (src, fs = 12) {
    const sh = (fs * 0.32).toFixed(2), small = (fs * 0.74).toFixed(1);
    let pend = 0, out = '';
    const open = extra => { const d = pend ? ' dy="' + (-pend).toFixed(2) + '"' : ''; pend = 0; return '<tspan' + d + (extra || '') + '>'; };
    lex(src).forEach(t => {
      if (t.k === 'ital') out += open(' class="mi"') + esc(t.s) + '</tspan>';
      else if (t.k === 'sub' || t.k === 'sup') {
        const d = t.k === 'sub' ? +sh : -sh * 1.25;
        const inner = lex(t.s).map(u => (u.k === 'ital' && (t.k === 'sup' || /[^A-Za-z]/.test(u.s)) ? '<tspan class="mi">' + esc(u.s) + '</tspan>' : esc(u.s))).join('');
        out += '<tspan dy="' + (d - pend).toFixed(2) + '" font-size="' + small + '">' + inner + '</tspan>';
        pend = d;
      } else out += open() + esc(t.s) + '</tspan>';
    });
    if (pend) out += '<tspan dy="' + (-pend).toFixed(2) + '">​</tspan>';
    return out;
  };

  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
