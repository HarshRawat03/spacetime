/* Shell: header search, index sidebar, routing between home, chapter and tool pages. */
(function (G) {
  'use strict';
  const ST = G.ST;
  const UI = ST.UI, B = ST.BOOK;
  const $ = (s, r = document) => r.querySelector(s);
  const esc = UI.esc;
  const tools = () => ST.TOOLS;
  const byCh = ch => tools().filter(t => t.ch === ch);
  const chTitle = ch => (B.ch[ch] ? B.ch[ch].t : '');

  function sideHTML() {
    let h = '';
    B.parts.forEach(p => {
      const list = p.ch.filter(c => byCh(c).length);
      if (!list.length) return;
      h += '<h3><span>' + esc(p.n) + '</span>' + esc(p.t) + '</h3>';
      list.forEach(c => {
        h += '<div class="chap"><em>' + (c === 31 ? 'A' : c) + '</em><a href="#ch' + c + '">' + esc(chTitle(c)) + '</a></div><ul>';
        byCh(c).forEach(t => { h += '<li><a href="#' + t.id + '" data-id="' + t.id + '"><span class="t">' + esc(t.tag) + '</span><span>' + esc(t.title) + '</span></a></li>'; });
        h += '</ul>';
      });
    });
    return h;
  }

  /* ---------- search ---------- */
  function norm(s) { return String(s || '').toLowerCase().normalize('NFKD').replace(/[̀-ͯ]/g, '').replace(/[^a-z0-9.\- ]+/g, ' '); }
  let INDEX = null;
  function buildIndex() {
    INDEX = tools().map(t => ({
      t,
      title: norm(t.title), tag: norm(t.tag),
      body: norm([t.lede, t.keys, chTitle(t.ch), 'chapter ' + t.ch, t.sec ? 'section ' + t.sec : '', (t.presets || []).map(p => p.name).join(' ')].join(' '))
    }));
  }
  function search(q) {
    if (!INDEX) buildIndex();
    const words = norm(q).split(/\s+/).filter(Boolean);
    if (!words.length) return [];
    const out = [];
    INDEX.forEach(e => {
      let score = 0;
      for (const w of words) {
        const inT = e.title.indexOf(w) >= 0, inTag = e.tag.indexOf(w) >= 0, inB = e.body.indexOf(w) >= 0;
        if (!inT && !inTag && !inB) return;
        score += (inTag ? 6 : 0) + (inT ? (e.title.startsWith(w) ? 5 : 3) : 0) + (inB ? 1 : 0);
      }
      out.push({ e, score });
    });
    return out.sort((a, b) => b.score - a.score || a.e.t.ch - b.e.t.ch).slice(0, 12).map(r => r.e.t);
  }
  function setupSearch() {
    const inp = $('#q'), box = $('#results');
    let sel = -1, list = [];
    const draw = () => {
      if (!inp.value.trim()) { box.hidden = true; return; }
      list = search(inp.value);
      box.innerHTML = list.length ? list.map((t, i) => '<a href="#' + t.id + '" role="option" aria-selected="' + (i === sel) + '"><span class="tag">' + esc(t.tag) + '</span><b>' + esc(t.title) + '</b><small>Ch. ' + t.ch + ' · ' + esc(chTitle(t.ch)) + '</small></a>').join('') : '<div class="none">Nothing matches. Try a chapter number, a reactor type, or “Example 5.4”.</div>';
      box.hidden = false;
    };
    inp.addEventListener('input', () => { sel = -1; draw(); });
    inp.addEventListener('keydown', e => {
      if (e.key === 'ArrowDown') { sel = Math.min(list.length - 1, sel + 1); draw(); e.preventDefault(); }
      else if (e.key === 'ArrowUp') { sel = Math.max(0, sel - 1); draw(); e.preventDefault(); }
      else if (e.key === 'Enter') { const t = list[sel >= 0 ? sel : 0]; if (t) { location.hash = t.id; inp.value = ''; box.hidden = true; inp.blur(); } }
      else if (e.key === 'Escape') { inp.value = ''; box.hidden = true; inp.blur(); }
    });
    box.addEventListener('click', e => { if (e.target.closest('a')) { inp.value = ''; box.hidden = true; } });
    document.addEventListener('click', e => { if (!e.target.closest('.search')) box.hidden = true; });
    document.addEventListener('keydown', e => {
      if (e.key === '/' && !/INPUT|TEXTAREA|SELECT/.test(document.activeElement.tagName)) { e.preventDefault(); inp.focus(); }
    });
  }

  /* ---------- pages ---------- */
  function partOfCh(ch) { return ST.partOf(ch); }
  function chapterHTML(ch) {
    const p = partOfCh(ch);
    const list = byCh(ch);
    let h = '<div class="chpage"><div class="th-row"><span class="tag tag-l">' + (p ? 'PART ' + esc(p.n) : '') + '</span><span class="th-ref">' + (p ? esc(p.t) : '') + ' · book p. ' + (B.ch[ch] ? B.ch[ch].p : '') + '</span></div>';
    h += '<h1>' + (ch === 31 ? 'Appendix' : 'Chapter ' + ch) + ' · ' + esc(chTitle(ch)) + '</h1>';
    h += '<table class="eqlist"><tbody>' + list.map(rowHTML).join('') + '</tbody></table></div>';
    return h;
  }
  function rowHTML(t) {
    const refs = (t.presets || []).filter(p => /^(Example|Ex\.|Problem|Fig|Table|Eq)/.test(p.name)).map(p => p.name.split(' — ')[0].replace(/^Examples?/, 'Ex.').replace(/^Problems?/, 'P').trim());
    const ex = refs.filter((r, i) => refs.indexOf(r) === i).slice(0, 4).map(r => r.replace(/ /g, '\u00a0')).join(', ');
    return '<tr><td><span class="tag">' + esc(t.tag) + '</span></td><td><a href="#' + t.id + '">' + esc(t.title) + '</a></td><td>' + UI.lab(t.lede || '') + '</td><td>' + esc(UI.secRef(t.sec)) + (ex ? '<span class="exs">' + esc(ex) + '</span>' : '') + '</td></tr>';
  }
  ST.App = {
    rowHTML, chapterHTML, byCh, chTitle,
    plantHTML() {
      let h = '';
      B.parts.forEach(p => {
        const chs = p.ch.filter(c => byCh(c).length);
        if (!chs.length) return;
        const n = chs.reduce((s, c) => s + byCh(c).length, 0);
        h += '<section class="part" id="part-' + p.code + '"><div class="sect-h"><span class="pn">' + esc(p.n) + '</span><h2>' + esc(p.t) + '</h2><p>' + n + ' tool' + (n === 1 ? '' : 's') + (p.code === 'A' ? ' in the appendix' : ' across chapter' + (chs.length > 1 ? 's ' : ' ') + chs.join(', ')) + '</p></div>';
        chs.forEach(c => {
          h += '<table class="eqlist"><caption><em>' + (c === 31 ? 'A' : String(c).padStart(2, '0')) + '</em>' + esc(chTitle(c)) + '</caption><tbody>' + byCh(c).map(rowHTML).join('') + '</tbody></table>';
        });
        h += '</section>';
      });
      return h;
    }
  };

  let current = null;
  function route() {
    const h = decodeURIComponent(location.hash.slice(1));
    const main = $('#main'), frame = $('.frame'), side = $('.side');
    side.classList.remove('open');
    let tool = null;
    if (!h || h === 'top') {
      frame.classList.add('home');
      main.innerHTML = '<div class="main-in">' + ST.Home.html() + '</div>' + footHTML();
      ST.Home.mount(main);
      document.title = 'Spacetime';
    } else if (/^ch\d+$/.test(h) && B.ch[+h.slice(2)]) {
      frame.classList.remove('home');
      main.innerHTML = '<div class="main-in">' + chapterHTML(+h.slice(2)) + '</div>' + footHTML();
    } else if (/^part-[A-Z]$/.test(h)) {
      frame.classList.add('home');
      main.innerHTML = '<div class="main-in">' + ST.Home.html() + '</div>' + footHTML();
      ST.Home.mount(main);
      const el = document.getElementById(h);
      if (el) setTimeout(() => el.scrollIntoView(), 0);
      return mark(null);
    } else if ((tool = tools().find(t => t.id === h))) {
      frame.classList.remove('home');
      main.innerHTML = '<div class="main-in"><article class="tool" id="tool"></article></div>' + footHTML();
      UI.mountTool($('#tool', main), tool);
    } else {
      location.hash = '';
      return;
    }
    mark(tool);
    if (current !== h) window.scrollTo(0, 0);
    current = h;
  }
  function mark(tool) {
    document.querySelectorAll('.side a[data-id]').forEach(a => {
      if (tool && a.dataset.id === tool.id) { a.setAttribute('aria-current', 'page'); a.scrollIntoView({ block: 'nearest' }); }
      else a.removeAttribute('aria-current');
    });
  }
  function footHTML() {
    return '<footer class="foot"><div>A project by Harsh Rawat, Chemical Engineering, IIT Jammu.</div><div>Companion to O. Levenspiel, <i>Chemical Reaction Engineering</i>, 3rd ed. (Wiley, 1999). Equation, figure and example numbers refer to that edition. Not affiliated with the author or the publisher.</div><div>Every calculation runs in your browser. Inputs are remembered on this device only.</div></footer>';
  }

  ST.start = function () {
    const root = document.getElementById('app') || document.body;
    root.innerHTML = '<header class="top"><div class="top-in"><a class="brand" href="#"><b>SPACE<i>·</i>TIME</b><span>Levenspiel’s CRE, computed</span></a>' +
      '<div class="search" role="search"><input id="q" type="search" placeholder="Search tools, chapters, examples" autocomplete="off" aria-label="Search tools" aria-controls="results"><kbd>/</kbd><div class="results" id="results" role="listbox" hidden></div></div>' +
      '<button class="idx-btn" type="button" aria-expanded="false" aria-controls="side">Index</button></div></header>' +
      '<div class="frame"><nav class="side" id="side" aria-label="Book index">' + sideHTML() + '</nav><main id="main"></main></div>';
    setupSearch();
    const btn = $('.idx-btn'), side = $('.side');
    btn.addEventListener('click', () => { const o = side.classList.toggle('open'); btn.setAttribute('aria-expanded', String(o)); });
    side.addEventListener('click', e => { if (e.target.closest('a')) { side.classList.remove('open'); btn.setAttribute('aria-expanded', 'false'); } });
    G.addEventListener('hashchange', route);
    route();
  };
})(window);
