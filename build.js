/* Build: one self-contained index.html with fonts, KaTeX styles and pre-rendered equations inlined. */
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const katex = require('katex');

const ROOT = __dirname;
const SRC = path.join(ROOT, 'src');
const OUT = path.join(ROOT, 'dist');
const NM = path.join(ROOT, 'node_modules');

const CORE_PRE = ['num', 'expr', 'net', 'kin', 'kit', 'tex', 'plot', 'therm', 'rtd', 'pellet', 'draw', 'book', 'ui'].map(n => path.join(SRC, 'core', n + '.js'));
const CORE_POST = ['app', 'home'].map(n => path.join(SRC, 'core', n + '.js'));
const TOOL_FILES = fs.readdirSync(path.join(SRC, 'tools')).filter(f => f.endsWith('.js')).sort().map(f => path.join(SRC, 'tools', f));

/* ---- load specs in node to collect equations ---- */
const ctx = { console, Math, module: undefined };
ctx.globalThis = ctx;
vm.createContext(ctx);
for (const f of CORE_PRE.concat(TOOL_FILES)) vm.runInContext(fs.readFileSync(f, 'utf8'), ctx, { filename: f });
const ST = ctx.ST;
const TEX = {};
let nEq = 0, bad = 0;
ST.TOOLS.forEach(t => {
  (t.eqs || []).forEach(e => {
    const tex = e.tex || e[0];
    const key = 'D:' + tex;
    if (TEX[key]) return;
    try { TEX[key] = katex.renderToString(tex, { displayMode: true, throwOnError: true, strict: false }); nEq++; }
    catch (err) { bad++; console.error('TeX error in', t.id, ':', err.message.slice(0, 160)); TEX[key] = katex.renderToString(tex, { displayMode: true, throwOnError: false }); }
  });
});
const ids = new Set();
ST.TOOLS.forEach(t => { if (ids.has(t.id)) console.error('DUPLICATE id', t.id); ids.add(t.id); });

/* ---- fonts ---- */
const b64 = p => fs.readFileSync(p).toString('base64');
const fsrc = (pkg, file) => path.join(NM, '@fontsource', pkg, 'files', file);
const faces = [
  ['Big Shoulders Display', 700, 'normal', fsrc('big-shoulders-display', 'big-shoulders-display-latin-700-normal.woff2')],
  ['Big Shoulders Display', 800, 'normal', fsrc('big-shoulders-display', 'big-shoulders-display-latin-800-normal.woff2')],
  ['Big Shoulders Stencil Display', 800, 'normal', fsrc('big-shoulders-stencil-display', 'big-shoulders-stencil-display-latin-800-normal.woff2')],
  ['Atkinson Hyperlegible Next', 400, 'normal', fsrc('atkinson-hyperlegible-next', 'atkinson-hyperlegible-next-latin-400-normal.woff2')],
  ['Atkinson Hyperlegible Next', 400, 'italic', fsrc('atkinson-hyperlegible-next', 'atkinson-hyperlegible-next-latin-400-italic.woff2')],
  ['Atkinson Hyperlegible Next', 700, 'normal', fsrc('atkinson-hyperlegible-next', 'atkinson-hyperlegible-next-latin-700-normal.woff2')],
  ['Atkinson Hyperlegible Next', 800, 'normal', fsrc('atkinson-hyperlegible-next', 'atkinson-hyperlegible-next-latin-800-normal.woff2')],
  ['Atkinson Hyperlegible Mono', 400, 'normal', fsrc('atkinson-hyperlegible-mono', 'atkinson-hyperlegible-mono-latin-400-normal.woff2')],
  ['Atkinson Hyperlegible Mono', 500, 'normal', fsrc('atkinson-hyperlegible-mono', 'atkinson-hyperlegible-mono-latin-500-normal.woff2')]
];
let fontCSS = faces.map(([fam, w, st, p]) => '@font-face{font-family:"' + fam + '";font-style:' + st + ';font-weight:' + w + ';font-display:swap;src:url(data:font/woff2;base64,' + b64(p) + ') format("woff2")}').join('\n');

/* KaTeX css with fonts inlined (woff2 only) */
let kcss = fs.readFileSync(path.join(NM, 'katex', 'dist', 'katex.min.css'), 'utf8');
kcss = kcss.replace(/src:url\(fonts\/(KaTeX_[A-Za-z0-9_-]+)\.woff2\) format\("woff2"\),url\(fonts\/[^)]+\) format\("woff"\),url\(fonts\/[^)]+\) format\("truetype"\)/g,
  (m, name) => 'src:url(data:font/woff2;base64,' + b64(path.join(NM, 'katex', 'dist', 'fonts', name + '.woff2')) + ') format("woff2")');
if (/url\(fonts\//.test(kcss)) console.error('KaTeX font URL left un-inlined');

const css = fs.readFileSync(path.join(SRC, 'css', 'site.css'), 'utf8');
const js = CORE_PRE.concat(TOOL_FILES, CORE_POST).map(f => '/* ' + path.relative(SRC, f) + ' */\n' + fs.readFileSync(f, 'utf8')).join('\n');

const html = '<title>Spacetime</title>\n<meta name="claude:viewport-fit-sides" content="cover">\n<meta name="description" content="Interactive calculators for every chapter of Levenspiel, Chemical Reaction Engineering, 3rd edition.">\n' +
  '<style>\n' + fontCSS + '\n' + kcss + '\n' + css + '\n</style>\n' +
  '<div id="app"></div>\n' +
  '<script>\n' + fs.readFileSync(path.join(NM, 'katex', 'dist', 'katex.min.js'), 'utf8').replace(/<\/script/gi, '<\\/script') + '\n</script>\n' +
  '<script>\n' + js.replace(/<\/script/gi, '<\\/script') + '\nST.start();\n</script>\n';

fs.mkdirSync(OUT, { recursive: true });
fs.writeFileSync(path.join(OUT, 'index.html'), html);
/* local preview wrapper that mimics the hosting skeleton */
fs.writeFileSync(path.join(OUT, 'preview.html'), '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover"><style>:root{color-scheme:light;padding-top:env(safe-area-inset-top);padding-bottom:env(safe-area-inset-bottom)}body{margin:0;font:14px system-ui}img{max-width:100%}[hidden]{display:none!important}</style></head><body>\n' + html + '</body></html>');
console.log('tools:', ST.TOOLS.length, '| equations:', nEq, (bad ? '| TeX errors: ' + bad : ''), '| size:', (Buffer.byteLength(html) / 1024).toFixed(0) + ' KB');
