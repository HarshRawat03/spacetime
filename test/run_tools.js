/* Runs every tool on its defaults and presets; prints book checks. Usage: node test/run_tools.js [filter] */
const fs = require('fs'), path = require('path'), vm = require('vm');
const SRC = path.join(__dirname, '..', 'src');
const ctx = { console, Math, module: undefined };
ctx.globalThis = ctx;
vm.createContext(ctx);
const core = ['num', 'expr', 'net', 'kin', 'kit', 'tex', 'plot', 'therm', 'rtd', 'pellet', 'draw', 'book', 'ui'];
core.forEach(n => vm.runInContext(fs.readFileSync(path.join(SRC, 'core', n + '.js'), 'utf8'), ctx, { filename: n }));
fs.readdirSync(path.join(SRC, 'tools')).filter(f => f.endsWith('.js')).sort().forEach(f => vm.runInContext(fs.readFileSync(path.join(SRC, 'tools', f), 'utf8'), ctx, { filename: f }));
const ST = ctx.ST, N = ST.N;
const filt = process.argv[2] || '';
let fails = 0, checks = 0, off = 0;
const fmt = v => (typeof v === 'number' ? N.fmt(v, 5) : String(v));
ST.TOOLS.filter(t => !filt || t.id.indexOf(filt) >= 0 || t.tag.indexOf(filt) >= 0).forEach(t => {
  const defaults = {};
  const seenK = new Set();
  t.inputs.forEach(f => { if (f.k) { if (seenK.has(f.k)) { console.log('  DUPKEY', t.tag, f.k); fails++; } seenK.add(f.k); defaults[f.k] = JSON.parse(JSON.stringify(f.v)); } });
  if (t.units) { defaults._cu = t.units.C || 'mol/L'; defaults._tu = t.units.t || 'min'; }
  const cases = [{ name: 'defaults', v: {} }].concat(t.presets || []);
  cases.forEach(p => {
    const v = Object.assign({}, JSON.parse(JSON.stringify(defaults)), JSON.parse(JSON.stringify(p.v || {})));
    const t0 = Date.now();
    let r;
    try { r = t.run(v); } catch (e) { fails++; console.log('  FAIL', t.tag, t.id, '|', p.name, '|', e.message); return; }
    const ms = Date.now() - t0;
    const k = (r.kpi || []).map(q => q.k + '=' + fmt(q.v)).join('  ');
    let line = '  ' + t.tag.padEnd(7) + (p.name).padEnd(30) + ms + 'ms  ' + k;
    const allChecks = Object.assign({}, p.book || {}, p.expect || {}, p.check || {});
    if (p.book || p.expect || p.check) {
      Object.keys(allChecks).forEach(key => {
        checks++;
        const b = allChecks[key], bv = typeof b === 'object' ? b.v : b, tol = (typeof b === 'object' && b.tol) || 0.03;
        const q = (r.kpi || []).find(z => z.k === key);
        if (!q) { off++; line += '\n      MISSING kpi ' + key; return; }
        if (typeof bv === 'string') { const okS = String(q.v) === bv; if (!okS) off++; line += '\n      ' + (okS ? 'ok  ' : 'OFF ') + key + ': book "' + bv + '" here "' + q.v + '"'; return; }
        const rel = Math.abs(q.v - bv) / Math.max(Math.abs(bv), 1e-12);
        if (!(rel <= tol)) off++;
        line += '\n      ' + (rel <= tol ? 'ok  ' : 'OFF ') + key + ': book ' + bv + ' here ' + fmt(q.v) + ' (' + (rel * 100).toFixed(2) + '%)';
      });
    }
    (r.plots || []).forEach((pl, i) => { (pl.series || []).forEach(s => { if (s.x && s.y && s.x.length !== s.y.length) { fails++; line += '\n      BAD series length in plot ' + i; } }); });
    if (process.env.V) console.log(line + (r.warn && r.warn.length ? '\n      warn: ' + r.warn.join(' | ') : ''));
    else console.log(line.split('\n')[0].slice(0, 220) + (line.indexOf('\n') > 0 ? '\n' + line.split('\n').slice(1).join('\n') : ''));
  });
});
console.log('\n' + ST.TOOLS.length + ' tools | failures: ' + fails + ' | book checks: ' + checks + ' (' + off + ' off)');
