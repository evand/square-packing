/* Overview grid: n = 1..324 coloured by a property of the record packing, with a year slider. */
(() => {
const $ = id => document.getElementById(id);
const NMAX = 324;
let IDX, TL, state = { mode: 'cat', year: 2026, layout: 'grid' };
const K = Math.round(Math.sqrt(NMAX));   // 18: rows of the triangle, one per side k = ceil(sqrt n)
// Where tile n goes.  'grid': 18 per row.  'tri': one row per k, holding n = (k-1)^2+1 .. k^2,
// centred (Ellsworth's triangular table).  'deficit': the same rows right-aligned, so each column
// is one value of k^2 - n and the families s(k^2 - c) = k line up.
function place(n) {
  if (state.layout === 'grid') return null;
  const k = Math.ceil(Math.sqrt(n)), col = state.layout === 'tri' ? (K - k) + (n - (k - 1) * (k - 1)) : 2 * K - 1 - (k * k - n);
  return { row: k, col };
}
// Dark or light numbers by the tile's own background: the tile colours are fixed (or follow the
// theme), so the page's ink colour is wrong on half of them in one theme or the other.
function inkFor(el) {
  const m = getComputedStyle(el).backgroundColor.match(/[\d.]+/g); if (!m) return;
  const L = ([r, g, b]) => [r, g, b].map(v => { v /= 255; return v <= 0.03928 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4; }).reduce((s, v, i) => s + v * [0.2126, 0.7152, 0.0722][i], 0);
  el.style.color = L(m.map(Number)) > 0.2 ? '#131A1E' : '#F2F4F1';
}

const CATS = {
  trivial: ['grid only', 'var(--sq-fill)'],
  diagonal: ['45° tilts only', '#9DC3E6'],
  one: ['one tilt angle', '#F4B183'],
  many: ['several tilt angles', '#C55A5A'],
};
const SYMS = {
  'D4 (four mirror axes)': ['4 mirror axes', '#1F77B4'], 'two mirror axes': ['2 mirror axes', '#6BAED6'], 'one mirror axis': ['1 mirror axis', '#BDD7E7'],
  'rotational (90°) only': ['90° rotation only', '#FDAE6B'], 'rotational (180°) only': ['180° rotation only', '#FDD0A2'], 'none': ['no symmetry', 'var(--sq-fill)'],
};
const proved = new Set();   // n known to be optimal, from Ellsworth's "Proved by" wording; refined by bounds data if present

function recordAsOf(n, year) {
  const L = TL[String(n)] || [];
  // packings for this n with a date <= year; best = smallest s; undated ones count as always known if they're the trivial one
  let best = null;
  for (const d of L) {
    if (!d.s) continue;
    const y = d.date ? d.date[0] : null;
    if (y != null && y > year) continue;
    if (y == null && !d.is_record) continue;   // undated alternatives: skip
    if (y == null && year < 2026) continue;
    // Optimiser starts and files with real errors (51_r2__invalid overlaps) are not packings;
    // "warn:" entries are numeric-precision notes on genuine ones (41_r4, 55_r2a).
    const f = IDX.files[d.file];
    if (f && (f.start || (f.errors || []).some(e => !e.startsWith('warn:')))) continue;
    if (!best || d.s < best.s - 1e-12) best = d;
  }
  return best;
}
function trivialS(n) { return Math.ceil(Math.sqrt(n)); }

function tileInfo(n) {
  const rec = recordAsOf(n, state.year);
  const triv = trivialS(n);
  const file = rec && rec.s < triv - 1e-9 ? rec.file : null;
  const f = file ? IDX.files[file] : null;
  const sum = f && f.summary;
  let cat = 'trivial';
  if (sum) cat = sum.category === 'trivial' ? 'trivial' : sum.category === 'diagonal' ? 'diagonal' : (sum.n_angles <= 2 ? 'one' : 'many');
  return { n, rec, file, sum, cat, s: file ? rec.s : triv, year: rec && rec.date ? rec.date[0] : null };
}

// Something Explore can actually show for this n: either a packing file of its own, or a catalogue
// record under exactly this n (which Explore resolves and annotates itself).
const openable = t => !!(t.file || IDX.records[String(t.n)]);

function colour(t) {
  if (state.mode === 'cat') return CATS[t.cat][1];
  if (state.mode === 'sym') { if (!t.file) return SYMS['D4 (four mirror axes)'][1]; return t.sum ? SYMS[t.sum.symmetry][1] : '#ccc'; }
  if (state.mode === 'free') { if (!t.file) return 'var(--sq-fill)'; if (!t.sum) return '#ccc'; if (t.sum.rigid) return '#2A6E4F'; const k = t.sum.free; return k === 0 ? '#DFEDE5' : k <= 2 ? '#F5EBDA' : k <= 6 ? '#EFC77A' : '#C9862B'; }
  if (state.mode === 'year') { if (!t.file || !t.year) return 'var(--sq-fill)'; const u = (t.year - 1979) / (2026 - 1979); return `hsl(${200 - 160 * u} 60% ${72 - 22 * u}%)`; }
  if (state.mode === 'waste') { const w = t.s * t.s - t.n, frac = w / (t.s * t.s); const u = Math.min(1, frac / 0.16); return `hsl(${120 - 120 * u} 50% ${80 - 25 * u}%)`; }
  return '#ccc';
}
function legend() {
  const L = $('legend'); L.innerHTML = '';
  let items = [];
  if (state.mode === 'cat') items = Object.values(CATS);
  else if (state.mode === 'sym') items = Object.values(SYMS);
  else if (state.mode === 'free') items = [['rigid', '#2A6E4F'], ['no free squares (but some move together)', '#DFEDE5'], ['1–2 free', '#F5EBDA'], ['3–6 free', '#EFC77A'], ['7+ free', '#C9862B'], ['grid only', 'var(--sq-fill)']];
  else if (state.mode === 'year') items = [[1979, colour({ file: 1, year: 1979 })], [1990, colour({ file: 1, year: 1990 })], [2005, colour({ file: 1, year: 2005 })], [2020, colour({ file: 1, year: 2020 })], [2026, colour({ file: 1, year: 2026 })]];
  else if (state.mode === 'waste') items = [['0%', colour({ s: 1, n: 1 })], ['4%', colour({ s: 1, n: .96 })], ['8%', colour({ s: 1, n: .92 })], ['16%+', colour({ s: 1, n: .84 })]];
  for (const [name, col] of items) { const sp = document.createElement('span'); sp.innerHTML = `<i style="background:${col}"></i>${name}`; L.appendChild(sp); }
}

function draw() {
  const G = $('grid'); G.innerHTML = '';
  G.className = 'grid' + (state.layout === 'grid' ? '' : ' tri');
  G.style.gridTemplateColumns = state.layout === 'grid' ? '' : `repeat(${2 * K - 1}, 1fr)`;
  $('layoutnote').textContent = {
    grid: '',
    tri: 'Row k holds every n with ⌈√n⌉ = k, from (k − 1)² + 1 to k², as in David Ellsworth\'s triangular table. Packings that beat the grid stop near the centre column, n ≈ k² − k; to its right nothing better than the grid is known.',
    deficit: 'Row k holds every n with ⌈√n⌉ = k, right-aligned, so each column is one deficit c = k² − n: the rightmost is the perfect squares, then k² − 1, k² − 2, … The five rightmost columns (c = 0 … 4) are proved grid-optimal for every k ≥ 6 (s(k² − c) = k; the families in the cards below).',
  }[state.layout];
  const counts = {};
  for (let n = 1; n <= NMAX; n++) {
    const t = tileInfo(n); counts[t.cat] = (counts[t.cat] || 0) + 1;
    // 135 of these 324 n have no catalogue entry at all. Linking them to explore.html?n=N made
    // gotoN snap silently to the nearest n that does, so the tile showed a different packing than
    // the one clicked; they are not links now, and the hover says why.
    const a = document.createElement(openable(t) ? 'a' : 'span');
    a.className = 'tile' + (Number.isInteger(Math.sqrt(n)) ? ' sq' : '') + (proved.has(n) ? ' proved' : '') + (openable(t) ? '' : ' nolink');
    a.style.background = colour(t); a.textContent = n;
    if (openable(t)) a.href = t.file ? `explore.html?p=${encodeURIComponent(t.file)}` : `explore.html?n=${n}`;
    a.onmouseenter = e => tip(e, t); a.onmousemove = e => move(e); a.onmouseleave = () => { $('tip').style.display = 'none'; };
    const p = place(n); if (p) { a.style.gridRow = p.row; a.style.gridColumn = p.col; }
    G.appendChild(a); inkFor(a);
  }
  legend();
  $('foot').textContent = `As of ${state.year === 2026 ? 'today' : state.year}: ${counts.trivial || 0} of ${NMAX} are best packed on the plain grid, ${counts.diagonal || 0} use only 45° tilts, ${counts.one || 0} have one other tilt angle and ${counts.many || 0} several. Dates are read from the catalogue's wording and rounded to the year.`;
}
function tip(e, t) {
  const T = $('tip'); const lines = [`n = ${t.n}`, `s = ${t.s.toFixed(t.file ? 6 : 0)}${t.file ? '…' : ' (grid)'}`];
  if (t.file) {
    const ev = (t.rec.events || []).filter(x => x.who).map(x => `${x.verb.toLowerCase()} by ${x.who} ${x.date.text}`);
    if (ev.length) lines.push(ev.slice(0, 2).join('; '));
    if (t.sum) lines.push(`${CATS[t.cat][0]} · ${t.sum.n_angles - 1} tilt angle${t.sum.n_angles === 2 ? '' : 's'} · ${t.sum.symmetry}`, t.sum.rigid ? 'rigid' : `${t.sum.free} free square${t.sum.free === 1 ? '' : 's'}`);
  } else if (proved.has(t.n)) lines.push('grid packing, proved optimal');
  else lines.push('no better packing than the grid is known');
  if (!openable(t)) lines.push('not in the catalogue — nothing to open');
  T.textContent = lines.join('\n'); T.style.display = 'block'; move(e);
}
function move(e) { const T = $('tip'); T.style.left = (e.clientX + 14) + 'px'; T.style.top = (e.clientY + 14) + 'px'; }

async function main() {
  [IDX, TL] = await Promise.all([fetch('data/index.json').then(r => r.json()), fetch('data/timeline.json').then(r => r.json())]);
  for (const [n, r] of Object.entries(IDX.records)) if (/Proved/.test(r.prose || '') || Number.isInteger(Math.sqrt(+n))) proved.add(+n);
  for (let n = 1; n <= NMAX; n++) if (Number.isInteger(Math.sqrt(n))) proved.add(n);
  try { const B = await fetch('data/lower_bounds.json').then(r => r.ok ? r.json() : null); if (B) for (const [n, b] of Object.entries(B)) { if (!b.best) continue; const t = tileInfo(+n); if (Math.abs(b.best.value - t.s) < 1e-9 && b.best.status === 'proved') proved.add(+n); } } catch (e) {}
  $('mode').onchange = () => { state.mode = $('mode').value; draw(); };
  const q = new URLSearchParams(location.search).get('layout'); if (['grid', 'tri', 'deficit'].includes(q)) { state.layout = q; $('layout').value = q; }
  $('layout').onchange = () => { state.layout = $('layout').value; history.replaceState(null, '', state.layout === 'grid' ? location.pathname : `?layout=${state.layout}`); draw(); };
  $('year').oninput = () => { state.year = +$('year').value; $('yearv').textContent = state.year === 2026 ? 'today' : state.year; draw(); };
  draw();
}
main();
})();
