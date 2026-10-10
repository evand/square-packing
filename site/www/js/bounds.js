/* Bounds page: gap chart (dumbbell per n), one-n timeline, table. Data: lower_bounds.json, timeline.json, index.json */
(() => {
const $ = id => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
const el = (t, a = {}, p = null) => { const e = document.createElementNS(NS, t); for (const [k, v] of Object.entries(a)) e.setAttribute(k, v); if (p) p.appendChild(e); return e; };
let LB, TL, IDX;
const trivial = n => Math.ceil(Math.sqrt(n));
const who = d => { const e = (d.events || []).find(x => x.who); return e ? e.who : ''; };

function upper(n) {   // {s, file, date, who}
  const L = (TL[String(n)] || []).filter(d => d.s && !(IDX.files[d.file] && (IDX.files[d.file].start || (IDX.files[d.file].errors || []).some(e => !e.startsWith('warn:')))));
  const rec = L.find(d => d.is_record) || (L.length ? L.reduce((b, d) => d.s < b.s ? d : b) : null);
  const t = trivial(n);
  if (!rec || rec.s >= t - 1e-9) return { s: t, trivial: true, who: '', date: null };
  const R = IDX.records[String(n)], pending = !!(R && R.credit && !R.credit.registered && rec.origin);
  return { s: rec.s, trivial: false, who: who(rec), date: rec.date, file: rec.file, pending };
}
function lower(n) {
  const b = LB && LB[String(n)] && LB[String(n)].best;
  if (b) return b;
  const t = Math.sqrt(n); return { value: t, source: 'area', year: null, status: 'proved' };
}
const isSettled = (u, l) => Math.abs(u.s - l.value) < 1e-9;
const SHORT = { 'area-bound': 'area', 'trivial-square': 'grid = area', 'Nagamochi2005': 'Nagamochi 2005', 'FriedmanDS7': 'Friedman', 'Gobel1979': 'Göbel 1979', 'Stromquist2003': 'Stromquist 2003', 'Green2000': 'Green 2000', 'Bentz2010': 'Bentz 2010', 'Bentz2016': 'Bentz 2016', 'KearneyShiu2002': 'Kearney–Shiu 2002', 'ElMoumni1999': 'El Moumni 1999', 'EvanDaniel2026': 'evand 2026', 'MiraAcc2026': 'Mira 2026', 'Burns2026': 'Burns 2026', 'Massaccesi2026': 'Massaccesi 2026', 'Fort2026': 'Fort 2026', 'jlevy2026': 'jlevy 2026', 'wand125_2026': 'wand125 2026', 'Kleddamag2026': 'Kleddamag 2026', 'Guzhou0806_2026': 'Guzhou0806 2026', 'tokoharu2026': 'tokoharu 2026', 'chelokot2026': 'chelokot 2026', 'Karakus2026': 'Karakuş 2026', 'Queuingtheorydotcom2026': 'Queuingtheorydotcom 2026', 'squarepacker2026': 'squarepacker 2026' };
const srcLabel = h => { if (!h || !h.source) return ''; if (SHORT[h.source]) return SHORT[h.source]; const S = LB && LB.sources && LB.sources[h.source]; return S ? `${S.authors} ${S.year || ''}`.trim() : h.source; };

// ---------- chart 1: gap per n ----------
function chart1() {
  const box = $('c1'); box.querySelectorAll('svg').forEach(s => s.remove());
  const nmax = Math.max(10, Math.min(100, +$('nmax').value || 60)), rel = $('rel').checked, drop = $('dropTrivial').checked;
  const rows = [];
  for (let n = 2; n <= nmax; n++) {
    const u = upper(n), l = lower(n);
    if (drop && u.trivial && isSettled(u, l)) continue;
    const t = trivial(n);
    rows.push({ n, u, l, settled: isSettled(u, l), yu: rel ? 0 : u.s - t, yl: rel ? -100 * (u.s - l.value) / u.s : l.value - t });
  }
  const W = 1000, H = 360, ml = 48, mr = 12, mt = 14, mb = 34;
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}` }, box);
  const xs = i => ml + (i + 0.5) / rows.length * (W - ml - mr);
  const ymin = rel ? Math.min(-1, ...rows.map(r => r.yl)) - 0.3 : Math.min(...rows.map(r => r.yl), -0.5) - 0.03;   // padding: the lowest dots sat on the edge
  const ymax = rel ? 0.5 : 0.05;
  const ys = v => mt + (ymax - v) / (ymax - ymin) * (H - mt - mb);
  // gridlines + axis labels
  const g = el('g', { class: 'grid' }, svg), ax = el('g', { class: 'axis' }, svg);
  const step = rel ? 2 : 0.1;
  for (let v = Math.ceil(ymin / step) * step; v <= ymax + 1e-9; v += step) {
    el('line', { x1: ml, x2: W - mr, y1: ys(v), y2: ys(v) }, g);
    const t = el('text', { x: ml - 6, y: ys(v) + 4, 'text-anchor': 'end' }, ax); t.textContent = rel ? `${v.toFixed(0)}%` : (Math.abs(v) < 1e-9 ? '⌈√n⌉' : v.toFixed(1));
  }
  el('line', { x1: ml, x2: W - mr, y1: ys(0), y2: ys(0), stroke: 'var(--muted)', 'stroke-width': 1.2 }, svg);
  const lab = el('text', { x: ml, y: H - 8, class: 'axis' }, ax); lab.textContent = rel ? 'how far below the best packing the floor is (%)' : 'side minus the grid side ⌈√n⌉  (0 = the plain grid)';
  const tt = box.querySelector('.tt');
  rows.forEach((r, i) => {
    const x = xs(i), y1 = ys(r.yu), y2 = ys(r.yl);
    const grp = el('g', { style: 'cursor:pointer' }, svg);
    el('rect', { x: x - 6, y: mt, width: 12, height: H - mt - mb, fill: 'transparent' }, grp);   // hit target
    if (!r.settled) el('line', { x1: x, x2: x, y1, y2, stroke: 'var(--rule)', 'stroke-width': 4, 'stroke-linecap': 'round' }, grp);
    if (r.settled) el('circle', { cx: x, cy: y1, r: 4.5, fill: 'var(--ink)', stroke: 'var(--surface)', 'stroke-width': 1.5 }, grp);
    else { el('circle', { cx: x, cy: y2, r: 4, fill: 'var(--chart-b)', stroke: 'var(--surface)', 'stroke-width': 1.5 }, grp); el('circle', { cx: x, cy: y1, r: 4, fill: 'var(--chart-a)', stroke: 'var(--surface)', 'stroke-width': 1.5 }, grp); }
    if (rows.length <= 70 || r.n % 5 === 0) { const t = el('text', { x, y: H - mb + 14, 'text-anchor': 'middle', class: 'axis', style: 'font-size:10px;fill:var(--muted)' }, ax); t.textContent = r.n; }
    grp.onmouseenter = () => { tt.textContent = tip(r); tt.style.display = 'block'; };
    grp.onmousemove = e => { const bb = box.getBoundingClientRect(); tt.style.left = Math.min(e.clientX - bb.left + 12, bb.width - 260) + 'px'; tt.style.top = (e.clientY - bb.top + 12) + 'px'; };
    grp.onmouseleave = () => { tt.style.display = 'none'; };
    grp.onclick = () => { $('npick').value = r.n; chart2(); $('c2').scrollIntoView({ behavior: 'smooth', block: 'center' }); };
  });
  const settled = rows.filter(r => r.settled).length, open = rows.length - settled;
  $('c1note').textContent = `${settled} of the ${rows.length} shown are settled; ${open} are open. Hover for details; click to see that n's history below, with a link to its packing.`;
}
function tip(r) {
  const lines = [`n = ${r.n}`, `best packing  ${r.u.s.toFixed(6)}${r.u.trivial ? ' (grid)' : ''}${r.u.who ? '  ' + r.u.who : ''}${r.u.date ? ' ' + r.u.date[0] : ''}`,
    `proven floor  ${r.l.value.toFixed(6)}  ${srcLabel(r.l)}${r.l.status === 'preprint' ? ' (unrefereed)' : ''}`];
  lines.push(r.settled ? 'settled' : `gap ${(r.u.s - r.l.value).toFixed(4)} (${(100 * (r.u.s - r.l.value) / r.u.s).toFixed(2)}%)`);
  return lines.join('\n');
}

// ---------- chart 2: one n over time ----------
function chart2() {
  const n = Math.max(2, Math.min(100, +$('npick').value || 17)); const box = $('c2'); box.querySelectorAll('svg').forEach(s => s.remove());
  const t = trivial(n);
  // upper-bound steps: dated packings, best-so-far
  // An undated record goes at the left edge (and says so) rather than being dropped: n = 5's record
  // is undated, and dropping it had this chart calling a proved value open while the table said settled.
  // Not packings: optimiser starts and files with real errors (51_r2__invalid overlaps), as on Overview.
  const usable = d => { const f = IDX.files[d.file]; return !f || !(f.start || (f.errors || []).some(e => !e.startsWith('warn:'))); };
  const packs = (TL[String(n)] || []).filter(d => d.s && (d.date || d.is_record) && usable(d)).map(d => ({ y: d.date ? d.date[0] + ((d.date[1] || 6) - 0.5) / 12 : 1979, undated: !d.date, s: d.s, who: who(d), file: d.file, d })).sort((a, b) => a.y - b.y);
  const up = [{ y: 1979, s: t, who: 'grid', label: 'plain grid' }]; let best = t;
  for (const p of packs) if (p.s < best - 1e-9) { best = p.s; up.push(p); }
  // lower-bound steps
  // `best` is normally repeated in `history`; merge it in anyway so a data slip cannot hide the current floor.
  const E = LB && LB[String(n)], hist = ((E && E.history) || []).slice();
  if (E && E.best && !hist.some(h => Math.abs(h.value - E.best.value) < 1e-12 && h.source === E.best.source && !h.reported)) hist.push(E.best);
  const lo = [{ y: 1979, v: Math.sqrt(n), src: 'area bound', when: '' }]; let bl = Math.sqrt(n);
  for (const h of hist.filter(h => !h.reported && h.status !== 'claimed' && h.status !== 'gap').sort((a, b) => hy(a) - hy(b))) if (h.value > bl + 1e-9) { bl = h.value; lo.push({ y: Math.max(1979, hy(h)), v: h.value, src: srcLabel(h), status: h.status, when: h.date || String(h.year || '') }); }
  // "since August 2026" zooms to the 2026 rush, when most floors moved; steps from before the
  // window start at its left edge.  The right edge follows the newest dated floor (plus a week),
  // so new results never fall off the end of the chart.
  const newest = Math.max(2026 + 9 / 12, ...Object.keys(LB || {}).filter(k => /^\d+$/.test(k)).flatMap(k => (LB[k].history || []).concat(LB[k].best || [])).map(hy));
  const zoom = $('zoom').checked, y0 = zoom ? 2026 + 7 / 12 : 1978, y1 = zoom ? newest + 7 / 365 : 2027;
  const inWin = (pts, key) => { const k = pts.findLastIndex(p => p.y < y0); return pts.filter((p, i) => i >= k).map(p => p[key]); };
  const allv = zoom ? [...inWin(up, 's'), ...inWin(lo, 'v')] : [...up.map(u => u.s), ...lo.map(l => l.v)]; let vmin = Math.min(...allv), vmax = Math.max(...allv); const pad = Math.max(0.02, (vmax - vmin) * 0.15); vmin -= pad; vmax += pad;
  const W = 1000, H = 300, ml = 62, mr = 16, mt = 14, mb = 30;
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}` }, box);
  const xs = y => ml + (y - y0) / (y1 - y0) * (W - ml - mr), ys = v => mt + (vmax - v) / (vmax - vmin) * (H - mt - mb);
  const g = el('g', { class: 'grid' }, svg), ax = el('g', { class: 'axis' }, svg);
  const MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  const ticks = zoom ? Array.from({ length: 2 * (Math.floor((y1 - 2026) * 12) - 7 + 1) }, (_, i) => [7 + (i >> 1), i % 2 ? 15 : 1]).map(([m, dd]) => [`${MON[m]} ${dd}`, 2026 + (m + (dd - 1) / 31) / 12]).filter(t => t[1] <= y1)
                     : Array.from({ length: 10 }, (_, i) => [String(1980 + 5 * i), 1980 + 5 * i]);
  for (const [t, y] of ticks) { el('line', { x1: xs(y), x2: xs(y), y1: mt, y2: H - mb }, g); const tx = el('text', { x: xs(y), y: H - mb + 16, 'text-anchor': 'middle' }, ax); tx.textContent = t; }
  const vstep = niceStep((vmax - vmin) / 5);
  for (let v = Math.ceil(vmin / vstep) * vstep; v <= vmax; v += vstep) { el('line', { x1: ml, x2: W - mr, y1: ys(v), y2: ys(v) }, g); const tx = el('text', { x: ml - 6, y: ys(v) + 4, 'text-anchor': 'end' }, ax); tx.textContent = v.toFixed(vstep < 0.01 ? 3 : 2); }
  const stepPath = (pts, key) => { let d = ''; pts.forEach((p, i) => { const x = xs(Math.max(p.y, y0)), y = ys(p[key]); d += i === 0 ? `M${x},${y}` : `H${x}V${y}`; }); d += `H${xs(y1)}`; return d; };
  el('path', { d: stepPath(up, 's'), fill: 'none', stroke: 'var(--chart-a)', 'stroke-width': 2 }, svg);
  el('path', { d: stepPath(lo, 'v'), fill: 'none', stroke: 'var(--chart-b)', 'stroke-width': 2 }, svg);
  const tt = box.querySelector('.tt');
  // An upper-bound dot is a packing: clicking it opens that packing in Explore.
  const dot = (x, y, col, text, href) => { const c = el('circle', { cx: x, cy: y, r: 5, fill: col, stroke: 'var(--surface)', 'stroke-width': 1.5, style: 'cursor:pointer' }, svg); const hit = el('circle', { cx: x, cy: y, r: 12, fill: 'transparent', style: 'cursor:pointer' }, svg); hit.onmouseenter = c.onmouseenter = () => { tt.textContent = text + (href ? '\nclick to open it in Explore' : ''); tt.style.display = 'block'; tt.style.left = Math.min(x / W * box.clientWidth + 12, box.clientWidth - 240) + 'px'; tt.style.top = (y / H * box.clientWidth * H / W - 10) + 'px'; }; hit.onmouseleave = c.onmouseleave = () => { tt.style.display = 'none'; }; if (href) hit.onclick = c.onclick = () => { location.href = href; }; };
  const R = IDX.records[String(n)] || {};
  const packHref = p => p.file ? (p.file === R.best ? `explore.html?n=${n}` : `explore.html?p=${encodeURIComponent(p.file)}`) : p.label === 'plain grid' && R.best_s != null && Math.abs(+R.best_s - t) < 1e-9 ? `explore.html?n=${n}` : null;   // the grid is drawn only where it is the record
  const shown = p => p.y >= y0;
  up.filter(shown).forEach(p => dot(xs(p.y), ys(p.s), 'var(--chart-a)', `${p.s.toFixed(6)}\n${p.label || p.who || ''} ${p.undated ? '(date not recorded)' : p.d && p.d.origin ? p.d.events[0].date.text : Math.floor(p.y)}${p.d && p.d.origin && R.credit && !R.credit.registered && p.file === R.best ? '\nnot yet on the register' : ''}`, packHref(p)));
  lo.filter(shown).forEach(p => dot(xs(p.y), ys(p.v), 'var(--chart-b)', `≥ ${p.v.toFixed(6)}\n${p.src || ''}${p.when && !String(p.src || '').includes(p.when) ? ' · ' + p.when : ''}${p.status === 'preprint' ? '\n(unrefereed)' : ''}`));
  // direct labels at the right end
  const lu = el('text', { x: W - mr - 4, y: ys(up[up.length - 1].s) - 7, 'text-anchor': 'end', style: 'font-size:11px;fill:var(--muted)' }, svg); lu.textContent = `best ${up[up.length - 1].s.toFixed(5)}`;
  const ll = el('text', { x: W - mr - 4, y: ys(lo[lo.length - 1].v) + 15, 'text-anchor': 'end', style: 'font-size:11px;fill:var(--muted)' }, svg); ll.textContent = `floor ${lo[lo.length - 1].v.toFixed(5)}`;
  const settled = Math.abs(up[up.length - 1].s - lo[lo.length - 1].v) < 1e-9;
  $('npickhint').textContent = settled ? `s(${n}) is settled` : `gap ${(up[up.length - 1].s - lo[lo.length - 1].v).toFixed(5)}`;
  // story cards
  const st = $('story'); st.innerHTML = '';
  const cards = [];
  const first = up[1]; if (first) cards.push(`First improvement on the grid: <b>${first.s.toFixed(5)}</b> by ${first.who || '?'}${first.undated ? ' (date not recorded)' : ' in ' + Math.floor(first.y)}.`);
  const last = up[up.length - 1]; if (up.length > 2) cards.push(`Current record <b>${last.s.toFixed(6)}</b> (${last.who || ''}${last.undated ? '' : ', ' + Math.floor(last.y)}) after ${up.length - 1} improvements.`);
  const ll2 = lo[lo.length - 1]; cards.push(`Proven floor <b>${ll2.v.toFixed(6)}</b> (${ll2.src || 'area'}${ll2.when && ll2.when.length > 4 ? ', ' + ll2.when : ''}${ll2.status === 'preprint' ? ', unrefereed' : ''}).`);
  if (!settled) cards.push(`Open by <b>${(last.s - ll2.v).toFixed(5)}</b> — ${(100 * (last.s - ll2.v) / last.s).toFixed(2)}% of the side.`);
  for (const c of cards) { const d = document.createElement('div'); d.innerHTML = c; st.appendChild(d); }
  const xb = $('npickexplore'); if (xb) { xb.href = `explore.html?n=${n}`; xb.textContent = `see the best packing for n = ${n} in Explore →`; }
  histTable(n, hist);
  history.replaceState(null, '', `?n=${n}`);
}
// Fractional year of a floor: its ISO date if it has one ("2026-09-04" or "2026-09"), else mid-year.
function hy(h) { const m = /^(\d{4})(?:-(\d{2}))?(?:-(\d{2}))?/.exec(h.date || ''); if (!m) return h.year || 1979; return +m[1] + ((m[2] ? +m[2] - 1 : 5.5) + (m[3] ? (+m[3] - 0.5) / 31 : 0.5)) / 12; }
// Every recorded floor for n, oldest first, with its source and note: the same data as the chart.
function histTable(n, hist) {
  const tb = $('histtbl').querySelector('tbody'); tb.innerHTML = '';
  $('histn').textContent = n;
  const esc = t => String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;');
  const rows = hist.filter(h => h.source !== 'area-bound').sort((a, b) => hy(a) - hy(b));
  for (const h of rows) {
    const S = LB.sources && LB.sources[h.source], url = h.url || (S && S.url);
    const src = url ? `<a href="${esc(url)}">${esc(srcLabel(h))}</a>` : esc(srcLabel(h));
    const tr = document.createElement('tr');
    tr.innerHTML = `<td class="num">${esc(h.date || h.year || '')}</td><td class="num">${h.value.toFixed(6)}</td><td>${src}${h.reported ? ' <span class="tag open">reported, not replayed</span>' : h.status === 'claimed' ? ' <span class="tag open">claimed, no published proof</span>' : h.status === 'gap' ? ' <span class="tag open">published proof incomplete</span>' : h.status === 'preprint' ? ' <span class="tag open">unrefereed</span>' : ''}</td><td>${esc(h.note || '')}</td>`;
    tb.appendChild(tr);
  }
  if (!rows.length) tb.innerHTML = '<tr><td colspan="4">Only the area bound √n is recorded for this n.</td></tr>';
}
function niceStep(x) { const p = Math.pow(10, Math.floor(Math.log10(x))); const m = x / p; return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; }

// ---------- table ----------
function table() {
  const tb = $('tbl').querySelector('tbody'); tb.innerHTML = '';
  for (let n = 1; n <= 100; n++) {
    const u = upper(n), l = lower(n), settled = isSettled(u, l);
    const tr = document.createElement('tr');
    tr.innerHTML = `<td class="num"><a href="explore.html?n=${n}" title="explore n = ${n}">${n}</a></td><td class="num">${u.s.toFixed(6)}${u.trivial ? ' (grid)' : ''}</td><td>${u.who}${u.date ? ' ' + u.date[0] : ''}${u.pending ? ' <span class="tag open" title="posted as a registration request on jlevy/squares; checked by our two exact certificate checkers">not yet registered</span>' : ''}</td><td class="num">${l.value.toFixed(6)}</td><td>${srcLabel(l)}${l.status === 'preprint' ? ' <span class="tag open">unrefereed</span>' : ''}</td><td class="num">${settled ? '—' : (u.s - l.value).toFixed(4)}</td><td>${settled ? '<span class="tag ok">settled</span>' : '<span class="tag open">open</span>'}</td>`;
    tb.appendChild(tr);
  }
}

async function main() {
  [TL, IDX] = await Promise.all([fetch('data/timeline.json').then(r => r.json()), fetch('data/index.json').then(r => r.json())]);
  try { LB = await fetch('data/lower_bounds.json').then(r => r.ok ? r.json() : null); } catch (e) { LB = null; }
  const q = Object.fromEntries(new URLSearchParams(location.search)); if (q.n) $('npick').value = q.n;
  for (const id of ['nmax', 'rel', 'dropTrivial']) $(id).onchange = chart1;
  $('npick').onchange = chart2; $('zoom').onchange = chart2;
  chart1(); chart2(); table();
}
main();
})();
