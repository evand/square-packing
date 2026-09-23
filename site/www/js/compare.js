/* Compare: all packings for one n on a timeline strip, and a morph between two of them. */
(() => {
const $ = id => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
const el = (t, a = {}, p = null) => { const e = document.createElementNS(NS, t); for (const [k, v] of Object.entries(a)) e.setAttribute(k, v); if (p) p.appendChild(e); return e; };
let IDX, TL, n, list = [], A = null, B = null, match = null, timer = null, gen = 0;
const cache = {};
const MONTH = ['', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
const dateStr = d => d && d.date ? (d.date[1] && !d.events.every(e => e.date.approx) ? `${MONTH[d.date[1]]} ${d.date[0]}` : `${d.date[0]}`) : 'undated';
const who = d => { const e = (d.events || []).find(x => x.who); return e ? e.who : ''; };

async function get(file) { if (!cache[file]) cache[file] = fetch(`data/p/${encodeURIComponent(file.replace('.svg', '.json'))}`).then(r => r.ok ? r.json() : null).then(r => { if (r) r.s = +r.s; return r; }); return cache[file]; }

// A file we can actually draw: parsed, no hard errors, so data/p/<file>.json exists.
const drawable = f => !!(f && IDX.files[f] && IDX.files[f].n_parsed && !(IDX.files[f].errors || []).some(e => !e.startsWith('warn')));
// The timeline entry carrying a given file, wherever it is filed. Ellsworth's "2, 3" box links
// square-2.svg but embeds square-3.svg, which does not exist, so n = 3's only drawable packing is
// filed under n = 2 and the dates that go with it have to be fetched from there.
const entryFor = file => { for (const L of Object.values(TL)) { const d = L.find(x => x.file === file); if (d) return d; } return null; };
// Above this many squares the correspondence is not attempted -- see setup().
const MATCH_MAX = 1200;
const setMsg = t => { const e = $('msg'); if (!e) return; e.textContent = t || ''; e.style.display = t ? '' : 'none'; };

function thumbSVG(rec, size = 150) {
  const s = rec.s, svg = el('svg', { viewBox: `-0.02 -0.02 ${s + 0.04} ${s + 0.04}` });
  const w = el('g', { transform: `matrix(1 0 0 -1 0 ${s})` }, svg);
  el('rect', { x: 0, y: 0, width: s, height: s, fill: 'none', stroke: 'var(--wall)', 'stroke-width': s / 150 }, w);
  for (const [cx, cy, th] of rec.squares) el('rect', { x: -.5, y: -.5, width: 1, height: 1, transform: `translate(${cx} ${cy}) rotate(${th})`, fill: th % 90 < 1e-9 ? 'var(--sq-fill)' : `hsl(${((th > 45 ? th - 90 : th) * 4 + 360) % 360} 55% 62%)`, stroke: 'var(--sq-stroke)', 'stroke-width': s / 250 }, w);
  return svg;
}

async function loadN(nn, asked = nn) {
  // loadN awaits a fetch per packing, so two calls in flight at once (mashing the arrows, or a slow
  // network) let the older one finish last and publish A from its own list -- against which the
  // newer render has already run, leaving A.rec unset and setup() throwing on A.rec.s. Only the
  // newest call is allowed to touch shared state past an await.
  const my = ++gen;
  n = nn; $('nbox').value = n;
  // Ten n are catalogued only inside another n's entry and are absent from timeline.json, so the
  // old `TL[q.n] ? +q.n : 17` sent ?n=147 to n = 17 and rewrote the URL to match. Any snap says so.
  const snap = asked !== n ? `Nothing is catalogued for n = ${asked}; showing the nearest n that is, n = ${n}. ` : '';
  const all = TL[String(n)] || [];
  let L = all.filter(d => d.s && drawable(d.file));
  // Nothing filed under this n is drawable. The catalogue may still picture this n inside another
  // n's entry (n = 3), or record an s with no drawing at all (n = 585 and seven more, where
  // Ellsworth's prose says the SVG has not been made) -- the two cases must not look the same.
  let borrowed = null;
  if (!L.length) {
    const r = IDX.records[String(n)];
    if (r && drawable(r.svg)) { borrowed = r.svg; L = [entryFor(r.svg) || { file: r.svg, s: +IDX.files[r.svg].s, date: null, events: [] }]; }
  }
  L.sort((a, b) => ((a.date || [0, 0])[0] - (b.date || [0, 0])[0]) || ((a.date || [0, 0])[1] - (b.date || [0, 0])[1]) || (b.s - a.s));
  history.replaceState(null, '', `?n=${n}`);
  if (!L.length) {
    // "0 packings on record" was untrue as well as blank: the record exists, the drawing does not.
    list = []; A = B = match = null;
    $('count').textContent = all.length ? `${all.length} on record, not drawn` : 'nothing on record';
    $('strip').innerHTML = ''; $('view').innerHTML = ''; $('view').removeAttribute('viewBox');
    $('infoA').innerHTML = ''; $('infoB').innerHTML = ''; $('diff').innerHTML = '';
    const meta = all.length && IDX.files[all[0].file];
    setMsg(snap + (meta
      ? `Nothing to compare for n = ${n}. The catalogue records s = ${meta.s_dec || '?'} for it, from ${all[0].file}, but carries no drawing of that packing: ${(meta.prose || '').replace(/\s+/g, ' ').trim()}`
      : `Nothing to compare for n = ${n}: no packing of it is in the catalogue.`));
    return;
  }
  list = L; $('count').textContent = `${L.length} packing${L.length === 1 ? '' : 's'} on record`;
  if (borrowed) {
    const pn = IDX.files[borrowed].n_parsed;
    setMsg(snap + `n = ${n} is not drawn on its own. The catalogue covers n = ${Math.min(n, pn)} and n = ${Math.max(n, pn)} with a single entry, because s is the same for both; below is the packing that entry draws, which holds ${pn} square${pn === 1 ? '' : 's'}.`);
  } else setMsg(snap);
  const strip = $('strip'); strip.innerHTML = '';
  const recs = await Promise.all(L.map(d => get(d.file)));
  if (my !== gen) return;   // superseded while fetching
  // Drop anything whose analysis did not load before anything indexes this list: the thumbnails
  // carry their position in it, and setup() reads A/B back out of it by that index.
  L.forEach((d, i) => { if (recs[i]) d.rec = recs[i]; });
  L = L.filter(d => d.rec); list = L;
  if (!L.length) {
    A = B = match = null;
    $('view').innerHTML = ''; $('infoA').innerHTML = ''; $('infoB').innerHTML = ''; $('diff').innerHTML = '';
    $('count').textContent = `${all.length} on record, not loadable`;
    setMsg(snap + `The packings recorded for n = ${n} could not be loaded: the exported analysis is missing or failed to parse.`);
    return;
  }
  $('count').textContent = `${L.length} packing${L.length === 1 ? '' : 's'} on record`;
  L.forEach((d, i) => {
    const rec = d.rec;
    const t = document.createElement('div'); t.className = 'thumb'; t.dataset.i = i;
    t.appendChild(thumbSVG(rec));
    // `t.innerHTML +=` here serialised the SVG just appended back to markup and reparsed it -- at
    // n = 9465 that is two 9465-rect round trips per thumbnail, and it froze the tab for ~45 s.
    const cap = document.createElement('div');
    cap.innerHTML = `<div class="s">${rec.s.toFixed(6)}${d.is_record ? ' <span class="rec">record</span>' : ''}</div><div class="d">${dateStr(d)}${who(d) ? ' · ' + who(d) : ''}</div>`;
    while (cap.firstChild) t.appendChild(cap.firstChild);
    t.onclick = e => { if (e.shiftKey) B = d; else A = d; if (A === B) B = null; setup(); };
    strip.appendChild(t);
  });
  // default: oldest non-trivial vs record, or the two most recent
  // default pair: the previous record (best packing clearly worse than the current record) vs the record
  const rec = L.find(d => d.is_record) || L[L.length - 1];
  const key = d => d.date ? d.date[0] * 12 + (d.date[1] || 6) : -1;
  let worse = L.filter(d => d !== rec && d.s > rec.s + 1e-3 && d.date && key(d) <= key(rec));   // dated no later than the record
  if (!worse.length) worse = L.filter(d => d !== rec && d.s > rec.s + 1e-3);
  const prev = worse.length ? worse.reduce((b, d) => d.s < b.s ? d : b) : L.find(d => d !== rec);
  A = prev || rec; B = prev ? rec : null; setup();
}

// Hungarian assignment on a cost matrix (n x n), returns col for each row
function hungarian(cost) {
  const N = cost.length, INF = 1e18; const u = new Array(N + 1).fill(0), v = new Array(N + 1).fill(0), p = new Array(N + 1).fill(0), way = new Array(N + 1).fill(0);
  for (let i = 1; i <= N; i++) {
    p[0] = i; let j0 = 0; const minv = new Array(N + 1).fill(INF), used = new Array(N + 1).fill(false);
    do {
      used[j0] = true; const i0 = p[j0]; let delta = INF, j1 = 0;
      for (let j = 1; j <= N; j++) if (!used[j]) { const cur = cost[i0 - 1][j - 1] - u[i0] - v[j]; if (cur < minv[j]) { minv[j] = cur; way[j] = j0; } if (minv[j] < delta) { delta = minv[j]; j1 = j; } }
      for (let j = 0; j <= N; j++) if (used[j]) { u[p[j]] += delta; v[j] -= delta; } else minv[j] -= delta;
      j0 = j1;
    } while (p[j0] !== 0);
    do { const j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0);
  }
  const res = new Array(N); for (let j = 1; j <= N; j++) if (p[j]) res[p[j] - 1] = j - 1; return res;
}
const angDiff = (a, b) => { let d = Math.abs(a - b) % 90; return Math.min(d, 90 - d); };

function setup() {
  document.querySelectorAll('.thumb').forEach(t => { const d = list[+t.dataset.i]; t.classList.toggle('a', d === A); t.classList.toggle('b', d === B); });
  if (!A || !A.rec) return;
  if (B && !B.rec) B = null;
  const ra = A.rec, rb = B && B.rec;
  $('infoA').innerHTML = `<div class="s">${ra.s.toFixed(8)}</div>${dateStr(A)}${who(A) ? '<br>' + who(A) : ''}<br><a href="explore.html?p=${encodeURIComponent(A.file)}">explore →</a>`;
  $('infoB').innerHTML = rb ? `<div class="s">${rb.s.toFixed(8)}</div>${dateStr(B)}${who(B) ? '<br>' + who(B) : ''}<br><a href="explore.html?p=${encodeURIComponent(B.file)}">explore →</a>` : '<span class="note">shift-click a packing to compare</span>';
  match = null; let why = '';
  // The correspondence costs an N x N matrix and an assignment, solved once per symmetry of the
  // square. Measured on the shipped data: n = 1037 takes 0.7 s for all eight, n = 2135 7.7 s, and
  // n = 9465 29 s with a 600 MB matrix -- which hangs, then kills, the tab. The n box accepts
  // anything in timeline.json, so the ceiling has to be here rather than on the input.
  if (rb && Math.max(ra.squares.length, rb.squares.length) > MATCH_MAX) {
    why = `Both packings are measured below, but they are too large to match square-for-square: that costs a ${MATCH_MAX}+ square assignment solved eight times over, which would lock the page up for tens of seconds. The morph is off above n = ${MATCH_MAX}.`;
  } else if (rb) {
    // normalise both to the same scale (unit squares; container sides differ slightly) and match centres; try the 8 symmetries of the square and keep the best
    const sa = ra.s, sb = rb.s; let best = null;
    const syms = [p => p, p => [sb - p[0], p[1], -p[2]], p => [p[0], sb - p[1], -p[2]], p => [sb - p[0], sb - p[1], p[2]], p => [p[1], p[0], 90 - p[2]], p => [sb - p[1], p[0], 90 + p[2]], p => [p[1], sb - p[0], 90 + p[2]], p => [sb - p[1], sb - p[0], 90 - p[2]]];
    for (let k = 0; k < syms.length; k++) {
      const qb = rb.squares.map(q => { const r = syms[k](q); return [r[0], r[1], ((r[2] % 90) + 90) % 90]; });
      const N = Math.max(ra.squares.length, qb.length); const cost = [];
      for (let i = 0; i < N; i++) { cost.push([]); for (let j = 0; j < N; j++) { const a = ra.squares[i], b = qb[j]; cost[i].push(a && b ? Math.hypot(a[0] - b[0], a[1] - b[1]) + angDiff(a[2], b[2]) / 60 : 3); } }
      const asg = hungarian(cost); const tot = asg.reduce((t, j, i) => t + cost[i][j], 0);
      if (!best || tot < best.tot) best = { tot, asg, qb, k };
    }
    match = best;
  }
  $('t').value = 0; draw();
  // No correspondence means no morph: leave the controls off rather than sliding the container
  // walls in and out around squares that cannot follow them.
  const canMorph = !!(rb && match);
  $('t').disabled = !canMorph; $('play').disabled = !canMorph;
  if (!canMorph && timer) { clearInterval(timer); timer = null; $('play').textContent = 'play'; }
  const moved = match ? match.asg.filter((j, i) => ra.squares[i] && match.qb[j] && (Math.hypot(ra.squares[i][0] - match.qb[j][0], ra.squares[i][1] - match.qb[j][1]) > 0.05 || angDiff(ra.squares[i][2], match.qb[j][2]) > 0.05)).length : 0;
  // Two packings can tie -- n = 9465's pair are alternatives with the same s -- and rendering that
  // as "+0.000000" reads as a change of nothing rather than as no change. Say which it is. The
  // exported values are doubles (see TODO, stage 1), so an exact tie is only claimed as such.
  const ds = rb ? ra.s - rb.s : 0, tied = rb && Math.abs(ds) < 5e-7;   // rb is null until a B is picked
  $('diff').innerHTML = !rb ? '' :
    (tied
      ? `<div>Both sides: <b>${ra.s.toFixed(6)}</b> — ${ds === 0 ? 'the same s in the exported values; these are alternative packings of one record' : 'equal to the 6 decimals shown, differing below that'}.</div>`
      : `<div>Side ${ra.s.toFixed(6)} → ${rb.s.toFixed(6)}: <b>${ds > 0 ? '−' : '+'}${Math.abs(ds).toFixed(6)}</b> (${(100 * Math.abs(ds) / ra.s).toFixed(3)}%)</div>`) +
    (match ? `<div>${moved} of ${ra.squares.length} squares move or turn by more than 0.05${match.k ? ' (after reflecting/rotating the newer packing to line up best)' : ''}</div>`
           : `<div class="note">${why}</div>`) +
    `<div>Distinct angles: ${ra.analysis.n_angles} → ${rb.analysis.n_angles}; free squares ${ra.analysis.free.length} → ${rb.analysis.free.length}</div>`;
}

function stopPlay() { if (timer) { clearInterval(timer); timer = null; $('play').textContent = 'play'; } }
function draw() {
  if (!A) { stopPlay(); return; }
  const svg = $('view'); svg.innerHTML = ''; const t = +$('t').value / 1000;
  const ra = A.rec, rb = B && B.rec; const s = (rb && match) ? ra.s + (rb.s - ra.s) * t : ra.s;
  svg.setAttribute('viewBox', `${-0.03 * s} ${-0.03 * s} ${1.06 * s} ${1.06 * s}`);
  const w = el('g', { transform: `matrix(1 0 0 -1 0 ${s})` }, svg);
  el('rect', { x: 0, y: 0, width: s, height: s, fill: 'none', stroke: 'var(--wall)', 'stroke-width': s / 300 }, w);
  ra.squares.forEach((q, i) => {
    let [cx, cy, th] = q, moved = false;
    if (rb && match) {
      const qb = match.qb[match.asg[i]];
      if (qb) {
        let dth = qb[2] - th; if (dth > 45) dth -= 90; if (dth < -45) dth += 90;
        moved = Math.hypot(qb[0] - cx, qb[1] - cy) > 0.05 || Math.abs(dth) > 0.05;
        cx += (qb[0] - cx) * t; cy += (qb[1] - cy) * t; th += dth * t;
      }
    }
    const tilt = ((th % 90) + 90) % 90; const tl = tilt > 45 ? tilt - 90 : tilt;
    el('rect', { x: -.5, y: -.5, width: 1, height: 1, transform: `translate(${cx} ${cy}) rotate(${th})`, fill: Math.abs(tl) < 1e-9 ? 'var(--sq-fill)' : `hsl(${(tl * 4 + 360) % 360} 55% 62%)`,
      stroke: moved ? 'var(--red)' : 'var(--sq-stroke)', 'stroke-width': moved ? s / 200 : s / 500 }, w);
  });
  if (rb && match && rb.squares.length > ra.squares.length) {
    // extra squares in B fade in
    match.qb.forEach((qb, j) => { if (match.asg.includes(j)) return; el('rect', { x: -.5, y: -.5, width: 1, height: 1, transform: `translate(${qb[0]} ${qb[1]}) rotate(${qb[2]})`, fill: 'var(--gold-soft)', stroke: 'var(--gold)', 'stroke-width': s / 300, opacity: t }, w); });
  }
}

async function main() {
  [IDX, TL] = await Promise.all([fetch('data/index.json').then(r => r.json()), fetch('data/timeline.json').then(r => r.json())]);
  const q = Object.fromEntries(new URLSearchParams(location.search));
  // Every n the page has something to say about: one filed under it, or a record entry pointing at
  // a packing filed elsewhere. timeline.json alone misses ten of them.
  const ns = [...new Set([...Object.keys(TL), ...Object.keys(IDX.records)].map(Number))].sort((a, b) => a - b);
  const nearest = v => ns.reduce((b, x) => Math.abs(x - v) < Math.abs(b - v) ? x : b, ns[0]);
  // A non-numeric n box must not leave a blank page: fall back to whatever is loaded, else 17.
  const go = v => { if (!Number.isFinite(v)) { $('nbox').value = n != null ? n : 17; return n != null ? null : loadN(17); }
                    return loadN(ns.includes(v) ? v : nearest(v), v); };
  $('nbox').onchange = () => go(+$('nbox').value);
  $('prev').onclick = () => { const c = ns.filter(x => x < n); if (c.length) loadN(c[c.length - 1]); };
  $('next').onclick = () => { const c = ns.filter(x => x > n); if (c.length) loadN(c[0]); };
  $('t').oninput = draw;
  $('swap').onclick = () => { if (B) { [A, B] = [B, A]; setup(); } };
  $('play').onclick = () => { if (timer) { clearInterval(timer); timer = null; $('play').textContent = 'play'; return; } let dir = 1; $('play').textContent = 'stop'; timer = setInterval(() => { let v = +$('t').value + 12 * dir; if (v >= 1000) { v = 1000; dir = -1; } if (v <= 0) { v = 0; dir = 1; } $('t').value = v; draw(); }, 40); };
  go(q.n != null && q.n !== '' ? +q.n : 17);
}
main();
})();
