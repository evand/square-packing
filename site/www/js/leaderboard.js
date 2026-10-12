/* Leaderboard: map of record ages, per-n credit chain, hall of fame, oldest / newest records.
   Data: data/credits.json (site/tools/credits.py), data/hist.json + data/p/*.json (a drawing per holder;
   site/tools/history_packings.py). */
(() => {
const $ = id => document.getElementById(id);
const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
let C, P, TODAY, HIST = {};
const NMAX = 324;

const dnum = d => d ? d[0] + ((d[1] || 6.5) - 0.5) / 12 : null;           // (year, month) -> fractional year
const dfmt = d => !d ? '' : d[1] ? `${['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'][d[1] - 1]} ${d[0]}` : String(d[0]);
const name = id => P[id] ? P[id].name : id;
const nameH = id => { const p = P[id]; if (!p) return esc(id); const h = p.handle && p.handle !== p.name ? ` <span class="handle">${esc(p.handle)}</span>` : ''; return esc(p.name) + h; };
const names = ids => (ids || []).map(nameH).join(', ');
const plain = t => esc(String(t || '').replace(/\$s\((\d+)\)\$/g, 's($1)').replace(/\$([^$]*)\$/g, '$1'));
const upper = n => C.n[n].upper;
const current = n => { const u = upper(n); return u[u.length - 1]; };
const lowerBest = n => { const l = C.n[n].lower; return l && l.length ? l[l.length - 1] : null; };
const isSquare = n => Math.round(Math.sqrt(n)) ** 2 === n;
function settled(n) {
  if (isSquare(n)) return true;
  const l = lowerBest(n), h = current(n);
  return !!(l && h.s && Math.abs(l.value - h.s) < 1e-9);
}
const credited = h => h.kind === 'construction' ? (h.origin || []) : (h.who || []);
function holderLabel(h) {
  if (h.kind === 'grid') return 'grid packing';
  if (h.kind === 'construction') return `construction on ${names(h.origin)}'s packing`;
  return names(h.who);
}

// ---------- age bins (sequential, one hue; darker = older) ----------
const BINS = [  // [label, test on year]
  ['2026', y => y >= 2026], ['2020–25', y => y >= 2020], ['2000–19', y => y >= 2000], ['1981–99', y => y >= 1981], ['1979–80', y => true]];
const bin = y => BINS.findIndex(b => b[1](y));

function scale() {
  $('scale').innerHTML = 'record set in: ' + BINS.map((b, i) => `<span><i style="background:var(--seq${i})"></i>${b[0]}</span>`).join('') +
    '<span><i style="box-shadow:inset 0 0 0 1px var(--rule)"></i>grid packing</span>';
}

// ---------- tallies (mirror credits.py held()) ----------
function tallies() {
  const T = {};
  const t = id => T[id] || (T[id] = { id, now: [], ever: new Set(), years: 0, built: [], lbNow: [], lbEver: new Set(), first: Infinity });
  for (let n = 1; n <= NMAX; n++) {
    const hs = upper(n);
    hs.forEach((h, i) => {
      const last = i === hs.length - 1;
      if (h.kind === 'construction') { if (last) (h.origin || []).forEach(w => t(w).built.push(n)); return; }
      if (!h.who || !h.who.length || !h.date) return;
      const y = Math.max(0, dnum(h.until || TODAY) - dnum(h.date)) / h.who.length;
      h.who.forEach(w => { const r = t(w); r.ever.add(n); r.years += y; r.first = Math.min(r.first, h.date[0]); if (last) r.now.push(n); });
    });
    for (const c of C.n[n].other || []) (c.who || []).forEach(w => { if (c.date) t(w).first = Math.min(t(w).first, c.date[0]); });
    const lc = C.n[n].lower || [];
    lc.forEach((l, i) => l.who.forEach(w => { const r = t(w); r.lbEver.add(n); r.first = Math.min(r.first, l.year || Infinity); if (i === lc.length - 1) r.lbNow.push(n); }));
  }
  return Object.values(T).map(r => ({ ...r, ever: [...r.ever].sort((a, b) => a - b), lbEver: [...r.lbEver].sort((a, b) => a - b),
    nowN: r.now.length, everN: r.ever.size, builtN: r.built.length, lbNowN: r.lbNow.length, lbEverN: r.lbEver.size }));
}

// ---------- map ----------
let selN = 11, hi = '', layout = 'deficit';   // the Overview's default too
const K = 18;
// 'tri': row k holds n = (k-1)^2+1 .. k^2, centred (Ellsworth's triangular table); 'deficit': the same rows
// right-aligned, one column per k^2 - n; 'grid': 18 per row.  Same placement as the Overview's grid.
function place(n) {
  if (layout === 'grid') return null;
  const k = Math.ceil(Math.sqrt(n) - 1e-12);
  return { row: k, col: layout === 'tri' ? (K - k) + (n - (k - 1) * (k - 1)) : 2 * K - 1 - (k * k - n) };
}
const LAYOUTNOTE = {
  tri: 'One cell per n; row k holds every n with ⌈√n⌉ = k, from (k − 1)² + 1 to k², as in David Ellsworth\'s triangular table.  Records that beat the grid stop near the centre column, n ≈ k² − k; to its right the grid packing is the best known.',
  deficit: 'One cell per n; row k holds every n with ⌈√n⌉ = k, right-aligned, so each column is one deficit c = k² − n: the rightmost is the perfect squares, then k² − 1, k² − 2, …',
  grid: 'One cell per n, 18 to a row.',
};
function map() {
  const g = $('grid'), tt = document.querySelector('#map .tt'), box = $('map');
  g.innerHTML = ''; g.className = 'mapgrid' + (layout === 'grid' ? '' : ' tri');
  $('layoutnote').textContent = LAYOUTNOTE[layout];
  for (let n = 1; n <= NMAX; n++) {
    const h = current(n), b = h.date ? bin(h.date[0]) : -1;
    const c = document.createElement('button');
    c.className = 'cell' + (b < 0 ? ' grid' : '') + (b >= 3 ? ' dark' : '');
    if (b >= 0) c.style.background = `var(--seq${b})`;
    c.textContent = n; c.dataset.n = n;
    c.setAttribute('aria-label', `n = ${n}: ${tipText(n)}`);
    const show = e => {
      tt.textContent = tipText(n); tt.style.display = 'block';
      const bb = box.getBoundingClientRect(), cb = c.getBoundingClientRect();
      tt.style.left = Math.min(cb.left - bb.left + cb.width + 6, bb.width - 300) + 'px';
      tt.style.top = (cb.top - bb.top + cb.height + 6) + 'px';
    };
    c.onmouseenter = show; c.onfocus = show;
    c.onmouseleave = c.onblur = () => { tt.style.display = 'none'; };
    c.onclick = () => { pick(n, true); };
    const p = place(n); if (p) { c.style.gridRow = p.row; c.style.gridColumn = p.col; }
    g.appendChild(c);
  }
  paint();
}
function tipText(n) {
  const h = current(n), lines = [`n = ${n}` + (settled(n) ? '  (settled)' : '')];
  if (h.s) lines.push(`side ${h.s.toFixed(6)}`);
  if (h.kind === 'grid') lines.push('grid packing' + (isSquare(n) ? ' (perfect square)' : ', never beaten'));
  else if (h.kind === 'construction') lines.push(`built on ${(h.origin || []).map(name).join(', ')}'s packing`);
  else lines.push(`${h.kind === 'improved' ? 'improved' : 'found'} by ${(h.who || []).map(name).join(', ')}`);
  if (h.date) lines.push(`since ${dfmt(h.date)} (${Math.max(0, dnum(TODAY) - dnum(h.date)).toFixed(1)} years)`);
  const prev = upper(n).length - 1;
  if (prev > 0) lines.push(`${prev} earlier record${prev > 1 ? 's' : ''}`);
  return lines.join('\n');
}
function paint() {
  const lab = $('labels').checked;
  const mine = new Set();
  if (hi) for (let n = 1; n <= NMAX; n++) if (credited(current(n)).includes(hi)) mine.add(n);
  for (const c of $('grid').children) {
    const n = +c.dataset.n;
    c.classList.toggle('lab', lab || n === selN || mine.has(n));
    c.classList.toggle('mine', mine.has(n));
    c.classList.toggle('dim', !!hi && !mine.has(n));
    c.classList.toggle('sel', n === selN);
  }
}

// ---------- per-n chain ----------
const KIND = { found: 'found', improved: 'improved', previous: 'earlier record', construction: 'construction', grid: 'grid' };
function chain() {
  const n = selN;
  $('cn').textContent = n; $('npick').value = n; $('cexplore').href = `explore.html?n=${n}`; $('cbounds').href = `bounds.html?n=${n}`; $('ccompare').href = `compare.html?n=${n}`;
  const hs = upper(n), other = C.n[n].other || [];
  const byFile = {};
  for (const c of other) (byFile[c.event ? 'ev:' + n : c.file] = byFile[c.event ? 'ev:' + n : c.file] || []).push(c);
  const li = hs.map((h, i) => {
    const last = i === hs.length - 1;
    const when = h.kind === 'grid' ? '' : `${dfmt(h.date) || '?'}${last ? ' – now' : h.until ? ' – ' + dfmt(h.until) : ''}`;
    const src = h.url ? `<a href="${esc(h.url)}">${esc(h.text)}</a>` : plain(h.text);
    const built = (h.based_on || []).filter(r => r.who && r.who.length).map(r =>
      r.n && r.n !== n ? `built on <a href="#n=${r.n}" data-n="${r.n}">s(${r.n})</a> by ${names(r.who)}` : r.pattern ? `continues a pattern of ${names(r.who)}` : `based on a packing by ${names(r.who)}`);
    let extra = [];
    if (last) {
      const groups = {};
      for (const c of other) if (c.kind !== 'nonrecord') (groups[c.kind] = groups[c.kind] || new Set()).add(...(c.who || []));
      const lab = { refined: 'refined by', improvement: 'improved by', alternative: 'alternatives by', proved: 'proved optimal by', invalid: 'withdrawn version by' };
      extra = Object.entries(groups).filter(([k, s]) => s.size && lab[k]).map(([k, s]) => `${lab[k]} ${names([...s].filter(Boolean))}`);
    }
    const side = h.s ? `<span class="sv">s = ${h.s.toFixed(9)}</span> · ` : '';
    return `<li class="${last ? 'now' : ''}"><div class="when">${when}</div><div><div><span class="who">${holderLabel(h)}</span> <span class="handle">${KIND[h.kind] || h.kind}</span></div>
      <div class="src">${side}${src}</div>${built.length || extra.length ? `<div class="credits">${[...built, ...extra].join(' · ')}</div>` : ''}</div><div class="th" data-i="${i}" title="show this packing below"></div></li>`;
  });
  $('chainU').innerHTML = li.join('');
  evolution(n);
  const lc = C.n[n].lower || [], lo = C.n[n].lower_other || [];
  const tag = s => s === 'proved' ? '<span class="tag ok">proved</span>' : '<span class="tag open">certificate, unrefereed</span>';
  let L = lc.slice().map((l, i) => {
    const last = i === lc.length - 1;
    const who = l.who.length ? names(l.who) : esc(l.authors || l.source);
    const link = l.url ? ` · <a href="${esc(l.url)}">source</a>` : '';
    return `<li class="${last ? 'now' : ''}"><div class="when">${l.year || ''}${last ? ' – now' : ''}</div><div><div><span class="who">${who}</span> ${tag(l.status)}</div><div class="src"><span class="sv">s ≥ ${l.value.toFixed(6)}</span>${link}</div></div></li>`;
  });
  if (!lc.length) L = [`<li><div class="when"></div><div class="src">${n > 100 ? 'Floors are collected for n ≤ 100 so far; here only the area bound √n is on record.' : 'Only the area bound √n.'}</div></li>`];
  if (lo.length) L.push(`<li><div class="when"></div><div class="src">also on record, not counted: ${lo.map(l => `${l.who.length ? l.who.map(name).join(', ') : esc(l.authors || l.source)} ${l.year || ''} (s ≥ ${l.value.toFixed(4)}, ${l.status === 'gap' ? 'proof has a gap' : l.status === 'claimed' ? 'claimed without proof' : 'reported'})`).join('; ')}</div></li>`);
  $('chainL').innerHTML = L.join('');
  for (const a of $('chainU').querySelectorAll('a[data-n]')) a.onclick = e => { e.preventDefault(); pick(+a.dataset.n, true); };
}
// ---------- packings: one drawing per holder ----------
const cache = {};
function packing(n, i) {
  const h = upper(n)[i], f = (HIST[n] || [])[i];
  if (!f) {   // grid holder (or no drawing): the plain grid, row by row
    const k = Math.ceil(Math.sqrt(n) - 1e-12), sq = [];
    for (let j = 0; j < n; j++) sq.push([(j % k) + 0.5, k - Math.floor(j / k) - 0.5, 0]);
    return Promise.resolve({ s: k, squares: sq, grid: true });
  }
  return cache[f] || (cache[f] = fetch(`data/p/${encodeURIComponent(f)}`).then(r => r.ok ? r.json() : null).then(r => { if (r) r.s = +r.s; return r; }).catch(() => null));
}
const SVGNS = 'http://www.w3.org/2000/svg';
function drawSVG(rec) {
  const s = rec.s, w = s / 160;
  const sq = rec.squares.map(([cx, cy, th]) => {
    const t = ((th % 90) + 90) % 90, tilt = t > 1e-6 && t < 90 - 1e-6;
    const fill = tilt ? `hsl(${(((t > 45 ? t - 90 : t) * 4) + 360) % 360} 55% 62%)` : 'var(--sq-fill)';
    return `<rect x="-.5" y="-.5" width="1" height="1" transform="translate(${cx} ${cy}) rotate(${th})" fill="${fill}" stroke="var(--sq-stroke)" stroke-width="${w * .7}"/>`;
  }).join('');
  return `<svg xmlns="${SVGNS}" viewBox="${-0.02 * s} ${-0.02 * s} ${1.04 * s} ${1.04 * s}"><g transform="matrix(1 0 0 -1 0 ${s})"><rect x="0" y="0" width="${s}" height="${s}" fill="none" stroke="var(--wall)" stroke-width="${w * 1.4}"/>${sq}</g></svg>`;
}
let stageI = 0, stageGen = 0, playT = null;
// The catalogue kept no drawing for some holders improved in place: they show the next holder's packing.
const later = (n, i) => { const H = HIST[n] || []; return !!H[i] && H[i] === H[i + 1] && !upper(n)[i].s; };
function sideOf(n, i, rec) { const h = upper(n)[i]; return h.s || (!later(n, i) && rec && rec.s) || null; }
async function evolution(n) {
  const my = ++stageGen, hs = upper(n);
  stopPlay();
  stageI = hs.length - 1;
  const recs = await Promise.all(hs.map((_, i) => packing(n, i)));
  if (my !== stageGen) return;
  const parts = [];
  hs.forEach((h, i) => {
    const s = sideOf(n, i, recs[i]), sp = i ? sideOf(n, i - 1, recs[i - 1]) : null;
    const d = s && sp ? `<span class="delta">−${(sp - s).toPrecision(3)}</span>` : '';
    const who = h.kind === 'grid' ? 'grid' : h.kind === 'construction' ? `on ${(h.origin || []).map(name).join(', ')}` : (h.who || []).map(name).join(', ');
    if (i) parts.push('<span class="arrow">→</span>');
    parts.push(`<figure data-i="${i}"><div class="pic">${recs[i] ? drawSVG(recs[i]) : ''}</div><figcaption><b>${esc(who)}</b><br>${h.kind === 'grid' ? '' : dfmt(h.date)}<br>${s ? s.toFixed(6) : 'side not recorded'} ${d}${later(n, i) ? '<br><i>drawing: the later packing</i>' : ''}</figcaption></figure>`);
  });
  $('evo').innerHTML = parts.join('');
  for (const f of $('evo').querySelectorAll('figure')) f.onclick = () => { stopPlay(); stage(n, +f.dataset.i, recs); };
  for (const t of $('chainU').querySelectorAll('.th')) {
    const i = +t.dataset.i; t.innerHTML = recs[i] ? drawSVG(recs[i]) : '';
    t.onclick = () => { stopPlay(); stage(n, i, recs); $('evo').scrollIntoView({ behavior: 'smooth', block: 'nearest' }); };
  }
  stage(n, stageI, recs);
  const e = $('evo'); e.scrollLeft = e.scrollWidth;
}
function stage(n, i, recs) {
  const hs = upper(n), h = hs[i], rec = recs[i];
  stageI = i;
  for (const f of $('evo').querySelectorAll('figure')) f.classList.toggle('on', +f.dataset.i === i);
  $('big').innerHTML = rec ? drawSVG(rec) : '<p class="note">no drawing</p>';
  const s = sideOf(n, i, rec), sp = i ? sideOf(n, i - 1, recs[i - 1]) : null;
  const tilted = rec ? rec.squares.filter(([, , t]) => { t = ((t % 90) + 90) % 90; return t > 1e-6 && t < 90 - 1e-6; }).length : 0;
  const extra = (later(n, i) ? `<div class="note">The catalogue improved this record in place and kept no drawing of this stage; shown is the next record's packing.</div>` : '') + (rec && rec.squares.length > n ? `<div class="note">The catalogue draws this record for n = ${rec.squares.length} at the same side; drop any one square for n = ${n}.</div>` : '');
  const when = h.kind === 'grid' ? '' : `${dfmt(h.date)}${i === hs.length - 1 ? ' – now' : h.until ? ' – ' + dfmt(h.until) : ''}`;
  $('bigcap').innerHTML = `<div class="eyebrow">record ${i + 1} of ${hs.length}</div>
    <div><span class="who" style="font-weight:600">${holderLabel(h)}</span> <span class="handle">${KIND[h.kind] || h.kind}</span></div>
    <div class="note">${when}</div>
    <div class="sv" style="font-family:'IBM Plex Mono',monospace">s = ${s ? s.toFixed(9) : 'not recorded'}${sp && s ? `<br><span class="delta">${(sp - s).toPrecision(3)} smaller than the record before</span>` : ''}</div>
    <div class="note">${tilted} of ${n} squares tilted${rec && rec.grid ? ' (the plain grid)' : ''}</div>${extra}
    <div class="btns"><button id="bprev">← earlier</button><button id="bplay">${playT ? 'pause' : 'play from the start'}</button><button id="bnext">later →</button></div>`;
  $('bprev').onclick = () => { stopPlay(); if (i > 0) stage(n, i - 1, recs); };
  $('bnext').onclick = () => { stopPlay(); if (i < hs.length - 1) stage(n, i + 1, recs); };
  $('bplay').onclick = () => { if (playT) { stopPlay(); stage(n, stageI, recs); return; } let j = 0; stage(n, 0, recs); playT = setInterval(() => { if (++j >= hs.length) { stopPlay(); stage(n, hs.length - 1, recs); return; } stage(n, j, recs); }, 1400); stage(n, 0, recs); };
}
function stopPlay() { if (playT) clearInterval(playT); playT = null; }

function pick(n, scroll) {
  selN = Math.max(1, Math.min(NMAX, n | 0)); chain(); paint();
  history.replaceState(null, '', `#n=${selN}${hi ? '&who=' + hi : ''}${layout !== 'deficit' ? '&layout=' + layout : ''}`);
  if (scroll) $('chainsec').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

// ---------- hall of fame ----------
let T, sortK = 'now', openId = '';
const KEYS = { now: 'nowN', ever: 'everN', years: 'years', first: 'first', built: 'builtN', lbNow: 'lbNowN', lbEver: 'lbEverN' };
function hof() {
  const k = KEYS[sortK], asc = sortK === 'first';
  const rows = T.slice().sort((a, b) => (asc ? a[k] - b[k] : b[k] - a[k]) || b.nowN - a.nowN || b.years - a.years || b.lbNowN - a.lbNowN || name(a.id).localeCompare(name(b.id)));
  const nl = ns => ns.map(n => `<a href="#n=${n}" data-n="${n}">${n}</a>`).join('');
  const z = v => v ? v : '<span style="color:var(--rule)">·</span>';
  $('hof').querySelector('tbody').innerHTML = rows.map(r => {
    const main = `<tr class="person${openId === r.id ? ' on' : ''}" data-id="${r.id}"><td>${nameH(r.id)}${P[r.id] && P[r.id].note ? ` <span class="handle">${esc(P[r.id].note)}</span>` : ''}</td><td class="num">${z(r.nowN)}</td><td class="num">${z(r.everN)}</td><td class="num">${r.years >= 0.05 ? r.years.toFixed(1) : z(0)}</td><td class="num">${isFinite(r.first) ? r.first : ''}</td><td class="num">${z(r.builtN)}</td><td class="num">${z(r.lbNowN)}</td><td class="num">${z(r.lbEverN)}</td></tr>`;
    if (openId !== r.id) return main;
    const parts = [];
    if (r.now.length) parts.push(`records now: <span class="nlist">${nl(r.now)}</span>`);
    const past = r.ever.filter(n => !r.now.includes(n));
    if (past.length) parts.push(`held before: <span class="nlist">${nl(past)}</span>`);
    if (r.built.length) parts.push(`built on their packings: <span class="nlist">${nl(r.built)}</span>`);
    if (r.lbNow.length) parts.push(`best floor now: <span class="nlist">${nl(r.lbNow)}</span>`);
    const lbPast = r.lbEver.filter(n => !r.lbNow.includes(n));
    if (lbPast.length) parts.push(`floors before: <span class="nlist">${nl(lbPast)}</span>`);
    return main + `<tr class="detail"><td colspan="8">${parts.join('<br>') || 'credited for refinements and alternatives; see each n'}</td></tr>`;
  }).join('');
  for (const th of $('hof').querySelectorAll('th.sort')) th.classList.toggle('on', th.dataset.k === sortK);
  for (const tr of $('hof').querySelectorAll('tr.person')) tr.onclick = e => {
    if (e.target.closest('a')) return;
    openId = openId === tr.dataset.id ? '' : tr.dataset.id;
    hi = openId; $('who').value = hi; hof(); paint(); pick(selN, false);
  };
  for (const a of $('hof').querySelectorAll('a[data-n]')) a.onclick = e => { e.preventDefault(); pick(+a.dataset.n, true); };
}

// ---------- oldest / newest ----------
function lists() {
  const open = [];
  for (let n = 1; n <= NMAX; n++) { const h = current(n); if (h.date && !settled(n)) open.push({ n, h }); }
  open.sort((a, b) => dnum(a.h.date) - dnum(b.h.date) || a.n - b.n);
  $('oldest').querySelector('tbody').innerHTML = open.slice(0, 60).map(({ n, h }) =>
    `<tr><td class="num"><a href="#n=${n}" data-n="${n}">${n}</a></td><td>${holderLabel(h)}</td><td class="num">${dfmt(h.date)}</td><td class="num">${(dnum(TODAY) - dnum(h.date)).toFixed(0)}</td></tr>`).join('');
  const recent = [];
  for (let n = 1; n <= NMAX; n++) {
    const hs = upper(n), h = hs[hs.length - 1];
    if (!h.day) continue;
    h.newest = h.last_day || h.day;
    const prev = hs[hs.length - 2];
    recent.push({ n, h, prev });
  }
  recent.sort((a, b) => b.h.newest.localeCompare(a.h.newest) || a.n - b.n);
  $('newest').querySelector('tbody').innerHTML = recent.map(({ n, h, prev }) => {
    const imp = (C.n[n].other || []).filter(c => c.kind === 'improvement').flatMap(c => c.who);
    const by = names(h.who) + (imp.length ? `<br><span class="handle">improved by ${imp.map(name).join(', ')}</span>` : '');
    const before = prev ? `${prev.kind === 'construction' ? 'construction, ' : prev.kind === 'grid' ? 'grid' : (prev.who || []).map(name).join(', ') + ', '}${prev.date ? prev.date[0] : ''}` : '';
    return `<tr><td class="num"><a href="#n=${n}" data-n="${n}">${n}</a></td><td>${by}</td><td class="num">${h.newest.slice(5)}</td><td>${esc(before)}</td></tr>`;
  }).join('');
  for (const a of document.querySelectorAll('#oldest a[data-n], #newest a[data-n]')) a.onclick = e => { e.preventDefault(); pick(+a.dataset.n, true); };
}

async function main() {
  [C, HIST] = await Promise.all([fetch('data/credits.json').then(r => r.json()), fetch('data/hist.json').then(r => r.ok ? r.json() : {}).catch(() => ({}))]);
  P = C.people;
  const g = C.generated.split('-').map(Number); TODAY = [g[0], g[1]];
  T = tallies();
  const sel = $('who');
  for (const r of T.slice().filter(r => r.nowN || r.builtN).sort((a, b) => name(a.id).localeCompare(name(b.id))))
    sel.insertAdjacentHTML('beforeend', `<option value="${r.id}">${esc(name(r.id))} (${r.nowN + r.builtN})</option>`);
  const q = new URLSearchParams(location.hash.slice(1));
  if (q.get('n')) selN = Math.max(1, Math.min(NMAX, +q.get('n') || 11));
  if (q.get('who') && P[q.get('who')]) { hi = openId = q.get('who'); sel.value = hi; }
  sel.onchange = () => { hi = openId = sel.value; hof(); paint(); pick(selN, false); };
  $('labels').onchange = paint;
  if (['grid', 'tri', 'deficit'].includes(q.get('layout'))) { layout = q.get('layout'); $('layout').value = layout; }
  $('layout').onchange = () => { layout = $('layout').value; map(); pick(selN, false); };
  $('npick').onchange = () => pick(+$('npick').value, false);
  for (const th of $('hof').querySelectorAll('th.sort')) th.onclick = () => { sortK = th.dataset.k; hof(); };
  scale(); map(); chain(); hof(); lists();
}
main();
})();
