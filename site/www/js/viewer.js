/* Explore a packing: renders one packing with analysis layers. Data from data/index.json + data/p/<file>.json */
(() => {
const $ = id => document.getElementById(id);
const svgNS = 'http://www.w3.org/2000/svg';
const S = { idx: null, rec: null, file: null, n: null, asked: null, tol: 0, colormode: 'angle', sel: null, hiGroup: null,
            vb: null, groups: [] };

// ---------- helpers ----------
const el = (tag, attrs = {}, parent = null) => {
  const e = document.createElementNS(svgNS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
};
const fmt = (x, d = 4) => (Math.abs(x) < 1e-12 ? '0' : (Math.abs(x) < 1e-3 ? x.toExponential(2) : x.toFixed(d)));
const fmtDeg = t => Math.abs(t) < 1e-9 ? '0°' : (Math.abs(t - Math.round(t)) < 1e-9 ? Math.round(t) + '°' : t.toFixed(4) + '°');
const tiltOf = th => { th = ((th % 90) + 90) % 90; return th > 45 ? th - 90 : th; };
const hueFill = (tilt) => `hsl(${(tilt * 4 + 360) % 360} 55% 62%)`;
const CAT = ['#5B9BD5','#ED7D31','#70AD47','#C9A227','#8E6BBF','#DB5B8A','#3FA9A9','#B0715A','#7A8B2F','#4C6FBF','#D9713C','#5FA35F',
             '#A4527F','#3F8E7A','#B85C5C','#6F7DB5','#C48E2A','#8062A8','#3C9CC7','#9E7A4E','#5C8EA6','#B0505E','#4C9E6E'];
const catFill = k => k < 0 ? 'var(--sq-fill)' : CAT[k % CAT.length];   // -1: a square whose merged group averages to ~0 tilt
const parseQuery = () => Object.fromEntries(new URLSearchParams(location.search));
// Seven catalogue files are named `square-NN_rM%27.svg`.  Dropped into a URL unescaped the
// server decodes the %27 back to an apostrophe and the fetch 404s, so encode the name.
const dataURL = file => `data/p/${encodeURIComponent(file.replace('.svg', '.json'))}`;
// Optional annotation slots: a browser holding a cached copy of an older explore.html must not take
// the rest of the sidebar down with it, so writing to one that is not there is a no-op.
const setText = (id, text) => { const e = $(id); if (e) e.textContent = text; return !!e; };
const GAP_MAX = 0.25;    // tools/analysis.py NEAR_MAX: pairs further apart than this are not exported
const GAP_SHOW = 12;     // how many of them the panel lists and the drawing marks
const IDX_FONT = 0.28;   // square-number height at fit, in units of one square's side

// Annotations -- gap markers and square numbers -- are sized for the whole-packing view.  Past that
// zoom they would go on growing with the drawing: at gap-level zoom a two-line label stands taller
// than the stage.  So hold them at the size they had when the packing was fitted.  Below fit they
// stay tied to the squares, which is what reads correctly when a square is a few pixels across.
const annScale = () => Math.min(1, S.vb.w / (S.rec.s * 1.06));   // 1.06: the padding fit() adds
// centred boxes in world units, for keeping labels off each other
const boxHit = (a, b) => Math.abs(a.x - b.x) * 2 < a.w + b.w && Math.abs(a.y - b.y) * 2 < a.h + b.h;
// does the segment touch the box?  (Liang-Barsky; a leader run across a number reads as a strikethrough)
function segBox(x1, y1, x2, y2, q) {
  const dx = x2 - x1, dy = y2 - y1;
  const P = [-dx, dx, -dy, dy];
  const Q = [x1 - (q.x - q.w / 2), (q.x + q.w / 2) - x1, y1 - (q.y - q.h / 2), (q.y + q.h / 2) - y1];
  let t0 = 0, t1 = 1;
  for (let i = 0; i < 4; i++) {
    if (P[i] === 0) { if (Q[i] < 0) return false; continue; }
    const r = Q[i] / P[i];
    if (P[i] < 0) { if (r > t1) return false; if (r > t0) t0 = r; }
    else { if (r < t0) return false; if (r < t1) t1 = r; }
  }
  return true;
}
// Elements whose stroke and marker sizes should hold their fitted size on screen.  Each value is
// [constant, multiple of lw]; the spec is kept on the node so a zoom can retune it in place instead
// of rebuilding the layer -- the free-region layer paints and serialises a canvas per square, which
// is not a thing to redo on every frame of a wheel gesture.
function annSet(e, k) {
  const lw = S.rec.s / 700;
  for (const [key, v] of Object.entries(e.__ann))
    e.setAttribute(key, Array.isArray(v[0]) ? v.map(([c, m]) => c + m * lw * k).join(' ') : v[0] + v[1] * lw * k);
}
const elS = (tag, attrs, ann, parent) => {
  const e = el(tag, { ...attrs, class: 'ann' }, parent);
  e.__ann = ann; annSet(e, annScale()); return e;
};
const rescaleAnn = () => {
  const k = annScale();
  for (const id of ['free', 'regions', 'contacts']) for (const e of $(id).querySelectorAll('.ann')) annSet(e, k);
};

// where a gap label may go, in order of preference; straight above, then the other axes, then corners
const LABEL_DIRS = [[0, 1], [0, -1], [1, 0], [-1, 0], [0.71, 0.71], [-0.71, 0.71], [0.71, -0.71], [-0.71, -0.71]];
const textBox = (x, y, lines, fs) =>
  ({ x, y, w: Math.max(...lines.map(t => t.length)) * fs * 0.62, h: fs * lines.length * 1.25 });
// Where the square numbers sit, whether or not they are switched on.  They are pinned to square
// centres and cannot yield, so the gap labels are the ones that have to go round them.
const idxBoxes = () => !$('L_idx').checked ? [] :
  S.rec.squares.map(([cx, cy], i) => textBox(cx, cy, [String(i)], IDX_FONT * annScale()));

// merge exact angle groups with a tolerance (degrees, circular mod 90)
function mergedGroups(tol) {
  const g = S.rec.analysis.angle_groups.map(x => ({ theta: x.theta, members: x.members.slice() }));
  if (!g.length) return [];
  const out = [];
  for (const x of g) {
    const last = out[out.length - 1];
    if (last && x.theta - last.max < tol + 1e-12) { last.members.push(...x.members); last.max = x.theta; last.thetas.push(x.theta); }
    else out.push({ min: x.theta, max: x.theta, thetas: [x.theta], members: x.members.slice() });
  }
  // wrap-around: last group near 90 merges with first near 0
  if (out.length > 1 && (out[0].min + 90 - out[out.length - 1].max) < tol + 1e-12) {
    const last = out.pop(); out[0].members.push(...last.members); out[0].thetas.push(...last.thetas.map(t => t - 90)); out[0].wrapped = true;
  }
  for (const o of out) { o.mean = o.thetas.reduce((a, b) => a + b, 0) / o.thetas.length; o.tilt = tiltOf(o.mean); o.exact = o.thetas.length; }
  // sort: unrotated first, then by |tilt|
  out.sort((a, b) => (Math.abs(a.tilt) < 1e-9 ? -1 : Math.abs(b.tilt) < 1e-9 ? 1 : Math.abs(a.tilt) - Math.abs(b.tilt)));
  return out;
}

function fillFor(i) {
  const a = S.rec.analysis;
  const tilt = a.tilt[i];
  if (S.colormode === 'none') return 'var(--sq-fill)';
  if (Math.abs(tilt) < 1e-9) return 'var(--sq-fill)';
  if (S.colormode === 'angle') return hueFill(tilt);
  const k = S.groups.findIndex(g => g.members.includes(i));
  const rotatedGroups = S.groups.filter(g => Math.abs(g.tilt) >= 1e-9);
  return catFill(rotatedGroups.findIndex(g => g.members.includes(i)));
}

// ---------- rendering ----------
function render() {
  const svg = $('view'); svg.innerHTML = '';
  const r = S.rec, a = r.analysis, s = r.s;
  const lw = s / 700;
  const world = el('g', { id: 'world', transform: `matrix(1 0 0 -1 0 ${s})` }, svg);
  el('rect', { x: 0, y: 0, width: s, height: s, fill: 'var(--surface)', stroke: 'var(--wall)', 'stroke-width': lw * 2.2 }, world);
  const sqG = el('g', { id: 'squares' }, world);
  r.squares.forEach((q, i) => {
    const [cx, cy, th] = q;
    const g = el('g', { transform: `translate(${cx} ${cy}) rotate(${th})`, class: 'sq', 'data-i': i }, sqG);
    el('rect', { x: -0.5, y: -0.5, width: 1, height: 1, fill: fillFor(i), stroke: 'var(--sq-stroke)', 'stroke-width': lw, 'stroke-linejoin': 'round' }, g);
  });
  el('g', { id: 'regions' }, world); el('g', { id: 'free' }, world); el('g', { id: 'contacts' }, world); el('g', { id: 'gapmarks' }, world); el('g', { id: 'sym' }, world); el('g', { id: 'selmark' }, world); el('g', { id: 'labels' }, world);
  if (!S.vb) fit();          // the annotation layers size themselves against the viewBox, so set it first
  drawFree(); drawContacts(); drawGaps(); drawSym(); drawLabels(); markSel();
  annK = annScale();
  applyVB();
}

function drawFree() {
  const a = S.rec.analysis, s = S.rec.s, lw = s / 700, on = $('L_free').checked;
  const g = $('free'), rg = $('regions'); g.innerHTML = ''; rg.innerHTML = '';
  if (!on) return;
  const free = new Set(a.free), wedged = new Set(a.wedged);
  S.rec.squares.forEach((q, i) => {
    if (!free.has(i) && !wedged.has(i)) return;   // squares that only move together with neighbours are reported in the tooltip
    const [cx, cy, th] = q;
    const cls = free.has(i) ? 'free' : 'wedged';
    const stroke = cls === 'wedged' ? 'var(--gold)' : 'var(--green)';
    const gg = el('g', { transform: `translate(${cx} ${cy}) rotate(${th})` }, g);
    elS('rect', { fill: 'none', stroke, ...(cls === 'wedged' ? {} : { 'stroke-dasharray': 'none' }) },
      { x: [-0.5, 2], y: [-0.5, 2], width: [1, -4], height: [1, -4], 'stroke-width': [0, 2.5],
        ...(cls === 'wedged' ? { 'stroke-dasharray': [[0, 2], [0, 4]] } : {}) }, gg);
    const reg = a.regions[i];
    if (reg && reg.rle && !reg.mask) {   // decode run-length rows
      reg.mask = reg.rle.split(';').map(row => { let out = '', ch = '0'; for (const n of row.split(',')) { out += ch.repeat(+n); ch = ch === '0' ? '1' : '0'; } return out; });
    }
    if (reg && reg.slide && reg.kind === 'slide') {
      // a channel of zero width: there is no area to shade, so show the far end and the direction
      const [ux, uy] = reg.slide.dir, t = reg.slide.len;
      const gh = el('g', { transform: `translate(${cx + ux * t} ${cy + uy * t}) rotate(${th})` }, rg);
      elS('rect', { x: -0.5, y: -0.5, width: 1, height: 1, fill: 'none', stroke: 'var(--green)', opacity: .8, 'stroke-linejoin': 'round' },
        { 'stroke-width': [0, 1.5], 'stroke-dasharray': [[0, 4], [0, 4]] }, gh);
      elS('line', { x1: cx, y1: cy, x2: cx + ux * t, y2: cy + uy * t, stroke: 'var(--green)', 'stroke-linecap': 'round' },
        { 'stroke-width': [0, 2.5] }, rg);
      elS('circle', { cx, cy, fill: 'var(--green)' }, { r: [0, 2.5] }, rg);
    }
    if (reg && reg.mask && reg.free) {
      // translation region of the centre, drawn from the sampled mask
      const N = reg.N, R = reg.R;
      const c = document.createElement('canvas'); c.width = N; c.height = N;
      const ctx = c.getContext('2d'); ctx.fillStyle = 'rgba(42,110,79,0.55)';
      for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) if (reg.mask[y][x] === '1') ctx.fillRect(x, y, 1, 1);
      // mask row 0 is dy = -R; the world group is y-flipped, so the image needs no extra flip
      el('image', { href: c.toDataURL(), x: cx - R, y: cy - R, width: 2 * R, height: 2 * R, preserveAspectRatio: 'none', style: 'image-rendering:pixelated' }, rg);
      elS('circle', { cx, cy, fill: 'var(--green)' }, { r: [0, 2.5] }, rg);
    }
  });
}

function drawContacts() {
  const a = S.rec.analysis, s = S.rec.s, lw = s / 700, g = $('contacts'); g.innerHTML = '';
  if (!$('L_contacts').checked) return;
  for (const c of a.contacts) {
    const col = c.type === 'corner-corner' ? 'var(--red)' : c.type.startsWith('edge') ? 'var(--accent)' : 'var(--gold)';
    if (c.type === 'edge-edge' || c.type === 'edge-wall') {
      if (c.points.length >= 2) elS('line', { x1: c.points[0][0], y1: c.points[0][1], x2: c.points[1][0], y2: c.points[1][1], stroke: col, 'stroke-linecap': 'round', opacity: .8 },
        { 'stroke-width': [0, 3] }, g);
    } else {
      elS('circle', { cx: c.points[0][0], cy: c.points[0][1], fill: col, stroke: 'var(--surface)' }, { r: [0, 3.5], 'stroke-width': [0, 1] }, g);
    }
  }
}

// nearest-point marker for a near miss, computed here from the squares
function sqCorners([cx, cy, th]) {
  const c = Math.cos(th * Math.PI / 180), s = Math.sin(th * Math.PI / 180);
  return [[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]].map(([dx, dy]) => [cx + c * dx - s * dy, cy + s * dx + c * dy]);
}
// separation of two convex quads along the best separating axis (negative = overlap)
function polyGap(A, B) {
  let best = -Infinity;
  for (const [P, Q] of [[A, B], [B, A]]) for (let k = 0; k < 4; k++) {
    const [x1, y1] = P[k], [x2, y2] = P[(k + 1) % 4];
    let nx = y2 - y1, ny = x1 - x2; const L = Math.hypot(nx, ny); nx /= L; ny /= L;
    const pP = P.map(p => p[0] * nx + p[1] * ny), pQ = Q.map(p => p[0] * nx + p[1] * ny);
    best = Math.max(best, Math.min(...pQ) - Math.max(...pP));
  }
  return best;
}
const WALLNAME = ['left wall', 'right wall', 'bottom wall', 'top wall'];
// what stops a square that can only slide: the gap that *closes*, not the neighbours it slides
// along (those stay in contact the whole way, so they would otherwise always win on final gap).
function slideTarget(i) {
  const r = S.rec.analysis.regions[i]; if (!r || !r.slide) return null;
  const [ux, uy] = r.slide.dir, t = r.slide.len, s = S.rec.s, q = S.rec.squares[i];
  const A0 = sqCorners(q), A1 = sqCorners([q[0] + ux * t, q[1] + uy * t, q[2]]);
  let best = null;
  const take = (g0, g1, label) => { if (g0 > 1e-9 && g1 < 1e-6 && (!best || g0 < best.g0)) best = { g0, label }; };
  S.rec.squares.forEach((p, j) => {
    if (j === i || Math.hypot(p[0] - q[0], p[1] - q[1]) > 1.5 + t) return;
    const B = sqCorners(p); take(polyGap(A0, B), polyGap(A1, B), `square ${j}`);
  });
  const walls = C => [Math.min(...C.map(p => p[0])), s - Math.max(...C.map(p => p[0])),
                      Math.min(...C.map(p => p[1])), s - Math.max(...C.map(p => p[1]))];
  const w0 = walls(A0), w1 = walls(A1);
  for (let w = 0; w < 4; w++) take(w0[w], w1[w], WALLNAME[w]);
  return best ? best.label : null;
}
// one sentence describing how square i is free to move
function freeText(i) {
  const r = S.rec.analysis.regions[i]; if (!r || !r.dx) return 'free';
  const rr = Math.max(...(r.rot || [0, 0]).map(Math.abs)), turn = rr > 1e-3 ? `, turns ${fmt(rr, 1)}°` : '';
  if (r.kind === 'slide' && r.slide) {
    const tgt = slideTarget(i);
    return `slides ${fmt(r.slide.len, 3)} — that one direction only${tgt ? `, up against ${tgt}` : ''}${turn}`;
  }
  if (r.kind === 'turn') return `turns ${fmt(rr, 2)}° on the spot`;
  return `slides ${fmt(-r.dx[0], 2)}←${fmt(r.dx[1], 2)}→ ${fmt(-r.dy[0], 2)}↓${fmt(r.dy[1], 2)}↑${turn}`;
}

// Where to draw a near miss: halfway across the gap, and centred on the part of the two squares
// that actually faces across it.  Taking the single closest vertex instead puts the marker at the
// end of a facing edge, which reads as measuring the wrong pair.
function nearPoint(m) {
  const s = S.rec.s;
  const A = sqCorners(S.rec.squares[m.i]);
  if (m.j == null) {   // wall
    const w = m.wall, ax = w < 2 ? 0 : 1, val = (w === 1 || w === 3) ? s : 0;
    const d = A.map(p => Math.abs(p[ax] - val)), dmin = Math.min(...d);
    const touching = A.filter((p, k) => d[k] < dmin + 1e-9);          // the corner, or the whole edge
    const along = touching.map(p => p[1 - ax]);
    const pt = []; pt[ax] = val + (val === 0 ? dmin / 2 : -dmin / 2);
    pt[1 - ax] = (Math.min(...along) + Math.max(...along)) / 2;
    return pt;
  }
  const B = sqCorners(S.rec.squares[m.j]);
  let best = null;
  for (const [P, Q] of [[A, B], [B, A]]) for (let k = 0; k < 4; k++) {
    const [x1, y1] = P[k], [x2, y2] = P[(k + 1) % 4];
    let nx = y2 - y1, ny = x1 - x2; const L = Math.hypot(nx, ny); nx /= L; ny /= L;
    const pP = P.map(p => p[0] * nx + p[1] * ny), pQ = Q.map(p => p[0] * nx + p[1] * ny);
    const gap = Math.min(...pQ) - Math.max(...pP);
    if (!best || gap > best.gap) best = { gap, nx, ny, P, Q, pmax: Math.max(...pP) };
  }
  const { gap, nx, ny, P, Q, pmax } = best;
  const tx = -ny, ty = nx;                                            // along the facing edge
  const tP = P.map(p => p[0] * tx + p[1] * ty), tQ = Q.map(p => p[0] * tx + p[1] * ty);
  const lo = Math.max(Math.min(...tP), Math.min(...tQ)), hi = Math.min(Math.max(...tP), Math.max(...tQ));
  // where the two footprints overlap along the edge; if they do not, meet in between
  const t = lo <= hi ? (lo + hi) / 2 : (Math.max(Math.min(...tP), Math.min(...tQ)) + Math.min(Math.max(...tP), Math.max(...tQ))) / 2;
  const n = pmax + gap / 2;
  return [nx * n + tx * t, ny * n + ty * t];
}
const gapLabel = m => `${m.i}–${m.j == null ? WALLNAME[m.wall] : m.j}`;
function drawGaps() {
  const a = S.rec.analysis, s = S.rec.s, lw = s / 700, g = $('gapmarks'); g.innerHTML = '';
  if (!$('L_gaps').checked) return;
  const k = annScale(), fs = lw * 11 * k;
  const placed = idxBoxes();   // several gaps can share a hole, and the square numbers are already there
  (a.near_misses || []).slice(0, GAP_SHOW).forEach(m => {
    const p = nearPoint(m);
    // the circle shows the gap at its true size; the floor only keeps it findable when zoomed out
    const r = Math.max(m.gap * 2, lw * 9 * k);
    el('circle', { cx: p[0], cy: p[1], r, fill: 'none', stroke: 'var(--red)', 'stroke-width': lw * 1.5 * k,
      'stroke-dasharray': `${lw * 3 * k} ${lw * 2 * k}` }, g);
    const lines = [gapLabel(m), fmt(m.gap, 4)];
    const base = r + lw * 16 * k, step = lw * 30 * k, lift = lw * 6 * k;
    // Somewhere to put the label.  Straight above is the default, but a gap inside a column of
    // squares has their numbers above it and below it, and nothing on the vertical is ever free --
    // so try the other directions too, then further out.  A candidate is refused if the label box
    // or its leader would touch anything already placed.  Only the neighbourhood can collide, and
    // a catalogue entry can hold 9465 numbers, so filter to it first.
    const R = base + 3 * step + 1;
    const near = placed.filter(q => Math.abs(q.x - p[0]) < R && Math.abs(q.y - p[1]) < R);
    let spot = null;
    for (let ring = 0; ring < 4 && !spot; ring++) for (const [ux, uy] of LABEL_DIRS) {
      const d = base + ring * step, b = textBox(p[0] + ux * d, p[1] + uy * d, lines, fs);
      if (near.some(q => boxHit(q, b))) continue;
      // leader from the rim of the circle to the edge of the label box
      const exit = Math.min(Math.abs(ux) > 1e-9 ? b.w / 2 / Math.abs(ux) : Infinity,
                            Math.abs(uy) > 1e-9 ? b.h / 2 / Math.abs(uy) : Infinity);
      const seg = [p[0] + ux * r, p[1] + uy * r, b.x - ux * exit, b.y - uy * exit];
      if (near.some(q => segBox(...seg, q))) continue;
      spot = { b, seg, lead: !(ring === 0 && ux === 0 && uy === 1) };
      break;
    }
    // nothing anywhere is clear: keep it next to its own circle rather than flung across the drawing
    if (!spot) spot = { b: textBox(p[0], p[1] + base, lines, fs), lead: false };
    const b = spot.b; placed.push(b);
    if (spot.lead)
      el('line', { x1: spot.seg[0], y1: spot.seg[1], x2: spot.seg[2], y2: spot.seg[3],
        stroke: 'var(--red)', 'stroke-width': lw * k, opacity: .5 }, g);
    // which pair, then how far: a bare number leaves you guessing what it measures
    // the label sits on top of the drawing, so give it a halo in the background colour: without one
    // a square's edge runs through the digits and reads as a strikethrough
    const lab = el('text', { 'font-size': fs, 'font-family': 'IBM Plex Mono, monospace', fill: 'var(--red)', 'text-anchor': 'middle',
      stroke: 'var(--surface)', 'stroke-width': lw * 2.6 * k, 'stroke-linejoin': 'round', 'paint-order': 'stroke fill',
      transform: `translate(${b.x} ${b.y - lift}) scale(1 -1)` }, g);
    el('tspan', { x: 0, y: -lw * 12 * k }, lab).textContent = lines[0];
    el('tspan', { x: 0, y: 0 }, lab).textContent = lines[1];
  });
}

function drawSym() {
  const a = S.rec.analysis, s = S.rec.s, lw = s / 700, g = $('sym'); g.innerHTML = '';
  if (!$('L_sym').checked) return;
  const k = annScale();
  const e = a.symmetry.elements, st = { stroke: 'var(--accent)', 'stroke-width': lw * 2 * k, 'stroke-dasharray': `${lw * 8 * k} ${lw * 6 * k}`, opacity: .9 };
  if (e.includes('mirror_v')) el('line', { x1: s / 2, y1: 0, x2: s / 2, y2: s, ...st }, g);
  if (e.includes('mirror_h')) el('line', { x1: 0, y1: s / 2, x2: s, y2: s / 2, ...st }, g);
  if (e.includes('mirror_d')) el('line', { x1: 0, y1: 0, x2: s, y2: s, ...st }, g);
  if (e.includes('mirror_a')) el('line', { x1: 0, y1: s, x2: s, y2: 0, ...st }, g);
  if (e.includes('rot180') || e.includes('rot90')) {
    el('circle', { cx: s / 2, cy: s / 2, r: lw * 10 * k, fill: 'none', stroke: 'var(--accent)', 'stroke-width': lw * 2 * k }, g);
    const t = el('text', { 'font-size': lw * 12 * k, fill: 'var(--accent)', 'text-anchor': 'middle', 'font-family': 'IBM Plex Mono, monospace',
      transform: `translate(${s / 2} ${s / 2 - lw * 24 * k}) scale(1 -1)` }, g);
    t.textContent = e.includes('rot90') ? '90°' : '180°';
  }
}

// The numbers are built at IDX_FONT and scaled by the group transform, so a zoom only has to
// retune `transform` (see rescaleIdx) -- a catalogue entry can hold 9465 of them, and rebuilding
// that many text nodes costs ~200 ms, which is not a thing to do on every frame of a wheel gesture.
function drawLabels() {
  const g = $('labels'); g.innerHTML = '';
  if (!$('L_idx').checked) return;
  const k = annScale();
  S.rec.squares.forEach(([cx, cy], i) => {
    const t = el('text', { 'font-size': IDX_FONT, 'text-anchor': 'middle', 'dominant-baseline': 'middle', fill: 'var(--ink)',
      'font-family': 'IBM Plex Mono, monospace', opacity: .8, transform: `translate(${cx} ${cy}) scale(${k} ${-k})` }, g);
    t.dataset.at = `translate(${cx} ${cy})`;
    t.textContent = i;
  });
}
function rescaleIdx() {
  const k = annScale();
  for (const t of $('labels').children) t.setAttribute('transform', `${t.dataset.at} scale(${k} ${-k})`);
}

function recolor() {
  S.groups = mergedGroups(S.tol);
  document.querySelectorAll('#squares .sq').forEach(g => {
    const i = +g.dataset.i; const rect = g.firstChild; rect.setAttribute('fill', fillFor(i));
    rect.setAttribute('opacity', S.hiGroup == null || S.hiGroup.includes(i) ? 1 : 0.25);
  });
  renderLegend();
}

// ---------- viewBox / zoom ----------
function fit() { const s = S.rec.s, pad = s * 0.03; S.vb = { x: -pad, y: -pad, w: s + 2 * pad, h: s + 2 * pad }; }
let annK = null, annQueued = false;
function applyVB() {
  const v = S.vb; $('view').setAttribute('viewBox', `${v.x} ${v.y} ${v.w} ${v.h}`);
  // annScale() is clamped at 1, so zooming out past fit changes nothing; panning never does.
  if (S.rec && annScale() !== annK) scheduleAnn();
}
// A wheel gesture fires many events per painted frame, so coalesce them into one.
function scheduleAnn() {
  if (annQueued) return;
  annQueued = true;
  requestAnimationFrame(() => {
    annQueued = false; if (!S.rec) return;
    annK = annScale();
    rescaleAnn();                          // free / regions / contacts, retuned in place
    drawGaps(); drawSym(); markSel();      // cheap layers, and their text offsets move with the scale
    rescaleIdx();
  });
}
function svgPoint(evt) {
  const svg = $('view'), pt = svg.createSVGMatrix ? svg.createSVGPoint() : null;
  pt.x = evt.clientX; pt.y = evt.clientY;
  const p = pt.matrixTransform(svg.getScreenCTM().inverse());
  return { x: p.x, y: S.rec.s - p.y };   // world coords (y-up)
}
function zoomTo(x, y, w) { S.vb = { x: x - w / 2, y: (S.rec.s - y) - w / 2, w, h: w }; applyVB(); }
function initZoom() {
  const svg = $('view'); let drag = null;
  svg.addEventListener('wheel', e => {
    e.preventDefault(); const v = S.vb, f = Math.exp(e.deltaY * 0.0015);
    const svgP = svg.createSVGPoint(); svgP.x = e.clientX; svgP.y = e.clientY; const p = svgP.matrixTransform(svg.getScreenCTM().inverse());
    S.vb = { x: p.x - (p.x - v.x) * f, y: p.y - (p.y - v.y) * f, w: v.w * f, h: v.h * f }; applyVB();
  }, { passive: false });
  // the square under the pointer is recorded on pointerdown: setPointerCapture retargets the
  // pointerup to the <svg>, so closest('.sq') on the up event would always miss.
  svg.addEventListener('pointerdown', e => { drag = { x: e.clientX, y: e.clientY, vb: { ...S.vb }, i: squareAt(e) }; svg.setPointerCapture(e.pointerId); svg.classList.add('dragging'); });
  svg.addEventListener('pointermove', e => {
    if (drag) { const sc = S.vb.w / svg.clientWidth; S.vb = { ...drag.vb, x: drag.vb.x - (e.clientX - drag.x) * sc, y: drag.vb.y - (e.clientY - drag.y) * sc }; applyVB(); }
    hover(e);
  });
  svg.addEventListener('pointerup', e => {
    const d = drag, moved = d && Math.hypot(e.clientX - d.x, e.clientY - d.y) > 3;
    drag = null; svg.classList.remove('dragging');
    if (d && !moved) select(d.i);
  });
  svg.addEventListener('pointerleave', () => { $('tip').style.display = 'none'; });
  $('reset').onclick = () => { fit(); applyVB(); };
}

// ---------- interaction ----------
function squareAt(e) { const t = e.target.closest && e.target.closest('.sq'); return t ? +t.dataset.i : null; }
function hover(e) {
  const i = squareAt(e), tip = $('tip');
  if (i == null) { tip.style.display = 'none'; return; }
  const a = S.rec.analysis, q = S.rec.squares[i];
  const gk = S.groups.findIndex(g => g.members.includes(i));
  const lines = [`square ${i}`, `tilt ${fmtDeg(a.tilt[i])}${S.groups.length ? `  (rotation group ${gk + 1} of ${S.groups.length})` : ''}`,
    `centre (${q[0].toFixed(4)}, ${q[1].toFixed(4)})`];
  if (a.free.includes(i)) lines.push('free: ' + freeText(i));
  else if (a.wedged.includes(i)) lines.push('wedged: has first-order slack, but no actual motion was found');
  else if (a.mobile[i]) lines.push('can move, but not on its own: only together with neighbours or by sliding along them');
  else lines.push('jammed');
  tip.textContent = lines.join('\n'); tip.style.display = 'block';
  const st = $('view').parentElement.getBoundingClientRect();
  tip.style.left = (e.clientX - st.left + 14) + 'px'; tip.style.top = (e.clientY - st.top + 14) + 'px';
}
function select(i) { S.sel = i; renderSel(); markSel(); }
// outline the selected square in the drawing
function markSel() {
  const s = S.rec.s, lw = s / 700, g = $('selmark'); if (!g) return;
  g.innerHTML = '';
  if (S.sel == null) return;
  const k = annScale();
  const [cx, cy, th] = S.rec.squares[S.sel];
  const gg = el('g', { transform: `translate(${cx} ${cy}) rotate(${th})` }, g);
  el('rect', { x: -0.5, y: -0.5, width: 1, height: 1, fill: 'none', stroke: 'var(--ink)', 'stroke-width': lw * 3.5 * k, 'stroke-linejoin': 'round' }, gg);
  el('rect', { x: -0.5, y: -0.5, width: 1, height: 1, fill: 'none', stroke: 'var(--surface)', 'stroke-width': lw * 1.5 * k, 'stroke-linejoin': 'round', 'stroke-dasharray': `${lw * 6 * k} ${lw * 5 * k}` }, gg);
}

function renderSel() {
  const box = $('sel'), i = S.sel;
  if (i == null) { box.textContent = 'Click a square.'; return; }
  const a = S.rec.analysis, q = S.rec.squares[i];
  const cs = a.contacts.filter(c => c.i === i || c.j === i);
  const wallName = WALLNAME;
  let h = `<div><b class="mono">square ${i}</b> · tilt ${fmtDeg(a.tilt[i])} · centre (${q[0].toFixed(6)}, ${q[1].toFixed(6)})</div>`;
  h += `<div class="note">${cs.length} contact${cs.length === 1 ? '' : 's'}:</div><div class="legend">`;
  for (const c of cs) {
    const other = c.j == null ? wallName[c.wall] : `square ${c.i === i ? c.j : c.i}`;
    h += `<div class="row">${c.type.replace('-', ' on ')} · ${other}${c.exact === false ? ` <span class="tag warn">gap ${fmt(c.gap_mp ?? 0, 2)}</span>` : ''}</div>`;
  }
  h += '</div>';
  box.innerHTML = h;
}

function renderLegend() {
  const L = $('legend'); L.innerHTML = '';
  const a = S.rec.analysis, n = S.rec.n;
  const rot = S.groups.filter(g => Math.abs(g.tilt) >= 1e-9);
  $('rotsum').textContent = rot.length === 0 ? 'No tilted squares.' :
    `${S.groups.length} distinct rotation${S.groups.length === 1 ? '' : 's'}` + (S.tol > 0 ? ` (merging angles within ${fmtDeg(S.tol)})` : '') + `; ${a.n_rotated} of ${n} squares tilted.`;
  S.groups.forEach((g, k) => {
    const row = document.createElement('div'); row.className = 'row';
    const sw = document.createElement('span'); sw.className = 'sw';
    sw.style.background = Math.abs(g.tilt) < 1e-9 ? 'var(--sq-fill)' : (S.colormode === 'group' ? catFill(rot.indexOf(g)) : S.colormode === 'angle' ? hueFill(g.tilt) : 'var(--sq-fill)');
    row.appendChild(sw);
    const lab = document.createElement('span');
    lab.textContent = fmtDeg(g.tilt) + (g.exact > 1 ? ` (${g.exact} nearly-equal angles, spread ${fmt(g.max - g.min, 4)}°)` : '');
    row.appendChild(lab);
    const cnt = document.createElement('span'); cnt.className = 'cnt'; cnt.textContent = g.members.length; row.appendChild(cnt);
    row.onmouseenter = () => { S.hiGroup = g.members; recolorOnly(); }; row.onmouseleave = () => { S.hiGroup = null; recolorOnly(); };
    L.appendChild(row);
  });
}
function recolorOnly() {
  document.querySelectorAll('#squares .sq').forEach(g => { const i = +g.dataset.i; g.firstChild.setAttribute('opacity', S.hiGroup == null || S.hiGroup.includes(i) ? 1 : 0.25); });
}

const TABLE_MAX = 324;   // Ellsworth's main table: every n up to here with a packing better than the grid
function renderSide() {
  const r = S.rec, a = r.analysis, meta = S.idx.files[S.file] || {}, recMeta = S.idx.records[S.n] || {};
  const isRecord = recMeta.svg === S.file, inTable = !!S.idx.records[S.n];
  document.querySelector('.side .eyebrow').textContent = isRecord ? 'best known packing' : inTable ? 'alternative / older packing' : 'catalogued packing (not in the main table)';
  $('title').textContent = `${r.n} squares`;
  // Two different mismatches, and conflating them misreports both: the catalogue may cover two n
  // with one entry (s equal for both, so the picture holds the other count), or the n asked for may
  // have no entry at all, in which case gotoN snapped to the nearest that has one.
  // Absence means different things either side of TABLE_MAX: up to it the main table lists every n
  // that beats the grid, but past it the catalogue covers only selected n.
  setText('shownfor',
    S.asked != null && S.asked !== S.n
      ? (S.asked <= TABLE_MAX
        ? `Ellsworth's table has no entry for n = ${S.asked}: no packing better than the trivial one — the ${S.asked} squares set square in a container of side ⌈√${S.asked}⌉ = ${Math.ceil(Math.sqrt(S.asked))} — is recorded. Showing the nearest n that does have an entry, n = ${S.n}.`
        : `The catalogue has no drawn packing for n = ${S.asked}; above n = ${TABLE_MAX} it covers only selected n. Showing the nearest n that has one, n = ${S.n}.`)
      : r.n === S.n ? ''
      : `You asked for n = ${S.n}. The catalogue covers n = ${Math.min(S.n, r.n)} and n = ${Math.max(S.n, r.n)} with one entry, because s is the same for both; this is the packing it pictures.`);
  const sShown = (+r.s).toFixed(12).replace(/0+$/, '').replace(/\.$/, '');
  $('sval').textContent = 's = ' + sShown + (+sShown === +r.s ? '' : '…');   // no ellipsis on an exact 2, 3, …
  // exact form.  The record's closed form and prose describe the record: never lend them to an alternative.
  const tex = meta.s_tex || (isRecord ? recMeta.s_tex : '') || '';
  const polys = meta.polys && meta.polys.length ? meta.polys : isRecord ? recMeta.polys : null;
  let closed = tex.replace(/\\Nn\{[^}]*\}/, '').replace(/\\(begin|end)\{aligned\}/g, '').replace(/\\\\\s*&?=\s*$/, '').replace(/&/g, '').replace(/=\s*$/, '').trim();
  const lock = closed.match(/\{\}\^\{(\d+)\}🔒/);
  const sf = $('sform'); sf.innerHTML = '';
  if (lock) {
    sf.innerHTML = `<span>root of a degree-${lock[1]} polynomial</span>`;
    if (polys && polys.length) {
      const d = document.createElement('details'); d.innerHTML = `<summary class="note">show polynomial</summary>`;
      const p = document.createElement('div'); p.style.overflowX = 'auto'; d.appendChild(p);
      try { katex.render(polys[0].replace(/=0$/, ' = 0'), p, { throwOnError: false, displayMode: false }); } catch (e) { p.textContent = polys[0]; }
      sf.appendChild(d);
    }
  } else if (closed && closed !== 's') {
    try { katex.render(closed, sf, { throwOnError: false }); } catch (e) { sf.textContent = closed; }
  }
  // tags
  const tags = [];
  tags.push(`<span class="tag neutral">${a.category === 'trivial' ? 'no tilts' : a.category === 'diagonal' ? '45° only' : 'tilted'}</span>`);
  tags.push(a.rigid ? `<span class="tag ok">rigid</span>` : `<span class="tag open">${a.free.length} free · ${a.mobile.filter(Boolean).length} movable</span>`);
  tags.push(`<span class="tag neutral">${a.symmetry.class}</span>`);
  const inexact = a.contacts.filter(c => c.exact === false).length;
  if (inexact) tags.push(`<span class="tag warn" title="some contacts are only numerically closed (packing not analytically optimised)">${inexact} numeric contacts</span>`);
  $('tags').innerHTML = tags.join(' ');
  const prose = (meta.prose || (isRecord ? recMeta.prose : '') || '').replace(/\n/g, ' ');
  $('prose').innerHTML = prose ? prose.replace(/\$([^$]+)\$/g, (m, t) => { try { return katex.renderToString(t, { throwOnError: false }); } catch (e) { return m; } }) : '';
  // Not encoded: seven names hold a literal %27, which Ellsworth's own links use as-is.
  $('srcline').innerHTML = `Geometry from <a href="https://kingbird.myphotos.cc/packing/${S.file}" target="_blank" rel="noopener">${S.file}</a> on David Ellsworth's <a href="https://kingbird.myphotos.cc/packing/squares_in_squares.html" target="_blank" rel="noopener">Squares in Squares</a> (after Erich Friedman).`;
  // where to go next for this n
  const nv = variantsFor(S.n).length, xl = [];
  if (nv > 1) xl.push(`<a href="compare.html?n=${S.n}">compare all ${nv} packings of ${S.n}</a>`);
  if (S.n <= 100) xl.push(`<a href="bounds.html?n=${S.n}">floor and history for n = ${S.n}</a>`);
  if (S.n === 12) xl.push(`<a href="s12/">our s(12) ≥ 3.9686</a>`);
  if (S.n === 13) xl.push(`<a href="proofs.html#bentz-proof">how s(13) = 4 is proved</a>`);
  $('xlinks').innerHTML = xl.join(' · ');
  // freedom
  const fl = $('freelist'); fl.innerHTML = '';
  const groupOnly = a.mobile.map((m, i) => m && !a.free.includes(i) && !a.wedged.includes(i) ? i : null).filter(x => x != null);
  // Ellsworth's rigid page says whether he considers this one rigid; our analysis is independent.
  const claim = meta.rigid_claim;
  setText('rigidcheck',
    claim === 'yes' ? (a.rigid ? 'Ellsworth calls this packing rigid too.' : 'Ellsworth calls this packing rigid; our analysis does not.') :
    claim === 'no' ? (a.rigid ? 'Ellsworth calls this one not fully rigid; our analysis finds it rigid.' : 'Ellsworth calls this one not fully rigid either.') : '');
  const slides = a.free.filter(i => (a.regions[i] || {}).kind === 'slide').length;
  $('freesum').textContent = a.rigid ? 'Rigid: no square can move at all (to first order).' :
    `${a.free.length} square${a.free.length === 1 ? '' : 's'} can move on ${a.free.length === 1 ? 'its' : 'their'} own` +
    (slides ? ` (${slides} of them only along a single line, wedged between parallel faces)` : '') +
    (groupOnly.length ? `, ${groupOnly.length} more only together with neighbours` : '') + (a.wedged.length ? `, ${a.wedged.length} wedged (first-order slack only; no actual motion found)` : '') + '.';
  const addRow = (i, txt, col) => { const row = document.createElement('div'); row.className = 'row'; row.innerHTML = `<span class="sw" style="background:${col}"></span><span>square ${i}: ${txt}</span>`; row.onclick = () => { select(i); const q = r.squares[i]; zoomTo(q[0], q[1], 3); }; fl.appendChild(row); };
  for (const i of a.free) addRow(i, freeText(i), 'var(--green)');
  for (const i of groupOnly) addRow(i, 'moves with or along neighbours', 'var(--green-soft)');
  for (const i of a.wedged) addRow(i, 'wedged', 'var(--gold)');
  // gaps
  const G = $('gaps'); G.innerHTML = '';
  const near = a.near_misses || [];
  setText('gapnote', `Closest non-touching pairs (square–square or square–wall), in side lengths. ` +
    (near.length > GAP_SHOW ? `The ${GAP_SHOW} smallest shown; only` : 'Only') +
    ` gaps below ${GAP_MAX} are measured. Click a row to zoom to it.`);
  near.slice(0, GAP_SHOW).forEach(m => {
    const row = document.createElement('div'); row.className = 'row';
    row.innerHTML = `<span>${gapLabel(m)}</span><span class="g">${fmt(m.gap_mp ?? m.gap, 5)}</span>`;
    row.onclick = () => { $('L_gaps').checked = true; drawGaps(); const p = nearPoint(m); zoomTo(p[0], p[1], Math.max(m.gap * 12, 0.6)); };
    G.appendChild(row);
  });
  if (!near.length) G.innerHTML = `<div class="note">Nothing within ${GAP_MAX} that is not already touching.</div>`;
  // stats
  const st = $('stats'); const cnt = t => a.contacts.filter(c => c.type === t).length;
  const exact = a.contacts.filter(c => c.exact !== false).length;
  st.innerHTML = [['wasted area', `${(r.s * r.s - r.n).toFixed(4)} (${(100 * (1 - r.n / (r.s * r.s))).toFixed(2)}%)`],
    ['contacts', `${a.contacts.length}: ${cnt('edge-edge')} edge–edge, ${cnt('corner-edge')} corner–edge, ${cnt('corner-corner')} corner–corner, ${cnt('edge-wall') + cnt('corner-wall')} on walls`],
    ['exact contacts', `${exact} of ${a.contacts.length}`],
    ['symmetry', a.symmetry.class], ['distinct angles', a.n_angles]].map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
}

// ---------- navigation ----------
// Other packings of the same n worth offering. Optimiser start configurations (f.start) are real
// packings but never competitive, so they are not listed -- a ?p= link to one still loads.
function variantsFor(n, keep) {
  const files = Object.values(S.idx.files).filter(f => f.n === n && f.n_parsed && !f.start && !(f.errors && f.errors.length)).map(f => f.svg);
  if (keep && !files.includes(keep)) files.push(keep);
  const rec = S.idx.records[n] && S.idx.records[n].svg;
  files.sort((a, b) => (a === rec ? -1 : b === rec ? 1 : a.localeCompare(b)));
  return files;
}
function recordNs() { return Object.keys(S.idx.records).map(Number).sort((a, b) => a - b); }
async function load(file, n) {
  // The tooltip is otherwise cleared only on pointerleave, so changing packings with the pointer
  // resting over the stage left one describing a square from the packing you just left.
  $('tip').style.display = 'none';
  // Newest call wins: two loads in flight could otherwise finish out of order and show the older one.
  const gen = S.gen = (S.gen || 0) + 1;
  const res = await fetch(dataURL(file));
  if (gen !== S.gen) return;
  if (!res.ok) {   // no exported analysis: say so everywhere, rather than leaving the last one up
    S.rec = null; S.file = file; S.n = n; $('nbox').value = n; $('title').textContent = 'No analysis for ' + file;
    $('sval').textContent = '—'; $('view').innerHTML = '';
    for (const id of ['sform', 'tags', 'prose', 'srcline', 'xlinks', 'rotsum', 'legend', 'freesum', 'freelist', 'gaps', 'stats', 'sel']) { const e = $(id); if (e) e.innerHTML = ''; }
    setText('panelerr', `${file} is in the catalogue but has no exported analysis — it failed the parser's checks, or the export skipped it.`);
    return;
  }
  const rec = await res.json(); if (gen !== S.gen) return;
  S.rec = rec; S.rec.s = +S.rec.s; S.file = file; S.n = n; S.sel = null; S.vb = null;
  $('nbox').value = n;
  const vs = $('variant'); vs.innerHTML = '';
  for (const f of variantsFor(n, file)) {
    const o = document.createElement('option'); o.value = f;
    o.textContent = f === (S.idx.records[n] || {}).svg ? `record (${f})` : f.replace('square-', '').replace('.svg', '');
    vs.appendChild(o);
  }
  vs.value = file; vs.style.display = vs.options.length > 1 ? '' : 'none';
  history.replaceState(null, '', `?p=${encodeURIComponent(file)}`);
  S.groups = mergedGroups(S.tol);
  render();
  // A throw here used to leave the panel half-built with no sign of why -- say so instead.
  try { setText('panelerr', ''); renderSide(); }
  catch (e) { console.error('renderSide failed', e); setText('panelerr', `This panel could not be built: ${e.message}. Reload with Ctrl-Shift-R; if it persists the page and its data are out of step.`); }
  renderLegend(); renderSel();
}
function gotoN(n) {
  const ns = recordNs(); if (!ns.length) return;
  if (!(Number.isInteger(n) && n >= 1)) {   // junk from the box or the URL: keep what is shown
    if (S.n != null) { $('nbox').value = S.n; return; }
    n = 17;
  }
  S.asked = n;   // reported in renderSide if the snap below moves us off it
  if (!S.idx.records[n]) {
    // Not in the main table, but the catalogue's sub-pages may still draw this n (up to n = 9465):
    // show the best of those before snapping anywhere.
    const own = variantsFor(n).sort((a, b) => +S.idx.files[a].s - +S.idx.files[b].s);
    if (own.length) { S.asked = null; load(own[0], n); return; }
    n = ns.reduce((best, x) => Math.abs(x - n) < Math.abs(best - n) ? x : best, ns[0]);
  }
  const rec = S.idx.records[n];
  const drawable = f => f && S.idx.files[f] && S.idx.files[f].n_parsed;
  const file = drawable(rec.svg) ? rec.svg : variantsFor(n)[0];
  if (file) load(file, n);
  else { S.rec = null; $('title').textContent = `no packing data for n = ${n}`; $('sval').textContent = '—';
         setText('shownfor', S.asked !== n ? `There is no catalogue entry for n = ${S.asked} either; the nearest is n = ${n}, and it has no drawing.` : ''); }
}
function initUI() {
  $('prev').onclick = () => { const ns = recordNs().filter(x => x < S.n); if (ns.length) gotoN(ns[ns.length - 1]); };
  $('next').onclick = () => { const ns = recordNs().filter(x => x > S.n); if (ns.length) gotoN(ns[0]); };
  $('nbox').onchange = () => gotoN(+$('nbox').value);
  $('variant').onchange = () => load($('variant').value, S.n);   // same n: S.asked still applies
  $('colormode').onchange = () => { S.colormode = $('colormode').value; recolor(); };
  $('tol').oninput = () => { const v = +$('tol').value; S.tol = v === 0 ? 0 : Math.pow(10, -5 + v / 100 * 5.7); $('tolv').textContent = S.tol === 0 ? '0' : S.tol < 0.01 ? S.tol.toExponential(1) : S.tol.toFixed(3); recolor(); };
  for (const id of ['L_free', 'L_contacts', 'L_gaps', 'L_sym', 'L_idx']) $(id).onchange = () => { drawFree(); drawContacts(); drawGaps(); drawSym(); drawLabels(); };
  initZoom();
}

async function main() {
  S.idx = await (await fetch('data/index.json')).json();
  initUI();
  const q = parseQuery();
  if (q.p) { S.asked = null; const f = S.idx.files[q.p]; load(q.p, f ? f.n : +(q.p.match(/\d+/) || [17])[0]); }
  else gotoN(+(q.n || 17));
}
main();
})();
