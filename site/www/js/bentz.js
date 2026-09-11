/* Bentz 2010 (s(13) = 4): clickable walk through the case tree.
   Renders into #bentz: a .tree column on the left, an SVG of the 4x4 container on the right,
   plus a caption carrying the node text and the budget line. Dependency-free. */
(() => {
const NS = 'http://www.w3.org/2000/svg';
const el = (t, a = {}, p = null) => { const e = document.createElementNS(NS, t); for (const [k, v] of Object.entries(a)) e.setAttribute(k, v); if (p) p.appendChild(e); return e; };
const h = (t, a = {}, p = null) => { const e = document.createElement(t); for (const [k, v] of Object.entries(a)) { if (k === 'text') e.textContent = v; else e.setAttribute(k, v); } if (p) p.appendChild(e); return e; };
const key = p => p[0].toFixed(4) + ',' + p[1].toFixed(4);
// text with a trailing "_x" rendered as a subscript
const txt = (p, x, y, s, size, fill, anchor) => {
  const t = el('text', { x, y, 'font-size': size, fill, 'font-family': 'IBM Plex Mono, monospace', ...(anchor ? { 'text-anchor': anchor } : {}) }, p);
  const i = s.indexOf('_');
  if (i < 0) { t.textContent = s; return t; }
  t.appendChild(document.createTextNode(s.slice(0, i)));
  el('tspan', { 'font-size': size * 0.72, dy: size * 0.22 }, t).textContent = s.slice(i + 1);
  return t;
};
const setOf = ps => new Set((ps || []).map(key));

const CSS = `
#bentz{display:grid;grid-template-columns:230px minmax(0,1fr);gap:18px;align-items:start;margin:4px 0}
@media (max-width:760px){#bentz{grid-template-columns:1fr}}
#bentz .tree{display:flex;flex-direction:column;gap:1px;font-size:.86rem;border:1px solid var(--rule);border-radius:3px;background:var(--surface);padding:6px}
#bentz .tree button{display:flex;gap:8px;align-items:baseline;width:100%;text-align:left;background:none;border:none;border-radius:3px;padding:5px 7px;line-height:1.25;color:var(--ink)}
#bentz .tree button:hover{background:var(--accent-soft)}
#bentz .tree button.on{background:var(--accent-soft);color:var(--accent);font-weight:600}
#bentz .tree button .b{margin-left:auto;font-family:"IBM Plex Mono",monospace;font-size:.72rem;color:var(--muted);flex:none}
#bentz .tree button.on .b{color:var(--accent)}
#bentz .tree .d1{padding-left:18px}
#bentz .tree .d2{padding-left:32px}
#bentz .stage2{display:flex;flex-direction:column;gap:10px;min-width:0}
#bentz svg{width:100%;aspect-ratio:1;background:var(--surface);border:1px solid var(--rule);border-radius:3px;display:block}
#bentz .cap{display:flex;flex-direction:column;gap:8px;font-size:.92rem}
#bentz .cap h3{font-size:1rem}
#bentz .cap .ref{font-family:"IBM Plex Mono",monospace;font-size:.72rem;color:var(--muted);letter-spacing:.02em}
#bentz .cap .nb{font-size:.84rem;color:var(--muted);border-left:2px solid var(--rule);padding-left:10px}
#bentz .budget{font-family:"IBM Plex Mono",monospace;font-size:.82rem;padding:6px 10px;border-radius:3px;background:var(--gold-soft);color:var(--gold);align-self:flex-start}
#bentz .budget.done{background:var(--green-soft);color:var(--green)}
#bentz .keys2{display:flex;flex-wrap:wrap;gap:12px;font-size:.78rem;color:var(--muted);font-family:"IBM Plex Mono",monospace}
#bentz .keys2 i{display:inline-block;width:10px;height:10px;border-radius:50%;vertical-align:-1px;margin-right:5px}
#bentz .keys2 i.sq{border-radius:2px}
`;

let D = null, byId = {}, order = [], cur = null, svg = null, cap = null, tree = null;

function build(host) {
  host.textContent = '';
  tree = h('div', { class: 'tree' }, host);
  const stage = h('div', { class: 'stage2' }, host);
  svg = el('svg', { viewBox: '-0.36 -0.36 4.72 4.72' }, stage);
  cap = h('div', { class: 'cap' }, stage);
  const keys = h('div', { class: 'keys2' }, stage);
  keys.innerHTML = '<span><i style="background:var(--red)"></i>new point</span>' +
    '<span><i style="background:var(--accent)"></i>free point</span>' +
    '<span><i style="background:var(--muted);opacity:.6"></i>spent (inside a pinned box)</span>' +
    '<span><i class="sq" style="background:var(--accent-soft);border:1px solid var(--accent)"></i>pinned box</span>' +
    '<span><i class="sq" style="background:var(--gold-soft);border:1px solid var(--gold)"></i>region a lemma misses</span>';
  order.forEach(n => {
    const b = h('button', { class: 'd' + (n.depth || 0), type: 'button' }, tree);
    h('span', { text: n.title }, b);
    if (n.budget && n.budget.placed) h('span', { class: 'b', text: (n.budget.free + n.budget.placed) }, b);
    b.onclick = () => show(n.id);
    n._btn = b;
  });
}

function show(id) {
  cur = id; const n = byId[id];
  order.forEach(m => m._btn.classList.toggle('on', m.id === id));
  draw(n);
  cap.textContent = '';
  h('h3', { text: n.title }, cap);
  h('div', { class: 'ref', text: n.ref }, cap);
  h('p', { text: n.text }, cap);
  if (n.note) h('p', { class: 'nb', text: n.note }, cap);
  const done = n.budget && n.budget.placed && (n.budget.free + n.budget.placed) < 13;
  h('div', { class: 'budget' + (done ? ' done' : ''), text: n.budget.line }, cap);
}

function draw(n) {
  svg.textContent = '';
  const w = el('g', { transform: 'matrix(1 0 0 -1 0 4)' }, svg);   // world: y up
  const lab = el('g', {}, svg);                                     // labels: (x, 4-y)
  const T = p => [p[0], 4 - p[1]];
  const path = poly => poly.map(p => p.join(' ')).join(' ');

  for (let i = 1; i < 4; i++) {
    el('line', { x1: i, y1: 0, x2: i, y2: 4, stroke: 'var(--rule)', 'stroke-width': 0.006 }, w);
    el('line', { x1: 0, y1: i, x2: 4, y2: i, stroke: 'var(--rule)', 'stroke-width': 0.006 }, w);
  }
  el('rect', { x: 0, y: 0, width: 4, height: 4, fill: 'none', stroke: 'var(--wall)', 'stroke-width': 0.02 }, w);

  (n.regions || []).forEach(r => {
    el('polygon', { points: path(r.poly), fill: 'var(--gold)', 'fill-opacity': 0.14, stroke: 'var(--gold)', 'stroke-width': 0.012, 'stroke-dasharray': '0.06 0.05' }, w);
    if (r.label) {
      const c = r.poly.reduce((a, p) => [a[0] + p[0] / r.poly.length, a[1] + p[1] / r.poly.length], [0, 0]);
      txt(lab, T(c)[0], T(c)[1] + 0.05, r.label, 0.14, 'var(--gold)', 'middle');
    }
  });
  (n.boxes || []).forEach(b => {
    const region = b.kind === 'region';
    el('polygon', {
      points: path(b.poly), fill: region ? 'var(--red)' : 'var(--accent)', 'fill-opacity': region ? 0.10 : 0.16,
      stroke: region ? 'var(--red)' : 'var(--accent)', 'stroke-width': 0.016,
      ...(region ? { 'stroke-dasharray': '0.08 0.05' } : {})
    }, w);
    const c = b.label_at || b.poly.reduce((a, p) => [a[0] + p[0] / b.poly.length, a[1] + p[1] / b.poly.length], [0, 0]);
    txt(lab, T(c)[0], T(c)[1] + 0.05, b.label, 0.15, region ? 'var(--red)' : 'var(--accent)', 'middle');
  });
  (n.segments || []).forEach(s => el('line', { x1: s.a[0], y1: s.a[1], x2: s.b[0], y2: s.b[1], stroke: 'var(--red)', 'stroke-width': 0.075, 'stroke-linecap': 'round', opacity: 0.32 }, w));

  const spent = setOf(n.spent), emph = setOf(n.emph);
  const par = n.parent ? byId[n.parent] : null;
  const seen = par ? setOf(par.points) : new Set();
  (n.alts || []).forEach(p => el('circle', { cx: p[0], cy: p[1], r: 0.05, fill: 'none', stroke: 'var(--muted)', 'stroke-width': 0.014, 'stroke-dasharray': '0.04 0.035' }, w));
  n.points.forEach(p => {
    const k = key(p), isSpent = spent.has(k), isNew = par && !seen.has(k);
    if (emph.has(k)) el('circle', { cx: p[0], cy: p[1], r: 0.105, fill: 'none', stroke: 'var(--gold)', 'stroke-width': 0.018 }, w);
    if (isNew && !isSpent) el('circle', { cx: p[0], cy: p[1], r: 0.1, fill: 'var(--red)', 'fill-opacity': 0.16 }, w);
    el('circle', {
      cx: p[0], cy: p[1], r: isSpent ? 0.042 : isNew ? 0.056 : 0.05,
      fill: isSpent ? 'var(--muted)' : isNew ? 'var(--red)' : 'var(--accent)',
      'fill-opacity': isSpent ? 0.55 : 1
    }, w);
  });
  (n.labels || []).forEach(L => {
    const q = T(L.p);
    txt(lab, q[0] + (L.dx !== undefined ? L.dx : 0.1), q[1] - (L.dy !== undefined ? L.dy : 0.08), L.t, 0.15, 'var(--ink)');
  });
}

async function init() {
  const host = document.getElementById('bentz');
  if (!host) return;
  if (!document.getElementById('bentz-css')) { const s = h('style', { id: 'bentz-css' }, document.head); s.textContent = CSS; }
  host.textContent = 'loading the case tree…';
  let d;
  try { d = await (await fetch('data/bentz13.json')).json(); }
  catch (e) { host.textContent = 'Could not load data/bentz13.json.'; return; }
  D = d; byId = {}; D.nodes.forEach(n => byId[n.id] = n);
  order = [];
  const walk = p => D.nodes.filter(n => n.parent === p).forEach(n => { order.push(n); walk(n.id); });
  walk(null);
  build(host);
  show(order[0].id);
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
