/* Unavoidable-points playground: drag/rotate a unit square over a certificate's points; count covered points;
   draw the count landscape at the current angle. */
(() => {
const $ = id => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
const el = (t, a = {}, p = null) => { const e = document.createElementNS(NS, t); for (const [k, v] of Object.entries(a)) e.setAttribute(k, v); if (p) p.appendChild(e); return e; };
let C = null;                 // certificate: {s, pts:[[x,y,w]], k}
let pose = { x: 1.5, y: 1.5, th: 0 };
const EPS = 1e-9;

async function loadCert(name) {
  const txt = await (await fetch('data/' + name)).text();
  const nums = txt.split(/\s+/).filter(Boolean).map(Number);
  const [sn, sd, D, W, m] = nums; const pts = [];
  for (let i = 0; i < m; i++) pts.push([nums[5 + 3 * i] / D, nums[6 + 3 * i] / D, nums[7 + 3 * i] / W]);
  const k = Math.round(1 / Math.min(...pts.map(p => p[2])));   // uniform weight 1/k
  C = { s: sn / sd, sTex: `${sn}/${sd}`, pts, k, m };
  pose = { x: C.s / 2, y: C.s / 2, th: 0 }; $('ang').value = 0;
  build(); update();
}

function covered(x, y, thDeg) {
  const t = thDeg * Math.PI / 180, c = Math.cos(t), s = Math.sin(t); const out = [];
  for (let i = 0; i < C.pts.length; i++) {
    const dx = C.pts[i][0] - x, dy = C.pts[i][1] - y;
    const u = c * dx + s * dy, v = -s * dx + c * dy;
    if (Math.abs(u) <= 0.5 + EPS && Math.abs(v) <= 0.5 + EPS) out.push(i);
  }
  return out;
}
function halfExtent(thDeg) { const t = thDeg * Math.PI / 180; return (Math.abs(Math.cos(t)) + Math.abs(Math.sin(t))) / 2; }
function clamp() { const h = halfExtent(pose.th); pose.x = Math.min(Math.max(pose.x, h), C.s - h); pose.y = Math.min(Math.max(pose.y, h), C.s - h); }

function build() {
  const svg = $('pl'); svg.innerHTML = ''; const s = C.s; svg.setAttribute('viewBox', `${-0.04 * s} ${-0.04 * s} ${1.08 * s} ${1.08 * s}`);
  const w = el('g', { transform: `matrix(1 0 0 -1 0 ${s})` }, svg);
  el('rect', { x: 0, y: 0, width: s, height: s, fill: 'none', stroke: 'var(--wall)', 'stroke-width': s / 250 }, w);
  // faint unit grid to give scale
  for (let i = 1; i < s; i++) { el('line', { x1: i, y1: 0, x2: i, y2: s, stroke: 'var(--rule)', 'stroke-width': s / 800 }, w); el('line', { x1: 0, y1: i, x2: s, y2: i, stroke: 'var(--rule)', 'stroke-width': s / 800 }, w); }
  const pg = el('g', { id: 'pts' }, w);
  C.pts.forEach((p, i) => el('circle', { cx: p[0], cy: p[1], r: s / 160, fill: 'var(--muted)', 'data-i': i }, pg));
  el('rect', { id: 'sq', x: -0.5, y: -0.5, width: 1, height: 1, fill: 'rgba(14,116,144,0.18)', stroke: 'var(--accent)', 'stroke-width': s / 250 }, w);
  const land = $('land'); land.width = 120; land.height = 120;
}
function update() {
  clamp();
  const cov = covered(pose.x, pose.y, pose.th); const set = new Set(cov);
  $('sq').setAttribute('transform', `translate(${pose.x} ${pose.y}) rotate(${pose.th})`);
  document.querySelectorAll('#pts circle').forEach(c => { const on = set.has(+c.dataset.i); c.setAttribute('fill', on ? 'var(--red)' : 'var(--muted)'); c.setAttribute('r', on ? C.s / 110 : C.s / 160); });
  $('cnt').textContent = cov.length; $('cnt').className = 'count' + (cov.length <= C.k ? ' low' : '');
  $('cntnote').textContent = cov.length === C.k ? `exactly ${C.k} — a tight spot` : `needs at least ${C.k}`;
  $('angv').textContent = (+$('ang').value).toFixed(1) + '°';
  if ($('showland').checked) landscape(); else $('land').getContext('2d').clearRect(0, 0, 120, 120);
}
function landscape() {
  const cv = $('land'), ctx = cv.getContext('2d'), N = cv.width, s = C.s, h = halfExtent(pose.th);
  const img = ctx.createImageData(N, N); let min = Infinity, tight = 0;
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    const x = (i + 0.5) / N * s, y = (j + 0.5) / N * s; const o = 4 * ((N - 1 - j) * N + i);
    if (x < h || x > s - h || y < h || y > s - h) { img.data[o + 3] = 0; continue; }
    const c = covered(x, y, pose.th).length; if (c < min) min = c; if (c === C.k) tight++;
    let col = c <= C.k ? [180, 50, 42] : c === C.k + 1 ? [245, 235, 218] : c <= C.k + 3 ? [220, 237, 241] : [14, 116, 144];
    if (c < C.k) col = [0, 0, 0];
    img.data[o] = col[0]; img.data[o + 1] = col[1]; img.data[o + 2] = col[2]; img.data[o + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);
  // marker for the current pose
  ctx.fillStyle = '#131A1E'; ctx.fillRect(Math.round(pose.x / s * N) - 1, N - 1 - Math.round(pose.y / s * N) - 1, 3, 3);
  $('landnote').textContent = `at ${(+$('ang').value).toFixed(1)}°: minimum ${min} over ${N}×${N} centre positions; ${(100 * tight / (N * N)).toFixed(1)}% of positions are tight`;
}
function scanAll() {
  let gmin = Infinity, worst = null; const N = 90, s = C.s;
  for (let a = 0; a <= 90; a += 0.5) {
    const h = halfExtent(a);
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      const x = (i + 0.5) / N * s, y = (j + 0.5) / N * s; if (x < h || x > s - h || y < h || y > s - h) continue;
      const c = covered(x, y, a).length; if (c < gmin) { gmin = c; worst = { x, y, a }; }
    }
  }
  $('scannote').textContent = `181 angles × ${N}×${N} positions: minimum ${gmin} (first found at ${worst.a}°, centre ${worst.x.toFixed(2)}, ${worst.y.toFixed(2)}). The exact verifier confirms ${C.k} over the whole continuum.`;
  pose = { x: worst.x, y: worst.y, th: worst.a }; $('ang').value = worst.a; update();
}

function initDrag() {
  const svg = $('pl'); let drag = null;
  const pt = e => { const p = svg.createSVGPoint(); p.x = e.clientX; p.y = e.clientY; const q = p.matrixTransform(svg.getScreenCTM().inverse()); return { x: q.x, y: C.s - q.y }; };
  svg.addEventListener('pointerdown', e => { const p = pt(e); drag = { dx: pose.x - p.x, dy: pose.y - p.y }; svg.setPointerCapture(e.pointerId); if (Math.hypot(pose.x - p.x, pose.y - p.y) > 1.2) { pose.x = p.x; pose.y = p.y; drag = { dx: 0, dy: 0 }; } update(); });
  svg.addEventListener('pointermove', e => { if (!drag) return; const p = pt(e); pose.x = p.x + drag.dx; pose.y = p.y + drag.dy; update(); });
  svg.addEventListener('pointerup', () => { drag = null; });
  svg.addEventListener('wheel', e => { e.preventDefault(); pose.th = ((pose.th + (e.deltaY > 0 ? 1 : -1) * 0.5) % 90 + 90) % 90; $('ang').value = pose.th; update(); }, { passive: false });
  window.addEventListener('keydown', e => { const st = { ArrowLeft: [-0.01, 0], ArrowRight: [0.01, 0], ArrowUp: [0, 0.01], ArrowDown: [0, -0.01] }[e.key]; if (st && document.activeElement.tagName !== 'INPUT') { pose.x += st[0]; pose.y += st[1]; update(); e.preventDefault(); } });
  $('ang').oninput = () => { pose.th = +$('ang').value; update(); };
  $('showland').onchange = update; $('scan').onclick = scanAll; $('cert').onchange = () => loadCert($('cert').value);
}
initDrag(); loadCert('cert_81.txt');
})();
