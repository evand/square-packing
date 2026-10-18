/* The packings behind the Lean bounds on c*(k), square by square, in floats for drawing.
 *
 * Each builder follows its Lean proof's layout exactly (same parameters, same placement), so the drawing is the
 * packing the theorem asserts.  Exactness lives in Lean; search/ryu/explorer_check.mjs re-checks these float
 * versions (count, walls, pairwise disjoint up to 1e-9) for every k the page offers.
 *
 *   fang(k)    k >= 15  lean/Sqpack/FangFinal.lean theorem2 with the best p, q (FangTheorem: chain, Excess.stepR/R')
 *   ryu(k)     k >= 6   lean/Sqpack/RyuBound.lean packs_sub_8sqrt          (Ryu's Theorem 1.4)
 *   layers(k)  k >= 7   lean/Sqpack/RyuLayers.lean packs_sub_4sqrt2_layers (filled ends, two top layers)
 *
 * A packing is {k, n, side, blocks, sq, parts, params}: blocks are aligned arrays of unit squares {x, y, w, h, g}
 * (lower-left corner, w x h squares), sq are single squares [cx, cy, cos, sin, g], g a group name for colouring.
 */
(function (root) {
'use strict';
let PIECES = null;   // Fang's pieces: {name: {n, W, H, sq: [[x, y, t], ...]}} (data/fang-pieces.json)

const isqrt = n => { let r = Math.floor(Math.sqrt(n)); while (r * r > n) r--; while ((r + 1) * (r + 1) <= n) r++; return r; };
const count = P => P.sq.length + P.blocks.reduce((s, B) => s + B.w * B.h, 0);

/* ---------------------------------------------------------------- Ryu's band (lean/Sqpack/RyuBand.lean) */
// Rows j < R of b squares at angle theta (cos g, sin s); centre ctr(i, j) = (s, 0) + (a g - (j+1/2) s, a s + (j+1/2) g)
// with a = i + 1/2 + j s/g; the whole band lifted by y0 (RyuEnds.band_lift).
function band(b, R, g, s, y0, grp, out) {
  for (let j = 0; j < R; j++) for (let i = 0; i < b; i++) {
    const a = i + 0.5 + j * (s / g);
    out.push([s + a * g - (j + 0.5) * s, a * s + (j + 0.5) * g + y0, g, s, grp]);
  }
}
// Axis layers in the end regions of a strip [0, w] x [0, L] (RyuEnds.axis_layers): the bottom region's layer j holds
// n_j squares [w - n_j, w] x [j, j + 1]; the top region is a bottom region turned by pi (PacksIn.rot180), so its
// layer j is [0, n_j] x [L - j - 1, L - j].
function bottomLayers(w, ns, grp, out) { ns.forEach((n, j) => { for (let i = 0; i < n; i++) out.push([w - n + i + 0.5, j + 0.5, 1, 0, grp]); }); }
function topLayers(L, ns, grp, out) { ns.forEach((n, j) => { for (let i = 0; i < n; i++) out.push([i + 0.5, L - j - 0.5, 1, 0, grp]); }); }
// The L-shape (RyuBand.lshape): m x m grid, strip 1 (w x (m + w)) at x = m, strip 2 (w x m) transposed at y = m.
function lshape(m, w, s1, s2) {
  const sq = [];
  for (const [x, y, c, s, g] of s1) sq.push([m + x, y, c, s, g]);
  for (const [x, y, c, s, g] of s2) sq.push([y, m + x, s, c, g]);   // transpose: (x, y) -> (y, x), angle pi/2 - theta
  return { blocks: [{ x: 0, y: 0, w: m, h: m, g: 'grid' }], sq, side: m + w };
}
const trig = t => [(1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)];

/* Ryu's Theorem 1.4 (RyuBound.packs_sub_8b with b = Nat.sqrt(k - 5) + 1): rows R1 = k - 4 and R2 = k - b - 4,
   t = 1/b + 1/((b^2 + 1)(6k + 4b)); k^2 - 8b squares. */
function ryu(k) {
  if (k < 6) return null;
  const b = isqrt(k - 5) + 1, m = k - b;
  const t = 1 / b + 1 / ((b * b + 1) * (6 * k + 4 * b));
  const [g, s] = trig(t), w = b * g + s;
  const s1 = [], s2 = [];
  band(b, k - 4, g, s, 0, 'band', s1);
  band(b, k - b - 4, g, s, 0, 'band', s2);
  const P = Object.assign(lshape(m, w, s1, s2), { k, series: 'ryu' });
  P.n = count(P);
  P.parts = [['grid', m * m, `${m} × ${m} axis grid`], ['band', b * (2 * k - b - 8), `${k - 4} + ${k - b - 4} tilted rows of ${b}`]];
  P.params = { b, m, t, theta: 2 * Math.atan(t), rows: [k - 4, k - b - 4] };
  P.claim = k * k - 8 * b;
  return P;
}

/* The filled L-shape with two top layers (RyuLayers.packs_layers, b = Nat.sqrt(2k - 8) + 1): rows b + 1 + r and
   1 + r (r = k - b - 4), t = 1/b + 1/(b D (6k + 4b)), D = b^2 + 1; ends: floor((b-1)/2) under each band,
   floor(b/2) over the right band, floor(b/2) + 1 and then 1 over the top band; k^2 - 4b squares. */
function layers(k) {
  if (k < 7) return null;
  const b = isqrt(2 * k - 8) + 1, m = k - b, r = k - b - 4;
  const D = b * b + 1, t = 1 / b + 1 / (b * D * (6 * k + 4 * b));
  const [g, s] = trig(t), w = b * g + s;
  const H = R => (R - 1) / g + b * s + g;     // a band of R rows is this long
  const nb = Math.floor((b - 1) / 2), nh = Math.floor(b / 2);
  const L1 = m + w, L2 = m, R1 = b + 1 + r, R2 = 1 + r;
  const s1 = [], s2 = [];
  band(b, R1, g, s, 0, 'band', s1); bottomLayers(w, [nb], 'end', s1); topLayers(L1, [nh], 'end', s1);
  band(b, R2, g, s, 0, 'band', s2); bottomLayers(w, [nb], 'end', s2); topLayers(L2, [nh + 1, 1], 'end', s2);
  const P = Object.assign(lshape(m, w, s1, s2), { k, series: 'layers' });
  P.n = count(P);
  P.parts = [['grid', m * m, `${m} × ${m} axis grid`], ['band', b * (R1 + R2), `${R1} + ${R2} tilted rows of ${b}`],
    ['end', 2 * nb + 2 * nh + 2, `band ends: ${nb} + ${nh} (right band), ${nb} + ${nh + 1} + 1 (top band)`]];
  P.params = { b, m, t, theta: 2 * Math.atan(t), rows: [R1, R2], ends: [nb, nh, nb, nh + 1, 1], slack: [L1 - H(R1), L2 - H(R2)] };
  P.claim = k * k - 4 * b;
  return P;
}

/* ---------------------------------------------------------------- Fang (lean/Sqpack/FangTheorem.lean) */
// A squeezable piece: N squares in (a - d) x (b - d).  From a certificate: square side S < m (m = ceil S), or a
// rectangle W x H in a x b (squeezable_of_packsIn: d = min(a - W, b - H)).
// keep: Lean's PacksIn.drop keeps the first N' squares (Fang's 308 in side < 18 serves as 307, FangFinal.base).
function piece(name, a, bb, grp, keep) {
  const C = PIECES[name];
  const sq = C.sq.slice(0, keep == null ? C.n : keep).map(([x, y, t]) => { const [c, s] = trig(t); return [x, y, c, s, grp]; });
  return { a, b: bb, d: Math.min(a - C.W, bb - C.H), blocks: [], sq, n: sq.length, tree: [name] };
}
// Squeezable.glue (FangGlue.glue): A (a x b) top left at (0, b'), B (a' x b') bottom right at (n - a', 0), aligned
// blocks D = (n - a) x b at (a - d, b' - d) and C = (n - a') x b' at (0, 0); n = b + b', d = min(dA, dB, 1).
function glue(A, B) {
  const n = A.b + B.b, d = Math.min(A.d, B.d, 1);
  const shift = (P, dx, dy) => ({ blocks: P.blocks.map(Q => Object.assign({}, Q, { x: Q.x + dx, y: Q.y + dy })), sq: P.sq.map(([x, y, c, s, g]) => [x + dx, y + dy, c, s, g]) });
  const PA = shift(A, 0, B.b), PB = shift(B, n - B.a, 0);
  const blocks = PA.blocks.concat(PB.blocks);
  if (n - A.a > 0) blocks.push({ x: A.a - d, y: B.b - d, w: n - A.a, h: A.b, g: 'block' });
  if (n - B.a > 0) blocks.push({ x: 0, y: 0, w: n - B.a, h: B.b, g: 'block' });
  return { a: n, b: n, d, blocks, sq: PA.sq.concat(PB.sq), tree: B.tree.concat(A.tree) };
}
const R257 = () => piece('R257', 17, 16, 'R257'), R443 = () => piece('R443', 22, 21, 'R443');
const base1 = r => r === 18 ? piece('P308', 18, 18, 'base', 307) : r <= 26 ? piece('P' + (r * r - r + 1), r, r, 'base') : glue(R257(), piece('P' + ((r - 16) * (r - 16) - (r - 16)), r - 16, r - 16, 'base'));

/* Theorem 2 with the best p, q: n = r + 16 p' + 21 q, r = 15 + t mod 16 in [15, 30], p' = p + t div 16 (here t div 16
   = 0, since p, q are optimal); base r (r <= 26: Fang's or a register piece; else R glued onto m^2 - m, m = r - 16),
   then p' copies of R, then q of R'.  Count n^2 - n + 1 + p + 2q. */
function fangPQ(k) {
  let best = null;
  for (let q = 0; 21 * q <= k - 15; q++) { const p = Math.floor((k - 15 - 21 * q) / 16); if (!best || p + 2 * q > best[0] + 2 * best[1]) best = [p, q]; }
  return best;
}
function fang(k) {
  if (k < 15 || !PIECES) return null;
  const [p, q] = fangPQ(k), t = k - 15 - 16 * p - 21 * q, r = 15 + (t % 16), pp = p + Math.floor(t / 16);
  let X = base1(r);
  for (let i = 0; i < pp; i++) X = glue(R257(), X);
  for (let i = 0; i < q; i++) X = glue(R443(), X);
  const P = { k, series: 'fang', blocks: X.blocks, sq: X.sq, side: k - X.d };
  P.n = count(P);
  const by = g => P.sq.filter(z => z[4] === g).length;
  const blk = P.blocks.reduce((s, B) => s + B.w * B.h, 0);
  P.parts = [['base', by('base'), r <= 26 ? `base piece: ${r * r - r + 1} in side < ${r}` : `base: ${(r - 16) * (r - 16) - (r - 16)} in side < ${r - 16}, glued to R`]];
  if (by('R257')) P.parts.push(['R257', by('R257'), `${pp + (r > 26 ? 1 : 0)} × R (257 in 17 × 16)`]);
  if (by('R443')) P.parts.push(['R443', by('R443'), `${q} × R′ (443 in 22 × 21)`]);
  P.parts.push(['block', blk, 'aligned blocks between the pieces']);
  P.params = { p, q, r, pp, d: X.d, tree: X.tree };
  P.claim = k * k - k + 1 + p + 2 * q;
  return P;
}

const SERIES = {
  fang: { kmin: 15, build: fang, bound: k => { const [p, q] = fangPQ(k); return k - 1 - p - 2 * q; } },
  ryu: { kmin: 6, build: ryu, bound: k => 8 * (isqrt(k - 5) + 1) },
  layers: { kmin: 7, build: layers, bound: k => 4 * (isqrt(2 * k - 8) + 1) },
};
// bound(k) = c: the series packs k^2 - c squares in side < k, so c*(k) <= c - 1.

// Every unit square of a packing, blocks expanded: [cx, cy, cos, sin, g].
function allSquares(P) {
  const out = P.sq.slice();
  for (const B of P.blocks) for (let i = 0; i < B.w; i++) for (let j = 0; j < B.h; j++) out.push([B.x + i + 0.5, B.y + j + 0.5, 1, 0, B.g]);
  return out;
}

const API = { setPieces: p => { PIECES = p; }, fang, ryu, layers, SERIES, allSquares, count, isqrt, fangPQ };
if (typeof module !== 'undefined' && module.exports) module.exports = API; else root.K2C = API;
})(typeof window !== 'undefined' ? window : globalThis);
