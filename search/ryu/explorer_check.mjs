// Float re-check of the explorer's drawings (site/www/k2c/constructions.js): for every series and k in range, the
// square count equals the Lean theorem's, side < k, every corner inside [0, side]^2 (tolerance 1e-7), and no two
// squares overlap (separating axis among the edge normals, tolerance 1e-7).  Exactness is Lean's; this checks that
// the page draws what the theorems say.
//
//   node search/ryu/explorer_check.mjs [KMAX]        (default 160)
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
const require = createRequire(import.meta.url);
const here = dirname(fileURLToPath(import.meta.url));
const www = join(here, '../../site/www/k2c');
const K2C = require(join(www, 'constructions.js'));
K2C.setPieces(JSON.parse(readFileSync(join(www, 'data/fang-pieces.json'), 'utf8')));
const KMAX = +(process.argv[2] || 160), EPS = 1e-7;   // display data: Fang's pieces are rounded to 9 decimals

const corners = ([x, y, c, s]) => [[.5, .5], [.5, -.5], [-.5, .5], [-.5, -.5]].map(([a, b]) => [x + c * a - s * b, y + s * a + c * b]);
function separated(p, q) {
  for (const [A, B] of [[p, q], [q, p]]) {
    const [x, y, c, s] = A, cb = corners(B);
    for (const [nx, ny] of [[c, s], [-s, c], [-c, -s], [s, -c]])
      if (cb.every(([X, Y]) => nx * (X - x) + ny * (Y - y) >= 0.5 - EPS)) return true;
  }
  return false;
}
function check(P) {
  const S = K2C.allSquares(P), side = P.side;
  if (S.length !== P.n || P.n !== P.claim) return `count ${S.length} / ${P.n} / claim ${P.claim}`;
  if (!(side < P.k)) return `side ${side}`;
  for (const q of S) for (const [X, Y] of corners(q)) if (X < -EPS || Y < -EPS || X > side + EPS || Y > side + EPS) return `wall ${q}`;
  const cell = new Map(), key = (i, j) => i * 100003 + j;
  S.forEach((q, i) => { const k = key(Math.floor(q[0] / 1.5), Math.floor(q[1] / 1.5)); (cell.get(k) || cell.set(k, []).get(k)).push(i); });
  let pairs = 0;
  for (let i = 0; i < S.length; i++) {
    const gi = Math.floor(S[i][0] / 1.5), gj = Math.floor(S[i][1] / 1.5);
    for (let di = -1; di <= 1; di++) for (let dj = -1; dj <= 1; dj++) for (const j of cell.get(key(gi + di, gj + dj)) || []) {
      if (j <= i) continue;
      const dx = S[i][0] - S[j][0], dy = S[i][1] - S[j][1];
      if (dx * dx + dy * dy >= 2 - EPS) continue;
      pairs++;
      if (!separated(S[i], S[j])) return `overlap ${S[i]} / ${S[j]}`;
    }
  }
  return null;
}
let bad = 0, total = 0;
for (const [name, ser] of Object.entries(K2C.SERIES)) {
  for (let k = ser.kmin; k <= KMAX; k++) {
    const P = ser.build(k), err = check(P);
    total++;
    if (k * k - P.n !== ser.bound(k)) { console.log(`${name} k=${k}: c = ${k * k - P.n} but bound() = ${ser.bound(k)}`); bad++; }
    if (err) { console.log(`${name} k=${k}: ${err}`); bad++; }
  }
  console.log(`${name}: k = ${ser.kmin}..${KMAX} checked`);
}
console.log(bad ? `NOT OK: ${bad} failures` : `OK: ${total} packings (count = Lean's, side < k, walls, pairwise disjoint)`);
process.exit(bad ? 1 : 0);
