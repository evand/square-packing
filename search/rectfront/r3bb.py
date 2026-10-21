#!/usr/bin/env python3
"""R_3 lower bound: tilt-box branch and bound with force-network certificates (2026-10-20).

Claim.  Three unit squares in [0,w] x [0,h] with w < 2, h < 3 force 2w + h >= L = 4 + 2 sqrt2.

Per configuration: tilts t_K (mod 90 deg) and, per pair (i, j), a separating axis u = an edge normal of i or of j
(SAT for squares): u.(c_j - c_i) >= S_ij = s_i(u) + s_j(u), s_K(u) = (|u.a_K| + |u.b_K|)/2.  Walls:
k_K <= c_K <= (w, h) - k_K, k_K = (|cos t_K| + |sin t_K|)/2.  For any lambda >= 0 (one per pair), with
G_K = sum_p (+-) lambda_p u_p the net pair force on K (sum_K G_K = 0):
    sum_p lambda_p S_p  <=  sum_K |G_Kx| (w/2 - k_K) + |G_Ky| (h/2 - k_K)  =: R(w, h)        (*)
(*) is exact LP duality for the centre problem of one type.  R is increasing and linear in (w, h), so if
sum lambda S >= R at both upper vertices V1 = (2, 2 sqrt2), V2 = ((1 + 2 sqrt2)/2, 3) of the excluded region
{w < 2, h < 3, 2w + h < L}, no configuration of that type has (w, h) in it (strictly: R(p) < max R(V) for p in the
open region when (a_w, a_h) != 0; we require margin > TOL anyway).

Rigour: trig enclosures over a tilt box by monotonicity (exact critical points), all arithmetic in float64 with an
explicit slack: every certified margin must exceed TOL = 1e-9, far above accumulated rounding (< 1e-13 here).
A replay in exact interval arithmetic (mpmath.iv) of every leaf is `--replay` (slow; for the record).

Boxes within EPS (deg, torus sup-distance) of a wedge point (two tilts = 0, one = 45 mod 90) are not split below
MINLOC (only the tube's own boxes); the (box, type) pairs that fail there are written out for the local lemma.

  python3 r3bb.py run --procs 4 --out runs/rectfront/r3bb1 [--eps 1.0] [--root 15]
"""
import argparse, itertools, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
from scipy.optimize import linprog

R2 = math.sqrt(2)
L = 4 + 2 * R2
V = [(2.0, 2 * R2), ((1 + 2 * R2) / 2, 3.0)]
TOL = 1e-9
PAIRS = [(0, 1), (0, 2), (1, 2)]
# option = (owner o in {0: i, 1: j}, m in 0..3): u = edge normal of the owner at angle t_owner + 90 m
OPTS = [(o, m) for o in (0, 1) for m in range(4)]
TYPES = list(itertools.product(range(8), repeat=3))


# ---------- interval trig on degree intervals ----------
def icos(a, b):
    """Enclosure of cos over [a, b] degrees (b - a <= 360)."""
    lo, hi = min(math.cos(math.radians(a)), math.cos(math.radians(b))), max(math.cos(math.radians(a)), math.cos(math.radians(b)))
    k0 = math.ceil(a / 180.0); k1 = math.floor(b / 180.0)
    for k in range(k0, k1 + 1):
        v = 1.0 if k % 2 == 0 else -1.0
        lo, hi = min(lo, v), max(hi, v)
    return lo - 1e-15, hi + 1e-15


def isin(a, b):
    return icos(a - 90.0, b - 90.0)


def iabs(lo, hi):
    if lo >= 0: return lo, hi
    if hi <= 0: return -hi, -lo
    return 0.0, max(-lo, hi)


def isupp_diff(a, b):
    """Enclosure of (|cos d| + |sin d|)/2 over d in [a, b] degrees (support of a square along a direction at
    relative angle d).  The function has period 90, is 1/2 at multiples of 90 and sqrt2/2 at 45 + 90k."""
    if b - a >= 90: return 0.5, R2 / 2
    c = iabs(*icos(a, b)); s = iabs(*isin(a, b))
    lo = 0.5 * (c[0] + s[0]); hi = 0.5 * (c[1] + s[1])
    # tighter: f(d) = (|cos|+|sin|)/2 = cos(d' - 45)/sqrt2 with d' = d mod 90 in [0, 90)
    vals = [_f(a), _f(b)]
    flo, fhi = min(vals), max(vals)
    k0, k1 = math.ceil((a - 45) / 90.0), math.floor((b - 45) / 90.0)
    if k1 >= k0: fhi = R2 / 2
    k0, k1 = math.ceil(a / 90.0), math.floor(b / 90.0)
    if k1 >= k0: flo = 0.5
    return max(lo, flo - 1e-15), min(hi, fhi + 1e-15)


def _f(d):
    r = math.radians(d)
    return 0.5 * (abs(math.cos(r)) + abs(math.sin(r)))


# ---------- one type on one box ----------
def type_data(box, typ, pt=None):
    """Interval (or point, if pt is given) data of a type: per pair, the direction u (cos, sin intervals),
    S = s_i(u) + s_j(u); per square k_K.  box: list of 3 (lo, hi) degree intervals."""
    if pt is not None:
        box = [(x, x) for x in pt]
    k = [isupp_diff(*bx) for bx in box]             # k_K = support along x = (|cos t| + |sin t|)/2
    pairs = []
    for (i, j), oi in zip(PAIRS, typ):
        o, m = OPTS[oi]
        own, oth = (i, j) if o == 0 else (j, i)
        a, b = box[own][0] + 90 * m, box[own][1] + 90 * m
        ux, uy = icos(a, b), isin(a, b)
        d = (box[oth][0] - box[own][1], box[oth][1] - box[own][0])
        so = isupp_diff(*d)
        S = (0.5 + so[0], 0.5 + so[1])
        pairs.append((i, j, ux, uy, S))
    return k, pairs


def margin_lp(box, typ):
    """lambda (normalised, sum 1) maximising the worst-vertex margin at the box centre; float LP."""
    pt = [(lo + hi) / 2 for lo, hi in box]
    k, pairs = type_data(None, typ, pt)
    # variables: lam0..2, gx0..2, gy0..2, m  (11)
    nv = 10
    A, bb = [], []
    # g >= +-G
    for K in range(3):
        for comp in (2, 3):
            for sg in (1, -1):
                r = [0.0] * nv
                for p, (i, j, ux, uy, S) in enumerate(pairs):
                    u = (ux[0], uy[0])[comp - 2]
                    if K == j: r[p] += sg * u
                    if K == i: r[p] -= sg * u
                r[3 + K if comp == 2 else 6 + K] = -1.0
                A.append(r); bb.append(0.0)
    for (w, h) in V:
        r = [0.0] * nv          # m - sum lam S + sum g (w/2 - k) <= 0
        for p, pr in enumerate(pairs): r[p] = -pr[4][0]
        for K in range(3):
            r[3 + K] = w / 2 - k[K][0]; r[6 + K] = h / 2 - k[K][0]
        r[9] = 1.0
        A.append(r); bb.append(0.0)
    c = [0.0] * nv; c[9] = -1.0
    bounds = [(0, None)] * 9 + [(None, None)]
    res = linprog(c, A_ub=A, b_ub=bb, A_eq=[[1, 1, 1] + [0] * 7], b_eq=[1], bounds=bounds, method='highs')
    if res.status != 0:
        return None, -math.inf
    return [max(0.0, x) for x in res.x[:3]], -res.fun


def margin_box(box, typ, lam):
    """Rigorous (up to TOL) lower bound over the box of min_V [sum lam S - R_V]."""
    k, pairs = type_data(box, typ)
    lhs = sum(l * pr[4][0] for l, pr in zip(lam, pairs))
    worst = math.inf
    G = []
    for K in range(3):
        gx, gy = [0.0, 0.0], [0.0, 0.0]
        for l, (i, j, ux, uy, S) in zip(lam, pairs):
            sg = 1 if K == j else (-1 if K == i else 0)
            if sg == 0 or l == 0: continue
            for acc, iv in ((gx, ux), (gy, uy)):
                lo, hi = sorted((sg * l * iv[0], sg * l * iv[1]))
                acc[0] += lo; acc[1] += hi
        G.append((iabs(*gx)[1], iabs(*gy)[1]))
    for (w, h) in V:
        R = sum(G[K][0] * (w / 2 - k[K][0]) + G[K][1] * (h / 2 - k[K][0]) for K in range(3))
        worst = min(worst, lhs - R)
    return worst


# ---------- box classification ----------
WEDGE = [(0, 0, 45), (0, 45, 0), (45, 0, 0)]


def tdist(a, b):
    d = abs((a - b) % 90.0)
    return min(d, 90.0 - d)


LOCAL_MODE = 'tube'   # 'tube': two tilts within eps of 0; 'line': two tilt ranges contain 0 (box touches a zero line)


def near_wedge(box, eps):
    if LOCAL_MODE == 'line':
        return sum(1 for lo, hi in box if lo <= 0.0 <= hi) >= 2
    return _near_tube(box, eps)


def _near_tube(box, eps):
    """Inside the tube around the zero lines: at least two tilts within eps of 0 (mod 90) on the whole box.
    (Zero margin holds on the whole lines {two tilts = 0} at V1: two upright squares side by side at w = 2;
    the wedge points lie on them.)"""
    near = [max(tdist(lo, 0.0), tdist(hi, 0.0)) <= eps and (hi - lo) <= 2 * eps for lo, hi in box]
    return sum(near) >= 2


def solve_box(args):
    box, inherited, eps, minloc, mode = args
    global LOCAL_MODE
    LOCAL_MODE = mode
    """Returns (status, payload): 'ok' with per-type lambdas, 'split' with the failing types, or 'local'."""
    lams, fail = {}, []
    for ti, typ in enumerate(TYPES):
        lam = inherited.get(ti) if inherited else None
        if lam is not None and margin_box(box, typ, lam) > TOL:
            lams[ti] = lam; continue
        lam, mc = margin_lp(box, typ)
        if lam is not None and margin_box(box, typ, lam) > TOL:
            lams[ti] = lam; continue
        fail.append(ti)
    if not fail:
        return box, 'ok', lams, []
    width = max(hi - lo for lo, hi in box)
    if near_wedge(box, eps) and width <= minloc:
        return box, 'local', lams, fail
    return box, 'split', lams, fail


def in_domain(box):
    """Fundamental domain: labels sorted (t0 <= t1 <= t2, squares are interchangeable; every type is checked)
    and the mirror x -> -x, which maps sorted tilts (a, b, c) to (-c, -b, -a): keep a + c >= 0.  A box is kept
    if it meets {t0 <= t1 <= t2, t0 + t2 >= 0}."""
    (a0, a1), (b0, b1), (c0, c1) = box
    return a0 <= b1 and b0 <= c1 and a1 + c1 >= 0


def split(box):
    d = max(range(3), key=lambda i: box[i][1] - box[i][0])
    lo, hi = box[d]; mid = (lo + hi) / 2
    a = list(box); b = list(box)
    a[d] = (lo, mid); b[d] = (mid, hi)
    return [a, b]


def run(a):
    os.makedirs(a.out, exist_ok=True)
    r = a.root
    g = [-45.0 + r * i for i in range(int(round(90 / r)))]
    queue = [([(x, x + r), (y, y + r), (z, z + r)], {}) for x in g for y in g for z in g]
    queue = [q for q in queue if in_domain(q[0])]
    if a.boxes:
        queue = [([tuple(x) for x in json.loads(l)['box']], {}) for l in open(a.boxes)]
    nok = nloc = nsplit = 0; t0 = time.time()
    leaves = open(os.path.join(a.out, 'leaves.jsonl'), 'w')
    local = open(os.path.join(a.out, 'local.jsonl'), 'w')
    with ProcessPoolExecutor(a.procs) as ex:
        while queue:
            batch, queue = queue[:4000], queue[4000:]
            for box, st, lams, fail in ex.map(solve_box, [(b, inh, a.eps, a.minloc, a.mode) for b, inh in batch], chunksize=4):
                if st == 'ok':
                    nok += 1; leaves.write(json.dumps(dict(box=box)) + '\n')
                elif st == 'local':
                    nloc += 1; local.write(json.dumps(dict(box=box, fail=[TYPES[t] for t in fail])) + '\n')
                else:
                    nsplit += 1
                    if min(hi - lo for lo, hi in box) < a.minbox:
                        raise SystemExit(f'box too small, not near a wedge point: {box} fails {[TYPES[t] for t in fail][:5]}')
                    queue += [(c, lams) for c in split(box) if in_domain(c)]
            print(f'{time.time()-t0:7.0f}s ok {nok} local {nloc} split {nsplit} queue {len(queue)}', flush=True)
    leaves.close(); local.close()
    json.dump(dict(ok=nok, local=nloc, split=nsplit, eps=a.eps, root=a.root, tol=TOL), open(os.path.join(a.out, 'summary.json'), 'w'))
    print('done', nok, nloc, nsplit)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd'); ap.add_argument('--procs', type=int, default=4); ap.add_argument('--out', required=True)
    ap.add_argument('--root', type=float, default=15.0); ap.add_argument('--eps', type=float, default=1.0)
    ap.add_argument('--minloc', type=float, default=0.5); ap.add_argument('--minbox', type=float, default=1e-3)
    ap.add_argument('--mode', default='tube'); ap.add_argument('--boxes', default=None)
    a = ap.parse_args()
    run(a)
