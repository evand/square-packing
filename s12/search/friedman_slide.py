"""Cut-and-slide descent on record packings (Friedman Conj. 1 experiment).

A packing of N unit squares in [0,s]^2, m = floor(s) < s.  Put r = s - m + EPS.  Each square is either removed or
translated by one of (0,0), (-r,0), (0,-r), (-r,-r); survivors must be pairwise interior-disjoint and lie in
[0,s-r]^2 (side < m).  Max survivors M by ILP.  Friedman's step needs M >= N - (2m+1)  (slack = M - N + 2m + 1).
Float geometry (SAT with tolerance), not a proof of anything.
"""
import json, math, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

TOL = 1e-9
EPS = 1e-5
STATUS = [None]
SHIFTS = [(0, 0), (1, 0), (0, 1), (1, 1)]  # multiples of -r


def corners(cx, cy, th):
    c, s = math.cos(th), math.sin(th)
    u, v = np.array([c, s]) / 2, np.array([-s, c]) / 2
    p = np.array([cx, cy])
    return np.array([p + u + v, p - u + v, p - u - v, p + u - v]), [np.array([c, s]), np.array([-s, c])]


def overlap(A, B):
    """interiors overlap by more than TOL along every SAT axis"""
    (pa, aa), (pb, ab) = A, B
    for ax in aa + ab:
        x, y = pa @ ax, pb @ ax
        if min(x.max(), y.max()) - max(x.min(), y.min()) <= TOL:
            return False
    return True


def load(name, d):
    v = d[name]
    s = float(v['s'])
    sq = [corners(float(a), float(b), math.radians(float(t))) for a, b, t in v['squares']]
    return s, sq


def check_packing(s, sq):
    bad = 0
    for i in range(len(sq)):
        p = sq[i][0]
        if p.min() < -1e-9 or p.max() > s + 1e-9:
            bad += 1
        for j in range(i):
            if overlap(sq[i], sq[j]):
                bad += 1
    return bad


def descend(s, sq, m=None):
    if m is None:
        m = math.floor(s)
    r = s - m + EPS
    side = s - r
    N = len(sq)
    moved = {}
    var = {}
    for i, (p, ax) in enumerate(sq):
        for q, (a, b) in enumerate(SHIFTS):
            pp = p - np.array([a * r, b * r])
            if pp.min() >= -TOL and pp.max() <= side + TOL:
                var[(i, q)] = len(var)
                moved[(i, q)] = (pp, ax)
    nv = len(var)
    rows = []
    for i in range(N):
        idx = [var[(i, q)] for q in range(len(SHIFTS)) if (i, q) in var]
        if len(idx) > 1:
            rows.append(idx)
    keys = list(var)
    cent = {k: moved[k][0].mean(0) for k in keys}
    for a in range(len(keys)):
        for b in range(a):
            ka, kb = keys[a], keys[b]
            if ka[0] == kb[0]:
                continue
            if np.linalg.norm(cent[ka] - cent[kb]) > math.sqrt(2) + 1e-6:
                continue
            if overlap(moved[ka], moved[kb]):
                rows.append([var[ka], var[kb]])
    if nv == 0:
        return 0, m, r, {}
    A = np.zeros((len(rows), nv))
    for k, row in enumerate(rows):
        A[k, row] = 1
    res = milp(c=-np.ones(nv), constraints=[LinearConstraint(A, -np.inf, 1)] if rows else [],
               integrality=np.ones(nv), bounds=Bounds(0, 1), options={'time_limit': 600})
    STATUS[0] = res.status
    x = np.round(res.x).astype(int)
    M = int(x.sum())
    assign = {k[0]: k[1] for k in keys if x[var[k]]}
    return M, m, r, assign


if __name__ == '__main__':
    d = json.load(open(sys.argv[1]))
    names = sys.argv[2:] or None
    best = {}
    for k, v in d.items():
        if v['errors']:
            continue
        n, s = v['n'], float(v['s'])
        if n not in best or s < best[n][0]:
            best[n] = (s, k)
    allp = names == ['--all']
    items = sorted(((v['n'], k) for k, v in d.items() if not v['errors']), key=lambda x: x[0]) if allp \
        else [(n, best[n][1]) for n in sorted(best)]
    for n, name in items:
        s = float(d[name]['s'])
        if names and not allp and name not in names:
            continue
        if n > 130 or abs(s - round(s)) < 1e-9 or s < 2:
            continue
        s, sq = load(name, d)
        bad = check_packing(s, sq)
        M, m, r, assign = descend(s, sq)
        c = (m + 1) ** 2 - n
        need = n - (2 * m + 1)
        print(f"{name:24s} n={n:4d} s={s:.6f} m={m} c={c:3d} r={r:.4f} bad={bad} survivors={M:4d} need={need:4d} "
              f"slack={M-need:+d}  => s({M}) < {m}", flush=True)
