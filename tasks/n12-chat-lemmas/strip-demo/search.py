"""Float-guided search for a wall-strip certificate (proves nothing; check.py does the proving).

Usage:  python3 search.py H [--a 457/500] [--step 1/100] [--max 400] [--out certs/h113.json]

Produces a certificate:  rational boxes covering the pose space, a menu of separation heights,
a rule Z[i][j] (the separation height used when box i is immediately left of box j), and the
table of node-value lower bounds the DP uses.  Box endpoints are exact rationals from the start
(bisection of rationals); only the bounds used to steer the search are floats.
"""
import argparse, json, math, time
from fractions import Fraction as Fr
import numpy as np
from leaf import Box, Dof

INF = 1e9
MARGIN = 1e-4   # float search target: DP value >= 4 + MARGIN (the exact check needs only >= 4)
NEGV = -1e9


def parse_fr(s):
    return Fr(s)


class Search:
    def __init__(self, h, a, step, menu=None):
        self.h, self.a = h, a
        lo, hi = h - Fr(1, 2), Fr(1)
        M = set([a])
        if menu:
            M |= set(menu)
        else:
            z = Fr(64, 100)
            while z < hi:
                if lo < z:
                    M.add(z)
                z += step
        self.M = sorted(m for m in M if lo < m < hi)
        self.nz = len(self.M)
        self.ia = self.M.index(a)
        self.Zs = self.M + [None]
        self.cache = {}

    # boxes are tuples of Fractions (t0, t1, c0, c1)
    def table(self, key):
        if key not in self.cache:
            fb = Box(*[float(x) for x in key], exact=False)
            eb = Box(*key, exact=True)
            T = np.full((self.nz + 1, self.nz + 1), NEGV)
            Mf = [float(z) for z in self.M] + [None]
            for i, z1 in enumerate(Mf):
                for j, z2 in enumerate(Mf):
                    if z1 is None and z2 is None:
                        continue
                    v = fb.V(z1, z2)
                    T[i, j] = NEGV if v is None else v
            self.cache[key] = (T, bool(eb.may_be_H(self.a)), eb)
        return self.cache[key]

    def fwd_bwd(self, T, H, Z):
        m, nz1, _ = T.shape
        nz = nz1 - 1
        TT = T[np.arange(m)[:, None, None], np.arange(nz1)[None, :, None], Z[:, None, :]]
        F = [None] * 5
        f = np.full((m, nz1, 2), INF)
        f[np.arange(m), nz, H] = 0.0
        F[1] = f
        for k in range(1, 4):
            g = np.full((m, nz1, 2), INF)
            for fl in range(2):
                val = np.where(f[:, :, fl][:, :, None] >= INF, INF, f[:, :, fl][:, :, None] + TT).min(1)
                for b2 in range(m):
                    np.minimum.at(g[b2, :, fl | H[b2]], Z[:, b2], val[:, b2])
            f = g
            F[k + 1] = f
        G = [None] * 5
        g = np.full((m, nz1, 2), INF)
        g[np.arange(m), :, H] = T[:, :, nz]
        G[4] = g
        for k in (3, 2, 1):
            gn = np.full((m, nz1, 2), INF)
            for fl in range(2):
                nxt = G[k + 1][np.arange(m)[None, :], Z, fl]
                val = np.where(nxt[:, None, :] >= INF, INF, TT + nxt[:, None, :]).min(2)
                for b in range(m):
                    gn[b, :, fl | H[b]] = np.minimum(gn[b, :, fl | H[b]], val[b])
            G[k] = gn
        tot = np.where(F[4][:, :, 1] >= INF, INF, F[4][:, :, 1] + T[:, :, nz])
        return tot.min(), F, G, tot

    def improve_rule(self, T, F, G):
        m, nz1, _ = T.shape
        nz = nz1 - 1
        score = np.full((m, m, nz), INF)
        for k in (1, 2, 3):
            f = F[k]
            A = np.where(f[:, :, None, :] >= INF, INF, f[:, :, None, :] + T[:, :, :nz, None]).min(1)
            g = G[k + 1][:, :nz, :]
            best = np.full((m, m, nz), INF)
            for f1 in range(2):
                for f2 in range(2):
                    if f1 | f2:
                        best = np.minimum(best, A[:, None, :, f1] + g[None, :, :, f2])
            score = np.minimum(score, best)
        return score.argmax(2)

    def solve(self, boxes, iters=8):
        tabs = [self.table(b) for b in boxes]
        T = np.array([t[0] for t in tabs])
        H = np.array([t[1] for t in tabs]).astype(int)
        m = len(boxes)
        Z = np.full((m, m), self.ia, int)
        bestv, bestZ = -INF, Z
        for _ in range(iters):
            v, F, G, tot = self.fwd_bwd(T, H, Z)
            if v > bestv + 1e-12:
                bestv, bestZ = v, Z.copy()
            Zn = self.improve_rule(T, F, G)
            if (Zn == Z).all():
                break
            Z = Zn
        v, F, G, tot = self.fwd_bwd(T, H, bestZ)
        return v, bestZ, F, G, T, H, tot

    def worst_chain(self, F, T, H, Z, tot):
        m = T.shape[0]
        b4, zl4 = np.unravel_index(tot.argmin(), tot.shape)
        chain = [(b4, zl4, 1)]
        for k in (3, 2, 1):
            b2, z2, fl2 = chain[-1]
            best = (3 * INF, None)
            for b in range(m):
                if Z[b, b2] != z2:
                    continue
                for fl in range(2):
                    if (fl | H[b2]) != fl2:
                        continue
                    v = np.where(F[k][b, :, fl] >= INF, 2 * INF, F[k][b, :, fl] + T[b, :, z2])
                    i = int(v.argmin())
                    if v[i] < best[0]:
                        best = (v[i], (b, i, fl))
            chain.append(best[1])
        return [c[0] for c in chain[::-1]]


def split(key, d=None):
    t0, t1, c0, c1 = key
    eb = Box(*key)
    cl = max(c0, eb.plb)
    # snap the lower c end to a short rational (cl is only a bound; boxes may contain c < p poses)
    dt = float(t1 - t0) * 1.6
    dc = float(c1 - cl)
    if d == 0 or (d is None and dt >= dc):
        tm = (t0 + t1) / 2
        return [(t0, tm, c0, c1), (tm, t1, c0, c1)]
    cm = Fr(math.floor(float((cl + c1) / 2) * 10 ** 6), 10 ** 6)
    if not (c0 < cm < c1):
        cm = (c0 + c1) / 2
    return [(t0, t1, c0, cm), (t0, t1, cm, c1)]


def nonempty(key):
    return Box(*key).nonempty()


def seed(h, a, nH):
    """Thin H strips (c from a + D_lb to h) plus the rest, with short rational t-cuts."""
    th = np.linspace(0, np.pi / 2, 400001)
    u = np.cos(th) + np.sin(th)
    D = (u - u * u + 1) / 2
    Hm = float(a) + D < float(h)
    if not Hm.any():
        return [(Fr(0), Fr(1), Fr(1, 2), h)]
    lo, hi = th[Hm].min() - 1e-4, th[Hm].max() + 1e-4
    cuts = [Fr(math.floor(math.tan(x / 2) * 10 ** 4), 10 ** 4) for x in np.linspace(lo, hi, nH + 1)]
    cuts[-1] = Fr(math.ceil(math.tan(hi / 2) * 10 ** 4), 10 ** 4)
    boxes = []
    for i in range(nH):
        eb = Box(cuts[i], cuts[i + 1], Fr(1, 2), h)
        cH = Fr(math.floor(float(a + eb.Dlb) * 10 ** 4), 10 ** 4)
        boxes += [(cuts[i], cuts[i + 1], cH, h), (cuts[i], cuts[i + 1], Fr(1, 2), cH)]
    boxes += [(Fr(0), cuts[0], Fr(1, 2), h), (cuts[-1], Fr(1), Fr(1, 2), h)]
    return [b for b in boxes if nonempty(b)]


def try_merge(S, boxes):
    """Undo splits where possible: merge two boxes forming a rectangle if the DP still reaches 4."""
    changed = True
    while changed:
        changed = False
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                A, B = boxes[i], boxes[j]
                if A[0] == B[0] and A[1] == B[1] and (A[3] == B[2] or B[3] == A[2]):
                    nb = (A[0], A[1], min(A[2], B[2]), max(A[3], B[3]))
                elif A[2] == B[2] and A[3] == B[3] and (A[1] == B[0] or B[1] == A[0]):
                    nb = (min(A[0], B[0]), max(A[1], B[1]), A[2], A[3])
                else:
                    continue
                trial = [b for k, b in enumerate(boxes) if k not in (i, j)] + [nb]
                if S.solve(trial)[0] >= 4 + MARGIN:
                    boxes = trial
                    changed = True
                    break
            if changed:
                break
    return boxes


def simplify_rule(S, boxes, Z, margin=None):
    margin = MARGIN if margin is None else margin
    """Reset rule entries to the default height a when the DP value stays >= 4 (fewer leaves)."""
    tabs = [S.table(b) for b in boxes]
    T = np.array([t[0] for t in tabs]); H = np.array([t[1] for t in tabs]).astype(int)
    m = len(boxes)
    for i in range(m):
        for j in range(m):
            if Z[i, j] == S.ia:
                continue
            old = Z[i, j]; Z[i, j] = S.ia
            if S.fwd_bwd(T, H, Z)[0] < 4 + margin:
                Z[i, j] = old
    # second pass: make each box's right/left heights as few as possible
    for i in range(m):
        for j in range(m):
            old = Z[i, j]
            if old == S.ia:
                continue
            for cand in sorted(set(Z[i, :]) | set(Z[:, j])):
                if cand == old:
                    continue
                Z[i, j] = cand
                if S.fwd_bwd(T, H, Z)[0] >= 4 + margin:
                    old = cand
                    break
                Z[i, j] = old
    return Z


def run(h, a, step, maxb, nH=1, merge=True, verbose=False, menu=None, simplify=True):
    S = Search(h, a, step, menu)
    boxes = seed(h, a, nH)
    t = time.time()
    while True:
        v, Z, F, G, T, H, tot = S.solve(boxes)
        if verbose:
            print(len(boxes), round(v, 5), flush=True)
        if v >= 4 + MARGIN or len(boxes) >= maxb:
            break
        ch = S.worst_chain(F, T, H, Z, tot)
        i = max(set(ch), key=lambda i: float(boxes[i][1] - boxes[i][0]) * 1.6
                + float(boxes[i][3] - max(boxes[i][2], Box(*boxes[i]).plb)))
        nb = [x for x in split(boxes[i]) if nonempty(x)]
        boxes = boxes[:i] + boxes[i + 1:] + nb
    ok = v >= 4 + MARGIN
    nraw = len(boxes)
    if ok and merge:
        boxes = try_merge(S, boxes)
        v, Z, F, G, T, H, tot = S.solve(boxes)
    if ok and simplify:
        Z = simplify_rule(S, boxes, Z)
        v = S.fwd_bwd(np.array([S.table(b)[0] for b in boxes]),
                      np.array([S.table(b)[1] for b in boxes]).astype(int), Z)[0]
    return S, boxes, Z, v, ok, nraw, time.time() - t


def needed_pairs(Z, i, m):
    """(zl, zr) index pairs box i can see in a chain (index -1 = END)."""
    L = sorted(set(int(Z[j, i]) for j in range(m))) + [-1]
    R = sorted(set(int(Z[i, j]) for j in range(m))) + [-1]
    return [(l, r) for l in L for r in R if not (l == -1 and r == -1)]


def export(S, boxes, Z, path, quantum=10 ** 6):
    """Write the certificate.  Leaf claims are the exact bounds rounded DOWN to multiples of 1/quantum."""
    fs = lambda x: str(x)
    m = len(boxes)
    leaves = []
    for i, key in enumerate(boxes):
        eb = Box(*key, exact=True)
        claims = []
        for l, r in needed_pairs(Z, i, m):
            z1 = None if l == -1 else S.M[l]
            z2 = None if r == -1 else S.M[r]
            v = eb.V(z1, z2)
            if v is None:
                continue          # no bound: the DP treats the pair as -infinity
            q = Fr(math.floor(v * quantum), quantum)
            claims.append([l, r, fs(q)])
        leaves.append(claims)
    cert = dict(h=fs(S.h), a=fs(S.a), menu=[fs(z) for z in S.M],
                boxes=[[fs(x) for x in b] for b in boxes],
                rule=[[int(Z[i, j]) for j in range(m)] for i in range(m)],
                leaves=leaves)
    json.dump(cert, open(path, 'w'), indent=0)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('h')
    ap.add_argument('--a', default='457/500')
    ap.add_argument('--step', default='1/100')
    ap.add_argument('--max', type=int, default=400)
    ap.add_argument('--nH', type=int, default=1)
    ap.add_argument('--out', default=None)
    ap.add_argument('-v', action='store_true')
    ap.add_argument('--nomerge', action='store_true')
    ap.add_argument('--nosimplify', action='store_true', help='skip the O(m^2) rule simplification')
    ap.add_argument('--menu', default=None, help='comma-separated rationals (a is always added)')
    args = ap.parse_args()
    h, a, step = Fr(args.h), Fr(args.a), Fr(args.step)
    menu = [Fr(x) for x in args.menu.split(',')] if args.menu else None
    S, boxes, Z, v, ok, nraw, dt = run(h, a, step, args.max, args.nH, not args.nomerge, args.v, menu,
                                       not args.nosimplify)
    print(f"h={args.h} ok={ok} boxes={len(boxes)} (before merge {nraw}) float-value={v:.6f} time={dt:.1f}s")
    if ok and args.out:
        export(S, boxes, Z, args.out)
        print("wrote", args.out)
