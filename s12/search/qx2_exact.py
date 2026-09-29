#!/usr/bin/env python3
"""qx2_exact.py -- exact (Fraction) tools for the qx2 family (task quadrant-exact-w2).

  resolve SOL            exact rational re-solve of a qx2_lp vertex: rebuild the LP on SOL's rows, take HiGHS's
                         simplex basis, recompute the tight rows EXACTLY (rational poses, exact chord fractions,
                         exact Lebesgue areas; theta = 0 grid rows as exact one-sided limits) and solve the basis
                         system in Fractions, sigma = 0 kept as an equality.  Writes SOL_exact.json.
  axis EXACT.json        exact check of the theta = 0 face (Lemma Z, QUADRANT_EXACT.md): every one-sided limit
                         corner of the breakpoint grid of the quadrant window, exactly.
  cover EXACT.json K OUT write the D4 box measure mu_K as a FORMAT.md v1 mixed cover (segments + the Lebesgue square
                         [a, K-a]^2 as one polygon), exact.

Everything that certifies is Fraction arithmetic.  The float LP only chooses the basis.
"""
import sys, os, math, json, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np

HALF = F(1, 2)


def fr(v, den=1000):
    """float model coordinate -> exact rational (model coordinates are multiples of 1/20)."""
    q = F(v).limit_denominator(den)
    assert abs(float(q) - v) < 1e-9, (v, q)
    return q


def trig(u):
    u = F(u); N = 1 + u * u
    return (1 - u * u) / N, 2 * u / N


# ------------------------------------------------------------------------------------------------ exact elements
class ExactModel:
    """exact version of a qx2 quadrant model: lines (orient, coordinate, t0, t1, var) and U = [a, inf)^2."""

    def __init__(self, m):
        self.R = m.R; self.w = fr(m.w); self.a = fr(m.a); self.nvar = m.nvar
        self.kappa = m.kappa; self.th0 = m.th0
        self.lines = []           # ('h', y, x0, x1, var)  /  ('v', x, y0, y1, var)
        self.fixed_mass = []
        for k, arr in m.A.items():
            if len(arr) == 0: continue
            if k in ('pts', 'cells'):
                raise ValueError(f"exact tools: {k} not supported (lines-only models)")
            for e in arr:
                var = int(e[-2])
                if var < 0: raise ValueError("fixed atoms not supported")
                self.lines.append((k[0], fr(e[0]), fr(e[1]), fr(e[2]), var))
        self.const_obj = fr(m.const_obj); self.const_sig = fr(m.const_sig)
        self.cobj = [fr(c) for c in m.cobj]; self.csig = [fr(c) for c in m.csig]

    # -- exact capture ----------------------------------------------------------------------------------------
    def coeffs(self, cx, cy, u, sx=0, sy=0):
        """exact (coefficient dict var -> fraction captured, Lebesgue area) at the closed square (cx, cy, u).
        u = 0 with sx/sy = +-1: one-sided limit as cx -> cx +- 0 (resp. cy)."""
        cx = F(cx); cy = F(cy); u = F(u)
        co = {}
        if u == 0:
            for (o, p, t0, t1, var) in self.lines:
                if o == 'h':     # y = p, x in [t0,t1]
                    if not _in1(p, cy, sy): continue
                    L = _ov(t0, t1, cx - HALF, cx + HALF)
                else:
                    if not _in1(p, cx, sx): continue
                    L = _ov(t0, t1, cy - HALF, cy + HALF)
                if L > 0: co[var] = co.get(var, 0) + L / (t1 - t0)
            ar = _ov(self.a, F(10 ** 6), cx - HALF, cx + HALF) * _ov(self.a, F(10 ** 6), cy - HALF, cy + HALF)
            return co, ar
        assert sx == 0 and sy == 0
        C, S = trig(u)
        for (o, p, t0, t1, var) in self.lines:
            if o == 'h': P0 = (F(0), p); d = (F(1), F(0))
            else: P0 = (p, F(0)); d = (F(0), F(1))
            iv = chord_iv(P0, d, cx, cy, C, S)
            if iv is None: continue
            lo, hi = iv
            L = min(hi, t1) - max(lo, t0)
            if L > 0: co[var] = co.get(var, 0) + L / (t1 - t0)
        Qv = square_vertices(cx, cy, C, S)
        ar = poly_area(clip_hp(clip_hp(Qv, -1, 0, -self.a), 0, -1, -self.a))     # x >= a, y >= a
        return co, ar

    def margin(self, u, ar):
        if self.kappa <= 0 or u == 0: return F(0)
        th = 2 * math.atan(float(u)); t = min(th % (math.pi / 2), math.pi / 2 - th % (math.pi / 2))
        return F(self.kappa * min(t, self.th0)).limit_denominator(10 ** 12) * max(F(0), 1 - ar)


def _in1(p, c, s):
    """p in [c-1/2, c+1/2] with one-sided limit semantics s (c -> c + s*0)."""
    lo, hi = c - HALF, c + HALF
    if s > 0: return lo < p <= hi
    if s < 0: return lo <= p < hi
    return lo <= p <= hi


def _ov(a0, a1, b0, b1):
    return max(F(0), min(a1, b1) - max(a0, b0))


def chord_iv(P0, d, cx, cy, C, S):
    """parameter interval {t : P0 + t d in closed square (c, theta)}, or None."""
    lo, hi = None, None
    ax, ay = P0[0] - cx, P0[1] - cy
    for a0, a1 in ((ax * C + ay * S, d[0] * C + d[1] * S), (-ax * S + ay * C, -d[0] * S + d[1] * C)):
        for sg in (1, -1):              # sg*(a0 + t a1) <= 1/2
            b0, b1 = sg * a0 - HALF, sg * a1
            if b1 > 0:
                v = -b0 / b1; hi = v if hi is None else min(hi, v)
            elif b1 < 0:
                v = -b0 / b1; lo = v if lo is None else max(lo, v)
            elif b0 > 0:
                return None
    if lo is None or hi is None or lo > hi: return None
    return lo, hi


def square_vertices(cx, cy, C, S):
    return [(cx + (sx * C - sy * S) / 2, cy + (sx * S + sy * C) / 2) for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]


def clip_hp(V, a, b, c):
    """V ∩ {a x + b y <= c}"""
    out = []
    n = len(V)
    for i in range(n):
        P, Qp = V[i], V[(i + 1) % n]
        sp = a * P[0] + b * P[1] - c; sq = a * Qp[0] + b * Qp[1] - c
        if sp <= 0: out.append(P)
        if (sp < 0 < sq) or (sq < 0 < sp):
            t = sp / (sp - sq)
            out.append((P[0] + t * (Qp[0] - P[0]), P[1] + t * (Qp[1] - P[1])))
    return out


def poly_area(V):
    if len(V) < 3: return F(0)
    s = F(0)
    for i in range(len(V)):
        x0, y0 = V[i]; x1, y1 = V[(i + 1) % len(V)]
        s += x0 * y1 - x1 * y0
    return abs(s) / 2


# ------------------------------------------------------------------------------------------------ exact re-solve
def snap_row(P, gridden=10, tol=1e-6):
    """float LP row pose -> (cx, cy, u, sx, sy) exact.  theta = 0 rows within tol of the 1/gridden grid become
    one-sided limits at the grid point (the LP used grid +- 1e-7 as a stand-in for the limit)."""
    cx, cy, th = (float(t) for t in P)
    if th == 0.0:
        out = []
        for c in (cx, cy):
            g = round(c * gridden) / gridden
            if abs(c - g) < tol and c != g:
                out += [F(round(c * gridden), gridden), 1 if c > g else -1]
            elif c == g:
                out += [F(round(c * gridden), gridden), 0]
            else:
                out += [F(c).limit_denominator(10 ** 9), 0]
        return out[0], out[2], F(0), out[1], out[3]
    u = F(math.tan(th / 2)).limit_denominator(10 ** 9)
    C, S = trig(u); hw = (abs(C) + abs(S)) / 2
    qx = max(F(cx).limit_denominator(10 ** 9), hw); qy = max(F(cy).limit_denominator(10 ** 9), hw)
    return qx, qy, u, 0, 0


def resolve(path, out=None):
    import highspy
    import qx2_tools as T
    import quadrant_lp as Q
    m, x, args, d = T.load(path)
    rows = d['rows']
    lp = Q.LP(m); lp.solver = 'simplex'
    t0 = time.time()
    lp.add(rows)
    xf, obj, st = lp.solve()
    print(f"re-solved float LP on {lp.nrows} rows: {st}, D = {m.R**2 - obj - m.const_obj:.9f} [{time.time()-t0:.1f}s]")
    h = lp.h
    basis = h.getBasis()
    cs = [s for s in basis.col_status]; rs = [s for s in basis.row_status]
    B = highspy.HighsBasisStatus.kBasic
    bcols = [j for j, s in enumerate(cs) if s == B]
    nbrows = [i for i, s in enumerate(rs) if s != B]
    print(f"basis: {len(bcols)} basic columns, {len(nbrows)} nonbasic rows (row 0 = sigma equality)")
    assert len(nbrows) == len(bcols), "basis not square?"
    EM = ExactModel(m)
    # row 0 of the HiGHS model is the sigma equality (added in LP.__init__); rows 1.. are lp.P in order
    sysA = []; sysb = []
    for i in nbrows:
        if i == 0:
            sysA.append({j: EM.csig[j] for j in bcols if EM.csig[j] != 0}); sysb.append(EM.w - EM.const_sig)
            continue
        P = lp.P[i - 1]
        cx, cy, u, sx, sy = snap_row(P)
        co, ar = EM.coeffs(cx, cy, u, sx, sy)
        rhs = 1 + EM.margin(u, ar) - ar
        sysA.append({j: co.get(j, F(0)) for j in bcols}); sysb.append(rhs)
    xs = gauss(sysA, sysb, bcols)
    xe = [F(0)] * m.nvar
    for j, v in zip(bcols, xs): xe[j] = v
    neg = [j for j in bcols if xe[j] < 0]
    D = m.R ** 2 - sum(EM.cobj[j] * xe[j] for j in range(m.nvar)) - EM.const_obj
    sig = EM.w - EM.const_sig - sum(EM.csig[j] * xe[j] for j in range(m.nvar))
    print(f"exact vertex: D = {D} = {float(D):.12f}; sigma = {sig}; negative entries: {len(neg)}; "
          f"max |x_exact - x_float| = {max(abs(float(xe[j]) - xf[j]) for j in range(m.nvar)):.3e}")
    res = dict(sol=path, args=vars(args), x=[str(v) for v in xe], D=str(D), Dfloat=float(D), sigma=str(sig),
               bcols=bcols, nbrows=nbrows, negative=neg)
    out = out or path.replace('.npz', '_exact.json')
    json.dump(res, open(out, 'w'))
    print("wrote", out)
    return res


def axis_corners(EM, R):
    """all (cx, cy, sx, sy) one-sided limit corners of the theta = 0 breakpoint grid in the quadrant window
    cx in [1/2, R + 11/4], cy in [1/2, min(cx, R + 7/4)] (a superset of the window of QUADRANT.md sec 2)."""
    ex = R + F(7, 4); ey = R + F(3, 4)
    g = breakgrid(EM, HALF, ex + 1)
    for cx in g:
        for cy in g:
            if cy > cx or cy > ey + 1: continue
            for sx in (-1, 1):
                for sy in (-1, 1):
                    if cx == HALF and sx < 0: continue
                    if cy == HALF and sy < 0: continue
                    yield cx, cy, sx, sy


def project(path, tol=1e-7, den=10 ** 12, out=None):
    """exact solution near the float one: every theta = 0 limit corner whose float value is < 1 + tol becomes an
    exact equality (value = 1), plus sigma = 0; support fixed; free variables rounded to 1/den, the pivots solved
    exactly.  (Tilted rows are not used: they carry the LP margin, and validity there is the checker's job.)"""
    import qx2_tools as T
    m, x, args, d = T.load(path)
    EM = ExactModel(m)
    sup = [j for j in range(m.nvar) if x[j] > 1e-13]
    eqs = []
    for cx, cy, sx, sy in axis_corners(EM, m.R):
        co, ar = EM.coeffs(cx, cy, 0, sx, sy)
        v = float(ar) + sum(float(c) * x[j] for j, c in co.items())
        if v < 1 + tol:
            if any(j not in sup and c != 0 for j, c in co.items()) and False: pass
            eqs.append(({j: c for j, c in co.items() if j in sup}, 1 - ar, v))
    eqs.append(({j: EM.csig[j] for j in sup if EM.csig[j] != 0}, EM.w - EM.const_sig, None))
    print(f"{len(sup)} support variables, {len(eqs)} equalities (tight theta = 0 limit corners + sigma)")
    # exact row reduction
    col = {j: k for k, j in enumerate(sup)}; n = len(sup)
    M = []
    for co, b, _ in eqs:
        r = [F(0)] * (n + 1)
        for j, c in co.items(): r[col[j]] = F(c)
        r[n] = F(b); M.append(r)
    piv = []; rk = 0
    for k in range(n):
        p = next((i for i in range(rk, len(M)) if M[i][k] != 0), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        pk = M[rk][k]; M[rk] = [v / pk for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][k] != 0:
                f = M[i][k]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        piv.append(k); rk += 1
    incons = [i for i in range(rk, len(M)) if M[i][n] != 0]
    print(f"rank {rk}; inconsistent rows after reduction: {len(incons)}")
    if incons:
        print("  max residual", max(abs(float(M[i][n])) for i in incons))
    free = [k for k in range(n) if k not in piv]
    xs = [None] * n
    for k in free: xs[k] = F(int(round(float(x[sup[k]]) * den)), den)
    for r_, k in enumerate(piv):
        v = M[r_][n] - sum(M[r_][f] * xs[f] for f in free)
        xs[k] = v
    xe = [F(0)] * m.nvar
    for k, j in enumerate(sup): xe[j] = xs[k]
    neg = [j for j in sup if xe[j] < 0]
    dev = max(abs(float(xe[j]) - x[j]) for j in range(m.nvar))
    D = m.R ** 2 - sum(EM.cobj[j] * xe[j] for j in range(m.nvar)) - EM.const_obj
    sig = EM.w - EM.const_sig - sum(EM.csig[j] * xe[j] for j in range(m.nvar))
    print(f"exact: D = {float(D):.12f} (float LP {float(d['D']):.12f}); sigma = {sig}; negative: {len(neg)}; "
          f"max |x_exact - x_float| = {dev:.3e}; denominators up to {max(v.denominator for v in xe)}")
    res = dict(sol=path, args=vars(args), x=[str(v) for v in xe], D=str(D), Dfloat=float(D), sigma=str(sig),
               n_eq=len(eqs), rank=rk, inconsistent=len(incons), negative=neg, dev=dev)
    out = out or path.replace('.npz', '_exact.json')
    json.dump(res, open(out, 'w')); print("wrote", out)
    return res


def gauss(Arows, b, cols):
    n = len(cols); idx = {c: k for k, c in enumerate(cols)}
    M = [[F(0)] * n + [F(bi)] for bi in b]
    for r, row in enumerate(Arows):
        for c, v in row.items(): M[r][idx[c]] = F(v)
    for k in range(n):
        piv = max(range(k, n), key=lambda r: (M[r][k] != 0, abs(M[r][k])))
        if M[piv][k] == 0: raise ValueError("singular basis system")
        M[k], M[piv] = M[piv], M[k]
        pk = M[k][k]
        for r in range(n):
            if r != k and M[r][k] != 0:
                f = M[r][k] / pk
                M[r] = [a - f * bb for a, bb in zip(M[r], M[k])]
    return [M[k][n] / M[k][k] for k in range(n)]


def load_exact(path):
    import qx2_lp as X
    d = json.load(open(path))
    args = argparse.Namespace(**d['args'])
    m = X.build(args)
    return m, ExactModel(m), [F(v) for v in d['x']], d


# ------------------------------------------------------------------------------------------------ theta = 0 face
def breakgrid(EM, lo, hi):
    xs = set([EM.a])
    for (o, p, t0, t1, var) in EM.lines: xs.update([p, t0, t1])
    out = set()
    for v in xs:
        for s in (-HALF, HALF):
            if lo <= v + s <= hi: out.add(v + s)
    return sorted(out)


def axis_check(path, verbose=True):
    """Lemma Z: at theta = 0 the mass is multilinear on each open cell of the grid (breakpoints +- 1/2), so its
    infimum over admissible poses of the quadrant window is the minimum over the one-sided limits at the cell
    corners (and the closed values are >= those limits).  Checked exactly for every corner and sign pattern."""
    m, EM, xe, d = load_exact(path)
    R = m.R; ex = R + F(7, 4); ey = R + F(3, 4)
    gx = [g for g in breakgrid(EM, HALF, ex + 1)]; gy = [g for g in breakgrid(EM, HALF, ex + 1)]
    worst = None; n = 0; tight = 0
    for cx in gx:
        for cy in gy:
            if cy > cx or cy > ey + 1 or cx > ex + 1: continue
            for sx in (-1, 1):
                for sy in (-1, 1):
                    if cx == HALF and sx < 0: continue
                    if cy == HALF and sy < 0: continue
                    co, ar = EM.coeffs(cx, cy, 0, sx, sy)
                    v = ar + sum(c * xe[j] for j, c in co.items())
                    n += 1
                    if v == 1: tight += 1
                    if worst is None or v < worst[0]: worst = (v, (cx, cy, sx, sy))
    if verbose:
        print(f"axis face: {n} one-sided limit corners, min = {worst[0]} = {float(worst[0]):.15f} at "
              f"{tuple(str(t) for t in worst[1])}; exactly tight: {tight}")
    return worst, n, tight


# ------------------------------------------------------------------------------------------------ the D4 box cover
def box_elements(m, EM, xe, k):
    """exact D4 box measure mu_k: list of segments (x0, y0, x1, y1, mass) and the Lebesgue square [a, k-a]^2.
    Built element by element: the corner module (KA == 0 elements, already diagonal-symmetric) at the 4 corners;
    the profile base atoms (m.base) at x = j + phase, j = R .. k-R, pieces inside [R, k-R], rotated to the 4 walls."""
    k = F(k); R = EM.R
    segs = []
    maps = [lambda a, b: (a, b), lambda a, b: (k - a, b), lambda a, b: (a, k - b), lambda a, b: (k - a, k - b)]
    for kk, arr in m.A.items():
        for e, t in zip(arr, m.KA[kk]):
            if t != 0: continue
            var = int(e[-2]); ms = xe[var]
            if ms == 0: continue
            p, t0, t1 = fr(e[0]), fr(e[1]), fr(e[2])
            for f in maps:
                if kk == 'hs': A, B = f(t0, p), f(t1, p)
                else: A, B = f(p, t0), f(p, t1)
                segs.append((A[0], A[1], B[0], B[1], ms))
    rots = [lambda a, b: (a, b), lambda a, b: (k - b, a), lambda a, b: (k - a, k - b), lambda a, b: (b, k - a)]
    for (typ, a, b, y0, y1, var, fm) in m.base:
        ms = xe[int(var)]
        if ms == 0: continue
        a, b, y0, y1 = fr(a), fr(b), fr(y0), fr(y1)
        for j in range(R, int(k) - R + 1):
            if j + b > k - R: continue
            for f in rots:
                if typ == 'h': A, B = f(j + a, y0), f(j + b, y0)
                elif typ == 'v': A, B = f(j + a, y0), f(j + a, y1)
                else: raise ValueError(typ)
                segs.append((A[0], A[1], B[0], B[1], ms))
    return segs


def write_cover(path_exact, k, out):
    import mixed_cover as MC
    m, EM, xe, d = load_exact(path_exact)
    segs = box_elements(m, EM, xe, k)
    a = EM.a; k = F(k)
    leb = (k - 2 * a) ** 2
    from math import lcm
    Dd = 1
    for s in segs:
        for v in s[:4]: Dd = lcm(Dd, v.denominator)
    Dd = lcm(Dd, a.denominator)
    Wd = 1
    for s in segs: Wd = lcm(Wd, s[4].denominator)
    Wd = lcm(Wd, leb.denominator)
    cover = dict(s_num=int(k), s_den=1, D=Dd, W=Wd, points=[],
                 segments=[(int(s[0] * Dd), int(s[1] * Dd), int(s[2] * Dd), int(s[3] * Dd), int(s[4] * Wd)) for s in segs],
                 polygons=[(int(leb * Wd), [(int(a * Dd), int(a * Dd)), (int((k - a) * Dd), int(a * Dd)),
                                              (int((k - a) * Dd), int((k - a) * Dd)), (int(a * Dd), int((k - a) * Dd))])])
    tot = sum(s[4] for s in segs) + leb
    D = F(d['D'])
    print(f"box k={k}: {len(segs)} segments + Lebesgue [{a},{k-a}]^2; total = {float(tot):.12f}; "
          f"k^2 - 4D = {float(k*k - 4*D):.12f}; equal: {tot == k*k - 4*D}")
    MC.write(out, cover, comment=f"qx2 box k={k} from {path_exact}; total {tot} = k^2 - 4D, D = {D}")
    MC.validate(MC.load(out))
    return cover


def write_family(path_exact, out):
    """human-readable exact description of the family (pi, nu) behind a sol_exact.json."""
    m, EM, xe, d = load_exact(path_exact)
    L = []
    L.append(f"# qx2 fixed-profile family, R = {m.R}, w = {m.w:g}, Lebesgue on U = [a, inf)^2 with a = {EM.a}; from {path_exact}")
    L.append(f"# D = {d['D']} = {float(F(d['D'])):.12f};  sigma = {d['sigma']}")
    L.append("# profile pi (one period, phases in [0,1), mirror orbits listed element by element; mass = total mass of the piece):")
    L.append("#   h  y  phase0 phase1  mass        (horizontal piece)      v  phase  y0 y1  mass   (vertical piece)")
    for (typ, a, b, y0, y1, var, fm) in m.base:
        ms = xe[int(var)]
        if ms == 0: continue
        if typ == 'h': L.append(f"h {fr(y0)} {fr(a)} {fr(b)} {ms}")
        else: L.append(f"v {fr(a)} {fr(y0)} {fr(y1)} {ms}")
    L.append(f"# plus Lebesgue (density 1) on [0,1) x [{EM.a}, {EM.w}] (per period mass {EM.w - EM.a})")
    L.append("# corner module nu on [0,R]^2 (diagonal orbits listed element by element): H y x0 x1 mass | V x y0 y1 mass")
    for kk in ('hs', 'vs'):
        for e, t in zip(m.A[kk], m.KA[kk]):
            if t != 0: continue
            ms = xe[int(e[-2])]
            if ms == 0: continue
            L.append(f"{'H' if kk == 'hs' else 'V'} {fr(e[0])} {fr(e[1])} {fr(e[2])} {ms}")
    L.append(f"# plus Lebesgue (density 1) on [{EM.a}, {m.R}]^2")
    open(out, 'w').write("\n".join(L) + "\n")
    print("wrote", out, len(L), "lines")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['resolve', 'project', 'axis', 'cover', 'family'])
    ap.add_argument('path'); ap.add_argument('--k', type=int, default=7); ap.add_argument('--out', default='')
    a = ap.parse_args()
    if a.cmd == 'resolve': resolve(a.path)
    elif a.cmd == 'project': project(a.path)
    elif a.cmd == 'axis': axis_check(a.path)
    elif a.cmd == 'cover': write_cover(a.path, a.k, a.out or a.path.replace('.json', f'_box{a.k}.txt'))
    elif a.cmd == 'family': write_family(a.path, a.out or a.path.replace('.json', '_family.txt'))
