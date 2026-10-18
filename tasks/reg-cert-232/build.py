#!/usr/bin/env python3
"""Exact certificate of the regularized (3 angle groups) drawing of the n = 232 record, over L = K(sqrt 2).

  1. reg.regularize(232, style='coincide'): the display regularization at 80 digits (angles 0, theta1, 45 mod 90)
  2. contacts of that drawing (wall and face-normal gaps < 1e-40 at 80 digits)
  3. the same contact equations with exact coefficients in L = Q[w]/g (field.py: w = t + sqrt 2), solved exactly
     for the 464 centre coordinates; every contact equation (not only the chosen independent ones) checked exactly
  4. written in minpoly.py's format (field f := g, u per square in L) for the unchanged verify_exact.py

  search/exact/.venv/bin/python3 build.py      (writes n-232-reg3.minpoly.json.gz, build.json)
"""
import gzip, json, os, sys
from fractions import Fraction as Fr
import numpy as np
import scipy.linalg
import mpmath as mpm
from mpmath import mp, mpf
from flint import fmpq_poly, fmpq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'search', 'regularize'))
import reg  # noqa: E402

mp.dps = 80
FLD = json.load(open(os.path.join(HERE, 'field.json')))
G = fmpq_poly([int(c) for c in FLD['g']])
TW = fmpq_poly([fmpq(Fr(x).numerator, Fr(x).denominator) for x in FLD['t_in_w']])
RW = fmpq_poly([fmpq(Fr(x).numerator, Fr(x).denominator) for x in FLD['sqrt2_in_w']])
WSTAR = mpm.findroot(lambda x: sum(int(c) * x ** i for i, c in enumerate(FLD['g'])), mpf(FLD['wstar']))
CERT = json.load(gzip.open(os.path.join(ROOT, 'search', 'exact', 'cur', 'data', 'n-232.minpoly.json.gz'), 'rt'))


def red(a):
    return a % G


def inv(a):
    d, s, _ = a.xgcd(G)          # s a + _ g = d (constant)
    assert d.degree() == 0
    return red(s / d[0])


def num(a):
    """numeric value at w* (mp)"""
    return sum(mpf(int(c.p)) / int(c.q) * WSTAR ** i for i, c in enumerate(a.coeffs()))


def from_t(coeffs):
    """K element given as coefficient list in t -> L (polynomial in w)"""
    r, tp = fmpq_poly([0]), fmpq_poly([1])
    for c in coeffs:
        q = Fr(c)
        r = red(r + tp * fmpq(q.numerator, q.denominator))
        tp = red(tp * TW)
    return r


ONE, ZERO, HALF = fmpq_poly([1]), fmpq_poly([0]), fmpq_poly([fmpq(1, 2)])
T = TW
den = inv(red(ONE + T * T))
FRAMES = {'A': (ONE, ZERO), 'B': (red((ONE - T * T) * den), red(2 * T * den)), 'C': (red(RW * HALF), red(RW * HALF))}
UVAL = {'A': ZERO, 'B': T, 'C': red(RW - 1)}


def main():
    P, R, rep = reg.regularize(232, style='coincide', log=lambda *a: None)
    assert rep['verify']['ok'] and rep['groups_after'] == 3, rep['verify']
    n, S = 232, R['S']
    Sx = from_t(CERT['S']['in_K'])
    assert abs(num(Sx) - S) < mpf('1e-60'), (num(Sx), S)
    # angle classes
    th1 = 2 * mpm.atan(num(T)) * 180 / mp.pi
    cls = []
    for tdeg in R['T']:
        a = reg.fold(tdeg)                       # mod 90 into (-45, 45]
        for k, ref in (('A', mpf(0)), ('B', th1), ('C', mpf(45)), ('C', mpf(-45))):
            if abs(reg.fold(a - ref)) < mpf('1e-40'):
                cls.append(k); break
        else:
            raise SystemExit(f'angle {a} in no class')
    from collections import Counter
    print('classes', Counter(cls))
    # numeric frames per square: the class frame (any rotation by 90 deg gives the same square)
    CF = [tuple(num(v) for v in FRAMES[k]) for k in cls]
    X, Y = R['X'], R['Y']
    # check the class frame really is the drawing's square (angle equal mod 90)
    for i in range(n):
        c, s = reg.mp_frame(R['T'][i])
        ok = any(abs(c * cc + s * ss - 1) < mpf('1e-60') for cc, ss in ((CF[i][0], CF[i][1]), (-CF[i][1], CF[i][0]),
                                                                       (-CF[i][0], -CF[i][1]), (CF[i][1], -CF[i][0])))
        assert ok, i

    def h(i, nx, ny):                            # support of square i in direction (nx, ny), numeric
        c, s = CF[i]
        return (abs(nx * c + ny * s) + abs(-nx * s + ny * c)) / 2

    TOL = mpf('1e-40')
    rows = []                                    # ('w', i, d, side) | ('p', i, j, own, k)  (normal k of square own)
    for i in range(n):
        for d, Z in ((0, X), (1, Y)):
            e = h(i, 1, 0) if d == 0 else h(i, 0, 1)
            if abs(Z[i] - e) < TOL:
                rows.append(('w', i, d, 0))
            if abs(S - e - Z[i]) < TOL:
                rows.append(('w', i, d, 1))
    minsep = mpf(10)
    for i in range(n):
        for j in range(i + 1, n):
            if (X[i] - X[j]) ** 2 + (Y[i] - Y[j]) ** 2 > mpf('2.01'):
                continue
            best = None
            for own in (i, j):
                c, s = CF[own]
                for k, (nx, ny) in enumerate(((c, s), (-s, c), (-c, -s), (s, -c))):
                    oth = j if own == i else i
                    gap = nx * (X[oth] - X[own]) + ny * (Y[oth] - Y[own]) - mpf(1) / 2 - h(oth, nx, ny)
                    if best is None or gap > best[0]:
                        best = (gap, own, k)
                    if abs(gap) < TOL:
                        rows.append(('p', i, j, own, k))
            minsep = min(minsep, best[0])
    print('rows', len(rows), 'min separation', mpm.nstr(minsep, 5))
    assert minsep > -TOL

    # exact rows: sum coef * z + const = 0, z index 2i+d
    def exact_row(rw):
        if rw[0] == 'w':
            _, i, d, side = rw
            c, s = FRAMES[cls[i]]
            cn, sn = CF[i]
            # support along e_d: (|c| + |s|)/2 with the numeric signs
            e = red((c * (1 if cn >= 0 else -1) + s * (1 if sn >= 0 else -1)) * HALF)
            if side == 0:
                return {2 * i + d: ONE}, red(-e)                          # z - e = 0
            return {2 * i + d: fmpq_poly([-1])}, red(from_t(CERT['S']['in_K']) - e)   # S - e - z = 0
        _, i, j, own, k = rw
        oth = j if own == i else i
        c, s = FRAMES[cls[own]]
        nx, ny = [(c, s), (-s, c), (-c, -s), (s, -c)][k]
        nx, ny = red(nx), red(ny)
        nxn, nyn = num(nx), num(ny)
        co, so = FRAMES[cls[oth]]
        a1, a2 = red(nx * co + ny * so), red(-nx * so + ny * co)
        sg = lambda v: (1 if num(v) > 0 else -1) if abs(num(v)) > mpf('1e-60') else 0
        ho = red((a1 * sg(a1) + a2 * sg(a2)) * HALF)
        coef = {}
        for q, sgn in ((oth, 1), (own, -1)):
            for d, v in ((0, nx), (1, ny)):
                coef[2 * q + d] = red(coef.get(2 * q + d, ZERO) + v * sgn)
        coef = {kk: v for kk, v in coef.items() if not v.is_zero()}
        return coef, red(-HALF - ho)
    ER = [exact_row(rw) for rw in rows]
    # numeric consistency of exact rows at the drawing
    Z = [None] * (2 * n)
    for i in range(n):
        Z[2 * i], Z[2 * i + 1] = X[i], Y[i]
    worst = max(abs(sum(num(v) * Z[k] for k, v in cf.items()) + num(c0)) for cf, c0 in ER)
    print('numeric residual of exact rows at the 80-digit drawing', mpm.nstr(worst, 5))
    assert worst < mpf('1e-50')
    # independent subset by QR with column pivoting on J^T
    J = np.zeros((len(ER), 2 * n))
    for r, (cf, _) in enumerate(ER):
        for k, v in cf.items():
            J[r, k] = float(num(v))
    _, Rq, piv = scipy.linalg.qr(J.T, pivoting=True, mode='economic')
    rank = int(np.sum(np.abs(np.diag(Rq)) > 1e-9))
    print('rank', rank, 'of', 2 * n)
    assert rank == 2 * n
    sel = list(piv[:rank])
    # exact sparse elimination (Markowitz-ish: process rows with fewest unknowns first)
    eqs = [(dict(ER[r][0]), ER[r][1]) for r in sel]
    sol = {}
    pivots = []                                   # (var, coef dict, const) with coef[var] = 1 after scaling
    active = eqs
    while active:
        active.sort(key=lambda e: len(e[0]))
        cf, c0 = active.pop(0)
        assert cf, 'dependent row'
        # pivot: the variable with largest |coef| numerically among those in the row
        var = max(cf, key=lambda k: abs(num(cf[k])))
        ip = inv(cf[var])
        cf = {k: red(v * ip) for k, v in cf.items()}
        c0 = red(c0 * ip)
        pivots.append((var, cf, c0))
        new = []
        for cf2, c2 in active:
            if var in cf2:
                m = cf2.pop(var)
                for k, v in cf.items():
                    if k == var:
                        continue
                    t = red(cf2.get(k, ZERO) - m * v)
                    if t.is_zero():
                        cf2.pop(k, None)
                    else:
                        cf2[k] = t
                c2 = red(c2 - m * c0)
            new.append((cf2, c2))
        active = new
    for var, cf, c0 in reversed(pivots):
        v = -c0
        for k, a in cf.items():
            if k != var:
                v = v - a * sol[k]
        sol[var] = red(v)
    # every contact row exactly
    for cf, c0 in ER:
        tot = c0
        for k, v in cf.items():
            tot = red(tot + v * sol[k])
        assert tot.is_zero(), 'a contact row fails exactly'
    print('all', len(ER), 'contact rows hold exactly in L')
    dev = max(abs(num(sol[k]) - Z[k]) for k in range(2 * n))
    print('max |exact - 80-digit drawing|', mpm.nstr(dev, 5))
    assert dev < mpf('1e-50')
    # write minpoly format
    st = lambda a: [str(Fr(int(c.p), int(c.q))) for c in a.coeffs()] or ['0']
    lo = mpm.mpf(WSTAR) - mpf('1e-45'); hi = mpm.mpf(WSTAR) + mpf('1e-45')
    tl = Fr(int(mpm.floor(lo * 10 ** 50)), 10 ** 50); th = Fr(int(mpm.ceil(hi * 10 ** 50)), 10 ** 50)
    D = dict(n=n, name='n-232-reg3', field=dict(f=[int(c) for c in FLD['g']], t_interval=[str(tl), str(th)],
                                                  t_approx=mpm.nstr(WSTAR, 70),
                                                  note='w = t + sqrt 2, t the field generator of n-232.minpoly.json'),
             S=dict(poly=CERT['S']['poly'], interval=CERT['S']['interval'], value=CERT['S']['value'], in_K=st(Sx)),
             classes=[dict(name='A', theta_deg='0', u=st(UVAL['A'])),
                      dict(name='B', theta_deg=mpm.nstr(th1, 40), u=st(UVAL['B'])),
                      dict(name='C', theta_deg='45', u=st(UVAL['C']))],
             squares=[dict(x=st(sol[2 * i]), y=st(sol[2 * i + 1]), u=st(UVAL[cls[i]]), m=0, cls=cls[i]) for i in range(n)],
             stats=dict(contacts=len(ER), classes=dict(Counter(cls))))
    out = os.path.join(HERE, 'n-232-reg3.minpoly.json.gz')
    with gzip.open(out, 'wt') as fh:
        json.dump(D, fh)
    json.dump(dict(rep={k: rep[k] for k in ('g', 'groups_after', 'verify')}, rows=len(ER), classes=dict(Counter(cls)),
                   max_dev=mpm.nstr(dev, 5)), open(os.path.join(HERE, 'build.json'), 'w'), indent=1, default=str)
    print('wrote', out)


if __name__ == '__main__':
    main()
