"""Exponent calculator for c*(k) upper-bound constructions (orders only; constants and logs dropped).

Model (from Ryu, k2-minus-c-upper v1.1, Prop. 6.3, §7 Lemma stair, §8 Lemma W; see tasks/k2c-exponent/README.md):
  c*(k) ~ 4k/b + E(b),   E(b) = T(b, 1/b)     (band end: thin trapezoid, length b, slope ~2/b, height free)
T(L, z) is the cost of filling a thin trapezoid of length L whose ceiling has slope z against the fill's own
orientation, height free.  Fillers:
  * trivial (axis layers):  L          if L z >= 1;   L^2 z  if L z <= 1;
  * chain (Ryu's ZC; Kb = 1 is the §7 staircase): columns of height X at tilt t, Kb = max(1, L z / (t^2 X)) blocks;
      terms  L t                (column triangles)          L z X t^3      (ceiling-slope part of c_col)
             Kb X t             (aligned closings)          Kb / t^2       (unaligned closings)
             Kb h_T w_T         (transition slivers)        L z            (the rise; right end)
             wedge slivers w0, w_e (lower order)            2 T(X, t)      (the two ends: walls of height X, slope t)
      constraints  X >= 1/t^2 (W3), t <= 1, X >= 1;  if Kb > 1: X >= L t and L <= t X^2 (W2), t >= z (W1);
                   if Kb = 1: t^2 X >= L z (one staircase spans the length).
Every term is a monomial, so for a fixed choice of filler and case at each level, minimizing the largest exponent
is a linear program in (log b, log X_i, log t_i) / log k.  `solve` enumerates the cases up to a depth.

  python3 search/ryu/exponents.py [DMAX]      best exponent per depth, with the optimal parameters
  python3 search/ryu/exponents.py --drop LABELS [DMAX]   what-if: omit terms whose label contains a substring
  python3 search/ryu/exponents.py --scale 'LABEL:E,...' [DMAX]   what-if: multiply terms by t^E at their level
  python3 search/ryu/exponents.py --terms D CASES   e.g. --terms 2 S,C : binding terms of one case vector
"""
import sys, itertools
import numpy as np
from scipy.optimize import linprog


class Lin:
    """affine form in the LP variables: dict var -> coef, plus constant."""
    def __init__(self, d=None, c=0.0):
        self.d = dict(d or {}); self.c = c

    def __add__(self, o):
        o = o if isinstance(o, Lin) else Lin(c=o)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, 0) + v
        return Lin(d, self.c + o.c)
    __radd__ = __add__

    def __neg__(self):
        return Lin({k: -v for k, v in self.d.items()}, -self.c)

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Lin) else Lin(c=o)))

    def __rsub__(self, o):
        return Lin(c=o) - self

    def __rmul__(self, a):
        return Lin({k: a * v for k, v in self.d.items()}, a * self.c)


def V(name):
    return Lin({name: 1.0})


def build(cases):
    """cases: tuple over levels, each 'S' (staircase, Kb=1) or 'C' (chain, Kb>1); then a leaf 'T+' (Lz>=1) or
    'T-' (Lz<=1).  Returns (terms, constraints) as lists of (label, Lin); constraint Lin >= 0."""
    terms, cons = [], []
    B = V('B')
    terms.append(('body 4k/b', 1 - B))
    cons += [('b>=1', B), ('b<=k', 1 - B)]
    l, s = B, -B                     # level-1 region: length b, slope 1/b
    for i, cs in enumerate(cases[:-1]):
        X, t = V(f'X{i}'), V(f't{i}')
        start = len(terms)
        p = f'L{i} '
        cons += [(p + 't<=1', -t), (p + 'X>=1', X), (p + 'X>=1/t^2 (W3)', X + 2 * t)]
        if cs == 'C':
            kb = l + s - 2 * t - X
            cons += [(p + 'Kb>=1', kb), (p + 'X>=Lt (W2 drift)', X - l - t), (p + 'L<=tX^2 (W2)', t + 2 * X - l),
                     (p + 't>=z (W1)', t - s)]
        else:
            kb = Lin()
            cons += [(p + 't^2X>=Lz (one block)', 2 * t + X - l - s)]
        terms += [(p + 'column triangles Lt', l + t), (p + 'ceiling c_col LzXt^3', l + s + X + 3 * t),
                  (p + 'closings-aligned KbXt', kb + X + t),
                  (p + 'rise Lz', l + s)]
        hT = [Lin(), s - 2 * t, s + X + t]
        wT = [Lin(), X + 4 * t - s, l + 4 * t]
        for a, h in enumerate(hT):
            for c, w in enumerate(wT):
                terms.append((p + f'slivers Kb h{a} w{c}', kb + h + w))
        w0 = [t, X + 3 * t]
        we = [t, X + t, l + s + t, X + 3 * t]
        for a, w in enumerate(w0 + we):
            terms.append((p + f'wedge {a}', w))
            terms.append((p + f'wedge {a} zw^2', s + 2 * w))
        for j in range(start, len(terms)):        # what-if: extra power of t on chosen terms (--scale)
            lab, f = terms[j]
            for sub, e in SCALE.items():
                if sub in lab:
                    f = f + e * t
            terms[j] = (lab, f)
        if cs == 'C':                # in the one-block staircase (E38) the closing lies inside the right wall
            terms.append((p + 'closings-unaligned Kb/t^2', kb - 2 * t))
        l, s = X, t                  # the two ends: walls of length X, slope t
    leaf = cases[-1]
    if leaf == 'T+':
        terms.append(('leaf trivial L', l)); cons.append(('leaf Lz>=1', l + s))
    else:
        terms.append(('leaf trivial L^2 z', 2 * l + s)); cons.append(('leaf Lz<=1', -(l + s)))
    return terms, cons


SCALE = {}       # substring -> extra power of the level's tilt t (what-if experiments: --scale)
DROP = ()        # substrings of term labels to omit (what-if experiments: --drop)


def lp(cases):
    terms, cons = build(cases)
    terms = [(lab, f) for lab, f in terms if not any(d in lab for d in DROP)]
    names = sorted({k for _, f in terms + cons for k in f.d})
    idx = {n: i for i, n in enumerate(names)}
    nv = len(names) + 1               # last: z
    A, bnd = [], []
    for _, f in terms:                # f <= z
        row = np.zeros(nv); row[-1] = -1
        for k, v in f.d.items():
            row[idx[k]] += v
        A.append(row); bnd.append(-f.c)
    for _, f in cons:                 # f >= 0  ->  -f <= 0
        row = np.zeros(nv)
        for k, v in f.d.items():
            row[idx[k]] -= v
        A.append(row); bnd.append(f.c)
    cvec = np.zeros(nv); cvec[-1] = 1
    r = linprog(cvec, A_ub=np.array(A), b_ub=np.array(bnd), bounds=[(-5, 5)] * nv, method='highs')
    if r.status != 0:
        return None
    sol = dict(zip(names, r.x[:-1]))
    val = lambda f: sum(v * sol[k] for k, v in f.d.items()) + f.c
    binding = [lab for lab, f in terms if val(f) > r.x[-1] - 1e-7]
    return r.x[-1], sol, binding


def solve(dmax):
    out = {}
    for d in range(dmax + 1):
        best = None
        for cs in itertools.product('SC', repeat=d):
            for leaf in ('T+', 'T-'):
                res = lp(cs + (leaf,))
                if res and (best is None or res[0] < best[0] - 1e-9):
                    best = (res[0], cs + (leaf,), res[1], res[2])
        out[d] = best
    return out


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--terms':
        d = int(sys.argv[2]); cs = tuple(sys.argv[3].split(',')) if d else ()
        for leaf in ('T+', 'T-'):
            res = lp(cs + (leaf,))
            if res:
                print(leaf, 'exponent %.6f' % res[0], {k: round(v, 4) for k, v in res[1].items()})
                print('   binding:', res[2])
        sys.exit()
    if len(sys.argv) > 2 and sys.argv[1] == '--scale':
        SCALE = {a.split(':')[0]: float(a.split(':')[1]) for a in sys.argv[2].split(',')}
        sys.argv = sys.argv[2:]
    if len(sys.argv) > 2 and sys.argv[1] == '--drop':
        DROP = tuple(sys.argv[2].split(','))
        sys.argv = sys.argv[2:]
    dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    from fractions import Fraction
    for d, (e, cs, sol, binding) in solve(dmax).items():
        fr = Fraction(e).limit_denominator(200)
        print(f'depth {d}: exponent {e:.6f} (~{fr}) cases {cs}  b = k^{sol["B"]:.4f}')
        print('   params:', {k: round(v, 4) for k, v in sol.items()})
        print('   binding:', binding)
