"""Exact form of a packing from its contact structure: number field, exact configuration, minimal polynomial of S.

Input: exactsolve's NAME.exact.txt (configuration, ~70 digits) and NAME.contacts.json (closed contacts used as
equations, free squares).  Output: NAME.minpoly.json, checked independently by verify_exact.py.

Method.
* Angle classes: squares in the system whose angles agree mod 90 deg, or are opposite mod 90 deg, to 1e-40 share one
  parameter t = tan(theta/2), theta in (-45, 45] deg: square i has (c, s) = R^m (c(u), s(u)) with u = +-t,
  c(u) = (1-u^2)/(1+u^2), s(u) = 2u/(1+u^2), R = rotation by 90 deg.  Axis-parallel squares: u = 0.
* With the angles fixed, every contact equation is linear in the centres and S.
* Flat directions (null space of the contact Jacobian in (x, y, S, t)) are pinned at nearby rationals, centres
  preferred over class angles.  S in the null space = S not fixed by the contacts (force-determined): not handled.
* Sparse elimination over Q(t) (pivots: nonzero at t*).  Rows left without a usable pivot vanish at t*; f = the
  irreducible factor of the gcd of their numerators with root t*.  Then back-substitution in K = Q[t]/f, and every
  equation is checked to be 0 in K.  The minimal polynomial of S is the first linear dependence among 1, S, S^2, ...
* Free squares (rattlers): centres and t rounded to rationals (they have clearance >= ~1e-6).

Only one free (unpinned) tilted class for now; more classes need multivariate elimination (task doc, roadmap 2).

  python3 minpoly.py DIR/NAME   (reads DIR/NAME.exact.txt, DIR/NAME.contacts.json; writes DIR/NAME.minpoly.json)
"""
import sys, json, time
from fractions import Fraction
import math
import flint
from mpmath import mp, mpf
import geom

Q = flint.fmpq_poly
T = Q([0, 1])


class Ring:
    """Q[t] (k <= 1: fmpq_poly) or Q[t1..tk] (fmpq_mpoly); the ring of the free class parameters."""
    def __init__(self, k):
        self.k = k
        if k <= 1:
            self.gens = [Q([0, 1])]
            self.one, self.zero = Q([1]), Q([0])
        else:
            self.ctx = flint.fmpq_mpoly_ctx.get(tuple(f't{i + 1}' for i in range(k)), 'lex')
            self.gens = list(self.ctx.gens())
            self.one, self.zero = self.ctx.from_dict({(0,) * k: 1}), self.ctx.from_dict({})

    def const(self, x):
        x = Fraction(x)
        q = flint.fmpq(x.numerator, x.denominator)
        return Q([q]) if self.k <= 1 else self.ctx.from_dict({(0,) * self.k: q} if x else {})

    def ev(self, p, tv):
        if self.k <= 1:
            return ev_poly(p, tv[0])
        v = mpf(0)
        for e, c in p.to_dict().items():
            m = mpf(int(c.p)) / int(c.q)
            for ti, ei in zip(tv, e):
                if ei:
                    m *= ti ** int(ei)
            v += m
        return v


R = Ring(1)


class NotHandled(Exception):
    pass


# ---------------------------------------------------------------- rational functions in t, linear forms over them
class RF:
    """num/den over the ring R, den with leading coefficient 1 (gcd-normalised)."""
    __slots__ = ('n', 'd')

    def __init__(self, n, d=None, norm=True):
        d = R.one if d is None else d
        if norm and not d.is_one():
            if n.is_zero():
                d = R.one
            else:
                g = n.gcd(d)
                if not g.is_one():
                    n, d = n / g, d / g     # exact division (fmpq_poly '/' by a divisor)
                lc = d.leading_coefficient()
                n, d = n / lc, d / lc
        self.n, self.d = n, d

    @staticmethod
    def c(x):
        return RF(R.const(x))

    def __add__(s, o):
        o = o if isinstance(o, RF) else RF.c(o)
        if s.d == o.d:
            return RF(s.n + o.n, s.d)
        return RF(s.n * o.d + o.n * s.d, s.d * o.d)
    __radd__ = __add__

    def __neg__(s):
        return RF(-s.n, s.d, False)

    def __sub__(s, o):
        return s + (-(o if isinstance(o, RF) else RF.c(o)))

    def __rsub__(s, o):
        return (-s) + o

    def __mul__(s, o):
        if isinstance(o, Lin):
            return o * s
        o = o if isinstance(o, RF) else RF.c(o)
        return RF(s.n * o.n, s.d * o.d)
    __rmul__ = __mul__

    def __truediv__(s, o):
        return RF(s.n * o.d, s.d * o.n)

    def zero(s):
        return s.n.is_zero()

    def ev(s, tv):
        return R.ev(s.n, tv) / R.ev(s.d, tv)


def ev_poly(p, t):
    c = p.coeffs()
    v = mpf(0)
    for a in reversed(c):
        v = v * t + mpf(int(a.p)) / int(a.q)
    return v


class Lin:
    """Affine form sum_v coef[v] * v + const, coefficients RF."""
    __slots__ = ('c', 'k')

    def __init__(self, c=None, k=None):
        self.c = c or {}
        self.k = k if k is not None else RF(R.zero)

    @staticmethod
    def var(v):
        return Lin({v: RF(R.one)})

    def __add__(s, o):
        if not isinstance(o, Lin):
            return Lin(dict(s.c), s.k + o)
        c = dict(s.c)
        for v, a in o.c.items():
            c[v] = c[v] + a if v in c else a
        return Lin({v: a for v, a in c.items() if not a.zero()}, s.k + o.k)
    __radd__ = __add__

    def __neg__(s):
        return Lin({v: -a for v, a in s.c.items()}, -s.k)

    def __sub__(s, o):
        return s + (-o)

    def __rsub__(s, o):
        return (-s) + o

    def __mul__(s, o):
        if isinstance(o, Lin):
            raise TypeError('product of two linear forms')
        o = o if isinstance(o, RF) else RF.c(o)
        if o.zero():
            return Lin()
        return Lin({v: a * o for v, a in s.c.items()}, s.k * o)
    __rmul__ = __mul__


# ---------------------------------------------------------------- K = Q[t]/f
class K:
    """K = Q[T]/f; ring elements (polynomials in the free class parameters) map in by t_i -> tmap[i]."""
    def __init__(self, f, tmap):
        self.f, self.tmap = f, tmap
        self.pw = [[Q([1])] for _ in tmap]

    def tpow(self, i, e):
        while len(self.pw[i]) <= e:
            self.pw[i].append(self.red(self.pw[i][-1] * self.tmap[i]))
        return self.pw[i][e]

    def to_K(self, p):
        if R.k <= 1:
            return p % self.f                                      # t = T (or a constant when k = 0)
        acc = Q([0])
        for e, c in p.to_dict().items():
            m = Q([c])
            for i, ei in enumerate(e):
                if ei:
                    m = self.red(m * self.tpow(i, ei))
            acc += m
        return self.red(acc)

    def red(self, p):
        return p % self.f

    def inv(self, p):
        g, a, b = (p % self.f).xgcd(self.f)
        assert g.degree() == 0, 'not invertible mod f'
        return self.red(a / g.coeffs()[0])

    def rf(self, r):
        return self.red(self.to_K(r.n) * self.inv(self.to_K(r.d)))


def cs_rat(u):
    u = Fraction(u)
    return (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)


def rot90(c, s, m):
    for _ in range(m % 4):
        c, s = -s, c
    return c, s


def frac_near(x, digits):
    return Fraction(int(mp.nint(x * mpf(10) ** digits)), 10 ** digits)


def run(base, quiet=False):
    t0 = time.time()
    log = (lambda *a: None) if quiet else (lambda *a: print(*a, file=sys.stderr))
    ct = json.load(open(base + '.contacts.json'))
    try:
        rep = json.load(open(base + '.json'))                      # exactsolve's report
    except FileNotFoundError:
        rep = {}
    if rep and not rep.get('cert_valid'):
        raise NotHandled(f"exactsolve's point is not a valid packing ({len(rep.get('violations') or [])} overlapping "
                         f"pairs, certificate invalid): re-solve first")
    L = [l.split() for l in open(base + '.exact.txt') if l.strip()]
    n = int(L[0][0])
    mp.dps = max(30, len(L[0][1]) - 2)
    S0 = mpf(L[0][1])
    X = [mpf(l[0]) for l in L[1:n + 1]]
    Y = [mpf(l[1]) for l in L[1:n + 1]]
    TH = [mpf(l[2]) for l in L[1:n + 1]]                          # degrees
    eqs = [tuple(c) for c in ct['equations']]
    # keep only equations that hold at the exact point (exactsolve's list can contain active-set incidences of free
    # squares that were moved afterwards for clearance)
    Cm0 = [mp.cos(th * mp.pi / 180) for th in TH]
    Sm0 = [mp.sin(th * mp.pi / 180) for th in TH]
    held = [c for c in eqs if abs(geom.eval_contact(c, X, Y, Cm0, Sm0, S0, mpf(1) / 2, order=1)[0]) < mpf(10) ** -40]
    if len(held) < len(eqs):
        log(f'  dropped {len(eqs) - len(held)} listed equations that do not hold at the point: '
            f'{[c for c in eqs if c not in held][:4]}')
    eqs = held
    sysq = sorted({c[1] for c in eqs} | {c[3] for c in eqs if c[0] == 'C'})
    free = [i for i in range(n) if i not in sysq]
    tolang = mpf(10) ** -40
    # corner-corner touches (exactsolve sets them aside as disjunctive): every incidence of a touching pair that holds
    # at the exact point (|g| < 1e-50, corner within the side's closed extent) is an equation too
    Cm = [mp.cos(th * mp.pi / 180) for th in TH]
    Sm = [mp.sin(th * mp.pi / 180) for th in TH]
    have = set(eqs)
    nvv = 0
    for i, j in ct.get('vv', []):
        if i not in sysq or j not in sysq:
            continue
        for own, oth in ((i, j), (j, i)):
            for k in range(4):
                for a in range(4):
                    c = ('C', oth, a, own, k)
                    g = geom.eval_contact(c, X, Y, Cm, Sm, S0, mpf(1) / 2, order=1)[0]
                    if abs(g) < mpf(10) ** -50 and abs(geom.tangential(c, X, Y, Cm, Sm, mpf(1) / 2)) <= \
                            mpf(1) / 2 + mpf(10) ** -50 and c not in have:
                        eqs.append(c)
                        have.add(c)
                        nvv += 1

    # angle classes over the squares in the system: class 0 = axis (u = 0)
    reps = [mpf(0)]                                                # representative angle, degrees, in (-45, 45]
    memb = {}                                                      # i -> (class, sign, m)
    for i in sysq:
        th = TH[i]
        r = th - 90 * mp.floor(th / 90)
        if r > 45:
            r -= 90                                                # r in (-45, 45]
        for k, rk in enumerate(reps):
            hit = None
            for sg in (1, -1):
                d = (r - sg * rk) / 90
                if abs(d - mp.nint(d)) < tolang:
                    hit = sg
                    break
            if hit:
                break
        else:
            k, hit = len(reps), 1
            reps.append(r)
        m = int(mp.nint((th - hit * reps[k]) / 90)) % 4
        memb[i] = (k, hit if k else 1, m)
    tstar = [mp.tan(r * mp.pi / 360) for r in reps]
    ncls = len(reps) - 1
    log(f'n = {n}: {len(sysq)} squares in the system, {len(free)} free; {ncls} tilted angle classes; '
        f'{len(eqs)} contact equations ({nvv} from corner-corner touches)')

    # class angles that exactsolve froze as exact flat directions: pinned at a nearby rational t
    frozen_sq = {v // 3 for v in ct.get('frozen_vars', []) if v % 3 == 2 and v // 3 in memb}
    pins = {}
    for i in frozen_sq:
        k = memb[i][0]
        if k and f't{k}' not in pins:
            if any(memb[j][0] == k for j in sysq if j not in frozen_sq):
                raise NotHandled(f'square {i}: flat rotation, but parallel to squares that do not rotate')
            pins[f't{k}'] = frac_near(tstar[k], 30)
    kts = [k for k in range(1, ncls + 1) if f't{k}' not in pins]     # free tilted classes
    global R
    R = Ring(len(kts))
    tv = [tstar[k] for k in kts] if kts else [mpf(0)]

    def u_cls(k):
        if k == 0:
            return Fraction(0)
        if k in kts:
            return RF(R.gens[kts.index(k)])
        return pins[f't{k}']
    cls_cs = {}
    for k in range(ncls + 1):
        u = u_cls(k)
        if isinstance(u, RF):
            d = RF(R.one) + u * u
            cls_cs[k] = ((RF(R.one) - u * u) / d, (RF.c(2) * u) / d)
        else:
            c, s = cs_rat(u)
            cls_cs[k] = (RF.c(c), RF.c(s))
    Cs, Ss = {}, {}
    for i in sysq:
        k, sg, m = memb[i]
        c, s = cls_cs[k]
        Cs[i], Ss[i] = rot90(c, s if sg > 0 else -s, m)
    half = RF.c(Fraction(1, 2))

    def build_rows(eqlist):
        Xs = {i: Lin(k=RF.c(pins[f'x{i}'])) if f'x{i}' in pins else Lin.var(f'x{i}') for i in sysq}
        Ys = {i: Lin(k=RF.c(pins[f'y{i}'])) if f'y{i}' in pins else Lin.var(f'y{i}') for i in sysq}
        return [contact_expr(c, Xs, Ys, Cs, Ss, Lin.var('S'), half) for c in eqlist]

    allv = {f'{a}{i}' for i in sysq for a in 'xy'} | {'S'}

    def free_vars(rows, piv_rows):
        # every unpinned centre coordinate counts, including ones that occur in no equation (their coefficient is
        # identically 0: e.g. a square held only by side contacts with horizontal normals is free vertically)
        return sorted((allv - set(pins)) - {v for v, _ in piv_rows})

    def consistent(left):
        thr = mpf(10) ** -50
        for r in left:
            if not kts:
                if r.c or not r.k.zero():
                    return False
            elif abs(r.k.ev(tv)) > thr:
                return False
        return True

    # sparse elimination over Q(t), pivots nonzero at t*.  Variables left without a pivot are free (given t the
    # system is linear).  Settle them first: add near-contacts (|g| < 1e-9 at the point, not yet equations) one at a
    # time when that removes a free variable and keeps the system consistent at t*; pin what is left at nearby
    # rationals and eliminate again.
    rows = build_rows(eqs)
    piv_rows, leftover = eliminate(rows, tv)
    freev = free_vars(rows, piv_rows)
    settled = []
    if freev:
        if 'S' in freev:
            raise NotHandled('S is not fixed by the contacts (force-determined S): not handled yet')
        dep = s_dependence(piv_rows, freev, tv)
        if dep:
            raise NotHandled(f'S depends on centres not fixed by the contacts ({dep[:4]}): force-determined, not handled')
        cands = near_incidences(sysq, X, Y, Cm, Sm, S0, mpf(10) ** -9, set(eqs))
        Xs0 = {i: Lin.var(f'x{i}') for i in sysq}
        Ys0 = {i: Lin.var(f'y{i}') for i in sysq}
        progress = True
        while freev and progress:
            progress = False
            fs = set(freev)
            for gabs, c in cands:
                if c in settled:
                    continue
                g = contact_expr(c, Xs0, Ys0, Cs, Ss, Lin.var('S'), half)
                if not (set(g.c) & fs):
                    continue
                rows2 = rows + [g]
                p2, l2 = eliminate(rows2, tv)
                f2 = free_vars(rows2, p2)
                if len(f2) < len(freev) and 'S' not in f2 and consistent(l2):
                    eqs.append(c)
                    settled.append(c)
                    rows, piv_rows, leftover, freev = rows2, p2, l2, f2
                    progress = True
                    break
        if settled:
            log(f'  settled {len(settled)} flat directions on near-contacts (|g| <= '
                f'{mp.nstr(max(g for g, c in cands if c in settled), 2)} at the point)')
    if freev:
        for v in freev:
            pins[v] = frac_near((X if v[0] == 'x' else Y)[int(v[1:])], 30)
        rows = build_rows(eqs)
        piv_rows, leftover = eliminate(rows, tv)
        if free_vars(rows, piv_rows):
            raise RuntimeError('variables still free after pinning')
    if pins:
        log(f'  pinned (flat directions): {", ".join(f"{v} = {float(q):.12g}" for v, q in pins.items())}')
    log(f'  elimination: {len(piv_rows)} pivots, {len(leftover)} consistency rows')
    stationary = False
    if not kts:
        if any(r.c or not r.k.zero() for r in leftover):
            raise RuntimeError('inconsistent rational system')
        f, tmap, troot = Q([0, 1]), [], mpf(0)                     # K = Q (T := 0)
    elif len(kts) == 1:
        g = R.zero
        for r in leftover:
            for a in list(r.c.values()) + [r.k]:
                if not a.zero():
                    g = a.n if g.is_zero() else g.gcd(a.n)
        if g.is_zero():
            # the contacts leave t free: S(t) along the family, and t* is a stationary point (force balance)
            vrf = {}
            for v, pr in reversed(piv_rows):
                acc = pr.k
                for w, a in pr.c.items():
                    if w != v:
                        acc = acc + a * vrf[w]
                vrf[v] = -acc / pr.c[v]
            Sr = vrf['S']
            g = Sr.n.derivative() * Sr.d - Sr.n * Sr.d.derivative()
            stationary = True
            log(f'  t not fixed by the contacts: S(t) of degree {Sr.n.degree()}/{Sr.d.degree()}; t* from dS/dt = 0')
        if g.is_zero() or g.degree() < 1:
            raise NotHandled('no consistency polynomial and S constant in t')
        gz = flint.fmpz_poly([int(x) for x in g.numer().coeffs()])
        cand = []
        for fac, mult in gz.factor()[1]:
            vals = abs(ev_poly(Q(fac), tv[0])) / max(1, max(abs(int(x)) for x in fac.coeffs()))
            cand.append((vals, fac))
        cand.sort(key=lambda z: z[0])
        if not (cand[0][0] < mpf(10) ** -(mp.dps // 2)) or (len(cand) > 1 and cand[1][0] < mpf(10) ** -10):
            raise RuntimeError(f'factor choice ambiguous: {[mp.nstr(c[0], 3) for c in cand[:3]]}')
        f, tmap, troot = Q(cand[0][1]), [T], tv[0]
        log(f'  f(t): degree {f.degree()}, from gcd of degree {g.degree()} ({len(cand)} factors)')
    else:
        f, tmap, troot, stationary = multi_field(leftover, tv, base, log, piv_rows)
    Kf = K(f, tmap)
    # back-substitution in K
    val = {}
    for v, pr in reversed(piv_rows):
        acc = Kf.rf(pr.k)
        for w, a in pr.c.items():
            if w != v:
                acc = Kf.red(acc + Kf.rf(a) * val[w])
        val[v] = Kf.red(-acc * Kf.inv(Kf.rf(pr.c[v])))
    # exact check of every equation in K
    for c, r in zip(eqs, rows):
        acc = Kf.rf(r.k)
        for w, a in r.c.items():
            acc = Kf.red(acc + Kf.rf(a) * val[w])
        if not acc.is_zero():
            raise RuntimeError(f'equation {c} not satisfied in K')
    Sk = val['S']
    # minimal polynomial of S over Q: first dependence among powers
    d = f.degree()
    pw = [Q([1])]
    p_S = None
    for e in range(1, d + 1):
        pw.append(Kf.red(pw[-1] * Sk))
        M = flint.fmpq_mat(e + 1, d, [pw[i].coeffs()[j] if j < len(pw[i].coeffs()) else 0
                                      for i in range(e + 1) for j in range(d)])
        if M.rank() <= e:
            # solve sum a_i pw_i = 0 with a_e = 1
            A = flint.fmpq_mat(d, e, [M[i, j] for j in range(d) for i in range(e)])
            b = flint.fmpq_mat(d, 1, [-M[e, j] for j in range(d)])
            # least squares-free: pick e independent rows
            sol = solve_consistent(A, b)
            co = [sol[i] for i in range(e)] + [Fraction(1)]
            den = 1
            for x in co:
                den = den * Fraction(x).denominator // math.gcd(den, Fraction(x).denominator)
            p_S = [int(Fraction(x) * den) for x in co]
            break
    pz = flint.fmpz_poly(p_S)
    g0 = 0
    for x in p_S:
        g0 = math.gcd(g0, abs(x))
    pz = flint.fmpz_poly([x // g0 for x in p_S])
    assert len(pz.factor()[1]) == 1 and pz.factor()[1][0][1] == 1, 'minimal polynomial not irreducible'
    Sv_num = ev_poly(Sk, troot)
    log(f'  S: minimal polynomial degree {pz.degree()}, height {max(abs(int(x)) for x in pz.coeffs())}; '
        f'S(t*) - S_file = {mp.nstr(Sv_num - S0, 3)}')
    if abs(Sv_num - S0) > mpf(10) ** -(mp.dps // 2):
        raise RuntimeError('S(t*) does not match the file')

    pS_list = [int(x) for x in pz.coeffs()]

    def kstr(p):
        return [str(Fraction(int(x.p), int(x.q))) for x in p.coeffs()]
    aS, bS = isolate(pS_list, Sv_num)
    out = {
        'n': n, 'name': base.split('/')[-1],
        'field': {'f': [int(x * f.denom()) for x in f.coeffs()] if kts else [0, 1],
                  't_interval': [str(x) for x in (isolate([int(x * f.denom()) for x in f.coeffs()], troot) if kts else
                                                  (Fraction(0), Fraction(0)))],
                  't_approx': mp.nstr(troot, 50)},
        'S': {'poly': pS_list, 'interval': [str(aS), str(bS)], 'value': mp.nstr(Sv_num, 50),
              'in_K': [str(Fraction(int(x.p), int(x.q))) for x in Sk.coeffs()]},
        'classes': [{'theta_deg': mp.nstr(r, 40),
                     'u': kstr(tmap[kts.index(k)]) if k in kts else [str(u_cls(k))]} for k, r in enumerate(reps)],
        'stationary': stationary,
        'squares': [],
        'pins': {v: str(q) for v, q in pins.items()},
        'equations': [list(c) for c in eqs],
        'stats': {'n_sys': len(sysq), 'n_free': len(free), 'n_classes_tilted': ncls, 'field_degree': f.degree() if kts else 1, 'n_classes_free': len(kts),
                  'S_degree': pz.degree(), 'seconds': None},
    }

    def coord(v):
        return [str(pins[v])] if v in pins else kstr(val[v])
    for i in range(n):
        if i in memb:
            k, sg, m = memb[i]
            out['squares'].append({'x': coord(f'x{i}'), 'y': coord(f'y{i}'), 'm': m,
                                   'u': kstr(tmap[kts.index(k)] * sg) if k in kts else [str(sg * Fraction(u_cls(k)))]})
        else:
            th = TH[i] * mp.pi / 180
            u = frac_near(mp.tan(th / 2), 15) if abs(mp.tan(th / 2)) <= 1 else None
            m = 0
            if u is None:                                          # rotate by 90 so that |tan(theta/2)| <= 1
                th2 = th - mp.pi / 2
                u, m = frac_near(mp.tan(th2 / 2), 15), 1
            out['squares'].append({'x': [str(frac_near(X[i], 15))], 'y': [str(frac_near(Y[i], 15))],
                                   'u': [str(u)], 'm': m, 'free': True})
    out['stats']['seconds'] = round(time.time() - t0, 2)
    log(f'  done in {out["stats"]["seconds"]} s')
    with open(base + '.minpoly.json', 'w') as fh:
        json.dump(out, fh)
    return out


def msolve_param(polys, k, tag):
    """Rational parametrization of the zero-dimensional system polys (fmpq_mpoly in t1..tk) by msolve:
    returns (f, den, [(g_i, c_i)]) with t_i = -g_i(A) / (c_i den(A)), f(A) = 0 (fmpq_poly)."""
    import subprocess, tempfile, os, re
    names = [f't{i + 1}' for i in range(k)]

    def mstr(p):
        d = p.to_dict()
        den = 1
        for c in d.values():
            den = den * int(c.q) // math.gcd(den, int(c.q))
        out = ''
        for e, c in d.items():                                     # plain syntax only: msolve misreads "(c)*..."
            cz = int(c * den)
            mon = '*'.join(f'{nm}^{ei}' if ei > 1 else nm for nm, ei in zip(names, e) if ei)
            out += ('-' if cz < 0 else '+') + str(abs(cz)) + (f'*{mon}' if mon else '')
        return out.lstrip('+') or '0'
    with tempfile.TemporaryDirectory() as td:
        fin, fout = os.path.join(td, 'in.ms'), os.path.join(td, 'out.ms')
        with open(fin, 'w') as fh:
            fh.write(','.join(names) + '\n0\n' + ',\n'.join(mstr(p) for p in polys) + '\n')
        r = subprocess.run(['msolve', '-P', '1', '-t', '1', '-f', fin, '-o', fout], capture_output=True, text=True,
                           timeout=3600)
        txt = open(fout).read() if os.path.exists(fout) else ''
    if not txt.strip():
        raise RuntimeError(f'msolve failed: {r.stderr[-300:]}')
    txt = re.sub(r'(-?\d+)\s*/\s*2\^(\d+)', r'Fraction(\1, 2**\2)', txt.strip().rstrip(':'))
    d = eval(txt, {'Fraction': Fraction})
    if d[0] != 0:
        raise NotHandled(f'msolve: system not zero-dimensional (code {d[0]})')
    body = d[1]
    vars_ = body[3]
    par = body[5][1]
    P = lambda c: Q([int(x) for x in c[1]]) if c[0] >= 0 else Q([0])
    f, den = P(par[0]), P(par[1])
    params = [(P(g), int(c)) for g, c in par[2]]
    if len(vars_) - 1 != len(params):
        raise RuntimeError('unexpected msolve output shape')
    # the last variable is the separating element A (msolve may add it, or reorder t1..tk); the others are
    # parametrized in the order listed
    if sorted(v for v in vars_ if v in names) != sorted(names):
        raise RuntimeError(f'unexpected msolve variables {vars_}')
    byname = dict(zip(vars_[:-1], params))
    byname[vars_[-1]] = (-(Q([0, 1]) * den), 1)                    # A itself: -g/(c den) = A
    return f, den, [byname[nm] for nm in names]


def back_rf(piv_rows):
    """Back-substitution over the rational function field: every pivot variable as an RF."""
    vrf = {}
    for v, pr in reversed(piv_rows):
        acc = pr.k
        for w, a in pr.c.items():
            if w != v:
                acc = acc + a * vrf[w]
        vrf[v] = -acc / pr.c[v]
    return vrf


def mdet(M):
    """Determinant of a small square matrix of fmpq_mpoly (cofactor expansion)."""
    if len(M) == 1:
        return M[0][0]
    acc = None
    for j in range(len(M)):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = M[0][j] * mdet(minor)
        acc = term if acc is None else (acc + term if j % 2 == 0 else acc - term)
    return acc


def multi_field(leftover, tv, base, log, piv_rows):
    """K for k >= 2 free class parameters: the consistency polynomials (numerators of the rows left after
    elimination), their irreducible factors that vanish at t*, msolve's parametrization, and the factor of f and root
    that reproduce t*.  Returns (f1, tmap, A*)."""
    k = R.k
    keep = []
    thr = mpf(10) ** -(mp.dps // 2)
    for r in leftover:
        for a in list(r.c.values()) + [r.k]:
            if a.zero():
                continue
            _, facs = a.n.factor()
            for fac, _m in facs:
                if fac.is_constant():
                    continue
                h = max(abs(float(c)) for c in fac.coeffs())
                if abs(R.ev(fac, tv)) / max(1.0, h) < thr and not any((fac - q).is_zero() or (fac + q).is_zero()
                                                                        for q in keep):
                    keep.append(fac)
    log(f'  {len(leftover)} consistency rows -> {len(keep)} distinct irreducible factors vanishing at t*: degrees '
        f'{[q.total_degree() for q in keep]}')
    stationary = False
    if len(keep) == k - 1:
        # the contacts leave a curve of angles; S is stationary along it at t* (force balance): Lagrange condition
        # det [grad C_1; ...; grad C_{k-1}; grad S] = 0, with the numerator of grad S
        Sr = back_rf(piv_rows)['S']
        names = [f't{i + 1}' for i in range(k)]
        gS = [Sr.n.derivative(nm) * Sr.d - Sr.n * Sr.d.derivative(nm) for nm in names]
        M = [[q.derivative(nm) for nm in names] for q in keep] + [gS]
        det = mdet(M)
        if det.is_zero():
            raise NotHandled('Lagrange determinant vanishes identically')
        _, facs = det.factor()
        new = []
        for fac, _m in facs:
            if fac.is_constant():
                continue
            h = max(abs(float(c)) for c in fac.coeffs())
            if abs(R.ev(fac, tv)) / max(1.0, h) < thr:
                new.append(fac)
        log(f'  angles fixed by force balance along a curve: Lagrange determinant of total degree {det.total_degree()}, '
            f'{len(new)} factors vanishing at t*: degrees {[q.total_degree() for q in new]}')
        if not new:
            raise RuntimeError('no factor of the Lagrange determinant vanishes at t*')
        keep += new
        stationary = True
    if len(keep) < k:
        raise NotHandled(f'{len(keep)} consistency polynomials for {k} class angles: more than one angle direction left '
                         f'to force balance, not handled')
    t0 = time.time()
    f, den, params = msolve_param(keep, k, base)
    log(f'  msolve: eliminating polynomial of degree {f.degree()} ({time.time() - t0:.1f} s)')
    fz = flint.fmpz_poly([int(x * f.denom()) for x in f.coeffs()])
    best = []
    prec0 = flint.ctx.prec
    flint.ctx.prec = int(3.4 * mp.dps) + 64                       # certified root balls at the working precision
    roots = [(fac, rt) for fac, _m in fz.factor()[1] for rt, _ in fac.complex_roots()]
    flint.ctx.prec = prec0
    for fac, rt in roots:
        if True:
            if abs(float(rt.imag.mid())) > 1e-6:
                continue
            A = mpf(rt.real.mid().str(mp.dps, radius=False, more=True))
            ts = [-ev_poly(g, A) / (c * ev_poly(den, A)) for g, c in params]
            err = max(abs(a - b) for a, b in zip(ts, tv))
            best.append((err, fac, A))
    best.sort(key=lambda z: z[0])
    if not best or best[0][0] > mpf(10) ** -30 or (len(best) > 1 and best[1][0] < mpf(10) ** -10):
        raise RuntimeError(f'root choice ambiguous: {[mp.nstr(b[0], 3) for b in best[:3]]}')
    _, fac, A = best[0]
    f1 = Q(fac)
    Kf = K(f1, [])
    dinv = Kf.inv(den % f1)
    tmap = [Kf.red(-g * dinv / c) for g, c in params]
    log(f'  field: degree {f1.degree()} (factor of the eliminating polynomial with root t*)')
    return f1, tmap, A, stationary


def near_incidences(sysq, X, Y, Cm, Sm, S0, tol, have):
    """Incidences (corner on a side line within the side's closed extent, or corner on a wall) with |g| < tol at the
    numerical point, not in `have`; sorted by |g|."""
    H = mpf(1) / 2
    out = []
    ss = set(sysq)
    for j in sysq:
        for a in range(4):
            for w in 'LRBT':
                c = ('W', j, a, w)
                g = geom.eval_contact(c, X, Y, Cm, Sm, S0, H, order=1)[0]
                if abs(g) < tol and c not in have:
                    out.append((abs(g), c))
    for i in sysq:
        for j in sysq:
            if i == j or (X[i] - X[j]) ** 2 + (Y[i] - Y[j]) ** 2 > 2.5:
                continue
            for k in range(4):
                for a in range(4):
                    c = ('C', j, a, i, k)
                    if c in have:
                        continue
                    g = geom.eval_contact(c, X, Y, Cm, Sm, S0, H, order=1)[0]
                    if abs(g) < tol and abs(geom.tangential(c, X, Y, Cm, Sm, H)) <= H + tol:
                        out.append((abs(g), c))
    out.sort(key=lambda z: z[0])
    return out


def s_dependence(piv_rows, freev, tv):
    """Free variables that S depends on (numerically at t = tv): back-substitution of unit vectors in f64."""
    dep = []
    for fv in freev:
        val = {fv: 1.0}
        for v, pr in reversed(piv_rows):
            acc = 0.0
            for w, a in pr.c.items():
                if w != v:
                    acc += float(a.ev(tv)) * val.get(w, 0.0)
            val[v] = -acc / float(pr.c[v].ev(tv))
        if abs(val.get('S', 0.0)) > 1e-9:
            dep.append(fv)
    return dep


def eliminate(rows, tv):
    """Sparse Gaussian elimination of affine forms over Q(t); a pivot must be nonzero at t = tv.  Returns the pivot
    rows (var, row) in order and the rows left without a usable pivot (all their coefficients vanish at tv)."""
    thr = mpf(10) ** -(mp.dps // 2)
    piv_rows, active = [], [r for r in rows if r.c or not r.k.zero()]
    while active:
        cnt = {}
        for r in active:
            for v in r.c:
                cnt[v] = cnt.get(v, 0) + 1
        best = None
        for ri, r in enumerate(active):
            for v, a in r.c.items():
                cost = (len(r.c) - 1) * (cnt[v] - 1)
                if best is not None and cost >= best[0]:
                    continue
                if abs(a.ev(tv)) > thr:
                    best = (cost, ri, v)
        if best is None:
            return piv_rows, active
        _, ri, v = best
        pr = active.pop(ri)
        a = pr.c[v]
        piv_rows.append((v, pr))
        nxt = []
        for r in active:
            if v in r.c:
                r = r - pr * (r.c[v] / a)
                r.c.pop(v, None)
            if r.c or not r.k.zero():
                nxt.append(r)
        active = nxt
    return piv_rows, []


def isolate(p, x0, w=Fraction(1, 10 ** 25)):
    """Rational [a, b] containing exactly one real root of the squarefree integer polynomial p, near x0: exact sign
    change p(a) p(b) < 0, and every other root (flint's certified complex root balls) farther than w from x0."""
    mid = frac_near(x0, 40)
    a, b = mid - w, mid + w
    pv = lambda x: sum(Fraction(c) * x ** i for i, c in enumerate(p))
    if not pv(a) * pv(b) < 0:
        raise RuntimeError('no sign change on the isolating interval')
    far = 0
    for r, _ in flint.fmpz_poly(p).complex_roots():
        d = abs(r - flint.acb(flint.arb(str(mid.numerator)) / flint.arb(str(mid.denominator))))
        if d > flint.arb(str(2 * w.numerator)) / flint.arb(str(w.denominator)):
            far += 1
    if far != len(p) - 2:
        raise RuntimeError(f'isolation failed: {len(p) - 1 - far} roots within 2w')
    return a, b


def solve_consistent(A, b):
    """A x = b, A (d x e) of full column rank, consistent: solve on e independent rows."""
    d, e = A.nrows(), A.ncols()
    rows = []
    for i in range(d):
        trial = rows + [i]
        M = flint.fmpq_mat(len(trial), e, [A[r, j] for r in trial for j in range(e)])
        if M.rank() == len(trial):
            rows = trial
            if len(rows) == e:
                break
    As = flint.fmpq_mat(e, e, [A[r, j] for r in rows for j in range(e)])
    bs = flint.fmpq_mat(e, 1, [b[r, 0] for r in rows])
    x = As.solve(bs)
    xs = [Fraction(int(x[i, 0].p), int(x[i, 0].q)) for i in range(e)]
    for i in range(d):                                             # consistency
        assert sum(Fraction(int(A[i, j].p), int(A[i, j].q)) * xs[j] for j in range(e)) == \
            Fraction(int(b[i, 0].p), int(b[i, 0].q))
    return xs


def contact_expr(c, X, Y, C, Sn, side, H):
    """g of a contact as Lin (geom's formula, specialised)."""
    if c[0] == 'W':
        _, i, a, w = c
        ox, oy = geom.offset(C[i], Sn[i], a, H)
        px, py = X[i] + ox, Y[i] + oy
        return {'L': px, 'R': side - px, 'B': py, 'T': side - py}[w]
    _, j, a, i, k = c
    nx, ny = geom.normal(C[i], Sn[i], k)
    ox, oy = geom.offset(C[j], Sn[j], a, H)
    dx, dy = X[j] + ox - X[i], Y[j] + oy - Y[i]
    return dx * nx + dy * ny - H


if __name__ == '__main__':
    for b in sys.argv[1:]:
        try:
            run(b)
        except NotHandled as ex:
            print(f'{b}: not handled: {ex}', file=sys.stderr)
            sys.exit(2)
