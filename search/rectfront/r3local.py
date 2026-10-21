#!/usr/bin/env python3
"""R_3 lower bound, the tube: two tilts near 0 (mod 90), where r3bb.py cannot close boxes (2026-10-20).

On the zero lines (two tilts exactly 0) some types have margin exactly 0 at a corner of the excluded region, so no
fixed lambda certifies a box that meets them.  Here lambda depends on the tilts: lambda(t) keeps the support and
the zero force components of the LP optimum at a probe point near the line ("cancellation rule": the net force
components that vanish at the probe are kept at exactly 0 by solving for the secondary weights).  Then, on a box
whose near-zero tilts r range over one sign each (a quadrant, so |sin r| is smooth there):

    m(t) >= m(t_line) + sum_i r_i * [inf over the box of d m / d r_i]        (mean value along the segment)

with d m / d r_i by interval forward-mode AD (iad.py; Clarke hulls where a |.| changes sign).  m(t_line) >= 0 is
checked by intervals on the line segment when it is > TOL, else it must be one of the analytic classes:
  SBS  lambda on one pair, separated along an axis that is horizontal on the line, both squares upright: the
       margin at V1 is exactly 0 (S = 1, R = w - 1), at V2 it is 2 - (1 + 2 sqrt2)/2 > 0.
  WEDGE / STACK: recorded for a separate analytic lemma (closed forms in RECTFRONT.md).
Rigour: the derivative test needs sign-definite bounds (> 0 or < 0 strictly), padding as in iad.py.
"""
import itertools, json, math, os, sys
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from iad import I, D, dcos, dsin, dabs, dfsq, ddiv
import r3bb as B

R2 = math.sqrt(2)
V = B.V
RAD = math.pi / 180


def model(tilts, typ, lam_fn):
    """m_V for V1, V2 as D numbers.  tilts: list of D (radians).  lam_fn(u, S) -> list of 3 D weights."""
    n = len(tilts[0].g)
    us, Ss = [], []
    for (i, j), oi in zip(B.PAIRS, typ):
        o, m = B.OPTS[oi]
        own, oth = (i, j) if o == 0 else (j, i)
        c, sn = dcos(tilts[own]), dsin(tilts[own])          # quarter turns exactly: no sin(pi) != 0
        us.append([(c, sn), (-sn, c), (-c, -sn), (sn, -c)][m])
        Ss.append(dfsq(tilts[oth] - tilts[own]) + 0.5)
    k = [dfsq(t) for t in tilts]
    lam = lam_fn(us, Ss)
    for l in lam:                                   # weights must be >= 0 on the whole box
        if l is not None and l.enc().lo < 0:
            raise ZeroDivisionError('negative weight')
    out = []
    for (w, h) in V:
        m = D.const(0, n)
        for l, S in zip(lam, Ss):
            if l is not None: m = m + l * S
        for K in range(3):
            gx = D.const(0, n); gy = D.const(0, n)
            for p, ((i, j), (ux, uy)) in enumerate(zip(B.PAIRS, us)):
                if lam[p] is None: continue
                sg = 1 if K == j else (-1 if K == i else 0)
                if sg: gx = gx + lam[p] * ux * sg; gy = gy + lam[p] * uy * sg
            m = m - dabs(gx) * (k[K] * -1 + w / 2) - dabs(gy) * (k[K] * -1 + h / 2)
        out.append(m)
    return out


def G_coef(us, K, comp):
    """coefficients a_p with G_K,comp = sum_p lam_p a_p."""
    a = []
    for (i, j), u in zip(B.PAIRS, us):
        sg = 1 if K == j else (-1 if K == i else 0)
        a.append(u[comp] * sg if sg else None)
    return a


def probe_rule(pt_deg, typ):
    """LP at a probe point: support, main pair, zero force components."""
    lam, mar = B.margin_lp([(x, x) for x in pt_deg], typ)
    if lam is None: return None
    k, pairs = B.type_data(None, typ, pt_deg)
    sup = [p for p in range(3) if lam[p] > 1e-7]
    main = max(sup, key=lambda p: lam[p])
    zeros = []
    for K in range(3):
        for comp in (0, 1):
            g = sum((lam[p] * (pairs[p][2][0] if comp == 0 else pairs[p][3][0]) * (1 if K == pairs[p][1] else -1))
                    for p in sup if K in pairs[p][:2])
            involved = [p for p in sup if K in pairs[p][:2]]
            if abs(g) < 1e-8 and involved:
                zeros.append((K, comp))
    return dict(lam=lam, sup=sup, main=main, zeros=zeros, margin=mar)


def make_lam_fn(rule, rho=None, rho_probe=None):
    sup, main = rule['sup'], rule['main']
    others = [p for p in sup if p != main]
    lam0 = rule['lam']

    def const(q, n):
        c = lam0[q] / lam0[main]
        if rho is None:
            return D.const(c, n)
        return rho * (c / rho_probe)                  # secondary weight proportional to the distance from the line

    def fn(us, Ss):
        n = len(us[0][0].g)
        lam = [None] * 3
        lam[main] = D.const(1.0, n)
        if not others:
            return lam
        # equations G_{K,c} = 0 for the recorded zeros that involve an unknown weight
        eqs = []
        for (K, comp) in rule['zeros']:
            a = G_coef(us, K, comp)
            if any(a[p] is not None for p in others):
                eqs.append(a)
        if len(others) == 1 and eqs:
            q = others[0]; a = eqs[0]
            if a[q] is None:
                lam[q] = const(q, n)
            else:
                rhs = -a[main] if a[main] is not None else D.const(0.0, n)
                lam[q] = ddiv(rhs, a[q])
            return lam
        if len(others) == 2 and len(eqs) >= 2:
            a, b = eqs[0], eqs[1]
            q1, q2 = others
            z = D.const(0.0, n)
            a1, a2 = a[q1] or z, a[q2] or z; b1, b2 = b[q1] or z, b[q2] or z
            ra = -(a[main] or z); rb = -(b[main] or z)
            det = a1 * b2 - a2 * b1
            lam[q1] = ddiv(ra * b2 - a2 * rb, det); lam[q2] = ddiv(a1 * rb - ra * b1, det)
            return lam
        for q in others:                      # fallback: ratios from the probe, scaled by rho / rho_probe
            lam[q] = const(q, n)
        return lam
    return fn


def check_quadrant(box_deg, typ, radial, signs, rho_deg):
    """box_deg: 3 degree intervals (the radial ones are replaced by [0, rho] * sign).  Returns (ok, info)."""
    nr = len(radial)
    # probe: radial tilts at 0.3 * rho * sign, others at the box centre
    pt = [(lo + hi) / 2 for lo, hi in box_deg]
    for r, s in zip(radial, signs): pt[r] = 0.3 * rho_deg * s
    rule = probe_rule(pt, typ)
    if rule is None or rule['margin'] <= 0:
        return False, dict(why='probe LP not positive', rule=rule)
    fn = make_lam_fn(rule)
    # line segment value
    tl = []
    for K in range(3):
        if K in radial:
            tl.append(D.const(I(0.0), nr))
        else:
            tl.append(D.const(I(box_deg[K][0] * RAD, box_deg[K][1] * RAD), nr))
    try:
        mline = model(tl, typ, fn)
    except ZeroDivisionError:
        return False, dict(why='division by 0 on the line', rule=rule)
    line_lo = min(x.v.lo for x in mline)
    # derivative over the quadrant box
    tb = []
    for K in range(3):
        if K in radial:
            idx = radial.index(K); s = signs[idx]
            iv = I(0.0, rho_deg * RAD) if s > 0 else I(-rho_deg * RAD, 0.0)
            g = [I(0)] * nr; g[idx] = I(1)
            tb.append(D(iv, g))
        else:
            tb.append(D.const(I(box_deg[K][0] * RAD, box_deg[K][1] * RAD), nr))
    try:
        mb = model(tb, typ, fn)
    except ZeroDivisionError:
        return False, dict(why='division by 0 on the box', rule=rule)
    der_ok = all((x.g[i].lo > 0 if signs[i] > 0 else x.g[i].hi < 0) for x in mb for i in range(nr))
    # also allow: line value strictly positive and absorbs a negative derivative
    slack = []
    for x in mb:
        worst = 0.0
        for i in range(nr):
            d = x.g[i].lo if signs[i] > 0 else -x.g[i].hi       # directional derivative lower bound (per unit |r_i|)
            worst += min(0.0, d) * rho_deg * RAD
        slack.append(worst)
    info = dict(line_lo=line_lo, der=[[(x.g[i].lo, x.g[i].hi) for i in range(nr)] for x in mb],
                sup=rule['sup'], zeros=rule['zeros'])
    if line_lo > B.TOL and line_lo + min(slack) > B.TOL:
        return True, dict(info, how='line>0')
    if der_ok:
        return True, dict(info, how='derivative')   # needs m(line) >= 0: interval or analytic (see caller)
    return False, dict(info, why='derivative not sign-definite')


def check_sector(box_deg, typ, radial, phi, rho_deg, probe_frac=0.3):
    """Rays from the line: radial tilts = rho * (cos phi, sin phi), rho in [0, rho_deg], phi in the degree interval
    `phi` (a sector on which the |.| signs are fixed).  One AD variable (rho)."""
    from iad import icos, isin
    nr = len(radial)
    pm = math.radians((phi[0] + phi[1]) / 2)
    pt = [(lo + hi) / 2 for lo, hi in box_deg]
    dirs_pt = [math.cos(pm), math.sin(pm)]
    for r, d in zip(radial, dirs_pt): pt[r] = probe_frac * rho_deg * d
    rule = probe_rule(pt, typ)
    if rule is None or rule['margin'] <= 0:
        return False, dict(why='probe LP not positive', rule=rule)
    rule['rho_probe'] = probe_frac * rho_deg * RAD
    ph = I(math.radians(phi[0]), math.radians(phi[1]))
    for variant in ('cancel', 'const'):
        r = _sector_eval(box_deg, typ, radial, rho_deg, rule, ph, variant)
        if r[0] or r[1].get('why') not in ('division by 0',):
            return r
    return r


def _sector_eval(box_deg, typ, radial, rho_deg, rule, ph, variant):
    from iad import icos, isin
    if variant == 'const':
        rule = dict(rule, zeros=[])                    # ratios from the probe, scaled with rho
    rho_d_line = D(I(0.0), [I(0)], I(0.0))
    rho_d_box = D(I(0.0, rho_deg * RAD), [I(1)], I(0.0))
    rp = rule.get('rho_probe', 0.3 * rho_deg * RAD)
    dvec = [icos(ph), isin(ph)]
    def tilts(rho_iv, deriv):
        out = []
        for K in range(3):
            if K in radial:
                d = dvec[radial.index(K)]
                out.append(D(rho_iv * d, [d if deriv else I(0)], I(0.0)))
            else:
                out.append(D.const(I(box_deg[K][0] * RAD, box_deg[K][1] * RAD), 1))
        return out
    import iad
    try:
        iad.RHO = None
        mline = model(tilts(I(0.0), False), typ, make_lam_fn(rule, rho_d_line, rp) if variant == 'const' else make_lam_fn(rule))
        iad.RHO = I(0.0, rho_deg * RAD)
        mb = model(tilts(I(0.0, rho_deg * RAD), True), typ, make_lam_fn(rule, rho_d_box, rp) if variant == 'const' else make_lam_fn(rule))
    except ZeroDivisionError:
        return False, dict(why='division by 0', rule=rule)
    line_lo = min(x.v.lo for x in mline)
    info = dict(line_lo=line_lo, line=[x.v.lo for x in mline], sup=rule['sup'], zeros=rule['zeros'],
                der=[(x.g[0].lo, x.g[0].hi) for x in mb])
    hows = []
    for ml, x in zip(mline, mb):                       # per vertex of the excluded region
        lv, d = ml.v.lo, x.g[0].lo
        if lv > B.TOL and lv + min(0.0, d) * rho_deg * RAD > B.TOL:
            hows.append('line>0')
        elif lv >= 0 and d > 0:
            hows.append('derivative')                  # lv >= 0 must hold exactly: see the caller / analytic classes
        else:
            return False, dict(info, why='vertex not certified')
    return True, dict(info, how=hows)


if __name__ == '__main__':
    box = [(0, 0.5), (19, 21), (0, 0.5)]
    for typ in [(5, 2, 3), (1, 0, 3)]:
        for ph in [(0, 45), (45, 90), (90, 180), (180, 225), (225, 270), (270, 360)]:
            ok, info = check_sector(box, typ, [0, 2], ph, 0.5)
            print(typ, ph, ok, info.get('how') or info.get('why'), 'line %.3g' % info.get('line_lo', float('nan')),
                  'der', [tuple(round(a, 4) for a in d) for d in info.get('der', [])])
    raise SystemExit

    box = [(0, 0.5), (19, 21), (0, 0.5)]
    for typ in [(5, 2, 3), (5, 2, 7), (1, 0, 3)]:
        for signs in itertools.product((1, -1), repeat=2):
            ok, info = check_quadrant(box, typ, [0, 2], signs, 0.5)
            print(typ, signs, ok, info.get('how') or info.get('why'), 'line', '%.3g' % info.get('line_lo', float('nan')),
                  'der', [[round(a, 4) for a in d[0]] for d in info.get('der', [])][:2])


# ---------------- driver over r3bb's local boxes ----------------
def sectors_for(box, radial):
    """45-degree sectors of the (r1, r2) plane meeting the box's rectangle with positive area."""
    (a0, a1), (b0, b1) = box[radial[0]], box[radial[1]]
    out = []
    for k in range(8):
        lo, hi = math.radians(k * 45), math.radians((k + 1) * 45)
        # sample the sector: does any small ray point lie strictly inside the rectangle?
        hit = False
        for f in (0.1, 0.5, 0.9):
            ph = lo + f * (hi - lo)
            for r in (1e-3, 0.3, 0.7, 1.0, 1.4):
                x, y = r * math.cos(ph), r * math.sin(ph)
                if a0 < x < a1 and b0 < y < b1: hit = True
        if hit: out.append((k * 45.0, (k + 1) * 45.0))
    return out


def certify_box_type(box, typ, eps, beta_step=0.5, min_sector=0.125):
    near = [K for K in range(3) if max(abs(box[K][0]), abs(box[K][1])) <= eps + 1e-12]
    if len(near) != 2:
        return False, dict(why=f'{len(near)} near-zero tilts (needs the 3-D origin lemma)')
    radial = near; T = [K for K in range(3) if K not in near][0]
    rho0 = math.hypot(*[max(abs(box[K][0]), abs(box[K][1])) for K in radial])
    lo, hi = box[T]
    nb = max(1, math.ceil((hi - lo) / beta_step - 1e-9))
    stats = dict(checks=0, hows={})
    for ib in range(nb):
        b = (lo + (hi - lo) * ib / nb, lo + (hi - lo) * (ib + 1) / nb)
        sbox = list(box); sbox[T] = b
        stack = sectors_for(box, radial)
        while stack:
            ph = stack.pop()
            ok, info = check_sector(sbox, typ, radial, ph, rho0)
            stats['checks'] += 1
            if ok:
                for hw in info['how']: stats['hows'][hw] = stats['hows'].get(hw, 0) + 1
                continue
            if ph[1] - ph[0] <= min_sector:
                return False, dict(why=info.get('why'), beta=b, phi=ph, info=info, stats=stats)
            m = (ph[0] + ph[1]) / 2
            stack += [(ph[0], m), (m, ph[1])]
    return True, stats


def run_local(path, eps, out, limit=None):
    rows = [json.loads(l) for l in open(path)]
    if limit: rows = rows[:limit]
    res = []
    for r in rows:
        box = [tuple(x) for x in r['box']]
        for typ in r['fail']:
            ok, info = certify_box_type(box, tuple(typ), eps)
            res.append(dict(box=box, typ=typ, ok=ok, info={k: v for k, v in info.items() if k != 'info'},
                            detail=(info.get('info') or {}).get('der') if not ok else None))
            print(box, typ, ok, info.get('why', ''), info.get('stats', info).get('checks') if isinstance(info, dict) else '', flush=True)
    with open(out, 'w') as f:
        for x in res: f.write(json.dumps(x, default=str) + '\n')
    print('total', len(res), 'failed', sum(not x['ok'] for x in res))


def polar_cell_box(sbox, radial, ph, rh):
    """Cartesian bounding box (degrees) of the polar cell rho in rh (deg), phi in ph (deg)."""
    from iad import icos, isin
    p = I(math.radians(ph[0]), math.radians(ph[1]))
    c, s = icos(p), isin(p)
    r = I(rh[0], rh[1])
    xa, ya = r * c, r * s
    b = list(sbox)
    b[radial[0]] = (xa.lo, xa.hi); b[radial[1]] = (ya.lo, ya.hi)
    return b


def certify_box_type2(box, typ, eps, beta_step=0.5, min_sector=0.02, min_rho=1e-5, max_checks=20000):
    near = [K for K in range(3) if max(abs(box[K][0]), abs(box[K][1])) <= eps + 1e-12]
    if len(near) != 2:
        return False, dict(why=f'{len(near)} near-zero tilts')
    radial = near; T = [K for K in range(3) if K not in near][0]
    rho0 = math.hypot(*[max(abs(box[K][0]), abs(box[K][1])) for K in radial])
    lo, hi = box[T]
    nb = max(1, math.ceil((hi - lo) / beta_step - 1e-9))
    stats = dict(checks=0, deriv=0, cell=0)
    stack = []
    for ib in range(nb):
        b = (lo + (hi - lo) * ib / nb, lo + (hi - lo) * (ib + 1) / nb)
        for ph in sectors_for(box, radial):
            stack.append((b, ph, (0.0, rho0)))
    while stack:
        b, ph, rh = stack.pop()
        stats['checks'] += 1
        if stats['checks'] > max_checks:
            return False, dict(why='too many checks', stats=stats, beta=b, phi=ph, rho=rh)
        sbox = list(box); sbox[T] = b
        if rh[0] == 0.0:
            ok, info = check_sector(sbox, typ, radial, ph, rh[1])
            if ok: stats['deriv'] += 1; continue
        else:
            cb = polar_cell_box(sbox, radial, ph, rh)
            lam, mc = B.margin_lp(cb, typ)
            if lam is not None and B.margin_box(cb, typ, lam) > B.TOL:
                stats['cell'] += 1; continue
            info = dict(why='cell margin')
        # split: the longest of (sector, rho) relative to their scale; beta last
        wph = ph[1] - ph[0]; wr = rh[1] - rh[0]
        if wph <= min_sector and wr <= min_rho:
            if b[1] - b[0] > 0.01:
                m = (b[0] + b[1]) / 2
                stack += [((b[0], m), ph, rh), ((m, b[1]), ph, rh)]; continue
            return False, dict(why=info.get('why'), beta=b, phi=ph, rho=rh, info=info, stats=stats)
        if rh[0] == 0.0 and wr > min_rho and (wph <= min_sector or wr * RAD * 50 > wph * math.pi / 180 * max(rh[1], 1e-9)):
            m = rh[1] / 4
            stack += [(b, ph, (0.0, m)), (b, ph, (m, rh[1]))]
        elif wph > min_sector:
            m = (ph[0] + ph[1]) / 2
            stack += [(b, (ph[0], m), rh), (b, (m, ph[1]), rh)]
        else:
            m = (rh[0] + rh[1]) / 2
            stack += [(b, ph, (rh[0], m)), (b, ph, (m, rh[1]))]
    return True, stats
