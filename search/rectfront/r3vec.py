#!/usr/bin/env python3
"""Batched tube checks for R_3 (numpy version of r3local.check_sector; 2026-10-21).

Cells: (beta interval of the third tilt, sector phi interval, rho interval [ra, rb]) for one (box, partial type) and
one weight rule (support / main pair / zero force components / variant).  Same tests as r3local.check_sector:
per vertex V of the excluded region, (line value > TOL and absorbs a negative slope over the segment) or
(line value >= 0 and slope > 0); weights >= 0; no division by an interval containing 0.  Cells failing every rule
go back to the scalar driver.
"""
import math
import numpy as np
import viad as V
import r3bb as B

RAD = math.pi / 180


def _vdiv_masked(x, y, bad):
    """ddiv with a per-cell mask: cells where y contains 0 are marked bad (their denominator set to 1)."""
    z0 = (y.v.lo <= 0) & (y.v.hi >= 0)
    z1 = (y.v0.lo <= 0) & (y.v0.hi >= 0)
    e = y.enc(); z2 = (e.lo <= 0) & (e.hi >= 0)
    badc = z0 | z1 | z2
    bad |= badc
    one = V.VI.c(1.0, len(badc))
    def fix(iv):
        return V.VI(np.where(badc, 1.0, iv.lo), np.where(badc, 1.0, iv.hi))
    y2 = V.VD(fix(y.v), y.g, fix(y.v0), y.R)
    return V.ddiv(x, y2)


def model_vec(tilts, typ, rule, variant, rho_d, rho_probe):
    """tilts: list of 3 VD (radians).  Returns ([m_V1, m_V2] as VD, bad mask)."""
    n = len(tilts[0].v.lo)
    bad = np.zeros(n, bool)
    us, Ss = [], []
    for (i, j), oi in zip(B.PAIRS, typ):
        if oi is None:
            z = V.VD.c(0.0, n); us.append((z, z)); Ss.append(z); continue
        o, m = B.OPTS[oi]
        own, oth = (i, j) if o == 0 else (j, i)
        c, s = V.dcos(tilts[own]), V.dsin(tilts[own])
        us.append([(c, s), (-s, c), (-c, -s), (s, -c)][m])
        Ss.append(V.dfsq(tilts[oth] - tilts[own]) + 0.5)
    k = [V.dfsq(t) for t in tilts]
    sup, main = rule['sup'], rule['main']
    others = [p for p in sup if p != main]
    lam0 = rule['lam']
    lam = [None] * 3
    lam[main] = V.VD.c(1.0, n, tilts[0].R)

    def const(q):
        c = lam0[q] / lam0[main]
        if variant == 'const':
            return rho_d * (c / rho_probe)
        return V.VD.c(c, n, tilts[0].R)

    def coef(K, comp):
        a = []
        for (i, j), u in zip(B.PAIRS, us):
            sg = 1 if K == j else (-1 if K == i else 0)
            a.append(u[comp] * sg if sg else None)
        return a
    zeros = rule['zeros'] if variant == 'cancel' else []
    eqs = [coef(K, comp) for (K, comp) in zeros]
    eqs = [a for a in eqs if any(a[p] is not None for p in others)]
    if others:
        if len(others) == 1 and eqs:
            q = others[0]; a = eqs[0]
            if a[q] is None:
                lam[q] = const(q)
            else:
                rhs = -a[main] if a[main] is not None else V.VD.c(0.0, n, tilts[0].R)
                lam[q] = _vdiv_masked(rhs, a[q], bad)
        elif len(others) == 2 and len(eqs) >= 2:
            a, b = eqs[0], eqs[1]; q1, q2 = others
            z = V.VD.c(0.0, n, tilts[0].R)
            a1, a2 = a[q1] or z, a[q2] or z; b1, b2 = b[q1] or z, b[q2] or z
            ra = -(a[main] or z); rb = -(b[main] or z)
            det = a1 * b2 - a2 * b1
            lam[q1] = _vdiv_masked(ra * b2 - a2 * rb, det, bad); lam[q2] = _vdiv_masked(a1 * rb - ra * b1, det, bad)
        else:
            for q in others: lam[q] = const(q)
    lam = [None if typ[p] is None else l for p, l in enumerate(lam)]
    for l in lam:
        if l is not None:
            bad |= l.enc().lo < 0
    out = []
    for (w, h) in B.V:
        m = V.VD.c(0.0, n, tilts[0].R)
        for l, S in zip(lam, Ss):
            if l is not None: m = m + l * S
        for K in range(3):
            gx = V.VD.c(0.0, n, tilts[0].R); gy = V.VD.c(0.0, n, tilts[0].R)
            for p, ((i, j), (ux, uy)) in enumerate(zip(B.PAIRS, us)):
                if lam[p] is None: continue
                sg = 1 if K == j else (-1 if K == i else 0)
                if sg: gx = gx + lam[p] * ux * sg; gy = gy + lam[p] * uy * sg
            m = m - V.dabs(gx) * (k[K] * -1 + w / 2) - V.dabs(gy) * (k[K] * -1 + h / 2)
        out.append(m)
    return out, bad


def check_cells(box, typ, radial, T, cells, rule, variant):
    """cells: array (N, 6): beta lo, beta hi, phi lo, phi hi (deg), rho lo, rho hi (deg).  Returns ok mask."""
    n = len(cells)
    blo, bhi, plo, phi_, ra, rb = (cells[:, i] for i in range(6))
    P = V.VI(np.radians(plo), np.radians(phi_))
    dvec = [V.vcos(P), V.vsin(P)]
    Rrel = V.VI(np.zeros(n), (rb - ra) * RAD)
    Ra = V.VI(ra * RAD, ra * RAD)
    Rab = V.VI(ra * RAD, rb * RAD)

    def tilts(line):
        out = []
        for K in range(3):
            if K in radial:
                d = dvec[radial.index(K)]
                if line:
                    out.append(V.VD(Ra * d, V.VI.c(0.0, n), Ra * d, None))
                else:
                    out.append(V.VD(Rab * d, d, Ra * d, Rrel))
            else:
                iv = V.VI(blo * RAD, bhi * RAD)
                out.append(V.VD(iv, V.VI.c(0.0, n), iv, None if line else Rrel))
        return out
    rp = rule['rho_probe']
    rho_line = V.VD(Ra, V.VI.c(0.0, n), Ra, None)
    rho_box = V.VD(Rab, V.VI.c(1.0, n), Ra, Rrel)
    mline, bad1 = model_vec(tilts(True), typ, rule, variant, rho_line, rp)
    mb, bad2 = model_vec(tilts(False), typ, rule, variant, rho_box, rp)
    ok = ~(bad1 | bad2)
    span = (rb - ra) * RAD
    for ml, x in zip(mline, mb):
        lv = ml.v.lo; d = x.g.lo
        okv = ((lv > B.TOL) & (lv + np.minimum(0.0, d) * span > B.TOL)) | ((lv >= 0) & (d > 0))
        ok &= okv
    return ok


def _scalar_cells(box, typ, radial, T, cells, st, max_checks, min_sector=0.25, min_rho=0.02):
    """The adaptive scalar loop of r3local.certify_box_type3, started from the given cells."""
    import r3local as L
    stack = [((c[0], c[1]), (c[2], c[3]), (c[4], c[5])) for c in cells]
    while stack:
        b, ph, rh = stack.pop()
        st['checks'] += 1
        if st['checks'] > max_checks:
            return False, dict(why='too many checks', stats=st, beta=b, phi=ph, rho=rh)
        sbox = list(box); sbox[T] = b
        if L.polar_cell_box(sbox, radial, ph, rh) is None:
            continue
        ok, info = L.check_sector_multi(sbox, typ, radial, ph, rh[1], rho_lo=rh[0])
        if ok: st['scalar'] += 1; continue
        cb = L.polar_cell_box(sbox, radial, ph, rh)
        lam, mc = B.margin_lp(cb, typ)
        if lam is not None and B.margin_box(cb, typ, lam) > B.TOL:
            st['cell'] += 1; continue
        if ph[1] - ph[0] > min_sector:
            m = (ph[0] + ph[1]) / 2
            stack += [(b, (ph[0], m), rh), (b, (m, ph[1]), rh)]; continue
        if b[1] - b[0] > 0.01:
            m = (b[0] + b[1]) / 2
            stack += [((b[0], m), ph, rh), ((m, b[1]), ph, rh)]; continue
        if rh[1] - rh[0] > min_rho:
            m = rh[0] + (rh[1] - rh[0]) / (4 if rh[0] == 0 else 2)
            stack += [(b, ph, (rh[0], m)), (b, ph, (m, rh[1]))]; continue
        return False, dict(why='cell unresolved', stats=st, beta=b, phi=ph, rho=rh)
    return True, st


def certify_box_type_vec(box, typ, eps, beta_step=0.1, sector0=5.0, rhos=(0.0, 0.3, 0.6, 1.0), max_checks=8000,
                         force_T=None):
    import r3local as L
    near = [K for K in range(3) if max(abs(box[K][0]), abs(box[K][1])) <= eps + 1e-12]
    if force_T is not None:
        near = [K for K in range(3) if K != force_T]
    elif len(near) == 3 or (len(near) == 2 and all(box[K][0] <= 0.0 <= box[K][1] for K in range(3))):
        return certify_origin_vec(box, typ, eps, max_checks)
    if len(near) != 2:
        return L.certify_box_type3(box, typ, eps, max_checks=max_checks)
    radial = near; T = [K for K in range(3) if K not in near][0]
    rho0 = math.hypot(*[max(abs(box[K][0]), abs(box[K][1])) for K in radial])
    rr = [r for r in rhos if r < rho0] + [rho0]
    lo, hi = box[T]
    nb = max(1, math.ceil((hi - lo) / beta_step - 1e-9))
    bs = [(lo + (hi - lo) * i / nb, lo + (hi - lo) * (i + 1) / nb) for i in range(nb)]
    secs = []
    for q in L.sectors_for(box, radial):
        k = q[0]
        while k < q[1] - 1e-9:
            secs.append((k, min(q[1], k + sector0))); k += sector0
    cells = np.array([[b[0], b[1], s[0], s[1], rr[i], rr[i + 1]] for b in bs for s in secs for i in range(len(rr) - 1)])
    # drop cells that miss the box
    keep = []
    for c in cells:
        sbox = list(box); sbox[T] = (c[0], c[1])
        keep.append(L.polar_cell_box(sbox, radial, (c[2], c[3]), (c[4], c[5])) is not None)
    cells = cells[np.array(keep)]
    st = dict(checks=0, vec=0, scalar=0, cell=0, cells=len(cells))
    done = np.zeros(len(cells), bool)
    bmid = (lo + hi) / 2
    for s in secs:
        sel = np.where((cells[:, 2] == s[0]) & ~done)[0]
        if not len(sel): continue
        pm = math.radians((s[0] + s[1]) / 2)
        for rp in (0.09, 0.45 * rho0, 0.8 * rho0):
            pt = [0.0, 0.0, 0.0]; pt[T] = bmid
            pt[radial[0]] = rp * math.cos(pm); pt[radial[1]] = rp * math.sin(pm)
            rule = L.probe_rule(pt, typ)
            if rule is None or rule['margin'] <= 0: continue
            rule['rho_probe'] = rp * RAD
            for variant in ('cancel', 'const', 'fixed'):
                sel = np.where((cells[:, 2] == s[0]) & ~done)[0]
                if not len(sel): break
                ok = check_cells(box, typ, radial, T, cells[sel], rule, variant)
                done[sel[ok]] = True
    st['vec'] = int(done.sum())
    rest = cells[~done]
    for rnd in range(12):                                  # batched refinement of the leftovers
        if not len(rest) or len(rest) > 20000: break
        nxt = []
        for c in rest:
            bl, bh, pl, ph_, ra, rb = c
            pm_, bm = (pl + ph_) / 2, (bl + bh) / 2
            if ra == 0.0:
                rsplit = [(0.0, rb / 4), (rb / 4, rb)]
            else:
                rsplit = [(ra, (ra + rb) / 2), ((ra + rb) / 2, rb)]
            kids = np.array([[b0, b1, p0, p1, r0, r1] for (b0, b1) in ((bl, bm), (bm, bh))
                             for (p0, p1) in ((pl, pm_), (pm_, ph_)) for (r0, r1) in rsplit])
            kd = np.zeros(len(kids), bool)
            for fr in (0.3, 0.1, 0.7):
                rp = ra + fr * (rb - ra)
                pt = [0.0, 0.0, 0.0]; pt[T] = bm
                pt[radial[0]] = rp * math.cos(math.radians(pm_)); pt[radial[1]] = rp * math.sin(math.radians(pm_))
                rule = L.probe_rule(pt, typ)
                if rule is None or rule['margin'] <= 0 or rp <= 0: continue
                rule['rho_probe'] = rp * RAD
                for variant in ('cancel', 'const', 'fixed'):
                    sel = np.where(~kd)[0]
                    if not len(sel): break
                    kd[sel[check_cells(box, typ, radial, T, kids[sel], rule, variant)]] = True
                if kd.all(): break
            st['vec'] += int(kd.sum())
            nxt += list(kids[~kd])
        rest = np.array(nxt) if nxt else np.zeros((0, 6))
        st['round%d' % rnd] = len(rest)
    ok, info = _scalar_cells(box, typ, radial, T, rest, st, max_checks)
    return ok, info


def certify_origin_vec(box, typ, eps, max_checks):
    """The all-upright corner.  If the LP certificate at the origin sits on one pair (P, Q) (side by side, tight on
    the line t_P = t_Q = 0), use rays from that line, with the third tilt's range cut into pieces shrinking
    geometrically towards 0; otherwise (stacks, tight only at the origin) rays from the origin."""
    import r3local as L
    lam, m = B.margin_lp([(0.0, 0.0)] * 3, typ)
    sup = [p for p in range(3) if lam is not None and lam[p] > 1e-7]
    order = []
    if len(sup) == 1:
        i, j = B.PAIRS[sup[0]]
        order.append([K for K in range(3) if K not in (i, j)][0])
    order += [T for T in range(3) if T not in order]
    tried = []
    for T in order[:1] if len(sup) == 1 else []:
        lo, hi = box[T]
        cuts = sorted({lo, hi, 0.0} | {sg * 10.0 ** -e for e in range(1, 3) for sg in (1, -1)})
        pieces = [(a, b) for a, b in zip(cuts, cuts[1:]) if a >= lo - 1e-12 and b <= hi + 1e-12 and b > a]
        allok = True
        for pc in pieces:
            bx = list(box); bx[T] = pc
            ok, info = certify_box_type_vec(bx, typ, eps, force_T=T, max_checks=max_checks)
            if not ok:
                allok = False; tried.append((T, pc, info.get('why'))); break
        if allok:
            return True, dict(method=f'line T={T}', pieces=len(pieces))
    r0 = math.sqrt(sum(max(abs(lo), abs(hi)) ** 2 for lo, hi in box))
    ok, info = L.certify_origin_type(typ, r0, max_checks=max_checks)
    if ok: return True, dict(info, method='origin')
    return False, dict(why=f'origin box: lines {tried}, origin {info.get("why")}')
