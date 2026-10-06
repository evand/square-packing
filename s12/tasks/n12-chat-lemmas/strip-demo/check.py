"""Exact checker for a wall-strip certificate (all arithmetic in fractions.Fraction).

    python3 check.py certs/h113.json            # verify
    python3 check.py certs/h113.json --brute    # also re-check the chain condition by brute force

What is verified (PROOF.md section 4 numbers the items):
  (K1) parameters: a <= sqrt2 - 1/2, every menu height z satisfies h - 1/2 < z < 1, h < 3/2;
  (K2) coverage: the boxes cover every pose (t, c) with t in [0,1], p(t) <= c <= h;
  (K3) every leaf claim  "V_box(z_l, z_r) >= q"  holds: q <= the exact rational lower bound of
       leaf.Box.V over the whole box (interval / closed-form bounds, no sampling);
       and for boxes touching theta = 0 or 90 deg, q <= 1 (the value at the axis-parallel pose);
  (K4) the H-flags: a box that may contain a pose with c > a + D(theta) is flagged (recomputed);
  (K5) chain condition: for every sequence of 4 boxes containing a flagged box, the sum of the
       claimed node values along the rule's separation heights is >= 4 (exact min-plus DP; with
       --brute also by enumerating all m^4 sequences).
Exit status 0 iff everything passes.
"""
import json, sys, itertools
from fractions import Fraction as Fr
from leaf import Box, self_test


class CertError(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise CertError(msg)


def load(path_or_dict):
    c = json.load(open(path_or_dict)) if isinstance(path_or_dict, str) else path_or_dict
    cert = dict(h=Fr(c['h']), a=Fr(c['a']), menu=[Fr(z) for z in c['menu']],
                boxes=[tuple(Fr(x) for x in b) for b in c['boxes']], rule=c['rule'],
                leaves=[[(l, r, Fr(q)) for l, r, q in L] for L in c['leaves']])
    return cert


def check_params(C):
    h, a, M = C['h'], C['a'], C['menu']
    need(0 < h < Fr(3, 2), "h out of range")
    need(a + Fr(1, 2) > 0 and (a + Fr(1, 2)) ** 2 <= 2, "a > sqrt2 - 1/2: the one-line case fails")
    need(h - Fr(1, 2) < a < 1, "the ordering line y = a must lie strictly in (h - 1/2, 1)")
    need(len(set(M)) == len(M), "menu has repeats")
    for z in M:
        need(h - Fr(1, 2) < z < 1, f"menu height {z} not strictly inside (h-1/2, 1)")
    return True


def p_lb(t0, t1):
    """Exact lower bound of p = (C+S)/2 over t in [t0, t1] (C+S is unimodal: min at an endpoint)."""
    return Box(t0, t1, Fr(0), Fr(1)).plb


def check_coverage(C):
    """Every feasible pose (p(t) <= c <= h) lies in some box."""
    h, boxes = C['h'], C['boxes']
    for b in boxes:
        need(0 <= b[0] < b[1] <= 1 and b[2] < b[3], f"bad box {b}")
    ts = sorted(set([Fr(0), Fr(1)] + [b[0] for b in boxes] + [b[1] for b in boxes]))
    for u0, u1 in zip(ts, ts[1:]):
        ivs = sorted((b[2], b[3]) for b in boxes if b[0] <= u0 and u1 <= b[1])
        lo = p_lb(u0, u1)          # every pose with t in [u0,u1] has c >= p >= lo
        reach = lo
        for c0, c1 in ivs:
            if c0 <= reach:
                reach = max(reach, c1)
        need(reach >= h, f"coverage gap for t in [{u0},{u1}]: covered only up to c = {reach} < h")
    return len(ts) - 1


def check_leaves(C):
    """Each claimed leaf value is <= the exact lower bound over the box. Returns the claim tables."""
    M, boxes = C['menu'], C['boxes']
    Zs = lambda i: None if i == -1 else M[i]
    tables, nleaf = [], 0
    for bi, (key, claims) in enumerate(zip(boxes, C['leaves'])):
        eb = Box(*key, exact=True)
        tab = {}
        for l, r, q in claims:
            need(not (l == -1 and r == -1), "END-END leaf")
            v = eb.V(Zs(l), Zs(r))
            need(v is not None and q <= v,
                 f"leaf FAILS: box {bi} {key} V({Zs(l)},{Zs(r)}) claimed {q} > proven bound {v}")
            if key[0] == 0 or key[1] == 1:
                need(q <= 1, f"leaf FAILS at the axis-parallel endpoint pose: box {bi}, claim {q} > 1")
            tab[(l, r)] = q
            nleaf += 1
        tables.append(tab)
    return tables, nleaf


def H_flags(C):
    return [Box(*k, exact=True).may_be_H(C['a']) for k in C['boxes']]


def chain_min(C, tables, flags):
    """Exact min over 4-sequences (with >= 1 flagged box) of the claimed chain sum.
    State after position k: (box, left height index) -> (min partial sum, per flag)."""
    m, rule = len(C['boxes']), C['rule']
    NEGINF = None

    def val(i, l, r):
        return tables[i].get((l, r), NEGINF)      # missing claim = no bound = -infinity

    # f[(b, zl, fl)] = min sum of node values of positions < current, current box b with left zl
    f = {}
    for b in range(m):
        f[(b, -1, int(flags[b]))] = Fr(0)
    for pos in range(3):
        g = {}
        for (b, zl, fl), s in f.items():
            for b2 in range(m):
                z = rule[b][b2]
                v = val(b, zl, z)
                if v is NEGINF:
                    return NEGINF, ("missing leaf", b, zl, z)
                key = (b2, z, fl | int(flags[b2]))
                if key not in g or s + v < g[key][0]:
                    g[key] = (s + v, None)
        f = {k: v[0] for k, v in g.items()}
    best, arg = None, None
    for (b, zl, fl), s in f.items():
        if not fl:
            continue
        v = val(b, zl, -1)
        if v is NEGINF:
            return NEGINF, ("missing end leaf", b, zl)
        if best is None or s + v < best:
            best, arg = s + v, (b, zl)
    return best, arg


def chain_brute(C, tables, flags):
    m, rule = len(C['boxes']), C['rule']
    best = None
    for seq in itertools.product(range(m), repeat=4):
        if not any(flags[b] for b in seq):
            continue
        zs = [-1] + [rule[seq[k]][seq[k + 1]] for k in range(3)] + [-1]
        s = Fr(0)
        for k, b in enumerate(seq):
            v = tables[b].get((zs[k], zs[k + 1]))
            if v is None:
                return None
            s += v
        if best is None or s < best:
            best = s
    return best


def check(path_or_dict, brute=False, quiet=False):
    self_test()
    C = load(path_or_dict)
    check_params(C)
    nstrip = check_coverage(C)
    tables, nleaf = check_leaves(C)
    flags = H_flags(C)
    need(any(flags), "no box may contain an H pose: nothing to prove beyond the one-line case")
    mn, arg = chain_min(C, tables, flags)
    need(mn is not None, f"chain condition: missing leaf {arg}")
    need(mn >= 4, f"chain condition FAILS: min chain value {mn} = {float(mn):.6f} < 4")
    out = dict(h=C['h'], boxes=len(C['boxes']), H_boxes=sum(flags), leaf_inequalities=nleaf,
               menu=C['menu'], min_chain=mn, strips=nstrip)
    if brute:
        mb = chain_brute(C, tables, flags)
        need(mb == mn, f"brute force disagrees: {mb} vs {mn}")
        out['brute'] = 'agrees'
    if not quiet:
        print(f"PASS  h = {C['h']} = {float(C['h'])}:  {out['boxes']} boxes ({out['H_boxes']} may hold an H pose), "
              f"{nleaf} leaf inequalities, menu heights {[str(z) for z in C['menu']]}")
        print(f"      min chain value = {float(mn):.7f} >= 4 (exact: {mn.numerator}/{mn.denominator})"
              + ("; brute force agrees" if brute else ""))
    return out


if __name__ == '__main__':
    path = sys.argv[1]
    try:
        check(path, brute='--brute' in sys.argv)
    except CertError as e:
        print("REJECT:", e)
        sys.exit(1)
