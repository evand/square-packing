#!/usr/bin/env python3
"""t4_cycles45: the uniform-tilt line at T = 4, n = 12, and the 45 deg certificate exactly.

Task tasks/t4-cycles/README.md item 3.  Imports t3_chain / t4_cycles / s6skel; modifies nothing.

Three questions.
  (a) What is delta*(T = 4, n = 12) at a common tilt t, and what shape is its dual?
  (b) At 45 deg: is the dual a cycle, and does the closed form of notes/t3-chain.md sec 5.3
      reproduce it?  The sec 3.4 GUESS there is that a T = 4 cycle needs k >= 2 sqrt2 (4 - sqrt2)
      = 7.31 links, i.e. a cycle of eight, because it assumes q = 4 turns.  The general condition
      is  k >= q sin(alpha/2) (T - u) , and q is what the configuration chooses.
  (c) Does a cycle of eight exist at 45 deg among twelve squares, and what does it bound?
"""
import argparse
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s6skel                                                     # noqa: E402
import t3_chain as TC                                             # noqa: E402
import t3_chain_shape as TS                                       # noqa: E402
import t4_cycles as TZ                                            # noqa: E402


def dual_at(z, n, T, lvl, wtol=1e-9, tol=1e-7):
    A, b, tags = TC.all_rows(z, n, T, lvl, tol=tol)
    fv, w = TC.dual_bound(A, b, range(len(b)), want_w=True)
    TH = np.array(z[3::3])
    supp = [(float(w[q]), tags[q], TC.row_axis(tags[q], TH)) for q in range(len(b))
            if w[q] > wtol]
    supp.sort(key=lambda s: -s[0])
    return A, b, tags, fv, w, supp


def cycle_from_support(supp, wtol=1e-9):
    """if the PAIR rows of a dual support form exactly one simple cycle, return it."""
    edges = [(min(tg[1], tg[2]), max(tg[1], tg[2])) for (w, tg, ax) in supp
             if tg[0] == 'pair' and w > wtol]
    if len(edges) != len(set(edges)):
        return None
    deg = defaultdict(int)
    for (i, j) in edges:
        deg[i] += 1; deg[j] += 1
    if not edges or any(d != 2 for d in deg.values()) or len(edges) != len(deg):
        return None
    cy = TZ.simple_cycles(edges, kmax=len(edges))
    return cy[0] if len(cy) == 1 and len(cy[0]) == len(edges) else None


def sec53(z, T, cyc, pairrows, tags, b, A):
    """the sec 5.3 closed form on a cycle, with its hypotheses checked (t4_cycles)."""
    return TZ.cycle_closed_form(z, len(z) // 3, T, cyc, pairrows, tags, b, A=A)


def cmd_uniform(a):
    rng = np.random.default_rng(a.seed)
    print(f'# uniform tilt at T = {a.T}, n = {a.n}.  delta* by multistart '
          f'(limit {a.limit} x 3 + {a.extra} random starts) -- a LOWER bound.')
    print(f'# {"t(deg)":>8s} {"delta*":>14s} {"fulldual":>14s} {"mu":>3s} {"E":>3s} {"V":>3s} '
          f'{"walls":>18s} {"cls":>10s} {"H":>13s} {"cyc(LP)":>13s} {"cyc len":>8s}')
    out = []
    for td in a.tilts:
        th = [math.radians(td)] * a.n
        v, z, _ = TC.delta_star(th, a.T, rng, limit=a.limit, extra=a.extra)
        A, b, tags, fv, w, supp = dual_at(z, a.n, a.T, v)
        g = TS.support_graph(supp, wtol=1e-9)
        gg = dict(g, walls={f'{k[0]}-{k[1]}': float(x) for k, x in g['walls'].items()})
        tight, _ = TZ.tight_indices(A, b, z, v)
        pg = TZ.tight_pair_graph(tags, tight)
        cv, cyc, cw, ncyc = TZ.best_cycle_bound(A, b, tags, pg, kmax=a.kmax)
        r = TC.best_h(z, a.n, a.T, v, modes=('H',), rows=(A, b, tags))
        print(f'  {td:8.4f} {v:+14.9f} {fv:+14.9f} {g["mu"]:3d} {g["nedge"]:3d} '
              f'{g["nvert"]:3d} {",".join(sorted(gg["walls"])):>18s} '
              f'{TZ.classify(gg):>10s} {r["H"]["bound"]:+13.6e} {cv:+13.6e} '
              f'{len(cyc) if cyc else 0:8d}')
        out.append(dict(tdeg=td, delta=float(v), fulldual=float(fv), mu=g['mu'],
                        E=g['nedge'], V=g['nvert'], walls=gg['walls'],
                        cls=TZ.classify(gg), H=float(r['H']['bound']),
                        cycLP=float(cv), cyc=list(cyc) if cyc else None, ncycles=ncyc,
                        z=[float(q) for q in z],
                        supp=[[s[0], list(map(TC._js, s[1])), s[2]] for s in supp]))
    if a.out:
        json.dump(out, open(a.out, 'w'), indent=1)


def cmd_deep45(a):
    """the 45 deg point in full: the dual, the sec 5.3 arithmetic, and cycles by length."""
    rng = np.random.default_rng(a.seed)
    th = [math.radians(a.tdeg)] * a.n
    v, z, tops = TC.delta_star(th, a.T, rng, limit=a.limit, extra=a.extra, keep=a.keep)
    print(f'# t = {a.tdeg} deg, T = {a.T}, n = {a.n}')
    print(f'# delta* (multistart, LOWER bound) = {v!r}')
    u = math.cos(math.radians(a.tdeg)) + math.sin(math.radians(a.tdeg))
    for zz in tops:
        vv = float(zz[0])
        A, b, tags, fv, w, supp = dual_at(zz, a.n, a.T, vv)
        g = TS.support_graph(supp, wtol=1e-9)
        gg = dict(g, walls={f'{k[0]}-{k[1]}': float(x) for k, x in g['walls'].items()})
        print(f'\n## local optimum delta = {vv:+.15f}, dual = {fv:+.15f}')
        print(f'   support graph: mu={g["mu"]} E={g["nedge"]} V={g["nvert"]} '
              f'walls={gg["walls"]} cls={TZ.classify(gg)}')
        for (ww, tg, ax) in supp:
            print(f'     {ww:.12f}  {tg}  {ax}')
        omega = sum(ww for (ww, tg, ax) in supp if tg[0] == 'lo')
        print(f'   Lemma W: omega = {omega:.12f},  omega(T+2-u)-1 = '
              f'{omega*(a.T+2-u)-1:+.15f}')
        cyc = cycle_from_support(supp)
        tight, _ = TZ.tight_indices(A, b, zz, vv)
        pg = TZ.tight_pair_graph(tags, tight)
        if cyc is not None:
            print(f'   the pair rows ARE one simple cycle: {cyc}')
            cf = TZ.cycle_closed_form(zz, a.n, a.T, cyc, pg, tags, b, A=A)
            print(f'   sec 5.3 closed form: {json.dumps(cf["best"])}')
            print(f'   naive (worst) form : {json.dumps(cf["naive"])}')
        else:
            print('   the pair rows are NOT a single simple cycle')
        # cycles by length
        print(f'   tight pair edges {len(pg)}')
        by = defaultdict(lambda: [math.inf, None, 0])
        cycs = TZ.simple_cycles(list(pg), kmax=a.kmax)
        rng2 = np.random.default_rng(11)
        use = cycs if len(cycs) <= a.maxcyc else [cycs[i] for i in
                                                  rng2.choice(len(cycs), a.maxcyc,
                                                              replace=False)]
        for c in use:
            keep = TZ.cycle_support(tags, c, pg)
            val, _ = TC.dual_bound(A, b, keep)
            by[len(c)][2] += 1
            if val < by[len(c)][0]:
                by[len(c)][0] = val
                by[len(c)][1] = list(c)
        print(f'   {len(cycs)} simple cycles in the tight graph '
              f'({len(use)} evaluated); best restricted dual by cycle length:')
        for L in sorted(by):
            val, c, cnt = by[L]
            cf = None
            if c is not None:
                cf = TZ.cycle_closed_form(zz, a.n, a.T, tuple(c), pg, tags, b, A=A)
            bb = (cf or {}).get('best')
            nn = (cf or {}).get('naive')
            print(f'     k={L:2d}  n={cnt:6d}  best LP {val:+.9f}  '
                  f'sec5.3 conforming {"none" if bb is None else "%+.9f (q=%d)" % (bb["value"], bb["q"])}'
                  f'  naive {"none" if nn is None else "%+.9f (q=%d, axis_ok=%s)" % (nn["value"], nn["q"], nn["axis_ok"])}')
        if a.first_only:
            break


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('uniform')
    p.add_argument('--T', type=float, default=4.0)
    p.add_argument('--n', type=int, default=12)
    p.add_argument('--tilts', type=float, nargs='*',
                   default=[1, 5, 10, 15, 20, 28.0725, 30, 35, 40, 44, 45])
    p.add_argument('--limit', type=int, default=300)
    p.add_argument('--extra', type=int, default=600)
    p.add_argument('--kmax', type=int, default=9)
    p.add_argument('--seed', type=int, default=20260921)
    p.add_argument('--out', default='')
    p.set_defaults(fn=cmd_uniform)
    p = sub.add_parser('deep45')
    p.add_argument('--T', type=float, default=4.0)
    p.add_argument('--n', type=int, default=12)
    p.add_argument('--tdeg', type=float, default=45.0)
    p.add_argument('--limit', type=int, default=500)
    p.add_argument('--extra', type=int, default=1500)
    p.add_argument('--keep', type=int, default=3)
    p.add_argument('--kmax', type=int, default=9)
    p.add_argument('--maxcyc', type=int, default=6000)
    p.add_argument('--first-only', action='store_true')
    p.add_argument('--seed', type=int, default=20260921)
    p.set_defaults(fn=cmd_deep45)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
