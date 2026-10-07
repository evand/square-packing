#!/usr/bin/env python3
"""EXACT certification of an anchor-clique-feasible packing measure at t = 399/100 (task G).

Takes a float support file (`runs/cq_*_support.txt` from search/clique_continuum.py), snaps every
pose to a rational rotation `theta = 2 arctan(p/q)` and a rational centre clamped exactly into the
closed container, and then

  * enumerates the arrangement vertices of the support exactly (search/dual_exact.py machinery)
    and gets the exact maximum coverage M;
  * separates anchor cliques K(p, A) (search/clique_family.py) on the current measure, re-solving
    the LP with the cut rows, until no violated anchor clique is found;
  * rounds the masses DOWN to multiples of 1/DM, rescales if necessary so that M <= 1 and every
    cut is <= 1 exactly, and re-checks everything in integers.

Output: mass (exact rational), M (exact), the exact value of every cut, and the exact support
file.  Because a measure that satisfies every clique constraint has mass at most the clique-LP
value, an exactly certified measure of mass >= 12 that satisfies every anchor clique is a
statement about how much the anchor family can ever be worth (see search/CLIQUE_CONTINUUM.md).

Usage
    python3 search/clique_exact.py runs/cq_A99sw1_support.txt A99sw1 [--rounds 8] [--procs 2]
"""
import sys, os, math, time, json, argparse
from fractions import Fraction as Fr
from math import gcd, lcm

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, 'runs')
sys.path.insert(0, HERE)
import numpy as np
import dual_exact as DE
import clique_family as CF

T = DE.T                      # 399/100


# ---------------------------------------------------------------- exact predicates
def sq_frame(sqrec):
    """(cx, cy, cos, sin) as Fractions for a dual_exact square record"""
    CX, CY, Dk, a, b, r = sqrec[:6]
    return Fr(CX, Dk), Fr(CY, Dk), Fr(a, r), Fr(b, r)


def pt_in_sq(q, sqrec):
    cx, cy, c, s = sq_frame(sqrec)
    dx = q[0] - cx
    dy = q[1] - cy
    return abs(dx * c + dy * s) <= Fr(1, 2) and abs(-dx * s + dy * c) <= Fr(1, 2)


def seg_meets_sq(a0, a1, sqrec):
    """closed segment [a0,a1] meets the closed unit square (exact SAT, three axes)"""
    cx, cy, c, s = sq_frame(sqrec)
    # the square's two normals
    for (ex, ey) in ((c, s), (-s, c)):
        cen = cx * ex + cy * ey
        v0 = a0[0] * ex + a0[1] * ey
        v1 = a1[0] * ex + a1[1] * ey
        lo, hi = min(v0, v1), max(v0, v1)
        if cen - Fr(1, 2) > hi or cen + Fr(1, 2) < lo:
            return False
    # the segment's normal
    dx = a1[0] - a0[0]
    dy = a1[1] - a0[1]
    ex, ey = -dy, dx
    if ex != 0 or ey != 0:
        h = a0[0] * ex + a0[1] * ey
        cen = cx * ex + cy * ey
        rad = (abs(c * ex + s * ey) + abs(-s * ex + c * ey)) / 2
        if cen - rad > h or cen + rad < h:
            return False
    return True


def cut_members(cut, squares):
    """mask over the images: membership in K(p, A) (exact)"""
    p = cut['p']
    A = cut['A']
    a0, a1 = A[0], A[-1]
    out = []
    for srec in squares:
        if pt_in_sq(p, srec):
            out.append(seg_meets_sq(a0, a1, srec))
        else:
            out.append(all(pt_in_sq(v, srec) for v in A))
    return out


# ---------------------------------------------------------------- main
def read_float_support(path):
    poses = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        w = line.split()
        if w[0] == 'pose':
            poses.append((float(w[1]), float(w[2]), math.radians(float(w[3])), float(w[4])))
    return poses


def images_float(poses, masses):
    out = []
    for k, (p, q, cx, cy) in enumerate(poses):
        th = 2 * math.atan2(p, q)
        for (x, y, pp, qq) in DE.images(p, q, cx, cy):
            thh = 2 * math.atan2(pp, qq)
            out.append((float(x), float(y), thh, masses[k] / 8.0))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('support')
    ap.add_argument('tag')
    ap.add_argument('--Q', type=int, default=100000)
    ap.add_argument('--Dc', type=int, default=1000000)
    ap.add_argument('--DM', type=int, default=10 ** 9)
    ap.add_argument('--procs', type=int, default=2)
    ap.add_argument('--rounds', type=int, default=8)
    ap.add_argument('--sep-pitch', type=float, default=0.02)
    ap.add_argument('--sep-thr', type=float, default=0.95)
    ap.add_argument('--sep-pts', type=int, default=120)
    ap.add_argument('--nseed', type=int, default=6)
    ap.add_argument('--maxstep', type=int, default=20)
    ap.add_argument('--no-cuts', action='store_true', help='point rows only (pure certification)')
    args = ap.parse_args()
    sys.set_int_max_str_digits(0)
    log = open(os.path.join(RUNS, f'cqx_{args.tag}.log'), 'w')

    def say(m):
        print(m, flush=True)
        log.write(m + '\n')
        log.flush()

    src = read_float_support(args.support)
    say(f"# {args.support}: {len(src)} float poses; snap Q={args.Q}, Dc={args.Dc}, t={T}")
    poses = []
    for (cx, cy, th, mu) in src:
        poses.append(DE.snap_pose(cx, cy, th, args.Q, args.Dc))
    squares, pose_of = DE.build_squares(poses)
    say(f"#   {len(squares)} images, all admissible (exact)")
    DE.enumerate_vertices(squares, full=False, procs=args.procs)
    inc = DE.incidences(args.procs)
    pat = DE.patterns(inc, pose_of)
    keys = list(pat.keys())
    nP = len(poses)
    say(f"#   {len(DE.VERTS)} vertices in F, {len(keys)} distinct incidence patterns")

    import scipy.sparse as sp
    from scipy.optimize import linprog
    rows, cols, vals = [], [], []
    for i, key in enumerate(keys):
        for k, c in key:
            rows.append(i)
            cols.append(k)
            vals.append(c / 8.0)
    Apt = sp.csr_matrix((vals, (rows, cols)), shape=(len(keys), nP))

    cuts = []
    Ccols = []
    mu = None
    for rnd in range(args.rounds):
        A = Apt if not Ccols else sp.vstack([Apt, sp.csr_matrix(np.array(Ccols))], format='csr')
        res = linprog(-np.ones(nP), A_ub=A, b_ub=np.ones(A.shape[0]), bounds=(0, None),
                      method='highs')
        assert res.status == 0, res.message
        mu = np.maximum(res.x, 0.0)
        say(f"  round {rnd}: LP mass {-res.fun:.9f} with {len(cuts)} anchor cuts")
        if args.no_cuts:
            break
        IM = images_float(poses, mu)
        IM = IM[IM[:, 3] > 1e-13]
        if len(IM) == 0:
            break
        imgs = [(a, b, c, d) for (a, b, c, d) in IM]
        cand = CF.tight_points(imgs, float(T), args.sep_pitch, args.sep_thr)
        cand.sort(key=lambda z: -z[1])
        if len(cand) > args.sep_pts:
            cand = cand[::max(1, len(cand) // args.sep_pts)][:args.sep_pts]
        found = []
        for (p, cv) in cand:
            r = CF.anchor_separate_all(IM, p, float(T), nseed=args.nseed, maxstep=args.maxstep)
            if r is not None and r['mass'] > 1.0 + 1e-7:
                found.append(r)
        found.sort(key=lambda r: -r['mass'])
        say(f"    separation: {len(found)} violated anchor cliques, worst "
            f"{found[0]['mass']:.6f}" if found else "    separation: none violated")
        if not found:
            break
        for r in found[:20]:
            # snap the anchor to rationals: any (p, A) is a valid clique
            p = (Fr(round(r['p'][0] * args.Dc), args.Dc), Fr(round(r['p'][1] * args.Dc), args.Dc))
            Av = [(Fr(round(v[0] * args.Dc), args.Dc), Fr(round(v[1] * args.Dc), args.Dc))
                  for v in r['A']]
            cut = dict(p=p, A=Av)
            mem = cut_members(cut, squares)
            col = np.zeros(nP)
            for si, m in enumerate(mem):
                if m:
                    col[pose_of[si]] += 1.0 / 8.0
            cuts.append(cut)
            Ccols.append(col)

    # ---------------- final re-solve with every cut in the model, then exact rounding
    A = Apt if not Ccols else sp.vstack([Apt, sp.csr_matrix(np.array(Ccols))], format='csr')
    res = linprog(-np.ones(nP), A_ub=A, b_ub=np.ones(A.shape[0]), bounds=(0, None),
                  method='highs')
    assert res.status == 0, res.message
    mu = np.maximum(res.x, 0.0)
    say(f"  final LP mass {-res.fun:.9f} with all {len(cuts)} anchor cuts")
    DM = args.DM
    MU = [int(math.floor(float(m) * DM)) for m in mu]
    cutmem = [cut_members(c, squares) for c in cuts]

    def exact_cut_values(MU):
        out = []
        for mem in cutmem:
            num = 0
            for si, m in enumerate(mem):
                if m:
                    num += MU[pose_of[si]]
            out.append(Fr(num, 8 * DM))
        return out

    M, key = DE.exact_max(pat, MU, DM)
    cv = exact_cut_values(MU)
    worst = max([M] + cv) if cv else M
    say(f"  rounded down to /{DM}: mass {sum(MU)/DM:.9f}, exact M = {float(M):.12f}, "
        f"max cut = {float(max(cv)) if cv else 0:.12f}")
    if worst > 1:
        MU = [(m * worst.denominator) // worst.numerator for m in MU]
        M, key = DE.exact_max(pat, MU, DM)
        cv = exact_cut_values(MU)
        worst = max([M] + cv) if cv else M
        say(f"  scaled by 1/{float(worst):.12f} and rounded down: mass {sum(MU)/DM:.9f}, "
            f"exact M = {float(M):.12f}, max cut = {float(max(cv)) if cv else 0:.12f}")
    assert worst <= 1, worst
    mass = Fr(sum(MU), DM)
    say(f"# EXACT: mass = {mass} = {float(mass):.12f};  M = {M} = {float(M):.12f};  "
        f"{len(cuts)} anchor cuts, max value {float(max(cv)) if cv else 0:.12f};  "
        f"mass {'>=' if mass >= 12 else '<'} 12")
    out = os.path.join(RUNS, f'cqx_{args.tag}_support.txt')
    DE.write_exact_support(out, poses, MU, DM, M, mass / M if M > 0 else 0,
                           f"anchor-clique-feasible, snapped from {os.path.basename(args.support)}")
    json.dump(dict(tag=args.tag, mass=str(mass), mass_float=float(mass), M=str(M),
                   M_float=float(M), cuts=len(cuts),
                   max_cut=str(max(cv)) if cv else '0',
                   max_cut_float=float(max(cv)) if cv else 0.0,
                   poses=sum(1 for m in MU if m > 0)),
              open(os.path.join(RUNS, f'cqx_{args.tag}.json'), 'w'), indent=1)
    say(f"# wrote {out}")


if __name__ == '__main__':
    main()
