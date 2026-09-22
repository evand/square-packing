#!/usr/bin/env python3
"""t3_exist_master: the wall-free form of the Farkas dual, and the corrected cycle lemma.

Task tasks/t3-existence/README.md.  Nothing in search/ is modified; t3_chain is imported.

MASTER FORMULA (proved in notes/t3-existence.md sec 2).  At any configuration z of n unit
squares in [0,T]^2 with margin delta, every wall row is satisfied, so a Farkas certificate is
determined by its PAIR-row weights alone: give each satisfied pair row e (a directed link
tail(e) -> head(e) with unit normal n_e and rhs m_e) a weight w_e >= 0, let

    r_i  =  sum_{e: head(e)=i} w_e n_e  -  sum_{e: tail(e)=i} w_e n_e        (the residual at i)

and cancel r_i with wall rows at i (cost |r_i|_1, the L1 norm -- walls carry +-e_x, +-e_y).  Then

    delta  <=  [ sum_i |r_i|_1 (T - u_i)/2  -  sum_e w_e m_e ]
               / [ sum_e w_e  +  sum_i |r_i|_1 ]  ,       u_i = |cos th_i| + |sin th_i| ,

and this wall completion is optimal, so minimising the right-hand side over w >= 0 reproduces the
full dual exactly.  In particular

    delta <= 0    iff    sum_i |r_i|_1 (T - u_i) / 2  <=  sum_e w_e m_e .

usage:
    python3 search/t3_exist_master.py verify [scan.jsonl] [nmax]   # master formula vs the LP
    python3 search/t3_exist_master.py cycle                        # the corrected cycle lemma
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t3_chain as TC                                             # noqa: E402

N, T = 6, 3


def link_of(tag, TH):
    """(tail, head, n, m) for a pair row tag = ('pair', i, j, o, kind, sg)."""
    _, i, j, o, kind, sg = tag
    c, s = math.cos(TH[o]), math.sin(TH[o])
    n = np.array([c, s]) if kind == 0 else np.array([-s, c])
    n = sg * n
    D = TH[j] - TH[i]
    m = 0.5 + 0.5 * (abs(math.cos(D)) + abs(math.sin(D)))
    # row reads  sg n.(c_j - c_i) >= m + delta  ->  link i -> j with normal sg*n
    return i, j, n, m


def master(wpair, tags, TH, n=N, Tc=T):
    """the master formula for a weighting of the pair rows.  wpair: dict row-index -> weight."""
    u = np.abs(np.cos(TH)) + np.abs(np.sin(TH))
    r = np.zeros((n, 2))
    Wm = 0.0
    Wt = 0.0
    for q, w in wpair.items():
        if w <= 0:
            continue
        i, j, nrm, m = link_of(tags[q], TH)
        r[j] += w * nrm
        r[i] -= w * nrm
        Wm += w * m
        Wt += w
    L1 = np.abs(r).sum(axis=1)
    cost = float((L1 * (Tc - u) / 2.0).sum())
    den = Wt + float(L1.sum())
    return (cost - Wm) / den if den > 0 else math.inf, cost, Wm, float(L1.sum())


def cmd_verify(path, nmax):
    rows = [json.loads(l) for l in open(path)][:nmax]
    worst = 0.0
    print(f"{'fam':>12} {'delta':>12} {'F (LP)':>13} {'master(w_F)':>13} {'diff':>10}")
    for r in rows:
        z = np.array(r['z'])
        A, b, tags = TC.all_rows(z, N, T, r['delta'])
        v, w = TC.dual_bound(A, b, range(len(b)), want_w=True)
        TH = np.array(z[3::3])
        wp = {q: w[q] for q in range(len(b)) if tags[q][0] == 'pair' and w[q] > 0}
        mv, _, _, _ = master(wp, tags, TH)
        d = abs(mv - v)
        worst = max(worst, d)
        if d > 1e-7 or rows.index(r) < 15:
            print(f"{r['fam']:>12} {r['delta']:12.6f} {v:13.9f} {mv:13.9f} {d:10.2e}")
    print(f"\nworst |master(w_pair of the LP dual) - LP value| over {len(rows)} : {worst:.3e}")


def cmd_cycle():
    """the corrected cycle lemma: turn cost is the L1 norm of the normal jump, not 2 sin(a/2)."""
    print("corrected cycle lemma:  a cycle of k links at common tilt t with q turns certifies")
    print("  delta <= 0   iff   k m  >=  q |n_in - n_out|_1 (T - u)/2 .")
    print("At a common tilt t the two edge normals are n0 = (C,S), n1 = (-S,C):")
    print("  |n0 - n1|_1 = (C+S) + |S-C| = 2C  for t <= 45 deg,   vs  2 sin(45 deg) = sqrt2.")
    print(f"\n{'t deg':>7} {'u':>8} {'2C':>8} {'2sin(a/2)':>10} {'k* (L1)':>9} {'k* (L2)':>9}")
    for td in (0, 5, 10, 20, 28.8, 30, 40, 45):
        t = math.radians(td)
        C, S = math.cos(t), math.sin(t)
        u = C + S
        l1 = 2 * C
        l2 = math.sqrt(2.0)
        print(f"{td:7.1f} {u:8.5f} {l1:8.5f} {l2:10.5f} "
              f"{4 * l1 * (T - u) / 2:9.4f} {4 * l2 * (T - u) / 2:9.4f}")
    print("\n(k* = the number of links a q = 4 cycle needs;  L2 is t3-chain.md sec 5.3's formula,")
    print(" valid only when the jump is axis-aligned, i.e. only at t = 45 deg where 2C = sqrt2.)")


def cmd_near():
    """the bare chain of three and the chain+2 legs at small tilts, exactly and to first order."""
    import t3_chain_hform as HF
    print("bare wall-to-wall chain of three at a COMMON tilt t (no legs), T = 3:")
    print(f"{'t deg':>8} {'exact H(3,0,t)':>15} {'t/4':>10} {'master formula':>15}")
    for td in (0.5, 1, 2, 5, 10, 20):
        t = math.radians(td)
        C, S = math.cos(t), math.sin(t)
        u = C + S
        # links with the common x-normal n = (C, S), both at weight 1
        r1 = np.array([-C, -S]); r2 = np.zeros(2); r3 = np.array([C, S])
        L = [abs(r).sum() for r in (r1, r2, r3)]
        cost = (L[0] + L[2]) * (3 - u) / 2.0
        val = (cost - 2.0) / (2.0 + sum(L))
        print(f"{td:8.2f} {HF.H(3, 0, t):15.8f} {t/4:10.6f} {val:15.8f}")
    print("\nchain of three + L legs at a common tilt (Lemma H), by the master formula:")
    print(f"{'t deg':>8} {'L=0':>12} {'L=1':>12} {'L=2':>12} {'L=3':>12}   (= H(3,L,t))")
    for td in (1, 5, 20, 28.8, 45):
        t = math.radians(td)
        row = f"{td:8.2f}"
        for Lg in (0, 1, 2, 3):
            row += f" {HF.H(3, Lg, t):12.6f}"
        print(row)
    print("\nso the bare chain is +t/4 at first order and the legs are worth exactly the tilt")
    print("they cancel: H(3,2,t) <= 0 for t <= 28.7959 deg, and Corollary P3 kills t >= 25.8431.")


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'verify'
    if cmd == 'cycle':
        cmd_cycle()
    elif cmd == 'near':
        cmd_near()
    else:
        p = sys.argv[2] if len(sys.argv) > 2 else 'runs/t3_chain_scan1.jsonl'
        nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 200
        cmd_verify(p, nmax)
