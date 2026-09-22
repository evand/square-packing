#!/usr/bin/env python3
"""t3_exist_w: Lemma W without the uniform-tilt hypothesis, and the scope audit of the
existence clause (E3).

Task tasks/t3-existence/README.md.  Nothing in search/ is modified; t3_chain is imported.

Lemma W' (proved here, T-uniform, arbitrary tilts).  Let w >= 0 be a normalised Farkas
certificate (sum_r w_r a_r = 0, sum_r w_r = 1) for the wall-carrying LP of
search/BANDCUT_K.md sec 1.1.  Write

    omega_x = w(lo-x rows) = w(hi-x rows),   omega_y = w(lo-y) = w(hi-y),
    Omega   = 2 (omega_x + omega_y)                      ("wall fraction"),
    B       = sum_{wall rows} w_r (u_{i(r)} - 1)/2  +  sum_{pair rows} w_r (W(D_ij) - 1)/2 >= 0
                                                          ("tilt bonus"),
    u_i = |cos th_i| + |sin th_i| in [1, sqrt2],  W(D) = |cos D| + |sin D| in [1, sqrt2].

Then   delta  <=  Omega (T+1)/2  -  1  -  B ,  and in particular

    delta <= 0   whenever   Omega  <=  2 (1 + B) / (T + 1).

usage:
    python3 search/t3_exist_w.py verify [scan.jsonl]      # identity on every dual of the scan
    python3 search/t3_exist_w.py scope  [scan.jsonl]      # which samples are in (E3)'s scope
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t3_chain as TC                                             # noqa: E402


def decompose(w, tags, TH, T):
    """-> dict with omega_x, omega_y, Omega, B, predicted bound, and the direct -sum w b."""
    ox = oy = hx = hy = 0.0
    B = 0.0
    tot = 0.0
    sumwb = 0.0
    u = np.abs(np.cos(TH)) + np.abs(np.sin(TH))
    for q, tg in enumerate(tags):
        if w[q] <= 0:
            continue
        tot += w[q]
        if tg[0] in ('lo', 'hi'):
            i = tg[2]
            P = u[i] / 2.0
            b = P if tg[0] == 'lo' else P - T
            B += w[q] * (u[i] - 1.0) / 2.0
            if tg[1] == 'x':
                if tg[0] == 'lo':
                    ox += w[q]
                else:
                    hx += w[q]
            else:
                if tg[0] == 'lo':
                    oy += w[q]
                else:
                    hy += w[q]
        else:
            _, i, j, o, kind, sg = tg
            D = TH[j] - TH[i]
            W = abs(math.cos(D)) + abs(math.sin(D))
            b = 0.5 + 0.5 * W
            B += w[q] * (W - 1.0) / 2.0
        sumwb += w[q] * b
    Om = ox + hx + oy + hy
    return dict(tot=tot, omega_x=ox, omega_x_hi=hx, omega_y=oy, omega_y_hi=hy,
                Omega=Om, B=B, pred=Om * (T + 1) / 2.0 - 1.0 - B, direct=-sumwb)


def cmd_verify(path):
    rows = [json.loads(l) for l in open(path)]
    worst_sym = worst_id = 0.0
    n, T = 6, 3
    print(f"{'fam':>12} {'delta':>12} {'Omega':>9} {'B':>9} {'W-pred':>12} "
          f"{'-sum w b':>12} {'|diff|':>9} {'|ox-hx|':>9}")
    shown = 0
    for r in rows:
        z = np.array(r['z'])
        A, b, tags = TC.all_rows(z, n, T, r['delta'])
        v, w = TC.dual_bound(A, b, range(len(b)), want_w=True)
        TH = np.array(z[3::3])
        d = decompose(w, tags, TH, T)
        diff = abs(d['pred'] - d['direct'])
        sym = max(abs(d['omega_x'] - d['omega_x_hi']), abs(d['omega_y'] - d['omega_y_hi']))
        worst_sym = max(worst_sym, sym)
        worst_id = max(worst_id, diff)
        if shown < 25 or diff > 1e-9:
            print(f"{r['fam']:>12} {r['delta']:12.6f} {d['Omega']:9.5f} {d['B']:9.5f} "
                  f"{d['pred']:12.8f} {d['direct']:12.8f} {diff:9.2e} {sym:9.2e}")
            shown += 1
    print(f"\nworst |W'-identity - direct| over {len(rows)} duals : {worst_id:.3e}")
    print(f"worst |w(lo-ax) - w(hi-ax)|      over {len(rows)} duals : {worst_sym:.3e}")


def cmd_scope(path):
    rows = [json.loads(l) for l in open(path)]
    inscope = [r for r in rows if r['delta'] > -1e-9]
    hfail = [r for r in rows if r['H'] > 1e-9]
    print(f"samples                      : {len(rows)}")
    print(f"in (E3)'s scope (delta >= 0) : {len(inscope)}")
    print(f"H > 0 anywhere               : {len(hfail)}")
    print(f"H > 0 AND delta >= 0         : {sum(1 for r in hfail if r['delta'] > -1e-9)}")
    print(f"max delta over H-failures    : {max(r['delta'] for r in hfail):.4e}")
    print(f"max delta over all samples   : {max(r['delta'] for r in rows):.4e}")
    fams = {}
    for r in inscope:
        fams.setdefault(r['fam'], []).append(r)
    print("\nin-scope samples by family (all have delta = 0):")
    for f in sorted(fams):
        rs = fams[f]
        th = max(max(min(abs(t % (math.pi / 2)), math.pi / 2 - abs(t % (math.pi / 2)))
                     for t in r['theta']) for r in rs)
        print(f"  {f:>6}  n={len(rs):3d}  max tilt over the family = {math.degrees(th):7.3f} deg"
              f"  H<=0 at {sum(1 for r in rs if r['H'] <= 1e-9)}/{len(rs)}"
              f"  CW<=0 at {sum(1 for r in rs if r['CW'] <= 1e-9)}/{len(rs)}")


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'verify'
    p = sys.argv[2] if len(sys.argv) > 2 else 'runs/t3_chain_scan1.jsonl'
    {'verify': cmd_verify, 'scope': cmd_scope}[cmd](p)
