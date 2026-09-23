#!/usr/bin/env python3
"""t3_chain_repair: what is the SMALLEST repair that turns a failing H into a certificate?

For every sampled optimum at which the H-restricted dual is > 0, this adds back one extra class
of rows at a time and reports which one first makes the bound <= 0:

  +1rung   the H plus the main-direction (ax-type) pair rows of ONE further pair, minimised over
           the 15 pairs.  A single "rung" across the H.  The support then has cyclomatic number 1.
  +allrung the H plus every ax-type pair row (i.e. only the "one chain" clause is dropped).
  full     everything (= delta*).
"""
import json
import math
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s6skel                                                    # noqa: E402
import t3_chain as TC                                            # noqa: E402


def repair(z, n, T, lvl):
    A, b, tags = TC.all_rows(z, n, T, lvl)
    TH = np.array(z[3::3])
    chs = TC.chains_of(z, n, lvl, int(round(T)))
    best = dict(H=math.inf, rung=math.inf, allr=math.inf, rungpair=None, path=None, axis=None)
    for (ax, path) in chs:
        keep = TC.h_support(tags, TH, ax, path, 'H')
        v, _ = TC.dual_bound(A, b, keep)
        if v < best['H']:
            best['H'], best['path'], best['axis'] = v, list(path), ax
        # every ax-type pair row
        extra = {q for q, tg in enumerate(tags)
                 if tg[0] == 'pair' and TC.row_axis(tg, TH) == ax}
        v, _ = TC.dual_bound(A, b, keep | extra)
        best['allr'] = min(best['allr'], v)
        # one pair at a time
        byp = {}
        for q in extra:
            byp.setdefault((tags[q][1], tags[q][2]), set()).add(q)
        for p, qs in byp.items():
            v, _ = TC.dual_bound(A, b, keep | qs)
            if v < best['rung']:
                best['rung'], best['rungpair'] = v, p
    return best


def main():
    recs = []
    for f in sys.argv[1:]:
        recs += [json.loads(line) for line in open(f)]
    bad = [r for r in recs if r['H'] > 1e-9]
    print(f'# {len(bad)} of {len(recs)} sampled optima at which the H-restricted dual is > 0.')
    print(f'# {"family":11s} {"delta*":>12s} {"H":>12s} {"H+1rung":>12s} {"H+allrungs":>12s} '
          f'{"rung pair":>10s}  max tilt')
    cnt = Counter()
    for r in sorted(bad, key=lambda r: r['fam']):
        z = np.array(r['z'])
        d = repair(z, 6, 3.0, r['delta'])
        tm = math.degrees(max(abs(s6skel.norm_tilt(t)) for t in r['theta']))
        cnt['n'] += 1
        cnt['rung_ok'] += (d['rung'] <= 1e-9)
        cnt['allr_ok'] += (d['allr'] <= 1e-9)
        print(f"  {r['fam']:11s} {r['delta']:+12.4e} {d['H']:+12.4e} {d['rung']:+12.4e} "
              f"{d['allr']:+12.4e} {str(d['rungpair']):>10s}  {tm:6.2f}")
    print(f"\n# one extra rung suffices at {cnt['rung_ok']}/{cnt['n']}; "
          f"all main-direction rows suffice at {cnt['allr_ok']}/{cnt['n']}")


if __name__ == '__main__':
    main()
