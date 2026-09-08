#!/usr/bin/env python3
"""Collect the `runs/t4_*.json` histories of `search/t4screen.py` into the table of
`search/T4SCREEN.md`: one row per (tag, pattern), the value trajectory over refinement stages,
and the convergence flags that say whether a number is a bound at all.

    python3 search/t4screen_table.py runs/t4_*.json
"""
import json, sys, os


def main():
    paths = sys.argv[1:]
    print(f'{"tag":8} {"t":6} {"corners":8} {"pattern":9} {"cq":3} '
          f'{"stage":5} {"LP":>10} {"M":>10} {"kmax":>8} {"LP_pure":>10} {"gain":>9} '
          f'{"poses":>6} {"rows":>6} {"cliq":>5}  conv')
    for p in sorted(paths):
        try:
            js = json.load(open(p))
        except Exception as e:
            print(f'# {p}: {e}'); continue
        a = js['args']; tag = js['tag']; t = js['t']
        last = {}
        for rec in js['hist']:
            last[(rec['pattern'], rec['stage'])] = rec
        for (pat, st) in sorted(last, key=lambda k: (k[0], k[1])):
            r = last[(pat, st)]
            conv = ('rows' if r['M'] <= 1 + 1e-9 else '    ') + '+' + \
                   ('cliq' if r['kmax'] <= 1 + 1e-6 else '    ')
            print(f'{tag:8} {t:<6} {a["corners"]:8} {pat:9} {int(bool(a["cliques"])):<3} '
                  f'{st:<5} {r["LP"]:10.6f} {r["M"]:10.7f} {r["kmax"]:8.5f} '
                  f'{r.get("LP_pure", float("nan")):10.6f} {r.get("gain", float("nan")):9.6f} '
                  f'{r["cols"]:6} {r["rows"]:6} {r["cq"]:5}  {conv}')


if __name__ == '__main__':
    main()
