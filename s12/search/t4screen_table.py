#!/usr/bin/env python3
"""Collect the `runs/t4_*.log` histories of `search/t4screen.py` into the table of
`search/T4SCREEN.md`.

The logs are APPENDED across the chunks of `search/t4screen_loop.sh` (the JSON is rewritten per
chunk, so it only ever holds the last chunk), and each chunk restarts its stage counter at 0.  So
the trajectory of a tag is the sequence of `STAGE` lines in its log, renumbered here as
`step = 0, 1, 2, ...` per (tag, pattern).

    python3 search/t4screen_table.py runs/t4_*.log
    python3 search/t4screen_table.py --inner runs/t4_L01010101.log     # the inner-loop trace too
"""
import re, sys, os

STAGE = re.compile(
    r'STAGE (\d+) (\S+): LP=([-\d.]+) M=([-\d.]+) kmax=([-\d.]+) L=([-\d.]+) '
    r'cols=(\d+) rows=(\d+) cq=(\d+)(?: \| matched pure LP=([-\d.]+) \(M=([-\d.]+)\) gain=([-+\d.]+))?')
HEAD = re.compile(r"# t4screen t=(\S+) r=\S+ corners=(\S+) patterns=(\S+) cliques=(\S+)")


def rows_of(path):
    t = corners = None
    cq = False
    out = []
    for line in open(path):
        h = HEAD.search(line)
        if h:
            t, corners, cq = h.group(1), h.group(2), h.group(4) == 'True'
            continue
        m = STAGE.search(line)
        if m:
            out.append(dict(t=t, corners=corners, cq=cq, pattern=m.group(2),
                            LP=float(m.group(3)), M=float(m.group(4)), kmax=float(m.group(5)),
                            cols=int(m.group(7)), rows=int(m.group(8)), cliq=int(m.group(9)),
                            LP_pure=float(m.group(10)) if m.group(10) else None,
                            M_pure=float(m.group(11)) if m.group(11) else None,
                            gain=float(m.group(12)) if m.group(12) else None))
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    print(f'{"tag":11} {"t":5} {"corn":5} {"pattern":9} {"cq":3} {"step":4} '
          f'{"LP":>10} {"M":>9} {"kmax":>8} {"LP_pure":>9} {"gain":>8} '
          f'{"poses":>6} {"rows":>6} {"cliq":>5}  converged')
    for p in sorted(args):
        tag = os.path.basename(p)[3:-4]
        step = {}
        for r in rows_of(p):
            k = r['pattern']
            step[k] = step.get(k, -1) + 1
            conv = ('rows ' if r['M'] <= 1 + 1e-9 else '  .  ') + \
                   ('cliq' if r['kmax'] <= 1 + 1e-6 else ' .  ')
            lp_pure = f'{r["LP_pure"]:9.6f}' if r['LP_pure'] is not None else '        -'
            gain = f'{r["gain"]:+8.6f}' if r['gain'] is not None else '       -'
            print(f'{tag:11} {r["t"]:5} {r["corners"]:5} {r["pattern"]:9} {int(r["cq"]):<3} {step[k]:<4} '
                  f'{r["LP"]:10.6f} {r["M"]:9.6f} {r["kmax"]:8.5f} {lp_pure} {gain} '
                  f'{r["cols"]:6} {r["rows"]:6} {r["cliq"]:5}  {conv}')


if __name__ == '__main__':
    main()
