#!/usr/bin/env python3
"""exactsolve, and while the result is not a local minimum because a corner-corner touch can slip (MILP descent, or no
equilibrium with smooth contacts), kick + slp2 (descend.py) and solve again.

  pipeline.py in.txt [--out DIR] [--cycles 5]      -> DIR/NAME.cK.{log,json,exact.txt,cert} per cycle K
"""
import os, sys, argparse, contextlib, io
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exactsolve as E
import descend as D
import slp2
from slp import save

ap = argparse.ArgumentParser()
ap.add_argument('input')
ap.add_argument('--out', default=os.path.join(HERE, 'descended'))
ap.add_argument('--cycles', type=int, default=5)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
name = os.path.splitext(os.path.basename(a.input))[0]
cur = a.input
for k in range(a.cycles):
    tag = f'{name}.c{k}'
    src = os.path.join(a.out, tag + '.txt')
    if cur != src:
        with open(cur) as f, open(src, 'w') as g:
            g.write(f.read())
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rep, P, kk = E.run(src, outdir=a.out)
    open(os.path.join(a.out, tag + '.log'), 'w').write(buf.getvalue())
    jammed = rep.get('S_exact') is not None and rep.get('vv_milp_dS') is not None and rep['vv_milp_dS'] > -1e-9
    print(f'cycle {k}: input S = {rep["S_input"]}, exact S = {(rep.get("S_exact") or "-")[:25]}, '
          f'MILP dS* = {rep.get("vv_milp_dS")}, {"jammed in every branch: stop" if jammed else "descend"}', flush=True)
    if jammed:
        break
    start = os.path.join(a.out, tag + '.exact.txt') if rep.get('S_exact') else src
    s, sq = slp2.load(start)
    s2, sq2 = D.descend(s, sq, tol=1e-10 if rep.get('S_exact') else 1e-6, log=lambda *x: print('   ', *x, flush=True))
    if s2 >= s - 1e-13:
        print('    no decrease; stop')
        break
    cur = os.path.join(a.out, f'{name}.c{k + 1}.in.txt')
    save(cur, s2, sq2)
