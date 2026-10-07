#!/usr/bin/env python3
"""Collect the nu_f batch (search/ceiling_launch.sh) into the table for search/CEILING.md.

  * L(s): best rigorous packing lower bound from runs/nuf_e*.json (exact-mode runs).
  * U(s'): for every certificate produced (runs/nuf_e*_best.txt exact-mode, runs/cert_c*.txt
    cell-mode), rescale to each s' of the grid, run the exact verifier at N (default 6000), and
    take  U(s') = min over certificates of  total / min_covered(s').

Usage:  python3 search/ceiling_collect.py [--N 6000] [--threads 8] [--par 4]
"""
import sys, os, glob, json, math
from fractions import Fraction
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_f

GRID = [Fraction(3920, 997), Fraction('3.94'), Fraction('3.945'), Fraction('3.95'), Fraction('3.9525'),
        Fraction('3.955'), Fraction('3.9575'), Fraction('3.96'), Fraction('3.965'), Fraction('3.97'),
        Fraction('3.98'), Fraction('3.99'), Fraction(4)]


def eval_cert(args):
    cert, N, th = args
    out = []
    tmp = cert + f'.scaled.{os.getpid()}.tmp'
    for st in GRID:
        try:
            tot, sp_ = nu_f.rescale_cert(cert, float(st), tmp)
            mv = nu_f.run_verifier(tmp, N=N, threads=th)
            U = tot / float(mv) if mv > 0 else float('inf')
            out.append((str(st), float(st), sp_, tot, float(mv), U))
        except Exception as e:
            out.append((str(st), float(st), None, None, None, float('inf')))
    if os.path.exists(tmp): os.remove(tmp)
    return cert, out


if __name__ == '__main__':
    a = sys.argv[1:]; N = 6000; th = 8; par = 4
    for i in range(0, len(a), 2):
        if a[i] == '--N': N = int(a[i + 1])
        if a[i] == '--threads': th = int(a[i + 1])
        if a[i] == '--par': par = int(a[i + 1])
    certs = ['certificates/s12_lower_3.931795.txt'] + sorted(glob.glob('runs/nuf_e*_best.txt')) + sorted(glob.glob('runs/cert_c*.txt'))
    with Pool(par) as p: res = p.map(eval_cert, [(c, N, th) for c in certs])
    U = {}; detail = {}
    for cert, out in res:
        detail[cert] = out
        for key, sf, sp_, tot, mv, u in out:
            if u < U.get(key, (float('inf'),))[0]: U[key] = (u, cert, tot, mv, sp_)
    Lb = {}
    for jf in sorted(glob.glob('runs/nuf_e*.json')):
        d = json.load(open(jf)); Lb[str(Fraction(d['s']))] = (d['best']['L'], d['best']['U'], jf, len(d['hist']), d['hist'][-1]['t'])
    json.dump(dict(N=N, U={k: v for k, v in U.items()}, L=Lb, detail=detail), open('runs/ceiling_table.json', 'w'), indent=1)
    print(f"| s | L(s) = rigorous lower bound on nu_f | U(s) = rigorous upper bound on COVER | from |")
    print("|---|---|---|---|")
    for st in GRID:
        k = str(st); u = U.get(k); l = Lb.get(k)
        ls = f"{l[0]:.4f}" if l else "-"
        us = f"{u[0]:.4f}" if u and u[0] < 1e9 else "-"
        src = os.path.basename(u[1]) if u else "-"
        print(f"| {float(st):.6f} | {ls} | {us} | {src} (min {u[3]:.6f}, total {u[2]:.5f}) |" if u else f"| {float(st):.6f} | {ls} | - | - |")
