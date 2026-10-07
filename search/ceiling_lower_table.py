#!/usr/bin/env python3
"""Best rigorous lower bound L(s) per s from runs/lower_*.json (post-hoc step of search/nu_f.py)."""
import glob, json
from fractions import Fraction
best = {}
for f in sorted(glob.glob('runs/lower_*.json')):
    d = json.load(open(f)); s = Fraction(d['s']); b = d['best']
    mass = max(h['mass'] for h in d['hist']) if d['hist'] else None
    if s not in best or b['L'] > best[s][0]: best[s] = (b['L'], b.get('M'), mass, d['coverLP'], f)
print("| s | L(s) rigorous | M | grid packing mass (heuristic) | cover-LP over the cert's points (heuristic) | source |")
print("|---|---|---|---|---|---|")
for s in sorted(best):
    L, M, mass, cov, f = best[s]
    print(f"| {float(s):.6f} | {L:.4f} | {M:.4f} | {mass:.4f} | {cov:.4f} | {f} |")
