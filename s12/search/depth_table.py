#!/usr/bin/env python3
"""per-stage table of the depth-test runs (runs/tl_D*.json)"""
import glob, json, os, sys
rows = []
for p in sorted(glob.glob(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', 'tl_D*.json'))):
    h = json.load(open(p))
    tag = h['tag']; a = h['args']
    for r in h['hist']:
        if r.get('it') != 'final':
            continue
        q = r.get('quads')
        rows.append((tag, a.get('quad'), a.get('quad_mode'), r['stage'], r['LP'], r['M'], r['kmax'],
                     r['interior'], q, r['cols'], r['rows'], r['cq'], r.get('LP_pure'), r.get('LP_nochord'),
                     r['bnd'], r['converged'], r['secs']))
print(f"{'tag':9s} {'q':5s} {'mode':4s} st {'LP':>10s} {'M':>11s} {'kmax':>8s} {'int':>8s} {'quads':>30s} {'cols':>5s} {'rows':>6s} {'cq':>4s} {'pure':>10s} {'nochord':>10s} {'bnd':>6s} conv {'secs':>6s}")
for (tag, q, m, st, LP, M, K, I, qs, c, rr, cq, pu, nc, b, cv, s) in rows:
    qs = '[' + ' '.join(f'{v:.4f}' for v in qs) + ']' if qs else '-'
    print(f"{tag:9s} {q or '-':5s} {m or '-':4s} {st:2d} {LP:10.6f} {M:11.9f} {K:8.6f} {I:8.6f} {qs:>30s} {c:5d} {rr:6d} {cq:4d} "
          f"{(f'{pu:10.6f}' if pu is not None else '-'):>10s} {(f'{nc:10.6f}' if nc is not None else '-'):>10s} {b:6.4f} {'YES' if cv else 'no ':4s} {s:6.0f}")
