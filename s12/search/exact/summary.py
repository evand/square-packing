#!/usr/bin/env python3
"""Markdown table of exactsolve results:  summary.py results/*.json"""
import json, sys

print('| input | n | S input | S exact (KKT point) | certified S\' | free squares | min lambda | 2nd order | corner-corner MILP |')
print('|---|---|---|---|---|---|---|---|---|')
for p in sys.argv[1:]:
    r = json.load(open(p))
    name = p.rsplit('/', 1)[-1][:-5]
    if r.get('S_exact') is None:
        why = r.get('status', '?')
        if r.get('vv_milp_dS') is not None:
            why += f' (first-order dS* = {r["vv_milp_dS"]:.1e})'
        print(f'| {name} | {r["n"]} | {r["S_input"]} | {why} | - | - | - | - | - |')
        continue
    so = r['second_order']
    st = so['status']
    if st.startswith('local minimum'):
        st = f'PSD, {so["n_zero"]} exact flat modes'
    elif st.startswith('PD'):
        st = 'PD (strict)'
    elif st.startswith('rigid'):
        st = 'rigid (null J_A = 0)'
    lam = r.get('lambda_A_maxmin')
    lam = float('nan') if lam is None else lam
    vv = r.get('vv_milp_dS')
    print(f'| {name} | {r["n"]} | {r["S_input"]} | {r["S_exact"][:32]} | {r["S_cert_decimal"][:32]} '
          f'{"(valid)" if r["cert_valid"] else "(INVALID)"} | {r["free_squares"]} | '
          f'{lam:.2e} | {st} | {"jammed" if vv is not None and vv > -1e-9 else vv} |')
