"""All 8 roots of the s(11) octic: solve the reduced system at each root -> runs/explain/phantom11.json."""
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'runs', 'explain', 'phantom11.json')
import sympy as sp, json, symsys
S, eqs, syms, cs, angs = symsys.build(11)
c, s = cs[0]
G = sp.groebner(eqs + [c**2 + s**2 - 1], *syms, c, s, S, order='lex')
oct_ = sp.Poly(S**8 - 20*S**7 + 178*S**6 - 842*S**5 + 1923*S**4 - 496*S**3 - 6754*S**2 + 12420*S - 6865, S)
roots = oct_.nroots(n=30)
print('roots:', [sp.N(r, 10) for r in roots])
order = list(syms) + [c, s]
out = []
for r in roots:
    sub = {S: r}
    # solve triangular: go through variables from last to first
    for v in reversed(order):
        cands = []
        for g in G.exprs:
            gs = sp.expand(g.subs(sub))
            if gs.free_symbols == {v}:
                cands.append(sp.Poly(gs, v))
        if not cands: raise SystemExit(f'no eq for {v}')
        # pick lowest degree nonzero poly; take roots consistent with all
        p = min(cands, key=lambda q: q.degree())
        sols = [z for z in p.nroots(n=30) if all(abs(complex(q.eval(z))) < 1e-15 * (1 + max(abs(complex(cc)) for cc in q.all_coeffs())) for q in cands)]
        if len(sols) != 1: print('ambiguous', v, len(sols), [sp.N(z, 6) for z in p.nroots(n=20)])
        sub[v] = sols[0] if sols else p.nroots(n=30)[0]
    out.append({str(k): [float(sp.re(v_)), float(sp.im(v_))] for k, v_ in sub.items()})
os.makedirs(os.path.dirname(OUT), exist_ok=True); json.dump(out, open(OUT, 'w'))
for o in out: print(round(o['S'][0], 6), round(o['S'][1], 6), 'c', [round(x, 4) for x in o['c0']], 's', [round(x, 4) for x in o['s0']])
