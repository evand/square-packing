import sys, sympy as sp
txt = open(sys.argv[1]).read().split('\n', 2)
names = txt[0].split(','); syms = {s: sp.Symbol(s) for s in names}
polys = [sp.sympify(p.replace('^', '**'), locals=syms) for p in txt[2].replace('\n', '').split(',') if p.strip()]
S = syms['S']; elim = {}
changed = True
while changed:
    changed = False
    for p in polys:
        for v in sorted(p.free_symbols - {S}, key=str):
            pp = sp.Poly(p, v)
            if pp.degree() == 1 and pp.coeffs()[0].is_number:
                sol = sp.solve(p, v)[0]
                elim[v] = sol
                polys = [sp.expand(q.subs(v, sol)) for q in polys]
                polys = [q for q in polys if q != 0]
                # numerators only
                polys = [sp.numer(sp.together(q)) for q in polys]
                changed = True; break
        if changed: break
uniq = []
for p in polys:
    p = sp.expand(p)
    if p != 0 and not any(sp.expand(p - q) == 0 or sp.expand(p + q) == 0 for q in uniq): uniq.append(p)
left = sorted({s for p in uniq for s in p.free_symbols} - {S}, key=str) + [S]
print(len(uniq), 'equations,', len(left), 'unknowns:', left, file=sys.stderr)
for p in uniq: print('  deg', sp.Poly(p, *left).total_degree(), str(p)[:150], file=sys.stderr)
with open(sys.argv[2], 'w') as fh:
    fh.write(','.join(map(str, left)) + '\n' + txt[1] + '\n')
    out = []
    for p in uniq:
        den = sp.lcm([sp.fraction(c)[1] for c in sp.Poly(p, *left).coeffs()])
        out.append(str(sp.expand(p * den)).replace('**', '^'))
    fh.write(',\n'.join(out) + '\n')
