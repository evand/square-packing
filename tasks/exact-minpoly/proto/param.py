import sys, re, flint
from fractions import Fraction
from mpmath import mp, mpf
mp.dps = 300
s = open(sys.argv[1]).read().strip().rstrip(':')
s = re.sub(r'(\d+)\s*/\s*2\^(\d+)', r'Fraction(\1, 2**\2)', s)
d = eval(s)
p = d[1]; vars_ = p[3]; lin = p[4]
elim, den, params = p[5][1][0], p[5][1][1], p[5][1][2]
P = lambda c: flint.fmpz_poly([int(x) for x in c[1]]) if c[0] >= 0 else flint.fmpz_poly([])
f, g = P(elim), P(den)
print('separating:', [v for v, c in zip(vars_, lin) if c], 'degree', f.degree())
iS = vars_.index('S')
pS, cS = params[iS]       # parametrised vars exclude the last (separating) one?
num = P(pS)
# S = -num(t) / (cS * g(t));  min poly of S: resultant_t( f(t), cS*g(t)*z + num(t) ) in z
R = flint.fmpz_mpoly_ctx.get(('t', 'z'), 'lex') if hasattr(flint, 'fmpz_mpoly_ctx') else None
import sympy as sp
t, z = sp.symbols('t z')
fs = sp.Poly(list(reversed([int(x) for x in f.coeffs()])), t)
hs = sp.Poly(cS * sp.Poly(list(reversed([int(x) for x in g.coeffs()])), t).as_expr() * z + sp.Poly(list(reversed([int(x) for x in num.coeffs()])), t).as_expr(), t, z)
res = sp.resultant(fs.as_expr(), hs.as_expr(), t)
fac = sp.factor_list(sp.Poly(res, z))
S0 = mpf(sys.argv[2])
for q, m in fac[1]:
    c = [int(x) for x in q.all_coeffs()]
    v = sum(ci * S0**(len(c)-1-i) for i, ci in enumerate(c))
    print('deg', q.degree(), 'mult', m, 'value at S*', mp.nstr(v, 3))
    if abs(v) < mpf(10)**-45: print('MINPOLY', c)
