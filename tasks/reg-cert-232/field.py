"""L = K(sqrt 2), K = Q[t]/f (the n = 232 record's field): primitive element w = t + sqrt 2, its minimal polynomial g,
and t, sqrt 2 as polynomials in w.  Prints a JSON summary."""
import json, sympy as sp
t, w, y = sp.symbols('t w y')
fc = [1, 3, -1, 3, -21, -7, 5, 1]                     # low -> high, from search/exact/cur/data/n-232.minpoly.json.gz
f = sum(c * t**i for i, c in enumerate(fc))
g = sp.Poly(sp.resultant(f.subs(t, w - y), y**2 - 2, y), w)
print('deg g', g.degree(), 'irreducible', g.is_irreducible)
# express t and sqrt2 in w: w^k = sum a_ij t^i r^j (r = sqrt 2), 14 x 14 rational matrix
r = sp.sqrt(2)
fP = sp.Poly(f, t)
def coords(expr):
    p = sp.Poly(sp.expand(expr), t, y)
    v = [0] * 14
    for (i, j), c in p.terms():
        # reduce y^j with y^2 = 2, t^i mod f
        c2 = c * 2 ** (j // 2); j %= 2
        rem = sp.Poly(t**i, t).rem(fP)
        for (ii,), cc in rem.terms():
            v[ii + 7 * j] += c2 * cc
    return v
M = sp.Matrix([coords((t + y)**k) for k in range(14)]).T     # column k = coords of w^k
tv = sp.Matrix(coords(t)); rv = sp.Matrix(coords(y))
A = M.LUsolve(tv); B = M.LUsolve(rv)
# numeric root
import mpmath as mpm
mpm.mp.dps = 60
tstar = mpm.findroot(lambda x: sum(c * x**i for i, c in enumerate(fc)), mpm.mpf('-0.26248080531744977935'))
wstar = tstar + mpm.sqrt(2)
gl = [int(c) for c in reversed(g.all_coeffs())]
print('g(w*) =', mpm.nstr(sum(c * wstar**i for i, c in enumerate(gl)), 5))
At = [sp.Rational(x) for x in A]; Bt = [sp.Rational(x) for x in B]
print('t(w*) - t* =', mpm.nstr(sum(mpm.mpf(x.p) / x.q * wstar**i for i, x in enumerate(At)) - tstar, 5))
print('r(w*) - sqrt2 =', mpm.nstr(sum(mpm.mpf(x.p) / x.q * wstar**i for i, x in enumerate(Bt)) - mpm.sqrt(2), 5))
json.dump(dict(g=gl, t_in_w=[str(x) for x in At], sqrt2_in_w=[str(x) for x in Bt], wstar=mpm.nstr(wstar, 55)),
          open('field.json', 'w'), indent=1)
print(gl)
