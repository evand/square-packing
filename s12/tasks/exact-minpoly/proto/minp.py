import sys, flint
from mpmath import mp, mpf
mp.dps = 5000
f = sys.argv[1]; maxd = int(sys.argv[2])
s = open(f).readline().split()[1]; S = mpf(s); nd = len(s) - 3
digits = int(nd*0.9); K = mpf(10)**digits
for d in range(1, maxd+1):
    rows = [[1 if j == i else 0 for j in range(d+1)] + [int(mp.nint(K*S**i))] for i in range(d+1)]
    M = flint.fmpz_mat(rows).lll()
    c = [int(M[0, j]) for j in range(d+1)]
    v = abs(sum(ci*S**i for i, ci in enumerate(c)))
    h = max(abs(x) for x in c)
    if c[d] and len(str(h)) < 0.6*digits/(d+1) and v < mpf(10)**(-nd+len(str(h))+30):
        p = flint.fmpz_poly(c)
        fac = p.factor()
        best = min(fac[1], key=lambda q: abs(sum(int(q[0][i])*S**i for i in range(q[0].degree()+1)))/max(1,max(abs(int(x)) for x in q[0].coeffs())))
        q = best[0]
        print('deg', q.degree(), 'resid', mp.nstr(abs(sum(int(q[i])*S**i for i in range(q.degree()+1))), 3), 'height', max(abs(int(x)) for x in q.coeffs()))
        print('poly', [int(x) for x in reversed(q.coeffs())])
        break
else: print('none up to', maxd, 'at', nd, 'digits')
