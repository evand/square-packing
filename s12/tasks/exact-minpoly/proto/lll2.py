import sys, flint, time
from mpmath import mp, mpf
mp.dps = 5000
f = sys.argv[1]; degs = [int(x) for x in sys.argv[2].split(',')]
s = open(f).readline().split()[1]; S = mpf(s); nd = len(s) - 3
digits = int(nd*0.9)
for d in degs:
    t = time.time(); K = mpf(10)**digits
    rows = [[1 if j == i else 0 for j in range(d+1)] + [int(mp.nint(K*S**i))] for i in range(d+1)]
    M = flint.fmpz_mat(rows).lll()
    c = [int(M[0, j]) for j in range(d+1)]
    v = abs(sum(ci*S**i for i, ci in enumerate(c))); h = max(abs(x) for x in c)
    print(d, 'resid', mp.nstr(v, 3), 'height digits', len(str(h)), 'true deg', max(i for i in range(d+1) if c[i]), 'time', round(time.time()-t, 1), flush=True)
