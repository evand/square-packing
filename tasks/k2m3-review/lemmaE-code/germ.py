"""Boxes with a corner at the corner germs (3/2,3/2), (23/10,3/2) (4 orientations), u0 = 0: sampled exact min m_s;
certify(tau = m_s + 1e-12) must fail; also sharpness certify(tau = m_s - 1e-6)."""
import sys, random, json, time
from fractions import Fraction as F
from mass import *
import qx2_zm as QZ
from sound import sample_min
cov, a, b = load(); M = Mass(cov, a, b); ex = QZ.Exact(cov, a, b)
rng = random.Random(int(sys.argv[1]))
pts = [(F(3,2),F(3,2)),(F(23,10),F(3,2)),(F(3,2),F(23,10))]
for (gx, gy) in pts[int(sys.argv[2])::2]:
    for sx in (1,-1):
        for sy in (1,-1):
            h = F(1,160); u1 = F(1,256)
            x0 = gx if sx > 0 else gx - h; y0 = gy if sy > 0 else gy - h
            box = (x0, x0+h, y0, y0+h, F(0), u1)
            t0 = time.time()
            ms = sample_min(M, box, rng, n=1000)[0]
            plus = ex.certify(box, tau=ms + F(1,10**12))[0]
            m6 = ex.certify(box, tau=ms - F(1,10**6))[0]
            rec = dict(box=[str(v) for v in box], ms=float(ms), plus=plus, m6=m6, reg=str(ex.u_regime(box))[:60], t=round(time.time()-t0,1))
            if plus: rec['BUG'] = True
            print(json.dumps(rec), flush=True)
