"""Targeted: sub-boxes of the run-B root [23/10, 12/5]^2 x [0, 1/16] (square touching both boundary lines of U at
theta -> 0), where the cap case split at d = 0 (cap_alts with two branches) and _gcdcert fire.  For each sub-box:
sampled exact min m_s; certify(tau = m_s + 1e-12) must fail; count split/gcd/sproc usage."""
import sys, random, json, time
from fractions import Fraction as F
from mass import *
import qx2_zm as QZ
from sound import sample_min

cov, a, b = load(); M = Mass(cov, a, b); ex = QZ.Exact(cov, a, b)
cnt = dict(split=0, gcd=0, sproc=0)
_ca = ex.cap_alts
def cap_alts(*A, **K):
    r = _ca(*A, **K)
    if len(r) > 1: cnt['split'] += 1
    return r
ex.cap_alts = cap_alts
_g = QZ.Exact._gcdcert
def gcdc(*A):
    r = _g(*A)
    if r: cnt['gcd'] += 1
    return r
ex._gcdcert = gcdc
_s = QZ.Exact._sproc
def sp(*A):
    r = _s(*A)
    if r: cnt['sproc'] += 1
    return r
ex._sproc = sp
rng = random.Random(int(sys.argv[1]))
boxes = []
for hx in (F(1, 80), F(1, 160), F(1, 640)):
    for u1 in (F(1, 16), F(1, 64), F(1, 512)):
        for k in range(2):
            x0 = F(23, 10) + hx * rng.randint(0, int(F(1, 10) / hx) - 1)
            y0 = F(23, 10) + hx * rng.randint(0, int(F(1, 10) / hx) - 1)
            if k == 0: x0 = y0 = F(23, 10)
            boxes.append((x0, x0 + hx, y0, y0 + hx, F(0), u1))
for box in boxes:
    t0 = time.time()
    for k in cnt: cnt[k] = 0
    sm = sample_min(M, box, rng, n=800)
    ms = sm[0]
    plus = ex.certify(box, tau=ms + F(1, 10 ** 12))[0]
    c1 = ex.certify(box)[0]
    rec = dict(box=[str(v) for v in box], ms=float(ms), plus=plus, cert1=c1, cnt=dict(cnt), t=round(time.time() - t0, 1))
    if plus: rec['BUG'] = True
    print(json.dumps(rec), flush=True)
