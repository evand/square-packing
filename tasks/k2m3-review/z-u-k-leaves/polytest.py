import sys, random
sys.path.insert(0, '/home/evand/math/square-packing/public/s12/search')
import qx2_zm as QZ, zm_mixed as ZM, mixed_cover as MC, zeromargin as zm
from ev import *
path = '/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt'
cov = ZM.Cover(MC.load(path)); ev = Cov(path)
chk = ZM.MixedChecker(cov, max_depth=18, use_chain=False, cert_mode=False)
m = cov.m; random.seed(int(sys.argv[1])); N = int(sys.argv[2])
def inQ(p, cx, cy, u):
    C, S, N_ = 1 - u * u, 2 * u, 1 + u * u
    X = ((p[0] - cx) * C + (p[1] - cy) * S) / N_; Y = (-(p[0] - cx) * S + (p[1] - cy) * C) / N_
    return abs(X) <= HALF and abs(Y) <= HALF
nb = 0; nk = 0; bad = 0; npose = 0
for t in range(N):
    k = random.randint(1, 8); wd = F(1, 10 * 2 ** k)
    x0 = F(random.randint(0, int(F(7, 2) / wd) - 1)) * wd; y0 = F(random.randint(0, int(F(7, 2) / wd) - 1)) * wd
    ku = random.randint(3, 12); du = F(1, 2 ** ku)
    u0 = du * random.randint(0, int(F(1, 2) / du) - 1); u1 = u0 + du
    box = (x0, x0 + wd, y0, y0 + wd, u0, u1)
    cu1 = zm.clip_bin(box, m)
    if cu1 < u1: box = box[:5] + (cu1,)
    B = zm.bin_data(box[4], box[5]); lo = B['wlo'] / 2
    if box[1] < lo or box[0] > m - lo or box[3] < lo or box[2] > m - lo: continue
    cx0, cx1, cy0, cy1, u0, u1 = box
    specs = chk.zc._adm_specs(box, B)
    r = ZM.REACH
    frame = [(cx0 - r, cy0 - r), (cx1 + r, cy0 - r), (cx1 + r, cy1 + r), (cx0 - r, cy1 + r)]
    K = frame
    for c in range(4):
        Kc = ZM.cond_region(specs, c, u0, u1 - u0, m, (cx0, cy0), frame)
        K = ZM.clip_convex(K, Kc) if Kc else []
        if not K: break
    nb += 1
    if not K: continue
    nk += 1
    for j in range(12):
        fx, fy, fu = [random.choice([0, 1, F(random.randint(0, 1000), 1000)]) for _ in range(3)]
        cx, cy, u = cx0 + fx * (cx1 - cx0), cy0 + fy * (cy1 - cy0), u0 + fu * (u1 - u0)
        if not ev.admissible(cx, cy, u): continue
        npose += 1
        for p in K:
            if not inQ(p, cx, cy, u):
                bad += 1
                if bad < 5: print('K vertex outside Q', box, (cx, cy, u), p)
                break
print('boxes', nb, 'nonempty K', nk, 'poses', npose, 'bad', bad)
