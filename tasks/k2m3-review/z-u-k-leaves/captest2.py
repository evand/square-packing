import sys, random
sys.path.insert(0, '/home/evand/math/square-packing/public/s12/search')
import qx2_zm as QZ, zm_mixed as ZM, mixed_cover as MC, zeromargin as zm
from ev import *
path = '/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt'
cov = ZM.Cover(MC.load(path)); ev = Cov(path)
a, b = QZ.u_square(cov); m = cov.m
seed = int(sys.argv[1]); N = int(sys.argv[2])
random.seed(seed)
# my own profiles for lines y=a and x=a
def myprof(segs):
    bs = sorted(set([p for p, q, w in segs] + [q for p, q, w in segs]))
    return bs, [sum((w / (q - p) for p, q, w in segs if p <= bs[i] and q >= bs[i + 1]), F(0)) for i in range(len(bs) - 1)]
PH = myprof(ev.H[a]); PV = myprof(ev.V[a])
def mindens(P, lo, hi):
    bs, r = P
    if lo < bs[0] or hi > bs[-1]: return F(0)
    vals = [r[i] for i in range(len(r)) if bs[i + 1] > lo and bs[i] < hi]
    if not vals and lo == hi:
        vals = [r[i] for i in range(len(r)) if bs[i] <= lo <= bs[i + 1]]
    return min(vals) if vals else F(0)
acc = {'CAP': 0, 'LEB': 0}; nposes = 0; worst = None; lemma_viol = 0
for trial in range(N):
    k = random.randint(2, 9); wdt = F(1, 10 * 2 ** k)
    ku = random.randint(6, 14); du = F(1, 2 ** ku)
    u0 = du * random.choice([0, 0, 1, 2, 3, random.randint(0, int(F(1, 16) / du) - 1)]); u1 = u0 + du
    # centre boxes near the boundary of U on the lower / left side, or near the corner
    def coord():
        t = random.random()
        base = a + F(1, 2) - F(random.randint(0, 300), 10000)
        if t < 0.85: return base + wdt * random.randint(-3, 3)
        return F(random.randint(int(a * 10) - 5, 35 * 10 - 1), 10)
    cx0 = coord(); cy0 = coord()
    cx0 = (cx0 // wdt) * wdt; cy0 = (cy0 // wdt) * wdt
    box = (cx0, cx0 + wdt, cy0, cy0 + wdt, u0, u1)
    if box[1] > F(7, 2) or box[3] > F(7, 2): continue
    cu1 = zm.clip_bin(box, m)
    if cu1 < u1: box = box[:5] + (cu1,)
    if box[5] <= box[4]: continue
    B = zm.bin_data(box[4], box[5])
    kind = None
    if box[4] ** 2 + 2 * box[4] - 1 >= 0: continue
    if QZ.cert_leb(box, B, a, b): kind = 'LEB'
    elif QZ.cert_cap(box, B, cov, a, b): kind = 'CAP'
    if kind is None: continue
    acc[kind] += 1
    X0, X1, Y0, Y1, U0, U1 = box
    samples = []
    for fx in (0, 1, F(1, 2), F(random.random())):
        for fy in (0, 1, F(1, 2), F(random.random())):
            for fu in (0, 1, F(1, 1000), F(random.random())):
                samples.append((X0 + fx * (X1 - X0), Y0 + fy * (Y1 - Y0), U0 + fu * (U1 - U0)))
    for cx, cy, u in samples:
        cx, cy, u = F(cx).limit_denominator(10**12), F(cy).limit_denominator(10**12), F(u).limit_denominator(10**12)
        if not (X0 <= cx <= X1 and Y0 <= cy <= Y1 and U0 <= u <= U1): continue
        if not ev.admissible(cx, cy, u): continue
        nposes += 1
        l, s = ev.parts(cx, cy, u)
        mm = l + s
        if worst is None or mm < worst[0]: worst = (mm, kind, box, (cx, cy, u))
        # the lemma's own hypotheses, checked at this pose
        vs = ev.verts(cx, cy, u)
        ok = all(v[0] <= b and v[1] <= b for v in vs)
        dy = a - min(v[1] for v in vs); dx = a - min(v[0] for v in vs)
        if dy > 0:
            ok &= cy >= a
            ch = ev.chord_h(cx, cy, u, a); ok &= ch is not None and dy <= mindens(PH, *ch)
        if dx > 0:
            ok &= cx >= a
            ch = ev.chord_v(cx, cy, u, a); ok &= ch is not None and dx <= mindens(PV, *ch)
        if not ok:
            lemma_viol += 1
            if lemma_viol < 5: print('LEMMA HYPOTHESIS FAILS', kind, box, (cx, cy, u), float(mm))
print('seed', seed, 'accepted', acc, 'poses', nposes, 'lemma-hyp failures', lemma_viol)
print('worst', float(worst[0]), worst[1], [str(v) for v in worst[2]], [float(v) for v in worst[3]])
