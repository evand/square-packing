import sys, json, gzip, random
sys.path.insert(0, '/home/evand/math/square-packing/public/s12/search')
import qx2_zm as QZ, zm_mixed as ZM, mixed_cover as MC, zeromargin as zm
from ev import *
path = '/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt'
cov = ZM.Cover(MC.load(path)); ev = Cov(path); m = F(7)
a, b = QZ.u_square(cov)
rec = {tuple(d['root']): d['st'] for d in (json.loads(l) for l in gzip.open('/home/evand/math/square-packing/public/s12/search/qx2_data/cert/runV2_cdade4b6_full.jsonl.gz', 'rt')) if 'root' in d}
roots = json.load(open('cheaproots.json'))
part, nparts = int(sys.argv[1]), int(sys.argv[2]); roots = roots[part::nparts]
chk = QZ.QXChecker(cov, exact_umax=F(1, 2), exact_from=3, use_exact=True, max_depth=18, use_chain=False, cert_mode=False, dump=True)
random.seed(part)
def w(u):
    c, s = zm.trig(u); return c + s
cnt = {}; bad = 0; mism = 0; minmass = {}
for r in roots:
    root = tuple(F(v) for v in r)
    st, unc, leaves = chk.run_box(root)
    if any(st[k] != rec[tuple(r)][k] for k in st if k != 'cpu'): mism += 1
    for box, kind, wit in leaves:
        cnt[kind] = cnt.get(kind, 0) + 1
        cx0, cx1, cy0, cy1, u0, u1 = box
        if kind == 'AXIS':
            K = 2 * min(cx1, m - cx0, cy1, m - cy0)
            if not (u0 == 0 and u1 == 0 and K <= 1): bad += 1; print('AXIS?', box, K)
            continue
        if kind == 'EMPTY':
            wl = min(w(u0), w(u1))
            if not (cx1 < wl / 2 or cx0 > m - wl / 2 or cy1 < wl / 2 or cy0 > m - wl / 2): bad += 1; print('EMPTY?', box)
            continue
        if kind == 'SYM' and not (u0 * u0 + 2 * u0 - 1 >= 0): bad += 1; print('SYM?', box)
        for j in range(10):
            fx, fy, fu = [random.choice([0, 1, F(random.randint(0, 1000), 1000)]) for _ in range(3)]
            cx, cy, u = cx0 + fx * (cx1 - cx0), cy0 + fy * (cy1 - cy0), u0 + fu * (u1 - u0)
            if kind in ('SYM',):
                # the reflected pose: (cy, cx), theta' = 90 - theta, i.e. u' = (1-u)/(1+u)
                cx, cy, u = cy, cx, (1 - u) / (1 + u)
            if kind == 'EXACT45' and u * u + 2 * u - 1 > 0:
                cx, cy, u = cy, cx, (1 - u) / (1 + u)
            if not ev.admissible(cx, cy, u): continue
            mm = ev.mass(cx, cy, u)
            if kind not in minmass or mm < minmass[kind]: minmass[kind] = mm
            if mm < 1: bad += 1; print('MASS < 1', kind, box, (cx, cy, u), float(mm))
print('part', part, 'roots', len(roots), 'stat mismatches vs V2', mism, 'leaves', cnt, 'bad', bad)
print({k: float(v) for k, v in minmass.items()})
