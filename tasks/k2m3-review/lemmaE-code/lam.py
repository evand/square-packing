"""Zero-margin / slope test at the tight theta = 0 points.
1. find the tight one-sided corners of the theta = 0 face (own enumeration, exact);
2. for boxes with a tight corner on a face/corner (u0 = 0): sampled min of (mass - 1)/u -> lam_s (exact at a pose);
   certify(box, tau=1, lam=lam_s + 1e-9) must FAIL; certify(tau=1, lam=0) and (lam = lam_s/2) recorded;
   certify(tau = 1 + 1e-12) must FAIL (the inf over u > 0 is 1).
usage: lam.py PART NPART"""
import sys, random, time, json
from fractions import Fraction as F
from mass import *
import qx2_zm as QZ
from sound import sample_min, ue_of

cov, a, b = load(); M = Mass(cov, a, b); ex = QZ.Exact(cov, a, b)
H_ = F(1, 2)


def tight_points():
    xs = set([a, b]); prof = {}
    for key, L in cov.lines.items():
        P = QZ.Profile(L['segs']); prof[key] = P; xs.update(P.b); xs.add(key[1])
    lo, hi = H_, F(7, 2)
    grid = sorted(set(v + s for v in xs for s in (-H_, H_) if lo <= v + s <= hi) | {lo, hi})
    out = []
    for cx in grid:
        for cy in grid:
            if cy > cx: continue
            for sx in (-1, 1):
                for sy in (-1, 1):
                    if (cx == lo and sx < 0) or (cx == hi and sx > 0) or (cy == lo and sy < 0) or (cy == hi and sy > 0):
                        continue
                    # exact mass of the square slightly moved by (sx, sy) * eps, eps -> 0: evaluate at eps = 1e-30
                    e = F(1, 10 ** 30)
                    v = M.exact(cx + sx * e, cy + sy * e, F(0)) if False else None
                    out.append((cx, cy, sx, sy))
    return grid, out


def mass0(cx, cy):
    # theta = 0 exact mass of the closed square (independent: via the polygon evaluator at u = 0)
    return M.exact(cx, cy, F(0))


def main():
    part, npart = int(sys.argv[1]), int(sys.argv[2])
    grid, pts = tight_points()
    e = F(1, 10 ** 20)
    tight = set()
    for cx, cy, sx, sy in pts:
        if mass0(cx + sx * e, cy + sy * e) < 1 + F(1, 10 ** 15):
            tight.add((cx, cy, sx, sy))
    tight = sorted(t for t in tight if min(t[0], t[1]) - H_ < a)
    print(len(tight), 'tight one-sided corners (cy <= cx)', flush=True)
    rng = random.Random(5)
    rng.shuffle(tight)
    for i, (cx, cy, sx, sy) in enumerate(tight):
        if i % npart != part: continue
        for h, u1 in ((F(1, 80), F(1, 64)), (F(1, 640), F(1, 1024))):
            x0 = cx if sx > 0 else cx - h; y0 = cy if sy > 0 else cy - h
            box = (x0, x0 + h, y0, y0 + h, F(0), u1)
            if x0 < 0 or y0 < 0: continue
            t0 = time.time()
            # sampled min of (m - 1)/u
            best = None
            for k in range(1500):
                u = u1 * F(rng.choice([1, 10, 100, 1000, 10 ** 4, 10 ** 6]), 10 ** 6) * F(rng.randint(1, 1000), 1000)
                c, s = QZ.zm.trig(u); w2 = (c + s) / 2
                def pk(lo, hi):
                    r = rng.random()
                    if r < 0.3: return lo
                    if r < 0.6: return hi
                    return lo + (hi - lo) * F(rng.randint(0, 1000), 1000)
                xl = max(box[0], w2); yl = max(box[2], w2)
                if xl > box[1] or yl > box[3]: continue
                px, py = pk(xl, box[1]), pk(yl, box[3])
                mf = M.flt(float(px), float(py), float(u))
                lam = (mf - 1) / float(u)
                if best is None or lam < best[0]: best = (lam, px, py, u)
            lam_f, px, py, u = best
            me = M.exact(px, py, u); lam_s = (me - 1) / u
            r_plus = ex.certify(box, tau=F(1), lam=lam_s + F(1, 10 ** 9))[0]
            r_t = ex.certify(box, tau=F(1) + F(1, 10 ** 12))[0]
            r_0 = ex.certify(box)[0]
            r_half = ex.certify(box, tau=F(1), lam=lam_s / 2)[0] if lam_s > 0 else None
            rec = dict(box=[str(v) for v in box], lam_s=float(lam_s), pose=[str(px), str(py), str(u)], m=float(me),
                       plus=r_plus, tau_plus=r_t, cert=r_0, half=r_half, t=round(time.time() - t0, 1))
            if r_plus or r_t: rec['BUG'] = True
            print(json.dumps(rec), flush=True)


if __name__ == '__main__':
    main()
