"""Soundness test of Exact.certify against the independent evaluator (mass.py).
For a box: sample admissible poses (u in (u0, ue]), float-minimise, evaluate the best exactly -> m_s (a true mass
value at an admissible pose).  Then certify(box, tau = m_s + 1e-12) MUST fail (else: BUG, with witness pose).
Sharpness: certify(box, tau = m_s - delta) for delta = 1e-6, 1e-3 (records how close the certificate gets).
Optional lambda test at u0 = 0 (claim mass >= tau + lam u).
usage: sound.py SEED N [mut]"""
import sys, random, time, math, json
from fractions import Fraction as F
from mass import *
import qx2_zm as QZ

T225 = F(41421356, 10 ** 8)          # < tan 22.5deg


def ue_of(u1):
    return u1 if u1 * u1 + 2 * u1 - 1 <= 0 else T225


def rnd_box(rng):
    r = rng.random()
    if r < 0.4:   # near the corner of U
        cx = rng.uniform(1.0, 2.6); cy = rng.uniform(1.0, 2.6)
    elif r < 0.7:  # near a wall
        cx = rng.uniform(0.5, 3.5); cy = rng.uniform(0.5, 0.75)
        if rng.random() < 0.5: cx, cy = cy, cx
    else:
        cx = rng.uniform(0.5, 3.5); cy = rng.uniform(0.5, 3.5)
    h = F(1, rng.choice([10, 20, 40, 80, 160, 320]))
    cx0 = F(round(cx / float(h))) * h; cy0 = F(round(cy / float(h))) * h
    hy = h if rng.random() < 0.7 else h / 2
    t = rng.random()
    if t < 0.4:
        u0 = F(0); u1 = F(1, rng.choice([16, 64, 256, 2048, 2 ** 14]))
    elif t < 0.85:
        k = rng.choice([16, 32, 64, 128, 256]); i = rng.randint(1, int(0.41 * k) - 1)
        u0 = F(i, k); u1 = F(i + 1, k)
    else:
        u0 = F(rng.randint(385, 412), 1000); u1 = u0 + F(rng.choice([5, 10, 20, 40]), 1000)
    return (cx0, cx0 + h, cy0, cy0 + hy, u0, u1)


def sample_min(M, box, rng, n=1500):
    cx0, cx1, cy0, cy1, u0, u1 = box
    ue = ue_of(u1)
    us = set()
    for f in (1e-9, 1e-6, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0):
        us.add(float(u0) + f * float(ue - u0))
    us = sorted(us)
    cand = []

    def lo_edge(v0, u):
        c = (1 - u * u) / (1 + u * u); s = 2 * u / (1 + u * u)
        return max(float(v0), (c + s) / 2 + 1e-15)
    for i in range(n):
        u = rng.choice(us) if rng.random() < 0.6 else float(u0) + rng.random() * float(ue - u0)
        if u <= 0: continue
        xl, yl = lo_edge(cx0, u), lo_edge(cy0, u)
        if xl > float(cx1) or yl > float(cy1): continue
        def pick(lo, hi):
            r = rng.random()
            if r < 0.25: return lo
            if r < 0.5: return hi
            return lo + rng.random() * (hi - lo)
        cx, cy = pick(xl, float(cx1)), pick(yl, float(cy1))
        cand.append((M.flt(cx, cy, u), cx, cy, u))
    if not cand: return None
    cand.sort()
    # local refine of the best few (float), then exact
    best = []
    for m, cx, cy, u in cand[:8]:
        step = [float(cx1 - cx0) / 8, float(cy1 - cy0) / 8, float(ue - u0) / 8]
        for it in range(60):
            improved = False
            for k in range(3):
                for sg in (-1, 1):
                    p = [cx, cy, u]; p[k] += sg * step[k]
                    p[2] = min(max(p[2], float(u0) + 1e-12 * float(ue - u0) + 1e-300), float(ue))
                    if p[2] <= 0: continue
                    xl, yl = lo_edge(cx0, p[2]), lo_edge(cy0, p[2])
                    p[0] = min(max(p[0], xl), float(cx1)); p[1] = min(max(p[1], yl), float(cy1))
                    mm = M.flt(*p)
                    if mm < m - 1e-15: m, cx, cy, u = mm, p[0], p[1], p[2]; improved = True
            if not improved:
                step = [s / 2 for s in step]
        best.append((m, cx, cy, u))
    out = None
    for m, cx, cy, u in best + cand[:4]:
        uq = F(u).limit_denominator(10 ** 15)
        if not (u0 < uq <= ue): uq = ue if uq > ue else u0 + (ue - u0) / 10 ** 9
        c, s = QZ.zm.trig(uq); w2 = (c + s) / 2
        cxq = min(max(F(cx).limit_denominator(10 ** 15), cx0, w2), cx1)
        cyq = min(max(F(cy).limit_denominator(10 ** 15), cy0, w2), cy1)
        if cxq < w2 or cyq < w2: continue
        me = M.exact(cxq, cyq, uq)
        if out is None or me < out[0]: out = (me, cxq, cyq, uq)
    return out


def main():
    seed, N = int(sys.argv[1]), int(sys.argv[2])
    mut = len(sys.argv) > 3 and sys.argv[3] == 'mut'
    rng = random.Random(seed)
    cov, a, b = load()
    if mut:
        for key, L in cov.lines.items():
            f = rng.choice([1, 1, F(97, 100), F(103, 100), F(9, 10)])
            L['segs'] = [(t0, t1, d * f) for t0, t1, d in L['segs']]
    M = Mass(cov, a, b); ex = QZ.Exact(cov, a, b)
    out = open(f'sound_{seed}{"_mut" if mut else ""}.jsonl', 'w')
    for it in range(N):
        box = rnd_box(rng)
        t0 = time.time()
        sm = sample_min(M, box, rng)
        if sm is None: continue
        ms = sm[0]
        rec = dict(box=[str(v) for v in box], ms=float(ms), pose=[str(v) for v in sm[1:]])
        ok, why = ex.certify(box, tau=ms + F(1, 10 ** 12))
        rec['plus'] = ok
        if ok:
            rec['BUG'] = True
            print('*** BUG', rec, flush=True)
        ok6, _ = ex.certify(box, tau=ms - F(1, 10 ** 6)); rec['m6'] = ok6
        if not ok6:
            ok3, why3 = ex.certify(box, tau=ms - F(1, 10 ** 3)); rec['m3'] = ok3; rec['why3'] = why3
        rec['reg'] = str(ex.u_regime(box))[:80]
        rec['t'] = round(time.time() - t0, 1)
        out.write(json.dumps(rec) + '\n'); out.flush()
        print(it, rec['box'][:1], rec['box'][4:], f"ms={float(ms):.9f}", 'plus', ok, 'm6', ok6, rec.get('m3'), rec['reg'][:40], rec['t'], flush=True)


if __name__ == '__main__':
    main()
