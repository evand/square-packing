"""Evaluate a fixed w-band profile with the quadrant_lp float oracle.  Elements: list of (kind, args, fm)."""
import sys, os, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import quadrant_lp as Q
def build(w, elems):
    m = Q.Model(0, w); m.band = True; m.G = ('band', float(w)); m.hlat = 0.1
    for kind, args, fm in elems: getattr(m, kind)(*args, fm=fm)
    m.finalize(); return m
def check(w, elems, seed=0, nrand=300000, show=8):
    m = build(w, elems); x = np.zeros(1)
    Q._init(m, x); rng = np.random.default_rng(seed)
    mn, fp, fv, worst = Q.oracle(Q.Serial(), m, x, rng, nrand=nrand, npol=600, nproc=1)
    print('m_v=%.5f  sigma-mass=%.5f  oracle min=%.6f at (%.5f, %.5f, %.3f deg)' % (m.const_obj, m.const_sig, mn, worst[0], worst[1], math.degrees(worst[2])))
    if len(fv):
        o = np.argsort(fv)[:show]
        for p, v in zip(fp[o], fv[o]): print('   %.5f  (%.4f, %.4f, %.3f deg)' % (v, p[0], p[1], math.degrees(p[2])))
    return mn, fp, fv
def seamline(y0, y1, dens, n=50):          # piecewise-uniform seam density, as n vsegs
    ys = np.linspace(y0, y1, n + 1); out = []
    for a, b in zip(ys[:-1], ys[1:]):
        d = dens((a + b) / 2) if callable(dens) else dens
        out.append(('prof_vseg', (0.0, a, b), d * (b - a)))
    return out
if __name__ == '__main__':
    A = seamline(0.5, 1.0, 1.5) + [('prof_hseg', (1.0, 0.0, 0.5), 0.125)]
    check(1.0, A)
