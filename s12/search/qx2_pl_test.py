#!/usr/bin/env python3
"""qx2_pl_test.py -- tests for qx2_pl.py (hat line elements).  Run: python3 qx2_pl_test.py [--quick]

1. brute force: every ramp of a hat model subdivided into N >= 2000 uniform pieces of mass x*(F(t_k+1) - F(t_k)),
   evaluated by the existing uniform code (an XModel with fixed hs/vs pieces), vs HatModel.contrib, at >= 10^4 poses
   (random, tiny theta, wall/corner-resting, near-lattice +-1e-7, tile-germ poses); also LP rows A x + const vs the
   oracle path contrib(P, x).
2. bookkeeping: sigma (LP equality LHS) and D = R^2 - cobj.x - const_obj vs direct geometric summation over the
   element arrays: profile mass in the half-open period cell [R+1, R+2) x [0,w], m_v = mass on the line x = R+1,
   corner-box mass mu([0,R]^2), D = R^2 - mu([0,R]^2) + m_v  (QUADRANT.md 1.2/1.4).
3. symmetry: diagonal (cx,cy,th) <-> (cy,cx,90-th) and profile mirror x -> 2R+5-x in the band.
4. all-equal hats (c) == the uniform model with all pieces 2c, at random poses, and the same D.
"""
import sys, os, math, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import qx2_pl as PL
Q = PL.Q; X = PL.X


def ns(**k):
    d = dict(R=2, w=2.0, delta=0.5, kappa=0.0, th0=3.0, hc=0.5, hp=0.5, elem='hat'); d.update(k)
    return argparse.Namespace(**d)


def test_poses(m, n, rng):
    G = m.G
    P = [Q.random_poses(G, n, rng)]
    k = n // 4
    T = Q.random_poses(G, k, rng); T[:, 2] = 10.0 ** rng.uniform(-9, -2, k)          # tiny theta
    T[k // 2:, 2] = math.pi / 2 - T[k // 2:, 2]
    P.append(T)
    W = Q.random_poses(G, k, rng); hw = Q.halfwidth(W[:, 2]); W[:, 1] = hw + 1e-9    # wall-resting, some in the corner
    W[: k // 3, 0] = hw[: k // 3] + 1e-9; P.append(W)
    NL = Q.near_lattice_poses(G, m.hlat, rng, k // 2, full=True); P.append(NL[rng.choice(len(NL), min(k, len(NL)), replace=False)])
    GP = Q.germ_poses(G, m.hlat); P.append(GP[rng.choice(len(GP), min(k, len(GP)), replace=False)])
    return np.concatenate(P)


def brute_model(m, x, N):
    """XModel with the same Lebesgue layer/margin, each ramp -> N fixed uniform pieces."""
    U = X.XModel(m.R, m.w, m.a, m.kappa, m.th0)
    for k, kk in (('rh', 'hs'), ('rv', 'vs')):
        for (c, p, q, d, var, fm) in m.A[k]:
            mass = x[int(var)] if var >= 0 else fm
            if mass <= 0: continue
            t = np.linspace(p, q, N + 1); s = (t - p) / (q - p)
            F = s ** 2 if d > 0 else 1 - (1 - s) ** 2
            for i in range(N): U._add(kk, (c, t[i], t[i + 1], -1, mass * (F[i + 1] - F[i])))
    U.finalize()
    return U


def t_brute(cfg, npose, N, seed, chunk=2000):
    rng = np.random.default_rng(seed)
    m = PL.build(cfg)
    x = rng.uniform(0, 0.3, m.nvar) * (rng.random(m.nvar) < 0.8)
    t0 = time.time(); U = brute_model(m, x, N)
    P = test_poses(m, npose, rng)
    err = 0.0; err_lp = 0.0; worst = None
    for i in range(0, len(P), chunk):
        B = P[i:i + chunk]
        v = m.contrib(B, xval=x)
        A, c = m.contrib(B); vl = A @ x + c
        u = U.contrib(B, xval=np.zeros(1))
        e = np.abs(v - u); j = int(e.argmax())
        if e[j] > err: err = float(e[j]); worst = B[j]
        err_lp = max(err_lp, float(np.abs(vl - v).max()))
    print(f"[brute] R={cfg.R} w={cfg.w} h={cfg.hc}/{cfg.hp} kappa={cfg.kappa}: {len(P)} poses, {len(m.A['rh'])+len(m.A['rv'])} "
          f"ramps x N={N}: max |hat - brute| = {err:.2e} (at {worst}), max |LP row - oracle| = {err_lp:.2e} "
          f"[{time.time()-t0:.0f}s]")
    # brute-force discretisation error <= sum over cut ramps of x/(4N^2) per chord end
    assert err < 2e-6, err
    assert err_lp < 1e-12, err_lp


def direct_mass(m, x, rect, halfopen_x=False):
    """mass of the ramp elements (corner + both walls) in the closed rectangle [x0,x1] x [y0,y1] (x1 open if
    halfopen_x: only matters for atoms on x = x1; ramps are absolutely continuous along their line)."""
    x0, x1, y0, y1 = rect; tot = 0.0
    F = lambda s, d: s * s if d > 0 else 1 - (1 - s) ** 2
    for k in ('rh', 'rv'):
        for (c, p, q, d, var, fm) in m.A[k]:
            mass = x[int(var)]
            if k == 'rh': inside = y0 <= c <= y1; lo, hi = x0, x1        # line y = c, param x
            else:
                inside = (x0 <= c < x1) if halfopen_x else (x0 <= c <= x1); lo, hi = y0, y1
            if not inside: continue
            s0 = min(max((max(lo, p) - p) / (q - p), 0), 1); s1 = min(max((min(hi, q) - p) / (q - p), 0), 1)
            if s1 > s0: tot += mass * (F(s1, d) - F(s0, d))
    return tot


def t_bookkeeping(cfg, seed):
    rng = np.random.default_rng(seed)
    m = PL.build(cfg); R, w, a = m.R, m.w, m.a
    x = rng.uniform(0, 1, m.nvar)
    sig_lp = float(m.csig @ x) + m.const_sig
    sig_dir = direct_mass(m, x, (R + 1, R + 2, 0, w), halfopen_x=True) + (w - a) * 1.0
    mv_dir = direct_mass(m, x, (R + 1, R + 1, 0, w))
    box = direct_mass(m, x, (0, R, 0, R)) + (R - a) ** 2
    D_dir = R ** 2 - box + mv_dir
    D_lp = m.D(x)
    mv_lp = float(sum(m.cobj[i] * x[i] for i in range(m.nvar) if m.label[i][0] == 'Hv'))
    print(f"[book] R={R} w={w} h={cfg.hc}: sigma-row LHS {sig_lp:.12f} vs direct {sig_dir:.12f}; m_v {mv_lp:.12f} vs "
          f"{mv_dir:.12f}; D {D_lp:.12f} vs direct R^2 - mu([0,R]^2) + m_v = {D_dir:.12f}")
    assert abs(sig_lp - sig_dir) < 1e-9 and abs(mv_lp - mv_dir) < 1e-9 and abs(D_lp - D_dir) < 1e-9


def t_symmetry(cfg, n, seed):
    rng = np.random.default_rng(seed)
    m = PL.build(cfg); R = m.R
    x = rng.uniform(0, 1, m.nvar)
    P = Q.random_poses(m.G, n, rng)
    Pd = np.c_[P[:, 1], P[:, 0], np.mod(math.pi / 2 - P[:, 2], math.pi / 2)]
    e1 = np.abs(m.contrib(P, xval=x) - m.contrib(Pd, xval=x)).max()
    # band mirror x -> 2R+5-x, poses well inside the band copies and away from the corner/left wall
    B = np.c_[rng.uniform(R + 1.5, R + 3.5, n), rng.uniform(0, R + 0.75, n), rng.uniform(0, math.pi / 2, n)]
    B[:, 1] = np.maximum(B[:, 1], Q.halfwidth(B[:, 2]) + 1e-7)
    Bm = np.c_[2 * R + 5 - B[:, 0], B[:, 1], np.mod(-B[:, 2], math.pi / 2)]
    e2 = np.abs(m.contrib(B, xval=x) - m.contrib(Bm, xval=x)).max()
    print(f"[sym] diagonal max diff {e1:.2e}; band mirror max diff {e2:.2e}")
    assert e1 < 1e-10 and e2 < 1e-10


def t_equal(cfg, n, seed):
    rng = np.random.default_rng(seed)
    mh = PL.build(cfg); mu = PL.build(argparse.Namespace(**{**vars(cfg), 'elem': 'uniform'}))
    c = 0.137
    xh = np.full(mh.nvar, c); xu = np.full(mu.nvar, 2 * c)
    P = test_poses(mh, n, rng)
    e = np.abs(mh.contrib(P, xval=xh) - mu.contrib(P, xval=xu)).max()
    sh = float(mh.csig @ xh) + mh.const_sig; su = float(mu.csig @ xu) + mu.const_sig
    print(f"[equal] all hats {c} vs uniform {2*c}: max mass diff {e:.2e} over {len(P)} poses; D {mh.D(xh):.12f} vs "
          f"{mu.D(xu):.12f}; sigma-LHS {sh:.12f} vs {su:.12f}")
    assert e < 1e-10 and abs(mh.D(xh) - mu.D(xu)) < 1e-10 and abs(sh - su) < 1e-10


if __name__ == '__main__':
    quick = '--quick' in sys.argv
    small = ns()                                             # R=2, w=2, a=1.5, h=0.5 (phase node 0.5 self-mirror)
    small_off = ns(delta=0.3, hc=0.5, hp=0.5)                # a = 1.7 off the node grid
    prod = ns(R=3, w=3.0, delta=0.2, hc=0.2, hp=0.2, kappa=0.1, th0=3.0)
    t_equal(small, 2000, 1); t_equal(prod, 4000, 2)
    t_symmetry(small, 3000, 3); t_symmetry(prod, 3000, 4)
    for c in (small, small_off, prod): t_bookkeeping(c, 5)
    t_brute(small, 2000 if quick else 10000, 2000, 6)
    t_brute(small_off, 2000 if quick else 10000, 2000, 7)
    t_brute(prod, 500 if quick else 2000, 2000, 8, chunk=100)
    print("all tests passed")
