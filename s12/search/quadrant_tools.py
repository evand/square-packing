#!/usr/bin/env python3
"""quadrant_tools.py -- helpers for quadrant_lp.py (task quadrant-lp): pose breakdown, deep local search,
box accounting check.  HEURISTIC float tools."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import quadrant_lp as Q


def breakdown(m, x, pose):
    """list of (kind, element tuple, fraction, mass) captured by the closed square at pose, plus Lebesgue."""
    P = np.array([pose], float); out = []
    for k, arr in m.A.items():
        for e in arr:
            sub = Q.Model(m.R, m.w, m.wleb); sub.E = {kk: [] for kk in m.E}; sub.E[k] = [tuple(e[:-2]) + (-1, 1.0)]
            sub.nvar = 0; sub.finalize()
            A, c = sub.contrib(P)
            fr = c[0] - sub.lebesgue(P)[0]
            var = int(e[-2]); mass = x[var] if var >= 0 else e[-1]
            if fr > 1e-12 and mass > 0: out.append((k, tuple(np.round(e[:-2], 6)), float(fr), float(fr * mass)))
    leb = float(m.lebesgue(P)[0])
    return out, leb


def deep(m, x, seeds, rounds=14, nper=200, rng=None):
    if rng is None: rng = np.random.default_rng(0)
    best = []
    for s in seeds:
        cur = np.array(s, float)[None]; cv = m.contrib(cur, xval=x)[0]
        for r in range(rounds):
            rad = 0.03 * 0.5 ** r; arad = math.radians(1.0) * 0.5 ** r
            Qp = cur + rng.normal(size=(nper, 3)) * np.array([rad, rad, arad])
            k = rng.random(nper) < 0.3; Qp[k, 1] = Q.halfwidth(Qp[k, 2]) + Q.WALL
            Qp[:, 2] = np.mod(Qp[:, 2], math.pi / 2); Qp = Q.clipQ(Qp, m.R)
            v = m.contrib(Qp, xval=x); j = int(v.argmin())
            if v[j] < cv: cv = v[j]; cur = Qp[j:j + 1]
        best.append((cv, cur[0]))
    best.sort(key=lambda t: t[0])
    return best


def box_saving(m, x, k):
    """independent box accounting: k^2 - mu_k([0,k]^2) for the D4 box built from the quadrant family
    (4 corner modules, profile on [R, k-R] along each wall, Lebesgue on [w,k-w]^2 minus the 4 corner squares).
    Returns (saving, 4D, sigma deficit, m_v)."""
    def mass(var, fm): return x[int(var)] if var >= 0 else fm
    MC_pts_cells = sum(mass(e[-2], e[-1]) for kk in ('pts', 'cells') for e, t in zip(m.A[kk], m.KA[kk]) if t == 0)
    MC_segs = sum(mass(e[-2], e[-1]) for kk in ('hs', 'vs') for e, t in zip(m.A[kk], m.KA[kk]) if t == 0)
    MC = MC_pts_cells + MC_segs
    wall = 0.0; cell = 0.0; mv = 0.0
    for (typ, a, b, y0, y1, var, fm) in m.base:
        ms = mass(var, fm); cell += ms
        if typ in 'pv' and abs(a) < 1e-12: mv += ms
        for j in range(m.R, k - m.R + 1):
            if j + b <= k - m.R + 1e-12: wall += ms
    leb = (k - 2 * m.wleb) ** 2 - 4 * (m.R - m.wleb) ** 2
    total = 4 * MC + 4 * wall + leb
    Dform = m.R ** 2 - MC - mv
    return k * k - total, 4 * Dform, m.w - cell, mv


# ============================================================================ independent box evaluation
class BoxModel(Q.Model):
    """the D4 box measure mu_k built element by element from a quadrant solution (all elements fixed)."""

    def __init__(self, m, x, k):
        super().__init__(m.R, m.w, m.wleb)
        self.k = k; self.J = 0
        def mass(var, fm): return float(x[int(var)]) if var >= 0 else float(fm)
        maps = [lambda a, b: (a, b), lambda a, b: (k - a, b), lambda a, b: (a, k - b), lambda a, b: (k - a, k - b)]
        # corner elements: D4 corner images (the module is diagonal-symmetric already)
        for kk, arr in m.A.items():
            for e, t in zip(arr, m.KA[kk]):
                if t != 0: continue
                ms = mass(e[-2], e[-1])
                if ms <= 0: continue
                for f in maps:
                    if kk == 'pts':
                        a, b = f(e[0], e[1]); self.E['pts'].append((a, b, -1, ms))
                    elif kk == 'hs':          # y = e0, x in [e1, e2]
                        (a0, y), (a1, _) = f(e[1], e[0]), f(e[2], e[0])
                        self.E['hs'].append((y, min(a0, a1), max(a0, a1), -1, ms))
                    elif kk == 'vs':          # x = e0, y in [e1, e2]
                        (x0, b0), (_, b1) = f(e[0], e[1]), f(e[0], e[2])
                        self.E['vs'].append((x0, min(b0, b1), max(b0, b1), -1, ms))
                    else:
                        (a0, b0), (a1, b1) = f(e[0], e[2]), f(e[1], e[3])
                        self.E['cells'].append((min(a0, a1), max(a0, a1), min(b0, b1), max(b0, b1), -1, ms))
        # profile: base cell atoms at j + phase on [R, k-R] along the bottom wall, then the 4 walls by rotation
        R = m.R
        rots = [lambda a, b: (a, b), lambda a, b: (k - b, a), lambda a, b: (k - a, k - b), lambda a, b: (b, k - a)]
        for (typ, a, b, y0, y1, var, fm) in m.base:
            ms = mass(var, fm)
            if ms <= 0: continue
            for j in range(R, k - R + 1):
                if j + b > k - R + 1e-12: continue
                for f in rots:
                    if typ == 'p':
                        u, v = f(j + a, y0); self.E['pts'].append((u, v, -1, ms))
                    elif typ == 'h':          # horizontal at height y0, x in [j+a, j+b]
                        p0, p1 = f(j + a, y0), f(j + b, y0)
                        if abs(p0[1] - p1[1]) < 1e-12: self.E['hs'].append((p0[1], min(p0[0], p1[0]), max(p0[0], p1[0]), -1, ms))
                        else: self.E['vs'].append((p0[0], min(p0[1], p1[1]), max(p0[1], p1[1]), -1, ms))
                    elif typ == 'v':          # vertical at x = j+a, y in [y0, y1]
                        p0, p1 = f(j + a, y0), f(j + a, y1)
                        if abs(p0[0] - p1[0]) < 1e-12: self.E['vs'].append((p0[0], min(p0[1], p1[1]), max(p0[1], p1[1]), -1, ms))
                        else: self.E['hs'].append((p0[1], min(p0[0], p1[0]), max(p0[0], p1[0]), -1, ms))
                    else:
                        p0, p1 = f(j + a, y0), f(j + b, y1)
                        self.E['cells'].append((min(p0[0], p1[0]), max(p0[0], p1[0]), min(p0[1], p1[1]), max(p0[1], p1[1]), -1, ms))
        for kk in self.E: self.K[kk] = [0] * len(self.E[kk])
        self.finalize()

    def total(self):
        t = sum(float(self.A[kk][:, -1].sum()) for kk in self.A)
        return t + (self.k - 2 * self.wleb) ** 2 - 4 * (self.R - self.wleb) ** 2

    def lebesgue(self, P):
        V = Q.corners(P); n = len(P); k = float(self.k); wl = self.wleb; R = float(self.R)
        f = lambda a, b, c, d: Q.area_rect(V, np.full(n, a), np.full(n, b), np.full(n, c), np.full(n, d))
        out = f(wl, k - wl, wl, k - wl)
        for (a, b, c, d) in [(wl, R, wl, R), (k - R, k - wl, wl, R), (wl, R, k - R, k - wl), (k - R, k - wl, k - R, k - wl)]:
            out -= f(a, b, c, d)
        return out


def box_clip(P, k):
    hw = Q.halfwidth(P[:, 2])
    P[:, 0] = np.clip(P[:, 0], hw + Q.WALL, k - hw - Q.WALL); P[:, 1] = np.clip(P[:, 1], hw + Q.WALL, k - hw - Q.WALL)
    return P


def box_oracle(B, nrand=2000000, h=0.1, npol=3000, rounds=10, seed=0):
    """float oracle on the whole box: random poses (30% wall-resting), near-lattice poses, local polish."""
    k = B.k; rng = np.random.default_rng(seed); x = np.zeros(1)
    P = np.c_[rng.uniform(0, k, nrand), rng.uniform(0, k, nrand), rng.uniform(0, math.pi / 2, nrand)]
    r = rng.random(nrand); hw = Q.halfwidth(P[:, 2])
    P[r < 0.15, 1] = hw[r < 0.15] + Q.WALL
    s = (r >= 0.15) & (r < 0.3); P[s, 0] = hw[s] + Q.WALL
    s = rng.random(nrand) < 0.2; P[s, 2] = rng.choice([0.0, 1e-6, 1e-3, math.pi / 4], size=int(s.sum()))
    g = np.arange(0, k + 1e-9, h); out = [P]
    for d in (1e-7, 1e-5):
        for sx in (-1, 1):
            for sy in (-1, 1):
                X, Y = np.meshgrid(g + 0.5 + sx * d, g + 0.5 + sy * d, indexing='ij')
                for th in (0.0, 1e-6):
                    out.append(np.c_[X.ravel(), Y.ravel(), np.full(X.size, th)])
    P = box_clip(np.concatenate(out), k)
    v = np.concatenate([B.contrib(c, xval=x) for c in np.array_split(P, max(1, len(P) // 20000))])
    o = np.argsort(v)[:npol]; cur = P[o].copy(); curv = v[o].copy()
    for rr in range(rounds):
        rad = 0.03 * 0.45 ** rr; arad = math.radians(1.5) * 0.45 ** rr
        Qp = np.repeat(cur, 16, 0); Qp = Qp + rng.normal(size=Qp.shape) * np.array([rad, rad, arad])
        Qp[:, 2] = np.mod(Qp[:, 2], math.pi / 2); Qp = box_clip(Qp, k)
        vv = B.contrib(Qp, xval=x).reshape(len(cur), 16); j = vv.argmin(1)
        bt = vv[np.arange(len(cur)), j] < curv
        cur[bt] = Qp.reshape(len(cur), 16, 3)[np.arange(len(cur)), j][bt]; curv[bt] = vv[np.arange(len(cur)), j][bt]
    j = int(curv.argmin())
    return float(min(v.min(), curv[j])), cur[j], float(v.min())


# ============================================================================ heavy check CLI
def heavy_check(m, x, seed=1, nrand=4000000, npol=6000, rounds=12, log=print):
    """a heavier, differently seeded oracle than the LP loop's: random poses (x3 wall-resting mix), the lattice at
    pitch 0.013 over 40 angles, near-lattice at d in {1e-8,1e-6,1e-4,1e-3}, germ family at 6 tiny angles, and a deep
    polish of the worst 6000 (diverse) poses.  Returns (min, pose, n_below_1-1e-6)."""
    rng = np.random.default_rng(seed); G = m.G
    parts = [Q.random_poses(G, nrand, rng),
             Q.lattice_poses(G, 0.013, list(np.linspace(0, 88.5, 40)) + [1e-6, 1e-3, 89.99])]
    ex, ey = Q._ext(G); band = isinstance(G, tuple); h = m.hlat
    g = np.arange(0, ex + 1e-9, h); gy = np.arange(0, ey + 1e-9, h)
    for d in (1e-8, 1e-6, 1e-4, 1e-3):
        for sx in (-1, 1):
            for sy in (-1, 1):
                X, Y = np.meshgrid(g + 0.5 + sx * d, gy + 0.5 + sy * d, indexing='ij')
                X = X.ravel(); Y = Y.ravel(); k = (Y <= ey) & ((Y <= X + 1e-6) | band)
                for th in (0.0, 1e-7, math.pi / 2 - 1e-7):
                    parts.append(Q.clipQ(np.c_[X[k], Y[k], np.full(k.sum(), th)], G))
    parts.append(Q.germ_poses(G, h, thetas=(1e-7, -1e-7, 1e-4, -1e-4, 1e-2, -1e-2)))
    parts.append(Q.germ_poses(G, h / 2, thetas=(1e-5, -1e-5)))
    P = np.concatenate(parts)
    v = np.concatenate([m.contrib(c, xval=x) for c in np.array_split(P, max(1, len(P) // 20000))])
    log(f"  heavy: {len(P)} poses, raw min {v.min():.7f}, below 1-1e-6: {(v < 1 - 1e-6).sum()}")
    o = np.argsort(v)
    key = np.floor(P[o, :2] / 0.02).astype(np.int64); kk = key[:, 0] * 100000 + key[:, 1]
    _, first = np.unique(kk, return_index=True)
    seeds = P[o[np.sort(first)][:npol]]
    cur = seeds.copy(); curv = m.contrib(cur, xval=x)
    for r in range(rounds):
        rad = 0.04 * 0.5 ** r; arad = math.radians(2.0) * 0.5 ** r
        Qp = np.repeat(cur, 20, 0) + rng.normal(size=(len(cur) * 20, 3)) * np.array([rad, rad, arad])
        s = rng.random(len(Qp)) < 0.25; Qp[s, 1] = Q.halfwidth(Qp[s, 2]) + Q.WALL
        s = rng.random(len(Qp)) < 0.2; Qp[s, 2] = 0.0
        Qp[:, 2] = np.mod(Qp[:, 2], math.pi / 2); Qp = Q.clipQ(Qp, G)
        vv = m.contrib(Qp, xval=x).reshape(len(cur), 20); j = vv.argmin(1)
        bt = vv[np.arange(len(cur)), j] < curv
        cur[bt] = Qp.reshape(len(cur), 20, 3)[np.arange(len(cur)), j][bt]; curv[bt] = vv[np.arange(len(cur)), j][bt]
    j = int(curv.argmin())
    mn = min(float(v.min()), float(curv[j])); pose = cur[j] if curv[j] <= v.min() else P[int(v.argmin())]
    log(f"  heavy: after polish min {mn:.7f} at ({pose[0]:.6f},{pose[1]:.6f},{math.degrees(pose[2]):.5f}deg); "
        f"polished below 1-1e-6: {(curv < 1 - 1e-6).sum()}")
    return mn, pose, int((v < 1 - 1e-6).sum() + (curv < 1 - 1e-6).sum())


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['check', 'box'])
    ap.add_argument('mode'); ap.add_argument('R', type=int); ap.add_argument('w', type=float); ap.add_argument('sol')
    ap.add_argument('--hc', type=float, default=0.1); ap.add_argument('--hp', type=float, default=0.1)
    ap.add_argument('--hca', type=float, default=0.25); ap.add_argument('--hpa', type=float, default=0.25)
    ap.add_argument('--k', type=int, default=0); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--nrand', type=int, default=4000000)
    a = ap.parse_args()
    m = Q.build(a); d = np.load(a.sol); x = d['x']
    print(f"{a.sol}: D_LP = {float(d['D']):.6f} (round {int(d['round'])})")
    if a.cmd == 'check':
        heavy_check(m, x, seed=a.seed, nrand=a.nrand)
    else:
        k = a.k or 2 * m.R + 2
        B = BoxModel(m, x, k)
        sv = box_saving(m, x, k)
        print(f"  box k={k}: saving (element sum) = {k*k - B.total():.6f}; box_saving = {sv[0]:.6f}; 4D = {sv[1]:.6f}; "
              f"sigma = {sv[2]:.2e}; m_v = {sv[3]:.6f}")
        mn, pose, raw = box_oracle(B, nrand=a.nrand // 2, seed=a.seed)
        print(f"  box oracle min {mn:.7f} (raw {raw:.7f}) at ({pose[0]:.5f},{pose[1]:.5f},{math.degrees(pose[2]):.4f}deg)")
