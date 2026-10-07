"""y-invariant seam LP (FRIEDMAN.md section 8): a 1-periodic x-profile rho (measure on [0,1), extended in y as rho x Leb_y)
must give every closed unit square (angle th, centre offset cx) mass >= 1.  rho = g*delta_0 + rest, total 1+e.
max g given e;  kappa = e/g is the excess price of one unit of seam in the bulk.  Float LP, discretised.
Usage: python3 search/seam_1d.py [N=400] [e1,e2,...] [--save DIR]
  --save DIR: also write DIR/1d_N{N}_e{e}.npy (x[i] = weight at i/N) for each e, for search/seam_1d_exact.py.
Importable: V, build_rows, solve, check (no work is done at import)."""
import sys, os, numpy as np, highspy


def V(th, d):                               # closed vertical chord of a unit square at angle th, horizontal offset d
    c, s = np.cos(th), np.sin(th); a = (c + s) / 2; d = np.abs(d)
    if s < 1e-12: return np.where(d <= 0.5 + 1e-12, 1.0, 0.0)
    return np.where(d <= (c - s) / 2 + 1e-12, 1 / c, np.clip((a - d) / (s * c), 0, None))


THS = np.radians(np.r_[np.arange(0, 45.01, 0.5), [1e-4, 0.05, 0.1, 0.2]])


def build_rows(N, ths=THS):
    xs = np.arange(N) / N                   # atom positions; xs[0] = 0 is the seam
    cxs0 = np.arange(2 * N) / (2 * N)
    cxs = np.r_[cxs0 + 1e-7, cxs0 - 1e-7, cxs0 + 0.5 / (2 * N) * 0.37]   # off-grid: no free endpoint atoms (eps-gap)
    rows = []
    for th in ths:
        for cx in cxs:
            d = (xs - cx + 0.5) % 1 - 0.5       # nearest-image offsets; chord support <= 0.71 so also add +/-1 images
            rows.append(V(th, d) + V(th, d + 1) + V(th, d - 1))
    return np.array(rows)


def solve(e, A, extra=None):
    """max x[0] s.t. sum x = 1+e, A x >= 1 (and extra x >= 1 if given), x >= 0."""
    N = A.shape[1]
    h = highspy.Highs(); h.setOptionValue('output_flag', False); inf = highspy.kHighsInf
    h.addVars(N, np.zeros(N), np.full(N, inf))
    cost = np.zeros(N); cost[0] = -1.0; h.changeColsCost(N, np.arange(N, dtype=np.int32), cost)
    h.addRow(1 + e, 1 + e, N, np.arange(N, dtype=np.int32), np.ones(N))
    for r in (A if extra is None else np.vstack([A, extra])):
        nz = np.nonzero(r)[0].astype(np.int32); h.addRow(1.0, inf, len(nz), nz, r[nz])
    h.run()
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:      # seen once (N=200, e=0.001, cut round 5)
        h.clearSolver(); h.setOptionValue('solver', 'ipm'); h.run()
    x = np.array(h.getSolution().col_value)
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal or abs(x.sum() - 1 - e) > 1e-6 or x.min() < -1e-7:
        raise RuntimeError(f'LP failed: status {h.modelStatusToString(h.getModelStatus())}, sum {x.sum()}, min {x.min()}')
    return x


def check(x, dth=0.01, ncx=20000):
    """dense float re-check of a 1D solution: all th on a dth-degree grid (plus tiny tilts), cx on a fine grid plus
    the critical centres where an atom sits just outside a chord end (eps-gap poses).  Returns (min mass, th, cx)."""
    xs = np.arange(len(x)) / len(x)
    sup = np.nonzero(x > 1e-12)[0]; X = xs[sup]; W = x[sup]
    worst = (9.0, None, None)
    thd = np.r_[np.arange(0, 45.0001, dth), [1e-6, 1e-4, 1e-3]]
    for t in np.radians(thd):
        c, s = np.cos(t), np.sin(t); a = (c + s) / 2; b = (c - s) / 2
        crit = np.r_[np.arange(ncx) / ncx, ((X[:, None] + np.array([a, -a, b, -b, .5, -.5])[None, :]).ravel()[:, None]
                                             + np.array([1e-9, -1e-9])[None, :]).ravel() % 1]
        tot = np.zeros(len(crit))
        for sh in (-1, 0, 1):
            d = X[None, :] + sh - crit[:, None]
            tot += (V(t, d) * W[None, :]).sum(1)
        i = tot.argmin()
        if tot[i] < worst[0]: worst = (tot[i], np.degrees(t), crit[i])
    return worst


def main(argv):
    save = None
    if '--save' in argv:
        i = argv.index('--save'); save = argv[i + 1]; argv = argv[:i] + argv[i + 2:]
    N = int(argv[0]) if len(argv) > 0 else 400
    ES = [float(t) for t in argv[1].split(',')] if len(argv) > 1 else [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
    xs = np.arange(N) / N
    A = build_rows(N)
    for e in ES:
        x = solve(e, A); g = x[0]
        mid = x[(xs >= .25) & (xs <= .75)].sum(); other = x[1:].sum() - mid
        print(f"e={e:<6} g={g:.4f}  kappa=e/g={e/g if g>1e-9 else float('inf'):.4f}   mass in [1/4,3/4]: {mid:.3f}  elsewhere(off-seam): {other:.3f}", flush=True)
        if save:
            os.makedirs(save, exist_ok=True); fn = os.path.join(save, f"1d_N{N}_e{e:g}.npy"); np.save(fn, x)
            print('   saved', fn, flush=True)
        if e == 0.05:
            sup = [(round(xs[i], 4), round(x[i], 4)) for i in np.nonzero(x > 1e-4)[0]]
            print('   support (x, mass) at e=0.05:', sup[:40], '...' if len(sup) > 40 else '')


if __name__ == '__main__':
    main(sys.argv[1:])
