"""Band LP restricted to the x-uniform + seam class: mu = lambda(dy) x Leb_x + rho(dy) x delta_Z(x).
Columns: full-width horizontal lines at y_i, full-width slabs [y_j,y_j+1], seam segments, seam points."""
import sys, os, math, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import quadrant_lp as Q
w = float(sys.argv[1]); hp = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 25
m = Q.Model(0, w); m.band = True; m.G = ('band', w); m.hlat = hp
ys = np.round(np.arange(0, w + 1e-9, hp), 10)
for y in ys[1:]: m.prof_hseg(y, 0.0, 0.5, var=-2)          # orbit {[0,.5],[.5,1]} = full line
for a, b in zip(ys[:-1], ys[1:]): m.prof_cell(0.0, 0.5, a, b, var=-2); m.prof_vseg(0.0, a, b, var=-2)
for y in ys[1:]: m.prof_point(0.0, y, var=-2)
m.cobj = [-c for c in m.cobj]; m.finalize()
lp = Q.LP(m); lp.solver = 'ipm'; lp.crossover = 'on'
lp.add(Q.lattice_poses(m.G, 0.05, list(range(0, 90, 3)) + [1e-4, 89.99]))
rng = np.random.default_rng(1); t0 = time.time()
for r in range(rounds):
    x, obj, st = lp.solve()
    Q._init(m, x)
    mn, fp, fv, worst = Q.oracle(Q.Serial(), m, x, rng, nrand=200000, npol=400, nproc=1)
    add = lp.add(Q.pick_rows(fp, fv, 3000)) if len(fp) else 0
    print('round %d  m_v=%.5f  oracle min %.6f at (%.4f,%.4f,%.2fdeg) add %d  [%.0fs]' % (r, -obj, mn, *worst[:2], math.degrees(worst[2]), add, time.time() - t0), flush=True)
    if mn >= 1 - 1e-6 and r >= 2: break
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs', 'seam_xu_w%g.npz' % w), x=x, obj=obj)
for j in np.nonzero(x > 1e-6)[0]: print('   ', m.label[j], '%.5f' % x[j])
