"""Re-solve a finished band LP (quadrant_lp.py band mode) from its saved rows and dump the row duals.

The duals are a fractional periodic packing of closed unit squares in the half-plane that certifies
m_v <= LP value on that row set (SEAM.md).  Output: runs/seam_dual/w<w>.npz with poses, duals, slacks.
Usage: python3 search/seam_dual.py <run_dir> <w>
"""
import sys, os, time, argparse, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import quadrant_lp as Q

run, w = sys.argv[1], float(sys.argv[2])
a = argparse.Namespace(mode='band', R=2, w=w, hc=0.1, hca=0.25, hp=0.1, hpa=0.25)
m = Q.build(a)
z = np.load(os.path.join(run, 'sol.npz'))
lp = Q.LP(m)
n = lp.add(z['rows'])
h = lp.h
h.setOptionValue('solver', 'ipm'); h.setOptionValue('run_crossover', 'on')
t = time.time(); h.run()
print('status', h.getModelStatus(), 'obj', h.getInfo().objective_function_value, 'rows', n, '%.0fs' % (time.time() - t))
sol = h.getSolution()
rd = np.array(sol.row_dual); x = np.array(sol.col_value)
# row 0 is the sigma equality (added in LP.__init__); rows 1.. are poses in lp.P order
out = os.path.join(os.path.dirname(run), 'seam_dual', 'w%g.npz' % w)
np.savez(out, P=lp.P, ydual=rd[1:], lam=rd[0], x=x, obj=h.getInfo().objective_function_value)
print('saved', out, 'lambda', rd[0], 'sum dual', rd[1:].sum(), 'nnz', (np.abs(rd[1:]) > 1e-9).sum())
