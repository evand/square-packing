import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
TREE, OUT = sys.argv[1], sys.argv[2]
sys.path.insert(0, TREE); os.chdir(TREE)
import numpy as np
from pk.engine import FQ, QOpt, SCREEN, POLISH
from pk.packing import Packing
import mcmin
f = FQ(); R = {}
for n in (40, 110, 237):
    s, sq = mcmin.load_deg(f'{HERE}/data/n{n}.txt')
    p = Packing(s, np.array(sq) + np.random.default_rng(n).normal(0, 0.01, (n, 3)))
    for name, o in (('screen', SCREEN.but(loosen=1.02)), ('polish', POLISH), ('rect', QOpt(loosen=1.0, rect_dx=0.2, pit=8, flip_top=0))):
        r = f.quench(p, o)
        R[(n, name)] = (r.p.s, r.p.sq.tobytes(), r.p.meta, {k: v for k, v in r.stats.items() if k not in ('sec', 't_alm', 'trace')})
pickle.dump(R, open(OUT, 'wb'))
