"""run_chain keeps the best feasible packing it saw (pk/stages.py, 10-26).  Fake stages, no fq.
  python3 pk/tests/test_keep_best.py"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pk import stages as S
from pk.packing import Packing

N = 4


def grid(s):
    """2 x 2 upright unit squares spread over side s (feasible for s >= 2)"""
    c = (0.5, s - 0.5)
    return Packing(s, np.array([(x, y, 0.0) for x in c for y in c]))


NAMES = []


def _reg(name, kind, f):
    S.stage(name, kind)(f); NAMES.append(name)


_reg('gen.tk_ok', 'gen', lambda st, rng, ctx: (S.State(grid(2.2)), {}))              # feasible, side 2.2
_reg('gen.tk_overlap', 'gen', lambda st, rng, ctx: (S.State(grid(1.5)), {}))         # overlapping, side 1.5
_reg('relax.tk_worse', 'relax', lambda st, rng, ctx: (st.but(p=grid(st.p.s + 0.3)), {}))
_reg('relax.tk_better', 'relax', lambda st, rng, ctx: (st.but(p=grid(2.05)), {}))
_reg('relax.tk_fail', 'relax', lambda st, rng, ctx: (st, dict(failed='infeasible')))
_reg('shrink.tk_new', 'shrink', lambda st, rng, ctx: (S.State(Packing(st.p.s, st.p.sq)), {}))   # a fresh State, as many stages build


def test_keep_best():
    ctx = dict(n=N)
    try:
        st, lg = S.run_chain('gen.tk_ok | relax.tk_worse', ctx, 1)          # the relax made it worse: keep the build
        assert st.p.s == 2.2 and lg[-1]['stage'] == 'keep.best' and lg[-1]['info']['why'] == 'worse'
        assert lg[-1]['info']['src'] == 'gen.tk_ok' and not S.failed(lg)

        st, lg = S.run_chain('gen.tk_ok | relax.tk_better', ctx, 1)         # it improved: the final state, no row
        assert st.p.s == 2.05 and lg[-1]['stage'] == 'relax.tk_better' and not S.failed(lg)

        st, lg = S.run_chain('gen.tk_ok | relax.tk_fail', ctx, 1)           # failed relax: the feasible build counts
        assert st.p.s == 2.2 and lg[-1]['info']['why'] == 'failed' and not S.failed(lg)

        st, lg = S.run_chain('gen.tk_overlap | relax.tk_fail', ctx, 1)      # nothing feasible: still a failure
        assert S.failed(lg) and lg[-1]['stage'] == 'relax.tk_fail'

        st, lg = S.run_chain('gen.tk_overlap | relax.tk_worse', ctx, 1)     # infeasible all the way: final as is
        assert lg[-1]['stage'] == 'relax.tk_worse' and st.best is None

        b, lg1 = S.run_chain('gen.tk_ok | shrink.tk_new', ctx, 1, keep_best=False)   # split chain: build half tracks only,
        assert lg1[-1]['stage'] == 'shrink.tk_new' and b.best[0] == 2.2              # through a stage that drops fields
        st, lg2 = S.run_chain('relax.tk_worse', ctx, 2, st=b)               # relax half: the build's best carried over
        assert st.p.s == 2.2 and lg2[-1]['stage'] == 'keep.best'

        st, lg = S.run_chain('relax.tk_worse', ctx, 1, st=S.State(grid(2.2)))   # the start itself is not a candidate
        assert st.p.s == 2.5 and lg[-1]['stage'] == 'relax.tk_worse'

        st, lg = S.run_chain('gen.tk_ok | relax.tk_worse', dict(n=N + 1), 1)   # wrong count: never tracked
        assert st.p.s == 2.5 and st.best is None
    finally:
        for k in NAMES:
            S.STAGES.pop(k, None)


if __name__ == '__main__':
    test_keep_best(); print('ok')
