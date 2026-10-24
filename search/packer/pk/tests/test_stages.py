"""Contract of pk/stages.py (10-22): every registered stage logs all its parameters, is deterministic for a seed, and
leaves frozen squares where they were.  A new stage is covered by registering it (no test to write; a plugin stage that
needs a particular start names its chain in the plugin's TEST_STARTS).  Plugins load as usual (PK_PLUGINS= skips them).

  python3 pk/tests/test_stages.py            (needs target/release/fq and anneal; PK_FQ overrides fq)
"""
import os, random, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pk import stages as S
from pk.packing import Packing

N = 40                                                     # 32 frozen + 8 to place: small and fast
TEST_STARTS = {k: v for m in S.PLUGINS for k, v in getattr(m, 'TEST_STARTS', {}).items()}
SLOW = {'recreate.exact', 'recreate.greedy', 'perturb.cross'}  # fillfast scans; cross needs mates


def start(seed=1):
    """6 x 6 upright lattice in a 6.9 box, the middle 2 x 2 cells removed (their unit boxes = the region), rest frozen"""
    K, F = 6, 0.9
    cells = [(i, j) for i in range(K) for j in range(K) if not (2 <= i < 4 and 2 <= j < 4)]
    sq = np.array([(i + 0.5 + F * i / K, j + 0.5 + F * j / K, 0.0) for i, j in cells])
    rects = [(i + F * i / K, j + F * j / K, i + 1 + F * i / K, j + 1 + F * j / K) for i in (2, 3) for j in (2, 3)]
    st = S.State(Packing(K + F, sq), region=rects)
    return S.run_stage('constrain.freeze', {}, st, random.Random(seed), dict(n=N))[0]


def filled(seed=1):
    st, _ = S.run_chain('recreate.bestofK(K=4) | relax.round.15c | constrain.release', dict(n=N), seed, st=start(seed))
    return st


def check(name):
    f = S.STAGES[name]
    if name in SLOW:
        return 'skip'
    if f.kind == 'gen':
        st0 = None
        if name == 'gen.file':
            return 'skip'                                  # needs a path; covered by trial.py
    elif f.kind in ('perturb',):
        st0 = filled()
    elif name in TEST_STARTS:                               # plugin stages that need their own start (meta)
        st0 = S.run_chain(TEST_STARTS[name], dict(n=N), 1)[0]
    else:
        st0 = start()
    fixed = {}
    runs = []
    for _ in range(2):
        st, row = S.run_stage(name, fixed, st0, random.Random(5), dict(n=N))
        runs.append((st, row))
    (a, ra), (b, rb) = runs
    assert set(f.space) <= set(ra['params']), (name, 'parameters not logged', sorted(set(f.space) - set(ra['params'])))
    assert ra['params'] == rb['params'], (name, 'params differ for one seed')
    assert np.array_equal(a.p.sq, b.p.sq) and a.p.s == b.p.s, (name, 'not deterministic for one seed')
    if st0 is not None and st0.frozen.any() and f.kind in ('recreate', 'relax', 'constrain') and name != 'constrain.release':
        k = int(st0.frozen.sum())
        assert np.allclose(a.p.sq[:k], st0.p.sq[:k], atol=1e-9) or row_failed(ra), (name, 'moved frozen squares')
    return 'ok'


def row_failed(r):
    return bool(r['info'].get('failed'))


def main():
    res = {}
    for name in sorted(S.STAGES):
        try:
            res[name] = check(name)
        except AssertionError as e:
            res[name] = f'FAIL {e}'
        except Exception as e:                             # e.g. perturb on frozen squares must raise, others must not
            res[name] = f'ERROR {e!r}'[:200]
    for k, v in res.items():
        print(f'{k:24s} {v}')
    bad = [k for k, v in res.items() if v.startswith(('FAIL', 'ERROR'))]
    # the frozen guard: perturb stages refuse frozen squares
    try:
        S.run_stage('perturb.kick', {}, start(), random.Random(1), dict(n=N)); bad.append('perturb.kick accepted frozen squares')
    except ValueError:
        pass
    print('FAILED: ' + ', '.join(bad) if bad else f'all {len(res)} stages ok')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
