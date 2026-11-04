"""Tests for pk/families.py (cached family codes, lazy fill, donor / best queries).  python pk/tests/test_families.py"""
import os, sqlite3, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pk import families as F

N = 5
REAL = F.classify


def _store(d, rows):
    """rows: [(s, family code the stub classifier returns, grid, added)]"""
    from pk.store import Store
    p = os.path.join(d, 's.sqlite')
    Store(p).close()
    assert F.migrate(p) and not F.migrate(p)
    db = sqlite3.connect(p)
    for i, (s, code, grid, added) in enumerate(rows):
        sq = np.full((N, 3), float(i))
        db.execute("INSERT INTO packing (hash, n, s, coords, status, pen, grid, added, tags) VALUES (?,?,?,?,?,?,?,?,?)",
                   (f'h{i}', N, s, sq.astype('<f8').tobytes(), 'raw', 0.0, grid, added, f'["{code}"]'))
    db.commit()
    return p, db


def _stub(calls):
    def classify(n, s, sq):
        calls.append(s)
        return {0: 'G', 1: 'J', 2: 'D', 3: 'J', 4: 'V', 5: 'J', 6: 'J'}[int(sq[0, 0])]
    return classify


def test_near_and_cache():
    with tempfile.TemporaryDirectory() as d:
        p, db = _store(d, [(3.00, 'G', 0, '2026-01-01'), (3.01, 'J', 0, '2026-01-01'), (3.02, 'D', 1, '2026-01-01'),
                           (3.0100000001, 'J', 0, '2026-01-01'), (3.03, 'V', 0, '2026-02-01'),
                           (3.20, 'J', 0, '2026-01-01'), (3.04, 'J', 0, '2026-01-01')])
        calls = []; F.classify = _stub(calls)
        r = F.near(db, N, 4.0, exclude=['G'], path=p, window=0.05)
        assert [x[1] for x in r] == [3.01, 3.03, 3.04]          # G excluded, grid excluded, side dedupe, window 0.05
        assert [x[3] for x in r] == ['J', 'V', 'J']
        n1 = len(calls)
        r2 = F.near(db, N, 4.0, exclude=['G'], path=p, window=0.05)
        assert [x[0] for x in r2] == [x[0] for x in r] and len(calls) == n1     # second read: all from the cache
        assert db.execute('SELECT count(*) FROM packing WHERE family IS NOT NULL').fetchone()[0] == n1
        assert F.best(db, N, 4.0, ['G'], path=p) == (3.01, r[0][0])
        assert F.best(db, N, 4.0, ['G'], before='2026-01-15', path=p)[0] == 3.01
        assert [x[1] for x in F.near(db, N, 4.0, ['G', 'J'], path=p, before='2026-01-15')] == []
        assert [x[1] for x in F.near(db, N, 4.0, ['G', 'J'], path=p)] == [3.03]
        assert F.best(db, N, 3.005, ['G'], path=p) == (None, None)


def test_fill_readonly():
    with tempfile.TemporaryDirectory() as d:
        p, db = _store(d, [(3.0, 'J', 0, '2026-01-01')])
        F.classify = _stub([])
        ro = sqlite3.connect(f'file:{p}?mode=ro', uri=True)
        r = F.near(ro, N, 4.0, path=p + '.missing-dir/x')       # write-back fails quietly: result still right
        assert [x[3] for x in r] == ['G']


def test_shape_dedupe():
    """near(shape=k): up to k packings per rounded side when they are different arrangements (differs())"""
    rng = np.random.default_rng(0)
    g = np.array([[x + .5, y + .5, 0.0] for x in range(4) for y in range(4)])        # 16 squares, side 4
    noise = g + np.r_[rng.normal(0, 1e-3, (16, 2)).T, np.zeros((1, 16))].T              # one packing, relax noise
    mirror = g.copy(); mirror[:, 0] = 4 - mirror[:, 0]                                  # a box symmetry of g
    moved = g.copy(); moved[:6, 0] += 0.4                                               # six squares rearranged
    assert not F.differs(4.0, g, noise) and not F.differs(4.0, g, mirror) and F.differs(4.0, g, moved)
    few = g.copy(); few[:2, 0] += 0.4
    assert not F.differs(4.0, g, few)                                                   # <= SHAPE_K wanderers
    with tempfile.TemporaryDirectory() as d:
        from pk.store import Store
        p = os.path.join(d, 's.sqlite'); Store(p).close(); F.migrate(p)
        db = sqlite3.connect(p)
        for i, (s, X) in enumerate([(4.0, g), (4.0000000001, noise), (4.0000000002, moved), (4.0000000003, moved + 1e-4),
                                    (4.1, g)]):
            db.execute("INSERT INTO packing (hash, n, s, coords, status, pen, grid, added, tags) VALUES (?,?,?,?,?,?,?,?,?)",
                       (f'h{i}', 16, s, X.astype('<f8').tobytes(), 'raw', 0.0, 0, '2026-01-01', '["J"]'))
        db.commit()
        F.classify = lambda n, s, sq: 'J'
        assert [round(x[1], 10) for x in F.near(db, 16, 5.0, path=p)] == [4.0, 4.1]
        assert [round(x[1], 10) for x in F.near(db, 16, 5.0, path=p, shape=3)] == [4.0, 4.0000000002]   # 4.1 = g again
        assert [round(x[1], 10) for x in F.near(db, 16, 5.0, path=p, shape=1)] == [4.0]
        assert [round(x[1], 10) for x in F.near(db, 16, 5.0, path=p, shape=3, window=0.05)] == [4.0, 4.0000000002]


def teardown_function(_):
    F.classify = REAL


if __name__ == '__main__':
    test_near_and_cache(); test_fill_readonly(); test_shape_dedupe(); print('ok')
