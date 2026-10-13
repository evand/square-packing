"""python -m pytest pk/tests  (or run directly)."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pk.packing import Packing, read, parse

HERE = os.path.dirname(os.path.abspath(__file__))
REC = os.path.join(HERE, '..', '..', 'seeds', 'rec110.txt')


def sym(p, sw, fx, fy):
    x, y, a = p.sq[:, 0].copy(), p.sq[:, 1].copy(), p.sq[:, 2].copy()
    if sw: x, y, a = y, x, 90 - a
    if fx: x, a = p.s - x, 90 - a
    if fy: y, a = p.s - y, 90 - a
    return Packing(p.s, np.c_[x, y, a])


def test_hash_invariant():
    p = read(REC)
    h = p.canon_hash()
    rng = np.random.default_rng(0)
    for sw in (0, 1):
        for fx in (0, 1):
            for fy in (0, 1):
                q = sym(p, sw, fx, fy)
                q = Packing(q.s, q.sq[rng.permutation(q.n)])
                assert q.canon_hash() == h
                assert abs(q.max_pen() - p.max_pen()) < 1e-12


def test_roundtrip_formats():
    p = read(REC)
    for fmt in ('ours', 'couzo', 'ellsworth', 'json'):
        q = parse(p.text(fmt))
        assert q.n == p.n and abs(q.s - p.s) < 1e-12 and np.allclose(q.sq, p.sq, atol=1e-9), fmt


if __name__ == '__main__':
    test_hash_invariant(); test_roundtrip_formats(); print('ok')


def test_rect_container_pen_and_hash():
    from pk.packing import Packing
    two = [[0.5, 0.5, 0], [1.5, 0.5, 0]]
    assert Packing(1.0, two, dx=1.0).max_pen() == 0.0          # 2 x 1 box holds both
    assert Packing(1.0, two).max_pen() > 0.9                    # 1 x 1 box does not
    a = Packing(1.0, two, dx=1.0)
    b = Packing(1.0, [[1.5, 0.5, 0], [0.5, 0.5, 0]], dx=1.0)   # relabelled
    assert a.canon_hash() == b.canon_hash()
    assert a.canon_hash() != Packing(1.0, two, dx=1.5).canon_hash()
