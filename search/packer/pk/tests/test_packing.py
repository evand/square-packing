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
