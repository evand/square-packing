"""pk.distinct: matching over box symmetries, move counts, verdicts (synthetic work dirs; no exactsolve).

  python pk/tests/test_distinct.py
"""
import json, os, random, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from pk import distinct as D

S = 3.2071067811865475


def config(shift=0.0):
    """8 squares in a 3.207 box: 7 axis-parallel + one 45-degree square; `shift` moves one held square"""
    sq = [[0.5, 0.5, 0.0], [1.5, 0.5, 0.0], [2.707, 0.5, 0.0], [0.5, 1.5, 0.0], [0.5, 2.707, 0.0],
          [1.6 + shift, 1.6, 45.0], [2.707, 2.707, 0.0], [2.707, 1.55, 12.5]]
    return np.array(sq)


EQS = [['W', 0, 2, 'B'], ['C', 0, 1, 1, 3], ['C', 1, 2, 3, 0]]     # exactsolve format: ['C', i, corner, j, side]


def workdir(root, name, sq, free=(), S_exact=None, side=S, order=None):
    """order[r] = config() row of row r (from sym_perm); contact equations are relabelled to match"""
    w = os.path.join(root, name); os.makedirs(w)
    S_exact = S_exact or f'{side:.40f}'
    with open(os.path.join(w, 'cand.exact.txt'), 'w') as f:
        f.write(f'{len(sq)} {side!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(a)!r}\n' for x, y, a in sq))
    json.dump(dict(S_exact=S_exact, S_cert_decimal=S_exact, second_order=dict(status='local minimum of S, strict', n_zero=3),
                   n_vv=1, vv_milp_dS=0.0), open(os.path.join(w, 'cand.json'), 'w'))
    new = {o: r for r, o in enumerate(order)} if order is not None else {i: i for i in range(len(sq))}
    eqs = [[e[0], new[e[1]]] + e[2:] if e[0] == 'W' else [e[0], new[e[1]], e[2], new[e[3]], e[4]] for e in EQS]
    json.dump(dict(equations=eqs, free=list(free)), open(os.path.join(w, 'cand.contacts.json'), 'w'))
    json.dump(dict(S_exact=S_exact, groups_after=3, orientation_g=0,
                   packing=dict(S=S_exact, sq=[[repr(float(v)) for v in r] for r in sq])),
              open(os.path.join(w, 'regularized.json'), 'w'))
    return w


def sym_perm(sq, k, seed=0, with_order=False):
    """apply box symmetry k and shuffle the rows (with_order: also return the row order)"""
    q = D._sym(sq, S, k)
    idx = list(range(len(q))); random.Random(seed).shuffle(idx)
    return (q[idx], idx) if with_order else q[idx]


def test_match():
    A = config()
    for k in range(8):
        B = sym_perm(A, k, seed=k)
        md, kk, perm, Bk = D.best_match(A, S, B, S)
        assert md < 1e-9, (k, md)
        mv, _ = D.moves(A, Bk[perm])
        assert mv.max() < 1e-9


def test_tilt():
    assert abs(D.tilt(89.9) - 0.1) < 1e-9 and abs(D.tilt(-44.0) - 44.0) < 1e-9 and D.tilt(90.0) < 1e-9


def test_verdicts():
    with tempfile.TemporaryDirectory() as t:
        a = workdir(t, 'a', config(), free=[7])
        q, order = sym_perm(config(), 5, 1, with_order=True)
        same = workdir(t, 'same', q, free=[order.index(7)], order=order)            # same rows, reoriented and shuffled
        C = D.compare(a, same)
        assert C['verdict'] == 'same minimum', D.summary_lines(C)
        moved = workdir(t, 'moved', config(shift=0.003), free=[7])     # same side; moves are descriptive only (10-19)
        K = json.load(open(os.path.join(moved, 'cand.contacts.json')))
        K['equations'][1][3] = 2; json.dump(K, open(os.path.join(moved, 'cand.contacts.json'), 'w'))
        C = D.compare(a, moved)
        assert C['verdict'] == 'same side, different contact graph' and C['difference']['moved_gt_1e4'] == 1, D.summary_lines(C)
        assert C['difference']['shared_angles'] == ['45.0000'], C['difference']
        rat = workdir(t, 'rat', np.vstack([config()[:7], [[2.6, 1.6, 30.0]]]), free=[7])  # only the rattler differs
        C = D.compare(a, rat)
        assert C['verdict'] == 'same minimum', D.summary_lines(C)
        far = workdir(t, 'far', config(), free=[7], side=S + 1e-8)
        C = D.compare(a, far)
        assert C['verdict'] == 'distinct minimum', D.summary_lines(C)
        assert '| exact side S* |' in D.table_md(C)


if __name__ == '__main__':
    test_match(); test_tilt(); test_verdicts()
    print('test_distinct: ok')
