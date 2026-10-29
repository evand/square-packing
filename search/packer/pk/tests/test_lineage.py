"""pk/lineage.py walk rules (10-28) on a scratch store: packing.parent, move-log and file-name links, same-side copies,
cross-n flags, and credit going to the nearest public ancestor found by someone else (our own posted claims don't count).
  python3 pk/tests/test_lineage.py"""
import os, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pk.store import Store
from pk.packing import Packing
from pk import lineage as LG


def grid(n, s, nudge=0.0):
    """n upright unit squares in rows, spread over side s (feasible when they fit); nudge moves the last one (a copy
    at the same side with different coordinates, as a re-imported start file differs by float noise)"""
    k = int(np.ceil(np.sqrt(n))); c = np.linspace(0.5, s - 0.5, k)
    sq = np.array([(c[i % k], c[i // k], 0.0) for i in range(n)]); sq[-1, 0] -= nudge
    return Packing(s, sq)


def main():
    d = tempfile.mkdtemp()
    st = Store(os.path.join(d, 's.sqlite'))
    add = lambda n, s, src, finder=None, ref=None, parent=None, nudge=0.0: st.add(grid(n, s, nudge), src, finder, ref, parent=parent)[0]
    # someone else's public packing at n = 8 (the root), our copy of it (same side, same n), a cross-n child at 9 via
    # packing.parent, a move-log child, a file-name child, then our posted claim at the tip
    pub = add(8, 3.30, 'register', 'Someone Else')
    copy = add(8, 3.30, 'run:runs/x', 'ours', nudge=0.05)
    a = add(9, 3.40, 'run:runs/x', 'ours', parent=copy)
    b = add(9, 3.39, 'run:runs/y', 'ours')
    st.db.execute("INSERT INTO expansion(packing, n, mechanism, cpu_s, added, child) VALUES (?, 9, 'move.test', 1, '2026', ?)", (a, b))
    c = add(9, 3.38, 'run:runs/z', 'ours', ref=f'runs/z/n9_p{b}_k1.txt')
    tip = add(9, 3.37, 'pending:#999', 'Evan Daniel (evand)')
    st.db.execute('UPDATE packing SET parent=? WHERE id=?', (c, tip))
    st.db.commit()

    steps = LG.walk(st.db, tip)
    ids = [s['id'] for s in steps]
    assert ids == [tip, c, b, a, copy, pub], ids
    hows = [s['how'] for s in steps]
    assert hows == ['parent', 'file-name', 'move-log', 'parent', 'same-side', None], hows
    assert 'cross-n' in steps[3]['flags'], steps[3]          # a (n = 9) <- copy (n = 8)
    assert 'root:public' in steps[-1]['flags']
    near, pubs = LG.credit(steps)
    assert near['id'] == pub and [s['id'] for s in pubs] == [pub], (near, pubs)   # tip is our own claim: not credit

    # an unrecorded root of ours: no credit start, flagged incomplete
    lone = add(10, 3.60, 'run:runs/gone', 'ours')
    s2 = LG.walk(st.db, lone)
    assert 'root:unrecorded' in s2[-1]['flags'] and 'run-dir-missing' in s2[-1]['flags'] and LG.credit(s2)[0] is None
    print('test_lineage: ok')


if __name__ == '__main__':
    main()
