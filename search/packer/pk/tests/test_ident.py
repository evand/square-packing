"""pk/ident: grouping, protection, apply, and the dedup-following hooks (core test faked: deviation = max centre move)."""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from pk import ident
from pk.store import Store
from pk.packing import Packing


def test_plan_apply(monkeypatch):
    monkeypatch.setattr(ident, 'core_mask', lambda a: (a[0], np.ones(3, bool)))
    monkeypatch.setattr(ident, 'core_dev', lambda a, b: float(np.abs(a[1][:, :2] - b[1][:, :2]).max()))
    monkeypatch.setattr(ident, 'census_protect', lambda n: set())
    from pk import policy
    monkeypatch.setattr(policy, 'protected_ids', lambda st, cfg: set())
    with tempfile.TemporaryDirectory() as d:
        st = Store(os.path.join(d, 's.sqlite'))
        base = np.array([[0.7, 0.7, 20.0], [2.0, 0.7, 20.0], [0.7, 2.0, 20.0]])   # tilted: no grid line
        def add(ds, dx, src, status='polished'):
            q = base.copy(); q[0, 0] += dx
            return st.add(Packing(2.5 + ds, q), src, status=status, check=False)[0]
        a = add(0, 0, 'A')                       # representative (lowest side)
        b = add(3e-12, 2e-5, 'B')                # copy: side within 1e-11, core within 1e-4
        c = add(6e-12, 3e-3, 'C')                # same side, core far: distinct
        e = add(5e-11, 1e-5, 'E')                # side outside the window: distinct
        f = add(2e-12, 1e-5, 'F', status='not-min')  # not a candidate
        child = st.add(Packing(2.6, base), 'K', parent=b, check=False)[0]
        P = ident.plan(st, 3, procs=1, log=lambda *x: None)
        assert {k: [o for o, *_ in v] for k, v in P['groups'].items()} == {a: [b]}, P['groups']
        ident.apply(st, P, log=lambda *x: None)
        assert st.db.execute('SELECT COUNT(*) FROM packing WHERE id=?', (b,)).fetchone()[0] == 0
        assert {r[0] for r in st.db.execute('SELECT source FROM sighting WHERE packing=?', (a,))} == {'A', 'B'}
        assert st.db.execute('SELECT parent FROM packing WHERE id=?', (child,)).fetchone()[0] == a
        assert st.get(b).meta['id'] == a and st.get(b).meta['dedup_of'] == b
        q = base.copy(); q[0, 0] += 2e-5
        pid, new = Store(os.path.join(d, 's.sqlite')).add(Packing(2.5 + 3e-12, q), 'again', check=False)
        assert (pid, new) == (a, False)          # the removed copy's hash resolves to its representative
        for x in (c, e, f):
            assert st.get(x).meta['id'] == x


def test_guard(monkeypatch):
    """two known exact minima within GUARD of the rows' side: no merge"""
    monkeypatch.setattr(ident, 'core_mask', lambda a: (a[0], np.ones(3, bool)))
    monkeypatch.setattr(ident, 'census_protect', lambda n: set())
    from pk import policy
    monkeypatch.setattr(policy, 'protected_ids', lambda st, cfg: set())
    with tempfile.TemporaryDirectory() as d:
        st = Store(os.path.join(d, 's.sqlite'))
        base = np.array([[0.7, 0.7, 20.0], [2.0, 0.7, 20.0], [0.7, 2.0, 20.0]])
        st.add(Packing(2.5, base), 'A', status='polished', check=False)
        q = base.copy(); q[0, 0] += 2e-5
        st.add(Packing(2.5 + 3e-12, q), 'B', status='polished', check=False)
        for S in ('2.5000000000001', '2.5000000000003'):
            st.db.execute('INSERT INTO minimum(n, key, S, cls) VALUES (3, ?, ?, ?)', (S, S, 'certified'))
        P = ident.plan(st, 3, procs=1, log=lambda *x: None)
        assert P['groups'] == {} and P['stats']['refused'] == {'near two known minima (guard band)': 1}, P['stats']
