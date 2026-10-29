"""ingest --fresh-later / --fresh-pending state (10-28): a deferred window stays pending until a fresh pass succeeds.
Fake runner (nothing fetched or run), temp state dir, empty scratch store.
  python3 pk/tests/test_ingest_later.py"""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pk import ingest as I
from pk.store import Store


def main():
    d = tempfile.mkdtemp(); root = os.path.join(d, 'ingest'); sp = os.path.join(d, 's.sqlite'); Store(sp).close()
    calls = []
    ok = lambda cmd, f: calls.append(cmd[1].split('/')[-1] + ' ' + cmd[2]) or 0
    q = lambda *a, **k: None
    I.run(store=sp, apply=True, since='2026-10-01T00:00:00', later=True, root=root, runner=ok, log=q)
    S = I.load_state(root)
    assert S['fresh_pending'] == '2026-10-01T00:00:00' and S['last'], S
    assert not any('fresh.py' in c for c in calls), calls                  # syncs only
    I.run(store=sp, apply=True, later=True, root=root, runner=ok, log=q)   # a second deferral keeps the earliest time
    assert I.load_state(root)['fresh_pending'] == '2026-10-01T00:00:00'
    assert I.run_pending(store=sp, root=root, runner=lambda c, f: 1, log=q) == 1
    assert I.load_state(root)['fresh_pending'] == '2026-10-01T00:00:00'    # failed: still pending
    calls.clear()
    assert I.run_pending(store=sp, root=root, runner=ok, log=q) == 0
    assert 'fresh_pending' not in I.load_state(root) and any('fresh.py' in c for c in calls), calls
    assert I.run_pending(store=sp, root=root, runner=ok, log=q) is None   # nothing left
    print('test_ingest_later: ok')


if __name__ == '__main__':
    main()
