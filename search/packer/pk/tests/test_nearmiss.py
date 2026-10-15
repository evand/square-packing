"""Near-miss queue, finisher and ingest logic on tiny synthetic temp stores (no fq / harden / network / candidate.py):
python3 pk/tests/test_nearmiss.py"""
import json, os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from pk.packing import Packing
from pk.store import Store
from pk import nearmiss as NM, ingest

F10 = 3.9999                                      # synthetic frontier at n = 10 (below the grid 4)


def pk10(s):
    """10 axis squares on a 4 x 4 lattice scaled to side s (feasible for s >= 4), at most 3 per row / column"""
    return Packing(s, [((c + 0.5) * s / 4, (r + 0.5) * s / 4, 0.0) for r, c in CELLS])


CELLS = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (1, 3), (2, 0), (2, 3), (3, 0), (3, 3)]


def fr(n):
    return F10 if n == 10 else float(int(np.ceil(np.sqrt(n))))


def mkstore(d):
    st = Store(os.path.join(d, 's.sqlite'))
    ids = {}
    for key, s, src, grid in (('a', 4.0001, 'run:runs/x', 0), ('a2', 4.00010005, 'run:runs/x', 0),   # a2: copy of a
                              ('b', 4.0003, 'battery:t', 0), ('c', 4.0005, 'run:runs/y', 0),
                              ('pub', 4.0002, 'register', 0), ('fin', 4.00015, 'finish:polish', 0),
                              ('grid', 4.00025, 'run:runs/z', 1), ('far', 4.01, 'run:runs/x', 0)):
        ids[key], _ = st.add(pk10(s), src, pre=dict(grid=grid))
    return st, ids


def test_lattice():
    assert pk10(4.0).max_pen() <= 1e-12 and pk10(4.0).grid_lines() == 0


def test_select():
    C = [dict(packing=i, n=10, s=4 + g, gap=g) for i, g in ((1, 3e-4), (2, 1e-4), (3, 1e-4 + 5e-8), (4, 5e-4))]
    out = NM.select(C, [])
    assert [c['packing'] for c in out] == [2, 1]                   # 3 = copy of 2; per n 2
    out = NM.select(C, [dict(n=10, s=4 + 2e-4, gap=2e-4), dict(n=10, s=4 + 4e-4, gap=4e-4)])
    assert [c['packing'] for c in out] == [2]                      # full n: only one better than every item gets in
    assert NM.select(C, [dict(n=11, s=1, gap=0)], per_n=1)[0]['packing'] == 2


def test_scan_enqueue():
    with tempfile.TemporaryDirectory() as d:
        st, ids = mkstore(d)
        c = NM.scan(st.db, fr, per_n=3)
        got = [x['packing'] for x in c]
        assert got == [ids['a'], ids['b'], ids['c']], (got, ids)  # public, finish:, grid, far, tie, copy excluded
        assert c[0]['origin'] == 'run:runs/x' and abs(c[0]['gap'] - (4.0001 - F10)) < 1e-12
        assert NM.scan(st.db, fr, per_n=3, include_public=True)[1]['packing'] == ids['pub']
        assert NM.scan(st.db, fr, per_n=3, min_gap=3e-4)[0]['packing'] == ids['b']   # floor: a (gap 2e-4) left out
        assert NM.enqueue(st, c[:2]) == 2 and NM.enqueue(st, c[:2]) == 0
        assert 'nearmiss' in st.get(ids['a']).meta['tags']
        assert [r['packing'] for r in NM.queue(st.db, 'queued')] == [ids['a'], ids['b']]
        assert NM.scan(st.db, fr, per_n=2) == []                   # n full, c is not better than a or b
        assert [x['packing'] for x in NM.scan(st.db, fr, per_n=3)] == [ids['c']]
        assert NM.scan(st.db, fr, ids=[]) == [] and NM.scan(st.db, fr, nmax=9) == []
        ro = NM.ReadOnly(os.path.join(d, 's.sqlite'))               # read-only view sees the queue
        assert len(NM.queue(ro.db)) == 2
        st.close()


def test_bar():
    it = dict(n=10, s=4.0, gap=-1e-4, frontier=4.0001)
    assert NM.bar(it, lambda n: 4.0) == 4.0001                    # the frontier is the hit itself
    assert NM.bar(it, lambda n: 3.99) == 3.99                     # someone beat it since
    assert NM.bar(dict(it, gap=2e-4), lambda n: 3.9998) == 3.9998


def fake_runner(delta, fail=False, calls=None):
    """stage runner: writes the input scaled to side s - delta (feasible for these lattices while s - delta >= 4)"""
    def run(inp, wd, n, target, procs):
        if calls is not None:
            calls.append(os.path.basename(wd))
        if fail:
            raise RuntimeError('boom')
        from pk.packing import read
        p = read(inp)
        q = pk10(p.s - delta)
        out = os.path.join(wd, 'out.txt'); q.write(out)
        return [(float(q.s), out)], 3
    return run


def test_finish():
    with tempfile.TemporaryDirectory() as d:
        st, ids = mkstore(d)
        NM.enqueue(st, NM.scan(st.db, fr, per_n=3))
        calls = []
        R = NM.finish(st, os.path.join(d, 'out'), limit=1, fnow=fr, log=lambda *a: None,
                      runners=dict(polish=fake_runner(2e-5, calls=calls), harden1=fake_runner(1e-5, calls=calls),
                                   harden=fake_runner(0.0, calls=calls)))
        r = R[0]
        assert r['packing'] == ids['a'] and r['err'] is None and calls == ['polish', 'harden1', 'harden'], (r, calls)
        assert abs(r['best_s'] - (4.0001 - 3e-5)) < 1e-12 and r['best_pid'] != ids['a']
        q = NM.queue(st.db)
        assert {x['packing']: x['status'] for x in q}[ids['a']] == 'done'
        row = st.db.execute('SELECT * FROM expansion WHERE packing=?', (ids['a'],)).fetchone()
        assert row['mechanism'] == 'finish' and row['nprop'] == 9 and abs(row['frontier_s'] - F10) < 1e-12
        assert [s['source'] for s in st.sightings(r['best_pid'])] == ['finish:harden1']
        assert st.get(r['best_pid']).meta['parent'] == ids['a']
        assert all(x['packing'] != r['best_pid'] for x in NM.scan(st.db, fr, per_n=9))   # outputs are not items
        # gate: harden skipped when the gap after harden1 is above nm_gate; failure -> failed
        calls.clear()
        R = NM.finish(st, os.path.join(d, 'out'), limit=1, fnow=fr, nm_gate=1e-4, log=lambda *a: None,
                      runners=dict(polish=fake_runner(0.0, calls=calls), harden1=fake_runner(0.0, calls=calls),
                                   harden=fake_runner(0.0, fail=True)))
        assert calls == ['polish', 'harden1'] and R[0]['stages'][-1]['skipped'].startswith('gap'), R
        R = NM.finish(st, os.path.join(d, 'out'), limit=1, fnow=fr, log=lambda *a: None,
                      runners=dict(polish=fake_runner(0.0, fail=True)), stages=['polish'])
        assert R[0]['err'] and {x['packing']: x['status'] for x in NM.queue(st.db)}[R[0]['packing']] == 'failed'
        # polish lands on the frontier (tie): later stages skipped
        NM.set_status(st, ids['a'], 'queued'); calls.clear()
        R = NM.finish(st, os.path.join(d, 'out'), limit=1, fnow=lambda n: 4.00007 if n == 10 else fr(n),
                      log=lambda *a: None, runners=dict(polish=fake_runner(3e-5, calls=calls),
                                                        harden1=fake_runner(0, calls=calls)))
        assert calls == ['polish'] and R[0]['stages'][-1]['skipped'].startswith('tie'), R
        st.close()


def test_finish_candidate_dry():
    with tempfile.TemporaryDirectory() as d:
        st, ids = mkstore(d)
        NM.enqueue(st, NM.scan(st.db, fr, per_n=1))
        f = lambda n: 4.00005 if n == 10 else fr(n)                 # frontier above a after one stage: a hit
        calls = []
        R = NM.finish(st, os.path.join(d, 'out'), fnow=f, log=lambda *a: None,
                      runners=dict(polish=fake_runner(1e-4, calls=calls), harden1=fake_runner(0, calls=calls)))
        r = R[0]
        assert calls == ['polish'] and r['candidate'][1].endswith('candidate.py') and 'candidate_rc' not in r
        q = NM.queue(st.db)[0]
        assert q['status'] == 'candidate' and q['note'].startswith('candidate (dry run)')
        # dry-run finish on a read-only view writes nothing
        NM.set_status(st, ids['a'], 'queued')
        st.db.commit()
        n0 = st.db.execute('SELECT count(*) FROM expansion').fetchone()[0]
        R = NM.finish(NM.ReadOnly(os.path.join(d, 's.sqlite')), None, dry=True, fnow=f, log=lambda *a: None)
        assert R[0]['dry'] and st.db.execute('SELECT count(*) FROM expansion').fetchone()[0] == n0
        st.close()


def write_outputs(d, n=10):
    """a compose2-style dir: near.jsonl (n + 2 squares at S_min, trimmed to n) + a harden-style results.jsonl"""
    os.makedirs(d, exist_ok=True)
    big = pk10(4.0004).sq.tolist() + [[0.5, 3.7, 0.0], [1.5, 3.7, 0.0]]
    with open(os.path.join(d, 'near.jsonl'), 'w') as f:
        f.write(json.dumps(dict(n=n, A='a', B='b', S_min=4.0004, need=4.0004 - F10, sq=big)) + '\n')
        f.write(json.dumps(dict(n=n, A='a', B='c', S_min=4.1, need=0.1, sq=big)) + '\n')        # too far
    pk10(4.0006).write(os.path.join(d, 'pol.txt'))
    with open(os.path.join(d, 'results.jsonl'), 'w') as f:
        f.write(json.dumps(dict(n=n, ok=True, s=4.0006, out='pol.txt')) + '\n')
        f.write(json.dumps(dict(n=n, ok=False, err='x')) + '\n')


def test_outputs():
    with tempfile.TemporaryDirectory() as d:
        st = Store(os.path.join(d, 's.sqlite'))
        od = os.path.join(d, 'compose_lo'); write_outputs(od)
        R = NM.read_outputs(od, frontier=fr)
        assert sorted(k for k, *_ in R) == ['harden', 'near'] and all(p.n == 10 for *_, p in R)
        items = NM.from_outputs(st, [od], fr)
        assert len(items) == 2 and items[0]['origin'].startswith('compose2:'), items
        assert NM.enqueue(st, items) == 2
        st.close()


def test_ingest():
    with tempfile.TemporaryDirectory() as d:
        sp, root = os.path.join(d, 's.sqlite'), os.path.join(d, 'ingest')
        st, ids = mkstore(d); st.close()
        logs = []
        try:
            ingest.run(store=sp, root=root, log=logs.append)
            raise AssertionError('first ingest without --since must fail')
        except SystemExit:
            pass
        r = ingest.run(store=sp, root=root, since='2000-01-01T00:00:00', log=logs.append)
        assert r['apply'] is False and not os.path.exists(root)                 # dry run writes nothing
        assert [n for n, _ in r['plan']] == ['sync-register', 'sync-pending', 'watch', 'fresh']
        assert '--nmax' in r['plan'][-1][1] and r['plan'][-1][1][r['plan'][-1][1].index('--since') + 1] == '2000-01-01T00:00:00'
        assert any('store already holds 1 packings' in l for l in logs), logs   # the 'register' row
        seen = []

        def runner(cmd, logf):
            seen.append(cmd)
            if 'fresh.py' in cmd[1]:
                write_outputs(os.path.join(cmd[cmd.index('--out') + 1], 'compose_lo'))
            return 0
        import pk.nearmiss as nm
        orig = nm.frontier_fn
        nm.frontier_fn = lambda path=None: fr                     # synthetic frontier (the real one reads events)
        try:
            r = ingest.run(store=sp, root=root, since='2000-01-01T00:00:00', apply=True, runner=runner, log=logs.append)
            assert r['ok'] and r['nearmiss_fresh'] == 2, r
            assert r['nearmiss_scan'] >= 1, r                      # store rows added since (a, b, c; per n 2 incl. fresh)
            S = ingest.load_state(root)
            assert S['last'] == r['start'] and len(seen) == 4
            r2 = ingest.run(store=sp, root=root, apply=True, runner=lambda c, f: 1 if 'watch' in c else 0,
                            log=logs.append)
            assert r2['since'] == r['start'] and not r2['ok'] and ingest.load_state(root)['last'] == r['start']
        finally:
            nm.frontier_fn = orig




def test_reclaim_stale():
    with tempfile.TemporaryDirectory() as d:
        st, ids = mkstore(d)
        NM.enqueue(st, NM.scan(st.db, fr, per_n=3))
        dead = 2 ** 22 + 12345                     # above the default pid_max on most systems: not a live process
        while NM._pid_alive(dead):
            dead += 1
        NM.set_status(st, ids['a'], 'running', note=f'pid {dead}')
        NM.set_status(st, ids['b'], 'running', note=f'pid {os.getpid()}')     # alive: left alone
        NM.set_status(st, ids['c'], 'running')                                 # no pid (pre-10-15): reported only
        msgs = []
        assert NM.reclaim_stale(st, log=msgs.append) == 1 and len(msgs) == 1
        S = {x['packing']: x['status'] for x in NM.queue(st.db)}
        assert (S[ids['a']], S[ids['b']], S[ids['c']]) == ('queued', 'running', 'running'), S
        st.close()


if __name__ == '__main__':
    for k, f in list(globals().items()):
        if k.startswith('test_'):
            f()
    print('ok')
