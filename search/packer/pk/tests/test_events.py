"""pk/events.py without search: the grid-only target rule, options (relax first, once-moves, window, eligibility),
claim balance across conditions, record (tries / wins, new parents, condition bests, gain), lease expiry.

  python3 pk/tests/test_events.py
"""
import math, os, random, sys, tempfile, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pk import events as E

GRID = lambda n: float(E.kgrid(n))


def test_eligible():
    # row 17: sub-17 known down to c = 16 (n = 273); row 16 has c = 15 (n = 241); row 18 down to c = 16 (308)
    side = {n: GRID(n) - 0.01 for n in range(226, 274)}          # rows 16 / 17: every n with c >= 16 (row 17) sub-grid
    side.update({n: GRID(n) - 0.01 for n in range(290, 309)})    # row 18: c >= 16
    side[241] = 15.99                                            # row 16, c = 15
    assert E.eligible(273, side)                                 # sub-grid
    assert E.eligible(274, side)                                 # first grid-only of row 17, row 16 has c = 15
    assert not E.eligible(275, side)                             # deeper
    assert not E.eligible(309, side)                             # first grid-only of row 18, but row 17 lacks c = 15
    side[274] = 16.999
    assert E.eligible(309, side) and E.eligible(275, side)       # 274 fell: next in row 17 and the one below open
    del side[242]
    assert not E.eligible(275, side)                             # row 16 lacks c = 14 (n = 256 - 14 = 242)
    print('eligible ok')


def campaign_file(d):
    p = os.path.join(d, 'tcamp.yaml')
    open(p, 'w').write("""goal: tcamp
budget_cpu_h: 1
predictions: [test]
events:
  ns: [100, 101]
  moves:
    relax: {chain: 'relax.round.05c.p5', once: true, relax: true}
    walk: {chain: 'perturb.kick | relax.round.05c.p5'}
    add: {chain: 'grow.add | relax.round.05c.p5', dn: 1}
  policies: [rule, uniform]
  window: 1.0e-3
  unit_cpu: 10
""")
    return p


def seed_db(db, path):
    cfg = E.load(path)
    t = time.time()
    db.execute("INSERT INTO campaign VALUES ('tcamp', ?, 1, 10, 'active', ?, ?)", (path, t, t))
    for cond in cfg['policies']:
        for pid, n, s, rel in [(1, 100, 9.95, 1), (2, 100, 9.95 + 0.5e-3 * 10, 0), (3, 100, 9.95 + 2e-3 * 10, 1), (4, 101, 10.6, 1)]:
            db.execute('INSERT INTO parent VALUES (?,?,?,?,?,?,?,?,?,NULL)', ('tcamp', cond, pid, n, s, 'V', t - 3600, 0, rel))
        for n, s in [(100, 9.95), (101, 10.6)]:
            db.execute('INSERT INTO nbest VALUES (?,?,?,?,?,?)', ('tcamp', cond, n, s, s, t - 7200))
    for n in E.row_ns([100, 101]):
        db.execute('INSERT INTO known VALUES (?,?,?)', ('tcamp', n, GRID(n) - 0.05 if n <= 101 else GRID(n)))
    return cfg


def test_options_claim_record():
    d = tempfile.mkdtemp()
    db = E.connect(os.path.join(d, 'ev.sqlite'))
    cfg = seed_db(db, campaign_file(d))
    O = E.options(db, cfg, 'rule', time.time())
    have = {(o['pid'], o['move']) for o in O}
    assert (2, 'relax') in have and (2, 'walk') not in have      # unrelaxed: relax only
    assert (1, 'walk') in have and (1, 'relax') not in have      # relaxed: everything but relax
    assert (1, 'add') in have                                     # 100 -> 101 (target in ns, sub-grid)
    assert not any(o['pid'] == 3 for o in O)                      # outside the window
    assert all(abs(o['n_age_h'] - 2.0) < 0.01 for o in O)          # the best at n arrived 2 h ago
    assert (4, 'add') not in have                                 # 102 is not a target
    rng = random.Random(1)
    u1 = E.claim(db, 1, rng); u2 = E.claim(db, 1, rng)
    assert {u1['cond'], u2['cond']} == {'rule', 'uniform'}       # equal CPU: the second claim goes to the other one
    # record a unit with one win and a child at n 100 below the condition's best
    u = u1
    E.record(db, u, dict(status='done', apps=[dict(s=9.949, cpu=1.0, win=1), dict(s=None, cpu=1.0, win=0)],
                         children=[(99, u['n'], 9.949, 'V')]), 2.0, 2.0)
    t, w = db.execute('SELECT tries, wins FROM ptry WHERE campaign=? AND cond=? AND pid=? AND move=?',
                      ('tcamp', u['cond'], u['pid'], u['move'])).fetchone()
    assert (t, w) == (2, 1), (t, w)
    assert db.execute('SELECT fresh FROM parent WHERE campaign=? AND cond=? AND pid=99', ('tcamp', u['cond'])).fetchone()[0] == 1
    g = db.execute('SELECT gain FROM unit WHERE id=?', (u['id'],)).fetchone()[0]
    if u['n'] == 100:
        assert abs(g - 0.001 / 10) < 1e-12, g
    # the other condition never sees the child
    other = 'uniform' if u['cond'] == 'rule' else 'rule'
    assert not db.execute('SELECT 1 FROM parent WHERE campaign=? AND cond=? AND pid=99', ('tcamp', other)).fetchone()
    # once-move: after a relax claim + record, the parent is relaxed and offers no relax
    db.execute("UPDATE parent SET relaxed=1 WHERE pid=2")
    assert not any(o['pid'] == 2 and o['move'] == 'relax' for o in E.options(db, cfg, 'rule', time.time()))
    # lease expiry
    E.expire(db, time.time() + 1e6)
    assert db.execute("SELECT status FROM unit WHERE id=?", (u2['id'],)).fetchone()[0] == 'lost'
    print('options / claim / record ok')


def test_dn_and_sources():
    from pk import stages as S
    if 'gen.testsrc' not in S.STAGES:
        @S.stage('gen.testsrc', 'gen')
        def _g(st, rng, ctx):
            raise RuntimeError('not run here')
        _g.version = lambda ctx=None: 'v1'
        _g.variants = 2
    d = tempfile.mkdtemp(); p = os.path.join(d, 'tsrc.yaml')
    open(p, 'w').write("""goal: tsrc
budget_cpu_h: 1
predictions: [test]
events:
  ns: [119]
  moves:
    lstep: {chain: 'relax.round.05c.p5', dn: '2*k-1'}
    src: {chain: 'gen.testsrc | relax.round.05c.p5', source: true}
  policies: [rule]
""")
    cfg = E.load(p)
    assert E.target(cfg, 'lstep', 100) == 119 and E.target(cfg, 'lstep', 101) == 122   # k = 10, 11
    assert E.parent_ns(cfg) == [100]
    try:
        E._dn_fn('__import__("os")'); assert False
    except ValueError:
        pass
    db = E.connect(os.path.join(d, 'ev.sqlite'))
    t = time.time()
    db.execute("INSERT INTO campaign VALUES ('tsrc', ?, 1, 10, 'active', ?, ?)", (p, t, t))
    for n in E.row_ns([119]):
        db.execute('INSERT INTO known VALUES (?,?,?)', ('tsrc', n, GRID(n) - 0.01))
    rng = random.Random(0)
    a, b = E.claim(db, 1, rng), E.claim(db, 1, rng)
    assert a['cond'] == b['cond'] == E.SOURCE and {a['variant'], b['variant']} == {0, 1}
    c = E.claim(db, 1, rng)
    assert c is None or c['cond'] != E.SOURCE                     # both variants built at this version
    print('dn / sources ok')


def test_open_space_and_chosen_params():
    sp = E.open_space("relax.round(r0=('log', 0.02, 0.3), steps=6)")
    assert sp == {'relax.round': {'r0': ('log', 0.02, 0.3)}}, sp          # fixed values and scalar defaults are not open
    assert E.open_space('grow.add | relax.round') == {'grow.add': {'f': ('uni', 1.0, 1.04), 'K': ('choice', [4, 16, 64])}}
    d = tempfile.mkdtemp()
    db = E.connect(os.path.join(d, 'ev.sqlite'))
    cfg = seed_db(db, campaign_file(d))
    db.execute("UPDATE parent SET relaxed=1")
    E.SCORERS['picker'] = lambda opts, ctx, rng: ([1.0] * len(opts), [{'grow.add': {'K': 64}}] * len(opts))
    cfg['policies'][:] = ['picker']                               # the cached config: one condition using it
    db.execute("INSERT INTO nbest SELECT campaign, 'picker', n, s0, s, t FROM nbest WHERE cond='rule'")
    db.execute("INSERT INTO parent SELECT campaign, 'picker', pid, n, s, fam, t_in, fresh, relaxed, unit FROM parent WHERE cond='rule'")
    u = E.claim(db, 1, random.Random(0))
    assert u['params'] == {'grow.add': {'K': 64}}, u['params']
    cfg['policies'][:] = ['rule', 'uniform']
    print('open space / chosen params ok')


def test_batch_queue():
    """batch claims: one fill queues several distinct parents; pops need no rescoring; a once-move tried meanwhile is
    skipped; a new condition best at n drops that n's queued entries"""
    d = tempfile.mkdtemp()
    db = E.connect(os.path.join(d, 'ev.sqlite'))
    cfg = seed_db(db, campaign_file(d))
    db.execute("UPDATE parent SET relaxed=1")
    for pid, n, s in [(5, 100, 9.951), (6, 100, 9.952), (7, 100, 9.953)] + [(p, 101, 10.6 + 1e-4 * p) for p in range(20, 40)]:
        for cond in cfg['policies']:
            db.execute('INSERT INTO parent VALUES (?,?,?,?,?,?,?,?,?,NULL)', ('tcamp', cond, pid, n, s, 'V', time.time(), 0, 1))
    calls = []
    def counting(opts, ctx, rng):
        calls.append(ctx['draws'])
        return np.array([[1.0 + 0.01 * o['pid'] + 0.001 * k for k in range(ctx['draws'])] for o in opts])
    E.SCORERS['counting'] = counting
    cfg['policies'][:] = ['counting']; cfg['batch'] = 16; cfg['draws'] = 3; cfg['batch_per_n'] = 16; cfg['window'] = 1.0
    db.execute("INSERT INTO nbest SELECT campaign, 'counting', n, s0, s, t FROM nbest WHERE cond='rule'")
    db.execute("INSERT INTO parent SELECT campaign, 'counting', pid, n, s, fam, t_in, fresh, relaxed, unit FROM parent WHERE cond='rule'")
    rng = random.Random(0)
    u1 = E.claim(db, 1, rng)
    assert calls == [3], calls                                    # one scorer call, all draws at once
    q = db.execute("SELECT pid FROM queue WHERE campaign='tcamp' AND cond='counting'").fetchall()
    pids = [u1['pid']] + [r[0] for r in q]
    assert len(pids) == len(set(pids)) == 16, pids                # one entry per parent
    u2 = E.claim(db, 1, rng)
    assert u2['pid'] != u1['pid'] and len(calls) == 1             # popped, not rescored
    E.record(db, u2, dict(status='done', apps=[dict(s=9.9, cpu=1.0, win=1)], children=[(98, 100, 9.9, 'V')]), 1.0, 1.0)
    assert not db.execute("SELECT 1 FROM queue WHERE cond='counting' AND (n0=100 OR n=100)").fetchone()   # dropped
    cfg['batch_gap'] = 1e6
    E.claim(db, 1, rng); assert len(calls) == 1                   # new best, but within batch_gap of the last fill
    cfg['batch_gap'] = 0.0
    E.claim(db, 1, rng); assert len(calls) == 2, calls            # new best since the last fill -> refilled
    E.claim(db, 1, rng); assert len(calls) == 2, calls            # flag cleared: no second refill
    u = E.claim(db, 1, rng)
    E.record(db, u, dict(status='done', apps=[dict(s=9.8, cpu=1.0, win=1)], children=[(97, 100, 9.8, 'V')]), 1.0, 1.0)
    E.record(db, u, dict(status='done', apps=[dict(s=9.7, cpu=1.0, win=1)], children=[(96, 100, 9.7, 'V')]), 1.0, 1.0)
    E.claim(db, 1, rng); E.claim(db, 1, rng)
    assert len(calls) == 3, calls                                 # two bests in a row: one refill (coalesced)
    cfg['policies'][:] = ['rule', 'uniform']; cfg['batch'] = 64; cfg['draws'] = 8; cfg['batch_gap'] = 10.0
    cfg['batch_per_n'] = 8; cfg['window'] = 1.0e-3
    print('batch queue ok')


def test_mixed_policy():
    rng = random.Random(0)
    assert {E.policy('rule+uniform', rng) for _ in range(50)} == {'rule', 'uniform'}
    assert E.policy('rule-frozen') == 'rule'
    print('mixed policy ok')


def test_stop_rule():
    d = tempfile.mkdtemp()
    db = E.connect(os.path.join(d, 'ev.sqlite'))
    cfg = dict(name='ts', stop=dict(conds=['a', 'b'], p=0.95, min_cpu_h=1.0))
    def units(cond, n, gains, cpu):
        for i in range(n):
            db.execute("INSERT INTO unit(campaign, cond, status, cpu, gain) VALUES ('ts', ?, 'done', ?, ?)", (cond, cpu, 1e-6 if i < gains else 0.0))
    units('a', 100, 30, 36.0); units('b', 50, 2, 36.0)          # a: 1 CPU-h, 30 hits; b: 0.5 CPU-h -> not yet
    assert E.stop_state(db, cfg) is None
    units('b', 50, 1, 36.0)                                      # b: 1 CPU-h, 3 hits
    z = E.stop_state(db, cfg)
    assert z['decided'] and z['p'] > 0.99, z
    print('stop rule ok')


if __name__ == '__main__':
    test_open_space_and_chosen_params()
    test_batch_queue()
    test_mixed_policy()
    test_stop_rule()
    test_eligible()
    test_options_claim_record()
    test_dn_and_sources()
