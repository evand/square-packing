"""Tests for pk/units.py claim / record (temp database; no store, no compute).  python pk/tests/test_units.py"""
import os, random, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pk import units as U
from pk.goal import Posterior

ARMS = [dict(id='10:a', n=10, chain='x | relax.y'), dict(id='10:b', n=10, chain='x | relax.y'), dict(id='20:a', n=20, chain='x | relax.y')]
VAL = dict(frontier=1e5, bar=100.0)


def _db(d, budget=1.0, unit=60.0):
    db = U.connect(os.path.join(d, 'u.sqlite'))
    db.execute("INSERT INTO goal VALUES ('g', 'x', ?, ?, 0, 1e-3, 'active', '2026-01-01', 0)", (budget, unit))
    for n, s in ((10, 4.0), (20, 5.0)):
        db.execute("INSERT INTO nstate(goal, n, sbest0, sbest, f0) VALUES ('g', ?, ?, ?, ?)", (n, s, s, s))
    return db


def _claim(db, now=100.0, rng=None):
    return U.claim(db, 1, rng or random.Random(0), arms_of=lambda g: ARMS, posterior=Posterior, now=now)


def test_arm_reference_polish_lineage():
    with tempfile.TemporaryDirectory() as d:
        db = _db(d)
        u = _claim(db)
        assert u['kind'] == 'construct'
        n = u['n']; s0 = db.execute("SELECT sbest FROM nstate WHERE n=?", (n,)).fetchone()[0]
        r = U.record(db, u, dict(jbest=s0 - 0.01, cands=[('polish', 7, s0 - 0.01), ('polish', 8, s0 + 0.5)]), 50.0, 55.0, VAL, {}, now=110)
        assert abs(r['store_gain'] - 0.01) < 1e-12
        u2 = _claim(db, 120)                                       # a candidate at n, no reference yet
        assert (u2['kind'], u2['n']) == ('reference', n)
        u3 = _claim(db, 121)                                       # reference running: nothing to polish yet -> arm
        assert u3['kind'] == 'construct'
        U.record(db, u2, dict(f0p=s0 - 0.005), 10.0, 10.0, VAL, {}, now=130)
        u4 = _claim(db, 140)                                       # the far candidate is skipped (margin 1e-3)
        assert (u4['kind'], u4['parent'], u4['arm']) == ('polish', u['id'], u['arm'])
        assert db.execute("SELECT status FROM cand WHERE packing=8").fetchone()[0] == 'skipped'
        r = U.record(db, u4, dict(fbest=s0 - 0.012, jbest=s0 - 0.012), 30.0, 30.0, VAL, {}, now=150)
        assert abs(r['gain'] - 0.007) < 1e-12 and abs(r['store_gain'] - 0.002) < 1e-12
        w, c = db.execute("SELECT SUM(reward), SUM(cpu) FROM unit WHERE arm=? AND status='done'", (u['arm'],)).fetchone()
        assert abs(w - 1e5 * 0.019) < 1e-6 and c == 80.0               # the arm pays for and earns its polish


def test_budget_and_lease():
    with tempfile.TemporaryDirectory() as d:
        db = _db(d, budget=3 * 60 / 3600)                          # room for 3 units of 60 s
        us = [_claim(db) for _ in range(3)]
        assert all(u and u['kind'] == 'construct' for u in us) and _claim(db) is None
        lost = U.expire(db, 100 + 601)
        assert len(lost) == 3
        assert _claim(db, 702)['kind'] == 'construct'                    # lost units cost nothing recorded: budget back


def test_bar_retires():
    with tempfile.TemporaryDirectory() as d:
        db = _db(d)
        u = _claim(db)
        r = U.record(db, u, dict(jbest=1.0), 1.0, 1.0, VAL, {u['n']: 2.0}, now=110)
        assert r['bar'] == 1
        for _ in range(5):
            assert _claim(db, 120)['n'] != u['n']



def test_walk_candidates_lineage():
    with tempfile.TemporaryDirectory() as d:
        db = _db(d)
        u = _claim(db)
        n = u['n']; s0 = db.execute("SELECT sbest FROM nstate WHERE n=?", (n,)).fetchone()[0]
        U.record(db, u, dict(jbest=s0 + 0.01, cands=[('walk', 5, s0 + 0.01), ('walk', 6, s0 + 0.2)]), 50.0, 50.0, VAL, {}, now=110)
        w = _claim(db, 120)                                        # the near start is walked; the far one is dropped
        assert (w['kind'], w['start'], w['arm'], w['parent']) == ('walk', 5, u['arm'], u['id'])
        assert db.execute("SELECT status FROM cand WHERE packing=6").fetchone()[0] == 'skipped'
        assert _claim(db, 121)['kind'] == 'construct'              # nothing else pending
        U.record(db, w, dict(jbest=s0 - 0.003), 60.0, 60.0, VAL, {}, now=200)
        r, c = db.execute("SELECT SUM(reward), SUM(cpu) FROM unit WHERE arm=? AND status='done'", (u['arm'],)).fetchone()
        assert abs(r - 300.0) < 1e-6 and c == 110.0                # the construct arm pays for and earns its walk


def test_walk_arm_kind():
    with tempfile.TemporaryDirectory() as d:
        db = _db(d)
        u = U.claim(db, 1, random.Random(0), arms_of=lambda g: [dict(id='10:walk', n=10, chain=None)], posterior=Posterior, now=100)
        assert (u['kind'], u.get('start')) == ('walk', None)


if __name__ == '__main__':
    test_arm_reference_polish_lineage(); test_budget_and_lease(); test_bar_retires(); test_walk_candidates_lineage(); test_walk_arm_kind(); print('ok')
