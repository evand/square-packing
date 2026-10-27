"""The unit loop (DESIGN.md "Next: the unit loop", step 2; 10-26): goals spend CPU in short single-threaded units,
claimed one at a time from a shared database by persistent workers.  No dispatcher, no quantum, no drain.

  python3 -m pk.units add GOAL [--unit-cpu 120] [--prio 0] [--polish-margin 1e-3]   register a goal (its YAML budget)
  python3 -m pk.units work --procs N [--max-idle S]   N workers: claim / run / record until stopped (or idle S seconds)
  python3 -m pk.units status | report GOAL | stop [GOAL]

Launch workers with jobs.py like any compute: `pk/jobs.py start units --procs N -- python3 -m pk.units work --procs N`.

Kinds (each one single-threaded task of at most ~`unit_cpu` CPU-s, on the goal's arms and the store):
  construct  a chain arm: build, screen, and fully relax the builds whose screen score is in the arm's best 1 / ratio
             (pk/screen.construct_serial; the cut comes from every earlier unit's screens at that arm), then ingest
             + converge.  Its best new starts become walk candidates, its best new wanted-family packings within
             `polish_margin` of the family best polish candidates.
  walk       explore.py from one start for the allowance, then ingest + converge (polish candidates as above).  The
             start is a construct unit's candidate (the walk is charged and credited to that arm: lineage), or, for a
             walk arm, the store packing in the donor window with the least walk CPU spent on it so far, near the
             family best first (store `expansion` rows, mechanism 'walk', record every walk).
  reference  the polished family best at n before the goal (goal.reference; clusters off: the 3 best), built once,
             only when a polish candidate at n needs it
  polish     one unit's polish candidates to a polish fixpoint (goal.polish_cached); charged and credited to its arm
  refresh    rebuild a generator's inputs (stage hook `refresh()`), claimed only when its `stale()` says so

Deterministic generators (stage hooks, set by the stage's module, 10-26).  A gen stage with `version(ctx)` (a key of
the inputs its output depends on besides its parameters and rng) and only finite ('choice') parameters is built once
per (version, parameter combination, variant), `variants` (default 1) rng draws each: a construct unit of such an arm
fully relaxes its not-yet-built combinations (no screen) and the arm is not offered again until the version changes.
`stale()` / `refresh()`: the inputs can be brought up to date; refresh is a unit kind, claimed when stale (no timers).

Claim (one BEGIN IMMEDIATE transaction): expire lost leases; pick a goal (priority, then least budget used); in it a
reference if one is needed, else the best pending polish, else the best pending walk, else an arm by Thompson sampling
(goal.Posterior) on reward per CPU-h from the unit rows.  Record (one transaction): gains against the per-n state
(store best, polished family best) read inside the transaction, so concurrent units never pay the same gain twice.
Accounting is SELECTs over the unit and worker rows (CPU per unit = getrusage of the worker, children included; claim
/ record wall per unit).  Requires `clusters: false` (the cluster census needs one serial counter).
"""
from __future__ import annotations
import argparse, json, math, os, random, resource, socket, sqlite3, subprocess, sys, tempfile, time, traceback

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = os.path.join(HERE, 'runs', 'units.sqlite')
SCHEMA = """
CREATE TABLE IF NOT EXISTS goal(name TEXT PRIMARY KEY, path TEXT, budget_cpu_h REAL, unit_cpu_s REAL, prio REAL DEFAULT 0,
  polish_margin REAL, status TEXT, t_goal TEXT, added REAL);
CREATE TABLE IF NOT EXISTS nstate(goal TEXT, n INT, sbest0 REAL, sbest REAL, f0 REAL, f0p REAL, fbest REAL,
  retired INT DEFAULT 0, PRIMARY KEY(goal, n));
CREATE TABLE IF NOT EXISTS unit(id INTEGER PRIMARY KEY, goal TEXT, kind TEXT, n INT, arm TEXT, variant TEXT, parent INT,
  start INT, worker INT, seed INT, status TEXT, lease_until REAL, started REAL, ended REAL, cpu REAL, wall REAL,
  claim_s REAL, record_s REAL, gain REAL DEFAULT 0, store_gain REAL DEFAULT 0, bar INT DEFAULT 0, reward REAL DEFAULT 0,
  jbest REAL, info TEXT);
CREATE INDEX IF NOT EXISTS unit_goal ON unit(goal, status, arm);
CREATE INDEX IF NOT EXISTS unit_kind ON unit(goal, kind, n, status);
CREATE TABLE IF NOT EXISTS cand(id INTEGER PRIMARY KEY, goal TEXT, n INT, kind TEXT, unit INT, packing INT, s REAL,
  status TEXT, taken_by INT);
CREATE INDEX IF NOT EXISTS cand_pending ON cand(goal, kind, status, n);
CREATE TABLE IF NOT EXISTS screen(arm TEXT, k INT, score REAL, unit INT);
CREATE INDEX IF NOT EXISTS screen_arm ON screen(arm, k);
CREATE TABLE IF NOT EXISTS worker(id INTEGER PRIMARY KEY, pid INT, host TEXT, job TEXT, started REAL, last_seen REAL,
  cpu REAL DEFAULT 0, unit_cpu REAL DEFAULT 0, idle_s REAL DEFAULT 0, claim_s REAL DEFAULT 0, record_s REAL DEFAULT 0,
  units INT DEFAULT 0, status TEXT);
CREATE TABLE IF NOT EXISTS ctl(key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS gen_done(arm TEXT, version TEXT, combo TEXT, variant INT, unit INT,
  PRIMARY KEY(arm, version, combo, variant));
"""
EXTRA = [('unit', 'stage', 'TEXT')]                       # columns added after a db was created
LEASE = dict(construct=3.0, walk=3.0, reference=20.0, polish=20.0, refresh=20.0)   # x unit_cpu_s (min 600 s): past it, a lost worker
EMPTY_RETRY = 1800                                        # an arm with no donors is retried after this many seconds
POLISH_TOP = 3                                            # polish candidates per unit (finalize's clusters-off default)
WALK_TOP = 2                                              # walk candidates per construct unit (its best new starts)
HIST = 400                                                # screen scores per arm the serial cut is taken from
REFRESH_GAP = 30                                          # s between refreshes of one stage (guards a hook that stays stale)
WALKED = 'walk'                                           # store expansion mechanism of a walk unit


def connect(path=DB):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    db = sqlite3.connect(path, timeout=120, isolation_level=None)   # autocommit; transactions are explicit
    db.execute('PRAGMA journal_mode=WAL'); db.execute('PRAGMA synchronous=NORMAL')
    db.executescript(SCHEMA)
    from . import queue as Q
    db.executescript(Q.SCHEMA)
    for tb, col, ty in EXTRA:
        if col not in {r[1] for r in db.execute(f'PRAGMA table_info({tb})')}:
            db.execute(f'ALTER TABLE {tb} ADD COLUMN {col} {ty}')
    return db


# ------------------------------------------------------------------ deterministic generators

_GEN = {}


def gen_stage(arm):
    """the arm's first stage function if it is a deterministic generator (has `version`), else None (memoized per chain)"""
    if not arm.get('chain'):
        return None
    if arm['chain'] not in _GEN:
        from . import stages as S
        try:
            name, fixed = S.parse(arm['chain'].split('|')[0])[0]
        except KeyError:                                          # not a registered stage: running the chain will say so
            name = None
        f = S.STAGES.get(name)
        _GEN[arm['chain']] = f if f is not None and hasattr(f, 'version') and combos(f, fixed) is not None else None
    return _GEN[arm['chain']]


def combos(f, fixed):
    """every distinct parameter combination of stage f given the fixed ones ([] if none vary); None if a parameter is
    not finite ('choice')"""
    import itertools
    keys, vals = [], []
    for k, v in f.space.items():
        if k in fixed:
            continue
        if not (isinstance(v, tuple) and v and v[0] == 'choice'):
            return None
        keys.append(k); vals.append(sorted(set(v[1]), key=repr))
    return [dict(zip(keys, c)) for c in itertools.product(*vals)]


def _todo(db, arm, f, version):
    """combinations x variants of arm not yet built at this version -> [(combo dict, variant)]"""
    from . import stages as S
    done = {(c, v) for c, v in db.execute('SELECT combo, variant FROM gen_done WHERE arm=? AND version=?', (arm['id'], version))}
    return [(c, v) for c in combos(f, S.parse(arm['chain'].split('|')[0])[0][1]) for v in range(int(getattr(f, 'variants', 1)))
            if (json.dumps(c, sort_keys=True), v) not in done]


def cpu():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


# ------------------------------------------------------------------ goals

_GOALS = {}


def load_goal(path):
    """goal.load, cached per (path, mtime); unit goals need clusters off"""
    from . import goal as G
    key = (path, os.path.getmtime(path))
    if key not in _GOALS:
        g = G.load(path)
        if g.clusters:
            raise ValueError(f'{g.goal}: the unit loop needs clusters: false')
        if g.queue:
            from . import queue as Q
            g.targets = Q.targets(g)
        _GOALS[key] = (g, G.plan(g))
    return _GOALS[key]


def store_ro():
    return sqlite3.connect(f'file:{os.path.join(HERE, "runs", "store.sqlite")}?mode=ro', uri=True, timeout=60)


def add_goal(db, spec, unit_cpu=120.0, prio=0.0, margin=1e-3):
    from . import goal as G, families as F
    path = G.resolve(spec); g, arms = load_goal(path)
    os.makedirs(g.dir, exist_ok=True)
    sf = os.path.join(g.dir, 'start.json')                # goal.reference reads the goal start from here
    if not os.path.exists(sf):
        json.dump(dict(t_goal=time.strftime('%Y-%m-%dT%H:%M:%S'), units=True), open(sf, 'w'))
    t_goal = json.load(open(sf))['t_goal']
    st = store_ro()
    try:
        ex = g.progress.get('exclude', [])
        rows = [(n, G.best_before(st, n, t_goal), F.best(st, n, G.kgrid(n) + 0.05, ex, before=t_goal)[0]) for n in g.targets]
        if g.queue:
            from . import queue as Q
            db.execute('BEGIN IMMEDIATE')
            if not db.execute('SELECT 1 FROM qbest WHERE goal=? LIMIT 1', (g.goal,)).fetchone():
                Q.setup(db, g, t_goal, st)
            db.execute('COMMIT')
    finally:
        st.close()
    db.execute('BEGIN IMMEDIATE')
    db.execute('INSERT OR REPLACE INTO goal VALUES (?,?,?,?,?,?,?,?,?)',
               (g.goal, path, g.budget_cpu_h, unit_cpu, prio, margin, 'active', t_goal, time.time()))
    for n, s, f0 in rows:
        db.execute('INSERT OR IGNORE INTO nstate(goal, n, sbest0, sbest, f0) VALUES (?,?,?,?,?)', (g.goal, n, s, s, f0))
        if g.retire and n in g.bars and s < g.bars[n]:
            db.execute('UPDATE nstate SET retired=1 WHERE goal=? AND n=?', (g.goal, n))
    db.execute('COMMIT')
    return g, len(arms)


# ------------------------------------------------------------------ claim

def expire(db, now):
    lost = [i for (i,) in db.execute("SELECT id FROM unit WHERE status='running' AND lease_until<?", (now,))]
    if lost:
        q = ','.join('?' * len(lost))
        db.execute(f"UPDATE unit SET status='lost' WHERE id IN ({q})", lost)
        db.execute(f"UPDATE cand SET status='pending', taken_by=NULL WHERE status='claimed' AND taken_by IN ({q})", lost)
    return lost


def spent_h(db, g):
    """CPU-h used (finished units) and committed (running units at their allowance)"""
    done = db.execute("SELECT COALESCE(SUM(cpu), 0) FROM unit WHERE goal=? AND status!='running'", (g['name'],)).fetchone()[0]
    run = db.execute("SELECT COUNT(*) FROM unit WHERE goal=? AND status='running'", (g['name'],)).fetchone()[0]
    return done / 3600, run * g['unit_cpu_s'] / 3600


def pick(db, g, arms, rng, now, posterior, window=0.05):
    """the next option in goal g -> dict(kind, n, arm, parent, cand) or None.  window: walk candidates more than this
    above the family best are dropped (the goal's donor window)"""
    name = g['name']
    from . import stages as S
    hooks = {}
    for a in arms:
        f = gen_stage(a)
        if f is not None and hasattr(f, 'stale'):
            hooks[f.stage_name] = f
    for sname, f in hooks.items():
        if not db.execute("SELECT 1 FROM unit WHERE kind='refresh' AND stage=? AND (status='running' OR ended>?)",
                          (sname, now - REFRESH_GAP)).fetchone() and f.stale():     # a hook that stays stale cannot spin
            return dict(kind='refresh', n=None, arm=None, parent=None, stage=sname)
    fam = "(SELECT COALESCE(fbest, f0, sbest) FROM nstate WHERE nstate.goal=cand.goal AND nstate.n=cand.n)"
    db.execute(f"UPDATE cand SET status='skipped' WHERE goal=? AND kind='polish' AND status='pending' AND s > {fam} + ?",
               (name, g['polish_margin']))
    db.execute(f"UPDATE cand SET status='skipped' WHERE goal=? AND kind='walk' AND status='pending' AND s > {fam} + ?",
               (name, window))
    db.execute("UPDATE cand SET status='skipped' WHERE goal=? AND status='pending' AND n IN "
               "(SELECT n FROM nstate WHERE goal=? AND retired=1)", (name, name))
    r = db.execute("SELECT DISTINCT c.n FROM cand c JOIN nstate s ON s.goal=c.goal AND s.n=c.n WHERE c.goal=? AND "
                   "c.kind='polish' AND c.status='pending' AND s.f0p IS NULL AND NOT EXISTS (SELECT 1 FROM unit u WHERE "
                   "u.goal=c.goal AND u.n=c.n AND u.kind='reference' AND u.status='running') LIMIT 1", (name,)).fetchone()
    if r:
        return dict(kind='reference', n=r[0], arm=None, parent=None)
    r = db.execute("SELECT c.unit, c.n, MIN(c.s - s.fbest) d FROM cand c JOIN nstate s ON s.goal=c.goal AND s.n=c.n "
                   "WHERE c.goal=? AND c.kind='polish' AND c.status='pending' AND s.f0p IS NOT NULL GROUP BY c.unit "
                   "ORDER BY d LIMIT 1", (name,)).fetchone()
    if r:
        arm = db.execute('SELECT arm FROM unit WHERE id=?', (r[0],)).fetchone()[0]
        return dict(kind='polish', n=r[1], arm=arm, parent=r[0])
    r = db.execute("SELECT c.id, c.unit, c.n, c.packing FROM cand c JOIN nstate s ON s.goal=c.goal AND s.n=c.n WHERE "
                   "c.goal=? AND c.kind='walk' AND c.status='pending' ORDER BY c.s - COALESCE(s.fbest, s.f0, s.sbest) LIMIT 1",
                   (name,)).fetchone()
    if r:
        arm = db.execute('SELECT arm FROM unit WHERE id=?', (r[1],)).fetchone()[0]
        return dict(kind='walk', n=r[2], arm=arm, parent=r[1], cand=r[0], start=r[3])
    stats = [dict(arm=a, reward=w, cpu_h=c / 3600) for a, w, c in db.execute(
        "SELECT arm, SUM(reward), SUM(cpu) FROM unit WHERE goal=? AND arm IS NOT NULL AND status NOT IN ('running', 'lost') "
        "GROUP BY arm", (name,))]
    pend = dict(db.execute("SELECT arm, COUNT(*) FROM unit WHERE goal=? AND status='running' AND arm IS NOT NULL GROUP BY arm",
                           (name,)))
    empty = {a for (a,) in db.execute("SELECT DISTINCT arm FROM unit WHERE goal=? AND status='empty' AND ended>?",
                                      (name, now - EMPTY_RETRY))}
    retired = {n for (n,) in db.execute('SELECT n FROM nstate WHERE goal=? AND retired=1', (name,))}
    elig = [a for a in arms if a['id'] not in empty and a['n'] not in retired]
    vers = {}
    def built(a):                                                 # a deterministic generator with nothing left to build
        f = gen_stage(a)
        if f is None:
            return False
        if pend.get(a['id']):
            return True
        v = vers.setdefault(f.stage_name, f.version())
        return not _todo(db, a, f, v)
    elig = [a for a in elig if not built(a)]
    if not elig:
        return None
    post = posterior(stats); uh = g['unit_cpu_s'] / 3600
    a = max(elig, key=lambda a: post.sample(a['id'], rng, pend.get(a['id'], 0) * uh))
    return dict(kind='construct' if a.get('chain') else 'walk', n=a['n'], arm=a['id'], parent=None)


def claim(db, wid, rng, arms_of=None, posterior=None, now=None):
    """one transaction: the next unit for worker wid -> unit dict (with its goal row) or None (nothing to do)"""
    from . import goal as G
    arms_of = arms_of or (lambda g: load_goal(g['path'])[1])
    posterior = posterior or G.Posterior
    now = time.time() if now is None else now
    db.execute('BEGIN IMMEDIATE')
    try:
        expire(db, now)
        cols = [c[1] for c in db.execute('PRAGMA table_info(goal)')]
        goals = [dict(zip(cols, r)) for r in db.execute("SELECT * FROM goal WHERE status='active'")]
        cand = []
        for g in goals:
            done, run = spent_h(db, g)
            if done + run + g['unit_cpu_s'] / 3600 > g['budget_cpu_h'] + 1e-9:
                if run == 0:
                    db.execute("UPDATE goal SET status='done' WHERE name=?", (g['name'],))
                continue
            cand.append((-g['prio'], (done + run) / max(g['budget_cpu_h'], 1e-9), g))
        for _, _, g in sorted(cand, key=lambda t: t[:2]):
            gl = load_goal(g['path'])[0]
            if gl.queue:
                from . import queue as Q
                o = Q.pick(db, g, gl, arms_of(g), rng, now)
            else:
                o = pick(db, g, arms_of(g), rng, now, posterior)
            if o is None:
                continue
            lease = now + max(600.0, LEASE[o['kind']] * g['unit_cpu_s'])
            seed = rng.randrange(1 << 30)
            uid = db.execute("INSERT INTO unit(goal, kind, n, arm, parent, start, stage, variant, worker, seed, status, "
                             "lease_until, started) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                             (g['name'], o['kind'], o['n'], o['arm'], o['parent'], o.get('start'), o.get('stage'),
                              o.get('variant'), wid, seed, 'running', lease, now)).lastrowid
            if o.get('q'):
                Q.on_claim(db, g['name'], uid, o)
            if o['kind'] == 'polish':
                db.execute("UPDATE cand SET status='claimed', taken_by=? WHERE unit=? AND kind='polish' AND status='pending'",
                           (uid, o['parent']))
            if o.get('cand'):
                db.execute("UPDATE cand SET status='claimed', taken_by=? WHERE id=?", (uid, o['cand']))
            db.execute('COMMIT')
            return dict({k: v for k, v in o.items() if k != 'q'}, id=uid, seed=seed, goal=g)
        db.execute('COMMIT')
        return None
    except BaseException:
        db.execute('ROLLBACK')
        raise


# ------------------------------------------------------------------ record

def record(db, u, res, cpu_s, wall, value, bars, retire=True, now=None, goal=None):
    """one transaction: gains vs the per-n state, the unit row, candidates, screen history, the polished family best"""
    now = time.time() if now is None else now
    name, n = u['goal']['name'], u['n']
    db.execute('BEGIN IMMEDIATE')
    try:
        row = db.execute('SELECT sbest, fbest, f0p FROM nstate WHERE goal=? AND n=?', (name, n)).fetchone()
        sbest, fbest, f0p = row if row else (math.inf, None, None)     # refresh units have no n
        gain = sgain = 0.0; bar = 0
        status = res.get('status', 'done')
        if u['kind'] == 'reference' and res.get('f0p') is not None:
            f0p = res['f0p']; fbest = min(fbest if fbest is not None else f0p, f0p)
            db.execute('UPDATE nstate SET f0p=?, fbest=? WHERE goal=? AND n=?', (f0p, fbest, name, n))
        if u['kind'] == 'polish' and res.get('fbest') is not None and fbest is not None and res['fbest'] < fbest:
            gain = fbest - res['fbest']; fbest = res['fbest']
            db.execute('UPDATE nstate SET fbest=? WHERE goal=? AND n=?', (fbest, name, n))
        jb = res.get('jbest')
        if jb is not None and jb < sbest:
            sgain = sbest - jb
            b = bars.get(n)
            bar = int(b is not None and jb < b <= sbest)
            db.execute('UPDATE nstate SET sbest=?, retired=MAX(retired, ?) WHERE goal=? AND n=?', (jb, int(bar and retire), name, n))
        reward = value.get('frontier', 1e5) * (gain + sgain) + value.get('bar', 100.0) * bar
        if goal is not None and goal.queue:                          # own best of the unit's condition, score units
            from . import queue as Q
            gain, reward = Q.on_record(db, u, res, goal, now)
        db.execute("UPDATE unit SET status=?, ended=?, cpu=?, wall=?, gain=?, store_gain=?, bar=?, reward=?, jbest=?, info=? "
                   "WHERE id=?", (status, now, cpu_s, wall, gain, sgain, bar, reward, jb, json.dumps(res.get('info', {})), u['id']))
        for kind, pid, s in res.get('cands', []):
            db.execute("INSERT INTO cand(goal, n, kind, unit, packing, s, status) VALUES (?,?,?,?,?,?,'pending')",
                       (name, n, kind, u['id'], pid, s))
        if res.get('screens'):
            db.executemany('INSERT INTO screen VALUES (?,?,?,?)', [(u['arm'], k, s, u['id']) for k, s in res['screens']])
        db.execute("UPDATE cand SET status='done' WHERE taken_by=?", (u['id'],))
        for v, c, k in res.get('gen_done', []):
            db.execute('INSERT OR IGNORE INTO gen_done VALUES (?,?,?,?,?)', (u['arm'], v, c, k, u['id']))
        db.execute('COMMIT')
    except BaseException:
        db.execute('ROLLBACK')
        raise
    return dict(gain=gain, store_gain=sgain, bar=bar, reward=reward)


# ------------------------------------------------------------------ running a unit (in the worker)

def donors_fast(goal, n, m, rng, match=True, window=None):
    """goal.donors over the cached family column (pk/families.near): same rows, ms instead of seconds"""
    from . import families as F, goal as G
    window = goal.donor_window if window is None else window
    st = store_ro()
    try:
        ex = goal.progress.get('exclude', []) if match else []
        rows = [(i, s, sq) for i, s, sq, _ in F.near(st, n, G.kgrid(n) + 0.05, ex, 400, window=window)]
    finally:
        st.close()
    if len(rows) <= m:
        return rows
    picks = [rows[0]] + rng.sample(rows[1:], m - 1)
    return sorted(picks, key=lambda r: r[1])


def _unit_dir(goal, u):
    jid = f"u{u['id']:07d}"; out = os.path.join(goal.dir, jid)   # fixed width: source LIKE '<jid>%' is exact
    os.makedirs(os.path.join(out, 'starts'), exist_ok=True)
    return jid, out


def _family_ref(goal, n):
    db = connect(); r = db.execute('SELECT COALESCE(fbest, f0, sbest) FROM nstate WHERE goal=? AND n=?', (goal.goal, n)).fetchone()
    db.close()
    return r[0] if r else None


def _ingest(goal, u, jid, out, t0, J):
    """starts + walk archive into the store (goal.ingest_job), converge this unit's packings, and its new packings at n
    -> (jbest, polish candidates, rows [(id, s, family)] of new packings, best first)"""
    from . import goal as G, families as F, converge
    from .store import Store
    json.dump(J, open(os.path.join(out, 'job.json'), 'w'))
    n = u['n']; k = G.kgrid(n)
    st = Store()
    try:
        G.ingest_job(st, goal, jid, out)
        converge.run(st, lo=n, hi=n, procs=1, ids=G._job_ids(st.db, goal, jid, t0), log=lambda *a: None)
        ids = G._job_ids(st.db, goal, jid, t0)
        rows = []
        if ids:
            q = ','.join('?' * len(ids))
            rows = st.db.execute(f"SELECT id, s, coords, family FROM packing WHERE id IN ({q}) AND n=? AND {G.FEAS} AND s<? "
                                 "AND added>=? ORDER BY s", (*ids, n, k - 1e-9, t0)).fetchall()   # new packings only (a
    finally:                                                      # walk re-imports its start: an old id, nothing to gain)
        st.close()
    ex = goal.progress.get('exclude', [])
    fam = F.fill(rows[:4 * POLISH_TOP], n) if rows else {}
    fb = _family_ref(goal, n)
    pol = [('polish', i, s) for i, s, _, _ in rows[:4 * POLISH_TOP]
           if fam[i] not in ex and (fb is None or s <= fb + u['goal']['polish_margin'])][:POLISH_TOP]
    return (rows[0][1] if rows else None), pol, [(i, s, fam.get(i)) for i, s, _, _ in rows]


def run_construct(u, goal, arms):
    """a chain arm: serial screen-then-promote construction for the allowance (no walk: walks are units of their own)"""
    import numpy as np
    sys.path.insert(0, HERE)
    import mcmin
    from layout import full_lines
    from . import goal as G, stages as S, screen as SCR
    from .packing import Packing
    arm = next(a for a in arms if a['id'] == u['arm'])
    n = arm['n']; k = G.kgrid(n); rng = random.Random(u['seed'])
    jid, out = _unit_dir(goal, u)
    json.dump(arm, open(os.path.join(out, 'arm.json'), 'w'))
    t0 = time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(time.time() - 1))
    c0 = cpu()
    gen_led = str(arm.get('chain') or '').strip().startswith('gen.')
    if u.get('start') and not gen_led:                            # a queue claim: the parent is chosen
        st_ = store_ro(); r = st_.execute('SELECT id, s, coords FROM packing WHERE id=?', (u['start'],)).fetchone(); st_.close()
        D = [(r[0], r[1], np.frombuffer(r[2], '<f8').reshape(-1, 3).copy())] if r else []
    else:
        D = [(None, None, None)] if gen_led else donors_fast(goal, n + arm['offset'], goal.donors, rng, arm.get('match_donors', True))
    if not D:
        return dict(status='empty', info=dict(donors=0))
    cfg = dict(SCR.DEFAULT, **(goal.screen if isinstance(goal.screen, dict) else {}))
    db = connect()
    hist = [s for (s,) in db.execute('SELECT score FROM screen WHERE arm=? AND k=? ORDER BY rowid DESC LIMIT ?',
                                     (arm['id'], cfg['k'], HIST))]
    db.close()
    budget = 0.95 * u['goal']['unit_cpu_s'] - 3.0                 # ~3 s: ingest + converge
    rt = []

    def cap():
        return max(30.0, goal.relax_cap * float(np.median(rt))) if len(rt) >= 3 else max(60.0, budget / 2)

    def pick_start():
        pid, s, sq = D[rng.randrange(len(D))]
        return pid, (None if gen_led else S.State(Packing(s, sq)))

    def keep(st, lg):
        if st.p.n != n or st.p.s >= k + 0.05:
            return False
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, 'e.txt'); mcmin.write_deg(p, float(st.p.s), [tuple(map(float, q)) for q in st.p.sq])
            return not full_lines(p, k)

    def on_relax(lg):
        rt.extend(r['cpu'] for r in lg if r['kind'] == 'relax' and not r['info'].get('failed'))

    log = dict(arm=arm['id'], donors=[d[0] for d in D])
    fg = gen_stage(arm); done = []
    if fg is not None:                                            # deterministic generator: each combination once
        version = fg.version()
        db = connect(); todo = _todo(db, arm, fg, version); db.close()
        spec = S.parse(arm['chain']); kept, recs = [], []
        for c, var in todo:
            if cpu() - c0 > budget:
                break
            sp = [(spec[0][0], dict(spec[0][1], **c))] + spec[1:]
            t = cpu()
            try:
                st, lg = S.run_chain(sp, dict(n=n, fq_timeout=cap()), rng.randrange(1 << 30), st=None)
            except Exception as e:
                recs.append(dict(combo=c, variant=var, failed=repr(e)[:80])); done.append((version, json.dumps(c, sort_keys=True), var))
                continue
            on_relax(lg)
            ok = not S.failed(lg) and keep(st, lg)
            recs.append(dict(combo=c, variant=var, s=float(st.p.s), ok=ok, cpu=cpu() - t,
                             build_s=lg[0]['s'], relax_cpu=sum(r['cpu'] for r in lg[1:])))
            done.append((version, json.dumps(c, sort_keys=True), var))
            if ok:
                kept.append((float(st.p.s), st.meta.get('parent'), np.asarray(st.p.sq, float)))
        log.update(built=len(recs), generator=dict(version=version, todo=len(todo), done=len(done)))
    else:
        kept, recs = SCR.construct_serial(arm['chain'], dict(n=n, cap=cap, on_relax=on_relax), budget - (cpu() - c0), cpu, rng,
                                          pick_start, keep, cfg, hist, log)
    json.dump(recs, open(os.path.join(out, 'screen.json'), 'w'))
    kept.sort(key=lambda t: t[0]); seen = set(); parents = {}
    for s, pid, sq in kept:
        if round(s, 6) in seen:
            continue
        seen.add(round(s, 6)); f = os.path.join(out, 'starts', f'p{pid}_{len(parents)}.txt'); G._write(f, s, sq)
        parents[os.path.basename(f)] = pid
    J = dict(log, built=log.get('built', sum(1 for r in recs if 'traj' in r)), kept=len(parents), parents=parents, cpu=cpu() - c0)
    jbest, pol, rows = _ingest(goal, u, jid, out, t0, J)
    ex = goal.progress.get('exclude', [])
    walk = [('walk', i, s) for i, s, f in rows if f is not None and f not in ex][:WALK_TOP]
    screens = [(cfg['k'], SCR.score_at(r['traj'], cfg['k'])) for r in recs if r.get('traj')]
    return dict(jbest=jbest, cands=pol + walk, screens=[x for x in screens if x[1] is not None], gen_done=done, new=rows,
                info=dict(built=J['built'], screened=log.get('screened'), relaxed=log.get('full_relaxed'), kept=len(parents),
                          new=len(rows), **({'generator': log['generator']} if 'generator' in log else {})))


def _walk_start(goal, u):
    """a walk arm's start: in the donor window at n, the packing with the least walk CPU so far, near the best first
    (score = (s - best) / window + log(1 + walk CPU / allowance)), skipping starts other running walks hold"""
    from . import families as F, goal as G
    n = u['n']; st = store_ro()
    try:
        rows = F.near(st, n, G.kgrid(n) + 0.05, goal.progress.get('exclude', []), 400, window=goal.donor_window)
        if not rows:
            return None
        q = ','.join('?' * len(rows))
        spent = dict(st.execute(f"SELECT packing, SUM(cpu_s) FROM expansion WHERE mechanism=? AND packing IN ({q}) GROUP BY packing",
                                (WALKED, *[r[0] for r in rows])))
    finally:
        st.close()
    db = connect()
    busy = {i for (i,) in db.execute("SELECT start FROM unit WHERE status='running' AND kind='walk' AND n=? AND id!=?",
                                     (n, u['id'])) if i}
    db.close()
    s0, A = rows[0][1], u['goal']['unit_cpu_s']
    score = lambda r: (r[1] - s0) / goal.donor_window + math.log1p(spent.get(r[0], 0.0) / A)
    free = [r for r in rows if r[0] not in busy] or rows
    return min(free, key=score)


def run_walk(u, goal, arms):
    """explore.py from one start for the allowance; the start is a construct unit's candidate or (walk arm) chosen by
    _walk_start.  Every walk leaves an expansion row on its start."""
    import numpy as np
    from . import goal as G
    from .store import Store
    jid, out = _unit_dir(goal, u)
    t0 = time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(time.time() - 1))
    c0 = cpu()
    if u.get('start'):
        st = store_ro(); r = st.execute('SELECT id, s, coords FROM packing WHERE id=?', (u['start'],)).fetchone(); st.close()
        start = (r[0], r[1], np.frombuffer(r[2], '<f8').reshape(-1, 3)) if r else None
    else:
        start = _walk_start(goal, u)
        if start is not None:
            db = connect(); db.execute('UPDATE unit SET start=? WHERE id=?', (start[0], u['id'])); db.close()
    if start is None:
        return dict(status='empty', info=dict(donors=0))
    pid, s, sq = start[:3]
    f = os.path.join(out, 'starts', f'id{pid}.txt'); G._write(f, s, sq)
    left = 0.95 * u['goal']['unit_cpu_s'] - 4.0 - (cpu() - c0)    # ~4 s: explore start-up, ingest + converge
    with open(os.path.join(out, 'walk.log'), 'w') as fh:
        # CPU budget of the explore tree (a step in flight still finishes: up to one step over, charged to the unit);
        # minutes and the timeout are hang guards only
        try:
            subprocess.run([sys.executable, 'explore.py', '--n', str(u['n']), '--starts', f, '--procs', '1',
                            '--cpu-h', f'{max(10.0, left) / 3600:.6f}', '--minutes', f'{max(1.0, 3 * left / 60):.2f}',
                            '--out', os.path.join(out, 'ex'), '--seed', str(u['seed'] % 100000), '--polish-gate', '2e-3',
                            '--weights', 'melt=2,jmelt=1,jrefill=1'],
                           cwd=HERE, stdout=fh, stderr=subprocess.STDOUT, timeout=10 * u['goal']['unit_cpu_s'] + 300)
        except subprocess.TimeoutExpired:
            fh.write('\nunits: explore killed at the hang guard\n')
    wcpu = cpu() - c0
    J = dict(arm=u['arm'], donors=[pid], parents={os.path.basename(f): pid}, cpu=wcpu)
    jbest, pol, rows = _ingest(goal, u, jid, out, t0, J)
    st = Store()
    try:
        st.record_expansion(pid, u['n'], WALKED, wcpu, run=f'runs/goals/{goal.goal}/{jid}', nprop=len(rows),
                            best_s=jbest, frontier_s=None)
    finally:
        st.close()
    return dict(jbest=jbest, cands=pol, new=rows, info=dict(start=pid, start_s=s, new=len(rows), walk_cpu=round(wcpu, 1)))


def run_refresh(u, goal, arms):
    from . import stages as S
    return dict(info=dict(stage=u['stage'], **(S.STAGES[u['stage']].refresh() or {})))


def run_reference(u, goal, arms):
    from . import goal as G
    from .store import Store
    st = Store()
    try:
        R = G.reference(st, goal, u['n'], 1)
    finally:
        st.close()
    return dict(f0p=R['f0p'], info=dict(f0=R['f0'], minima=R.get('minima'), polish=R.get('polish')))


def run_polish(u, goal, arms):
    import numpy as np
    from . import goal as G
    from .store import Store
    db = connect()
    C = db.execute("SELECT packing FROM cand WHERE taken_by=?", (u['id'],)).fetchall(); db.close()
    n, k = u['n'], G.kgrid(u['n'])
    st = Store()
    try:
        items = []
        for (pid,) in C:
            r = st.db.execute('SELECT s, coords FROM packing WHERE id=?', (pid,)).fetchone()
            if r:
                items.append((pid, r[0], np.frombuffer(r[1], '<f8').reshape(-1, 3).copy()))
        P, stats = G.polish_cached(st.db, items, 1, goal.polish_reuse)
        G._add_polished(st, f"run:runs/goals/{goal.goal}/u{u['parent']:07d}/polish", P, {})
    finally:
        st.close()
    ok = [(v[0], v[1]) for v in P.values() if v is not None and v[0] < k - 1e-9]
    want = [s for s, sq in ok if G._wanted(goal, n, s, sq, k)]
    return dict(fbest=min(want) if want else None, jbest=min([s for s, _ in ok], default=None),
                info=dict(items=len(items), **stats))


RUN = dict(construct=run_construct, walk=run_walk, reference=run_reference, polish=run_polish, refresh=run_refresh)


def worker(i, path, max_idle, job):
    """one persistent worker: claim / run / record until the ctl stop flag, or `max_idle` seconds with nothing to do"""
    db = connect(path)
    rng = random.Random(os.getpid() * 1000003 + int(time.time()))
    wid = db.execute("INSERT INTO worker(pid, host, job, started, last_seen, status) VALUES (?,?,?,?,?,'up')",
                     (os.getpid(), socket.gethostname(), job, time.time(), time.time())).lastrowid
    idle0 = None
    try:
        while True:
            if (db.execute("SELECT value FROM ctl WHERE key='stop'").fetchone() or ['0'])[0] == '1':
                break
            w = time.time()
            u = claim(db, wid, rng)
            cs = time.time() - w
            if u is None:
                idle0 = idle0 or time.time()
                if max_idle is not None and time.time() - idle0 > max_idle:
                    break
                time.sleep(10)
                db.execute('UPDATE worker SET idle_s=idle_s+?, claim_s=claim_s+?, last_seen=?, cpu=? WHERE id=?',
                           (time.time() - w - cs, cs, time.time(), cpu(), wid))
                continue
            idle0 = None
            goal, arms = load_goal(u['goal']['path'])
            c0, w0 = cpu(), time.time()
            try:
                res = RUN[u['kind']](u, goal, arms)
            except Exception:
                res = dict(status='error', info=dict(error=traceback.format_exc()[-800:]))
            c, wl = cpu() - c0, time.time() - w0
            w1 = time.time()
            record(db, u, res, c, wl, goal.value, goal.bars, goal.retire, goal=goal)
            rs = time.time() - w1
            db.execute('UPDATE unit SET claim_s=?, record_s=? WHERE id=?', (cs, rs, u['id']))
            db.execute('UPDATE worker SET units=units+1, unit_cpu=unit_cpu+?, claim_s=claim_s+?, record_s=record_s+?, '
                       'last_seen=?, cpu=? WHERE id=?', (c, cs, rs, time.time(), cpu(), wid))
    finally:
        db.execute("UPDATE worker SET status='down', last_seen=?, cpu=? WHERE id=?", (time.time(), cpu(), wid))


# ------------------------------------------------------------------ reports

def report(db, name):
    g = db.execute('SELECT path, budget_cpu_h, unit_cpu_s, t_goal, status FROM goal WHERE name=?', (name,)).fetchone()
    if not g:
        return f'no unit goal {name}'
    med = lambda xs: sorted(xs)[len(xs) // 2] if xs else float('nan')
    out = [f'# unit goal {name} ({g[4]}): budget {g[1]} CPU-h, units {g[2]:.0f} CPU-s, start {g[3]}', '',
           '| kind | units | CPU-h | median CPU-s | median wall-s | claim + record ms (median) | gain: family / store | errors | empty | lost |',
           '|---|---|---|---|---|---|---|---|---|---|']
    tot = 0.0
    for kind, in db.execute('SELECT DISTINCT kind FROM unit WHERE goal=? ORDER BY kind', (name,)).fetchall():
        R = db.execute("SELECT cpu, wall, COALESCE(claim_s, 0) + COALESCE(record_s, 0), status, gain, store_gain FROM unit "
                       "WHERE goal=? AND kind=?", (name, kind)).fetchall()
        fin = [r for r in R if r[0] is not None]
        c = sum(r[0] for r in fin) / 3600; tot += c
        out.append(f"| {kind} | {len(R)} | {c:.3f} | {med([r[0] for r in fin]):.0f} | {med([r[1] for r in fin]):.0f} | "
                   f"{1000 * med([r[2] for r in fin]):.1f} | {sum(r[4] or 0 for r in fin):.2e} / {sum(r[5] or 0 for r in fin):.2e} | "
                   f"{sum(r[3] == 'error' for r in R)} | {sum(r[3] == 'empty' for r in R)} | {sum(r[3] == 'lost' for r in R)} |")
    W = db.execute("SELECT SUM(last_seen - started), SUM(cpu), SUM(unit_cpu), SUM(idle_s), SUM(claim_s), SUM(record_s), "
                   "COUNT(*) FROM worker WHERE id IN (SELECT DISTINCT worker FROM unit WHERE goal=?)", (name,)).fetchone()
    if W[0]:
        out += ['', f'Accounting: units CPU {tot:.3f} h.  Workers that ran it (all goals): {W[6]}, wall {W[0] / 3600:.3f} '
                f'worker-h, CPU {W[1] / 3600:.3f} h of which units {W[2] / 3600:.3f} h; utilization {W[2] / W[0]:.0%} '
                f'(units CPU / worker wall), idle {W[3] / 3600:.3f} h, claim + record {100 * (W[4] + W[5]) / max(W[2], 1e-9):.2f} % '
                f'of unit CPU, worker overhead CPU (imports, claims) {(W[1] - W[2]) / 3600:.3f} h']
    out += ['', '| arm | construct / walk / polish units | CPU-h (all its units) | gain: family / store | reward / CPU-h |',
            '|---|---|---|---|---|']
    for a, kc, kw, kp, c, gn, sg, w in db.execute(
            "SELECT arm, SUM(kind='construct'), SUM(kind='walk'), SUM(kind='polish'), SUM(cpu) / 3600, SUM(gain), "
            "SUM(store_gain), SUM(reward) FROM unit WHERE goal=? AND arm IS NOT NULL AND cpu IS NOT NULL GROUP BY arm "
            "ORDER BY SUM(reward) DESC, SUM(cpu) DESC LIMIT 25", (name,)):
        out.append(f'| {a} | {kc} / {kw} / {kp} | {c:.3f} | {gn:.2e} / {sg:.2e} | {w / max(c, 1e-9):.1f} |')
    S = db.execute('SELECT n, sbest0, sbest FROM nstate WHERE goal=? ORDER BY n', (name,)).fetchall()
    moved = [(n, a, b) for n, a, b in S if b < a - 1e-12]
    out += ['', f'Targets moved: {len(moved)} of {len(S)}; score gain {sum((a - b) / math.sqrt(n) for n, a, b in moved):.3e} '
            f'(store best vs goal start, side / sqrt n)'] + [f'  {n}: {a:.9f} -> {b:.9f}' for n, a, b in moved]
    gl = load_goal(g[0])[0]
    if gl.queue:
        from . import queue as Q
        out += ['', Q.report(db, name)]
    return '\n'.join(out)


def status(db):
    out = []
    for name, st, b in db.execute('SELECT name, status, budget_cpu_h FROM goal ORDER BY added'):
        c = db.execute("SELECT COALESCE(SUM(cpu), 0) / 3600, SUM(status='running') FROM unit WHERE goal=?", (name,)).fetchone()
        out.append(f'{name:24s} {st:7s} {c[0]:.2f} / {b} CPU-h, running {c[1] or 0}')
    up = db.execute("SELECT COUNT(*), MAX(last_seen) FROM worker WHERE status='up'").fetchone()
    out.append(f'workers up: {up[0]}' + (f' (last seen {time.time() - up[1]:.0f} s ago)' if up[1] else ''))
    return '\n'.join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='pk.units', description=__doc__.split('\n')[0])
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('add'); p.add_argument('goal'); p.add_argument('--unit-cpu', type=float, default=120.0)
    p.add_argument('--prio', type=float, default=0.0); p.add_argument('--polish-margin', type=float, default=1e-3)
    p = sp.add_parser('work'); p.add_argument('--procs', type=int, default=1); p.add_argument('--max-idle', type=float)
    sp.add_parser('status')
    p = sp.add_parser('report'); p.add_argument('goal')
    p = sp.add_parser('stop'); p.add_argument('goal', nargs='?')
    sp.add_parser('resume', help='clear the worker stop flag')
    a = ap.parse_args(argv)
    db = connect()
    if a.cmd == 'add':
        g, na = add_goal(db, a.goal, a.unit_cpu, a.prio, a.polish_margin)
        print(f'{g.goal}: {na} arms over {len(g.targets)} targets, budget {g.budget_cpu_h} CPU-h, units {a.unit_cpu:.0f} CPU-s')
    elif a.cmd == 'work':
        db.execute("INSERT OR REPLACE INTO ctl VALUES ('stop', '0')")
        job = os.environ.get('PK_JOB_NAME')
        if a.procs == 1:
            worker(0, DB, a.max_idle, job)
        else:
            import multiprocessing as mp
            ps = [mp.Process(target=worker, args=(i, DB, a.max_idle, job)) for i in range(a.procs)]
            for p in ps:
                p.start()
            for p in ps:
                p.join()
    elif a.cmd == 'status':
        print(status(db))
    elif a.cmd == 'report':
        print(report(db, a.goal))
    elif a.cmd == 'stop':
        if a.goal:
            db.execute("UPDATE goal SET status='stopped' WHERE name=?", (a.goal,))
        else:
            db.execute("INSERT OR REPLACE INTO ctl VALUES ('stop', '1')")
    elif a.cmd == 'resume':
        db.execute("INSERT OR REPLACE INTO ctl VALUES ('stop', '0')")


if __name__ == '__main__':
    main()
