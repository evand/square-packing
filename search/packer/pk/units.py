"""The unit loop (DESIGN.md "Next: the unit loop", step 2; 10-26): goals spend CPU in short single-threaded units,
claimed one at a time from a shared database by persistent workers.  No dispatcher, no quantum, no drain.

  python3 -m pk.units add GOAL [--unit-cpu 180] [--prio 0] [--polish-margin 1e-3]   register a goal (its YAML budget)
  python3 -m pk.units work --procs N [--max-idle S]   N workers: claim / run / record until stopped (or idle S seconds)
  python3 -m pk.units status | report GOAL | stop [GOAL]

Launch workers with jobs.py like any compute: `pk/jobs.py start units --procs N -- python3 -m pk.units work --procs N`.

Kinds (all on the goal's arms and store, as `pk.py goal run` does, minus its job-sized steps):
  arm        one arm for `unit_cpu` CPU-s: construct + walk (goal.run_job, finish=False), then ingest + converge its own
             output; its best few wanted-family packings within `polish_margin` of the family best become candidates
  reference  the polished family best at n before the goal (goal.reference; clusters off: the 3 best), built once,
             only when a candidate at n needs it
  polish     one arm unit's candidates to a polish fixpoint (goal.polish_cached); gains credited to that arm (lineage)

Claim (one BEGIN IMMEDIATE transaction): expire lost leases; pick a goal (priority, then least budget used); in it a
reference if one is needed, else the best pending polish, else an arm by Thompson sampling (goal.Posterior) on reward
per CPU-h from the unit rows.  Record (one transaction): gains against the per-n state (store best, polished family
best) read inside the transaction, so concurrent units never pay the same gain twice.  Accounting is SELECTs over the
unit and worker rows (CPU per unit = getrusage of the worker, children included; claim / record wall per unit).
Requires `clusters: false` (the cluster census needs one serial counter).
"""
from __future__ import annotations
import argparse, dataclasses, json, math, os, random, resource, socket, sqlite3, sys, time, traceback

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = os.path.join(HERE, 'runs', 'units.sqlite')
SCHEMA = """
CREATE TABLE IF NOT EXISTS goal(name TEXT PRIMARY KEY, path TEXT, budget_cpu_h REAL, unit_cpu_s REAL, prio REAL DEFAULT 0,
  polish_margin REAL, status TEXT, t_goal TEXT, added REAL);
CREATE TABLE IF NOT EXISTS nstate(goal TEXT, n INT, sbest0 REAL, sbest REAL, f0p REAL, fbest REAL, retired INT DEFAULT 0,
  PRIMARY KEY(goal, n));
CREATE TABLE IF NOT EXISTS unit(id INTEGER PRIMARY KEY, goal TEXT, kind TEXT, n INT, arm TEXT, variant TEXT, parent INT,
  worker INT, seed INT, status TEXT, lease_until REAL, started REAL, ended REAL, cpu REAL, wall REAL, claim_s REAL,
  record_s REAL, gain REAL DEFAULT 0, store_gain REAL DEFAULT 0, bar INT DEFAULT 0, reward REAL DEFAULT 0, jbest REAL,
  info TEXT);
CREATE INDEX IF NOT EXISTS unit_goal ON unit(goal, status, arm);
CREATE INDEX IF NOT EXISTS unit_kind ON unit(goal, kind, n, status);
CREATE TABLE IF NOT EXISTS cand(id INTEGER PRIMARY KEY, goal TEXT, n INT, unit INT, packing INT, s REAL, status TEXT,
  polish_unit INT);
CREATE INDEX IF NOT EXISTS cand_pending ON cand(goal, status, n);
CREATE TABLE IF NOT EXISTS worker(id INTEGER PRIMARY KEY, pid INT, host TEXT, job TEXT, started REAL, last_seen REAL,
  cpu REAL DEFAULT 0, unit_cpu REAL DEFAULT 0, idle_s REAL DEFAULT 0, claim_s REAL DEFAULT 0, record_s REAL DEFAULT 0,
  units INT DEFAULT 0, status TEXT);
CREATE TABLE IF NOT EXISTS ctl(key TEXT PRIMARY KEY, value TEXT);
"""
LEASE = dict(arm=3.0, reference=20.0, polish=20.0)       # x unit_cpu_s (min 600 s): a lease past this is a lost worker
EMPTY_RETRY = 1800                                        # an arm with no donors is retried after this many seconds
POLISH_TOP = 3                                            # candidates per arm unit (finalize's clusters-off default)


def connect(path=DB):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    db = sqlite3.connect(path, timeout=120, isolation_level=None)   # autocommit; transactions are explicit
    db.execute('PRAGMA journal_mode=WAL'); db.execute('PRAGMA synchronous=NORMAL')
    db.executescript(SCHEMA)
    return db


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
        _GOALS[key] = (g, G.plan(g))
    return _GOALS[key]


def store_ro():
    return sqlite3.connect(f'file:{os.path.join(HERE, "runs", "store.sqlite")}?mode=ro', uri=True, timeout=60)


def add_goal(db, spec, unit_cpu=180.0, prio=0.0, margin=1e-3):
    from . import goal as G
    path = G.resolve(spec); g, arms = load_goal(path)
    os.makedirs(g.dir, exist_ok=True)
    sf = os.path.join(g.dir, 'start.json')                # goal.reference reads the goal start from here
    if not os.path.exists(sf):
        json.dump(dict(t_goal=time.strftime('%Y-%m-%dT%H:%M:%S'), units=True), open(sf, 'w'))
    t_goal = json.load(open(sf))['t_goal']
    st = store_ro()
    rows = [(g.goal, n, G.best_before(st, n, t_goal)) for n in g.targets]
    st.close()
    db.execute('BEGIN IMMEDIATE')
    db.execute('INSERT OR REPLACE INTO goal VALUES (?,?,?,?,?,?,?,?,?)',
               (g.goal, path, g.budget_cpu_h, unit_cpu, prio, margin, 'active', t_goal, time.time()))
    for name, n, s in rows:
        db.execute('INSERT OR IGNORE INTO nstate(goal, n, sbest0, sbest) VALUES (?,?,?,?)', (name, n, s, s))
        if g.retire and n in g.bars and s < g.bars[n]:
            db.execute('UPDATE nstate SET retired=1 WHERE goal=? AND n=?', (name, n))
    db.execute('COMMIT')
    return g, len(arms)


# ------------------------------------------------------------------ claim

def expire(db, now):
    lost = [i for (i,) in db.execute("SELECT id FROM unit WHERE status='running' AND lease_until<?", (now,))]
    if lost:
        q = ','.join('?' * len(lost))
        db.execute(f"UPDATE unit SET status='lost' WHERE id IN ({q})", lost)
        db.execute(f"UPDATE cand SET status='pending', polish_unit=NULL WHERE status='claimed' AND polish_unit IN ({q})", lost)
    return lost


def spent_h(db, g):
    """CPU-h used (finished units) and committed (running units at their allowance)"""
    done = db.execute("SELECT COALESCE(SUM(cpu), 0) FROM unit WHERE goal=? AND status!='running'", (g['name'],)).fetchone()[0]
    run = db.execute("SELECT COUNT(*) FROM unit WHERE goal=? AND status='running'", (g['name'],)).fetchone()[0]
    return done / 3600, run * g['unit_cpu_s'] / 3600


def pick(db, g, arms, rng, now, posterior):
    """the next option in goal g -> dict(kind, n, arm, parent) or None"""
    name = g['name']
    db.execute("UPDATE cand SET status='skipped' WHERE goal=? AND status='pending' AND s > "
               "(SELECT fbest FROM nstate WHERE nstate.goal=cand.goal AND nstate.n=cand.n) + ?", (name, g['polish_margin']))
    db.execute("UPDATE cand SET status='skipped' WHERE goal=? AND status='pending' AND n IN "
               "(SELECT n FROM nstate WHERE goal=? AND retired=1)", (name, name))
    r = db.execute("SELECT DISTINCT c.n FROM cand c JOIN nstate s ON s.goal=c.goal AND s.n=c.n WHERE c.goal=? AND "
                   "c.status='pending' AND s.f0p IS NULL AND NOT EXISTS (SELECT 1 FROM unit u WHERE u.goal=c.goal AND "
                   "u.n=c.n AND u.kind='reference' AND u.status='running') LIMIT 1", (name,)).fetchone()
    if r:
        return dict(kind='reference', n=r[0], arm=None, parent=None)
    r = db.execute("SELECT c.unit, c.n, MIN(c.s - s.fbest) d FROM cand c JOIN nstate s ON s.goal=c.goal AND s.n=c.n "
                   "WHERE c.goal=? AND c.status='pending' AND s.f0p IS NOT NULL GROUP BY c.unit ORDER BY d LIMIT 1",
                   (name,)).fetchone()
    if r:
        arm = db.execute('SELECT arm FROM unit WHERE id=?', (r[0],)).fetchone()[0]
        return dict(kind='polish', n=r[1], arm=arm, parent=r[0])
    stats = [dict(arm=a, reward=w, cpu_h=c / 3600) for a, w, c in db.execute(
        "SELECT arm, SUM(reward), SUM(cpu) FROM unit WHERE goal=? AND arm IS NOT NULL AND status NOT IN ('running', 'lost') "
        "GROUP BY arm", (name,))]
    pend = dict(db.execute("SELECT arm, COUNT(*) FROM unit WHERE goal=? AND status='running' AND kind='arm' GROUP BY arm", (name,)))
    empty = {a for (a,) in db.execute("SELECT DISTINCT arm FROM unit WHERE goal=? AND status='empty' AND ended>?",
                                      (name, now - EMPTY_RETRY))}
    retired = {n for (n,) in db.execute('SELECT n FROM nstate WHERE goal=? AND retired=1', (name,))}
    elig = [a for a in arms if a['id'] not in empty and a['n'] not in retired]
    if not elig:
        return None
    post = posterior(stats); uh = g['unit_cpu_s'] / 3600
    a = max(elig, key=lambda a: post.sample(a['id'], rng, pend.get(a['id'], 0) * uh))
    return dict(kind='arm', n=a['n'], arm=a['id'], parent=None)


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
            o = pick(db, g, arms_of(g), rng, now, posterior)
            if o is None:
                continue
            lease = now + max(600.0, LEASE[o['kind']] * g['unit_cpu_s'])
            seed = rng.randrange(1 << 30)
            cur = db.execute("INSERT INTO unit(goal, kind, n, arm, parent, worker, seed, status, lease_until, started) "
                             "VALUES (?,?,?,?,?,?,?,?,?,?)", (g['name'], o['kind'], o['n'], o['arm'], o['parent'], wid,
                                                             seed, 'running', lease, now))
            uid = cur.lastrowid
            if o['kind'] == 'polish':
                db.execute("UPDATE cand SET status='claimed', polish_unit=? WHERE unit=? AND status='pending'", (uid, o['parent']))
            db.execute('COMMIT')
            return dict(o, id=uid, seed=seed, goal=g)
        db.execute('COMMIT')
        return None
    except BaseException:
        db.execute('ROLLBACK')
        raise


# ------------------------------------------------------------------ record

def record(db, u, res, cpu_s, wall, value, bars, retire=True, now=None):
    """one transaction: gains vs the per-n state, the unit row, candidates, the polished family best"""
    now = time.time() if now is None else now
    name, n = u['goal']['name'], u['n']
    db.execute('BEGIN IMMEDIATE')
    try:
        sbest, fbest, f0p = db.execute('SELECT sbest, fbest, f0p FROM nstate WHERE goal=? AND n=?', (name, n)).fetchone()
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
        db.execute("UPDATE unit SET status=?, ended=?, cpu=?, wall=?, gain=?, store_gain=?, bar=?, reward=?, jbest=?, info=? "
                   "WHERE id=?", (status, now, cpu_s, wall, gain, sgain, bar, reward, jb, json.dumps(res.get('info', {})), u['id']))
        for pid, s in res.get('cands', []):
            db.execute("INSERT INTO cand(goal, n, unit, packing, s, status) VALUES (?,?,?,?,?,'pending')", (name, n, u['id'], pid, s))
        if u['kind'] == 'polish':
            db.execute("UPDATE cand SET status='done' WHERE polish_unit=?", (u['id'],))
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


def run_arm(u, goal, arms):
    from . import goal as G, families as F, converge
    from .store import Store
    arm = next(a for a in arms if a['id'] == u['arm'])
    jid = f"u{u['id']:07d}"; out = os.path.join(goal.dir, jid)   # fixed width: source LIKE '<jid>%' is exact
    os.makedirs(out, exist_ok=True)
    json.dump(arm, open(os.path.join(out, 'arm.json'), 'w'))
    t0 = time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(time.time() - 1))
    G.run_job(dataclasses.replace(goal, score_share=0.05), arm, out, u['goal']['unit_cpu_s'], u['seed'], finish=False)
    J = json.load(open(os.path.join(out, 'job.json')))
    if J.get('empty'):
        return dict(status='empty', info=dict(donors=0))
    n = arm['n']; k = G.kgrid(n)
    st = Store()
    try:
        G.ingest_job(st, goal, jid, out)
        ids = G._job_ids(st.db, goal, jid, t0)
        converge.run(st, lo=n, hi=n, procs=1, ids=ids, log=lambda *a: None)
        ids = G._job_ids(st.db, goal, jid, t0)
        rows = []
        if ids:
            q = ','.join('?' * len(ids))
            rows = st.db.execute(f"SELECT id, s, coords, family, added>=? FROM packing WHERE id IN ({q}) AND n=? AND {G.FEAS} "
                                 "AND s<? ORDER BY s", (t0, *ids, n, k - 1e-9)).fetchall()
    finally:
        st.close()
    jbest = rows[0][1] if rows else None
    db = connect(); fb = db.execute('SELECT fbest FROM nstate WHERE goal=? AND n=?', (goal.goal, n)).fetchone(); db.close()
    fb = fb[0] if fb and fb[0] is not None else None
    if fb is None:                                                # no reference yet: the unpolished family best
        sro = store_ro()
        try:
            fb = F.best(sro, n, k + 0.05, goal.progress.get('exclude', []), before=u['goal']['t_goal'])[0]
        finally:
            sro.close()
    cands = []
    near = [r[:4] for r in rows if r[4] and (fb is None or r[1] <= fb + u['goal']['polish_margin'])]   # new packings only
                                                                  # (a walk re-imports its donors: old ids, nothing to gain)
    if near:
        fam = F.fill(near[:4 * POLISH_TOP], n)
        cands = [(i, s) for i, s, _, _ in near[:4 * POLISH_TOP] if fam[i] not in goal.progress.get('exclude', [])][:POLISH_TOP]
    return dict(jbest=jbest, cands=cands, info=dict(built=J.get('built', 0), kept=J.get('kept', 0), donors=len(J.get('donors', [])),
                                                     packings=len(rows), job_cpu=round(J.get('cpu', 0), 1)))


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
    C = db.execute("SELECT packing FROM cand WHERE polish_unit=?", (u['id'],)).fetchall(); db.close()
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


RUN = dict(arm=run_arm, reference=run_reference, polish=run_polish)


def worker(i, path, max_idle, job):
    """one persistent worker: claim / run / record until the ctl stop flag, or `max_idle` seconds with nothing to do"""
    from . import goal as G
    G.donors = donors_fast                                        # this process only; goal.py switches after 10-26's run
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
            record(db, u, res, c, wl, goal.value, goal.bars, goal.retire)
            rs = time.time() - w1
            db.execute('UPDATE unit SET claim_s=?, record_s=? WHERE id=?', (cs, rs, u['id']))
            db.execute('UPDATE worker SET units=units+1, unit_cpu=unit_cpu+?, claim_s=claim_s+?, record_s=record_s+?, '
                       'last_seen=?, cpu=? WHERE id=?', (c, cs, rs, time.time(), cpu(), wid))
    finally:
        db.execute("UPDATE worker SET status='down', last_seen=?, cpu=? WHERE id=?", (time.time(), cpu(), wid))


# ------------------------------------------------------------------ reports

def report(db, name):
    from . import goal as G
    g = db.execute('SELECT path, budget_cpu_h, unit_cpu_s, t_goal, status FROM goal WHERE name=?', (name,)).fetchone()
    if not g:
        return f'no unit goal {name}'
    out = [f'# unit goal {name} ({g[4]}): budget {g[1]} CPU-h, units {g[2]:.0f} CPU-s, start {g[3]}', '',
           '| kind | units | CPU-h | median CPU-s | median wall-s | claim + record ms (median) | errors | empty | lost |',
           '|---|---|---|---|---|---|---|---|---|']
    tot = 0.0
    for kind, in db.execute('SELECT DISTINCT kind FROM unit WHERE goal=? ORDER BY kind', (name,)).fetchall():
        R = db.execute("SELECT cpu, wall, COALESCE(claim_s, 0) + COALESCE(record_s, 0), status FROM unit WHERE goal=? AND kind=?",
                       (name, kind)).fetchall()
        fin = [r for r in R if r[0] is not None]
        med = lambda xs: sorted(xs)[len(xs) // 2] if xs else float('nan')
        c = sum(r[0] for r in fin) / 3600; tot += c
        out.append(f"| {kind} | {len(R)} | {c:.3f} | {med([r[0] for r in fin]):.0f} | {med([r[1] for r in fin]):.0f} | "
                   f"{1000 * med([r[2] for r in fin]):.1f} | {sum(r[3] == 'error' for r in R)} | {sum(r[3] == 'empty' for r in R)} | "
                   f"{sum(r[3] == 'lost' for r in R)} |")
    W = db.execute("SELECT SUM(last_seen - started), SUM(cpu), SUM(unit_cpu), SUM(idle_s), SUM(claim_s), SUM(record_s), "
                   "COUNT(*) FROM worker WHERE id IN (SELECT DISTINCT worker FROM unit WHERE goal=?)", (name,)).fetchone()
    if W[0]:
        out += ['', f'Accounting: units CPU {tot:.3f} h.  Workers that ran it (all goals): {W[6]}, wall {W[0] / 3600:.3f} '
                f'worker-h, CPU {W[1] / 3600:.3f} h of which units {W[2] / 3600:.3f} h; utilization {W[2] / W[0]:.0%} '
                f'(units CPU / worker wall), idle {W[3] / 3600:.3f} h, claim + record {100 * (W[4] + W[5]) / max(W[2], 1e-9):.2f} % '
                f'of unit CPU, worker overhead CPU (imports, claims) {(W[1] - W[2]) / 3600:.3f} h']
    out += ['', '| arm | units | CPU-h (incl. its polish) | gain: family / store | reward / CPU-h |', '|---|---|---|---|---|']
    for a, k, c, gn, sg, w in db.execute(
            "SELECT arm, SUM(kind='arm'), SUM(cpu) / 3600, SUM(gain), SUM(store_gain), SUM(reward) FROM unit WHERE goal=? AND "
            "arm IS NOT NULL AND cpu IS NOT NULL GROUP BY arm ORDER BY SUM(reward) DESC, SUM(cpu) DESC LIMIT 25", (name,)):
        out.append(f'| {a} | {k} | {c:.3f} | {gn:.2e} / {sg:.2e} | {w / max(c, 1e-9):.1f} |')
    S = db.execute('SELECT n, sbest0, sbest, f0p, fbest, retired FROM nstate WHERE goal=? ORDER BY n', (name,)).fetchall()
    moved = [(n, a, b) for n, a, b, _, _, _ in S if b < a - 1e-12]
    out += ['', f'Targets moved: {len(moved)} of {len(S)}; score gain {sum((a - b) / math.sqrt(n) for n, a, b in moved):.3e} '
            f'(store best vs goal start, side / sqrt n)'] + [f'  {n}: {a:.9f} -> {b:.9f}' for n, a, b in moved]
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
    p = sp.add_parser('add'); p.add_argument('goal'); p.add_argument('--unit-cpu', type=float, default=180.0)
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
