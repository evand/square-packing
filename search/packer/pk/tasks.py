#!/usr/bin/env python3
"""The compute pool (11-01): many short single-threaded units, highest priority first, no core counts.

A unit is one command: usually one application (pick a packing, apply moves, maybe store a child), or a short loop;
aim for 10-60 s, a few minutes at most.  A campaign is a stream of units with a priority: a list (`add`), or a command
repeated until its CPU budget is spent (`repeat`; the command picks its own work when it starts).  Nothing is preempted:
a new high-priority campaign takes each slot as a unit ends.

  python3 pk/tasks.py add CAMPAIGN [--prio high|normal|low] [--cwd DIR] -- CMD ARGS...    one unit
  python3 pk/tasks.py add CAMPAIGN [--prio ...] --file LINES                              one unit per line ('-': stdin)
  python3 pk/tasks.py repeat CAMPAIGN --budget CPU_H [--prio low] -- CMD ARGS...
  python3 pk/tasks.py repeat CAMPAIGN --every SECONDS [--prio ...] -- CMD   one at a time, no budget, keeps the pool up
  python3 pk/tasks.py ls                                   campaigns, slots, load
  python3 pk/tasks.py wait CAMPAIGN                        until it ends (exit 1 if paused by failures)
  python3 pk/tasks.py cancel CAMPAIGN [--kill]             drop what is pending; --kill also ends running units
  python3 pk/tasks.py resume CAMPAIGN                      after a failure pause
  python3 pk/tasks.py limit [N]                            slots (default 16); `limit 0` drains for timing runs

Rules.  A unit starts while running units + outside load < limit.  Outside load = runnable threads (/proc/loadavg)
minus running units and NOISE, smoothed: renders, long `jobs.py start` runs and everything else make the pool start fewer units.
Next unit: highest priority, then the campaign with the fewest running units, then fewest finished (equal priorities
share evenly, also when slots free one at a time).
A command exits 75 for "nothing to do now": its campaign backs off, and a repeat campaign ends after IDLE_END s of
that.  FAIL_STOP failures in a row pause a campaign.  Budgets are checked as units finish, so they can overrun slightly.
A script that runs calls in parallel uses tasks.Executor (no private pools, 11-02): it stays an ordinary sleeping
process and every call is a unit (the function must be importable, as with spawn; see the class).
Logs: runs/tasks/<campaign>/log (every unit, appended when it ends), u<id>.log (a failed unit's own output).
"""
from __future__ import annotations
import argparse, fcntl, importlib.util, json, os, pickle, resource, shlex, shutil, signal, sqlite3, subprocess, sys, threading, time, traceback
from concurrent.futures import Future

try:
    from .jobs import THREAD_VARS
except ImportError:
    from jobs import THREAD_VARS

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.abspath(os.environ.get('PK_TASKS', os.path.join(HERE, '..', 'runs', 'tasks')))
LIMIT = 16               # slots (~/math/CLAUDE.md: about 16 busy processes)
PRIO = {'high': 0, 'normal': 1, 'low': 2}
EX_IDLE = 75             # EX_TEMPFAIL: "nothing to do right now"
BACKOFF = 15.0           # s a campaign waits after an idle exit
IDLE_END = 600.0         # s of only idle exits that end a repeat campaign
FAIL_STOP = 5
NOISE = 2                # runnable threads above the units even with no outside work (units' helpers, desktop, monitoring)
LINGER = 60.0            # s the dispatcher stays after the last unit
SCHEMA = """
CREATE TABLE IF NOT EXISTS campaign(name TEXT PRIMARY KEY, prio INT, repeat TEXT, cwd TEXT, budget_h REAL, status TEXT,
  note TEXT, t0 REAL, next_t REAL DEFAULT 0, last_work REAL, units INT DEFAULT 0, fails INT DEFAULT 0,
  consec INT DEFAULT 0, cpu REAL DEFAULT 0, wall REAL DEFAULT 0, max_wall REAL DEFAULT 0, every REAL);
CREATE TABLE IF NOT EXISTS unit(id INTEGER PRIMARY KEY, campaign TEXT, cmd TEXT, cwd TEXT, state TEXT, pid INT, t0 REAL);
CREATE INDEX IF NOT EXISTS unit_state ON unit(state, campaign);
CREATE TABLE IF NOT EXISTS ctl(key TEXT PRIMARY KEY, value REAL);
"""


def connect():
    os.makedirs(DIR, exist_ok=True)
    db = sqlite3.connect(os.path.join(DIR, 'tasks.sqlite'), timeout=30, isolation_level=None)
    db.execute('PRAGMA journal_mode=WAL'); db.execute('PRAGMA synchronous=NORMAL')
    db.executescript(SCHEMA)
    if 'every' not in {r[1] for r in db.execute('PRAGMA table_info(campaign)')}:
        db.execute('ALTER TABLE campaign ADD COLUMN every REAL')
    return db


def limit(db):
    r = db.execute("SELECT value FROM ctl WHERE key='limit'").fetchone()
    return int(r[0]) if r else LIMIT


def runnable():
    """runnable threads on the machine, this reader excluded"""
    return int(open('/proc/loadavg').read().split()[3].split('/')[0]) - 1


# ------------------------------------------------------------------ the dispatcher

def _pick(db, now):
    """-> (campaign, unit id or None, repeat cmd, cwd) of the next unit, or None"""
    return db.execute(
        "SELECT name, (SELECT MIN(id) FROM unit WHERE campaign=name AND state='pending') u, repeat, cwd FROM campaign "
        "WHERE status='active' AND next_t<=? AND (u IS NOT NULL OR repeat IS NOT NULL) "
        "AND NOT (every IS NOT NULL AND EXISTS (SELECT 1 FROM unit WHERE campaign=name AND state='running')) "   # --every: one at a time, even when a unit hangs past the interval
        "ORDER BY prio, "
        "(SELECT COUNT(*) FROM unit WHERE campaign=name AND state='running'), units + fails, t0 LIMIT 1", (now,)).fetchone()


def _launch(db, pick):
    name, uid, rep, cwd = pick
    if uid is None:
        uid = db.execute("INSERT INTO unit(campaign, cmd, cwd, state) VALUES (?,?,?,'pending')", (name, rep, cwd)).lastrowid
    cmd, cwd = db.execute('SELECT cmd, cwd FROM unit WHERE id=?', (uid,)).fetchone()
    d = os.path.join(DIR, name); os.makedirs(d, exist_ok=True)
    env = dict(os.environ, PYTHONUNBUFFERED='1', PK_JOB_NAME=name)
    for v in THREAD_VARS:
        env.setdefault(v, '1')
    with open(os.path.join(d, f'u{uid}.log'), 'w') as out:
        try:
            p = subprocess.Popen(json.loads(cmd), cwd=cwd, stdout=out, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                 start_new_session=True, env=env)
        except OSError as e:                         # a command that cannot start is a failed unit
            out.write(f'cannot start: {e}\n'); _finish(db, uid, name, 127, 0.0, 0.0); return None
    db.execute("UPDATE unit SET state='running', pid=?, t0=? WHERE id=?", (p.pid, time.time(), uid))
    db.execute('UPDATE campaign SET next_t=?+every WHERE name=? AND every IS NOT NULL', (time.time(), name))
    return p, uid, name, time.time()


def _finish(db, uid, name, rc, cpu, wall):
    now = time.time()
    f = os.path.join(DIR, name, f'u{uid}.log')
    fail = rc not in (0, EX_IDLE)
    if not fail:                                     # kept only on failure; one log per campaign otherwise
        if os.path.getsize(f):
            with open(f, 'rb') as a, open(os.path.join(DIR, name, 'log'), 'ab') as b:
                b.write(a.read())
        os.remove(f)
    db.execute('BEGIN IMMEDIATE')
    db.execute('DELETE FROM unit WHERE id=?', (uid,))
    db.execute('UPDATE campaign SET cpu=cpu+?, wall=wall+?, max_wall=MAX(max_wall, ?), units=units+?, fails=fails+?, '
               'consec=CASE WHEN ? THEN consec+1 WHEN ? THEN consec ELSE 0 END, '
               'last_work=CASE WHEN ?=0 THEN ? ELSE last_work END, next_t=CASE WHEN ? THEN ? ELSE next_t END WHERE name=?',
               (cpu, wall * (rc != EX_IDLE), wall, rc == 0, fail, fail, rc == EX_IDLE, rc, now, rc == EX_IDLE, now + BACKOFF, name))
    db.execute('UPDATE campaign SET next_t=?+every WHERE name=? AND every IS NOT NULL', (now, name))
    st, rep, budget, used, consec, last, every = db.execute('SELECT status, repeat, budget_h, cpu, consec, last_work, every '
                                                            'FROM campaign WHERE name=?', (name,)).fetchone()
    left = db.execute("SELECT COUNT(*) FROM unit WHERE campaign=? AND state='pending'", (name,)).fetchone()[0]
    run = db.execute("SELECT COUNT(*) FROM unit WHERE campaign=? AND state='running'", (name,)).fetchone()[0]
    new, note = None, None
    if st == 'active':
        if consec >= FAIL_STOP:
            new, note = 'paused', f'{consec} failures in a row; see {f}'
        elif budget is not None and used >= budget * 3600:
            new, note = 'done', 'budget spent'
        elif rep and not every and now - last > IDLE_END:
            new, note = 'done', 'nothing to do'
        elif not rep and not left and not run:
            new, note = 'done', None
    if new:
        db.execute('UPDATE campaign SET status=?, note=? WHERE name=?', (new, note, name))
        if new == 'done':
            db.execute("DELETE FROM unit WHERE campaign=? AND state='pending'", (name,))
    db.execute('COMMIT')


def dispatch():
    lock = open(os.path.join(DIR, 'dispatcher.lock'), 'w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        return                                       # one is already running
    db = connect()
    n = db.execute("DELETE FROM unit WHERE state='running'").rowcount      # a dead dispatcher's units: they finish on
    if n:                                            # their own and count as outside load until then
        print(f'{n} units of an earlier dispatcher left to finish unrecorded', flush=True)
    live, ema, last_sample, last_busy = {}, 0.0, 0.0, time.time()
    signal.signal(signal.SIGTERM, lambda *a: sys.exit(0))    # units keep running, as after a crash
    while True:
        now = time.time()
        while live:                                  # reap
            try:
                pid, status, ru = os.wait4(-1, os.WNOHANG)
            except ChildProcessError:
                break
            if not pid:
                break
            p, uid, name, t0 = live.pop(pid)
            p.returncode = os.waitstatus_to_exitcode(status)
            _finish(db, uid, name, p.returncode, ru.ru_utime + ru.ru_stime, now - t0)
        if now - last_sample >= 1.0:
            ema, last_sample = 0.8 * ema + 0.2 * max(0, runnable() - 1 - len(live)), now     # - 1: the dispatcher itself
        cap = limit(db)
        outside = 0 if os.environ.get('PK_TASKS_NO_OUTSIDE') else max(0, round(ema - NOISE))     # (tests on a busy machine)
        while len(live) + outside < cap:
            p = _pick(db, now)
            if not p:
                break
            r = _launch(db, p)
            if r:
                live[r[0].pid] = r       # keep the Popen: a dropped one is reaped by the next Popen call, behind wait4
        if live:
            last_busy = now
        elif now - last_busy > LINGER and not _work_left(db):
            fcntl.flock(lock, fcntl.LOCK_UN)         # then look once more: an `add` that saw the lock held has
            if not _work_left(db):                   # committed by now, and one that comes later starts a new dispatcher
                return
            fcntl.flock(lock, fcntl.LOCK_EX)
        time.sleep(0.1)


def _work_left(db):
    return db.execute("SELECT 1 FROM campaign c WHERE status='active' AND (repeat IS NOT NULL OR EXISTS "
                      "(SELECT 1 FROM unit u WHERE u.campaign=c.name))").fetchone() is not None


def ensure_dispatcher():
    lock = open(os.path.join(DIR, 'dispatcher.lock'), 'w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        return
    lock.close()
    out = open(os.path.join(DIR, 'dispatcher.log'), 'a')
    subprocess.Popen([sys.executable, os.path.abspath(__file__), '_dispatch'], cwd=HERE, stdout=out, stderr=out,
                     stdin=subprocess.DEVNULL, start_new_session=True)


# ------------------------------------------------------------------ commands

def _campaign(db, name, prio, cwd, rep=None, budget=None, every=None):
    db.execute("INSERT INTO campaign(name, prio, repeat, cwd, budget_h, every, status, t0) VALUES (?,?,?,?,?,?,'active',?) "
               "ON CONFLICT(name) DO UPDATE SET prio=excluded.prio, repeat=excluded.repeat, cwd=excluded.cwd, "
               "budget_h=excluded.budget_h, every=excluded.every, next_t=0, "
               "status=CASE status WHEN 'paused' THEN 'paused' ELSE 'active' END, note=CASE status WHEN 'paused' THEN note END",
               (name, PRIO[prio], rep, cwd, budget, every, time.time()))
    db.execute("UPDATE campaign SET last_work=? WHERE name=?", (time.time(), name))


def add(name, cmds, prio='normal', cwd=None):
    """queue units (each a list of arguments) in a campaign -> how many"""
    db = connect(); cwd = os.path.abspath(cwd or os.getcwd())
    db.execute('BEGIN IMMEDIATE')
    _campaign(db, name, prio, cwd)
    db.executemany("INSERT INTO unit(campaign, cmd, cwd, state) VALUES (?,?,?,'pending')",
                   [(name, json.dumps(c), cwd) for c in cmds])
    db.execute('COMMIT')
    if db.execute("SELECT 1 FROM campaign WHERE name=? AND status='paused'", (name,)).fetchone():
        print(f'tasks: campaign {name} is paused by failures; units wait until `tasks.py resume {name}`', file=sys.stderr, flush=True)
    ensure_dispatcher()
    return len(cmds)


def repeat(name, cmd, budget=None, prio='low', cwd=None, every=None):
    db = connect()
    _campaign(db, name, prio, os.path.abspath(cwd or os.getcwd()), json.dumps(cmd), budget, every)
    ensure_dispatcher()


def wait(name, poll=1.0):
    """block until the campaign ends -> its status ('done', 'paused', 'cancelled')"""
    db = connect()
    while True:
        r = db.execute('SELECT status FROM campaign WHERE name=?', (name,)).fetchone()
        if r is None:
            sys.exit(f'no campaign {name}')
        running = db.execute("SELECT COUNT(*) FROM unit WHERE campaign=? AND state='running'", (name,)).fetchone()[0]
        if r[0] != 'active' and not running:
            return r[0]
        time.sleep(poll)


def cancel(name, kill=False):
    db = connect()
    db.execute("UPDATE campaign SET status='cancelled', consec=0 WHERE name=? AND status IN ('active', 'paused')", (name,))
    db.execute("DELETE FROM unit WHERE campaign=? AND state='pending'", (name,))
    if kill:
        for pid, in db.execute("SELECT pid FROM unit WHERE campaign=? AND state='running'", (name,)).fetchall():
            try:
                os.killpg(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass


def ls(show_done=5):
    db = connect()
    rows = db.execute("SELECT name, prio, status, note, units, fails, cpu, wall, max_wall, repeat IS NOT NULL, budget_h, "
                      "(SELECT COUNT(*) FROM unit u WHERE u.campaign=c.name AND state='running'), "
                      "(SELECT ?-MIN(t0) FROM unit u WHERE u.campaign=c.name AND state='running'), "
                      "(SELECT COUNT(*) FROM unit u WHERE u.campaign=c.name AND state='pending') FROM campaign c "
                      "ORDER BY status='done', status!='active', CASE WHEN status='active' THEN prio ELSE 0 END, t0 DESC",
                      (time.time(),)).fetchall()
    names = {v: k for k, v in PRIO.items()}
    print(f'{"campaign":24s} {"prio":6s} {"status":9s} {"run":>4s} {"wait":>6s} {"units":>7s} {"fail":>5s} {"CPU-h":>8s} '
          f'{"mean s":>7s} {"max s":>6s}')
    ended = 0
    for r in rows:
        ended += r[2] == 'done'
        if r[2] == 'done' and ended > show_done:
            continue
        n = r[4] + r[5]
        budget = f' of {r[10]:g}' if r[10] else ''
        flag = f'  LONG UNIT {r[12]:.0f} s' if (r[12] or 0) > 600 else ''        # a unit running now (max s: all time)
        print(f'{r[0]:24s} {names[r[1]]:6s} {r[2]:9s} {r[11]:4d} {"repeat" if r[9] else r[13]:>6} {r[4]:7d} {r[5]:5d} '
              f'{r[6] / 3600:8.2f} {r[7] / max(n, 1):7.1f} {r[8]:6.0f}{budget}{flag}{"  " + r[3] if r[3] else ""}')
    run = db.execute("SELECT COUNT(*) FROM unit WHERE state='running'").fetchone()[0]
    cap, outside = limit(db), max(0, runnable() - run)
    print(f'limit {cap}; running units {run}; outside load ~{outside}; idle ~{max(0, cap - run - outside)}')


# ------------------------------------------------------------------ calls from a running script

def default_campaign():
    """the running script's name (runs of one script may share a campaign: each drops only its own calls)"""
    main = getattr(sys.modules['__main__'], '__file__', None)
    return os.path.splitext(os.path.basename(main))[0] if main else 'calls'


class Executor:
    """submit / map like ProcessPoolExecutor, each call (or chunk of map calls) a unit of `campaign` (default: the
    script's name).  Every unit is a fresh process with the caller's environment: `initializer(*initargs)` runs in it
    before its calls, and `timeout` (s, per unit) kills an overdue unit even inside C code (its calls raise TimeoutError;
    not a campaign failure).  Leaving the block waits for the calls, as ProcessPoolExecutor does; leaving it by an
    exception (or shutdown(cancel_futures=True)) drops this executor's calls that have not started."""
    def __init__(self, campaign=None, prio='normal', initializer=None, initargs=(), timeout=None):
        self.campaign, self.prio, self.calls = campaign or default_campaign(), prio, {}
        self.job = dict(init=(initializer, tuple(initargs)) if initializer else None, timeout=timeout)
        self.dir = os.path.join(DIR, self.campaign, f'calls{os.getpid()}'); os.makedirs(self.dir, exist_ok=True)
        self.lock, self.seq, self.stop = threading.Lock(), 0, False
        db = connect()                               # a new run is new code: a pause left by an earlier one does not bind it
        if db.execute("UPDATE campaign SET status='active', consec=0, note=NULL WHERE name=? AND status='paused'",
                      (self.campaign,)).rowcount:
            print(f'tasks: campaign {self.campaign} was paused by an earlier run; resumed', file=sys.stderr, flush=True)
        self.poller = threading.Thread(target=self._poll, daemon=True); self.poller.start()

    def _put(self, fn, units, one):
        """queue units, each a list of (args, kwargs) -> their Futures (one: a unit's result is its single call's)"""
        main = getattr(sys.modules['__main__'], '__file__', None)
        paths = []
        for calls in units:
            with self.lock:
                self.seq += 1; p = os.path.join(self.dir, str(self.seq))
            with open(p + '.tmp', 'wb') as o:
                pickle.dump((os.path.abspath(main) if main else None, dict(os.environ)), o)
                pickle.dump(dict(self.job, fn=fn, calls=calls, one=one), o)
            os.replace(p + '.tmp', p + '.call'); paths.append(p)
        add(self.campaign, [[sys.executable, os.path.abspath(__file__), '_call', p] for p in paths], self.prio)
        fs = [Future() for _ in paths]
        with self.lock:                              # after the units exist: a call without one reads as lost
            self.calls.update(zip(paths, fs))
        return fs

    def submit(self, fn, /, *args, **kwargs):
        return self._put(fn, [[(args, kwargs)]], True)[0]

    def map(self, fn, *its, timeout=None, chunksize=1):
        """ordered results; chunksize calls per unit (batch calls much shorter than a unit's 10-60 s)"""
        calls = [(a, {}) for a in zip(*its)]
        fs = self._put(fn, [calls[i:i + chunksize] for i in range(0, len(calls), chunksize)], False)
        return (r for f in fs for r in f.result(timeout))

    def _poll(self):
        global UNIT_CPU
        db = connect()
        while not self.stop:
            with self.lock:
                calls = list(self.calls.items())
            if calls:                                # units of this executor still queued or running
                live = {json.loads(c)[-1] for c, in db.execute('SELECT cmd FROM unit WHERE campaign=? AND cmd LIKE ?',
                                                               (self.campaign, f'%{self.dir}/%'))}
                paused = db.execute("SELECT note FROM campaign WHERE name=? AND status='paused'", (self.campaign,)).fetchone()
            for p, f in calls:
                if os.path.exists(p + '.done'):
                    ok, v, cpu = pickle.load(open(p + '.done', 'rb'))
                    UNIT_CPU += cpu
                elif paused:                         # nothing more will run: fail now, so the script stops and says why
                    ok, v = False, RuntimeError(f'campaign {self.campaign} paused: {paused[0]}')
                elif p not in live:                  # its unit ended without a result (killed, cancelled)
                    ok, v = False, RuntimeError(f'unit of {p} ended without a result')
                else:
                    continue
                if not f.cancelled():
                    f.set_result(v) if ok else f.set_exception(v)
                with self.lock:
                    del self.calls[p]
            time.sleep(0.2)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, *exc):
        self.shutdown(cancel_futures=exc_type is not None)

    def shutdown(self, wait=True, cancel_futures=False):
        if self.stop:
            return
        if wait and not cancel_futures:
            while self.calls:
                time.sleep(0.2)
        db = connect()
        db.execute('BEGIN IMMEDIATE')
        db.execute("DELETE FROM unit WHERE campaign=? AND state='pending' AND cmd LIKE ?", (self.campaign, f'%{self.dir}/%'))
        db.execute("UPDATE campaign SET status='done' WHERE name=? AND status='active' AND repeat IS NULL "
                   "AND NOT EXISTS (SELECT 1 FROM unit WHERE campaign=name)", (self.campaign,))
        db.execute('COMMIT')
        self.stop = True; self.poller.join()
        shutil.rmtree(self.dir, ignore_errors=True)


def _calls(job):
    try:
        if job.get('init'):
            job['init'][0](*job['init'][1])
        v = [job['fn'](*a, **k) for a, k in job['calls']]
        return True, v[0] if job['one'] else v
    except Exception as e:
        traceback.print_exc(); return False, e


UNIT_CPU = 0.0       # CPU s of this process's finished Executor units (pk/cpubudget and pk/manifest add it)
LEGACY = False       # a call written by an Executor from before 11-02: its result is (ok, value), no CPU


def _done(p, r):
    u, c = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    r = r if LEGACY else (*r, u.ru_utime + u.ru_stime + c.ru_utime + c.ru_stime)
    try:
        b = pickle.dumps(r)
    except Exception as e:                           # an unpicklable result or exception
        b = pickle.dumps((False, RuntimeError(f'{r[1]!r} (unpicklable: {e})'), *r[2:]))
    with open(p + '.tmp', 'wb') as o:
        o.write(b)
    os.replace(p + '.tmp', p + '.done')
    return 0 if r[0] else 1


def _call(p):
    """run one Executor unit -> exit status"""
    with open(p + '.call', 'rb') as f:
        main, env = (lambda m: m if isinstance(m, tuple) else (m, None))(pickle.load(f))     # (path, env); a path before 11-02
        for k, v in (env or {}).items():             # the caller's environment, as a ProcessPoolExecutor worker has it
            if k not in THREAD_VARS and k != 'PK_JOB_NAME':
                os.environ[k] = v
        if main:                                     # functions of the calling script (as multiprocessing's spawn does)
            sys.path.insert(0, os.path.dirname(os.path.realpath(main)))
            spec = importlib.util.spec_from_file_location('__mp_main__', main)
            m = importlib.util.module_from_spec(spec); sys.modules['__mp_main__'] = sys.modules['__main__'] = m
            spec.loader.exec_module(m)
        try:
            job = pickle.load(f)
        except Exception as e:
            traceback.print_exc(); return _done(p, (False, e))
    if isinstance(job, tuple):                       # (fn, args): written by an Executor from before 11-02
        global LEGACY
        job, LEGACY = dict(fn=job[0], calls=[(job[1], {})], one=True), True
    t = job.get('timeout')
    if not t:
        return _done(p, _calls(job))
    pid = os.fork()                                  # the calls run in a child, so that a hang in C code can be killed
    if pid == 0:
        c = 1
        try:
            os.setpgid(0, 0); c = _done(p, _calls(job)); sys.stdout.flush(); sys.stderr.flush()
        finally:
            os._exit(c)
    try:
        os.setpgid(pid, pid)                         # (also here: the kill below may come before the child's own call)
    except OSError:
        pass
    end = time.time() + t
    while time.time() < end:
        q, st = os.waitpid(pid, os.WNOHANG)
        if q:
            return os.waitstatus_to_exitcode(st) if os.path.exists(p + '.done') else \
                _done(p, (False, RuntimeError(f'call died (exit {os.waitstatus_to_exitcode(st)})')))
        time.sleep(0.05)
    os.killpg(pid, signal.SIGKILL); os.waitpid(pid, 0)
    print(f'tasks: unit killed after {t:g} s', flush=True)
    _done(p, (False, TimeoutError(f'unit killed after {t:g} s')))
    return 0


def main():
    raw = sys.argv[1:]
    cut = raw.index('--') if '--' in raw else len(raw)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('add'); p.add_argument('name'); p.add_argument('--prio', default='normal', choices=PRIO)
    p.add_argument('--cwd'); p.add_argument('--file')
    p = sp.add_parser('repeat'); p.add_argument('name'); p.add_argument('--budget', type=float); p.add_argument('--every', type=float)
    p.add_argument('--prio', default='low', choices=PRIO); p.add_argument('--cwd')
    sp.add_parser('ls')
    p = sp.add_parser('wait'); p.add_argument('name')
    p = sp.add_parser('cancel'); p.add_argument('name'); p.add_argument('--kill', action='store_true')
    p = sp.add_parser('resume'); p.add_argument('name')
    p = sp.add_parser('limit'); p.add_argument('n', nargs='?', type=int)
    sp.add_parser('_dispatch'); sp.add_parser('up')
    p = sp.add_parser('_call'); p.add_argument('path')
    a = ap.parse_args(raw[:cut])
    argv = raw[cut + 1:]
    if a.cmd == 'add':
        if a.file:
            lines = (sys.stdin if a.file == '-' else open(a.file)).read().splitlines()
            cmds = [shlex.split(l) for l in lines if l.strip()]
        else:
            cmds = [argv] if argv else sys.exit('usage: tasks.py add CAMPAIGN -- CMD ARGS...')
        print(f'{a.name}: {add(a.name, cmds, a.prio, a.cwd)} units queued')
    elif a.cmd == 'repeat':
        if not argv or not (a.budget or a.every):
            sys.exit('usage: tasks.py repeat CAMPAIGN --budget CPU_H | --every SECONDS -- CMD ARGS...')
        repeat(a.name, argv, a.budget, a.prio, a.cwd, a.every)
    elif a.cmd == 'ls':
        ls()
    elif a.cmd == 'wait':
        st = wait(a.name); print(f'{a.name}: {st}'); sys.exit(st == 'paused')
    elif a.cmd == 'cancel':
        cancel(a.name, a.kill)
    elif a.cmd == 'resume':
        db = connect(); db.execute("UPDATE campaign SET status='active', consec=0, next_t=0, note=NULL WHERE name=?", (a.name,))
        ensure_dispatcher()
    elif a.cmd == '_call':
        sys.exit(_call(a.path))
    elif a.cmd == 'up':                          # after a reboot: resume what is active
        connect(); ensure_dispatcher()
    elif a.cmd == 'limit':
        db = connect()
        if a.n is not None:
            db.execute("INSERT OR REPLACE INTO ctl VALUES ('limit', ?)", (a.n,))
        print(limit(db))
    else:
        dispatch()


if __name__ == '__main__':
    main()
