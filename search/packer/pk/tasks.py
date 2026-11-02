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
A script that keeps state between tries stays an ordinary sleeping process and swaps ProcessPoolExecutor for
tasks.Executor(campaign, prio=...): every call is a unit (the function must be importable, as with spawn).
Logs: runs/tasks/<campaign>/log (every unit, appended when it ends), u<id>.log (a failed unit's own output).
"""
from __future__ import annotations
import argparse, fcntl, importlib.util, json, os, pickle, shlex, shutil, signal, sqlite3, subprocess, sys, threading, time, traceback
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
        "WHERE status='active' AND next_t<=? AND (u IS NOT NULL OR repeat IS NOT NULL) ORDER BY prio, "
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
    db.execute("UPDATE campaign SET status='cancelled' WHERE name=? AND status='active'", (name,))
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
                      "(SELECT COUNT(*) FROM unit u WHERE u.campaign=c.name AND state='pending') FROM campaign c "
                      "ORDER BY status='done', status!='active', CASE WHEN status='active' THEN prio ELSE 0 END, t0 DESC").fetchall()
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
        flag = '  LONG UNIT' if r[8] > 600 else ''
        print(f'{r[0]:24s} {names[r[1]]:6s} {r[2]:9s} {r[11]:4d} {"repeat" if r[9] else r[12]:>6} {r[4]:7d} {r[5]:5d} '
              f'{r[6] / 3600:8.2f} {r[7] / max(n, 1):7.1f} {r[8]:6.0f}{budget}{flag}{"  " + r[3] if r[3] else ""}')
    run = db.execute("SELECT COUNT(*) FROM unit WHERE state='running'").fetchone()[0]
    print(f'limit {limit(db)}; running units {run}; outside load ~{max(0, runnable() - run)}')


# ------------------------------------------------------------------ calls from a running script

class Executor:
    """submit / map like ProcessPoolExecutor, each call a unit of `campaign`.  Leaving the block drops what has not started."""
    def __init__(self, campaign, prio='normal'):
        self.campaign, self.prio, self.calls = campaign, prio, {}
        self.dir = os.path.join(DIR, campaign, f'calls{os.getpid()}'); os.makedirs(self.dir, exist_ok=True)
        self.lock, self.seq, self.stop = threading.Lock(), 0, False
        self.poller = threading.Thread(target=self._poll, daemon=True); self.poller.start()

    def submit(self, fn, *args):
        f = Future()
        with self.lock:
            self.seq += 1; p = os.path.join(self.dir, str(self.seq))
        main = getattr(sys.modules['__main__'], '__file__', None)
        with open(p + '.tmp', 'wb') as o:
            pickle.dump(os.path.abspath(main) if main else None, o); pickle.dump((fn, args), o)
        os.replace(p + '.tmp', p + '.call')
        add(self.campaign, [[sys.executable, os.path.abspath(__file__), '_call', p]], self.prio)
        with self.lock:                              # after the unit exists: a call without one reads as lost
            self.calls[p] = f
        return f

    def map(self, fn, *its):
        return (f.result() for f in [self.submit(fn, *a) for a in zip(*its)])

    def _poll(self):
        db = connect()
        while not self.stop:
            with self.lock:
                calls = list(self.calls.items())
            gone = calls and not db.execute('SELECT 1 FROM unit WHERE campaign=? LIMIT 1', (self.campaign,)).fetchone()
            for p, f in calls:
                if os.path.exists(p + '.done'):
                    ok, v = pickle.load(open(p + '.done', 'rb'))
                elif gone:                           # its unit ended without a result (killed, cancelled)
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

    def __exit__(self, *exc):
        cancel(self.campaign); self.stop = True; self.poller.join()
        shutil.rmtree(self.dir, ignore_errors=True)


def _call(p):
    """run one Executor call -> exit status"""
    with open(p + '.call', 'rb') as f:
        main = pickle.load(f)
        if main:                                     # functions of the calling script (as multiprocessing's spawn does)
            sys.path.insert(0, os.path.dirname(os.path.realpath(main)))
            spec = importlib.util.spec_from_file_location('__mp_main__', main)
            m = importlib.util.module_from_spec(spec); sys.modules['__mp_main__'] = sys.modules['__main__'] = m
            spec.loader.exec_module(m)
        try:
            fn, args = pickle.load(f)
            r = (True, fn(*args))
        except Exception as e:
            traceback.print_exc(); r = (False, e)
    with open(p + '.tmp', 'wb') as o:
        pickle.dump(r, o)
    os.replace(p + '.tmp', p + '.done')
    return 0 if r[0] else 1


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
