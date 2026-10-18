#!/usr/bin/env python3
"""Background jobs without pattern matching (10-10; after `pgrep -f` / `pkill -f` repeatedly matched the controlling
shell's own command line, `setsid` forked away from `$!`, and output vanished behind `| tail` pipes).

Every job runs in its own session (the session id is the job's pid, recorded at launch), always logs to a file, and is
stopped by signalling every process group in that session.  Nothing here ever matches command lines.

  python3 pk/jobs.py start NAME [--cwd DIR] [--procs N] [--threads T] [--strict] -- CMD ARGS...
                                                             launch (refuses if NAME is still running).  Math-library
                                                             threads (OMP / OpenBLAS / MKL / numexpr) default to 1 per
                                                             process (our jobs are process-parallel; 10-14: in-process
                                                             exactsolve x 12 workers ran ~47 threads each, load ~300).
                                                             --procs N = the job's declared processes: warns when the
                                                             running threads + N exceed the target (--strict: refuse)
  python3 pk/jobs.py free                                    suggested --procs now (target - running threads, >= 1)
  python3 pk/jobs.py status [NAME]                           alive, elapsed, processes / CPU in the group, last log line,
                                                             machine load vs the ~16-process target
  python3 pk/jobs.py stop NAME [--grace 10]                  SIGTERM the group, SIGKILL after the grace period
  python3 pk/jobs.py log NAME [-n 20]                        tail of the log
  python3 pk/jobs.py peek NAME                               run the job's --peek command (early results, e.g. a
                                                             bench `score` of its partial output)

Queue (10-17): jobs wait in runs/jobs/queue.json and a small scheduler (itself a job, `_queue`, started on demand,
exits when the queue is empty) starts them when cores are free (use smoothed over polls, reset after a launch).  Accounting: the declared --procs of every running
job (a job without one counts its measured running threads) + running threads outside our jobs; a job starts when
that + its procs <= TARGET - RESERVE (normal / low) or <= TARGET x SLACK (high: quick tests go ahead of, and on top
of, a background run).  Order: priority, then FIFO; `--backfill` lets a later job of the same priority start when an
earlier one doesn't fit yet (default: strict FIFO).  Elastic jobs: `--min-procs M` with `{procs}` in the command
starts the job with as many processes as are free (M ... --procs).
  python3 pk/jobs.py queue add NAME [--prio high|normal|low] [--procs N] [--min-procs M] [--backfill] [--peek CMD]
                               [--cwd DIR] -- CMD ARGS...    ({procs} in CMD / --peek is substituted at launch)
  python3 pk/jobs.py queue ls                                waiting jobs (+ free cores by the scheduler's accounting)
  python3 pk/jobs.py queue rm NAME

Job records: runs/jobs/<NAME>/{meta.json, log}.  meta.json: cmd, cwd, pid (= process group id), start time, git commit
(+ dirty flag), load at launch.
"""
from __future__ import annotations
import argparse, json, os, signal, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
JOBS = os.path.abspath(os.environ.get('PK_JOBS', os.path.join(HERE, '..', 'runs', 'jobs')))
TARGET = 16          # busy processes the machine is meant to carry (~/math/CLAUDE.md)
SLACK = 1.25         # warn above TARGET x SLACK ("somewhat above 16 is fine")
THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS')


def running_threads(pgid=None, samples=3, dt=0.3):
    """threads in state R (on a CPU or runnable), machine-wide or in one process group; averaged over a few samples
    (instantaneous, unlike the load average, which lags by minutes)"""
    tot = 0
    for i in range(samples):
        out = subprocess.run(['ps', '-L', '-eo', 'sid=,stat='], capture_output=True, text=True).stdout
        for line in out.splitlines():
            f = line.split()
            if len(f) == 2 and f[1].startswith('R') and (pgid is None or int(f[0]) == pgid):
                tot += 1
        if i + 1 < samples:
            time.sleep(dt)
    return tot / samples - (1 if pgid is None else 0)        # minus the ps itself


def free():
    return max(1, int(TARGET - running_threads()))


def _dir(name):
    return os.path.join(JOBS, name)


def _meta(name):
    p = os.path.join(_dir(name), 'meta.json')
    return json.load(open(p)) if os.path.exists(p) else None


def _group(sid):
    """Processes of a job: [(pid, pcpu, etimes, pgid)], by numeric SESSION id (= the job's pid, start_new_session).
    Session, not process group: `timeout`, `xargs`-spawned shells etc. may move children into their own process group,
    but they stay in the session (10-14: a `timeout` child was invisible to status and stop)."""
    out = subprocess.run(['ps', '-eo', 'pid=,sid=,pcpu=,etimes=,pgid=,stat='], capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        f = line.split()
        if len(f) == 6 and int(f[1]) == sid and not f[5].startswith('Z'):     # zombies: exited, not yet reaped
            rows.append((int(f[0]), float(f[2]), int(f[3]), int(f[4])))
    return rows


def _git(cwd):
    try:
        sha = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=cwd, capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=cwd,
                                    capture_output=True, text=True).stdout.strip())
        return sha + ('+dirty' if dirty else '')
    except OSError:
        return None


def alive(name):
    m = _meta(name)
    return bool(m) and bool(_group(m['pid']))


def start(name, cmd, cwd=None, procs=None, threads=1, strict=False, peek=None, prio=None):
    if alive(name):
        sys.exit(f'job {name} is still running (pid {_meta(name)["pid"]}); stop it or pick another name')
    run = running_threads()
    if procs and run + procs * threads > TARGET * SLACK:
        msg = (f'{name}: {run:.0f} threads already running + {procs} x {threads} declared > {TARGET} target '
               f'(x{SLACK}); suggested --procs {max(1, int(TARGET - run))}')
        if strict:
            sys.exit('refused: ' + msg)
        print('warning: ' + msg)
    cwd = os.path.abspath(cwd or os.getcwd())
    d = _dir(name)
    os.makedirs(d, exist_ok=True)
    log = open(os.path.join(d, 'log'), 'w')
    env = dict(os.environ, PYTHONUNBUFFERED='1')
    for v in THREAD_VARS:                            # the caller's explicit setting wins
        env.setdefault(v, str(threads))
    # start_new_session: setsid() in the child before exec, no extra fork, so p.pid is the job and its group id
    p = subprocess.Popen(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                         start_new_session=True, env=env)
    meta = dict(name=name, pid=p.pid, cmd=cmd, cwd=cwd, start=time.strftime('%Y-%m-%dT%H:%M:%S'), t0=time.time(),
                git=_git(cwd), load=os.getloadavg()[0], running_at_start=run, procs=procs,
                threads={v: env[v] for v in THREAD_VARS}, peek=peek, prio=prio)
    json.dump(meta, open(os.path.join(d, 'meta.json'), 'w'), indent=1)
    print(f'started {name}: pid/pgid {p.pid}, log {os.path.join(d, "log")}')
    return p.pid


def stop(name, grace=10.0):
    m = _meta(name)
    if not m:
        sys.exit(f'no job {name}')
    pg = m['pid']
    if not _group(pg):
        print(f'{name}: not running'); return
    def sig(sg):                                     # every process group in the session, then stragglers by pid
        for g in {r[3] for r in _group(pg)}:
            try:
                os.killpg(g, sg)
            except ProcessLookupError:
                pass
        for r in _group(pg):
            try:
                os.kill(r[0], sg)
            except ProcessLookupError:
                pass
    sig(signal.SIGTERM)
    t = time.time()
    while _group(pg) and time.time() - t < grace:
        time.sleep(0.5)
    if _group(pg):
        sig(signal.SIGKILL); time.sleep(0.5)
    print(f'{name}: stopped' if not _group(pg) else f'{name}: still {len(_group(pg))} processes after SIGKILL')


def _last_line(name):
    p = os.path.join(_dir(name), 'log')
    if not os.path.exists(p):
        return ''
    with open(p, 'rb') as f:
        f.seek(0, 2); f.seek(max(0, f.tell() - 4096))
        lines = [l for l in f.read().decode(errors='replace').splitlines() if l.strip()]
    return lines[-1] if lines else ''


def status(name=None):
    names = [name] if name else sorted(os.listdir(JOBS)) if os.path.isdir(JOBS) else []
    busy = 0
    for nm in names:
        m = _meta(nm)
        if not m:
            continue
        g = _group(m['pid'])
        cpu = sum(r[1] for r in g)
        el = time.time() - m.get('t0', time.time())
        rt = running_threads(m['pid'], samples=2) if g else 0
        state = f'RUNNING {len(g)} procs, {rt:.0f} running threads' if g else 'done'
        flag = ''
        if g:
            busy += rt
            cap = max(len(g), m.get('procs') or 0)
            if rt > 1.5 * cap + 1:
                flag = f'  OVERSUBSCRIBED ({rt:.0f} running threads for {cap} processes: math-library threads?)'
        if g or name:
            print(f'{nm:24s} {state:36s} {el / 60:7.1f} min  pid {m["pid"]}  git {m.get("git")}{flag}')
            print(f'{"":24s} | {_last_line(nm)[:150]}')
    l1, l5, _ = os.getloadavg()
    run = running_threads()
    warn = '  OVER TARGET' if run > TARGET * SLACK else ''
    print(f'running threads now {run:.0f} vs target ~{TARGET} (these jobs {busy:.0f}); free ~{max(0, int(TARGET - run))}; '
          f'load {l1:.1f} (5 min {l5:.1f}, lags){warn}')


# ------------------------------------------------------------------ queue

QFILE = os.path.join(JOBS, 'queue.json')
RESERVE = 2              # cores normal / low jobs leave for high-priority quick tests
PRIO = {'high': 0, 'normal': 1, 'low': 2}
QNAME = '_queue'


class _qlock:
    def __enter__(self):
        import fcntl
        os.makedirs(JOBS, exist_ok=True)
        self.f = open(QFILE + '.lock', 'w'); fcntl.flock(self.f, fcntl.LOCK_EX)
        self.q = json.load(open(QFILE)) if os.path.exists(QFILE) else []
        return self

    def __exit__(self, *exc):
        if exc[0] is None:
            tmp = QFILE + '.tmp'; json.dump(self.q, open(tmp, 'w'), indent=1); os.replace(tmp, QFILE)
        self.f.close()


def committed():
    """(cores held by our running jobs by declaration, running threads outside them)"""
    held = jobs_rt = 0.0
    for nm in (os.listdir(JOBS) if os.path.isdir(JOBS) else []):
        m = _meta(nm)
        if not m or nm == QNAME or not _group(m['pid']):
            continue
        rt = running_threads(m['pid'], samples=4)
        young = time.time() - m.get('t0', 0) < 300           # a job still spinning up holds what it declared
        held += max(rt, m.get('procs') or 0) if young else rt
        jobs_rt += rt
    return held, max(0.0, running_threads(samples=4) - jobs_rt)


def q_add(name, cmd, cwd, procs, min_procs, prio, backfill, peek):
    if prio not in PRIO:
        sys.exit(f'--prio: one of {list(PRIO)}')
    if alive(name):
        sys.exit(f'job {name} is running')
    with _qlock() as L:
        if any(e['name'] == name for e in L.q):
            sys.exit(f'{name} is already queued')
        L.q.append(dict(name=name, cmd=cmd, cwd=os.path.abspath(cwd or os.getcwd()), procs=procs or 1,
                        min_procs=min_procs, prio=prio, backfill=backfill, peek=peek, added=time.time()))
        print(f'queued {name} (prio {prio}, procs {procs or 1}{f", min {min_procs}" if min_procs else ""}; '
              f'{len(L.q)} waiting)')
    if not alive(QNAME):
        start(QNAME, [sys.executable, os.path.abspath(__file__), 'queue', 'run'], cwd=HERE)


def q_ls():
    q = json.load(open(QFILE)) if os.path.exists(QFILE) else []
    held, ext = committed()
    for e in sorted(q, key=lambda e: (PRIO[e['prio']], e['added'])):
        print(f"{e['name']:24s} {e['prio']:6s} procs {e['procs']}{'/min ' + str(e['min_procs']) if e['min_procs'] else '':8s} "
              f"waiting {(time.time() - e['added']) / 60:6.1f} min  {' '.join(e['cmd'])[:90]}")
    print(f'{len(q)} waiting; held by jobs {held:.0f} + outside {ext:.0f} of target {TARGET} '
          f'(normal / low start up to {TARGET - RESERVE}, high up to {TARGET * SLACK:.0f}); scheduler '
          f'{"running" if alive(QNAME) else "not running"}')


def q_rm(name):
    with _qlock() as L:
        n0 = len(L.q); L.q = [e for e in L.q if e['name'] != name]
        print(f'removed {name}' if len(L.q) < n0 else f'{name} not queued')


def _sub(x, procs):
    return None if x is None else x.replace('{procs}', str(procs))


def q_run(poll=15.0, settle=45.0):
    """the scheduler loop: start what fits, in priority / FIFO order; exit when the queue is empty"""
    ema, started_last = None, False                        # smoothed use: thread counts swing +-3 between samples
    while True:
        started = False
        try:                                               # reap jobs this scheduler started (else they stay zombies)
            while os.waitpid(-1, os.WNOHANG)[0]:
                pass
        except ChildProcessError:
            pass
        with _qlock() as L:
            if not L.q:
                print(time.strftime('%H:%M:%S'), 'queue empty: exiting', flush=True)
                return
            held, ext = committed()
            ema = held + ext if ema is None or started_last else 0.5 * ema + 0.5 * (held + ext)
            used = ema
            blocked = set()                                # priorities with a strict-FIFO job waiting
            for e in sorted(L.q, key=lambda e: (PRIO[e['prio']], e['added'])):
                if any(PRIO[b] < PRIO[e['prio']] for b in blocked) or (e['prio'] in blocked and not e['backfill']):
                    continue
                cap = TARGET * SLACK if e['prio'] == 'high' else TARGET - RESERVE
                room = int(cap - used)
                n = min(e['procs'], room) if e['min_procs'] else e['procs']
                if (room < (e['min_procs'] or e['procs'])) and used > 0.5:  # an idle machine takes anything
                    blocked.add(e['prio']); continue
                n = max(1, n)
                if alive(e['name']):
                    print(f"{e['name']}: a job of that name is running; dropped from the queue", flush=True)
                else:
                    cmd = [_sub(c, n) for c in e['cmd']]
                    start(e['name'], cmd, e['cwd'], n, peek=_sub(e['peek'], n), prio=e['prio'])
                    print(time.strftime('%H:%M:%S'), f"started {e['name']} with {n} procs (used {used:.0f} before)", flush=True)
                    used += n; started = True
                L.q = [x for x in L.q if x['name'] != e['name']]
        started_last = started
        time.sleep(settle if started else poll)


def peek(name):
    m = _meta(name)
    if not m:
        sys.exit(f'no job {name}')
    if not m.get('peek'):
        sys.exit(f'{name} has no --peek command')
    subprocess.run(m['peek'], shell=True, cwd=m['cwd'])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('start'); p.add_argument('name'); p.add_argument('--cwd')
    p.add_argument('--procs', type=int, help='processes this job will run (checked against the target)')
    p.add_argument('--threads', type=int, default=1, help='math-library threads per process (default 1)')
    p.add_argument('--strict', action='store_true', help='refuse instead of warning when over the target')
    p.add_argument('--peek', help='shell command that prints early results (run by `jobs.py peek NAME`)')
    p.add_argument('argv', nargs='*')
    sub.add_parser('free')
    p = sub.add_parser('peek'); p.add_argument('name')
    p = sub.add_parser('queue'); p.add_argument('qcmd', choices=['add', 'ls', 'rm', 'run']); p.add_argument('name', nargs='?')
    p.add_argument('--prio', default='normal'); p.add_argument('--procs', type=int); p.add_argument('--min-procs', type=int)
    p.add_argument('--backfill', action='store_true'); p.add_argument('--peek'); p.add_argument('--cwd')
    p = sub.add_parser('stop'); p.add_argument('name'); p.add_argument('--grace', type=float, default=10.0)
    p = sub.add_parser('status'); p.add_argument('name', nargs='?')
    p = sub.add_parser('log'); p.add_argument('name'); p.add_argument('-n', type=int, default=20)
    raw = sys.argv[1:]
    cut = raw.index('--') if '--' in raw else len(raw)          # options before `--`, the command after it
    a = ap.parse_args(raw[:cut])
    if a.cmd == 'start':
        argv = raw[cut + 1:] or a.argv
        if not argv:
            sys.exit('usage: jobs.py start NAME -- CMD ARGS...')
        start(a.name, argv, a.cwd, a.procs, a.threads, a.strict, peek=a.peek)
    elif a.cmd == 'peek':
        peek(a.name)
    elif a.cmd == 'queue':
        if a.qcmd == 'add':
            argv = raw[cut + 1:]
            if not a.name or not argv:
                sys.exit('usage: jobs.py queue add NAME [options] -- CMD ARGS...')
            q_add(a.name, argv, a.cwd, a.procs, a.min_procs, a.prio, a.backfill, a.peek)
        elif a.qcmd == 'ls':
            q_ls()
        elif a.qcmd == 'rm':
            q_rm(a.name)
        else:
            q_run()
    elif a.cmd == 'free':
        print(free())
    elif a.cmd == 'stop':
        stop(a.name, a.grace)
    elif a.cmd == 'status':
        status(a.name)
    elif a.cmd == 'log':
        subprocess.run(['tail', '-n', str(a.n), os.path.join(_dir(a.name), 'log')])


if __name__ == '__main__':
    main()
