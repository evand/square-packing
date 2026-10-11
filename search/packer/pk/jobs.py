#!/usr/bin/env python3
"""Background jobs without pattern matching (10-10; after `pgrep -f` / `pkill -f` repeatedly matched the controlling
shell's own command line, `setsid` forked away from `$!`, and output vanished behind `| tail` pipes).

Every job runs in its own session / process group (the group id is the job's pid, recorded at launch), always logs to a
file, and is stopped by killing that group.  Nothing here ever matches command lines.

  python3 pk/jobs.py start NAME [--cwd DIR] -- CMD ARGS...   launch (refuses if NAME is still running)
  python3 pk/jobs.py status [NAME]                           alive, elapsed, processes / CPU in the group, last log line,
                                                             machine load vs the ~16-process target
  python3 pk/jobs.py stop NAME [--grace 10]                  SIGTERM the group, SIGKILL after the grace period
  python3 pk/jobs.py log NAME [-n 20]                        tail of the log

Job records: runs/jobs/<NAME>/{meta.json, log}.  meta.json: cmd, cwd, pid (= process group id), start time, git commit
(+ dirty flag), load at launch.
"""
from __future__ import annotations
import argparse, json, os, signal, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
JOBS = os.path.abspath(os.environ.get('PK_JOBS', os.path.join(HERE, '..', 'runs', 'jobs')))
TARGET = 16          # busy processes the machine is meant to carry (~/math/CLAUDE.md)


def _dir(name):
    return os.path.join(JOBS, name)


def _meta(name):
    p = os.path.join(_dir(name), 'meta.json')
    return json.load(open(p)) if os.path.exists(p) else None


def _group(pgid):
    """Processes in a process group: [(pid, pcpu, etimes)], read from ps by numeric group id only."""
    out = subprocess.run(['ps', '-eo', 'pid=,pgid=,pcpu=,etimes='], capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        f = line.split()
        if len(f) == 4 and int(f[1]) == pgid:
            rows.append((int(f[0]), float(f[2]), int(f[3])))
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


def start(name, cmd, cwd=None):
    if alive(name):
        sys.exit(f'job {name} is still running (pid {_meta(name)["pid"]}); stop it or pick another name')
    cwd = os.path.abspath(cwd or os.getcwd())
    d = _dir(name)
    os.makedirs(d, exist_ok=True)
    log = open(os.path.join(d, 'log'), 'w')
    env = dict(os.environ, PYTHONUNBUFFERED='1')
    # start_new_session: setsid() in the child before exec, no extra fork, so p.pid is the job and its group id
    p = subprocess.Popen(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                         start_new_session=True, env=env)
    meta = dict(name=name, pid=p.pid, cmd=cmd, cwd=cwd, start=time.strftime('%Y-%m-%dT%H:%M:%S'), t0=time.time(),
                git=_git(cwd), load=os.getloadavg()[0])
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
    os.killpg(pg, signal.SIGTERM)
    t = time.time()
    while _group(pg) and time.time() - t < grace:
        time.sleep(0.5)
    if _group(pg):
        os.killpg(pg, signal.SIGKILL); time.sleep(0.5)
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
        cpu = sum(c for _, c, _ in g)
        el = time.time() - m.get('t0', time.time())
        state = f'RUNNING {len(g)} procs, ~{cpu / 100:.1f} cores' if g else 'done'
        if g:
            busy += cpu / 100
        if g or name:
            print(f'{nm:24s} {state:32s} {el / 60:7.1f} min  pid {m["pid"]}  git {m.get("git")}')
            print(f'{"":24s} | {_last_line(nm)[:150]}')
    l1, l5, _ = os.getloadavg()
    print(f'load {l1:.1f} (5 min {l5:.1f}) vs target ~{TARGET}; these jobs ~{busy:.1f} cores; '
          f'headroom ~{TARGET - l1:+.1f}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('start'); p.add_argument('name'); p.add_argument('--cwd'); p.add_argument('argv', nargs=argparse.REMAINDER)
    p = sub.add_parser('stop'); p.add_argument('name'); p.add_argument('--grace', type=float, default=10.0)
    p = sub.add_parser('status'); p.add_argument('name', nargs='?')
    p = sub.add_parser('log'); p.add_argument('name'); p.add_argument('-n', type=int, default=20)
    a = ap.parse_args()
    if a.cmd == 'start':
        argv = a.argv[1:] if a.argv[:1] == ['--'] else a.argv
        if not argv:
            sys.exit('usage: jobs.py start NAME -- CMD ARGS...')
        start(a.name, argv, a.cwd)
    elif a.cmd == 'stop':
        stop(a.name, a.grace)
    elif a.cmd == 'status':
        status(a.name)
    elif a.cmd == 'log':
        subprocess.run(['tail', '-n', str(a.n), os.path.join(_dir(a.name), 'log')])


if __name__ == '__main__':
    main()
