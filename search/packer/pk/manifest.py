"""Run manifests (10-11 review, REVIEW-2026-10-11.md section 5 item 2).

    from pk import manifest
    m = manifest.start(out_dir, 'explore', seed=a.seed, args=vars(a))   # writes <out>/manifest.json at once
    ...
    m.finish(summary={...})                                              # optional: status 'ok' + summary

Recorded at start: name, argv, cwd, args, seed, git commit + dirty files (repo-relative, tracked only), sha256 of the Rust
binaries the packer calls (fq, anneal, packer in target/release), python version, host, load average, start time.
At the end (finish(), or atexit / SIGTERM / SIGHUP as status 'incomplete'): end time, wall seconds, load average, and
measured CPU seconds = this process + all waited-for descendants (getrusage SELF + CHILDREN; worker pools count once
they are shut down, subprocess.run children as they exit).  'cpu_s' is measured CPU, 'worker_wall_s' = wall x procs
if the caller passes procs: earlier logs call the latter "CPU-h" / "proc-h" (PACKER.md 10-11).

A resumed run with an existing manifest.json keeps it as manifest.<k>.json (one per segment).
"""
from __future__ import annotations
import atexit, hashlib, json, os, platform, resource, signal, socket, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
PACKER = os.path.dirname(HERE)
BINARIES = ('fq', 'anneal', 'packer')
_SHA = {}


def _sha256(path):
    if path not in _SHA:
        try:
            h = hashlib.sha256()
            with open(path, 'rb') as f:
                for b in iter(lambda: f.read(1 << 20), b''):
                    h.update(b)
            st = os.stat(path)
            _SHA[path] = dict(sha256=h.hexdigest(), mtime=time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(st.st_mtime)))
        except OSError:
            _SHA[path] = None
    return _SHA[path]


def git_state(cwd=PACKER):
    def g(*args):
        return subprocess.run(['git', *args], cwd=cwd, capture_output=True, text=True).stdout
    commit = g('rev-parse', 'HEAD').strip()
    dirty = [l[3:] for l in g('status', '--porcelain', '--untracked-files=no').splitlines() if l.strip()]   # 'XY path'
    return dict(commit=commit, dirty=dirty)


def cpu_seconds():
    """Measured CPU of this process + waited-for descendants (user + system)."""
    s = resource.getrusage(resource.RUSAGE_SELF); c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return s.ru_utime + s.ru_stime + c.ru_utime + c.ru_stime


class Manifest:
    def __init__(self, out, name, seed=None, args=None, procs=None, extra=None, fname='manifest.json'):
        os.makedirs(out, exist_ok=True)
        self.path = os.path.join(out, fname)
        if os.path.exists(self.path):
            stem, ext = os.path.splitext(fname)
            k = 1
            while os.path.exists(os.path.join(out, f'{stem}.{k}{ext}')):
                k += 1
            os.replace(self.path, os.path.join(out, f'{stem}.{k}{ext}'))
        self.t0, self.cpu0, self.pid, self.procs = time.time(), cpu_seconds(), os.getpid(), procs
        self.d = dict(name=name, status='running', argv=sys.argv, cwd=os.getcwd(), seed=seed,
                      args={k: v for k, v in (args or {}).items() if _jsonable(v)}, procs=procs,
                      git=git_state(), binaries={b: _sha256(os.path.join(PACKER, 'target', 'release', b)) for b in BINARIES},
                      python=platform.python_version(), host=socket.gethostname(), cpus=os.cpu_count(),
                      start=time.strftime('%Y-%m-%dT%H:%M:%S'), load_start=list(os.getloadavg()), **(extra or {}))
        self._write()
        self.done = False
        atexit.register(self._atexit)
        for sig in (signal.SIGTERM, signal.SIGHUP):
            if signal.getsignal(sig) in (signal.SIG_DFL, None):
                signal.signal(sig, _exit_on_signal)

    def _write(self):
        tmp = self.path + '.tmp'
        with open(tmp, 'w') as f:
            json.dump(self.d, f, indent=1)
        os.replace(tmp, self.path)

    def _close(self, status, summary=None):
        wall = time.time() - self.t0
        self.d.update(status=status, end=time.strftime('%Y-%m-%dT%H:%M:%S'), wall_s=round(wall, 1),
                      cpu_s=round(cpu_seconds() - self.cpu0, 1), load_end=list(os.getloadavg()))
        if self.procs:
            self.d['worker_wall_s'] = round(wall * self.procs, 1)
        if summary:
            self.d['summary'] = summary
        self._write()
        self.done = True

    def finish(self, summary=None):
        if not self.done:
            self._close('ok', summary)

    def _atexit(self):
        if not self.done and os.getpid() == self.pid:      # forked workers inherit atexit handlers
            self._close('incomplete')


def _exit_on_signal(sig, frame):
    sys.exit(128 + sig)


def _jsonable(v):
    try:
        json.dumps(v); return True
    except (TypeError, ValueError):
        return False


def file_sha256(path):
    r = _sha256(os.path.abspath(path))
    return r and r['sha256']


def start(out, name, seed=None, args=None, procs=None, extra=None, fname='manifest.json') -> Manifest:
    return Manifest(out, name, seed=seed, args=args, procs=procs, extra=extra, fname=fname)
