"""Engine: persistent fq worker + explicit quench options + the proposal pipeline (screen -> classify -> polish).

  fq = FQ()                                   # one per process; restarts itself if fq dies
  r = fq.quench(packing, QOpt(loosen=1.02, pit=8, flip_top=0))   -> QResult(packing, stats)

Knobs (all in QOpt; None = fq default):
  loosen   dilation before the ALM screen (pressure release: 1.0 none, 1.02 / 1.05 looser)
  mu0, mu_max, tol, skin, outer      ALM penalty schedule (mu0 = initial stiffness: low = soft, exploring)
  pit, flip_top, stag_tol, stag_w, rmax, r0, rmin, cc_tol   SLP polish (pit = iterations; flip_top = corner-corner
           branch search width, 0 = none; stag_tol = stop when the side stalls)
  no_alm   polish only (input already feasible);  no_polish   ALM only
  kick, seed   gaussian kick inside fq (x, y: sigma; angle: 20 sigma degrees)
"""
from __future__ import annotations
import atexit, threading, json, math, os, subprocess, time
from dataclasses import dataclass, field, asdict, replace
import numpy as np
from .packing import Packing

HERE = os.path.dirname(os.path.abspath(__file__))
FQ_BIN = os.environ.get('PK_FQ', os.path.join(HERE, '..', 'target', 'release', 'fq'))


@dataclass(frozen=True)
class QOpt:
    loosen: float | None = None
    mu0: float | None = None
    mu_max: float | None = None
    tol: float | None = None
    skin: float | None = None
    outer: int | None = None
    pit: int | None = None
    flip_top: int | None = None
    stag_tol: float | None = None
    stag_w: int | None = None
    rmax: float | None = None
    r0: float | None = None
    rmin: float | None = None
    cc_tol: float | None = None
    kick: float | None = None
    seed: int | None = None
    no_alm: bool = False
    no_polish: bool = False
    shrink: bool = False
    rect_dx: float | None = None     # container [0, s + dx] x [0, s + dy] (10-12)
    rect_dy: float | None = None
    round_r0: float | None = None    # soft-to-rigid pre-stage (fq --round-*)
    round_compress: bool = False
    round_steps: int | None = None
    soft_a0: float | None = None     # soft squares (log aspect) pre-stage (fq --soft-*; 10-12)
    soft_compress: bool = False
    soft_steps: int | None = None
    soft_kick: float | None = None
    hold_f: float | None = None

    def args(self):
        out = []
        for k, v in asdict(self).items():
            flag = '--' + k.replace('_', '-')
            if isinstance(v, bool):
                if v: out.append(flag)
            elif v is not None:
                out += [flag, repr(v) if isinstance(v, float) else str(v)]
        return out

    def but(self, **kw):
        return replace(self, **kw)


# the explorer's two stages (explore.py 10-08 recipe)
SCREEN = QOpt(pit=8, flip_top=0)            # + loosen drawn per proposal
POLISH = QOpt(no_alm=True, loosen=1.0)


@dataclass
class QResult:
    p: Packing | None
    stats: dict
    sec: float

    @property
    def ok(self):
        return self.p is not None and self.stats.get('min_gap', -1) >= 0 and self.stats.get('min_wall', -1) >= 0

    @property
    def s(self):
        return self.p.s if self.p is not None else math.inf

    @property
    def fp(self):
        return self.stats.get('fp')


class FQ:
    """Persistent `fq serve` child.  Not thread-safe; use one per process."""

    def __init__(self, binary=FQ_BIN):
        self.binary = binary
        self.proc = None
        self.calls = 0
        atexit.register(self.close)

    def _start(self):
        self.proc = subprocess.Popen([self.binary, 'serve'], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.DEVNULL, text=True, bufsize=1 << 16)

    def close(self):
        if self.proc is not None:
            try:
                self.proc.stdin.close(); self.proc.wait(timeout=5)
            except Exception:
                self.proc.kill()
            self.proc = None

    def quench(self, p: Packing, opt: QOpt = QOpt(), timeout: float | None = 3600.0) -> QResult:
        """One quench.  `timeout` (s, None = none) kills a hung fq child (it is restarted on the next call); the result is
        then an error QResult with error='timeout'.  The 1 h default is a hang guard, not a budget."""
        if self.proc is None or self.proc.poll() is not None:
            self._start()
        t0 = time.time()
        proc = self.proc
        timer = threading.Timer(timeout, proc.kill) if timeout else None
        if timer:
            timer.daemon = True; timer.start()
        req = ['quench ' + ' '.join(opt.args()), f'{p.n} {float(p.s)!r}']
        req += [f'{float(x)!r} {float(y)!r} {float(a)!r}' for x, y, a in p.sq]
        try:
            self.proc.stdin.write('\n'.join(req) + '\n'); self.proc.stdin.flush()
            head = self.proc.stdout.readline()
            if not head:
                raise EOFError('fq died')
            stats = json.loads(head)
            if 'error' in stats:
                self.proc.stdout.readline()
                return QResult(None, stats, time.time() - t0)
            n, s = self.proc.stdout.readline().split()
            rows = [self.proc.stdout.readline().split() for _ in range(int(n))]
            end = self.proc.stdout.readline().strip()
            assert end == 'END', end
        except (EOFError, BrokenPipeError, json.JSONDecodeError, AssertionError, ValueError) as e:
            timed_out = timer is not None and not timer.is_alive() and proc.poll() is not None
            self.close()
            err = f'timeout ({timeout:g} s)' if timed_out else repr(e)[:200]
            return QResult(None, dict(error=err), time.time() - t0)
        finally:
            if timer:
                timer.cancel()
        self.calls += 1
        q = Packing(float(s), np.array(rows, float), meta=dict(parent=p.meta.get('id')))
        return QResult(q, stats, time.time() - t0)


_FQ = None


def fq() -> FQ:
    """Process-wide worker (lazily started; safe in ProcessPoolExecutor workers)."""
    global _FQ
    if _FQ is None:
        _FQ = FQ()
    return _FQ
