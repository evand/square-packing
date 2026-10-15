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

One fq call path (10-14): every production fq quench goes through `fq_quench(inp, opt, mode=...)` with a QOpt, and
`QOpt.args()` is the one place fq flags are built (`QOpt.from_args` parses flag lists, fq's last-wins semantics).
mode='subprocess' (`fq quench --in A --out B FLAGS`; production: spawn cost ~0.1 %) or 'serve' (the persistent
worker above; the battery).  Same flags -> byte-identical packing output in both modes (checked 10-14,
tasks/packer-arch-2026-10-14/refactor-check.md).
"""
from __future__ import annotations
import atexit, threading, json, math, os, subprocess, time
from dataclasses import dataclass, field, fields, asdict, replace
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
    # the rest of fq's quench flags (10-14: QOpt covers every flag `fq quench` reads)
    round_kick: float | None = None
    round_p: float | None = None
    round_loosen: float | None = None
    soft_p: float | None = None
    d0: float | None = None          # --shrink schedule
    dmin: float | None = None
    hops: int | None = None          # --shrink hopping (10-15): kicks per failed shrink step
    hop_sigma: float | None = None
    hop_k: int | None = None
    vtol: float | None = None
    tr_grow: float | None = None
    tr_shrink: float | None = None
    sqp: bool = False
    finish: bool = False
    ident: bool = False
    v: bool = False                  # verbose (stderr)
    dump_lp: str | None = None       # subprocess mode only (fq main)

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

    @classmethod
    def from_args(cls, tokens, base=None):
        """fq flag list (e.g. ['--loosen', '1.0', '--no-alm'] or a --polish-extra string's split) -> QOpt, on top of
        `base`.  Last occurrence wins (fq's `arg` uses rposition).  Unknown flags raise."""
        kw, toks, i = {}, list(tokens), 0
        types = {f.name: str(f.type) for f in fields(cls)}
        while i < len(toks):
            t = toks[i]
            name = t[2:].replace('-', '_') if t.startswith('--') else None
            if name not in types:
                raise ValueError(f'unknown fq flag {t!r} in {toks}')
            ty = types[name]
            if ty == 'bool':
                kw[name] = True; i += 1; continue
            if i + 1 >= len(toks):
                raise ValueError(f'fq flag {t} needs a value')
            v = toks[i + 1]
            kw[name] = int(v) if ty.startswith('int') else float(v) if ty.startswith('float') else str(v)
            i += 2
        return replace(base if base is not None else cls(), **kw)

    def with_args(self, tokens):
        """self with fq flags `tokens` applied on top (they win, as when appended to fq's argv)."""
        return QOpt.from_args(tokens, base=self)


# the explorer's two stages (explore.py 10-08 recipe)
SCREEN = QOpt(pit=8, flip_top=0)            # + loosen drawn per proposal
POLISH = QOpt(no_alm=True, loosen=1.0)
# deep polish of a feasible layout (compose2 / rectlib seam polish; fresh re-polish with pit=2000)
POLISH_DEEP = POLISH.but(flip_top=64, pit=600)


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
        t0 = time.time()
        stats, head, rows, _ = self.request(opt, p.s, p.sq, timeout)
        if head is None:
            return QResult(None, stats, time.time() - t0)
        n, s = head.split()
        q = Packing(float(s), np.array([r.split() for r in rows], float), meta=dict(parent=p.meta.get('id')))
        return QResult(q, stats, time.time() - t0)

    def request(self, opt: QOpt, s, sq, timeout: float | None = 3600.0):
        """Raw request: -> (stats, 'n s' header line, n row lines, the JSON line) as fq printed them ({:.17} / {:.15},
        the same text as `fq quench --out` / its stdout); (stats with 'error', None, None, None) on failure."""
        if self.proc is None or self.proc.poll() is not None:
            self._start()
        proc = self.proc
        timer = threading.Timer(timeout, proc.kill) if timeout else None
        if timer:
            timer.daemon = True; timer.start()
        req = ['quench ' + ' '.join(opt.args()), f'{len(sq)} {float(s)!r}']
        req += [f'{float(x)!r} {float(y)!r} {float(a)!r}' for x, y, a in sq]
        try:
            self.proc.stdin.write('\n'.join(req) + '\n'); self.proc.stdin.flush()
            head = self.proc.stdout.readline()
            if not head:
                raise EOFError('fq died')
            stats = json.loads(head)
            if 'error' in stats:
                self.proc.stdout.readline()
                return stats, None, None, None
            nl = self.proc.stdout.readline()
            n, _ = nl.split()
            rows = [self.proc.stdout.readline() for _ in range(int(n))]
            end = self.proc.stdout.readline().strip()
            assert end == 'END', end
        except (EOFError, BrokenPipeError, json.JSONDecodeError, AssertionError, ValueError) as e:
            timed_out = timer is not None and not timer.is_alive() and proc.poll() is not None
            self.close()
            err = f'timeout ({timeout:g} s)' if timed_out else repr(e)[:200]
            return dict(error=err), None, None, None
        finally:
            if timer:
                timer.cancel()
        self.calls += 1
        return stats, nl, rows, head


_FQ = None


def fq() -> FQ:
    """Process-wide worker (lazily started; safe in ProcessPoolExecutor workers)."""
    global _FQ
    if _FQ is None:
        _FQ = FQ()
    return _FQ


# ---------------------------------------------------------------- the one fq call path (10-14)

def write_cfg(path, s, sq, mod90=False):
    """fq input file: `n s` then `x y deg` rows, every number repr(float) (exact).  mod90: angles written mod 90 (as
    mcmin.write_deg / hop.quench always did; fq reads either)."""
    with open(path, 'w') as f:
        f.write(f'{len(sq)} {float(s)!r}\n' + ''.join(
            f'{float(x)!r} {float(y)!r} {float(a % 90 if mod90 else a)!r}\n' for x, y, a in sq))


def read_cfg(path):
    """fq output file -> (s, [(x, y, deg)]) as Python floats (n/d fractions accepted, as mcmin.load_deg)."""
    from fractions import Fraction
    L = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    f = lambda v: float(Fraction(v)) if '/' in v else float(v)
    n, s = int(L[0][0]), f(L[0][1])
    return s, [(f(a), f(b), f(c)) for a, b, c, *_ in L[1:n + 1]]


@dataclass
class FQRun:
    """One `fq quench` (either mode).  stdout = fq's stdout (the JSON line; serve mode: the same line + newline),
    out = the output file (fq's --out text, byte-identical in both modes) or None if fq wrote none."""
    argv: list
    returncode: int
    stdout: str
    stderr: str
    out: str | None

    def stats(self, whole=False):
        """fq's JSON (last stdout line; whole=True: json.loads(stdout), as hop.quench always parsed it)."""
        return json.loads(self.stdout) if whole else json.loads(self.stdout.strip().splitlines()[-1])

    def read(self):
        return read_cfg(self.out)


def fq_cmd(inp, out, opt: QOpt, binary=None):
    """`fq quench` argv: [fq, quench, --in, INP, --out, OUT, *opt.args()] (fqreplay's recorder relies on this order)."""
    return [binary or FQ_BIN, 'quench', '--in', inp] + (['--out', out] if out else []) + opt.args()


def fq_quench(inp, opt: QOpt = QOpt(), *, out=None, tmp=None, mode='subprocess', binary=None, timeout=None,
              mod90=False) -> FQRun:
    """The one fq quench call.  inp = input file path, a Packing, or (s, sq) (written to tmp/a.txt; mod90 = angles
    mod 90 in that file).  out = output file (default tmp/b.txt).  mode 'subprocess' (production) runs
    `fq quench --in --out FLAGS` (a timeout raises subprocess.TimeoutExpired, as subprocess.run does); 'serve' sends
    the same packing + flags to this process's `fq serve` worker and writes fq's output text to `out` (on error /
    timeout: returncode 1, stderr = the error, out None).  Same flags + same input -> byte-identical output file."""
    if isinstance(inp, Packing):
        inp = (inp.s, inp.sq)
    if out is None:
        out = os.path.join(tmp, 'b.txt')
    if mode == 'subprocess':
        if not isinstance(inp, str):
            path = os.path.join(tmp, 'a.txt')
            write_cfg(path, inp[0], inp[1], mod90=mod90)
            inp = path
        argv = fq_cmd(inp, out, opt, binary)
        r = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
        return FQRun(argv, r.returncode, r.stdout, r.stderr, out if os.path.exists(out) else None)
    if mode != 'serve':
        raise ValueError(mode)
    if opt.dump_lp:
        raise ValueError('dump_lp needs mode=subprocess')
    if isinstance(inp, str):
        s, sq = read_cfg(inp)
    else:
        s, sq = inp
        if mod90:
            sq = [(x, y, a % 90) for x, y, a in sq]
    w = fq() if binary is None or binary == FQ_BIN else _worker(binary)
    stats, head, rows, js = w.request(opt, s, sq, timeout if timeout is not None else 3600.0)
    argv = ['serve:quench'] + opt.args()
    if head is None:
        return FQRun(argv, 1, '', stats.get('error', ''), None)
    with open(out, 'w') as f:
        f.write(head + ''.join(rows))
    return FQRun(argv, 0, js, '', out)


_WORKERS = {}


def _worker(binary):
    if binary not in _WORKERS:
        _WORKERS[binary] = FQ(binary)
    return _WORKERS[binary]
