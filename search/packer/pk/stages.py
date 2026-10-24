"""Stages and chains (10-22): one contract for every generator, ruin, fill, constraint, relax and perturbation, so any
idea composes with any other and a test is a chain spec, not a new script.  How to add one: ../DESIGN.md "Adding an idea".

A stage is f(st: State, rng: random.Random, ctx: dict, **params) -> (State, info dict), registered with a kind and a
declared parameter space (pk.moves specs: unspecified parameters are drawn, every value used is logged):

  gen        make a start            gen.file (and plugin generators)
  ruin       open a region           (plugin stages; drop() is the helper)
  recreate   add squares to it       recreate.<placer> for every pk.fill placer (target count ctx['n'])
  constrain  freeze / release        constrain.freeze (every current square), constrain.release
  relax      T = 0 descent           relax.<schedule> for every pk.schedules entry (frozen squares passed as fq --freeze)
  perturb    one stochastic move     perturb.<move> for every pk.moves move (melt, anneal, kick, ...); refuses frozen squares

A chain is a string of stages joined by '|', e.g.
  'gen.file(path="runs/x/n{n}.txt") | constrain.freeze | recreate.bestofK(K=16) | relax.round.15c | constrain.release'
run_chain(spec, ctx, seed) -> (State, log): one log row per stage (name, kind, params, measured CPU, side and count after).
Evaluate a chain with ../trial.py (suite tail + baseline), never on its raw output (CLAUDE.md rule).

TODO (10-22, see TODO.md): loops (repeat / accept combinators for ruin-and-recreate and basin hopping); explore.py,
chain.py and hop.py still wire their stages by hand.

Plugins (10-23): stages, placers and trial suites can live outside this repo.  Every directory on PK_PLUGINS (os.pathsep
separated; default: the sibling checkout ../plugins/packer when it exists, set PK_PLUGINS= to disable) is put on
sys.path and its pk_plugin.py imported before the registries are adapted; a plugin registers with @stage, @pk.fill.placer
or pk.moves as usual.  PLUGINS lists the loaded plugin modules (trial.py reads their TRIAL_SUITES).
"""
from __future__ import annotations
import ast, math, os, random, re, resource, tempfile, time
from dataclasses import dataclass, field, replace
import numpy as np

from .packing import Packing
from . import moves as MV

STAGES: dict = {}


@dataclass
class State:
    p: Packing
    frozen: np.ndarray = None                # bool per square (None = nothing frozen)
    region: list = None                      # rectangles (x0, y0, x1, y1) the next recreate stage fills
    meta: dict = field(default_factory=dict)  # generator facts later stages may read (e.g. lattice cells)

    def __post_init__(self):
        if self.frozen is None:
            self.frozen = np.zeros(self.p.n, bool)

    def but(self, **kw):
        return replace(self, **kw)

    def freeze_spec(self):
        """fq --freeze ranges for the frozen mask ('' = none)"""
        idx = np.flatnonzero(self.frozen); out = []
        for i in idx:
            if out and out[-1][1] == i - 1: out[-1][1] = i
            else: out.append([i, i])
        return ','.join(f'{a}-{b}' if a < b else f'{a}' for a, b in out)


def stage(_name, _kind, **space):
    def deco(f):
        f.stage_name, f.kind, f.space = _name, _kind, space
        STAGES[_name] = f
        return f
    return deco


def _cpu():
    ch = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time() + ch.ru_utime + ch.ru_stime


# ================================================================== gen

@stage('gen.file', 'gen', path=None)
def gen_file(st, rng, ctx, path):
    """a packing file; '{n}' in path is replaced by ctx n (e.g. runs/x/n{n}.txt)"""
    from .packing import read
    return State(read(path.format(n=ctx['n']))), {}


# ================================================================== ruin

def drop(st, keep, rects):
    """ruin helper: keep squares where keep is True; rects join the region (meta 'cell' filtered too)"""
    m = st.meta.copy()
    if 'cell' in m:
        m['cell'] = [c for c, k in zip(m['cell'], keep) if k]
    return State(Packing(st.p.s, st.p.sq[keep]), st.frozen[keep], (st.region or []) + rects, m)


# ================================================================== recreate (every pk.fill placer)

def _recreate(pname):
    from . import fill as FL
    f = FL.PLACERS[pname]

    def run(st, rng, ctx, **P):
        m = ctx['n'] - st.p.n
        X, extra = f(st.p.s, st.p.sq, m, rng, st.region, **P)
        new = len(X) - st.p.n
        return State(Packing(st.p.s, X), np.r_[st.frozen, np.zeros(new, bool)], None, st.meta), \
            dict(added=int(new), pen=round(FL.packing_pen(st.p.s, X), 6), **{k: v for k, v in extra.items() if k != 'angles'})
    run.__doc__ = f.__doc__
    return stage(f'recreate.{pname}', 'recreate', **f.space)(run)


# ================================================================== constrain

@stage('constrain.freeze', 'constrain')
def freeze(st, rng, ctx):
    """freeze every current square (typically after ruin, before recreate: the kept lattice)"""
    return st.but(frozen=np.ones(st.p.n, bool)), {}


@stage('constrain.release', 'constrain')
def release(st, rng, ctx):
    """unfreeze every square"""
    return st.but(frozen=np.zeros(st.p.n, bool)), {}


# ================================================================== relax (every pk.schedules entry)

def _relax(name):
    from .schedules import get

    def run(st, rng, ctx):
        from .engine import fq_quench
        opt = get(name)
        if st.frozen.any():
            opt = opt.but(freeze=st.freeze_spec())
        if opt.kick or opt.round_kick or opt.soft_kick:
            opt = opt.but(seed=rng.randrange(1 << 30))
        with tempfile.TemporaryDirectory(prefix='stg_') as td:
            r = fq_quench(st.p, opt, tmp=td, mod90=True, timeout=ctx.get('fq_timeout', 1800))
            try:
                d = r.stats(whole=True); s2, sq2 = r.read()
            except Exception as e:
                return st, dict(failed=repr(e)[:120])
        ok = d.get('min_gap', -1) >= 0 and d.get('min_wall', -1) >= 0
        info = {k: d.get(k) for k in ('s', 'viol', 's_alm', 'frozen_dmax') if d.get(k) is not None}
        if d.get('round'): info['viol_hard'] = d['round'].get('viol_hard')
        if not ok:
            return st, dict(info, failed='infeasible')
        if st.frozen.any() and (info.get('frozen_dmax') or 0) > 1e-9:
            return st, dict(info, failed='frozen moved')   # fq's repair had to rescale: the free squares did not fit
        return st.but(p=Packing(s2, sq2)), info
    run.__doc__ = 'fq quench ' + ' '.join(get(name).args())
    return stage(f'relax.{name}', 'relax')(run)


# ================================================================== perturb (every pk.moves move)

def _perturb(mname):
    f = MV.MOVES[mname]

    def run(st, rng, ctx, **P):
        if st.frozen.any():
            raise ValueError(f'perturb.{mname} moves every square: release frozen squares first')
        sq, info = f(st.p, rng, **P)
        s = info.get('s_new', st.p.s)
        return st.but(p=Packing(s, sq)), {k: v for k, v in info.items() if isinstance(v, (int, float, str, bool))}
    run.__doc__ = f.__doc__
    return stage(f'perturb.{mname}', 'perturb', **f.space)(run)


def _register_adapters():
    from . import fill as FL
    from .schedules import SCHEDULES
    for p in FL.PLACERS: _recreate(p)
    for s in SCHEDULES: _relax(s)
    for m in MV.MOVES:
        if m != 'cross':                                   # needs mates chosen by the caller (battery)
            _perturb(m)



PLUGINS: list = []


def _load_plugins():
    import importlib.util, sys
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))         # search/packer
    env = os.environ.get('PK_PLUGINS')
    dirs = env.split(os.pathsep) if env is not None else [os.path.join(here, '..', '..', '..', 'plugins', 'packer')]
    for d in filter(None, dirs):
        d = os.path.realpath(d)
        if not os.path.isfile(os.path.join(d, 'pk_plugin.py')):
            if env:
                raise ImportError(f'PK_PLUGINS: no pk_plugin.py in {d}')
            continue
        if d not in sys.path:
            sys.path.append(d)                                                   # the plugin's own modules
        spec = importlib.util.spec_from_file_location(f'pk_plugin_{len(PLUGINS)}', os.path.join(d, 'pk_plugin.py'))
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        PLUGINS.append(mod)


_load_plugins()
_register_adapters()


# ================================================================== chains

_TOK = re.compile(r'^\s*([A-Za-z_][\w.]*)\s*(?:\((.*)\))?\s*$', re.S)


def parse(spec):
    """'a(k=1, x="s") | b' -> [('a', {'k': 1, 'x': 's'}), ('b', {})]; values are Python literals"""
    out = []
    for part in spec.split('|'):
        m = _TOK.match(part)
        if not m:
            raise ValueError(f'bad stage {part!r}')
        name, args = m.group(1), m.group(2)
        if name not in STAGES:
            raise KeyError(f'unknown stage {name!r} (pk.py stages lists them)')
        kw = {}
        if args and args.strip():
            call = ast.parse(f'f({args})', mode='eval').body
            kw = {k.arg: ast.literal_eval(k.value) for k in call.keywords}
        unknown = set(kw) - set(STAGES[name].space)
        if unknown:
            raise TypeError(f'{name}: unknown parameters {sorted(unknown)}; space {sorted(STAGES[name].space)}')
        out.append((name, kw))
    return out


def run_stage(name, fixed, st, rng, ctx):
    f = STAGES[name]
    P = {k: (fixed[k] if k in fixed else MV.draw(v, rng)) for k, v in f.space.items()}
    t0 = _cpu()
    st2, info = f(st, rng, ctx, **P)
    row = dict(stage=name, kind=f.kind, params={k: v for k, v in P.items() if isinstance(v, (int, float, str, bool, type(None)))},
               cpu=round(_cpu() - t0, 2), s=float(st2.p.s), count=int(st2.p.n), frozen=int(st2.frozen.sum()), info=info)
    return st2, row


def run_chain(spec, ctx, seed, st=None):
    """run every stage of spec in order from st (None: the first stage must be a gen) -> (State, log rows)"""
    rng = random.Random(seed)
    log = []
    for name, kw in (parse(spec) if isinstance(spec, str) else spec):
        st, row = run_stage(name, kw, st, rng, ctx)
        log.append(row)
    return st, log


def table():
    """(name, kind, space, first doc line) for every stage"""
    return [(k, f.kind, f.space, ((f.__doc__ or '').strip().splitlines() or [''])[0]) for k, f in sorted(STAGES.items())]
