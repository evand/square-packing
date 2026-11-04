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
  shrink     change n downward      shrink.drop (remove squares), shrink.lcut (inexact L removal); construction only
  grow       change n upward         (plugin stages; 10-24) a bigger box and/or added squares toward ctx['n'], left
                                     unrelaxed, with st.region set for a following recreate; refuses frozen squares

A chain is a string of stages joined by '|', e.g.
  'gen.file(path="runs/x/n{n}.txt") | constrain.freeze | recreate.bestofK(K=16) | relax.round.15c | constrain.release'
run_chain(spec, ctx, seed) -> (State, log): one log row per stage (name, kind, params, measured CPU, side and count after).
A chain keeps the best feasible packing it saw (10-26, Evan: "keep the best feasible packing found, regardless of where
in the chain it comes from"): after every stage, a packing with ctx['n'] squares and Packing.max_pen <= FEAS_PEN that
is smaller than the best so far becomes State.best, which travels with the State (so a chain split into a build call and
a relax call keeps an exact build too).  If the final state is worse, infeasible or failed, run_chain returns the best
instead and logs a 'keep.best' row; use failed(log), not the rows' info, to decide whether a chain produced something.
(Why: relax schedules are quenches, not descents; from a feasible start, e.g. an exact cut, round.15c.p5 returned a
larger side 16 times in 68 on 10-26, up to +0.011.)
Evaluate a chain with ../trial.py (suite tail + baseline), never on its raw output (CLAUDE.md rule).

TODO (10-22, see TODO.md): loops (repeat / accept combinators for ruin-and-recreate and basin hopping); explore.py,
chain.py and hop.py still wire their stages by hand.

Plugins (10-23): stages, placers and trial suites can live outside this repo.  Every directory on PK_PLUGINS (os.pathsep
separated; default: the sibling checkout ../plugins/packer when it exists, set PK_PLUGINS= to disable) is put on
sys.path and its pk_plugin.py imported before the registries are adapted; a plugin registers with @stage, @pk.fill.placer
or pk.moves as usual.  PLUGINS lists the loaded plugin modules (trial.py reads their TRIAL_SUITES).
"""
from __future__ import annotations
import ast, math, os, random, re, resource, subprocess, tempfile, time
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
    best: tuple = None                       # (s, Packing, stage): best feasible packing at ctx n so far (run_chain)

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


STATUS = ('lab', 'default', 'retired')


def stage(_name, _kind, status='default', evidence='', log_info=(), **space):
    """register f as stage _name of kind _kind over parameter space `space`.  status (11-02): 'lab' (being tried:
    `pk.py try`), 'default' (in use), 'retired' (kept so old chains replay); evidence: a short pointer (CLAIMS line, goal
    file, `try:NAME`).  `pk.py stages` lists both with the move log's use of every stage.  log_info (11-04): keys of the
    stage's info dict that the move log keeps with its params (small JSON values: a random choice the stage made)"""
    assert status in STATUS, (_name, status)

    def deco(f):
        f.stage_name, f.kind, f.space, f.status, f.evidence = _name, _kind, space, status, evidence
        f.log_info = tuple(log_info)
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


@stage('gen.random', 'gen', c=('uni', 1.2, 1.6))
def gen_random(st, rng, ctx, c):
    """ctx n squares at uniform random centres and angles in a box of side c * sqrt(n) (overlaps allowed: a rounded-
    square compress resolves them; random starts for the disk -> square pressure path, 10-24)"""
    n = int(ctx['n']); S = c * math.sqrt(n)
    sq = np.array([(rng.uniform(0.5, S - 0.5), rng.uniform(0.5, S - 0.5), rng.uniform(0.0, 90.0)) for _ in range(n)])
    return State(Packing(S, sq)), dict(side0=round(S, 4))


# ================================================================== n-changing moves (10-24, explore-loop plan)
# The elementary closures as moves: drop squares (n + d -> n), cut an inexact L (row k + 1 -> row k), add squares
# (n - d -> n).  Construction only: a relax stage follows in the chain.  All refuse frozen squares.

def _free(st):
    if st.frozen.any():
        raise ValueError('n-changing stages move or remove squares: release frozen squares first')


@stage('shrink.drop', 'shrink', pick=('choice', ['extreme', 'random']))
def shrink_drop(st, rng, ctx, pick):
    """remove st.n - ctx n squares: 'extreme' = among those touching the bounding box (removal can shrink the side),
    'random' = any; the box is kept (the following relax compacts)"""
    _free(st)
    d = st.p.n - int(ctx['n'])
    if d < 1:
        raise ValueError(f"shrink.drop: ctx n {ctx['n']} must be below the packing's {st.p.n}")
    sq = np.asarray(st.p.sq, float).copy()
    a = np.radians(sq[:, 2]); h = 0.5 * (np.abs(np.cos(a)) + np.abs(np.sin(a)))
    lo, hi = sq[:, :2] - h[:, None], sq[:, :2] + h[:, None]
    keep = np.ones(len(sq), bool)
    for _ in range(d):
        idx = np.flatnonzero(keep)
        if pick == 'extreme':
            L, H = lo[idx].min(0), hi[idx].max(0)
            ext = idx[(np.abs(lo[idx] - L) < 1e-6).any(1) | (np.abs(hi[idx] - H) < 1e-6).any(1)]
            idx = ext if len(ext) else idx
        keep[int(idx[rng.randrange(len(idx))])] = False
    return State(Packing(st.p.s, sq[keep])), dict(dropped=d)


@stage('shrink.lcut', 'shrink', cut=('choice', ['tilt', 'random']))
def shrink_lcut(st, rng, ctx, cut):
    """inexact L removal (lcut.py): a horizontal and a vertical unit strip at (cx, cy); the st.n - ctx n squares most
    inside them are removed, squares beyond a strip shift by -1, box side - 1.  cut 'tilt' centres the cross near a
    tilted square (the irregular region), 'random' anywhere"""
    _free(st)
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import lcut as LC
    d = st.p.n - int(ctx['n'])
    if d < 1:
        raise ValueError(f"shrink.lcut: ctx n {ctx['n']} must be below the packing's {st.p.n}")
    s, sq = st.p.s, np.asarray(st.p.sq, float)
    tl = np.abs(((sq[:, 2] + 45) % 90) - 45) > 1.0
    if cut == 'tilt' and tl.any():
        T = sq[tl]; cx, cy = T[rng.randrange(len(T)), :2] + np.array([rng.gauss(0, 1), rng.gauss(0, 1)])
    else:
        cx, cy = rng.uniform(0.5, s - 1.5), rng.uniform(0.5, s - 1.5)
    cx, cy = float(np.clip(cx - 0.5, 0, s - 1)), float(np.clip(cy - 0.5, 0, s - 1))
    C = LC.apply_cut([tuple(q) for q in sq], cx, cy, d)
    s1 = s - 1
    X = np.array([(min(max(x, 0.5), s1 - 0.5), min(max(y, 0.5), s1 - 0.5), a) for x, y, a in C])
    return State(Packing(s1, X)), dict(cx=round(cx, 3), cy=round(cy, 3), removed=d)


@stage('grow.add', 'grow', f=('uni', 1.0, 1.04), K=('choice', [4, 16, 64]))
def grow_add(st, rng, ctx, f, K):
    """add ctx n - st.n squares: centres dilated by f, then pk.fill bestofK (K random poses per square) anywhere"""
    _free(st)
    from . import fill as FL
    m = int(ctx['n']) - st.p.n
    if m < 1:
        raise ValueError(f"grow.add: ctx n {ctx['n']} must exceed the packing's {st.p.n}")
    sq = np.asarray(st.p.sq, float).copy(); S = st.p.s * f
    sq[:, :2] *= f
    X, info = FL.fill(S, sq, rng, placer='bestofK', m=m, K=K)
    return State(Packing(S, np.asarray(X, float))), dict(added=int(info.get('added', 0)), pen=info.get('pen'))


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


@stage('constrain.best', 'constrain', status='lab', evidence='keep-if-better cycles (11-02 move registry)')
def keep_best(st, rng, ctx):
    """continue from the best feasible packing at ctx n seen so far in this chain (State.best; nothing frozen, no
    region), so 'ruin | recreate | relax | constrain.best | ruin | ...' repeats a move on the best so far"""
    if st.best is None or (st.p is st.best[1]):
        return st, dict(reverted=False)
    return st.but(p=st.best[1], frozen=np.zeros(st.best[1].n, bool), region=None), dict(reverted=True, src=st.best[2])


for _f in (freeze, release, keep_best):
    _f.deterministic = True


# ================================================================== relax (every pk.schedules entry)

def _quench(st, rng, ctx, opt):
    """one fq quench of st under QOpt opt (frozen squares, seeds for kicks, the ctx deadline) -> (State, info)"""
    from .engine import fq_quench
    if st.frozen.any():
        opt = opt.but(freeze=st.freeze_spec())
    if opt.kick or opt.round_kick or opt.soft_kick:
        opt = opt.but(seed=rng.randrange(1 << 30))
    cap = ctx.get('fq_timeout')
    if cap:                                                    # 10-25: fq stops itself at the cap and returns
        opt = opt.but(deadline=float(cap))                     # the repaired state (best so far)
    with tempfile.TemporaryDirectory(prefix='stg_') as td:
        try:
            r = fq_quench(st.p, opt, tmp=td, mod90=True, timeout=(cap + 60) if cap else 1800)
        except subprocess.TimeoutExpired:                       # hang guard only (the deadline came first)
            return st, dict(failed='timeout', timeout=cap)
        try:
            d = r.stats(whole=True); s2, sq2 = r.read()
        except Exception as e:
            return st, dict(failed=repr(e)[:120])
    if d.get('diverged'):                                      # the compress collapsed (pair list exploded)
        return st, dict(failed='diverged', sec=d.get('sec'))
    ok = (d.get('min_gap') is None or d['min_gap'] >= 0) and (d.get('min_wall') is None or d['min_wall'] >= 0)
    info = {k: d.get(k) for k in ('s', 'viol', 's_alm', 'frozen_dmax', 'deadline') if d.get(k) is not None}
    if d.get('round'): info['viol_hard'] = d['round'].get('viol_hard')
    if not ok:
        return st, dict(info, failed='infeasible')
    if st.frozen.any() and (info.get('frozen_dmax') or 0) > 1e-9:
        return st, dict(info, failed='frozen moved')   # fq's repair had to rescale: the free squares did not fit
    return st.but(p=Packing(s2, sq2)), info


def _relax(name):
    from .schedules import get

    def run(st, rng, ctx):
        return _quench(st, rng, ctx, get(name))
    run.__doc__ = 'fq quench ' + ' '.join(get(name).args())
    o = get(name)
    run.deterministic = not (o.kick or o.round_kick or o.soft_kick)     # fq draws a seed only for kicks
    return stage(f'relax.{name}', 'relax')(run)


@stage('relax.round', 'relax', r0=0.15, steps=5, compress=True, press=True, kick=0.0, pit=8, loosen=1.0, press_outer=30)
def relax_round(st, rng, ctx, r0, steps, compress, press, kick, pit, loosen, press_outer):
    """soft-to-rigid relax with its knobs as parameters (10-27; defaults = relax.round.15c.p5): rounding radius r0,
    hardening steps, compress at r0, side free while hardening (press), kick per step, SLP iterations, start loosen,
    outer ALM iterations per press step"""
    from .engine import SCREEN
    opt = SCREEN.but(loosen=float(loosen), round_r0=float(r0), round_compress=bool(compress), round_press=bool(press),
                     round_steps=int(steps), pit=int(pit), **({'round_kick': float(kick)} if kick else {}),
                     **({'round_press_outer': int(press_outer)} if press else {}))
    return _quench(st, rng, ctx, opt)


relax_round.deterministic = lambda P: not P.get('kick')


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
    """'a(k=1, x="s") | b' -> [('a', {'k': 1, 'x': 's'}), ('b', {})]; values are Python literals.  A value that is a
    parameter space, e.g. r=('uni', 1.0, 6.0), is drawn per application (run_stage), like an unspecified parameter
    but over the chain's range"""
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


SPACE_KINDS = ('log', 'uni', 'int', 'choice')


def is_space(v):
    """a chain value that is a parameter space, e.g. ('log', 0.02, 0.3): drawn per application (10-27)"""
    return isinstance(v, tuple) and len(v) >= 2 and v[0] in SPACE_KINDS


def run_stage(name, fixed, st, rng, ctx):
    f = STAGES[name]
    P = {k: (MV.draw(fixed[k], rng) if is_space(fixed.get(k)) else fixed[k] if k in fixed else MV.draw(v, rng))
         for k, v in f.space.items()}
    t0 = _cpu()
    st2, info = f(st, rng, ctx, **P)
    row = dict(stage=name, kind=f.kind, params={k: v for k, v in P.items() if isinstance(v, (int, float, str, bool, type(None)))},
               cpu=round(_cpu() - t0, 2), s=float(st2.p.s), count=int(st2.p.n), frozen=int(st2.frozen.sum()), info=info)
    return st2, row


FEAS_PEN = 1e-9          # feasible = max SAT penetration at most this (the store's bar: pen <= 1e-9)
TRACK_MAX_N = 3000       # Packing.max_pen builds an n x n distance matrix; larger packings are not tracked


def _track(best, p, n, name):
    """best (s, Packing, stage) or None, updated with packing p if it has n squares, is smaller and feasible"""
    if p.n != n or p.n > TRACK_MAX_N or (best is not None and p.s >= best[0]):
        return best
    return (float(p.s), p, name) if p.max_pen() <= FEAS_PEN else best


def failed(log):
    """did this run_chain call fail to produce a packing?  (a 'keep.best' row means it returned an earlier feasible one)"""
    if log and log[-1]['stage'] == 'keep.best':
        return False
    return any(r['info'].get('failed') for r in log)


def run_chain(spec, ctx, seed, st=None, keep_best=True):
    """run every stage of spec in order from st (None: the first stage must be a gen) -> (State, log rows).
    Tracks the best feasible packing at ctx['n'] after every stage (State.best, carried across calls; the start state's
    own packing is not a candidate).  keep_best: if the final state is not that best (worse, infeasible, wrong count) or
    a stage failed, return the best instead, with a 'keep.best' log row.  keep_best=False (the build half of a split
    chain, which is relaxed next) only tracks."""
    rng = random.Random(seed)
    log = []
    n = int(ctx['n']) if ctx.get('n') is not None else None
    best = st.best if st is not None else None
    for name, kw in (parse(spec) if isinstance(spec, str) else spec):
        st, row = run_stage(name, kw, st, rng, ctx)
        log.append(row)
        if n is not None:
            best = _track(best, st.p, n, name)
            st = st.but(best=best)
    if keep_best and best is not None and (st.p is not best[1] or any(r['info'].get('failed') for r in log)):
        why = ('failed' if any(r['info'].get('failed') for r in log) else 'count' if st.p.n != n
               else 'worse' if st.p.s >= best[0] else 'infeasible')
        final_s = float(st.p.s)
        st = st.but(p=best[1], frozen=np.zeros(best[1].n, bool), region=None)
        log.append(dict(stage='keep.best', kind='keep', params={}, cpu=0.0, s=best[0], count=int(best[1].n), frozen=0,
                        info=dict(src=best[2], why=why, final_s=final_s)))
    return st, log


def deterministic(spec):
    """does the chain give the same output for every seed?  True only if every stage declares it (attribute
    `deterministic`: bool, or a function of the fixed parameters) and no parameter is left to draw"""
    for name, kw in (parse(spec) if isinstance(spec, str) else spec):
        f = STAGES[name]; d = getattr(f, 'deterministic', False)
        if any(is_space(kw.get(k, v)) if k in kw else isinstance(v, tuple) for k, v in f.space.items()):
            return False
        if not (d(dict({k: v for k, v in f.space.items()}, **kw)) if callable(d) else d):
            return False
    return True


def table():
    """(name, kind, space, first doc line) for every stage"""
    return [(k, f.kind, f.space, ((f.__doc__ or '').strip().splitlines() or [''])[0]) for k, f in sorted(STAGES.items())]


def names(spec):
    """stage names in a chain string or a move-log params key list (unknown tokens dropped; no parse errors)"""
    out = set()
    for part in str(spec or '').split('|'):
        m = re.match(r'\s*([A-Za-z_][\w.]*)', part)
        if m and m.group(1) in STAGES:
            out.add(m.group(1))
    return out
