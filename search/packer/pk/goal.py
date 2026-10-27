"""Goals: targeted exploration as one loop (10-24; tasks/explore-loop/DESIGN.md).

A goal is data (a YAML file); one scheduler spends CPU on arms; every arm is a stage chain run on donors selected from
the store; every output goes back to the store; reward is measured from the store.

  pk.py goal plan  NAME|FILE     arms and their donor pool sizes, nothing run
  pk.py goal run   NAME|FILE     the loop (launch it with pk/jobs.py; it runs `procs` jobs at a time)
  pk.py goal report NAME|FILE    per arm: jobs, CPU, reward, posterior; per target: frontier, matches

Where goal files live (10-24): a goal's arms, results and predictions are research results, so goal files belong to a
plugin (`<plugin dir>/goals/NAME.yaml`, found by name); this repo keeps only `goals/example.yaml`.  Run output goes to
runs/goals/NAME/ here either way (gitignored).

Goal file (every key but `goal`, `targets`, `predictions` has a default):
  goal: name                      targets: [n, ...]
  progress: {exclude: [G, D, S], within: 0.02}   a match = a new packing at a target whose familymap code is not
                                                 excluded, within `within` of the best such packing at job start
                                                 (family-relative near-record; Evan 10-24, the 148 case)
  value: {match: 1.0, frontier: 1.0e5, bar: 100}  reward = match * new clusters + frontier * (family-best gain + store-
                                                 best gain), side units
  bars: {n: side}                                a bonus (value.bar) for the job that first gets below side at n
  budget_cpu_h: 15   quantum_cpu_h: auto (budget / (procs x 8), >= 0.05)   procs: 4   relax: round.15c.p5   donor_window: 0.05   donors: 6
  arms: auto | [ {kind: walk|chain, offset: 0, chain: "...", match_donors: true}, ... ]
  arm_kinds: [walk, add]         (optional) with arms auto, keep only these kinds
  predictions: {...}  (required: the scheduler refuses a goal without them)   baseline: [...]   notes: "..."

Arms (`auto`, per target n in row k = ceil sqrt n): walk (donors at n), add (n - 1, grow.add), drop (n + 1, n + 2,
shrink.drop), lcut (n + 2k + 1 + e, e = 0..2, shrink.lcut), and every plugin template (`ARM_TEMPLATES` in a plugin's
pk_plugin.py: f(goal, n) -> [arm dict]).  Construction arms spend `construct_share` of the quantum left after scoring (1 - score_share) building starts
(chain + relax), then walk the best four: they are judged on refined outcomes like everything else.
"""
from __future__ import annotations
import json, math, os, random, re, resource, sqlite3, subprocess, sys, tempfile, time
from dataclasses import dataclass, field
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # search/packer
RUNS = os.path.join(HERE, 'runs', 'goals')
WAVES = 8                                                         # auto quantum: the run is about this many job lengths
FEAS = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9) AND (grid IS NULL OR grid=0)"


# ------------------------------------------------------------------ goal spec

@dataclass
class Goal:
    goal: str
    targets: list
    predictions: dict
    progress: dict = field(default_factory=lambda: dict(exclude=['G', 'D', 'S'], within=0.02))
    value: dict = field(default_factory=lambda: dict(match=1.0, frontier=1.0e5, bar=100.0))
    retire: bool = True                                         # stop scheduling a target once its bar is crossed
    ref_cap: int = 3000                                         # reference: at most this many pre-goal packings per n
    bars: dict = field(default_factory=dict)                    # {n: side}: a bonus (value.bar) for the job that first gets below
    budget_cpu_h: float = 15.0
    quantum_cpu_h: object = 'auto'                              # 'auto': budget / (procs x WAVES)
    procs: int = 4
    relax: str = 'round.15c.p5'
    donor_window: float = 0.05
    donors: int = 6
    score_share: float = 0.25                                   # of the quantum, kept for converge + polish
    construct_share: float = 0.3
    polish_reuse: bool = True                                   # reuse trusted polish results (polish_cached)
    clusters: bool = True                                       # False (10-25): no cluster-novelty reward; the reference
                                                                # polishes only the REF_TOP best pre-goal packings (f0p)
    polish_top: object = None                                   # finalize polishes only this many best candidates (None: all
                                                                # in the window; default 3 when clusters is False)
    screen: object = None                                       # pk/screen.py: {k, kmax, ratio, audit} or true
    relax_cap: float = 4.0                                      # relax stage timeout: x the job's median relax
    arms: object = 'auto'
    arm_kinds: object = None                                    # auto arms: keep only these kinds (walk, add, drop, lcut, plugin)
    queue: object = None                                        # unit loop: priority-queue claims over an n range (pk/queue.py)
    baseline: list = field(default_factory=list)
    notes: str = ''
    result: str = ''                 # written after the run (the reading of the result)
    path: str = ''

    @property
    def dir(self):
        return os.path.join(RUNS, self.goal)


def goal_dirs():
    """where goals are found by name: each loaded plugin's goals/, then this repo's goals/ (the example)"""
    from . import stages as S
    return [os.path.join(os.path.dirname(os.path.abspath(m.__file__)), 'goals') for m in S.PLUGINS] + \
        [os.path.join(HERE, 'goals')]


def resolve(spec):
    """a goal file path, or a goal name looked up in goal_dirs()"""
    if os.path.exists(spec):
        return os.path.abspath(spec)
    name = os.path.basename(spec)
    name = name if name.endswith(('.yaml', '.yml')) else name + '.yaml'
    for d in goal_dirs():
        if os.path.exists(os.path.join(d, name)):
            return os.path.join(d, name)
    raise FileNotFoundError(f'goal {spec}: not a file, and not in {goal_dirs()}')


def load(path):
    import yaml
    path = resolve(path)
    d = yaml.safe_load(open(path))
    g = Goal(**{k: v for k, v in d.items() if k in Goal.__dataclass_fields__}, path=os.path.abspath(path))
    unknown = set(d) - set(Goal.__dataclass_fields__)
    if unknown:
        raise ValueError(f'{path}: unknown keys {sorted(unknown)}')
    if not g.predictions:
        raise ValueError(f'{path}: no predictions (pre-register before running)')
    num = lambda x: float(x) if isinstance(x, (int, float, str)) and not isinstance(x, bool) else x
    g.value = {k: num(v) for k, v in g.value.items()}           # PyYAML reads 1.0e4 (no exponent sign) as a string
    g.progress = {k: (num(v) if k != 'exclude' else v) for k, v in g.progress.items()}
    if g.quantum_cpu_h == 'auto':                                 # 10-24: two waves of 0.5 CPU-h left one straggler
        g.quantum_cpu_h = max(0.05, float(g.budget_cpu_h) / (int(g.procs) * WAVES))   # alone for 30 min, 15 cores idle
    for f in ('budget_cpu_h', 'quantum_cpu_h', 'donor_window', 'construct_share', 'score_share', 'relax_cap'):
        setattr(g, f, float(getattr(g, f)))
    g.procs, g.donors = int(g.procs), int(g.donors)
    g.bars = {int(k): float(v) for k, v in (g.bars or {}).items()}
    return g


def kgrid(n):
    return math.ceil(math.sqrt(n) - 1e-12)


# ------------------------------------------------------------------ selector

def family(n, s, sq):
    sys.path.insert(0, HERE)
    from familymap import classify
    return classify(n, s, sq)['code']


def best_match(db, goal, n, before='9999', limit=20000):
    """lowest side at n whose family the goal wants (Evan 10-24: near-record within the wanted families, not to the overall
    frontier, which at 148 is a diamond 0.06 below every non-diamond packing) -> (s, id) or (None, None)"""
    k = kgrid(n); ex = goal.progress.get('exclude', []); seen = set()
    for i, s, c in db.execute(f"SELECT id, s, coords FROM packing WHERE n=? AND {FEAS} AND s<? AND added<? ORDER BY s LIMIT ?",
                              (n, k + 0.05, before, limit)):
        if round(s, 6) in seen:
            continue
        seen.add(round(s, 6))
        if family(n, s, np.frombuffer(c, '<f8').reshape(-1, 3)) not in ex:
            return s, i
    return None, None


def donors(goal, n, m, rng, match=True, window=None):
    """up to m store packings at n.  match: the goal's families only, within `donor_window` of the best such packing
    (family-relative near-record); spread over the window by side rank, best first.  Not match: within the window of
    the overall best.  -> [(id, s, sq)]"""
    window = goal.donor_window if window is None else window
    db = sqlite3.connect(f'file:{os.path.join(HERE, "runs/store.sqlite")}?mode=ro', uri=True)
    k = kgrid(n); ex = goal.progress.get('exclude', []) if match else []
    rows, seen, lo = [], set(), None
    for i, s, c in db.execute(f"SELECT id, s, coords FROM packing WHERE n=? AND {FEAS} AND s<? ORDER BY s LIMIT 50000", (n, k + 0.05)):
        if lo is not None and s > lo + window:
            break
        if round(s, 6) in seen:
            continue
        seen.add(round(s, 6))
        sq = np.frombuffer(c, '<f8').reshape(-1, 3)
        if ex and family(n, s, sq) in ex:
            continue
        lo = s if lo is None else lo
        rows.append((i, s, sq.copy()))
        if len(rows) >= 400:
            break
    if len(rows) <= m:
        return rows
    picks = [rows[0]] + rng.sample(rows[1:], m - 1)         # the family best, then a spread of the window
    return sorted(picks, key=lambda r: r[1])


# ------------------------------------------------------------------ planner

def plan(goal):
    """arms: dicts with id, n, offset (donor n = n + offset), chain (None = walk), match_donors"""
    if goal.arms != 'auto':
        A = []
        for n in goal.targets:
            for a in goal.arms:
                A.append(dict(a, n=n, id=f"{n}:{a.get('name', a['kind'])}:{a.get('offset', 0)}"))
        return A
    R = goal.relax; A = []
    for n in goal.targets:
        k = kgrid(n)
        A.append(dict(id=f'{n}:walk', n=n, offset=0, chain=None, match_donors=True))
        A.append(dict(id=f'{n}:add', n=n, offset=-1, chain=f'grow.add | relax.{R}', match_donors=True))
        for o in (1, 2):
            A.append(dict(id=f'{n}:drop+{o}', n=n, offset=o, chain=f'shrink.drop | relax.{R}', match_donors=True))
        for e in (0, 1, 2):
            A.append(dict(id=f'{n}:lcut+{2 * k + 1 + e}', n=n, offset=2 * k + 1 + e, chain=f'shrink.lcut | relax.{R}',
                          match_donors=True))
        from . import stages as S
        for m in S.PLUGINS:
            for tpl in getattr(m, 'ARM_TEMPLATES', []):
                for a in tpl(goal, n):
                    A.append(dict(dict(match_donors=True, chain=None), **a, n=n))
    if goal.arm_kinds:
        A = [a for a in A if arm_kind(a['id']) in goal.arm_kinds]
    return A


def arm_kind(aid):
    """'335:drop+1' -> 'drop'"""
    return re.match(r'[a-z_]*', aid.split(':')[1]).group(0)


def _pool_size(args):
    goal, a, seed = args
    if str(a.get('chain') or '').strip().startswith('gen.'):
        return 1
    return len(donors(goal, a['n'] + a['offset'], 3, random.Random(seed), a.get('match_donors', True)))


def pool_sizes(goal, arms, seed=0, prev=None, procs=1):
    """{arm id: donors available (<= 3)}.  One store scan per arm (0.1-1 s each; 560 arms ~10 min serial, 10-25), so:
    procs > 1 runs the scans in a process pool, and with prev only the arms that had no donors are rechecked (the store
    only grows, so a nonempty pool stays nonempty)."""
    out = dict(prev or {})
    todo = [a for a in arms if not out.get(a['id'])]
    args = [(goal, a, seed + i) for i, a in enumerate(todo)]
    if procs > 1 and len(todo) > procs:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(procs) as ex:
            v = list(ex.map(_pool_size, args, chunksize=4))
    else:
        v = [_pool_size(x) for x in args]
    out.update({a['id']: k for a, k in zip(todo, v)})
    return out


# ------------------------------------------------------------------ job (runs in its own process)

def _cpu():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def _write(path, s, sq):
    with open(path, 'w') as fh:
        fh.write(f'{len(sq)} {float(s)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in sq))


def score_share(goal, n):
    """the walk leaves this fraction of the quantum for ingest + converge + polish: `score_share` until two jobs at n
    have finished, then 1.3 x their median measured share (clamped to 0.03-0.6)"""
    jl = os.path.join(goal.dir, 'jobs.jsonl')
    f = []
    for l in open(jl) if os.path.exists(jl) else []:
        r = json.loads(l); c = r.get('cpu_split')
        if r.get('n') == n and c and sum(c.values()) > 0:
            f.append((sum(c.values()) - c.get('walk', 0)) / sum(c.values()))
    if len(f) < 2:
        return goal.score_share
    return min(0.6, max(0.03, 1.3 * sorted(f)[len(f) // 2]))


def run_job(goal, arm, out, cpu_s, seed, finish=True):
    """one arm x one quantum: donors from the store, construct (chain arms), walk; files only (the scheduler imports).
    finish=False (pk/units.py): stop after the walk; the caller ingests and polishes as units of their own"""
    sys.path.insert(0, HERE)
    import mcmin
    from layout import full_lines
    from . import stages as S
    from .packing import Packing
    rng = random.Random(seed); c0 = _cpu()
    n = arm['n']; k = kgrid(n)
    os.makedirs(os.path.join(out, 'starts'), exist_ok=True)
    gen_led = str(arm.get('chain') or '').strip().startswith('gen.')
    D = [(None, None, None)] if gen_led else donors(goal, n + arm['offset'], goal.donors, rng, arm.get('match_donors', True))
    log = dict(arm=arm['id'], donors=[d[0] for d in D], built=0, kept=0)
    parents = {}
    if not D:
        json.dump(dict(log, empty=True, cpu=_cpu() - c0), open(os.path.join(out, 'job.json'), 'w')); return
    from . import screen as SCR
    if arm.get('chain') and goal.screen and SCR.split_chain(arm['chain']):
        # 10-25: screen-then-promote (pk/screen.py): build many, screen each (a few compress iterations), fully relax
        # the best by screen score (+ a random audit share); every record goes to job.json for backtests
        rt = []

        def cap():
            return max(30.0, goal.relax_cap * float(np.median(rt))) if len(rt) >= 3 else max(60.0, goal.construct_share * cpu_s * (1 - score_share(goal, n)) / 4)

        def on_relax(lg):
            rt.extend(r['cpu'] for r in lg if r['kind'] == 'relax' and not r['info'].get('failed'))

        def pick_start():
            pid, s, sq = D[rng.randrange(len(D))]
            return pid, (None if gen_led else S.State(Packing(s, sq)))

        def keep(st, lg):
            if st.p.n != n or st.p.s >= k + 0.05:
                return False
            with tempfile.TemporaryDirectory() as td:
                p = os.path.join(td, 'e.txt'); mcmin.write_deg(p, float(st.p.s), [tuple(map(float, q)) for q in st.p.sq])
                return not full_lines(p, k)
        kept, recs = SCR.construct(arm['chain'], dict(n=n, cap=cap, on_relax=on_relax), goal.construct_share * cpu_s * (1 - score_share(goal, n)),
                                   _cpu, rng, pick_start, keep, goal.screen if isinstance(goal.screen, dict) else None, log)
        log['built'] = sum(1 for r in recs if 'traj' in r)
        json.dump(recs, open(os.path.join(out, 'screen.json'), 'w'))
    if arm.get('chain') and not (goal.screen and SCR.split_chain(arm['chain'])):
        kept = []; rt = []                                        # rt: CPU of each finished relax stage
        while _cpu() - c0 < goal.construct_share * cpu_s * (1 - score_share(goal, n)):
            pid, s, sq = D[rng.randrange(len(D))]
            # 10-25: relax cap.  At n ~ 1000 a few relaxes run 10-70x the median and end at the grid (978: 8 of 60
            # took 63 % of the relax time, none useful); cap at relax_cap x the job's running median (first: 1/4 of
            # the construct budget)
            cap = max(30.0, goal.relax_cap * float(np.median(rt))) if len(rt) >= 3 else max(60.0, goal.construct_share * cpu_s * (1 - score_share(goal, n)) / 4)
            try:
                st, lg = S.run_chain(arm['chain'], dict(n=n, fq_timeout=cap), rng.randrange(1 << 30),
                                     st=None if gen_led else S.State(Packing(s, sq)))
            except Exception as e:
                log.setdefault('errors', []).append(repr(e)[:120]); continue
            log['built'] += 1
            rt += [r['cpu'] for r in lg if r['kind'] == 'relax' and not r['info'].get('failed')]
            log['relax_timeouts'] = log.get('relax_timeouts', 0) + sum(r['info'].get('failed') == 'timeout' for r in lg)
            if st.p.n != n or S.failed(lg) or st.p.s >= k + 0.05:
                continue
            with tempfile.TemporaryDirectory() as td:
                p = os.path.join(td, 'e.txt'); mcmin.write_deg(p, float(st.p.s), [tuple(map(float, q)) for q in st.p.sq])
                if full_lines(p, k):
                    continue
            kept.append((float(st.p.s), pid if pid is not None else st.meta.get('parent'),   # gen-led: the stage's
                         np.asarray(st.p.sq, float)))                               # meta names its parent
    if arm.get('chain'):
        kept.sort(key=lambda t: t[0]); seen = set(); files = []
        for s, pid, sq in kept:
            if round(s, 6) in seen:
                continue
            seen.add(round(s, 6)); f = os.path.join(out, 'starts', f'p{pid}_{len(files)}.txt'); _write(f, s, sq)
            if pid is not None:
                parents[os.path.basename(f)] = pid
            files.append(f)
        log['kept'] = len(files); walk_from = files[:4]
    else:
        walk_from = []
        for pid, s, sq in D[:4]:
            f = os.path.join(out, 'starts', f'id{pid}.txt'); _write(f, s, sq); parents[os.path.basename(f)] = pid
            walk_from.append(f)
    left = cpu_s * (1 - score_share(goal, n)) - (_cpu() - c0)
    if walk_from and left > 30:
        with open(os.path.join(out, 'walk.log'), 'w') as fh:
            subprocess.run([sys.executable, 'explore.py', '--n', str(n), '--starts', *walk_from, '--procs', '1',
                            '--minutes', f'{left / 60:.2f}', '--out', os.path.join(out, 'ex'), '--seed', str(seed % 100000),
                            '--polish-gate', '2e-3', '--weights', 'melt=2,jmelt=1,jrefill=1'],
                           cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
    json.dump(dict(log, parents=parents, cpu=_cpu() - c0), open(os.path.join(out, 'job.json'), 'w'))
    if finish:
        finish_job(goal, arm, out, cpu_s, c0, time.strftime('%Y-%m-%dT%H:%M:%S',
                                                            time.localtime(os.path.getmtime(os.path.join(out, 'arm.json')))))


def finish_job(goal, arm, out, cpu_s, c0, t0):
    """10-24 (Evan: score incrementally, inside the compute budget): the job imports its own output, converges it and
    polishes its candidates, all on its own CPU clock; job.json gets the timing split and finalized = True"""
    from .store import Store
    from . import converge
    jid = os.path.basename(out.rstrip('/')); T = {'walk': _cpu() - c0}; t = _cpu()
    st = Store(); ref_cpu = 0.0
    try:
        J = ingest_job(st, goal, jid, out); T['ingest'] = _cpu() - t; t = _cpu()
        ids = _job_ids(st.db, goal, jid, t0)
        converge.run(st, lo=arm['n'], hi=arm['n'], procs=1, ids=ids, log=lambda *a: None); T['converge'] = _cpu() - t; t = _cpu()
        reference(st, goal, arm['n'], 1); ref_cpu = _cpu() - t; t = _cpu()   # 10-25: on demand (0 when cached); the goal's
        F = finalize(st, goal, jid, t0, arm['n'], 1, max(60.0, cpu_s - (_cpu() - c0 - ref_cpu))); T['polish'] = _cpu() - t
    finally:                                                                 # budget counts it via ref_n*.json, not the job
        st.close()
    json.dump(dict(J, cpu=_cpu() - c0 - ref_cpu, ref_cpu=round(ref_cpu), cpu_split={k: round(v) for k, v in T.items()},
                   finalized=True, **F), open(os.path.join(out, 'job.json'), 'w'))


# ------------------------------------------------------------------ store side: import, reward

def ingest_job(st, goal, jid, out):
    """starts (with donor parents) and the walk's archive into the store; walk roots get their donor as parent"""
    from .packing import read
    J = json.load(open(os.path.join(out, 'job.json')))
    src = f'run:runs/goals/{goal.goal}/{jid}'
    P = J.get('parents', {})
    for f, pid in P.items():
        path = os.path.join(out, 'starts', f)
        if os.path.exists(path):
            st.add(read(path), src + '/starts', 'ours', ref=f, status='screened', parent=pid, commit=False)
    for f in sorted(os.listdir(os.path.join(out, 'starts'))) if os.path.isdir(os.path.join(out, 'starts')) else []:
        if f not in P:                                            # starts from gen-led arms (no donor)
            st.add(read(os.path.join(out, 'starts', f)), src + '/starts', 'ours', ref=f, status='screened', commit=False)
    st.db.commit()
    ex = os.path.join(out, 'ex')
    if os.path.isdir(ex):
        subprocess.run([sys.executable, 'pk.py', 'import-run', ex], cwd=HERE, capture_output=True)
        arch = os.path.join(ex, 'archive.jsonl')
        if os.path.exists(arch):
            for line in open(arch):
                e = json.loads(line)
                if str(e.get('kind', '')).startswith('start:'):
                    pid = P.get(e['kind'][6:])
                    if pid:
                        st.db.execute("UPDATE packing SET parent=? WHERE parent IS NULL AND id IN (SELECT packing FROM sighting "
                                      "WHERE source LIKE ? AND ref=?)", (pid, f'run:{os.path.relpath(ex, HERE)}%', str(e['i'])))
            st.db.commit()
    return J


def best_before(db, n, t):
    k = kgrid(n)
    r = db.execute(f"SELECT min(s) FROM packing WHERE n=? AND {FEAS} AND s<? AND added<?", (n, k - 1e-9, t)).fetchone()[0]
    return r if r is not None else float(k)


def _polish_one(a):
    i, s, sq = a
    sys.path.insert(0, HERE)
    import hop
    from .engine import POLISH
    with tempfile.TemporaryDirectory(prefix='gp_') as td:
        try:
            r = hop.quench(s, [tuple(map(float, q)) for q in sq], td, opt=POLISH)
        except Exception:
            r = None
    return i, ((float(r[0]), np.asarray(r[1], float)) if r else None)


def polish_many(items, procs):
    """[(id, s, sq)] -> {id: (s, sq) or None}: fq POLISH (no ALM, loosen 1.0), the finalize step before counting"""
    if not items:
        return {}
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max(1, min(procs, len(items)))) as ex:
        return dict(ex.map(_polish_one, items, chunksize=max(1, len(items) // (4 * procs))))   # few items: one per worker


# sources whose packings are fq POLISH outputs of polish_many (10-25): re-polishing them is a no-op (median 2e-13,
# max 7.5e-7 over a sample); explorer 'polished' packings are not (up to 2.4e-4 more on a re-polish), so not trusted
POLISH_SOURCES = ("goal-ref:%", "run:runs/goals/%/polish")


FIXED = 'polish-fixed'


def polish_fix(items, procs, rounds=3, tol=1e-10):
    """polish until a round changes the side by <= tol (at most `rounds`; Evan 10-25: stop once re-polishing stops
    changing things).  -> {id: (s, sq, fixed) or None}; fixed = the last round moved <= tol (a polish fixpoint)"""
    out, todo = {}, list(items)
    for r in range(rounds):
        P = polish_many(todo, procs); nxt = []
        for i, s, sq in todo:
            v = P.get(i)
            if v is None:
                out.setdefault(i, None); continue
            moved = s - v[0]
            out[i] = (v[0], v[1], abs(moved) <= tol)
            if abs(moved) > tol:
                nxt.append((i, v[0], v[1]))
        todo = nxt
        if not todo:
            break
    return out


def polish_cached(db, items, procs, reuse=True):
    """polish with reuse: a packing that is a polish fixpoint (tag polish-fixed) or a trusted polish output, or has
    such a polished child (the best one), is not polished again; the rest are polished to a fixpoint (polish_fix) and
    tagged when fixed.  reuse=False: one plain polish_many round (the pre-10-25 behaviour).
    -> ({id: (s, sq[, fixed]) or None}, dict(self=, child=, polished=))"""
    out, stats = {}, dict(self=0, child=0, polished=0)
    if not reuse:
        P = polish_many(items, procs); stats['polished'] = len(items)
        return P, stats
    if items:
        ids = [i for i, _, _ in items]; q = ','.join('?' * len(ids))
        tr = ' OR '.join('g.source LIKE ?' for _ in POLISH_SOURCES)
        selfp = {i for (i,) in db.execute(f"SELECT DISTINCT g.packing FROM sighting g WHERE g.packing IN ({q}) AND ({tr})",
                                         (*ids, *POLISH_SOURCES))}
        selfp |= {i for (i,) in db.execute(f"SELECT id FROM packing WHERE id IN ({q}) AND tags LIKE ?", (*ids, f'%"{FIXED}"%'))}
        for i, s, sq in items:
            if i in selfp:
                out[i] = (s, sq); stats['self'] += 1
        rest = [i for i in ids if i not in out]
        if rest:
            q = ','.join('?' * len(rest))
            for par, s, c in db.execute(f"SELECT p.parent, p.s, p.coords FROM packing p LEFT JOIN sighting g ON g.packing = p.id "
                                        f"WHERE p.parent IN ({q}) AND (({tr}) OR p.tags LIKE ?) AND {FEAS} ORDER BY p.s DESC",
                                        (*rest, *POLISH_SOURCES, f'%"{FIXED}"%')):
                out[par] = (s, np.frombuffer(c, '<f8').reshape(-1, 3).copy())       # ascending last: the best child wins
            stats['child'] = sum(1 for i in rest if i in out)
    todo = [x for x in items if x[0] not in out]
    out.update(polish_fix(todo, procs)); stats['polished'] = len(todo)
    return out, stats


def window(goal, ref_side):
    """the match window: `within_rel` x side if set (comparable across n), else the absolute `within`"""
    r = goal.progress.get('within_rel')
    return float(r) * ref_side if r is not None else float(goal.progress.get('within', 0.02))


def _add_polished(st, src, P, parents):
    """polished packings into the store (parent = the packing polished); -> {original id: polished store id}"""
    from .packing import Packing
    out = {}
    for i, v in P.items():
        if v is None:
            continue
        pid, _ = st.add(Packing(v[0], v[1]), src, 'ours', ref=str(i), status='polished', parent=parents.get(i, i), commit=False)
        if len(v) > 2 and v[2]:                                   # a polish fixpoint: never re-polished (polish_cached)
            st.tag(pid, FIXED)
        out[i] = pid
    st.db.commit()
    return out


def _wanted(goal, n, s, sq, k):
    return s < k - 1e-9 and family(n, s, sq) not in goal.progress.get('exclude', [])


def reference(st, goal, n, procs=4):
    """once per (goal, n): the polished wanted-family minima in the window before the goal started; they seed the
    cluster leaders (tag goal_<name>), so a job is paid only for clusters nobody had.  Cached in the goal dir."""
    import fcntl
    cf = os.path.join(goal.dir, f'ref_n{n}.json')
    if os.path.exists(cf):
        return json.load(open(cf))
    with open(cf[:-5] + '.lock', 'w') as lk:                     # 10-25: built on demand by the first job to need
        fcntl.flock(lk, fcntl.LOCK_EX)                            # it; concurrent jobs at n wait, then reuse it
        if os.path.exists(cf):
            return json.load(open(cf))
        return _build_reference(st, goal, n, procs, cf)


def goal_start(goal):
    """the goal's start: start.json (written by run() before any job, 10-25), else the old inference from finished
    jobs and job-directory times (goals run before 10-25)"""
    sf = os.path.join(goal.dir, 'start.json')
    if os.path.exists(sf):
        return json.load(open(sf))['t_goal']
    jl = os.path.join(goal.dir, 'jobs.jsonl')
    rows = [json.loads(l) for l in open(jl)] if os.path.exists(jl) else []
    dirs = [d for d in os.listdir(goal.dir) if d.startswith('j')] if os.path.isdir(goal.dir) else []
    return min([r['t0'] for r in rows] + [time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(os.path.getmtime(os.path.join(goal.dir, d))))
                                         for d in dirs] + [time.strftime('%Y-%m-%dT%H:%M:%S')])


REF_TOP = 3                                                       # clusters off: pre-goal packings polished per n


def _build_reference(st, goal, n, procs, cf):
    from . import cluster as CL
    db = st.db; k = kgrid(n); c0 = _cpu()                         # cpu_h covers selection + classification too
    t_goal = goal_start(goal)
    f0, _ = best_match(db, goal, n, before=t_goal)
    f0 = f0 if f0 is not None else float(k)
    w = window(goal, f0)
    items = []
    for i, s, c in db.execute(f"SELECT id, s, coords FROM packing WHERE n=? AND {FEAS} AND s<=? AND added<? ORDER BY s LIMIT ?",
                              (n, min(f0 + 1.5 * w, k - 1e-9), t_goal, int(goal.ref_cap))):
        sq = np.frombuffer(c, '<f8').reshape(-1, 3).copy()
        if _wanted(goal, n, s, sq, k):
            items.append((i, s, sq))
    if not goal.clusters:                                         # 10-25: no novelty reward -> only f0p is needed
        items = items[:REF_TOP]
    P, pst = polish_cached(db, items, procs, goal.polish_reuse)
    ids = _add_polished(st, f'goal-ref:{goal.goal}', P, {})
    keep = [(ids[i], v[0], v[1]) for i, v in P.items() if v and i in ids and _wanted(goal, n, v[0], v[1], k)]
    f0p = min([x[1] for x in keep] + [f0])
    keep = [x for x in keep if x[1] <= f0p + window(goal, f0p)]
    for f in (os.path.join(CL.CACHE, f'n{n}_goal_{goal.goal}.json'),):
        if os.path.exists(f):
            os.remove(f)
    L = CL.labels(n, keep, tag=f'goal_{goal.goal}') if goal.clusters else {}
    R = dict(t_goal=t_goal, f0=f0, f0p=f0p, window=window(goal, f0p), minima=len({round(x[1], 8) for x in keep}),
             clusters=len(set(L.values())) if goal.clusters else None, cpu_h=round((_cpu() - c0) / 3600, 3), polish=pst)   # cpu: in the budget
    json.dump(R, open(cf + '.tmp', 'w'), indent=1); os.replace(cf + '.tmp', cf)   # readers never see a partial file
    return R


def _job_ids(db, goal, jid, t0):
    """this job's packings: its sources (not its polish) + converge children added since t0"""
    like = f'run:runs/goals/{goal.goal}/{jid}%'
    ids = {i for (i,) in db.execute('SELECT DISTINCT packing FROM sighting WHERE source LIKE ? AND source NOT LIKE ?',
                                    (like, like.rstrip('%') + '/polish%'))}
    if ids:
        q = ','.join('?' * len(ids))
        ids |= {i for (i,) in db.execute(f"SELECT id FROM packing WHERE parent IN ({q}) AND added>=? AND id NOT IN "
                                         f"(SELECT packing FROM sighting WHERE source LIKE '%/polish')", (*ids, t0))}
    return ids


def finalize(st, goal, jid, t0, n, procs=1, cpu_left=None):
    """the expensive half of scoring (10-24: runs inside the job, on its CPU budget): fq POLISH this job's packings at
    target n inside 1.5 x the window, best first, until `cpu_left` CPU-seconds are spent; the polished packings enter
    the store under .../<jid>/polish.  -> dict(cands, polished, jbest)"""
    db = st.db; k = kgrid(n)
    ref = reference(st, goal, n, procs)
    f0p, w = ref['f0p'], ref['window']
    items, jbest = [], float('inf')
    for i in _job_ids(db, goal, jid, t0):
        r = db.execute(f"SELECT n, s, coords FROM packing WHERE id=? AND {FEAS}", (i,)).fetchone()
        if not r or r[0] != n:
            continue
        jbest = min(jbest, r[1])
        sq = np.frombuffer(r[2], '<f8').reshape(-1, 3).copy()
        if r[1] <= f0p + 1.5 * w and _wanted(goal, n, r[1], sq, k):
            items.append((i, r[1], sq))
    items.sort(key=lambda t: t[1]); c0 = _cpu(); done = 0; batch = max(4, 4 * procs)
    ncand = len(items)
    top = goal.polish_top if goal.polish_top is not None else (None if goal.clusters else 3)
    if top is not None:                                           # 10-25: gain needs only the best few polished
        items = items[:int(top)]
    while done < len(items) and (cpu_left is None or _cpu() - c0 < cpu_left):
        P, _ = polish_cached(db, items[done:done + batch], procs, goal.polish_reuse); done += len(P)
        _add_polished(st, f'run:runs/goals/{goal.goal}/{jid}/polish', P, {})
    return dict(cands=ncand, polished=done, jbest=jbest if jbest < float('inf') else None)


def count(st, goal, jid, t0, n, jbest=None, fprev=None):
    """the cheap half (the scheduler, one job at a time so cluster leaders never race): this job's polished packings,
    deduped at 1e-8, counted as NEW region clusters (pk/cluster leaders, tag goal_<name>) that no pre-goal minimum and
    no earlier job had, within the window of the polished family best at goal start; gain = drop of the polished
    family best; store_gain = drop of the overall store best"""
    from . import cluster as CL
    db = st.db; k = kgrid(n)
    ref = reference(st, goal, n)
    f0p, w = ref['f0p'], ref['window']
    b0 = best_before(db, n, t0)
    jbest = float('inf') if jbest is None else jbest
    keep, seen, raw = [], set(), 0
    for pid, s, c in db.execute(f"SELECT DISTINCT p.id, p.s, p.coords FROM packing p JOIN sighting g ON g.packing = p.id "
                                f"WHERE g.source = ? AND p.n = ?", (f'run:runs/goals/{goal.goal}/{jid}/polish', n)):
        jbest = min(jbest, s)
        sq = np.frombuffer(c, '<f8').reshape(-1, 3).copy()
        if s > f0p + w or not _wanted(goal, n, s, sq, k) or round(s, 8) in seen:
            continue
        seen.add(round(s, 8)); keep.append((pid, s, sq)); raw += 1
    L = CL.labels(n, keep, tag=f'goal_{goal.goal}') if keep and goal.clusters else {}
    mine = {x[0] for x in keep}
    new_clusters = {lab for lab in L.values() if int(lab.split(':')[1]) in mine}
    fbest = min([x[1] for x in keep] + [f0p])
    best = min([b0] + [x[1] for x in keep if x[1] < k - 1e-9] + ([jbest] if jbest < k - 1e-9 else []))
    fref = f0p if fprev is None else min(f0p, fprev)               # 10-24: vs the family best already recorded, so a
    gain, sgain = max(0.0, fref - fbest), max(0.0, b0 - best)      # second job reaching the same side is not paid again
    bar = goal.bars.get(n)                                        # 10-24: 1e4 x a 0.003 record gain paid ~30, less
    crossed = int(bar is not None and best < bar <= b0)           # than a busy cluster arm; gains now 1e5, bars a bonus
    v = goal.value
    return dict(matches=len(new_clusters), minima=len(keep), raw=raw, gain=gain, store_gain=sgain, fbest=fbest,
                best_gap=(jbest - b0) if jbest < float('inf') else None,
                bar=crossed, reward=v.get('match', 1.0) * len(new_clusters) + v.get('frontier', 1e5) * (gain + sgain)
                + v.get('bar', 100.0) * crossed)


def reward(st, goal, jid, t0, n, procs=4, fprev=None):
    """v2 (Evan 10-24): finalize (polish everything) + count; for rescore and for jobs that did not finalize"""
    F = finalize(st, goal, jid, t0, n, procs)
    return count(st, goal, jid, t0, n, F['jbest'], fprev)


def rescore(goal, procs=12, log=print):
    """recompute every recorded job's reward with the current reward(), in job order (novelty is chronological); the
    previous values are kept as v1_* fields; jobs.jsonl is rewritten (backup: jobs.v1.jsonl)"""
    from .store import Store
    from . import cluster as CL
    jl = os.path.join(goal.dir, 'jobs.jsonl')
    rows = [json.loads(l) for l in open(jl)]
    if not os.path.exists(jl.replace('.jsonl', '.v1.jsonl')):
        open(jl.replace('.jsonl', '.v1.jsonl'), 'w').write(''.join(json.dumps(r) + '\n' for r in rows))
    for f in [os.path.join(goal.dir, x) for x in os.listdir(goal.dir) if x.startswith('ref_n') and x.endswith('.json')] + \
             [os.path.join(CL.CACHE, f'n{n}_goal_{goal.goal}.json') for n in goal.targets]:
        if os.path.exists(f):
            os.remove(f)
    st = Store()
    try:
        fb = {}
        for r in sorted(rows, key=lambda r: r['job']):
            R = reward(st, goal, r['job'], r['t0'], r['n'], procs, fb.get(r['n']))
            fb[r['n']] = min(fb.get(r['n'], R['fbest']), R['fbest'])
            for kk in ('matches', 'reward', 'gain', 'store_gain'):
                r.setdefault('v1_' + kk, r.get(kk))
            r.update(R)
            log(f"{r['job']} {r['arm']}: v1 matches {r['v1_matches']} -> new clusters {R['matches']} (minima {R['minima']}, raw {R['raw']})")
    finally:
        st.close()
    open(jl, 'w').write(''.join(json.dumps(r) + '\n' for r in sorted(rows, key=lambda r: r['job'])))


# ------------------------------------------------------------------ scheduler

def arm_family(arm_id):
    """'179:walk' -> 'walk'; '17:disk.p20:0' -> 'disk.p20'; '145:lcut+27' -> 'lcut' (the shared-prior group)"""
    import re
    return re.sub(r'[-+]\d+$', '', arm_id.split(':')[1])


class Posterior:
    """Gamma-Poisson reward rate per CPU-h.  Prior mean per arm = its family's pooled rate across targets (cross
    construction that pays at 148 starts high at 147; 10-24), shrunk to the goal-wide rate when the family has little
    CPU; prior weight a0 = 0.5 pseudo-events either way"""
    def __init__(self, rows):
        self.r, self.c, fr, fc = {}, {}, {}, {}
        for x in rows:
            a, f = x['arm'], arm_family(x['arm'])
            self.r[a] = self.r.get(a, 0.0) + x['reward']; self.c[a] = self.c.get(a, 0.0) + x['cpu_h']
            fr[f] = fr.get(f, 0.0) + x['reward']; fc[f] = fc.get(f, 0.0) + x['cpu_h']
        tot_c = sum(self.c.values()); self.mean = (sum(self.r.values()) / tot_c) if tot_c > 0.5 else 2.0
        self.fam = {f: (fr[f] + self.mean * 1.0) / (fc[f] + 1.0) for f in fc}     # 1 CPU-h of the goal-wide rate as shrinkage

    def sample(self, arm, rng, pending=0.0):
        m = self.fam.get(arm_family(arm), self.mean)
        a0, b0 = 0.5, 0.5 / max(m, 1e-3)
        a = a0 + self.r.get(arm, 0.0) + m * pending; b = b0 + self.c.get(arm, 0.0) + pending
        return rng.gammavariate(a, 1.0 / b)


def run(goal, seed=1, log=print):
    from .store import Store
    os.makedirs(goal.dir, exist_ok=True)
    jl = os.path.join(goal.dir, 'jobs.jsonl')
    rows = [json.loads(l) for l in open(jl)] if os.path.exists(jl) else []
    rng = random.Random(seed + len(rows))
    arms = plan(goal); A = {a['id']: a for a in arms}
    sf = os.path.join(goal.dir, 'start.json')                    # 10-25: references are built on demand by the
    if not os.path.exists(sf):                                    # first job to finish at each n (no serial setup:
        json.dump(dict(t_goal=goal_start(goal)), open(sf, 'w'))   # ~1 h of a mostly idle machine for 40 targets)
    refs = lambda: sum(json.load(open(os.path.join(goal.dir, f))).get('cpu_h', 0)
                       for f in os.listdir(goal.dir) if f.startswith('ref_n') and f.endswith('.json'))
    # 10-25 accounting: the dispatcher's own CPU (planning, pool scans incl. their worker processes, count / cluster
    # labels, recovery) is in the budget too; run.json keeps one segment per run() call (wall + setup CPU)
    af = os.path.join(goal.dir, 'run.json')
    prev = json.load(open(af))['segments'] if os.path.exists(af) else []
    seg = dict(start=time.strftime('%Y-%m-%dT%H:%M:%S'), job=os.environ.get('PK_JOB_NAME'), wall_h=0.0, setup_cpu_h=0.0)
    w0 = time.time(); self0 = sum(resource.getrusage(resource.RUSAGE_SELF)[:2]); pool_cpu = [0.0]
    def setup_cpu():
        return (sum(resource.getrusage(resource.RUSAGE_SELF)[:2]) - self0 + pool_cpu[0]) / 3600
    def acct(final=False):
        seg.update(wall_h=round((time.time() - w0) / 3600, 4), setup_cpu_h=round(setup_cpu(), 4),
                   **(dict(end=time.strftime('%Y-%m-%dT%H:%M:%S')) if final else {}))
        json.dump(dict(segments=prev + [seg]), open(af + '.tmp', 'w'), indent=1); os.replace(af + '.tmp', af)
    overhead = lambda: refs() + sum(x['setup_cpu_h'] for x in prev) + setup_cpu()
    def eta():
        left = max(0.0, goal.budget_cpu_h - sum(r['cpu_h'] for r in rows) - overhead())
        return f"ETA {time.strftime('%H:%M', time.localtime(time.time() + 3600 * (left / goal.procs + 0.5 * goal.quantum_cpu_h)))}"
    log(f'goal {goal.goal}: {len(arms)} arms over {goal.targets}; budget {goal.budget_cpu_h} CPU-h, spent '
        f'{sum(r["cpu_h"] for r in rows) + overhead():.2f} (references {refs():.2f}); quantum {goal.quantum_cpu_h:.3f} '
        f'CPU-h x {goal.procs} procs; {eta()}')

    def score(name, out, arm, t0, recovered=False):
        """count (cheap) for a job that finalized itself; a job from before 10-24 or one that crashed: finish it here"""
        from . import converge
        J = json.load(open(os.path.join(out, 'job.json')))
        st = Store()
        try:
            if not J.get('finalized'):
                J = ingest_job(st, goal, name, out)
                converge.run(st, lo=arm['n'], hi=arm['n'], procs=goal.procs, ids=_job_ids(st.db, goal, name, t0),
                             log=lambda *a: None)
                J.update(finalize(st, goal, name, t0, arm['n'], goal.procs))
            fp = [r['fbest'] for r in rows if r['n'] == arm['n'] and r.get('fbest') is not None]
            R = count(st, goal, name, t0, arm['n'], J.get('jbest'), min(fp) if fp else None)
        finally:
            st.close()
        row = dict(job=name, arm=arm['id'], n=arm['n'], t0=t0, cpu_h=J.get('cpu', 0) / 3600,
                   donors=len(J.get('donors', [])), built=J.get('built', 0), kept=J.get('kept', 0),
                   cands=J.get('cands'), polished=J.get('polished'), cpu_split=J.get('cpu_split'),
                   **(dict(recovered=True) if recovered else {}), **R)
        rows.append(row); open(jl, 'a').write(json.dumps(row) + '\n')
        log(f"{name} {'recovered' if recovered else 'done'} {row['arm']}: {row['cpu_h']:.2f} CPU-h {row['cpu_split']}, "
            f"polished {row['polished']} / {row['cands']}, new clusters {R['matches']}, gain {R['gain']:.2e}, "
            f"reward {R['reward']:.2f}; {eta()}")

    done = {r['job'] for r in rows}                               # resume: finished jobs never recorded (a crash)
    for name in sorted(d for d in os.listdir(goal.dir) if d.startswith('j') and d not in done):
        out = os.path.join(goal.dir, name)
        if not os.path.exists(os.path.join(out, 'job.json')) or not os.path.exists(os.path.join(out, 'arm.json')):
            continue
        arm = json.load(open(os.path.join(out, 'arm.json')))
        t0 = time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(os.path.getmtime(os.path.join(out, 'arm.json'))))
        score(name, out, arm, t0, True)
    retired = set()
    if goal.retire and goal.bars:                                 # 10-24 (Evan): a fast-follow target stops once its
        st = Store()                                              # bar is crossed (or already is at the start)
        try:
            retired = {n for n, b in goal.bars.items() if best_before(st.db, n, '9999') < b}
        finally:
            st.close()
        if retired:
            log(f'retired at start (already below the bar): {sorted(retired)}')
    c = sum(resource.getrusage(resource.RUSAGE_CHILDREN)[:2])
    sizes = pool_sizes(goal, arms, procs=goal.procs); live = {}       # its worker processes: setup CPU
    pool_cpu[0] += sum(resource.getrusage(resource.RUSAGE_CHILDREN)[:2]) - c
    acct(); last_acct = time.time()
    jid = max([int(d[1:]) for d in os.listdir(goal.dir) if d.startswith('j') and d[1:].isdigit()] + [0])
    last_sizes = time.time()
    while True:
        spent = sum(r['cpu_h'] for r in rows) + overhead() + goal.quantum_cpu_h * len(live)
        if time.time() - last_acct > 60:
            acct(); last_acct = time.time()
        if time.time() - last_sizes > 600:
            sizes = pool_sizes(goal, arms, seed=jid, prev=sizes); last_sizes = time.time()
        while len(live) < goal.procs and spent + goal.quantum_cpu_h <= goal.budget_cpu_h + 1e-9:
            post = Posterior(rows)
            pend = {}
            for x in live.values():
                pend[x['arm']] = pend.get(x['arm'], 0.0) + goal.quantum_cpu_h
            retired |= {r['n'] for r in rows if goal.retire and r.get('bar')}
            elig = [a for a in arms if sizes.get(a['id'], 0) > 0 and a['n'] not in retired]
            if not elig:
                if not getattr(run, '_noarm', False):
                    log('no arm left (no donors, or every target retired)'); run._noarm = True
                break
            pick = max(elig, key=lambda a: post.sample(a['id'], rng, pend.get(a['id'], 0.0)))
            jid += 1; name = f'j{jid:04d}'; out = os.path.join(goal.dir, name); os.makedirs(out, exist_ok=True)
            t0 = time.strftime('%Y-%m-%dT%H:%M:%S')
            json.dump(pick, open(os.path.join(out, 'arm.json'), 'w'))
            p = subprocess.Popen([sys.executable, 'pk.py', 'goal', 'job', goal.path, out, str(goal.quantum_cpu_h * 3600),
                                  str(seed * 100003 + jid)], cwd=HERE, stdout=open(os.path.join(out, 'job.log'), 'w'),
                                 stderr=subprocess.STDOUT)
            live[name] = dict(p=p, arm=pick['id'], t0=t0, out=out)
            spent += goal.quantum_cpu_h
            log(f'{name} start {pick["id"]} (pool {sizes[pick["id"]]})')
        if not live:
            break
        time.sleep(5)
        for name in [x for x, v in live.items() if v['p'].poll() is not None]:
            v = live.pop(name)
            score(name, v['out'], A[v['arm']], v['t0'])
    acct(final=True)
    log(report(goal))
    acct(final=True)                                              # again: the report's own CPU / wall (~minutes)


def accounting(goal, rows):
    """one line: wall and CPU of the whole run (10-25): jobs + references + dispatcher setup (run.json), and the
    jobs.py wrapper's own measurement (every waited-for process, plus the sink + converge pass after the run)"""
    refs = sum(json.load(open(os.path.join(goal.dir, f))).get('cpu_h', 0)
               for f in os.listdir(goal.dir) if f.startswith('ref_n') and f.endswith('.json'))
    jobs = sum(r['cpu_h'] for r in rows)
    af = os.path.join(goal.dir, 'run.json')
    S = json.load(open(af))['segments'] if os.path.exists(af) else []
    setup = sum(x['setup_cpu_h'] for x in S); wall = sum(x['wall_h'] for x in S)
    tot = jobs + refs + setup
    line = (f'Accounting: CPU {tot:.2f} h = jobs {jobs:.2f} + references {refs:.2f} + dispatcher / setup {setup:.2f}; '
            f'wall {wall:.2f} h' + (f' (utilization {tot / (wall * goal.procs):.0%} of {goal.procs} procs)' if wall else
                                    ' (no run.json: a run before 10-25)'))
    names = [x['job'] for x in S if x.get('job')]
    if names:
        from . import jobs as J
        ac = J.accounting(names[-1])
        if ac and ac.get('cpu_h') is not None:
            line += (f'.  jobs.py {names[-1]}: wall {ac["wall_h"]:.2f} h, CPU {ac["cpu_h"]:.2f} h (run {ac["main_cpu_h"]:.2f}, '
                     f'sink + converge {ac["sink_cpu_h"]:.2f} in {ac["sink_wall_h"]:.2f} h)')
    return line


def report(goal):
    jl = os.path.join(goal.dir, 'jobs.jsonl')
    rows = [json.loads(l) for l in open(jl)] if os.path.exists(jl) else []
    post = Posterior(rows); by = {}
    for r in rows:
        by.setdefault(r['arm'], []).append(r)
    out = [f'# goal {goal.goal}: {len(rows)} jobs, {sum(r["cpu_h"] for r in rows):.2f} CPU-h of {goal.budget_cpu_h}',
           '', accounting(goal, rows), '', '| arm | jobs | CPU-h | new clusters (v1 matches) | new minima | gain: family / store | best gap to store best | reward / CPU-h |',
           '|---|---|---|---|---|---|---|---|']
    for a, R in sorted(by.items(), key=lambda t: -sum(r['reward'] for r in t[1])):
        c = sum(r['cpu_h'] for r in R); w = sum(r['reward'] for r in R)
        g = [r['best_gap'] for r in R if r.get('best_gap') is not None]
        v1 = sum(r.get('v1_matches', r['matches']) or 0 for r in R)
        out.append(f"| {a} | {len(R)} | {c:.2f} | {sum(r['matches'] for r in R)} ({v1}) | {sum(r.get('minima', 0) for r in R)} | "
                   f"{max(r['gain'] for r in R):.2e} / {max(r.get('store_gain', 0) for r in R):.2e} | {min(g) if g else float('nan'):+.2e} | {w / max(c, 1e-9):.2f} |")
    db = sqlite3.connect(f'file:{os.path.join(HERE, "runs/store.sqlite")}?mode=ro', uri=True)
    out += ['', '| target | store best now | at goal start | wanted-family best now | at goal start |', '|---|---|---|---|---|']
    t0 = min((r['t0'] for r in rows), default=time.strftime('%Y-%m-%dT%H:%M:%S'))
    fm = lambda x: f'{x:.9f}' if x is not None else '-'
    for n in goal.targets:
        out.append(f'| {n} | {best_before(db, n, "9999"):.9f} | {best_before(db, n, t0):.9f} | '
                   f'{fm(best_match(db, goal, n)[0])} | {fm(best_match(db, goal, n, before=t0)[0])} |')
    out += ['', 'Library: distinct (1e-6) wanted-family packings within +d of the wanted-family best at goal start '
            '(start -> now; the reference stays fixed so the counts compare)', '',
            '| target | +0.001 | +0.005 | +0.02 |', '|---|---|---|---|']
    for n in goal.targets:
        f0 = best_match(db, goal, n, before=t0)[0]
        if f0 is None:
            out.append(f'| {n} | - | - | - |'); continue
        a, b = library(db, goal, n, f0, t0), library(db, goal, n, f0, '9999')
        out.append(f'| {n} | ' + ' | '.join(f'{a[d]} -> {b[d]}' for d in (0.001, 0.005, 0.02)) + ' |')
    return '\n'.join(out)


def library(db, goal, n, ref, before, ds=(0.001, 0.005, 0.02)):
    """distinct (1e-6) wanted-family packings at n with side <= ref + d, added before `before` -> {d: count}"""
    ex = goal.progress.get('exclude', []); seen = set(); keep = []
    for s, c in db.execute(f"SELECT s, coords FROM packing WHERE n=? AND {FEAS} AND s<=? AND added<? ORDER BY s",
                           (n, ref + max(ds), before)):
        if round(s, 6) in seen:
            continue
        seen.add(round(s, 6))
        if family(n, s, np.frombuffer(c, '<f8').reshape(-1, 3)) not in ex:
            keep.append(s)
    return {d: sum(1 for s in keep if s <= ref + d) for d in ds}
