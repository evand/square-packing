"""Goals: targeted exploration as one loop (10-24; tasks/explore-loop/DESIGN.md).

A goal is data (a YAML file); one scheduler spends CPU on arms; every arm is a stage chain run on donors selected from
the store; every output goes back to the store; reward is measured from the store.

  pk.py goal plan  goals/fam145-150.yaml      arms and their donor pool sizes, nothing run
  pk.py goal run   goals/fam145-150.yaml      the loop (launch it with pk/jobs.py; it runs `procs` jobs at a time)
  pk.py goal report fam145-150                per arm: jobs, CPU, reward, posterior; per target: frontier, matches

Goal file (every key but `goal`, `targets`, `predictions` has a default):
  goal: name                      targets: [n, ...]
  progress: {exclude: [G, D, S], within: 0.02}   a match = a new packing at a target whose familymap code is not
                                                 excluded, within `within` of the best such packing at job start
                                                 (family-relative near-record; Evan 10-24, the 148 case)
  value: {match: 1.0, frontier: 1.0e4}           reward = match * new matches + frontier * (family-best gain + store-
                                                 best gain), side units
  budget_cpu_h: 15   quantum_cpu_h: 0.5   procs: 4   relax: round.15c.p5   donor_window: 0.05   donors: 6
  arms: auto | [ {kind: walk|chain, offset: 0, chain: "...", match_donors: true}, ... ]
  predictions: {...}  (required: the scheduler refuses a goal without them)   baseline: [...]   notes: "..."

Arms (`auto`, per target n in row k = ceil sqrt n): walk (donors at n), add (n - 1, grow.add), drop (n + 1, n + 2,
shrink.drop), lcut (n + 2k + 1 + e, e = 0..2, shrink.lcut), and every plugin template (`ARM_TEMPLATES` in a plugin's
pk_plugin.py: f(goal, n) -> [arm dict]).  Construction arms spend `construct_share` of the quantum building starts
(chain + relax), then walk the best four: they are judged on refined outcomes like everything else.
"""
from __future__ import annotations
import json, math, os, random, resource, sqlite3, subprocess, sys, tempfile, time
from dataclasses import dataclass, field
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # search/packer
RUNS = os.path.join(HERE, 'runs', 'goals')
FEAS = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9) AND (grid IS NULL OR grid=0)"


# ------------------------------------------------------------------ goal spec

@dataclass
class Goal:
    goal: str
    targets: list
    predictions: dict
    progress: dict = field(default_factory=lambda: dict(exclude=['G', 'D', 'S'], within=0.02))
    value: dict = field(default_factory=lambda: dict(match=1.0, frontier=1.0e4))
    budget_cpu_h: float = 15.0
    quantum_cpu_h: float = 0.5
    procs: int = 4
    relax: str = 'round.15c.p5'
    donor_window: float = 0.05
    donors: int = 6
    construct_share: float = 0.3
    arms: object = 'auto'
    baseline: list = field(default_factory=list)
    notes: str = ''
    path: str = ''

    @property
    def dir(self):
        return os.path.join(RUNS, self.goal)


def load(path):
    import yaml
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
    for f in ('budget_cpu_h', 'quantum_cpu_h', 'donor_window', 'construct_share'):
        setattr(g, f, float(getattr(g, f)))
    g.procs, g.donors = int(g.procs), int(g.donors)
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
    return A


def pool_sizes(goal, arms, seed=0):
    rng = random.Random(seed); out = {}
    for a in arms:
        if str(a.get('chain') or '').strip().startswith('gen.'):
            out[a['id']] = 1; continue
        out[a['id']] = len(donors(goal, a['n'] + a['offset'], 3, rng, a.get('match_donors', True)))
    return out


# ------------------------------------------------------------------ job (runs in its own process)

def _cpu():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def _write(path, s, sq):
    with open(path, 'w') as fh:
        fh.write(f'{len(sq)} {float(s)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in sq))


def run_job(goal, arm, out, cpu_s, seed):
    """one arm x one quantum: donors from the store, construct (chain arms), walk; files only (the scheduler imports)"""
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
    if arm.get('chain'):
        kept = []
        while _cpu() - c0 < goal.construct_share * cpu_s:
            pid, s, sq = D[rng.randrange(len(D))]
            try:
                st, lg = S.run_chain(arm['chain'], dict(n=n), rng.randrange(1 << 30),
                                     st=None if gen_led else S.State(Packing(s, sq)))
            except Exception as e:
                log.setdefault('errors', []).append(repr(e)[:120]); continue
            log['built'] += 1
            if st.p.n != n or any(r['info'].get('failed') for r in lg) or st.p.s >= k + 0.05:
                continue
            with tempfile.TemporaryDirectory() as td:
                p = os.path.join(td, 'e.txt'); mcmin.write_deg(p, float(st.p.s), [tuple(map(float, q)) for q in st.p.sq])
                if full_lines(p, k):
                    continue
            kept.append((float(st.p.s), pid, np.asarray(st.p.sq, float)))
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
    left = cpu_s - (_cpu() - c0)
    if walk_from and left > 30:
        with open(os.path.join(out, 'walk.log'), 'w') as fh:
            subprocess.run([sys.executable, 'explore.py', '--n', str(n), '--starts', *walk_from, '--procs', '1',
                            '--minutes', f'{left / 60:.2f}', '--out', os.path.join(out, 'ex'), '--seed', str(seed % 100000),
                            '--polish-gate', '2e-3', '--weights', 'melt=2,jmelt=1,jrefill=1'],
                           cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
    json.dump(dict(log, parents=parents, cpu=_cpu() - c0), open(os.path.join(out, 'job.json'), 'w'))


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


def reward(st, goal, jid, t0, n):
    """from the store: this job's packings (its sources + converge children of them) at target n added since t0.
    matches = new distinct (1e-6) wanted-family packings within `within` of the best wanted-family packing at t0;
    gain = how much that family best dropped; store_gain = how much the overall store best dropped (both side units)"""
    db = st.db
    like = f'run:runs/goals/{goal.goal}/{jid}%'
    ids = {i for (i,) in db.execute('SELECT DISTINCT packing FROM sighting WHERE source LIKE ?', (like,))}
    if ids:
        q = ','.join('?' * len(ids))
        ids |= {i for (i,) in db.execute(f"SELECT id FROM packing WHERE parent IN ({q}) AND added>=?", (*ids, t0))}
    b0 = best_before(db, n, t0); f0, _ = best_match(db, goal, n, before=t0)
    k = kgrid(n)
    if f0 is None:                         # no wanted-family packing at n yet: any new one below the grid is a match
        f0 = float(k) - goal.progress.get('within', 0.02)
    within = goal.progress.get('within', 0.02); ex = goal.progress.get('exclude', [])
    new_keys, best, fbest, jbest = set(), b0, f0, float('inf')
    for i in ids:
        r = db.execute(f"SELECT n, s, coords FROM packing WHERE id=? AND {FEAS}", (i,)).fetchone()
        if not r or r[0] != n or r[1] >= k + 0.05:
            continue
        s = r[1]; jbest = min(jbest, s)
        if s < k - 1e-9:
            best = min(best, s)
        if family(n, s, np.frombuffer(r[2], '<f8').reshape(-1, 3)) in ex:
            continue
        fbest = min(fbest, s)
        if s > f0 + within or s >= k - 1e-9 or round(s, 6) in new_keys:
            continue
        if db.execute(f"SELECT 1 FROM packing WHERE n=? AND {FEAS} AND abs(s-?)<5e-7 AND added<? LIMIT 1", (n, s, t0)).fetchone():
            continue
        new_keys.add(round(s, 6))
    gain, sgain = max(0.0, f0 - fbest), max(0.0, b0 - best)
    v = goal.value
    return dict(matches=len(new_keys), gain=gain, store_gain=sgain, best_gap=(jbest - b0) if jbest < float('inf') else None,
                reward=v.get('match', 1.0) * len(new_keys) + v.get('frontier', 1e4) * (gain + sgain))


# ------------------------------------------------------------------ scheduler

class Posterior:
    """Gamma-Poisson reward rate per CPU-h with a shared empirical prior (all arms' pooled rate)"""
    def __init__(self, rows):
        self.r, self.c = {}, {}
        for x in rows:
            self.r[x['arm']] = self.r.get(x['arm'], 0.0) + x['reward']; self.c[x['arm']] = self.c.get(x['arm'], 0.0) + x['cpu_h']
        tot_c = sum(self.c.values()); self.mean = (sum(self.r.values()) / tot_c) if tot_c > 0.5 else 2.0

    def sample(self, arm, rng, pending=0.0):
        a0, b0 = 0.5, 0.5 / max(self.mean, 1e-3)
        a = a0 + self.r.get(arm, 0.0) + self.mean * pending; b = b0 + self.c.get(arm, 0.0) + pending
        return rng.gammavariate(a, 1.0 / b)


def run(goal, seed=1, log=print):
    from .store import Store
    from . import converge
    os.makedirs(goal.dir, exist_ok=True)
    jl = os.path.join(goal.dir, 'jobs.jsonl')
    rows = [json.loads(l) for l in open(jl)] if os.path.exists(jl) else []
    rng = random.Random(seed + len(rows))
    arms = plan(goal); A = {a['id']: a for a in arms}
    log(f'goal {goal.goal}: {len(arms)} arms over {goal.targets}; budget {goal.budget_cpu_h} CPU-h, spent '
        f'{sum(r["cpu_h"] for r in rows):.2f}')
    done = {r['job'] for r in rows}                               # resume: finished jobs never recorded (a crash)
    for name in sorted(d for d in os.listdir(goal.dir) if d.startswith('j') and d not in done):
        out = os.path.join(goal.dir, name)
        if not os.path.exists(os.path.join(out, 'job.json')) or not os.path.exists(os.path.join(out, 'arm.json')):
            continue
        arm = json.load(open(os.path.join(out, 'arm.json')))
        t0 = time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(os.path.getmtime(os.path.join(out, 'arm.json'))))
        st = Store()
        try:
            J = ingest_job(st, goal, name, out)
            converge.run(st, since=t0, procs=1, log=lambda *a: None)
            R = reward(st, goal, name, t0, arm['n'])
        finally:
            st.close()
        row = dict(job=name, arm=arm['id'], n=arm['n'], t0=t0, cpu_h=J.get('cpu', 0) / 3600, donors=len(J.get('donors', [])),
                   built=J.get('built', 0), kept=J.get('kept', 0), recovered=True, **R)
        rows.append(row); open(jl, 'a').write(json.dumps(row) + '\n')
        log(f"{name} recovered {arm['id']}: {row['cpu_h']:.2f} CPU-h, matches {R['matches']}, reward {R['reward']:.2f}")
    sizes = pool_sizes(goal, arms); live = {}
    jid = max([int(d[1:]) for d in os.listdir(goal.dir) if d.startswith('j') and d[1:].isdigit()] + [0])
    last_sizes = time.time()
    while True:
        spent = sum(r['cpu_h'] for r in rows) + goal.quantum_cpu_h * len(live)
        if time.time() - last_sizes > 600:
            sizes = pool_sizes(goal, arms, seed=jid); last_sizes = time.time()
        while len(live) < goal.procs and spent + goal.quantum_cpu_h <= goal.budget_cpu_h + 1e-9:
            post = Posterior(rows)
            pend = {}
            for x in live.values():
                pend[x['arm']] = pend.get(x['arm'], 0.0) + goal.quantum_cpu_h
            elig = [a for a in arms if sizes.get(a['id'], 0) > 0]
            if not elig:
                log('no arm has donors'); break
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
            st = Store()
            try:
                J = ingest_job(st, goal, name, v['out'])
                converge.run(st, since=v['t0'], procs=1, log=lambda *a: None)
                R = reward(st, goal, name, v['t0'], A[v['arm']]['n'])
            finally:
                st.close()
            row = dict(job=name, arm=v['arm'], n=A[v['arm']]['n'], t0=v['t0'], cpu_h=J.get('cpu', 0) / 3600,
                       donors=len(J.get('donors', [])), built=J.get('built', 0), kept=J.get('kept', 0), **R)
            rows.append(row)
            open(jl, 'a').write(json.dumps(row) + '\n')
            log(f"{name} done {v['arm']}: {row['cpu_h']:.2f} CPU-h, matches {R['matches']}, gain {R['gain']:.2e}, reward {R['reward']:.2f}")
    log(report(goal))


def report(goal):
    jl = os.path.join(goal.dir, 'jobs.jsonl')
    rows = [json.loads(l) for l in open(jl)] if os.path.exists(jl) else []
    post = Posterior(rows); by = {}
    for r in rows:
        by.setdefault(r['arm'], []).append(r)
    out = [f'# goal {goal.goal}: {len(rows)} jobs, {sum(r["cpu_h"] for r in rows):.2f} CPU-h of {goal.budget_cpu_h}',
           '', '| arm | jobs | CPU-h | matches | gain: family / store | best gap to store best | reward / CPU-h |', '|---|---|---|---|---|---|---|']
    for a, R in sorted(by.items(), key=lambda t: -sum(r['reward'] for r in t[1])):
        c = sum(r['cpu_h'] for r in R); w = sum(r['reward'] for r in R)
        g = [r['best_gap'] for r in R if r.get('best_gap') is not None]
        out.append(f"| {a} | {len(R)} | {c:.2f} | {sum(r['matches'] for r in R)} | {max(r['gain'] for r in R):.2e} / {max(r.get('store_gain', 0) for r in R):.2e} | {min(g) if g else float('nan'):+.2e} | {w / max(c, 1e-9):.2f} |")
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
