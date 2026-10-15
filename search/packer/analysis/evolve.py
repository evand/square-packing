#!/usr/bin/env python3
"""LLM program evolution for square-packing constructors (pilot, 10-15; family J = Xu #432's LLM-evolved programs).

AlphaEvolve-style: an LLM writes `construct(n, rng)` programs that emit structured STARTING layouts; our fixed backend
(fq ALM squeeze + SLP polish) does the numerics.  Proposals come from headless Claude Code (`claude -p`, subscription
tokens), run isolated: empty cwd outside ~/math, no tools, no settings / CLAUDE.md / MCP, no session persistence; the
prompt carries only the system text (evolve_system.md), parent programs and their evaluation reports.

Time travel: targets = frontier_at(n, --cutoff) (default 2026-10-08T00:00Z, just before #432); the LLM never sees
anything later.  Training n (fitness) = --train; `report` evaluates the best programs on --holdout too and only then
compares with the posted values (Xu's), which the loop never reads.

Fitness per program = mean over train n of g = (s_best - target) / (grid - target), grid = ceil(sqrt n): 0 ties the
target, 1 = grid, < 0 beats it.  Population: --islands islands; a proposal picks an island round-robin, a parent by
tournament (3) inside it and an inspiration (best of another island).  Models: haiku / sonnet by --mix; every
--director proposals an opus call sees the top programs and is asked for a substantially new design (new island
member at the worst island).

  python3 analysis/evolve.py run --out runs/evolve/p1 --proposals 300 [--llm-par 6 --procs 16]
  python3 analysis/evolve.py report --out runs/evolve/p1 [--top 5]
"""
import argparse, concurrent.futures as cf, json, math, os, random, re, shutil, subprocess, sys, tempfile, threading, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK)
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
CUTOFF = '2026-10-08T00:00'
TRAIN = '51,84,103,105,126,127,146,175'
HOLDOUT = '70,86,102,108,123,129,131,261,267,295'
LLM_CWD = '/tmp/claude-1000/evolve-cwd'
BANNED = re.compile(r'\b(os|sys|subprocess|socket|shutil|pathlib|importlib|ctypes|multiprocessing|threading|signal|'
                    r'builtins|urllib|requests|http|pickle|marshal|glob|io)\b|__import__|\bopen\s*\(|\beval\s*\(|\bexec\s*\(|'
                    r'\bcompile\s*\(|__builtins__|__subclasses__|__globals__')
IMPORT_OK = {'math', 'numpy', 'random', 'itertools', 'functools'}


# ------------------------------------------------------------------ evaluation (worker processes)

_ENG = None


def static_check(code):
    for m in re.finditer(r'^\s*(?:from\s+([\w.]+)\s+import|import\s+([\w., ]+?)(?:\s+as\s+\w+)?\s*$)', code, re.M):
        mods = [x.strip().split('.')[0].split(' as ')[0] for x in (m.group(1) or m.group(2)).split(',')]
        bad = [x for x in mods if x not in IMPORT_OK]
        if bad:
            return f'forbidden import {bad}'
    m = BANNED.search(code)
    return f'forbidden token {m.group(0)!r}' if m else None


def describe(p):
    """short text description of an optimised layout for the LLM"""
    a = p.sq[:, 2] % 90.0
    t = np.minimum(a, 90 - a) > 1.0
    s = p.s
    out = [f'{int(t.sum())} tilted']
    if t.any():
        ang = np.round(np.where(a[t] > 45, a[t] - 90, a[t])).astype(int)      # signed tilt in (-45, 45]
        vals, cnt = np.unique(ang, return_counts=True)
        top = sorted(zip(cnt, vals), reverse=True)[:4]
        out.append('angles ' + ', '.join(f'{v} deg x{c}' for c, v in top))
        x, y = p.sq[t, 0], p.sq[t, 1]
        out.append(f'tilted region x {x.min():.1f}-{x.max():.1f}, y {y.min():.1f}-{y.max():.1f} of {s:.3f}')
    return '; '.join(out)


def ascii_map(p):
    """1 char per half unit: '#' axis-aligned centre, '/' tilted centre, '.' empty"""
    k = int(math.ceil(p.s * 2))
    G = [['.'] * k for _ in range(k)]
    for x, y, a in p.sq:
        a %= 90.0
        i, j = min(k - 1, int(x * 2)), min(k - 1, int(y * 2))
        G[k - 1 - j][i] = '/' if min(a, 90 - a) > 1.0 else '#'
    return '\n'.join(''.join(r) for r in G)


def eval_job(args):
    """(program path, n, seed, maxc, loosens) -> dict(n, s, sq, desc, err, cands)"""
    global _ENG
    prog, n, seed, maxc, loosens = args[:5]
    fill = len(args) > 5 and args[5]
    from pk.engine import FQ, QOpt, POLISH, SCREEN
    from pk.packing import Packing
    if _ENG is None:
        _ENG = FQ()
    td = tempfile.mkdtemp(prefix='evo_')
    try:
        r = subprocess.run(['unshare', '-rn', sys.executable, '-I', os.path.join(HERE, 'evolve_runner.py'), prog, str(n),
                            str(seed), str(maxc), '1' if fill else '0'], cwd=td, capture_output=True, text=True, timeout=60)
        line = (r.stdout.strip().splitlines() or [''])[-1]
        try:
            res = json.loads(line)
        except Exception:
            return dict(n=n, err=f'runner failed: {(r.stderr or r.stdout)[-300:]}')
    except subprocess.TimeoutExpired:
        return dict(n=n, err='timeout (60 s wall)')
    finally:
        shutil.rmtree(td, ignore_errors=True)
    if 'error' in res:
        return dict(n=n, err=res['error'])
    best = None
    fills = []
    mins = []
    tf = tq = 0.0
    impl = os.environ.get('PK_FILL_IMPL', 'fast')
    g_ = os.environ.get('PK_EVAL_GATE', '2e-3')           # polish only screens within this of the best (C42 gating)
    gate = None if g_ == 'off' else float(g_)
    screens = []

    def _keep(pp):
        nonlocal best
        if pp is None or pp.max_pen() > 1e-7:
            return
        if pp.ntilted() >= 5 and all(abs(pp.s - m[0]) > 1e-9 for m in mins):   # distinct non-grid minima
            mins.append((float(pp.s), pp.sq.tolist()))
        if best is None or pp.s < best.s:
            best = pp        # FP1 (10-15): fillfast = same counts, in-process, no fork
    for s, sq in res['cands']:
        if fill:                                          # the given squares stay; the exact MILP adds uprights
            import mondrian as M
            import fillfast
            T = np.array(sq, float).reshape(-1, 3)
            got = None
            t0 = time.time()
            for f in (1.0, 1.01, 1.02, 1.04):
                if impl == 'fast':
                    L = fillfast.fill(s * f, T * [f, f, 1], 'exact')
                else:
                    L = M.exact_fill(T * [f, f, 1], s * f, time_limit=10)
                if f == 1.0:
                    fills.append(len(L))
                if len(L) >= n:
                    got = (s * f, M.trim(L, n)); break
            tf += time.time() - t0
            if got is None:
                continue
            s, sq = got
        p = Packing(float(s), np.array(sq, float))
        t0 = time.time()
        if gate is None:                                  # ungated (gen 0 / gen 1 behaviour)
            for L in loosens:
                q = _ENG.quench(p, QOpt(loosen=L), timeout=60)
                if q.p is None:
                    continue
                q2 = _ENG.quench(q.p, POLISH, timeout=60)
                _keep(q2.p if q2.p is not None else q.p)
        else:                                             # gated (FP1): screen every start, polish only the near-best
            for L in loosens:
                q = _ENG.quench(p, SCREEN.but(loosen=L), timeout=60)
                if q.p is not None:
                    screens.append(q.p)
        tq += time.time() - t0
    if gate is not None and screens:
        t0 = time.time()
        b = min(x.s for x in screens)
        for x in sorted(screens, key=lambda x: x.s):
            if x.s > b + gate:
                break
            q2 = _ENG.quench(x, POLISH, timeout=60)
            _keep(q2.p if q2.p is not None else x)
        tq += time.time() - t0
    fx = f'; exact fill at the given s reached {max(fills)} / {n} squares (best candidate)' if fills else ''
    if best is None:
        return dict(n=n, err='no candidate optimised to a feasible packing' + fx, cands=len(res['cands']), t_fill=tf, t_quench=tq)
    return dict(n=n, s=float(best.s), sq=best.sq.tolist(), desc=describe(best) + fx, cands=len(res['cands']),
                minima=sorted(mins)[:40], t_fill=tf, t_quench=tq)


# ------------------------------------------------------------------ LLM proposals (threads)

def call_llm(model, prompt, system, timeout=600):
    os.makedirs(LLM_CWD, exist_ok=True)
    cmd = ['claude', '-p', '--model', model, '--tools', '', '--setting-sources', '', '--strict-mcp-config',
           '--no-session-persistence', '--output-format', 'json', '--system-prompt', system]
    env = {k: v for k, v in os.environ.items() if not k.startswith('CLAUDE_CODE_') and k not in ('CLAUDECODE',)}
    t0 = time.time()
    try:
        r = subprocess.run(cmd, input=prompt, capture_output=True, text=True, cwd=LLM_CWD, timeout=timeout, env=env)
        d = json.loads(r.stdout)
    except Exception as e:
        return dict(err=f'llm call failed: {e!r}'[:300], sec=time.time() - t0)
    text = d.get('result') or ''
    mu = d.get('modelUsage') or {}
    tok = sum(v.get('inputTokens', 0) + v.get('outputTokens', 0) + v.get('cacheReadInputTokens', 0) +
              v.get('cacheCreationInputTokens', 0) for v in mu.values())
    m = re.findall(r'```python\s*\n(.*?)```', text, re.S)
    return dict(text=text, code=m[-1] if m else None, usd=d.get('total_cost_usd'), tokens=tok, sec=time.time() - t0,
                idea=text.split('```')[0].strip()[:600])


# ------------------------------------------------------------------ driver

class Evo:
    def __init__(self, a):
        a.out = os.path.abspath(a.out)
        self.a, self.out = a, a.out
        os.makedirs(f'{self.out}/progs', exist_ok=True)
        from pk import eventbench as EB
        self.train = [int(x) for x in a.train.split(',')]
        self.target = {n: EB.frontier_at(n, a.cutoff) for n in self.train}
        self.system = open(os.path.join(HERE, 'evolve_system.md')).read()
        self.db = []                                      # program records
        self.lock = threading.Lock()
        self.rng = random.Random(a.seed)
        self.best_n = {}                                  # n -> (s, pid)
        dbf = f'{self.out}/programs.jsonl'
        if os.path.exists(dbf):
            self.db = [json.loads(l) for l in open(dbf)]
            for r in self.db:
                for n, e in r.get('per_n', {}).items():
                    if 's' in e and (int(n) not in self.best_n or e['s'] < self.best_n[int(n)][0]):
                        self.best_n[int(n)] = (e['s'], r['id'])

    def norm_gap(self, n, s):
        g = math.ceil(math.sqrt(n) - 1e-12)
        return (s - self.target[n]) / max(g - self.target[n], 1e-9)

    def evaluate(self, pool, path):
        jobs = [(path, n, self.a.seed * 1000 + n, self.a.maxc, (1.0, 1.03)) for n in self.train]
        res = list(pool.map(eval_job, jobs))
        per = {}
        for r in res:
            n = r['n']
            if 's' in r:
                per[n] = dict(s=r['s'], g=self.norm_gap(n, r['s']), desc=r['desc'], sq=r['sq'], cands=r['cands'])
            else:
                per[n] = dict(err=r['err'], g=1.5)
        fit = float(np.mean([min(per[n]['g'], 1.5) for n in self.train]))
        return fit, per

    def report_text(self, rec, with_map=True):
        lines = [f"fitness {rec['fit']:.4f}"]
        for n in self.train:
            e = rec['per_n'][str(n)] if str(n) in rec['per_n'] else rec['per_n'][n]
            if 's' in e:
                lines.append(f"n={n}: side {e['s']:.6f}, best known {self.target[n]:.6f}, grid {math.ceil(math.sqrt(n) - 1e-12)}, "
                             f"normalised gap {e['g']:+.3f}; {e['desc']}")
            else:
                lines.append(f"n={n}: ERROR {e['err']}")
        if with_map:
            ok = [(e['g'], n) for n in self.train for e in [rec['per_n'][str(n)]] if 's' in e]
            if ok:
                g, n = min(ok)
                from pk.packing import Packing
                e = rec['per_n'][str(n)]
                lines.append(f"\nMap of the optimised n={n} layout (2 chars per unit, '#' axis-aligned centre, '/' tilted "
                             f"centre, '.' empty; top row = top of the box):\n" + ascii_map(Packing(e['s'], np.array(e['sq']))))
        return '\n'.join(lines)

    def pick(self, k):
        """parent (island k, tournament 3) + inspiration (best of another island)"""
        isl = [r for r in self.db if r.get('island') == k and r['fit'] < 1.4]
        others = [r for r in self.db if r.get('island') != k and r['fit'] < 1.4]
        par = min(self.rng.sample(isl, min(3, len(isl))), key=lambda r: r['fit']) if isl else None
        ins = min(others, key=lambda r: r['fit']) if others else None
        return par, ins

    def prompt_for(self, kind, k):
        targets = ', '.join(f'n={n} (best known {self.target[n]:.6f})' for n in self.train)
        if kind == 'seed':
            return (f'Target n: {targets}.\n\nWrite a first constructor program.  Pick one clear structural idea '
                    f'(idea #{k + 1} of several independent attempts: make it different from the most obvious design).'), []
        if kind == 'director':
            top = sorted(self.db, key=lambda r: r['fit'])[:5]
            body = '\n\n'.join(f"--- program {r['id']} (fitness {r['fit']:.4f}) idea: {r.get('idea', '')[:400]}\n"
                               f"{self.report_text(r, with_map=False)}" for r in top)
            return (f'Target n: {targets}.\n\nThe current best programs and their results are summarised below.  They '
                    f'may all be stuck in one kind of structure.  Propose a SUBSTANTIALLY DIFFERENT structural design (not a '
                    f'parameter tweak) that could beat them, and write it as a complete program.\n\n{body}'), [r['id'] for r in top]
        par, ins = self.pick(k)
        if par is None:
            return self.prompt_for('seed', k)
        txt = (f'Target n: {targets}.\n\nParent program (fitness {par["fit"]:.4f}):\n```python\n{par["code"]}```\n\n'
               f'Its evaluation:\n{self.report_text(par)}\n')
        ids = [par['id']]
        if ins is not None and ins['id'] != par['id']:
            txt += (f'\nAnother good program for inspiration (fitness {ins["fit"]:.4f}; idea: {ins.get("idea", "")[:300]}):\n'
                    f'```python\n{ins["code"]}```\n')
            ids.append(ins['id'])
        txt += ('\nImprove the parent: lower fitness.  You may change parameters, structure, or combine ideas from the other '
                'program; keep what works.  Errors above must be fixed.')
        return txt, ids

    def one(self, pool, i):
        a = self.a
        k = i % a.islands
        if i < a.seeds:
            kind, model = 'seed', 'sonnet'
        elif a.director and i % a.director == 0:
            kind, model = 'director', 'opus'
            k = max(range(a.islands), key=lambda j: min([r['fit'] for r in self.db if r.get('island') == j] or [9]))
        else:
            kind = 'mutate'
            model = 'sonnet' if self.rng.random() < a.mix else 'haiku'
        with self.lock:
            prompt, parents = self.prompt_for(kind, k)
        res = call_llm(model, prompt, self.system)
        rec = dict(id=i, island=k, kind=kind, model=model, parents=parents, usd=res.get('usd'), tokens=res.get('tokens'),
                   llm_sec=round(res.get('sec', 0), 1), idea=res.get('idea', ''), t=time.strftime('%H:%M:%S'))
        code = res.get('code')
        if not code:
            rec.update(code='', fit=1.5, per_n={}, err=res.get('err') or 'no python block')
        else:
            bad = static_check(code)
            path = f'{self.out}/progs/p{i:04d}.py'
            open(path, 'w').write(code)
            if bad:
                rec.update(code=code, fit=1.5, per_n={}, err=bad)
            else:
                t0 = time.time()
                fit, per = self.evaluate(pool, path)
                rec.update(code=code, fit=fit, per_n={str(n): e for n, e in per.items()}, eval_sec=round(time.time() - t0, 1))
        with self.lock:
            self.db.append(rec)
            for n, e in rec['per_n'].items():
                if 's' in e and (int(n) not in self.best_n or e['s'] < self.best_n[int(n)][0]):
                    self.best_n[int(n)] = (e['s'], i)
            with open(f'{self.out}/programs.jsonl', 'a') as f:
                f.write(json.dumps(rec) + '\n')
            best = min(r['fit'] for r in self.db)
            print(f"[{rec['t']}] #{i} {kind:8s} {model:6s} isl {k} fit {rec['fit']:.4f} (best {best:.4f}) "
                  f"llm {rec['llm_sec']}s {rec.get('err', '')[:80]}", flush=True)

    def run(self):
        from pk import manifest
        a = self.a
        man = manifest.start(self.out, 'evolve', args=vars(a), procs=a.procs)
        start = len(self.db)
        with cf.ProcessPoolExecutor(a.procs) as pool, cf.ThreadPoolExecutor(a.llm_par) as tp:
            futs = []
            for i in range(start, a.proposals):
                while len([f for f in futs if not f.done()]) >= a.llm_par:
                    time.sleep(1)
                if i >= a.seeds and len(self.db) < min(a.seeds, 4):
                    while len(self.db) < min(a.seeds, 4):    # mutations need parents
                        time.sleep(2)
                futs.append(tp.submit(self.one, pool, i))
            for f in futs:
                f.result()
        man.finish()
        print('best per n:', {n: (round(s, 6), round(s - self.target[n], 6), pid) for n, (s, pid) in sorted(self.best_n.items())})


def report(a):
    """best --top programs on train + holdout, more seeds; tiers vs the cutoff frontier and the posted J value"""
    from pk import eventbench as EB
    a.out = os.path.abspath(a.out)
    db = [json.loads(l) for l in open(f'{a.out}/programs.jsonl')]
    top = sorted([r for r in db if r['fit'] < 1.4], key=lambda r: r['fit'])[:a.top]
    E = json.load(open(os.path.join(PK, 'bench', 'events.json'))); E = E['events'] if isinstance(E, dict) else E
    J = {e['n']: e['new'] for e in E if e['family'] == 'J'}
    ns = [int(x) for x in (a.train + ',' + a.holdout).split(',')]
    best = {}
    with cf.ProcessPoolExecutor(a.procs) as pool:
        jobs = [(f"{a.out}/progs/p{r['id']:04d}.py", n, sd, a.maxc, (1.0, 1.03, 1.06)) for r in top for n in ns
                for sd in range(a.seeds_report)]
        for (path, n, sd, *_), res in zip(jobs, pool.map(eval_job, jobs)):
            if 's' in res and (n not in best or res['s'] < best[n][0]):
                best[n] = (res['s'], path, res['sq'])
    os.makedirs(f'{a.out}/report', exist_ok=True)
    from collections import Counter
    C = Counter()
    for n in ns:
        T = EB.frontier_at(n, a.cutoff)
        if n not in best:
            print(f'n={n:3d} none'); continue
        s, path, sq = best[n]
        with open(f'{a.out}/report/n{n}.txt', 'w') as f:
            f.write(f'{n} {s!r}\n' + ''.join(f'{x!r} {y!r} {t!r}\n' for x, y, t in sq))
        tier = 'REPLICATED' if n in J and s <= J[n] + 1e-9 else 'BEAT-ANTE' if s < T - 1e-9 else 'no'
        C[tier] += 1
        print(f"n={n:3d} {'train' if str(n) in a.train.split(',') else 'hold '} ours {s:.8f} ante {T:.8f} xu {J.get(n, float('nan')):.8f} "
              f"gap {s - T:+.2e}  frac of Xu's gain {((T - s) / (T - J[n])) if n in J else float('nan'):+.2f} {tier}  ({os.path.basename(path)})")
    print(dict(C))
    print('LLM use:', sum(1 for r in db), 'proposals,', sum(r.get('tokens') or 0 for r in db), 'tokens, list-price equiv $',
          round(sum(r.get('usd') or 0 for r in db), 2), {m: sum(1 for r in db if r['model'] == m) for m in ('haiku', 'sonnet', 'opus')})


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    for name in ('run', 'report'):
        p = sp.add_parser(name); p.add_argument('--out', required=True)
        p.add_argument('--cutoff', default=CUTOFF); p.add_argument('--train', default=TRAIN); p.add_argument('--holdout', default=HOLDOUT)
        p.add_argument('--procs', type=int, default=16); p.add_argument('--maxc', type=int, default=24)
        p.add_argument('--seed', type=int, default=1)
    r = sp.choices['run']
    r.add_argument('--proposals', type=int, default=300); r.add_argument('--llm-par', type=int, default=6)
    r.add_argument('--islands', type=int, default=4); r.add_argument('--seeds', type=int, default=8)
    r.add_argument('--mix', type=float, default=0.3, help='share of sonnet among mutations (rest haiku)')
    r.add_argument('--director', type=int, default=30, help='every k-th proposal is an opus redesign (0 = never)')
    q = sp.choices['report']; q.add_argument('--top', type=int, default=5); q.add_argument('--seeds-report', type=int, default=4)
    a = ap.parse_args()
    if a.cmd == 'run':
        Evo(a).run()
    else:
        report(a)
