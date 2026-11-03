"""Interactive experiments as logged moves (11-02): `pk.py try` and `pk.py stages`.

  pk.py try CHAIN [CHAIN ...] (--ns 791,853 [--best K] [--window W] | --ids 1,2 | --where SQL) --tries T
            [--control CHAIN | --no-control] [--prio high] [--name NAME] [--seeds S] [--hold] [--admit better|window]
  pk.py try --report NAME
  pk.py stages [KIND] [--status lab] [--md]

A try applies every arm (each CHAIN, plus the control: the standard post-construction relax) to the same parents with
the same seeds: T paired applications per (parent, arm), all of them, in shuffled order.  Every application is one pool
unit (pk/tasks.py Executor, campaign try-NAME) and one move-log row (expansion: mechanism = the chain, run = try:NAME,
params = the stage params as drawn + seed, arm, rep, admit); an output below its parent (--admit window: within the
window of the store best at its n) goes to the store with the parent link
(tags try:NAME, and react:hold with --hold), losses are rows with child NULL.  The report reads only those rows, so it
works while the run goes on, after it, and after the driver died.  Same application step as the event loop
(events.apply_chain, Store.record).
"""
from __future__ import annotations
import json, math, os, random, sqlite3, sys, time

CONTROL = 'relax.round.15c.p5'          # DESIGN.md "Stages": the default relax after constructions (CLAIMS C50)
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEAS = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9)"


def store_path():
    from .store import DEFAULT
    return DEFAULT


def ro(path=None):
    return sqlite3.connect(f'file:{path or store_path()}?mode=ro', uri=True, timeout=60)


# ------------------------------------------------------------------ parents

def parents(ns=None, best=4, window=2e-3, exclude=('G', 'D', 'S'), ids=None, where=None):
    """store ids: explicit ids, a WHERE clause on packing (feasible rows only), or per n the `best` best distinct
    packings (side to 6 decimals) within window * sqrt(n) of the first, regular families excluded (as an event-loop pool)"""
    from . import families as F
    db = ro()
    try:
        if ids:
            q = ','.join('?' * len(ids))
            got = [i for (i,) in db.execute(f'SELECT id FROM packing WHERE id IN ({q}) AND {FEAS}', list(ids))]
            if len(got) < len(ids):
                raise SystemExit(f'not feasible store packings: {sorted(set(ids) - set(got))}')
            return list(ids)
        if where:
            return [i for (i,) in db.execute(f'SELECT id FROM packing WHERE ({where}) AND {FEAS} ORDER BY n, s')]
        out = []
        for n in ns:
            out += [i for i, *_ in F.near(db, n, F.kgrid(n) - 1e-9, exclude, best, dedupe=6, window=window * math.sqrt(n),
                                          path=store_path())]
        return out
    finally:
        db.close()


# ------------------------------------------------------------------ one application (a pool unit)

def apply(name, arm, chain, pid, rep, seed, fq_timeout=300.0, tags=(), lim=None):
    """apply chain to store packing pid once -> dict(arm, pid, rep, s0, s, cpu, child); logs the move-log row itself.
    lim None (admit better): the output is stored only below the parent; else below lim (admit window)"""
    import numpy as np
    from .events import apply_chain
    from .store import Store
    st = Store()
    try:
        n, s0, c = st.db.execute('SELECT n, s, coords FROM packing WHERE id=?', (pid,)).fetchone()
        p, s1, params, cpu = apply_chain(chain, n, s0, np.frombuffer(c, '<f8').reshape(-1, 3), seed, fq_timeout)
        from .families import kgrid
        keep = s1 is not None and s1 < min(s0 - 1e-12 if lim is None else lim, kgrid(n) - 1e-9)
        cid = st.record(pid, chain, dict(params, seed=seed, arm=arm, rep=rep, admit='better' if lim is None else 'window'),
                        p if keep else None, cpu, run=f'try:{name}',
                        n=n, s=s1, tags=tags)
        return dict(arm=arm, pid=pid, rep=rep, s0=s0, s=s1, cpu=cpu, child=cid)
    finally:
        st.close()


# ------------------------------------------------------------------ the driver

def run(chains, pids, tries, control=CONTROL, name=None, prio='high', seeds=0, hold=False, fq_timeout=300.0,
        every=30.0, out=sys.stdout, admit='better', window=2e-3):
    """submit every (parent, rep, arm) application as a pool unit and wait, printing the report every `every` s.
    admit 'window': outputs within window * sqrt(n) of the store best at n are stored too (score-neutral steps)"""
    from . import stages as S, tasks
    name = name or time.strftime('t%m%d-%H%M%S')
    arms = ([('ctl', control)] if control else []) + [(f'a{i + 1}', c) for i, c in enumerate(chains)]
    for _, c in arms:
        S.parse(c)                                               # unknown stages / parameters fail here, not in units
    db = ro()
    try:
        if db.execute('SELECT 1 FROM expansion WHERE run=? LIMIT 1', (f'try:{name}',)).fetchone():
            raise SystemExit(f'try {name}: rows exist already (pick another --name; read it with --report {name})')
        lim = {}
        if admit == 'window':
            for pid, n in db.execute(f"SELECT id, n FROM packing WHERE id IN ({','.join('?' * len(pids))})", pids).fetchall():
                b = db.execute(f'SELECT MIN(s) FROM packing WHERE n=? AND {FEAS}', (n,)).fetchone()[0]
                lim[pid] = b + window * math.sqrt(n)
    finally:
        db.close()
    rng = random.Random(f'{name}:{seeds}')
    plan = [(pid, j, random.Random(f'{seeds}:{pid}:{j}').randrange(1 << 30)) for pid in pids for j in range(tries)]
    units = [(arm, c, pid, j, sd) for pid, j, sd in plan for arm, c in arms]
    rng.shuffle(units)                                           # arms interleave in time (equal load on average)
    os.makedirs(os.path.join(HERE, 'runs', 'try'), exist_ok=True)
    json.dump(dict(name=name, arms=dict(arms), parents=pids, tries=tries, seeds=seeds, hold=hold, prio=prio, admit=admit,
                   lim=lim,
                   store=store_path(), argv=sys.argv, t=time.time()),
              open(os.path.join(HERE, 'runs', 'try', f'{name}.json'), 'w'), indent=1)
    tags = (f'try:{name}',) + (('react:hold',) if hold else ())
    print(f'try {name}: {len(arms)} arms x {len(pids)} parents x {tries} tries = {len(units)} units '
          f'(pool campaign try-{name}, prio {prio}, store {store_path()})', file=out, flush=True)
    t0, last, done, errs = time.time(), time.time(), 0, []
    with tasks.Executor(f'try-{name}', prio=prio, timeout=fq_timeout + 300) as ex:
        fs = [ex.submit(apply, name, arm, c, pid, j, sd, fq_timeout, tags, lim.get(pid)) for arm, c, pid, j, sd in units]
        while done < len(fs):
            time.sleep(1.0)
            done = sum(f.done() for f in fs)
            if time.time() - last > every:
                last = time.time()
                print(f'\n[{done}/{len(fs)} units, {(last - t0) / 60:.1f} min]\n' + report(name), file=out, flush=True)
        errs = [repr(f.exception())[:200] for f in fs if f.exception() is not None]
    print(f'\n[done: {len(fs)} units, {(time.time() - t0) / 60:.1f} min, {len(errs)} unit errors]\n' + report(name),
          file=out, flush=True)
    for e in errs[:5]:
        print('  error:', e, file=out)
    return name


# ------------------------------------------------------------------ report

def rows(name, path=None):
    """move-log rows of a try -> [dict(pid, n, chain, arm, rep, cpu, s0, s, child)]"""
    db = ro(path)
    try:
        R = db.execute('SELECT e.packing, e.n, e.mechanism, e.params, e.cpu_s, e.best_s, e.child, p.s FROM expansion e '
                       'JOIN packing p ON p.id=e.packing WHERE e.run=?', (f'try:{name}',)).fetchall()
    finally:
        db.close()
    out = []
    for pid, n, mech, params, cpu, s, child, s0 in R:
        P = json.loads(params or '{}')
        out.append(dict(pid=pid, n=n, chain=mech, arm=P.get('arm', '?'), rep=P.get('rep'), cpu=cpu, s0=s0, s=s, child=child))
    return out


def report(name, path=None):
    """per arm x n: tries, failures, wins vs parent, best gain, gain / CPU-h; paired vs the control (same parent, seed):
    sign counts and median of d = r(ctl) - r(arm), r = min(parent side, output side) (d > 0: the arm did better)"""
    R = rows(name, path)
    if not R:
        return f'try {name}: no rows yet'
    med = lambda xs: sorted(xs)[len(xs) // 2] if xs else float('nan')
    res = lambda r: min(r['s0'], r['s']) if r['s'] is not None else r['s0']
    chain = {r['arm']: r['chain'] for r in R}
    ctl = {(r['pid'], r['rep']): res(r) for r in R if r['arm'] == 'ctl'}
    out = [f'# try {name}: {len(R)} applications, {sum(r["cpu"] for r in R) / 3600:.2f} CPU-h', '',
           '| arm | n | tries | failed | wins | best gain | gain / CPU-h | pairs | + / - / = vs ctl | median paired diff |',
           '|---|---|---|---|---|---|---|---|---|---|']
    keys = sorted({(r['arm'], r['n']) for r in R}, key=lambda k: (k[0] != 'ctl', k))
    for arm, n in keys + [(a, None) for a in sorted({a for a, _ in keys}, key=lambda a: (a != 'ctl', a))]:
        G = [r for r in R if r['arm'] == arm and (n is None or r['n'] == n)]
        g = [max(0.0, r['s0'] - r['s']) for r in G if r['s'] is not None]
        h = sum(r['cpu'] for r in G) / 3600
        D = [ctl[(r['pid'], r['rep'])] - res(r) for r in G if arm != 'ctl' and (r['pid'], r['rep']) in ctl]
        pm = f'{sum(d > 1e-12 for d in D)} / {sum(d < -1e-12 for d in D)} / {sum(abs(d) <= 1e-12 for d in D)}' if D else ''
        out.append(f"| {arm} | {n if n is not None else 'all'} | {len(G)} | {sum(r['s'] is None for r in G)} | "
                   f"{sum(x > 1e-12 for x in g)} | {max(g, default=0):.3e} | {sum(g) / max(h, 1e-9):.3e} | {len(D) or ''} | "
                   f"{pm} | {f'{med(D):+.3e}' if D else ''} |")
    out += [''] + [f'  {a}: {c}' for a, c in sorted(chain.items(), key=lambda x: (x[0] != 'ctl', x[0]))]
    return '\n'.join(out)


# ------------------------------------------------------------------ stage registry x move log

def stage_rows(status=None, kind=None, path=None):
    """every registered stage (name, kind, status, params, evidence, doc) with its move-log use: applications, children,
    CPU-h and last use over expansion rows whose chain (mechanism) or logged params name the stage"""
    from . import stages as S
    use = {}
    try:
        db = ro(path)
        try:
            R = db.execute("SELECT mechanism, (SELECT group_concat(key, '|') FROM json_each(CASE WHEN json_valid(params) "
                           "THEN params ELSE '{}' END) WHERE type='object'), COUNT(*), SUM(child IS NOT NULL), SUM(cpu_s), "
                           "MAX(added) FROM expansion GROUP BY 1, 2").fetchall()
        finally:
            db.close()
    except sqlite3.Error:
        R = []
    for mech, keys, k, w, c, t in R:
        for nm in S.names(mech) | S.names(keys):
            a = use.setdefault(nm, [0, 0, 0.0, ''])
            a[0] += k; a[1] += w or 0; a[2] += c or 0.0; a[3] = max(a[3], t or '')
    out = []
    for nm, f in sorted(S.STAGES.items()):
        st_ = getattr(f, 'status', 'default')
        if (status and st_ != status) or (kind and f.kind != kind):
            continue
        a = use.get(nm, [0, 0, 0.0, ''])
        out.append(dict(name=nm, kind=f.kind, status=st_, evidence=getattr(f, 'evidence', ''),
                        params=', '.join(f'{k}={v}' for k, v in f.space.items()),
                        doc=((f.__doc__ or '').strip().splitlines() or [''])[0], apps=a[0], children=a[1],
                        cpu_h=a[2] / 3600, last=a[3][:10]))
    return out


def stage_table(rows_, md=False):
    if md:
        cell = lambda x: str(x).replace('|', '\\|')
        out = ['| stage | kind | status | params | evidence | what | apps | children | CPU-h | last used |',
               '|---|---|---|---|---|---|---|---|---|---|']
        out += [f"| `{r['name']}` | {r['kind']} | {r['status']} | {cell(r['params'])} | {cell(r['evidence'])} | "
                f"{cell(r['doc'])} | {r['apps']} | {r['children']} | {r['cpu_h']:.1f} | {r['last']} |" for r in rows_]
        return '\n'.join(out)
    out = []
    for r in rows_:
        out.append(f"{r['kind']:9s} {r['name']:22s} {r['status']:7s} {r['apps']:7d} apps {r['children']:6d} ch "
                   f"{r['cpu_h']:7.1f} CPU-h {r['last']:10s} {r['doc'][:70]}")
        if r['params'] or r['evidence']:
            out.append(f"{'':40s}{r['params']}" + (f"  [{r['evidence']}]" if r['evidence'] else ''))
    return '\n'.join(out)
