"""Interactive experiments as logged moves (11-02): `pk.py try` and `pk.py stages`.

  pk.py try ARM [ARM ...] (--ns 791,853 [--best K] [--window W] | --ids 1,2 | --where SQL) --tries T
            [--control ARM | --no-control] [--prio high] [--name NAME] [--seeds S] [--hold] [--admit better|window]
            [--requires F,F] [--asof ISO]
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

Arms (11-02): an ARM is a chain, 'dn=D: CHAIN' (an n-changing arm: parent at n0, ctx n = n0 + D), or '@NAME' (a named
move: a MOVES dict {name: arm} in a loaded plugin / lab module).  Every row is scored against a reference: the store
best at its target n when the plan is made (or among rows added before --asof, to replay a past win).  Wins = outputs
below that reference (events.below: relative tolerance WIN_RTOL); 'below parent' is shown for dn = 0 arms.  The paired
diff is on each application's actual output (the chain's final or best feasible packing at its target n, not capped at
the parent), in score units: u = (output side - reference at n0 + dn) / sqrt(n0 + dn), so arms of different dn share
one scale; d = u(ctl) - u(arm) > 0: the arm did better; pairs where either side failed (no feasible output) are left
out of the paired columns and counted in 'failed'.  The control is the default relax at the parent's n unless --control names another
arm (for an underfill arm, the same chain with dn=0 isolates the n change).  Outputs of dn != 0 arms are stored
within --window of the store best at their n (as the event loop does across n).  A control that is deterministic
(stages.deterministic) runs once per parent and pairs with every rep.  --requires keeps parents whose cached features
(pk/features.py) meet every requirement.
"""
from __future__ import annotations
import json, math, os, random, re, sqlite3, sys, time

CONTROL = 'relax.round.15c.p5'          # DESIGN.md "Stages": the default relax after constructions (CLAIMS C50)
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEAS = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9)"


def store_path():
    from .store import DEFAULT
    return DEFAULT


def ro(path=None):
    return sqlite3.connect(f'file:{path or store_path()}?mode=ro', uri=True, timeout=60)


# ------------------------------------------------------------------ parents

def parents(ns=None, best=4, window=2e-3, exclude=('G', 'D', 'S'), ids=None, where=None, shape=0):
    """store ids: explicit ids, a WHERE clause on packing (feasible rows only), or per n the `best` best distinct
    packings (side to 6 decimals; shape > 0: distinct arrangements, families.near) within window * sqrt(n) of the first,
    regular families excluded (as an event-loop pool)"""
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
                                          path=store_path(), shape=shape)]
        return out
    finally:
        db.close()


# ------------------------------------------------------------------ one application (a pool unit)

_ARM = re.compile(r'^\s*dn\s*=\s*([-+]?\d+)\s*:(.*)$', re.S)


def moves():
    """named arms: MOVES dicts of loaded plugins and lab modules"""
    from . import stages as S
    out = {}
    for m in list(S.PLUGINS) + [v for k, v in sorted(sys.modules.items()) if k.startswith('lab_')]:
        out.update(getattr(m, 'MOVES', None) or {})
    return out


def arm_spec(a):
    """'dn=-1: CHAIN' | '@NAME' | CHAIN -> (dn, chain)"""
    a = a.strip()
    if a.startswith('@'):
        M = moves()
        if a[1:] not in M:
            raise SystemExit(f'unknown move {a} (named moves: {sorted(M)})')
        a = M[a[1:]]
    m = _ARM.match(a)
    return (int(m.group(1)), m.group(2).strip()) if m else (0, a)


def best_at(db, n, asof=None):
    q = f'SELECT MIN(s) FROM packing WHERE n=? AND {FEAS}' + (' AND added<?' if asof else '')
    r = db.execute(q, (n, asof) if asof else (n,)).fetchone()[0]
    return float(r) if r is not None else None


def apply(name, arm, chain, pid, rep, seed, fq_timeout=300.0, tags=(), lim=None, dn=0):
    """apply chain to store packing pid once at target n = parent n + dn -> dict(arm, pid, rep, s0, s, cpu, child);
    logs the move-log row itself.  lim None (admit better, dn 0): the output is stored only below the parent; else
    below lim (admit window / cross n)"""
    import numpy as np
    from .events import apply_chain, WIN_RTOL
    from .store import Store
    st = Store()
    try:
        n, s0, c = st.db.execute('SELECT n, s, coords FROM packing WHERE id=?', (pid,)).fetchone()
        p, s1, params, cpu = apply_chain(chain, n + dn, s0, np.frombuffer(c, '<f8').reshape(-1, 3), seed, fq_timeout)
        from .families import kgrid
        keep = s1 is not None and s1 < min(s0 * (1 - WIN_RTOL) if lim is None else lim, kgrid(n + dn) - 1e-9)
        cid = st.record(pid, chain, dict(params, seed=seed, arm=arm, rep=rep, admit='better' if lim is None else 'window',
                                         **({'dn': dn} if dn else {})),
                        p if keep else None, cpu, run=f'try:{name}',
                        n=n, s=s1, tags=tags)
        return dict(arm=arm, pid=pid, rep=rep, s0=s0, s=s1, cpu=cpu, child=cid)
    finally:
        st.close()


# ------------------------------------------------------------------ the driver

def run(chains, pids, tries, control=CONTROL, name=None, prio='high', seeds=0, hold=False, fq_timeout=300.0,
        every=30.0, out=sys.stdout, admit='better', window=2e-3, asof=None, requires=()):
    """submit every (parent, rep, arm) application as a pool unit and wait, printing the report every `every` s.
    admit 'window': outputs within window * sqrt(n) of the store best at n are stored too (score-neutral steps).
    Arms: see the module doc (dn, @NAME); asof: reference bests from rows added before this ISO time"""
    from . import stages as S, tasks
    name = name or time.strftime('t%m%d-%H%M%S')
    arms = ([('ctl', control)] if control else []) + [(f'a{i + 1}', c) for i, c in enumerate(chains)]
    arms = [(a, *arm_spec(c)) for a, c in arms]                  # (arm, dn, chain)
    for _, _, c in arms:
        S.parse(c)                                               # unknown stages / parameters fail here, not in units
    if requires:
        from . import features as FT
        k = len(pids); pids = FT.filter_ids(pids, requires)
        print(f'requires {list(requires)}: {len(pids)} of {k} parents', file=out)
        if not pids:
            raise SystemExit('no parents meet --requires')
    db = ro()
    try:
        if db.execute('SELECT 1 FROM expansion WHERE run=? LIMIT 1', (f'try:{name}',)).fetchone():
            raise SystemExit(f'try {name}: rows exist already (pick another --name; read it with --report {name})')
        n0 = dict(db.execute(f"SELECT id, n FROM packing WHERE id IN ({','.join('?' * len(pids))})", pids).fetchall())
        T = sorted({n0[p] + dn for p in pids for _, dn, _ in arms} | set(n0.values()))
        ref = {n: best_at(db, n, asof) for n in T}               # the scoring reference per target n
        now = {n: best_at(db, n) for n in T}                     # admission: the store as it is
        lim = {}
        for a, dn, _ in arms:
            for p in pids:
                n = n0[p] + dn
                if (dn or admit == 'window') and now[n] is not None:
                    lim[f'{a}:{p}'] = now[n] + window * math.sqrt(n)
    finally:
        db.close()
    ctl_once = bool(control) and S.deterministic(arms[0][2])
    rng = random.Random(f'{name}:{seeds}')
    plan = [(pid, j, random.Random(f'{seeds}:{pid}:{j}').randrange(1 << 30)) for pid in pids for j in range(tries)]
    rng.shuffle(plan)                                            # (parent, rep) blocks in random order, arms shuffled
    units = []                                                   # within a block: pairs complete as the run goes
    for pid, j, sd in plan:
        blk = [(arm, dn, c, pid, j, sd) for arm, dn, c in arms if not (arm == 'ctl' and ctl_once and j > 0)]
        rng.shuffle(blk); units += blk                           # (a deterministic control: once per parent)
    os.makedirs(os.path.join(HERE, 'runs', 'try'), exist_ok=True)
    json.dump(dict(name=name, arms={a: c for a, _, c in arms}, dn={a: dn for a, dn, _ in arms}, parents=pids,
                   tries=tries, seeds=seeds, hold=hold, prio=prio, admit=admit, lim=lim, asof=asof,
                   ref={str(n): v for n, v in ref.items()}, ctl_once=ctl_once, requires=list(requires),
                   store=store_path(), argv=sys.argv, t=time.time()),
              open(os.path.join(HERE, 'runs', 'try', f'{name}.json'), 'w'), indent=1)
    tags = (f'try:{name}',) + (('react:hold',) if hold else ())
    print(f'try {name}: {len(arms)} arms x {len(pids)} parents x {tries} tries = {len(units)} units '
          f'(pool campaign try-{name}, prio {prio}, store {store_path()})', file=out, flush=True)
    t0, last, done, errs = time.time(), time.time(), 0, []
    with tasks.Executor(f'try-{name}', prio=prio, timeout=fq_timeout + 300) as ex:
        fs = [ex.submit(apply, name, arm, c, pid, j, sd, fq_timeout, tags, lim.get(f'{arm}:{pid}'), dn)
              for arm, dn, c, pid, j, sd in units]
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
    """move-log rows of a try -> [dict(pid, n, nt, dn, chain, arm, rep, cpu, s0, s, child, out, params)] (n: the parent's,
    nt: target; s: the side returned (keep.best may return an earlier stage's packing); out: the chain's own final side
    when keep.best logged it (11-04), else s; params: the logged params)"""
    db = ro(path)
    try:
        R = db.execute('SELECT e.packing, e.n, e.mechanism, e.params, e.cpu_s, e.best_s, e.child, p.s FROM expansion e '
                       'JOIN packing p ON p.id=e.packing WHERE e.run=?', (f'try:{name}',)).fetchall()
    finally:
        db.close()
    out = []
    for pid, n, mech, params, cpu, s, child, s0 in R:
        P = json.loads(params or '{}'); dn = int(P.get('dn', 0))
        out.append(dict(pid=pid, n=n, nt=n + dn, dn=dn, chain=mech, arm=P.get('arm', '?'), rep=P.get('rep'), cpu=cpu, s0=s0,
                        s=s, child=child, out=P.get('_final_s', s) if s is not None else None, params=P))
    return out


def plan_of(name):
    f = os.path.join(HERE, 'runs', 'try', f'{name}.json')
    return json.load(open(f)) if os.path.exists(f) else {}


def report(name, path=None):
    """per arm x target n: tries, failures, outputs below the parent (dn 0), wins (below the reference best at the
    target n), best output - reference, gain / CPU-h (sum of gains below the reference); paired vs the control (same
    parent, seed): sign counts and median of d = u(ctl) - u(arm), u = the output's gap to the reference in score units
    (module doc; d > 0: the arm did better; failed sides left out);
    median output gap: (output side - reference at its n) / sqrt n over feasible outputs (informative when the parent is
    the reference itself, where every u is 0).  Paired diffs and output gaps use the chain's own final side (rows: out),
    so a side-keeping stage early in a chain (its output equals the parent's side) does not mask what follows; wins,
    below parent and best - ref use the side returned (what is stored)"""
    from .events import below, WIN_RTOL
    R = rows(name, path)
    if not R:
        return f'try {name}: no rows yet'
    P = plan_of(name)
    ref = {int(k): v for k, v in (P.get('ref') or {}).items()}
    miss = sorted(({r['nt'] for r in R} | {r['n'] for r in R}) - {k for k, v in ref.items() if v is not None})
    if miss:                                                     # an old plan: the store best now (none: the grid side)
        from .families import kgrid
        db = ro(path)
        try:
            ref.update({n: best_at(db, n) or float(kgrid(n)) for n in miss})
        finally:
            db.close()
    med = lambda xs: sorted(xs)[len(xs) // 2] if xs else float('nan')
    gap = lambda s, n: (s - ref[n]) / math.sqrt(n)

    def u(r):                                                  # the output itself (11-02: not capped at the parent;
        return gap(r['out'], r['nt']) if r['out'] is not None else None   # 11-04: nor at a side-keeping stage's output)
    chain = {r['arm']: (r['dn'], r['chain']) for r in R}
    ctl = {(r['pid'], r['rep']): u(r) for r in R if r['arm'] == 'ctl'}
    ctl1 = {r['pid']: u(r) for r in R if r['arm'] == 'ctl'} if P.get('ctl_once') else {}
    tie = 2 * WIN_RTOL
    out = [f'# try {name}: {len(R)} applications, {sum(r["cpu"] for r in R) / 3600:.2f} CPU-h'
           + (f'; references as of {P["asof"]}' if P.get('asof') else '') + (' (control once per parent)' if ctl1 else ''),
           '', '| arm | n | ref | tries | failed | below parent | wins | best - ref | gain / CPU-h | pairs | '
           '+ / - / = vs ctl | median paired diff (score) | median output gap (score) |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    keys = sorted({(r['arm'], r['nt']) for r in R}, key=lambda k: (k[0] != 'ctl', k))
    for arm, n in keys + [(a, None) for a in sorted({a for a, _ in keys}, key=lambda a: (a != 'ctl', a))]:
        G = [r for r in R if r['arm'] == arm and (n is None or r['nt'] == n)]
        ok = [r for r in G if r['s'] is not None]
        g = [max(0.0, ref[r['nt']] - r['s']) for r in ok if below(r['s'], ref[r['nt']])]
        h = sum(r['cpu'] for r in G) / 3600
        D = [] if arm == 'ctl' else [(ctl.get((r['pid'], r['rep'])) if (r['pid'], r['rep']) in ctl else ctl1.get(r['pid']),
                                      u(r)) for r in G]
        D = [c - a for c, a in D if c is not None and a is not None]
        pm = f'{sum(d > tie for d in D)} / {sum(d < -tie for d in D)} / {sum(abs(d) <= tie for d in D)}' if D else ''
        bp = sum(below(r['s'], r['s0']) for r in ok if r['dn'] == 0) if any(r['dn'] == 0 for r in G) else ''
        bs = min((r['s'] - ref[r['nt']] for r in ok), default=None) if n is not None else None
        og = [gap(r['out'], r['nt']) for r in ok]
        out.append(f"| {arm} | {n if n is not None else 'all'} | {f'{ref[n]:.6f}' if n is not None else ''} | {len(G)} | "
                   f"{len(G) - len(ok)} | {bp} | {len(g)} | {f'{bs:+.3e}' if bs is not None else ''} | "
                   f"{sum(g) / max(h, 1e-9):.3e} | {len(D) or ''} | {pm} | {f'{med(D):+.3e}' if D else ''} | "
                   f"{f'{med(og):+.3e}' if og else ''} |")
    out += [''] + [f'  {a}: ' + (f'dn={d}: ' if d else '') + c for a, (d, c) in sorted(chain.items(), key=lambda x: (x[0] != 'ctl', x[0]))]
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
