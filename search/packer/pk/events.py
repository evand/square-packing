"""Event-driven search loop (10-27; replaces the goal loop, the unit loop's arm bandit and pk/queue.py).

A new packing is an event: it joins a pool of parents, and persistent workers claim (parent, move) options one at a
time, ranked by a scorer.  A move is a stage chain applied to the parent until the unit's CPU allowance is spent; every
application is logged in the store's move log (`expansion`, child = the stored output); its outputs within the window
of the best at their n become parents at once.  No dispatcher, no quantum, no drain.

  python3 -m pk.events add CAMPAIGN          register a campaign (a YAML file in a plugin's goals/, by name or path);
                                             again after editing it: new targets get pools, budget / unit edits apply
  python3 -m pk.events work --procs N [--max-idle S]
  python3 -m pk.events status | report CAMPAIGN | stop [CAMPAIGN] | resume
Launch through the pool: `pk/tasks.py repeat NAME --budget CPU_H --prio low -- python3 -m pk.events work --campaign NAME --max-units 3`.

Campaign file (the pre-registration: `predictions` is required; the reading of the result goes under `result`):
  goal: NAME
  budget_cpu_h: 20
  predictions: [...]
  events:
    ns: [...]                  # or n: [lo, hi]; targets.  Parents may sit at n - dn for cross-n moves
    weight: {default: 1.0}     # w(n): utility multiplier per target n (production: 1 = score units per CPU)
    moves:                     # name -> stage chain; dn: target n - parent n; once: one application per parent
      relax:  {chain: 'relax.round.05c.p5', once: true, relax: true}
      add:    {chain: 'grow.add | relax.round.05c.p5', dn: 1}
      lcut:   {chain: 'shrink.lcut | relax.round.05c.p5', dn: '-(2*k-3)'}   # dn in k = ceil sqrt(parent n): an outer L
      build:  {chain: 'gen.random | relax.round.15c.p5', source: true}     # a source: no parent (below)
    policies: [rule, uniform]  # conditions run side by side, equal CPU, each with its own pool and bests (A/B);
                               # 'rule+uniform' = one condition mixing scorers per claim (the scorer is in the unit's info)
    frozen: []                 # conditions whose outputs never join their pool (isolates the value of chaining)
    window: 2.0e-3             # parents within this score gap (side / sqrt n) of the condition's best at their n
    pool: 8                    # start pool per (condition, parent n): the best distinct store packings in the window
    dedupe: 6                  # distinct = side rounded to this many decimals (4: near-copies of one basin merge)
    exclude: [G, D, S]         # parent families never used (regular packings: perturbing them does not work)
    unit_cpu: 30               # CPU-s per unit (a once-move runs once)
    batch: 64                  # claims are queued in batches (10-29): one scoring pass fills up to `batch` entries,
    draws: 8                   #   round-robin over `draws` Thompson draws by rank, one entry per parent and at most
    batch_per_n: 8             #   batch_per_n per parent n; refilled at <= batch / 4 left or after batch_age s; a
    batch_age: 300             #   new condition best at n drops that n's entries and triggers a full refill (at most
    batch_gap: 10              #   one per batch_gap s; batch_age is a backstop).  Scorers see ctx['draws'] and may
                               #   return a (options, draws) array (one feature pass for every draw)
    stop: {conds: [a, b], p: 0.95, min_cpu_h: 2}   # optional sequential stop rule (below)

Sources (moves with `source: true`, chain starting with a gen stage): built once per (target n, generator version,
variant) (stage hooks `version(ctx)` (ctx = {'n': n, 'side': known side}; None = skip n) / `variants`, as for deterministic generators), refreshed by a unit of their own when
the hook `stale()` says so; claimed before any scored option.  Their outputs join every condition's pool and count as
nobody's gain, so conditions stay comparable.

Stop rule (pre-registered with the campaign): conditions a and b are compared on improvements of their own best per
CPU-h (units with a score gain; Gamma(0.5 + k, CPU-h) posteriors); once both have min_cpu_h, the campaign stops as soon
as P(rate a > rate b) leaves [1 - p, p].  The decision is kept in ctl (`stop:<campaign>`) and shown by `report`.

Rules, every claim:
  - targets: grid-only n are eligible only as the first grid-only entry of their row (k = ceil sqrt n, c = k^2 - n), and
    only once the row above has that c (`eligible`; Evan 10-27: monotonicity makes deeper entries no easier).  Known
    sides = the store best and the register at campaign start, then each condition's own bests.
  - a parent that is not relaxed offers only the campaign's relax move (marked `relax: true`), so every other move is
    credited against a relaxed parent; outputs of a chain ending in a relax stage enter relaxed.
  - relax options go first (nearest the best first); otherwise options are queued by sampled utility x w(n) (batch
    above); scorers (SCORERS, by the condition name before '-'):
    uniform, rule (Evan's: near the best, few tries, not regular) and plugin models.  A scorer gets the option rows and
    returns sampled utilities (Thompson draws or a weighted race), so exploration lives in the scorer; it may instead
    return (utilities, settings): settings[i] = {stage: {param: value}} for option i, chosen from `open_space(chain)`
    (the planner as optimizer), used by every application of the unit; else the chain draws per application.  ctx carries the
    campaign config, the condition, the time and {move: mean CPU-s per application so far}.
"""
from __future__ import annotations
import argparse, json, math, os, random, resource, socket, sqlite3, sys, time, traceback
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = os.path.join(HERE, 'runs', 'events.sqlite')
RUNS = os.path.join(HERE, 'runs', 'events')
FEAS = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9)"
SCHEMA = """
CREATE TABLE IF NOT EXISTS campaign(name TEXT PRIMARY KEY, path TEXT, budget_cpu_h REAL, unit_cpu REAL, status TEXT,
  t0 REAL, added REAL);
CREATE TABLE IF NOT EXISTS parent(campaign TEXT, cond TEXT, pid INT, n INT, s REAL, fam TEXT, t_in REAL, fresh INT,
  relaxed INT, unit INT, PRIMARY KEY(campaign, cond, pid));
CREATE INDEX IF NOT EXISTS parent_n ON parent(campaign, cond, n, s);
CREATE TABLE IF NOT EXISTS ptry(campaign TEXT, cond TEXT, pid INT, move TEXT, tries INT, wins INT,
  PRIMARY KEY(campaign, cond, pid, move));
CREATE TABLE IF NOT EXISTS nbest(campaign TEXT, cond TEXT, n INT, s0 REAL, s REAL, t REAL, PRIMARY KEY(campaign, cond, n));
CREATE TABLE IF NOT EXISTS known(campaign TEXT, n INT, s REAL, PRIMARY KEY(campaign, n));
CREATE TABLE IF NOT EXISTS unit(id INTEGER PRIMARY KEY, campaign TEXT, cond TEXT, pid INT, n0 INT, n INT, move TEXT,
  worker INT, seed INT, status TEXT, lease_until REAL, started REAL, ended REAL, cpu REAL, wall REAL, apps INT DEFAULT 0,
  wins INT DEFAULT 0, gain REAL DEFAULT 0, best_s REAL, claim_s REAL, record_s REAL, util REAL, options INT, info TEXT);
CREATE INDEX IF NOT EXISTS unit_c ON unit(campaign, status);
CREATE INDEX IF NOT EXISTS unit_n ON unit(campaign, cond, n, ended);
CREATE TABLE IF NOT EXISTS worker(id INTEGER PRIMARY KEY, pid INT, host TEXT, job TEXT, started REAL, last_seen REAL,
  cpu REAL DEFAULT 0, unit_cpu REAL DEFAULT 0, idle_s REAL DEFAULT 0, claim_s REAL DEFAULT 0, record_s REAL DEFAULT 0,
  units INT DEFAULT 0, status TEXT);
CREATE TABLE IF NOT EXISTS ctl(key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS queue(id INTEGER PRIMARY KEY, campaign TEXT, cond TEXT, pid INT, n0 INT, n INT, move TEXT,
  util REAL, params TEXT, scorer TEXT, nopt INT, t REAL, UNIQUE(campaign, cond, pid, move));
CREATE INDEX IF NOT EXISTS queue_c ON queue(campaign, cond, id);
CREATE TABLE IF NOT EXISTS gen_done(campaign TEXT, move TEXT, n INT, version TEXT, variant INT, unit INT,
  PRIMARY KEY(campaign, move, n, version, variant));
"""
DEFAULT = dict(weight={'default': 1.0}, moves={}, policies=['rule', 'uniform'], frozen=[], window=2.0e-3, pool=8, dedupe=6,
               exclude=['G', 'D', 'S'], unit_cpu=30.0, recent_h=2.0, fq_timeout=300, stop=None, parents=None,
               batch=64, draws=8, batch_per_n=8, batch_age=300.0, batch_gap=10.0)
STOP_EVERY = 60.0                               # s between checks of a campaign's stop rule
LEASE = 10.0                                   # x unit_cpu (min 600 s): past it, a running unit's worker is lost
STALE = 300.0                                  # s: a worker not seen for this long counts as down
FAIL_RUN = 3                                   # a unit stops after this many failed applications in a row
SOURCE = '_source'                             # the cond of source / refresh units (shared by every condition)
REFRESH_GAP = 30.0                             # s between stale() checks of a source's inputs
FILL_LEASE = 120.0                             # s: a claimed queue refill not finished by then may be taken over
FILL_WAIT = 30.0                               # s a claim waits on an empty queue that another worker is filling


from .families import kgrid  # noqa: E402


def cpu():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def connect(path=None):
    path = path or DB
    os.makedirs(os.path.dirname(path), exist_ok=True)
    db = sqlite3.connect(path, timeout=120, isolation_level=None)   # autocommit; transactions are explicit
    db.execute('PRAGMA journal_mode=WAL'); db.execute('PRAGMA synchronous=NORMAL')
    db.executescript(SCHEMA)
    if 't' not in {r[1] for r in db.execute('PRAGMA table_info(nbest)')}:      # 10-27: time of the last improvement
        db.execute('ALTER TABLE nbest ADD COLUMN t REAL')
    return db


def store_path():
    from .store import DEFAULT as P
    return P


def store_ro():
    return sqlite3.connect(f'file:{store_path()}?mode=ro', uri=True, timeout=60)


# ------------------------------------------------------------------ targets: the grid-only rule

def eligible(n, side):
    """may target n be worked?  side: {n: best known side} (missing = the grid side).  Sub-grid n: yes.  Grid-only n
    (best known = k): only the first grid-only entry of its row (c = c(k) - 1, c(k) = the smallest c with a known
    sub-k packing; 2k - 2 if none), and only once row k - 1 has a sub-(k-1) packing at the same c."""
    k = kgrid(n); c = k * k - n
    if side.get(n, k) < k - 1e-9:
        return True
    sub = [cc for cc in range(0, 2 * k - 1) if side.get(k * k - cc, k) < k - 1e-9]
    first = (min(sub) - 1) if sub else 2 * k - 2
    if c != first:
        return False
    k1 = k - 1
    if c > 2 * k1 - 2:                         # the row above has no entry at this c
        return True
    return side.get(k1 * k1 - c, k1) < k1 - 1e-9


def row_ns(ns):
    """every n in the triangle rows of ns and the rows above (what `eligible` reads)"""
    out = set()
    for n in ns:
        k = kgrid(n)
        out.update(range(max(2, (k - 2) ** 2 + 1), k * k + 1))
    return sorted(out)


# ------------------------------------------------------------------ campaigns

def goal_dirs():
    from . import stages as S
    return [os.path.join(os.path.dirname(os.path.abspath(m.__file__)), 'goals') for m in S.PLUGINS] + \
        [os.path.join(HERE, 'goals')]


def resolve(spec):
    """a campaign file path, or a name looked up in each plugin's goals/ (then this repo's goals/)"""
    if os.path.exists(spec):
        return os.path.abspath(spec)
    name = os.path.basename(spec)
    name = name if name.endswith(('.yaml', '.yml')) else name + '.yaml'
    for d in goal_dirs():
        if os.path.exists(os.path.join(d, name)):
            return os.path.join(d, name)
    raise FileNotFoundError(f'campaign {spec}: not a file, and not in {goal_dirs()}')


_CFG = {}


def load(spec):
    """campaign file -> dict(name, budget_cpu_h, path, **events config); cached per (path, mtime)"""
    import yaml
    path = resolve(spec)
    key = (path, os.path.getmtime(path))
    if key not in _CFG:
        d = yaml.safe_load(open(path))
        if not d.get('predictions'):
            raise ValueError(f'{path}: no predictions (pre-register before running)')
        if 'events' not in d:
            raise ValueError(f'{path}: no events block')
        e = dict(DEFAULT, **d['events'])
        unknown = set(d['events']) - set(DEFAULT) - {'ns', 'n'}
        if e.get('stop') and len(e['stop'].get('conds', [])) != 2:
            raise ValueError(f'{path}: stop.conds names two conditions')
        if unknown:
            raise ValueError(f'{path}: unknown events keys {sorted(unknown)}')
        ns = sorted(int(x) for x in e['ns']) if e.get('ns') else list(range(int(e['n'][0]), int(e['n'][1]) + 1))
        moves = {}
        for name, m in e['moves'].items():
            m = dict(m) if isinstance(m, dict) else dict(chain=m)
            dn = m.get('dn', 0); dn = int(dn) if isinstance(dn, (int, float)) else str(dn)
            src = bool(m.get('source', False))
            if src and not m['chain'].strip().startswith('gen.'):
                raise ValueError(f'{path}: source move {name} must start with a gen stage')
            moves[name] = dict(chain=m['chain'], dn=dn, dnf=_dn_fn(dn), once=bool(m.get('once', False)), source=src,
                               relax=bool(m.get('relax', False)), tags=m.get('tags', {}),
                               relaxed_out=m['chain'].split('|')[-1].strip().startswith('relax.'))
        w = {k if k == 'default' else int(k): float(v) for k, v in (e['weight'] or {}).items()}
        for k in ('window', 'unit_cpu', 'recent_h', 'fq_timeout'):
            e[k] = float(e[k])
        _CFG[key] = dict(e, name=d['goal'], budget_cpu_h=float(d['budget_cpu_h']), path=path, ns=ns, moves=moves, weight=w,
                         pool=int(e['pool']), dedupe=int(e['dedupe']))
    return _CFG[key]


def _dn_fn(spec):
    """dn: an int, or an arithmetic expression in k = ceil sqrt(parent n), e.g. '-(2*k-3)' (an outer L of a parent with
    side in (k - 1, k))"""
    if isinstance(spec, int):
        return lambda n0: spec
    import ast
    tree = ast.parse(str(spec), mode='eval')
    ok = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Load, ast.Add, ast.Sub, ast.Mult, ast.USub)
    if not all(isinstance(x, ok) for x in ast.walk(tree)) or {x.id for x in ast.walk(tree) if isinstance(x, ast.Name)} - {'k'}:
        raise ValueError(f'dn {spec!r}: integers, k, + - * only')
    code = compile(tree, '<dn>', 'eval')
    return lambda n0: int(eval(code, {'__builtins__': {}}, {'k': kgrid(n0)}))


def target(cfg, move, n0):
    return n0 + cfg['moves'][move]['dnf'](n0)


def weight(cfg, n):
    return cfg['weight'].get(n, cfg['weight'].get('default', 1.0))


def parent_ns(cfg):
    """every parent n some move maps into a target"""
    T = set(cfg['ns']); out = set()
    for name, m in cfg['moves'].items():
        if m['source']:
            continue
        if isinstance(m['dn'], int):
            out.update(n - m['dn'] for n in T)
            continue
        lo, hi = min(T), max(T)
        for n0 in range(2, 2 * hi + 50):
            if target(cfg, name, n0) in T:
                out.add(n0)
    return sorted(n for n in out if n >= 2)


def _relaxed(st, ids):
    """store ids that are outputs of a logged move ending in a relax, or had a relax move logged on them"""
    if not ids:
        return set()
    q = ','.join('?' * len(ids))
    a = {i for (i,) in st.execute(f"SELECT DISTINCT packing FROM expansion WHERE mechanism LIKE 'relax%' AND packing IN ({q})", ids)}
    b = {i for (i,) in st.execute(f'SELECT DISTINCT child FROM expansion WHERE child IN ({q})', ids)}
    return a | b


def add(db, spec):
    """register a campaign: start pools per condition, known sides, condition bests"""
    from . import families as F
    cfg = load(spec); name = cfg['name']
    old = db.execute('SELECT 1 FROM campaign WHERE name=?', (name,)).fetchone() is not None
    have = {n for (n,) in db.execute('SELECT DISTINCT n FROM nbest WHERE campaign=?', (name,))}
    import datetime
    ts = lambda x: datetime.datetime.fromisoformat(x[:19]).timestamp()
    st = store_ro()
    try:
        known, since = {}, {}
        for n in row_ns(set(cfg['ns']) | set(parent_ns(cfg))):
            r = st.execute(f'SELECT s, added FROM packing WHERE n=? AND {FEAS} ORDER BY s LIMIT 1', (n,)).fetchone()
            f = st.execute('SELECT MIN(s) FROM frontier WHERE n=?', (n,)).fetchone()[0]
            v = min(x for x in (r[0] if r else None, f, float(kgrid(n))) if x is not None)
            known[n] = v
            if r:
                since[n] = r[1]                      # when the store best at n arrived (its last improvement)
        rows = []
        if cfg['parents']:                           # explicit start parents (store ids) instead of the nearest pools
            ids = [int(i) for i in cfg['parents']]
            pn = set(parent_ns(cfg))
            for i, n, s, fam in st.execute(f"SELECT id, n, s, family FROM packing WHERE id IN ({','.join('?' * len(ids))}) "
                                           f"AND {FEAS}", ids):
                if n in pn and n not in have:
                    rows.append((n, i, s, fam))
        for n in ([] if cfg['parents'] else parent_ns(cfg)):
            if n < 2 or n in have:                  # an edited campaign file: only the new parent n get pools
                continue
            for i, s, _, fam in F.near(st, n, kgrid(n) - 1e-9, cfg['exclude'], cfg['pool'], dedupe=cfg['dedupe'],
                                       window=cfg['window'] * math.sqrt(n)):
                rows.append((n, i, s, fam))
        rel = _relaxed(st, [i for _, i, _, _ in rows])
        added = dict(st.execute(f"SELECT id, added FROM packing WHERE id IN ({','.join('?' * len(rows)) or 'NULL'})",
                                [i for _, i, _, _ in rows]))
    finally:
        st.close()
    t0 = time.time()
    db.execute('BEGIN IMMEDIATE')
    if old:                                          # re-add: re-activate, pick up budget / unit edits and new targets
        db.execute("UPDATE campaign SET status='active', path=?, budget_cpu_h=?, unit_cpu=? WHERE name=?",
                   (cfg['path'], cfg['budget_cpu_h'], cfg['unit_cpu'], name))
    else:
        db.execute('INSERT INTO campaign VALUES (?,?,?,?,?,?,?)', (name, cfg['path'], cfg['budget_cpu_h'], cfg['unit_cpu'], 'active', t0, t0))
    db.executemany('INSERT OR IGNORE INTO known VALUES (?,?,?)', [(name, n, s) for n, s in known.items()])
    ns = sorted(set(cfg['ns']) | set(parent_ns(cfg)))
    for cond in cfg['policies']:
        db.executemany('INSERT OR IGNORE INTO parent VALUES (?,?,?,?,?,?,?,0,?,NULL)',
                       [(name, cond, i, n, s, fam, ts(added[i]), int(i in rel)) for n, i, s, fam in rows])
        db.executemany('INSERT OR IGNORE INTO nbest VALUES (?,?,?,?,?,?)',
                       [(name, cond, n, known[n], known[n], ts(since[n]) if n in since else None) for n in ns if n in known])
    db.execute('COMMIT')
    return cfg, len(rows)


# ------------------------------------------------------------------ scorers

SCORERS = {}


def scorer(name):
    def deco(f):
        SCORERS[name] = f
        return f
    return deco


@scorer('uniform')
def _uniform(opts, ctx, rng):
    return np.array([rng.random() for _ in opts])


@scorer('rule')
def _rule(opts, ctx, rng):
    """Evan's rule (10-27): near the best, few tries, not regular (regular families are never parents); proportional
    sampling as an exponential race (argmax of w / Exp(1))"""
    w = np.array([1.0 / (1.0 + o['gap'] / 1e-5) / (1.0 + o['tries'] / 2.0) for o in opts])
    return w / np.array([rng.expovariate(1.0) for _ in opts])


def policy(cond, rng=None):
    """a condition's scorer: the name before '-'; 'a+b+...' mixes scorers, one drawn per claim (Evan 10-27: when a
    choice does not matter for the question at hand, mix the options and keep learning)"""
    p = cond.split('-')[0]
    if '+' in p:
        return (rng or random).choice(p.split('+'))
    return p


_OPEN = {}


def open_space(chain):
    """{stage: {param: space}} for the chain's parameters left to draw: not given, or given as a space (a range);
    what a scorer may choose for an option (10-27, planner as optimizer)"""
    if chain not in _OPEN:
        from . import stages as S
        out = {}
        for name, kw in S.parse(chain):
            sp = {k: (kw[k] if S.is_space(kw.get(k)) else v) for k, v in S.STAGES[name].space.items()
                  if (k not in kw or S.is_space(kw[k])) and isinstance(kw.get(k, v), tuple)}
            if sp:
                out[name] = sp
        _OPEN[chain] = out
    return _OPEN[chain]


# ------------------------------------------------------------------ claim

def expire(db, now):
    lost = db.execute("SELECT id, campaign, cond, pid, move FROM unit WHERE status='running' AND lease_until<?", (now,)).fetchall()
    for i, c, cond, pid, mv in lost:                       # a lost unit gives its try back (else a lost relax blocks
        db.execute("UPDATE unit SET status='lost' WHERE id=?", (i,))                  # its parent for good)
        db.execute('UPDATE ptry SET tries=tries-1 WHERE campaign=? AND cond=? AND pid=? AND move=?', (c, cond, pid, mv))
    lost = [r[0] for r in lost]
    db.execute("UPDATE worker SET status='down' WHERE status='up' AND last_seen<?", (now - STALE,))
    return lost


def spent(db, name):
    """CPU-s per condition: finished units at their CPU, running ones at the allowance"""
    uc = db.execute('SELECT unit_cpu FROM campaign WHERE name=?', (name,)).fetchone()[0]
    out = {}
    for cond, c in db.execute("SELECT cond, SUM(CASE WHEN status='running' THEN ? ELSE COALESCE(cpu, 0) END) FROM unit "
                              "WHERE campaign=? AND status!='lost' GROUP BY cond", (uc, name)):
        out[cond] = c or 0.0
    return out


def sides(db, name, cond):
    """known sides for the eligibility rule: campaign-start known sides, improved by the condition's own bests"""
    side = dict(db.execute('SELECT n, s FROM known WHERE campaign=?', (name,)))
    for n, s in db.execute('SELECT n, s FROM nbest WHERE campaign=? AND cond=?', (name, cond)):
        side[n] = min(side.get(n, s), s)
    return side


def options(db, cfg, cond, now):
    """every (parent, move) option of a condition, with the features scorers read"""
    name = cfg['name']
    best, tb = {}, {}
    for n, s_, t_ in db.execute('SELECT n, s, t FROM nbest WHERE campaign=? AND cond=?', (name, cond)):
        best[n], tb[n] = s_, t_
    side = sides(db, name, cond)
    ok = {n: eligible(n, side) for n in cfg['ns']}
    tries = {}
    for pid, mv, t, w in db.execute('SELECT pid, move, tries, wins FROM ptry WHERE campaign=? AND cond=?', (name, cond)):
        tries[(pid, mv)] = (t, w)
    tot = {}
    for (pid, _), (t, w) in tries.items():
        a, b = tot.get(pid, (0, 0)); tot[pid] = (a + t, b + w)
    rec = dict(db.execute("SELECT n, SUM(apps) + SUM(status='running') FROM unit WHERE campaign=? AND cond=? AND "
                          "(status='running' OR ended>?) GROUP BY n", (name, cond, now - 3600 * cfg['recent_h'])))
    relax = [m for m, v in cfg['moves'].items() if v['relax']]
    out = []
    for pid, n0, s, fam, t_in, fresh, relaxed in db.execute(
            'SELECT pid, n, s, fam, t_in, fresh, relaxed FROM parent WHERE campaign=? AND cond=?', (name, cond)):
        b = best.get(n0, s)
        gap = max(0.0, (s - b) / math.sqrt(n0))
        if gap > cfg['window']:
            continue
        mvs = relax if (relax and not relaxed) else [m for m, v in cfg['moves'].items() if m not in relax and not v['source']]
        for m in mvs:
            v = cfg['moves'][m]; n = n0 + v['dnf'](n0)
            if n not in ok or not ok[n]:
                continue
            t, w = tries.get((pid, m), (0, 0))
            if v['once'] and t:
                continue
            T, W = tot.get(pid, (0, 0))
            out.append(dict(pid=pid, n0=n0, n=n, s=s, fam=fam, move=m, tags=v['tags'], gap=gap, best=b,
                            known=side.get(n0, b), age_h=max(now - t_in, 1.0) / 3600, fresh=fresh, tries=T, wins=W,
                            n_age_h=(now - tb[n0]) / 3600 if tb.get(n0) else 1000.0,
                            mtries=t, mwins=w, recent=rec.get(n, 0) or 0))
    return out


_STOP_AT = {}


def stop_state(db, cfg):
    """the stop rule's current reading -> dict(p, k, h, decided) or None (no rule / not enough CPU yet)"""
    r = cfg.get('stop')
    if not r:
        return None
    a, b = r['conds']
    k, h = {}, {}
    for c in (a, b):
        k[c], h[c] = db.execute('SELECT COUNT(CASE WHEN gain > 0 THEN 1 END), COALESCE(SUM(cpu), 0) / 3600 FROM unit '
                                'WHERE campaign=? AND cond=? AND cpu IS NOT NULL', (cfg['name'], c)).fetchone()
    if min(h.values()) < float(r.get('min_cpu_h', 2.0)):
        return None
    g = np.random.default_rng(0)
    ra = g.gamma(0.5 + k[a], 1.0 / h[a], 20000); rb = g.gamma(0.5 + k[b], 1.0 / h[b], 20000)
    p = float(np.mean(ra > rb)); q = float(r.get('p', 0.95))
    return dict(p=p, k=k, h={c: round(v, 3) for c, v in h.items()}, decided=p > q or p < 1 - q)


_STALE_AT = {}


def _source_option(db, cfg, now):
    """the next source unit: a refresh when a source's inputs are stale (hook `stale()`, checked every REFRESH_GAP s,
    one refresh at a time), else a build of the first (move, n, variant) not built at the generator's current version
    for that n (hook `version(ctx)`, ctx = {'n': n, 'side': the campaign's known side at n or None}: a generator whose
    output at n depends on n-specific inputs returns a per-n key, so an input change rebuilds only the n it affects, and
    None when it has nothing worth building at n; it is called for every target n on every claim, so it must be cheap;
    `variants` per version; static generators: once) at an eligible target n -> option or None"""
    from . import stages as S
    name = cfg['name']
    srcs = [(m, v) for m, v in cfg['moves'].items() if v['source']]
    if not srcs:
        return None
    running = {(m, i) for m, i in db.execute("SELECT move, info FROM unit WHERE campaign=? AND cond=? AND status='running'",
                                              (name, SOURCE))}
    side = dict(db.execute('SELECT n, s FROM known WHERE campaign=?', (name,)))
    for m, v in srcs:
        f = S.STAGES[S.parse(v['chain'].split('|')[0])[0][0]]
        if hasattr(f, 'stale') and now - _STALE_AT.get((name, m), 0.0) > REFRESH_GAP:
            _STALE_AT[(name, m)] = now
            if not any(mm == m and '"refresh"' in (i or '') for mm, i in running) and f.stale():
                return dict(kind='refresh', move=m, n=None)
        done = {(n, v, k) for n, v, k in db.execute('SELECT n, version, variant FROM gen_done WHERE campaign=? AND move=?',
                                                    (name, m))}
        for n in cfg['ns']:
            if not eligible(n, side):
                continue
            ver = f.version(dict(n=n, side=side.get(n))) if hasattr(f, 'version') else 'static'
            if ver is None:                               # the generator has nothing to build at n
                continue
            ver = str(ver)
            for k in range(int(getattr(f, 'variants', 1))):
                if (n, ver, k) not in done:
                    return dict(kind='build', move=m, n=n, version=ver, variant=k)
    return None


def _open_campaigns(db, wid, rng, now, t0, only=None):
    """write transaction 1 (short): expire lost units, apply stop rules and budgets, claim a due source unit ->
    the source unit dict, or [(name, cfg, spent)] of the campaigns with scored options to look at"""
    db.execute('BEGIN IMMEDIATE')
    try:
        expire(db, now)
        live = []
        for (name, uc, budget) in db.execute("SELECT name, unit_cpu, budget_cpu_h FROM campaign WHERE status='active' "
                                             "AND (?1 IS NULL OR name=?1) ORDER BY added", (only,)).fetchall():
            cfg = load(db.execute('SELECT path FROM campaign WHERE name=?', (name,)).fetchone()[0])
            if cfg.get('stop') and now - _STOP_AT.get(name, 0.0) > STOP_EVERY:
                _STOP_AT[name] = now
                z = stop_state(db, cfg)
                if z and z['decided']:
                    db.execute("UPDATE campaign SET status='done' WHERE name=?", (name,))
                    db.execute('INSERT OR REPLACE INTO ctl VALUES (?,?)', (f'stop:{name}', json.dumps(dict(z, t=now))))
                    continue
            sp = spent(db, name)
            if sum(sp.values()) + uc > 3600 * budget:
                if not db.execute("SELECT 1 FROM unit WHERE campaign=? AND status='running'", (name,)).fetchone():
                    db.execute("UPDATE campaign SET status='done' WHERE name=?", (name,))
                continue
            o = _source_option(db, cfg, now)
            if o is not None:
                uid = db.execute('INSERT INTO unit(campaign, cond, pid, n0, n, move, worker, seed, status, lease_until, started, '
                                 'options, claim_s, info) VALUES (?,?,NULL,NULL,?,?,?,?,?,?,?,0,?,?)',
                                 (name, SOURCE, o['n'], o['move'], wid, rng.randrange(1 << 30), 'running',
                                  now + max(600.0, LEASE * uc), now, time.time() - t0, json.dumps(o))).lastrowid
                if o['kind'] == 'build':
                    db.execute('INSERT OR IGNORE INTO gen_done VALUES (?,?,?,?,?,?)', (name, o['move'], o['n'], o['version'], o['variant'], uid))
                db.execute('COMMIT')
                return dict(id=uid, campaign=name, cond=SOURCE, seed=uid, cfg=cfg, pid=None, n0=None, s=None, **o)
            live.append((name, cfg, sp))
        db.execute('COMMIT')
        return live
    except BaseException:
        db.execute('ROLLBACK')
        raise


def _jsonable(x):
    return x.item() if hasattr(x, 'item') else str(x)


def _draw_columns(opts, cfg, cond, name, db, now, rng):
    """scored utilities for the options: [(utilities, chosen settings or None, scorer name)] per Thompson draw.
    The scorer gets ctx['draws'] = D and may return a (len(opts), D) array (and settings as D lists): one feature pass
    for every draw; a scorer returning one column is called D times"""
    D = max(1, int(cfg['draws']))
    per = {m: c / a for m, c, a in db.execute('SELECT move, SUM(cpu), SUM(apps) FROM unit WHERE campaign=? '
                                              'AND apps>0 GROUP BY move', (name,))}
    w = np.array([weight(cfg, o['n']) for o in opts])
    cols = []
    while len(cols) < D:
        sc = policy(cond, rng)
        r = SCORERS[sc](opts, dict(cfg=cfg, cond=cond, now=now, cpu=per, draws=D - len(cols)), rng)
        u, ch = r if isinstance(r, tuple) else (r, None)
        u = np.asarray(u, float)
        if u.ndim == 1:
            cols.append((u * w, ch, sc))
        else:
            for d in range(u.shape[1]):
                cols.append((u[:, d] * w, ch[d] if ch is not None else None, sc))
    return cols[:D]


def _fill(db, name, cfg, cond, rng, now):
    """refill the condition's claim queue (Evan 10-29: near-best options are plentiful, so pick in batches): options and
    scores in a read transaction (no write lock), then the queue replaced in one short write.  Relax options first
    (nearest the best first, as before); else round-robin over D Thompson draws by rank, at most one entry per parent
    and batch_per_n per parent n, up to `batch` entries -> entries queued"""
    B = max(1, int(cfg['batch'])); cap = max(1, int(cfg['batch_per_n']))
    t_fill = time.time()                               # a new best after this moment keeps its refill flag
    db.execute('BEGIN')
    try:
        opts = options(db, cfg, cond, now)
        rows = []
        rel = [j for j, o in enumerate(opts) if cfg['moves'][o['move']]['relax']]
        if rel:
            for j in sorted(rel, key=lambda j: opts[j]['gap'])[:B]:
                rows.append((j, math.inf, None, None))
        elif opts:
            cols = _draw_columns(opts, cfg, cond, name, db, now, rng)
            order = [np.argsort(-u, kind='stable') for u, _, _ in cols]
            used, per_n = set(), {}
            for r in range(len(opts)):
                for d, (u, ch, sc) in enumerate(cols):
                    if len(rows) >= B:
                        break
                    j = int(order[d][r]); o = opts[j]
                    if o['pid'] in used or per_n.get(o['n0'], 0) >= cap or not u[j] > 0:
                        continue
                    used.add(o['pid']); per_n[o['n0']] = per_n.get(o['n0'], 0) + 1
                    rows.append((j, float(u[j]), ch[j] if ch is not None else None, sc))
                if len(rows) >= B:
                    break
            if not rows:                               # every utility <= 0: fall back to the first draw's argmax
                u, ch, sc = cols[0]; j = int(np.argmax(u))
                rows.append((j, float(u[j]), ch[j] if ch is not None else None, sc))
    finally:
        db.execute('COMMIT')
    db.execute('BEGIN IMMEDIATE')
    try:
        db.execute('DELETE FROM queue WHERE campaign=? AND cond=?', (name, cond))
        for j, ui, prm, sc in rows:
            o = opts[j]
            db.execute('INSERT OR IGNORE INTO queue(campaign, cond, pid, n0, n, move, util, params, scorer, nopt, t) '
                       'VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                       (name, cond, o['pid'], o['n0'], o['n'], o['move'], min(ui, 1e300),
                        None if prm is None else json.dumps(prm, default=_jsonable), sc, len(opts), now))
        db.execute('DELETE FROM ctl WHERE key=?', (f'fill:{name}:{cond}',))
        db.execute('INSERT OR REPLACE INTO ctl VALUES (?,?)', (f'filled:{name}:{cond}', str(t_fill)))
        db.execute('DELETE FROM ctl WHERE key=? AND CAST(value AS REAL)<=?', (f'due:{name}:{cond}', t_fill))
        db.execute('COMMIT')
    except BaseException:
        db.execute('ROLLBACK')
        raise
    return len(rows)


def _pop(db, wid, rng, name, cfg, cond, now, t0):
    """one write transaction: take the queue head of (campaign, cond) as a unit -> (unit or None, refill due).  A
    once-move tried since it was queued is dropped (another worker took it); refill is due when the queue is short
    (<= batch / 4), older than batch_age (a backstop), or a condition best moved since the last fill (record() sets the
    flag; at most one such refill per batch_gap s), and the caller then holds the refill lease"""
    B = max(1, int(cfg['batch']))
    db.execute('BEGIN IMMEDIATE')
    try:
        uc, budget, status = db.execute('SELECT unit_cpu, budget_cpu_h, status FROM campaign WHERE name=?', (name,)).fetchone()
        if status != 'active' or sum(spent(db, name).values()) + uc > 3600 * budget:
            db.execute('COMMIT')
            return None, False
        u = None
        while True:
            q = db.execute('SELECT id, pid, n0, n, move, util, params, scorer, nopt FROM queue WHERE campaign=? AND cond=? '
                           'ORDER BY id LIMIT 1', (name, cond)).fetchone()
            if q is None:
                break
            qid, pid, n0, n, mv, ui, prm, sc, nopt = q
            db.execute('DELETE FROM queue WHERE id=?', (qid,))
            if mv not in cfg['moves']:
                continue
            t = db.execute('SELECT tries FROM ptry WHERE campaign=? AND cond=? AND pid=? AND move=?', (name, cond, pid, mv)).fetchone()
            if cfg['moves'][mv]['once'] and t and t[0]:
                continue
            s = db.execute('SELECT s FROM parent WHERE campaign=? AND cond=? AND pid=?', (name, cond, pid)).fetchone()
            if s is None:
                continue
            seed = rng.randrange(1 << 30)
            uid = db.execute('INSERT INTO unit(campaign, cond, pid, n0, n, move, worker, seed, status, lease_until, started, '
                             'util, options, claim_s) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                             (name, cond, pid, n0, n, mv, wid, seed, 'running', now + max(600.0, LEASE * uc), now,
                              ui, nopt, time.time() - t0)).lastrowid
            db.execute('INSERT INTO ptry VALUES (?,?,?,?,1,0) ON CONFLICT(campaign, cond, pid, move) DO UPDATE SET '
                       'tries=tries+1', (name, cond, pid, mv))
            u = dict(id=uid, campaign=name, cond=cond, seed=seed, cfg=cfg, params=json.loads(prm) if prm else None,
                     scorer=sc, pid=pid, n0=n0, n=n, move=mv, s=s[0])
            break
        left, oldest = db.execute('SELECT COUNT(*), MIN(t) FROM queue WHERE campaign=? AND cond=?', (name, cond)).fetchone()
        due = u is None or left <= B // 4 or (oldest is not None and now - oldest > float(cfg['batch_age']))
        if not due:                                   # a new condition best since the last fill (coalesced: at most
            ev = db.execute('SELECT value FROM ctl WHERE key=?', (f'due:{name}:{cond}',)).fetchone()   # one per batch_gap s)
            if ev is not None:
                f = db.execute('SELECT value FROM ctl WHERE key=?', (f'filled:{name}:{cond}',)).fetchone()
                tf = float(f[0]) if f else 0.0
                due = float(ev[0]) > tf and time.time() - tf >= float(cfg['batch_gap'])
        if due and u is not None:                     # one refill at a time while the queue still serves claims
            key = f'fill:{name}:{cond}'
            r = db.execute('SELECT value FROM ctl WHERE key=?', (key,)).fetchone()
            if r is not None and time.time() - float(r[0]) < FILL_LEASE:
                due = False
            else:
                db.execute('INSERT OR REPLACE INTO ctl VALUES (?,?)', (key, str(time.time())))
        db.execute('COMMIT')
        return u, due
    except BaseException:
        db.execute('ROLLBACK')
        raise


def _take_lease(db, name, cond):
    """the refill lease of (campaign, cond) in ctl -> True if this caller now holds it"""
    key = f'fill:{name}:{cond}'; t = time.time()
    db.execute('BEGIN IMMEDIATE')
    try:
        r = db.execute('SELECT value FROM ctl WHERE key=?', (key,)).fetchone()
        ok = r is None or t - float(r[0]) >= FILL_LEASE
        if ok:
            db.execute('INSERT OR REPLACE INTO ctl VALUES (?,?)', (key, str(t)))
        db.execute('COMMIT')
        return ok
    except BaseException:
        db.execute('ROLLBACK')
        raise


def claim(db, wid, rng, now=None, only=None):
    """the next unit for worker wid -> unit dict or None.  A short write transaction for expiry, stop rules, budgets
    and sources; then the head of a condition's claim queue (conditions in spent-CPU order).  An empty queue is filled
    first by the claiming worker; a short or old one is refilled by the worker that noticed, after taking its unit
    (one refill at a time per condition: lease in ctl)."""
    now = time.time() if now is None else now
    t0 = time.time()
    live = _open_campaigns(db, wid, rng, now, t0, only)
    if isinstance(live, dict):
        return live
    for name, cfg, sp in live:
        for cond in sorted(cfg['policies'], key=lambda c: (sp.get(c, 0.0), rng.random())):
            u, due = _pop(db, wid, rng, name, cfg, cond, now, t0)
            if u is None:
                if not due:
                    continue                            # campaign closed or over budget
                w0 = time.time()                        # empty queue: one worker fills, the others wait for it
                while not _take_lease(db, name, cond) and time.time() - w0 < FILL_WAIT:
                    time.sleep(0.25)
                    u, due = _pop(db, wid, rng, name, cfg, cond, now, t0)
                    if u is not None or not due:
                        break
                if u is None and due:
                    if not _fill(db, name, cfg, cond, rng, now):
                        continue                        # no options for this condition
                    u, due = _pop(db, wid, rng, name, cfg, cond, now, t0)
                if u is None:
                    continue
            if due:
                _fill(db, name, cfg, cond, rng, now)
            return u
    return None


def drop_queued(db, name, cond, n):
    """queued entries whose parent or target is at n (its condition best moved: the gap features changed)"""
    db.execute('DELETE FROM queue WHERE campaign=? AND cond=? AND (n0=? OR n=?)', (name, cond, n, n))


# ------------------------------------------------------------------ run (in the worker)

def run_unit(u):
    """apply the move's chain to the parent until the allowance is spent (once-moves: one application).  Every
    application -> one move-log row; outputs below the parent (same n) or within the window of the condition's best
    (cross n), not of an excluded family -> stored with the parent link -> returned as new parents."""
    from . import stages as S, families as F
    from .packing import Packing
    from .store import Store
    cfg = u['cfg']; mv = cfg['moves'][u['move']]; n = u['n']
    if u['cond'] == SOURCE:
        return run_source(u)
    st0 = store_ro()
    try:
        r = st0.execute('SELECT s, coords FROM packing WHERE id=?', (u['pid'],)).fetchone()
    finally:
        st0.close()
    if r is None:
        return dict(status='empty', apps=[], children=[])
    s0, sq0 = float(r[0]), np.frombuffer(r[1], '<f8').reshape(-1, 3).copy()
    db = connect(); bref = db.execute('SELECT s FROM nbest WHERE campaign=? AND cond=? AND n=?',
                                      (u['campaign'], u['cond'], n)).fetchone(); db.close()
    bn = bref[0] if bref else float(kgrid(n))
    lim = (s0 - 1e-12) if mv['dn'] == 0 else bn + cfg['window'] * math.sqrt(n)
    rng = random.Random(u['seed']); c0 = cpu()
    run = f"runs/events/{u['campaign']}/u{u['id']:07d}"
    spec = [(name, dict(kw, **(u.get('params') or {}).get(name, {}))) for name, kw in S.parse(mv['chain'])]
    apps, children = [], []
    store = None
    try:
        while True:
            seed = rng.randrange(1 << 30); t = cpu()
            try:
                st, lg = S.run_chain(spec, dict(n=n, fq_timeout=cfg['fq_timeout']), seed, st=S.State(Packing(s0, sq0.copy())))
                p = st.p
                ok = p.n == n and not S.failed(lg) and float(p.max_pen()) <= 1e-9
                params = {r_['stage']: {k: v for k, v in r_.get('params', {}).items() if isinstance(v, (int, float, str, bool))}
                          for r_ in lg if r_.get('stage') != 'keep.best'}
            except Exception as e:
                ok, p, params = False, None, dict(error=repr(e)[:160])
            c = cpu() - t
            s1 = float(p.s) if ok else None
            cid = None
            if ok and s1 < min(lim, kgrid(n) - 1e-9):
                store = store or Store()
                cid, created = store.add(p, f'run:{run}', parent=u['pid'], status='raw')
                fam = F.fill([(cid, s1, p.sq.astype('<f8').tobytes(), None)], n, path=store_path()).get(cid)
                if fam not in cfg['exclude']:
                    children.append((cid, n, s1, fam))
            store = store or Store()
            store.record_expansion(u['pid'], u['n0'], u['move'], c, params=dict(params, seed=seed), run=run, nprop=1,
                                   best_s=s1, child=cid)
            # a win is a new best of the condition at the target n, for every move (Evan 10-27: within and across n alike;
            # below the parent is not a win: a drop child always is, and the scorer models new bests)
            apps.append(dict(s=s1, cpu=round(c, 2), win=int(s1 is not None and s1 < bn - 1e-12)))
            if mv['once'] or cpu() - c0 > cfg['unit_cpu'] - (c if apps else 0) * 0.5:
                break
            if len(apps) >= FAIL_RUN and all(a['s'] is None for a in apps[-FAIL_RUN:]):
                break                                    # this setting keeps failing: give the claim back
    finally:
        if store is not None:
            store.close()
    return dict(status='done', apps=apps, children=children,
                info={k: u[k] for k in ('params', 'scorer') if u.get(k)})


def run_source(u):
    """a refresh (the generator's `refresh()` hook) or one build of a source move at target n; outputs within the window
    of the best known side join every condition's pool (shared, so conditions stay comparable)"""
    from . import stages as S, families as F
    from .store import Store
    cfg = u['cfg']; mv = cfg['moves'][u['move']]
    f = S.STAGES[S.parse(mv['chain'].split('|')[0])[0][0]]
    if u['kind'] == 'refresh':
        return dict(status='done', apps=[], children=[], info=dict(refresh=f.refresh() or {}))
    n = u['n']
    db = connect(); r = db.execute('SELECT MIN(s) FROM nbest WHERE campaign=? AND n=?', (u['campaign'], n)).fetchone(); db.close()
    lim = (r[0] if r and r[0] is not None else float(kgrid(n))) + cfg['window'] * math.sqrt(n)
    t = cpu(); children = []
    try:
        st, lg = S.run_chain(mv['chain'], dict(n=n, fq_timeout=cfg['fq_timeout']), u['seed'])
        p = st.p; ok = p.n == n and not S.failed(lg) and float(p.max_pen()) <= 1e-9
    except Exception as e:
        ok, p, lg = False, None, [dict(error=repr(e)[:160])]
    s1 = float(p.s) if ok else None
    if ok and s1 < min(lim, kgrid(n) - 1e-9):
        store = Store()
        try:
            cid, _ = store.add(p, f"run:runs/events/{u['campaign']}/u{u['id']:07d}", status='raw')
        finally:
            store.close()
        fam = F.fill([(cid, s1, p.sq.astype('<f8').tobytes(), None)], n, path=store_path()).get(cid)
        if fam not in cfg['exclude']:
            children.append((cid, n, s1, fam))
    return dict(status='done', apps=[dict(s=s1, cpu=round(cpu() - t, 2), win=0)], children=children,
                info=dict(version=u.get('version'), variant=u.get('variant')))


# ------------------------------------------------------------------ record

def record(db, u, res, cpu_s, wall, now=None):
    """one transaction: unit row, tries / wins, new parents, the condition's bests -> score gain"""
    now = time.time() if now is None else now
    name, cond, cfg = u['campaign'], u['cond'], u['cfg']
    mv = cfg['moves'][u['move']]
    apps = res.get('apps', [])
    wins = sum(a['win'] for a in apps)
    db.execute('BEGIN IMMEDIATE')
    try:
        if apps and u['pid'] is not None:
            db.execute('UPDATE ptry SET tries=tries+?, wins=wins+? WHERE campaign=? AND cond=? AND pid=? AND move=?',
                       (len(apps) - 1, wins, name, cond, u['pid'], u['move']))
        if mv['relax']:
            db.execute('UPDATE parent SET relaxed=1 WHERE campaign=? AND pid=?', (name, u['pid']))
        gain = 0.0
        conds = cfg['policies'] if cond == SOURCE else [cond]          # source outputs: every condition, gain to none
        for c in conds:
            for cid, n, s, fam in res.get('children', []):
                if c not in cfg['frozen']:
                    db.execute('INSERT OR IGNORE INTO parent VALUES (?,?,?,?,?,?,?,1,?,?)',
                               (name, c, cid, n, s, fam, now, int(mv['relaxed_out']), u['id']))
                r = db.execute('SELECT s FROM nbest WHERE campaign=? AND cond=? AND n=?', (name, c, n)).fetchone()
                if r is None:
                    db.execute('INSERT INTO nbest VALUES (?,?,?,?,?,?)', (name, c, n, s, s, now))
                elif s < r[0]:
                    gain += (r[0] - s) / math.sqrt(n) if c == cond else 0.0
                    db.execute('UPDATE nbest SET s=?, t=? WHERE campaign=? AND cond=? AND n=?', (s, now, name, c, n))
                    drop_queued(db, name, c, n)
                    db.execute('INSERT OR REPLACE INTO ctl VALUES (?,?)', (f'due:{name}:{c}', str(time.time())))
        best = min([a['s'] for a in apps if a['s'] is not None], default=None)
        db.execute('UPDATE unit SET status=?, ended=?, cpu=?, wall=?, apps=?, wins=?, gain=?, best_s=?, info=? WHERE id=?',
                   (res.get('status', 'done'), now, cpu_s, wall, len(apps), wins, gain, best,
                    json.dumps(res.get('info', {}))[:2000], u['id']))
        db.execute('COMMIT')
    except BaseException:
        db.execute('ROLLBACK')
        raise
    return gain


# ------------------------------------------------------------------ workers

def worker(i, path, max_idle, job, max_units=None, only=None):
    """claim / run / record until the ctl stop flag, `max_idle` seconds with nothing to do, or `max_units` units (for
    the pool, pk/tasks.py: nothing to claim then ends at once) -> units run.  only: claim from this campaign alone (a
    pool campaign's priority is then its event campaign's)"""
    db = connect(path)
    rng = random.Random(os.getpid() * 1000003 + int(time.time() * 1000))
    wid = db.execute("INSERT INTO worker(pid, host, job, started, last_seen, status) VALUES (?,?,?,?,?,'up')",
                     (os.getpid(), socket.gethostname(), job, time.time(), time.time())).lastrowid
    idle0 = None; done = 0
    try:
        while (max_units is None or done < max_units) and (db.execute("SELECT value FROM ctl WHERE key='stop'").fetchone() or ['0'])[0] != '1':
            w = time.time()
            u = claim(db, wid, rng, only=only)
            cs = time.time() - w
            if u is None:
                if max_units is not None:
                    break
                idle0 = idle0 or time.time()
                if max_idle is not None and time.time() - idle0 > max_idle:
                    break
                time.sleep(10)
                db.execute("UPDATE worker SET idle_s=idle_s+?, claim_s=claim_s+?, last_seen=?, cpu=?, status='up' WHERE id=?",
                           (time.time() - w - cs, cs, time.time(), cpu(), wid))
                continue
            idle0 = None; done += 1
            c0, w0 = cpu(), time.time()
            try:
                res = run_unit(u)
            except Exception:
                res = dict(status='error', info=dict(error=traceback.format_exc()[-800:]))
            c, wl = cpu() - c0, time.time() - w0
            w1 = time.time()
            record(db, u, res, c, wl)
            rs = time.time() - w1
            db.execute('UPDATE unit SET record_s=? WHERE id=?', (rs, u['id']))
            db.execute("UPDATE worker SET units=units+1, unit_cpu=unit_cpu+?, claim_s=claim_s+?, record_s=record_s+?, "
                       "last_seen=?, cpu=?, status='up' WHERE id=?", (c, cs, rs, time.time(), cpu(), wid))
    finally:
        db.execute("UPDATE worker SET status='down', last_seen=?, cpu=? WHERE id=?", (time.time(), cpu(), wid))
    return done


# ------------------------------------------------------------------ reports

def report(db, name):
    r = db.execute('SELECT path, budget_cpu_h, unit_cpu, t0, status FROM campaign WHERE name=?', (name,)).fetchone()
    if not r:
        return f'no campaign {name}'
    med = lambda xs: sorted(xs)[len(xs) // 2] if xs else float('nan')
    out = [f'# campaign {name} ({r[4]}): budget {r[1]} CPU-h, units {r[2]:.0f} CPU-s, start '
           f'{time.strftime("%Y-%m-%d %H:%M", time.localtime(r[3]))}', '',
           '| condition | units | apps | CPU-h | n moved | score gain | score / CPU-h | parents (new) |',
           '|---|---|---|---|---|---|---|---|']
    for cond, in db.execute('SELECT DISTINCT cond FROM nbest WHERE campaign=? ORDER BY cond', (name,)).fetchall():
        U = db.execute('SELECT COUNT(*), COALESCE(SUM(apps), 0), COALESCE(SUM(cpu), 0) / 3600 FROM unit WHERE campaign=? AND '
                       'cond=? AND cpu IS NOT NULL', (name, cond)).fetchone()
        M = db.execute('SELECT n, s0, s FROM nbest WHERE campaign=? AND cond=? AND s < s0 - 1e-12', (name, cond)).fetchall()
        g = sum((a - b) / math.sqrt(n) for n, a, b in M)
        P = db.execute('SELECT COUNT(*), SUM(fresh) FROM parent WHERE campaign=? AND cond=?', (name, cond)).fetchone()
        out.append(f'| {cond} | {U[0]} | {U[1]} | {U[2]:.2f} | {len(M)} | {g:.3e} | {g / max(U[2], 1e-9):.2e} | {P[0]} ({P[1] or 0}) |')
    out += ['', '| condition | move | units | apps | CPU-h | apps below the best | score gain |', '|---|---|---|---|---|---|---|']
    for row in db.execute('SELECT cond, move, COUNT(*), SUM(apps), SUM(cpu) / 3600, SUM(wins), SUM(gain) FROM unit WHERE '
                          'campaign=? AND cpu IS NOT NULL GROUP BY cond, move ORDER BY cond, move', (name,)):
        out.append('| {} | {} | {} | {} | {:.2f} | {} | {:.2e} |'.format(*row))
    S = db.execute("SELECT COUNT(*), SUM(status='error'), SUM(status='lost'), SUM(status='empty') FROM unit WHERE campaign=?",
                   (name,)).fetchone()
    C = [x for (x,) in db.execute('SELECT COALESCE(claim_s, 0) + COALESCE(record_s, 0) FROM unit WHERE campaign=? AND cpu IS NOT NULL',
                                  (name,))]
    W = db.execute('SELECT SUM(last_seen - started), SUM(unit_cpu), SUM(idle_s), COUNT(*) FROM worker WHERE id IN '
                   '(SELECT DISTINCT worker FROM unit WHERE campaign=?)', (name,)).fetchone()
    out += ['', f'units {S[0]}: errors {S[1] or 0}, lost {S[2] or 0}, empty {S[3] or 0}; claim + record median {1000 * med(C):.1f} ms']
    if W[0]:
        out.append(f'workers {W[3]}: utilization {W[1] / W[0]:.0%} (unit CPU / worker wall, all campaigns), idle {W[2] / 3600:.2f} h')
    cfg = load(r[0]) if os.path.exists(r[0]) else {}
    z = stop_state(db, cfg) if cfg.get('stop') else None
    if z:
        d = db.execute('SELECT value FROM ctl WHERE key=?', (f'stop:{name}',)).fetchone()
        a, b = cfg['stop']['conds']
        out.append(f"stop rule {a} vs {b}: P({a} rate > {b} rate) = {z['p']:.3f}; improvements {z['k']}, CPU-h {z['h']}"
                   + (f"; STOPPED at {time.strftime('%H:%M', time.localtime(json.loads(d[0])['t']))}" if d else ''))
    moved = db.execute('SELECT cond, n, s0, s FROM nbest WHERE campaign=? AND s < s0 - 1e-12 ORDER BY n, cond', (name,)).fetchall()
    out += ['', 'moved:'] + [f'  {c} {n}: {a:.9f} -> {b:.9f} ({(a - b) / math.sqrt(n):.2e})' for c, n, a, b in moved]
    return '\n'.join(out)


def status(db):
    out = []
    for name, st, b in db.execute('SELECT name, status, budget_cpu_h FROM campaign ORDER BY added'):
        c = db.execute("SELECT COALESCE(SUM(cpu), 0) / 3600, SUM(status='running') FROM unit WHERE campaign=?", (name,)).fetchone()
        out.append(f'{name:24s} {st:7s} {c[0]:.2f} / {b} CPU-h, running {c[1] or 0}')
    up = db.execute("SELECT COUNT(*), MAX(last_seen) FROM worker WHERE status='up' AND last_seen>?", (time.time() - STALE,)).fetchone()
    out.append(f'workers up: {up[0]}' + (f' (last seen {time.time() - up[1]:.0f} s ago)' if up[1] else ''))
    return '\n'.join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='pk.events', description=__doc__.split('\n')[0])
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('add'); p.add_argument('campaign')
    p = sp.add_parser('work'); p.add_argument('--procs', type=int, default=1); p.add_argument('--max-idle', type=float)
    p.add_argument('--max-units', type=int, help='for the pool: run this many units, exit 75 if there was nothing to do')
    p.add_argument('--campaign', help='claim from this campaign only')
    sp.add_parser('status')
    p = sp.add_parser('report'); p.add_argument('campaign')
    p = sp.add_parser('stop'); p.add_argument('campaign', nargs='?')
    sp.add_parser('resume', help='clear the worker stop flag')
    a = ap.parse_args(argv)
    from . import stages  # noqa: F401  (loads plugins: their stages and scorers)
    db = connect()
    if a.cmd == 'add':
        cfg, k = add(db, a.campaign)
        print(f"{cfg['name']}: {len(cfg['ns'])} targets, {len(cfg['moves'])} moves, policies {cfg['policies']}, "
              f"{k} new start parents per condition, budget {cfg['budget_cpu_h']} CPU-h")
    elif a.cmd == 'work':
        if a.max_units:                          # a pool unit must not clear a stop that someone set
            sys.exit(0 if worker(0, DB, None, os.environ.get('PK_JOB_NAME'), a.max_units, a.campaign) else 75)
        db.execute("INSERT OR REPLACE INTO ctl VALUES ('stop', '0')")
        job = os.environ.get('PK_JOB_NAME')
        if a.procs == 1:
            worker(0, DB, a.max_idle, job, only=a.campaign)
        else:
            import multiprocessing as mp
            ps = [mp.Process(target=worker, args=(i, DB, a.max_idle, job, None, a.campaign)) for i in range(a.procs)]
            for p in ps:
                p.start()
            for p in ps:
                p.join()
    elif a.cmd == 'status':
        print(status(db))
    elif a.cmd == 'report':
        print(report(db, a.campaign))
    elif a.cmd == 'stop':
        if a.campaign:
            db.execute("UPDATE campaign SET status='stopped' WHERE name=?", (a.campaign,))
        else:
            db.execute("INSERT OR REPLACE INTO ctl VALUES ('stop', '1')")
    elif a.cmd == 'resume':
        db.execute("INSERT OR REPLACE INTO ctl VALUES ('stop', '0')")


if __name__ == '__main__':                     # `python3 -m pk.events`: run the imported module, so plugin scorers
    from pk import events as _self              # (registered into pk.events) are the ones claims read
    _self.main()
