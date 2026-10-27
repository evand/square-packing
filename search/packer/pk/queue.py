"""Priority-queue claims for the unit loop (10-26).  A goal with a `queue:` block has no fixed targets: every claim picks
a (parent packing, move) pair anywhere in an n range, and a new packing becomes a parent as soon as its unit records it.

Conditions.  `queue.conditions` lists parent policies that run side by side in one goal, each with a private pool of
parents (the goal-start store, then only its own outputs), so they never feed each other and each is scored by its own
best per n:
  queue   sample (parent, move kind) with probability proportional to a priority function (`queue.priority`, a name in
          PRIORITY; plugins may register their own) of the parent's gap to the condition's best at its n (score
          units), tries of that move kind on that parent, and its age
  base    a parent uniformly from the pool within `base_window` (side) of the best at its n, then a move kind uniformly
          (the donor-window picker without the arm bandit)
A condition's policy is the part of its name before '-' (queue-frozen runs the queue policy).  Conditions listed in
`frozen` never add their own outputs to their pool (parents stay the goal-start store): live vs frozen with the same
policy isolates the value of chaining from new entries.  Gaps are measured against the best in the condition's pool
(for a live condition, its own best).  `ns: [...]` lists the n explicitly instead of the range `n: [lo, hi]`.
Claims go to the condition with the least CPU so far (equal CPU per condition).

Moves.  The goal's auto arms are planned over every n in the range; a move at parent n0 is any arm whose donor n
(arm n + offset) is n0, grouped by kind (G.arm_kind); within a kind one arm is drawn uniformly.  Generator-led arms
(`gen.*`, whose input is the whole store) are not moves: they would leak between conditions.

Rows: qpool (parents), qbest (each condition's best per n), qclaim (every claim's features and weights, for fitting a
priority model later).  Goal block:
  queue: {n: [lo, hi], conditions: [queue, base], window: 2.0e-3, base_window: 0.05, pool: 8, priority: default}
"""
from __future__ import annotations
import json, math, time

SCHEMA = """
CREATE TABLE IF NOT EXISTS qpool(goal TEXT, cond TEXT, n INT, packing INT, s REAL, t REAL, unit INT,
  PRIMARY KEY(goal, cond, packing));
CREATE INDEX IF NOT EXISTS qpool_n ON qpool(goal, cond, n, s);
CREATE TABLE IF NOT EXISTS qbest(goal TEXT, cond TEXT, n INT, s0 REAL, s REAL, PRIMARY KEY(goal, cond, n));
CREATE TABLE IF NOT EXISTS qtry(goal TEXT, cond TEXT, packing INT, kind TEXT, tries INT, PRIMARY KEY(goal, cond, packing, kind));
CREATE TABLE IF NOT EXISTS qclaim(unit INT PRIMARY KEY, goal TEXT, cond TEXT, packing INT, n0 INT, n INT, kind TEXT,
  gap REAL, tries INT, age_h REAL, w REAL, wsum REAL, options INT);
"""
DEFAULT = dict(conditions=['queue', 'base'], frozen=[], window=2.0e-3, base_window=0.05, pool=8, priority='default')

PRIORITY = {}


def priority(name):
    def deco(f):
        PRIORITY[name] = f
        return f
    return deco


@priority('default')
def _default(gap, tries, age_h, kind):
    """10-26 lineage backtest (store parent links since 10-24): big-hit rate per child falls ~2x by gap 1e-5, ~7x by 1e-4, ~30x by
    1e-3 (score units); ~5x over the first ~10 tries of a near-best parent, then a floor.  Age is not used here (it is
    confounded with tries; logged in qclaim)."""
    return (1.0 / (1.0 + gap / 1e-5)) * (0.2 + 0.8 / (1.0 + tries / 2.0))


def cfg(goal):
    return dict(DEFAULT, **(goal.queue or {}))


def targets(goal):
    c = cfg(goal)
    if c.get('ns'):
        return sorted(int(n) for n in c['ns'])
    lo, hi = c['n']
    return list(range(int(lo), int(hi) + 1))


def policy(cond):
    return cond.split('-')[0]


def moves(arms):
    """{donor n: {kind: [arm, ...]}} over the planned arms, generator-led arms left out"""
    from . import goal as G
    M = {}
    for a in arms:
        if str(a.get('chain') or '').strip().startswith('gen.'):
            continue
        M.setdefault(a['n'] + a['offset'], {}).setdefault(G.arm_kind(a['id']), []).append(a)
    return M


def setup(db, goal, t_goal, store):
    """start pools: per condition and n in the range, up to `pool` distinct feasible packings within base_window of the
    best (families outside progress.exclude), added before the goal start"""
    from . import goal as G, families as F
    c = cfg(goal); ex = goal.progress.get('exclude', [])
    t0 = time.mktime(time.strptime(t_goal, '%Y-%m-%dT%H:%M:%S'))
    rows = []
    for n in targets(goal):
        for i, s, _, _ in F.near(store, n, G.kgrid(n) + 0.05, ex, int(c['pool']), window=float(c['base_window']), before=t_goal):
            rows.append((n, i, s))
    best = {}
    for n, i, s in rows:
        best[n] = min(best.get(n, math.inf), s)
    for cond in c['conditions']:
        db.executemany('INSERT OR IGNORE INTO qpool VALUES (?,?,?,?,?,?,NULL)', [(goal.goal, cond, n, i, s, t0) for n, i, s in rows])
        db.executemany('INSERT OR IGNORE INTO qbest VALUES (?,?,?,?,?)', [(goal.goal, cond, n, s, s) for n, s in best.items()])
    return len(rows)


def spent(db, name):
    """CPU-s per condition, finished units at their CPU and running ones at the unit allowance"""
    out = {}
    for cond, c in db.execute("SELECT variant, SUM(COALESCE(cpu, (SELECT unit_cpu_s FROM goal WHERE name=?))) FROM unit "
                              "WHERE goal=? AND variant IS NOT NULL AND status!='lost' GROUP BY variant", (name, name)):
        out[cond] = c or 0.0
    return out


def pick(db, g, goal, arms, rng, now):
    """the next queue option -> dict(kind, n, arm, start, variant, q) or None"""
    from . import goal as G
    c = cfg(goal); name = g['name']
    M = _moves(goal, arms)
    sp = spent(db, name)
    cond = min(c['conditions'], key=lambda x: (sp.get(x, 0.0), rng.random()))
    best = dict(db.execute('SELECT n, MIN(s) FROM qpool WHERE goal=? AND cond=? GROUP BY n', (name, cond)))
    pol = policy(cond)
    tries = {(p, k): t for p, k, t in db.execute('SELECT packing, kind, tries FROM qtry WHERE goal=? AND cond=?', (name, cond))}
    opts = []
    for n, pid, s, t in db.execute('SELECT n, packing, s, t FROM qpool WHERE goal=? AND cond=?', (name, cond)):
        if n not in M or n not in best:
            continue
        gap = max(0.0, (s - best[n]) / math.sqrt(n))
        if pol == 'base':
            if s > best[n] + float(c['base_window']):
                continue
        elif gap > float(c['window']):
            continue
        for kind in M[n]:
            opts.append((pid, n, kind, gap, tries.get((pid, kind), 0), (now - t) / 3600))
    if not opts:
        return None
    if pol == 'base':                                     # uniform parent, then uniform kind at its n
        parents = sorted({(o[0], o[1]) for o in opts})
        pid, n0 = parents[rng.randrange(len(parents))]
        mine = [o for o in opts if o[0] == pid]
        o = mine[rng.randrange(len(mine))]; w = 1.0; wsum = float(len(parents))
    else:
        f = PRIORITY[c['priority']]
        W = [f(gap, tr, age, kind) for _, _, kind, gap, tr, age in opts]
        wsum = sum(W); x = rng.random() * wsum; acc = 0.0
        for o, w in zip(opts, W):
            acc += w
            if acc >= x:
                break
    pid, n0, kind, gap, tr, age = o
    arm = M[n0][kind][rng.randrange(len(M[n0][kind]))]
    return dict(kind='construct' if arm.get('chain') else 'walk', n=arm['n'], arm=arm['id'], parent=None, start=pid,
                variant=cond, q=dict(cond=cond, packing=pid, n0=n0, n=arm['n'], kind=kind, gap=gap, tries=tr, age_h=age,
                                     w=w, wsum=wsum, options=len(opts)))


_M = {}


def _moves(goal, arms):
    key = (goal.goal, len(arms))
    if key not in _M:
        _M[key] = moves(arms)
    return _M[key]


def on_claim(db, name, uid, o):
    q = o['q']
    db.execute('INSERT INTO qtry VALUES (?,?,?,?,1) ON CONFLICT(goal, cond, packing, kind) DO UPDATE SET tries=tries+1',
               (name, q['cond'], q['packing'], q['kind']))
    db.execute('INSERT INTO qclaim VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)', (uid, name, q['cond'], q['packing'], q['n0'], q['n'],
               q['kind'], q['gap'], q['tries'], q['age_h'], q['w'], q['wsum'], q['options']))


def on_record(db, u, res, goal, now):
    """inside record's transaction: the unit's new packings at its n join its condition's pool (within base_window of
    that condition's best); the condition's best moves -> (side gain, score gain)"""
    c = cfg(goal); name, cond, n = u['goal']['name'], u.get('variant'), u['n']
    if not cond or n is None:
        return 0.0, 0.0
    ex = set(goal.progress.get('exclude', []))
    new = [(i, s) for i, s, f in res.get('new', []) if f not in ex]
    r = db.execute('SELECT s FROM qbest WHERE goal=? AND cond=? AND n=?', (name, cond, n)).fetchone()
    b0 = r[0] if r else math.inf
    b1 = min([b0] + [s for _, s in new])
    for i, s in ([] if cond in c.get('frozen', []) else new):
        if s <= b1 + float(c['base_window']):
            db.execute('INSERT OR IGNORE INTO qpool VALUES (?,?,?,?,?,?,?)', (name, cond, n, i, s, now, u['id']))
    if b1 < b0:
        db.execute('INSERT INTO qbest VALUES (?,?,?,?,?) ON CONFLICT(goal, cond, n) DO UPDATE SET s=excluded.s',
                   (name, cond, n, b1 if r is None else b0, b1))
    g = (b0 - b1) if math.isfinite(b0) and b1 < b0 else 0.0
    return g, g / math.sqrt(n)


def report(db, name):
    """per condition: CPU, units, own-best score gain (sum over n of (s0 - s) / sqrt n), n moved, gain curve by CPU"""
    out = [f'## queue conditions ({name})', '', '| condition | units | CPU-h | n moved | score gain | score / CPU-h |', '|---|---|---|---|---|---|']
    conds = [c for (c,) in db.execute('SELECT DISTINCT cond FROM qbest WHERE goal=?', (name,))]
    curves = {}
    for cond in conds:
        U = db.execute("SELECT COUNT(*), COALESCE(SUM(cpu), 0) / 3600 FROM unit WHERE goal=? AND variant=? AND cpu IS NOT NULL",
                       (name, cond)).fetchone()
        R = db.execute('SELECT n, s0, s FROM qbest WHERE goal=? AND cond=? AND s < s0', (name, cond)).fetchall()
        sg = sum((a - b) / math.sqrt(n) for n, a, b in R)
        out.append(f'| {cond} | {U[0]} | {U[1]:.2f} | {len(R)} | {sg:.3e} | {sg / max(U[1], 1e-9):.2e} |')
        cum_c = cum_g = 0.0; pts = []
        for cpu_s, rw in db.execute("SELECT cpu, reward FROM unit WHERE goal=? AND variant=? AND cpu IS NOT NULL ORDER BY ended",
                                    (name, cond)):
            cum_c += cpu_s / 3600; cum_g += rw or 0.0; pts.append((cum_c, cum_g))
        curves[cond] = pts
    out += ['', 'score gain at CPU-h (each condition\'s own best, summed over n):']
    for cond, pts in curves.items():
        marks = []
        for h in (0.25, 0.5, 1, 2, 4, 8, 16):
            v = [g for c, g in pts if c <= h]
            if pts and pts[-1][0] >= h:
                marks.append(f'{h:g} h {v[-1] if v else 0.0:.2e}')
        out.append(f'  {cond}: ' + ', '.join(marks))
    k = db.execute("SELECT q.cond, q.kind, COUNT(*), SUM(u.cpu) / 3600, SUM(u.reward), SUM(u.reward > 0) FROM qclaim q JOIN unit u "
                   "ON u.id = q.unit WHERE q.goal=? AND u.cpu IS NOT NULL GROUP BY q.cond, q.kind ORDER BY q.cond, SUM(u.reward) DESC",
                   (name,)).fetchall()
    out += ['', '| condition | move kind | units | CPU-h | score gain | units with a gain |', '|---|---|---|---|---|---|']
    out += [f'| {c} | {kd} | {m} | {h:.2f} | {r:.2e} | {w} |' for c, kd, m, h, r, w in k]
    return '\n'.join(out)
