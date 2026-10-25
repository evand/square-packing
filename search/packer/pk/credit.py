"""Lineage credit audit (10-24; Evan: "we're not accurately propagating credit upstream").

Store reads only.  An *event* is a packing that lowered the store best at its n when it was added (feasible, sub-grid).
Each event's gain g = previous best - its side (side units; also g / sqrt n, score units) is paid up its ancestry
(`packing.parent`): ancestor at generation d (the packing itself d = 0) gets weight lambda^d, weights normalized to
sum 1, then aggregated by generator key.  Shrink-only steps (polish, converge, goal-ref) pass credit through: they
are skipped when counting generations, and earn nothing.

Generator keys: `goal:<arm family>/<stage>` for explore-loop jobs (stage = build | walk), `run:<run dir>` for older
runs (pk.yieldrep.key_of), `public:<source>` for imports.  Per-use credit divides by how many store packings the key
made (a popular generator should not win by volume).

  pk.py credit [--since 2026-10-10] [--lam 0.5 0.8 1.0] [--top 25]
"""
from __future__ import annotations
import collections, json, math, os, re

from .yieldrep import key_of

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEAS = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9) AND (grid IS NULL OR grid=0)"
PASS = ('polish', 'converge', 'goal-ref')                        # shrink-only: credit passes through


def _arm_family(goal, job, cache={}):
    k = (goal, job)
    if k not in cache:
        try:
            a = json.load(open(os.path.join(HERE, 'runs', 'goals', goal, job, 'arm.json')))['id']
            cache[k] = re.sub(r'[-+]\d+$', '', a.split(':')[1])
        except Exception:
            cache[k] = '?'
    return cache[k]


def gen_key(source):
    m = re.match(r'run:runs/goals/([^/]+)/(j\d+)/?([^/]*)', source)
    if m:
        stage = {'starts': 'build', 'ex': 'walk', 'polish': 'polish'}.get(m.group(3), m.group(3) or '?')
        return 'polish' if stage == 'polish' else f'goal:{_arm_family(m.group(1), m.group(2))}/{stage}'
    if source.startswith('goal-ref:'):
        return 'goal-ref'
    return key_of(source)


def events(db, since='0000'):
    """[(id, n, s, gain, added)]: packings that lowered the store best at n when added (added >= since)"""
    out, cur_n, best = [], None, None
    for i, n, s, added in db.execute(f"SELECT id, n, s, added FROM packing WHERE {FEAS} ORDER BY n, id"):
        k = math.ceil(math.sqrt(n) - 1e-12)
        if n != cur_n:
            cur_n, best = n, float(k)
        if s < best - 1e-9:
            if added >= since and best < k:                       # the first sub-grid packing at n is not a gain
                out.append((i, n, s, best - s, added))
            best = s
    return out


def audit(st, since='2026-10-10', lams=(0.5, 0.8, 1.0), depth=60):
    db = st.db
    first = {}
    for pid, src in db.execute("SELECT packing, source FROM sighting s WHERE added = (SELECT min(added) FROM sighting "
                               "WHERE packing = s.packing)"):
        first.setdefault(pid, src)
    key = {}
    def k_of(pid):
        if pid not in key:
            key[pid] = gen_key(first.get(pid, '?'))
        return key[pid]
    parent = dict(db.execute("SELECT id, parent FROM packing WHERE parent IS NOT NULL"))
    made = collections.Counter(k_of(p) for p in first)
    E = events(db, since)
    E_rows = E
    res = {lam: collections.defaultdict(lambda: [0.0, 0.0, 0]) for lam in lams}   # key -> [side, score, events]
    chains = []
    for pid, n, s, g, added in E:
        chain, cur, seen = [], pid, set()
        while cur is not None and cur not in seen and len(chain) < depth:
            seen.add(cur)
            kk = k_of(cur)
            if kk not in PASS:
                chain.append(kk)
            cur = parent.get(cur)
        if not chain:
            chain = ['polish']
        chains.append((n, g, chain))
        for lam in lams:
            w = [lam ** d for d in range(len(chain))]; tw = sum(w)
            per = collections.defaultdict(float)
            for kk, wi in zip(chain, w):
                per[kk] += wi / tw
            for kk, f in per.items():
                r = res[lam][kk]; r[0] += f * g; r[1] += f * g / math.sqrt(n); r[2] += 1
    return dict(events=len(E), since=since, made=made, res=res, chains=chains, event_rows=E_rows)


def fmt(A, top=25):
    L = [f"credit audit: {A['events']} store-best events since {A['since']}; credit = gain paid up the lineage "
         f"(lambda per generation, shrink-only steps pass through)"]
    lams = list(A['res'])
    base = A['res'][lams[0]]
    tot = {lam: sum(v[1] for v in A['res'][lam].values()) for lam in lams}
    keys = sorted(set().union(*[A['res'][l].keys() for l in lams]), key=lambda k: -A['res'][lams[-1]].get(k, [0, 0, 0])[1])
    L.append('| key | made | events touched | ' + ' | '.join(f'score share l={l}' for l in lams) +
             ' | score / 1e4 made (l=' + str(lams[-1]) + ') |')
    L.append('|---|---|---|' + '---|' * len(lams) + '---|')
    for kk in keys[:top]:
        sh = [A['res'][l].get(kk, [0, 0, 0])[1] / max(tot[l], 1e-300) for l in lams]
        ev = A['res'][lams[-1]].get(kk, [0, 0, 0])[2]
        per = A['res'][lams[-1]].get(kk, [0, 0, 0])[1] / max(A['made'].get(kk, 0), 1) * 1e4
        L.append(f"| {kk} | {A['made'].get(kk, 0)} | {ev} | " + ' | '.join(f'{x:.1%}' for x in sh) + f' | {per:.2e} |')
    depths = collections.Counter(min(len(c), 10) for _, _, c in A['chains'])
    L.append('')
    L.append('lineage length (non-shrink generations, capped at 10): ' + ', '.join(f'{d}: {c}' for d, c in sorted(depths.items())))
    return '\n'.join(L)


# ------------------------------------------------------------------ online (10-24, Evan: "credit should be online and
# continuous").  Tables are created by rebuild() (the backfill); until then the store's insert hook is a no-op.

LAM = 0.8
TABLES = """
CREATE TABLE IF NOT EXISTS credit_best (n INTEGER PRIMARY KEY, s REAL NOT NULL, packing INTEGER);
CREATE TABLE IF NOT EXISTS credit (key TEXT PRIMARY KEY, made INTEGER NOT NULL DEFAULT 0, events INTEGER NOT NULL DEFAULT 0,
  side REAL NOT NULL DEFAULT 0, score REAL NOT NULL DEFAULT 0, cpu_s REAL NOT NULL DEFAULT 0, updated TEXT);
CREATE TABLE IF NOT EXISTS credit_event (packing INTEGER PRIMARY KEY, n INTEGER NOT NULL, gain REAL NOT NULL,
  chain TEXT NOT NULL, added TEXT NOT NULL);
"""


def ready(db):
    return db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='credit_best'").fetchone() is not None


def _now():
    import time
    return time.strftime('%Y-%m-%dT%H:%M:%S')


def _first_source(db, pid):
    r = db.execute('SELECT source FROM sighting WHERE packing=? ORDER BY added LIMIT 1', (pid,)).fetchone()
    return r[0] if r else '?'


def _chain(db, pid, source=None, depth=60):
    """non-shrink generator keys from pid up its ancestry"""
    out, cur, seen = [], pid, set()
    while cur is not None and cur not in seen and len(out) < depth:
        seen.add(cur)
        kk = gen_key(source if (cur == pid and source) else _first_source(db, cur))
        if kk not in PASS:
            out.append(kk)
        r = db.execute('SELECT parent FROM packing WHERE id=?', (cur,)).fetchone()
        cur = r[0] if r else None
    return out or ['polish']


def _pay(db, chain, n, gain, lam=LAM):
    w = [lam ** d for d in range(len(chain))]; tw = sum(w)
    per = collections.defaultdict(float)
    for kk, wi in zip(chain, w):
        per[kk] += wi / tw
    for kk, f in per.items():
        db.execute('INSERT INTO credit(key, events, side, score, updated) VALUES (?, 1, ?, ?, ?) ON CONFLICT(key) DO UPDATE '
                   'SET events = events + 1, side = side + excluded.side, score = score + excluded.score, '
                   'updated = excluded.updated', (kk, f * gain, f * gain / math.sqrt(n), _now()))


def on_add(db, pid, n, s, source, feasible):
    """the store's insert hook (a new packing): count it for its generator; if it lowers the best at n, pay the gain up
    its ancestry.  No-op until rebuild() has run."""
    if not ready(db):
        return
    kk = gen_key(source)
    db.execute('INSERT INTO credit(key, made, updated) VALUES (?, 1, ?) ON CONFLICT(key) DO UPDATE SET made = made + 1, '
               'updated = excluded.updated', (kk, _now()))
    k = math.ceil(math.sqrt(n) - 1e-12)
    if not feasible or s >= k - 1e-9:
        return
    r = db.execute('SELECT s FROM credit_best WHERE n=?', (n,)).fetchone()
    if r is None:
        db.execute('INSERT INTO credit_best VALUES (?, ?, ?)', (n, s, pid)); return   # first sub-grid packing: no gain
    if s < r[0] - 1e-9:
        chain = _chain(db, pid, source)
        _pay(db, chain, n, r[0] - s)
        db.execute('INSERT OR REPLACE INTO credit_event VALUES (?, ?, ?, ?, ?)', (pid, n, r[0] - s, json.dumps(chain), _now()))
        db.execute('UPDATE credit_best SET s=?, packing=? WHERE n=?', (s, pid, n))


def add_cpu(db, key, cpu_s):
    if ready(db):
        db.execute('INSERT INTO credit(key, cpu_s, updated) VALUES (?, ?, ?) ON CONFLICT(key) DO UPDATE SET '
                   'cpu_s = cpu_s + excluded.cpu_s, updated = excluded.updated', (key, cpu_s, _now()))
        db.commit()


def rebuild(st, lam=LAM, log=print):
    """backfill from scratch (also the consistency check): best per n, events since the store began, made per key"""
    db = st.db
    db.executescript('DROP TABLE IF EXISTS credit_best; DROP TABLE IF EXISTS credit; DROP TABLE IF EXISTS credit_event;' + TABLES)
    A = audit(st, since='0000', lams=(lam,))
    for kk, c in A['made'].items():
        db.execute('INSERT INTO credit(key, made, updated) VALUES (?, ?, ?)', (kk, c, _now()))
    for kk, (side, score, ev) in A['res'][lam].items():
        db.execute('INSERT INTO credit(key, events, side, score, updated) VALUES (?, ?, ?, ?, ?) ON CONFLICT(key) DO UPDATE '
                   'SET events = excluded.events, side = excluded.side, score = excluded.score', (kk, ev, side, score, _now()))
    for (pid, n, s, g, added), (_, _, chain) in zip(A['event_rows'], A['chains']):
        db.execute('INSERT OR REPLACE INTO credit_event VALUES (?, ?, ?, ?, ?)', (pid, n, g, json.dumps(chain), added))
    for n, s, i in db.execute(f"SELECT n, min(s), id FROM packing WHERE {FEAS} GROUP BY n").fetchall():
        k = math.ceil(math.sqrt(n) - 1e-12)
        if s < k - 1e-9:
            db.execute('INSERT OR REPLACE INTO credit_best VALUES (?, ?, ?)', (n, s, i))
    db.commit()
    log(f"credit rebuilt: {A['events']} events, {len(A['made'])} keys")
    return A


def table(db, top=25):
    rows = db.execute('SELECT key, made, events, side, score, cpu_s FROM credit ORDER BY score DESC LIMIT ?', (top,)).fetchall()
    L = ['| key | made | events | score credit | per 1e4 made | per CPU-h |', '|---|---|---|---|---|---|']
    for kk, m, ev, side, sc, cpu in rows:
        L.append(f"| {kk} | {m} | {ev} | {sc:.4f} | {sc / max(m, 1) * 1e4:.2e} | {(sc / (cpu / 3600)) if cpu else float('nan'):.2e} |")
    return '\n'.join(L)
