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
    for i, n, s, added in db.execute(f"SELECT id, n, s, added FROM packing WHERE {FEAS} ORDER BY n, added, s"):
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
    return dict(events=len(E), since=since, made=made, res=res, chains=chains)


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
