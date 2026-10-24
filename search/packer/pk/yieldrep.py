"""Lineage yield report (10-24, explore-loop plan step 1): which generators / sources feed the records and near-records?

For every sub-grid n: the store's best packing (record set R), and distinct (1e-7) packings within +1e-3 (N3) and +1e-2
(N2) of it, at most --cap per n.  Each packing's ancestry is walked through `packing.parent` to its root; every packing
is labelled by its first sighting's source, grouped to a generator key (run dir under runs/, or the source kind:
register, kaggle, watch, converge, ...).  Credit per key:
  made    the key produced the packing itself
  touched the key appears anywhere in its ancestry (packing included)
  origin  the key of its root (no parent)
Public keys are marked; "ours" = any packing first sighted from one of our sources.

  pk.py yield [--cap 50] [--json runs/yield.json]
"""
from __future__ import annotations
import collections, json, math, re

PUBLIC = ('register', 'kaggle', 'kingbird', 'watch', 'pending', 'pool:reg')


def key_of(source: str) -> str:
    if source.startswith('run:runs/'):
        parts = source[len('run:runs/'):].split('/')
        return 'run:' + re.sub(r'\d+$', '', parts[0])            # night1014 -> night, tt3 -> tt, hunt1012 -> hunt
    for p in PUBLIC:
        if source.startswith(p):
            return 'public:' + p.split(':')[0]
    return source.split(':')[0]


def report(st, cap=50, nmax=1000):
    db = st.db
    first = {}
    for pid, src in db.execute("SELECT packing, source FROM sighting s WHERE added = (SELECT min(added) FROM sighting "
                               "WHERE packing = s.packing)"):
        first.setdefault(pid, src)
    parent = dict(db.execute("SELECT id, parent FROM packing WHERE parent IS NOT NULL"))
    feas = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9) AND (grid IS NULL OR grid=0)"
    sets = {'R': [], 'N3': [], 'N2': []}
    for n in range(2, nmax + 1):
        k = math.ceil(math.sqrt(n) - 1e-12)
        rows = db.execute(f"SELECT id, s FROM packing WHERE n=? AND {feas} AND s<? ORDER BY s LIMIT 5000", (n, k - 1e-9)).fetchall()
        if not rows:
            continue
        b = rows[0][1]; sets['R'].append(rows[0][0])
        seen = set(); c3 = c2 = 0
        for pid, s in rows:
            if s > b + 1e-2:
                break
            if round(s, 7) in seen:
                continue
            seen.add(round(s, 7))
            if s <= b + 1e-3 and c3 < cap:
                sets['N3'].append(pid); c3 += 1
            if c2 < cap:
                sets['N2'].append(pid); c2 += 1
    out = {}
    for name, ids in sets.items():
        made, touched, origin = collections.Counter(), collections.Counter(), collections.Counter()
        for pid in ids:
            chain, cur, guard = [], pid, 0
            while cur is not None and guard < 100000:
                chain.append(cur); cur = parent.get(cur); guard += 1
                if cur in chain[-3:]:
                    break
            keys = [key_of(first.get(c, '?')) for c in chain]
            made[keys[0]] += 1; origin[keys[-1]] += 1
            for k_ in set(keys):
                touched[k_] += 1
        out[name] = dict(size=len(ids), made=made.most_common(), touched=touched.most_common(), origin=origin.most_common())
    return out


def fmt(R, top=25):
    lines = []
    for name, title in (('R', 'store bests (one per sub-grid n)'), ('N3', 'within +1e-3 of the best'), ('N2', 'within +1e-2 of the best')):
        d = R[name]
        lines.append(f'== {title}: {d["size"]} packings')
        tb = {k: v for k, v in d['touched']}
        lines.append(f'{"key":28s} {"made":>6s} {"touched":>8s} {"origin":>7s}')
        keys = sorted(tb, key=lambda k: -tb[k])[:top]
        md, og = dict(d['made']), dict(d['origin'])
        for k in keys:
            lines.append(f'{k:28s} {md.get(k, 0):6d} {tb[k]:8d} {og.get(k, 0):7d}')
    return '\n'.join(lines)
