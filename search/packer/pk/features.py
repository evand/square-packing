"""Per-packing features for parent filters (11-02): `requires:` on a campaign move (pk/events.py) and `pk.py try
--requires` (pk/lab.py).  A feature is f(s, sq) -> a JSON value, registered by a plugin with @feature(name, version);
values are cached per store id in the store table `feature` (and in-process), recomputed when the version changes.

A requirement is a string: 'name' (truthy), '!name' (falsy), 'name=v' (string equality), 'name>=x' / 'name<=x' /
'name>x' / 'name<x' (numeric).  ok(pid, reqs) is True when every one holds (a missing feature fails the requirement).
"""
from __future__ import annotations
import json, re, sqlite3
import numpy as np

FEATURES: dict = {}                    # name -> f, with f.version
_MEM: dict = {}                        # (pid, name, version) -> value
_REQ = re.compile(r'^\s*(!?)\s*([A-Za-z_][\w.]*)\s*(?:(>=|<=|=|>|<)\s*(.+?))?\s*$')
SCHEMA = 'CREATE TABLE IF NOT EXISTS feature(pid INTEGER, name TEXT, version TEXT, value TEXT, PRIMARY KEY(pid, name))'


def feature(name, version=1):
    def deco(f):
        f.feature_name, f.version = name, str(version)
        FEATURES[name] = f
        return f
    return deco


def _path():
    from .store import DEFAULT
    return DEFAULT


def get(pids, names, path=None):
    """{pid: {name: value}} for store ids and feature names (computed once, cached in the store)"""
    import pk.stages  # noqa: F401  (plugins register their features when the stage registry loads)
    path = path or _path()
    out = {p: {} for p in pids}
    todo = []
    for p in pids:
        for nm in names:
            f = FEATURES.get(nm)
            if f is None:
                raise KeyError(f'unknown feature {nm!r} (registered: {sorted(FEATURES)})')
            k = (p, nm, f.version)
            if k in _MEM:
                out[p][nm] = _MEM[k]
            else:
                todo.append(k)
    if not todo:
        return out
    db = sqlite3.connect(path, timeout=5)
    try:
        db.execute(SCHEMA)
        new = []
        for p, nm, v in todo:
            r = db.execute('SELECT version, value FROM feature WHERE pid=? AND name=?', (p, nm)).fetchone()
            if r and r[0] == v:
                val = json.loads(r[1])
            else:
                row = db.execute('SELECT s, coords FROM packing WHERE id=?', (p,)).fetchone()
                if row is None:
                    continue
                val = FEATURES[nm](float(row[0]), np.frombuffer(row[1], '<f8').reshape(-1, 3).copy())
                new.append((p, nm, v, json.dumps(val)))
            _MEM[(p, nm, v)] = val
            out[p][nm] = val
        if new:
            try:
                db.executemany('INSERT OR REPLACE INTO feature VALUES (?,?,?,?)', new); db.commit()
            except sqlite3.OperationalError:            # store busy: the in-process cache still holds them
                pass
    finally:
        db.close()
    return out


def parse(req):
    m = _REQ.match(str(req))
    if not m:
        raise ValueError(f'bad requirement {req!r}')
    neg, name, op, val = m.groups()
    return bool(neg), name, op, val


def holds(vals, req):
    neg, name, op, val = parse(req)
    if name not in vals:
        return False
    x = vals[name]
    if op is None:
        r = bool(x)
    elif op == '=':
        r = str(x) == val
    else:
        try:
            x, y = float(x), float(val)
        except (TypeError, ValueError):
            return False
        r = dict([('>=', x >= y), ('<=', x <= y), ('>', x > y), ('<', x < y)])[op]
    return r != neg


def ok(pid, reqs, path=None):
    names = sorted({parse(r)[1] for r in reqs})
    vals = get([pid], names, path)[pid]
    return all(holds(vals, r) for r in reqs)


def filter_ids(pids, reqs, path=None):
    if not reqs:
        return list(pids)
    names = sorted({parse(r)[1] for r in reqs})
    V = get(list(pids), names, path)
    return [p for p in pids if all(holds(V[p], r) for r in reqs)]
