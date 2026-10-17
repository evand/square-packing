"""Packing store: every packing we touch (register, pending issues, other people's repos, our runs and candidates) in one
sqlite file, deduplicated up to box symmetry, with provenance (sightings), status and tags.

Tables
  packing(id, hash UNIQUE, n, s, s_str, coords BLOB f64 (n,3), status, s_exact, fp, ntilted, grid, pen, parent, tags, added)
  sighting(packing, source, finder, ref, kind, added)   -- where a packing was seen (many per packing)
  frontier(n, s, s_str, finder, source, ref, pending, updated)  -- best known per n, one row per (n, source)
  basin(n, key, cls, S, status, added)  -- exactsolve verdict per screened f64 side key (f'{s:.9f}'): a lookup cache,
                                           not a minimum (pk/identity.py)
  minimum(n, key, S, cls, status, src, pid, added, updated)  -- one row per exact minimum (10-11): key = identity.min_key(S)
                                           (1e-17), cls certified | weak; census / "new minimum" counts read this table
  copy(n, side_key, min_key, via, added)  -- a screened side that descends / solves to a minimum (copies of one minimum)
  expansion(packing, n, mechanism, params, run, cpu_s, nprop, best_s, frontier_s, added)  -- ledger (10-14): CPU spent
                                           expanding a store packing (one row per lineage root x run) and the best side its
                                           descendants reached; every driver writes it, parent pickers read it
  nearmiss(packing, n, s, frontier, gap, origin, status, added, updated, best_s, best_pid, cpu_s, note)  -- near-miss
                                           queue (10-14, pk/nearmiss.py): our packings within ~1e-3 of frontier_now(n);
                                           status queued / running / done / candidate / failed; drained by `pk.py finish`
  packing.cluster  -- region label for parent sampling (10-14; definition from analysis/similarity110.md), NULL = unset

status: raw (imported, f64-feasible within tol) < screened < polished < certified; or not-min / infeasible / unresolved.
"""
from __future__ import annotations
import collections, datetime, json, math, os, sqlite3
import numpy as np
from .packing import Packing

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.environ.get('PK_STORE', os.path.join(HERE, '..', 'runs', 'store.sqlite'))
RANK = dict(infeasible=-1, raw=0, screened=1, polished=2, unresolved=2, **{'not-min': 3}, weak=4, certified=5)

SCHEMA = """
CREATE TABLE IF NOT EXISTS packing (
  id INTEGER PRIMARY KEY, hash TEXT UNIQUE NOT NULL, n INTEGER NOT NULL, s REAL NOT NULL, s_str TEXT,
  coords BLOB NOT NULL, status TEXT NOT NULL DEFAULT 'raw', s_exact TEXT, fp TEXT, ntilted INTEGER, grid INTEGER,
  pen REAL, parent INTEGER, tags TEXT NOT NULL DEFAULT '[]', added TEXT NOT NULL);
CREATE INDEX IF NOT EXISTS packing_ns ON packing(n, s);
CREATE TABLE IF NOT EXISTS sighting (
  packing INTEGER NOT NULL REFERENCES packing(id), source TEXT NOT NULL, finder TEXT, ref TEXT, kind TEXT, added TEXT,
  UNIQUE(packing, source, ref));
CREATE INDEX IF NOT EXISTS sighting_src ON sighting(source);
CREATE TABLE IF NOT EXISTS basin (
  n INTEGER NOT NULL, key TEXT NOT NULL, cls TEXT NOT NULL, S TEXT, status TEXT, added TEXT, UNIQUE(n, key));
CREATE TABLE IF NOT EXISTS minimum (
  n INTEGER NOT NULL, key TEXT NOT NULL, S TEXT NOT NULL, cls TEXT NOT NULL, status TEXT, src TEXT, pid INTEGER,
  added TEXT, updated TEXT, UNIQUE(n, key));
CREATE TABLE IF NOT EXISTS copy (
  n INTEGER NOT NULL, side_key TEXT NOT NULL, min_key TEXT NOT NULL, via TEXT, added TEXT, UNIQUE(n, side_key, min_key));
CREATE TABLE IF NOT EXISTS piece (
  id INTEGER PRIMARY KEY, hash TEXT UNIQUE NOT NULL, m INTEGER NOT NULL, d REAL NOT NULL, h REAL NOT NULL,
  coords BLOB NOT NULL, status TEXT NOT NULL DEFAULT 'raw', pen REAL, source TEXT NOT NULL, ref TEXT,
  tags TEXT NOT NULL DEFAULT '[]', added TEXT NOT NULL);
CREATE INDEX IF NOT EXISTS piece_mdh ON piece(m, d, h);
CREATE TABLE IF NOT EXISTS expansion (
  id INTEGER PRIMARY KEY, packing INTEGER NOT NULL REFERENCES packing(id), n INTEGER NOT NULL, mechanism TEXT NOT NULL,
  params TEXT, run TEXT, cpu_s REAL NOT NULL, nprop INTEGER, best_s REAL, frontier_s REAL, added TEXT NOT NULL);
CREATE INDEX IF NOT EXISTS expansion_pk ON expansion(packing);
CREATE INDEX IF NOT EXISTS expansion_n ON expansion(n);
CREATE TABLE IF NOT EXISTS nearmiss (
  packing INTEGER PRIMARY KEY REFERENCES packing(id), n INTEGER NOT NULL, s REAL NOT NULL, frontier REAL, gap REAL,
  origin TEXT, status TEXT NOT NULL DEFAULT 'queued', added TEXT NOT NULL, updated TEXT, best_s REAL, best_pid INTEGER,
  cpu_s REAL, note TEXT);
CREATE INDEX IF NOT EXISTS nearmiss_st ON nearmiss(status, gap);
CREATE TABLE IF NOT EXISTS sync (
  key TEXT PRIMARY KEY, value TEXT, at TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS frontier (
  n INTEGER NOT NULL, s REAL NOT NULL, s_str TEXT, finder TEXT, source TEXT NOT NULL, ref TEXT, pending INTEGER NOT NULL,
  packing INTEGER, updated TEXT, UNIQUE(n, source));
"""


# Sighting sources that mean "someone published this packing" (10-17: the one definition; fresh, nearmiss, ingest,
# eventbench and distinct used to carry their own copies, eventbench's without watch:, and every one let pool:ours, our
# own seed pool, count as public).
PUBLIC = ('register', 'pending:', 'kingbird', 'pool:', 'watch:')
PRIVATE = ('pool:ours',)


def is_public(source):
    return source.startswith(PUBLIC) and source not in PRIVATE


def public_sql(alias='g'):
    """SQL condition: sighting `alias` is a public source"""
    return ('((' + ' OR '.join(f"{alias}.source LIKE '{p}%'" for p in PUBLIC) + ')'
            + ''.join(f" AND {alias}.source != '{x}'" for x in PRIVATE) + ')')


def now():
    return datetime.datetime.now().isoformat(timespec='seconds')


class Store:
    def __init__(self, path=DEFAULT):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        self.db = sqlite3.connect(path, timeout=60)
        self.db.row_factory = sqlite3.Row
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.executescript(SCHEMA)
        if 'cluster' not in {r[1] for r in self.db.execute('PRAGMA table_info(packing)')}:   # 10-14 migration
            self.db.execute('ALTER TABLE packing ADD COLUMN cluster TEXT')
            self.db.commit()

    def close(self):
        self.db.commit(); self.db.close()

    # ---------- write ----------
    def add(self, p: Packing, source: str, finder=None, ref=None, kind=None, status='raw', parent=None, tags=(),
            check=True, tol=1e-9, fp=None, commit=True, pre=None, s_exact=None):
        """Insert (or find) a packing; record the sighting.  Returns (id, created).  Infeasible (pen > tol) packings are
        stored with status 'infeasible' so the rejection is remembered."""
        pre = pre or {}                                  # precomputed hash / grid / ntilted / pen (parallel importers)
        h = pre.get('hash') or p.canon_hash()
        row = self.db.execute('SELECT id, status, tags FROM packing WHERE hash=?', (h,)).fetchone()
        created = row is None
        if created:
            pen = pre['pen'] if 'pen' in pre else (p.max_pen() if check else None)
            st = 'infeasible' if (pen is not None and pen > tol) else status
            cur = self.db.execute(
                'INSERT INTO packing(hash,n,s,s_str,coords,status,s_exact,fp,ntilted,grid,pen,parent,tags,added) '
                'VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                (h, p.n, float(p.s), p.s_str, p.sq.astype('<f8').tobytes(), st, s_exact, fp,
                 pre.get('ntilted', p.ntilted()) if 'ntilted' in pre else p.ntilted(),
                 pre['grid'] if 'grid' in pre else p.grid_lines(), pen, parent, json.dumps(sorted(set(tags))), now()))
            pid = cur.lastrowid
        else:
            pid = row['id']
            if RANK.get(status, 0) > RANK.get(row['status'], 0) and row['status'] != 'infeasible':
                self.db.execute('UPDATE packing SET status=? WHERE id=?', (status, pid))
            if s_exact:
                self.db.execute('UPDATE packing SET s_exact=? WHERE id=? AND s_exact IS NULL', (s_exact, pid))
            if tags:
                t = sorted(set(json.loads(row['tags'])) | set(tags))
                self.db.execute('UPDATE packing SET tags=? WHERE id=?', (json.dumps(t), pid))
        self.db.execute('INSERT OR IGNORE INTO sighting VALUES (?,?,?,?,?,?)', (pid, source, finder, ref, kind, now()))
        if commit:
            self.db.commit()
        return pid, created

    def set(self, pid, **kw):
        if 'tags' in kw:
            kw['tags'] = json.dumps(sorted(set(kw['tags'])))
        cols = ', '.join(f'{k}=?' for k in kw)
        self.db.execute(f'UPDATE packing SET {cols} WHERE id=?', (*kw.values(), pid)); self.db.commit()

    def tag(self, pid, *tags, remove=()):
        t = set(json.loads(self.db.execute('SELECT tags FROM packing WHERE id=?', (pid,)).fetchone()[0]))
        self.set(pid, tags=(t | set(tags)) - set(remove))

    def set_frontier(self, n, s, s_str, finder, source, ref, pending, packing=None):
        self.db.execute('INSERT OR REPLACE INTO frontier VALUES (?,?,?,?,?,?,?,?,?)',
                        (n, float(s), s_str, finder, source, ref, int(pending), packing, now()))

    def set_basin(self, n, key, cls, S=None, status=None):
        """Certification of a basin, keyed by the polished side to 1e-9 (f'{s:.9f}')."""
        self.db.execute('INSERT OR REPLACE INTO basin VALUES (?,?,?,?,?,?)', (n, key, cls, S, status, now()))
        self.db.commit()

    def basin(self, n, key):
        return self.db.execute('SELECT * FROM basin WHERE n=? AND key=?', (n, key)).fetchone()

    # ---------- rectangle pieces: m squares in (h + d) x h (Packing with dx = d), for composition ----------
    def add_piece(self, p: Packing, source: str, ref=None, tags=(), tol=1e-9, commit=True):
        """Insert (or find) a piece; infeasible ones are stored as such.  Every distinct local minimum is worth keeping
        (seams differ; Evan 10-12).  Returns (id, created)."""
        h = p.canon_hash()
        row = self.db.execute('SELECT id, tags FROM piece WHERE hash=?', (h,)).fetchone()
        if row:
            if tags:
                t = sorted(set(json.loads(row['tags'])) | set(tags))
                self.db.execute('UPDATE piece SET tags=? WHERE id=?', (json.dumps(t), row['id']))
            if commit:
                self.db.commit()
            return row['id'], False
        pen = p.max_pen()
        cur = self.db.execute('INSERT INTO piece(hash,m,d,h,coords,status,pen,source,ref,tags,added) '
                              'VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                              (h, p.n, float(p.dx), float(p.s), p.sq.astype('<f8').tobytes(),
                               'infeasible' if pen > tol else 'raw', pen, source, ref, json.dumps(sorted(set(tags))), now()))
        if commit:
            self.db.commit()
        return cur.lastrowid, True

    def pieces(self, m=None, d=None, mmin=None, source=None, feasible=True):
        """Pieces as Packings (s = h, dx = d; meta: id, source, tags), ordered by (m, d, h)."""
        w, a = [], []
        for col, op, v in (('m', '=', m), ('d', '=', d), ('m', '>=', mmin)):
            if v is not None:
                w.append(f'{col}{op}?'); a.append(v)
        if source:
            w.append('source LIKE ?'); a.append(source)
        if feasible:
            w.append("status!='infeasible'")
        q = 'SELECT * FROM piece' + (' WHERE ' + ' AND '.join(w) if w else '') + ' ORDER BY m, d, h'
        out = []
        for r in self.db.execute(q, a):
            out.append(Packing(r['h'], np.frombuffer(r['coords'], '<f8').reshape(-1, 3), dx=r['d'],
                               meta=dict(id=r['id'], source=r['source'], ref=r['ref'], tags=json.loads(r['tags']),
                                         added=r['added'])))
        return out

    # ---------- minima (pk/identity.py) ----------
    def find_minimum(self, n, S):
        """The minimum row whose key is S's key or a neighbour (1e-17), or None."""
        from .identity import same_min_keys
        ks = same_min_keys(S)
        return self.db.execute(f'SELECT * FROM minimum WHERE n=? AND key IN (?,?,?)', (n, *ks)).fetchone()

    def add_minimum(self, n, S, cls, status=None, src=None, pid=None, commit=True):
        """Record an exact minimum (cls certified | weak).  Merges with an existing row within 1e-17: keeps the longer S
        string, upgrades weak -> certified (never down).  Returns (key, is_new)."""
        from .identity import min_key, RANK as IR, CERTIFIED, WEAK
        assert cls in (CERTIFIED, WEAK), cls
        r = self.find_minimum(n, S)
        t = now()
        if r is None:
            k = min_key(S)
            self.db.execute('INSERT INTO minimum VALUES (?,?,?,?,?,?,?,?,?)', (n, k, str(S), cls, status, src, pid, t, t))
            out = (k, True)
        else:
            S2 = str(S) if len(str(S)) > len(r['S']) else r['S']
            up = IR[cls] > IR[r['cls']]
            self.db.execute('UPDATE minimum SET S=?, cls=?, status=?, src=?, updated=? WHERE n=? AND key=?',
                            (S2, cls if up else r['cls'], status if up else r['status'], src if up else r['src'], t, n, r['key']))
            out = (r['key'], False)
        if commit:
            self.db.commit()
        return out

    def record_verdict(self, n, side_key, r, via, pid=None):
        """One exactsolve verdict r (census_exact.classify dict) for a screened side: basin cache row (side key) and,
        for certified / weak, the minimum row + the copy link side_key -> minimum.  Returns the minimum key or None."""
        from .identity import CERTIFIED, WEAK
        cls = r.get('cls', 'unresolved')
        self.db.execute('INSERT OR REPLACE INTO basin VALUES (?,?,?,?,?,?)', (n, side_key, cls, r.get('S'), r.get('status'), now()))
        mk = None
        if cls in (CERTIFIED, WEAK) and r.get('S'):
            mk, _ = self.add_minimum(n, r['S'], cls, r.get('status'), src=via, pid=pid, commit=False)
            self.add_copy(n, side_key, mk, via, commit=False)
        self.db.commit()
        return mk

    def add_copy(self, n, side_key, mkey, via=None, commit=True):
        self.db.execute('INSERT OR IGNORE INTO copy VALUES (?,?,?,?,?)', (n, side_key, mkey, via, now()))
        if commit:
            self.db.commit()

    def side_minima(self, n):
        """{side_key: (min_key, cls, S float)} via the copy table (a side key linked to two minima, a 1e-9 collision, maps
        to the lower one).  Reports count distinct minima through this instead of side keys (10-11)."""
        M = {r['key']: (r['cls'], float(r['S'][:20])) for r in self.db.execute('SELECT key, cls, S FROM minimum WHERE n=?', (n,))}
        out = {}
        for sk, mk in self.db.execute('SELECT side_key, min_key FROM copy WHERE n=?', (n,)):
            if mk in M and (sk not in out or M[mk][1] < out[sk][2]):
                out[sk] = (mk, M[mk][0], M[mk][1])
        return out

    def minima(self, n, cls=None):
        q, a = 'SELECT * FROM minimum WHERE n=?', [n]
        if cls:
            cl = [cls] if isinstance(cls, str) else list(cls)
            q += f' AND cls IN ({",".join("?" * len(cl))})'; a += cl
        return self.db.execute(q + ' ORDER BY key', a).fetchall()

    # ---------- read ----------
    def get(self, pid) -> Packing:
        r = self.db.execute('SELECT * FROM packing WHERE id=?', (pid,)).fetchone()
        if r is None:
            raise KeyError(pid)
        return self._pk(r)

    def _pk(self, r) -> Packing:
        sq = np.frombuffer(r['coords'], '<f8').reshape(-1, 3)
        return Packing(r['s'], sq, r['s_str'], meta=dict(id=r['id'], status=r['status'], n=r['n'], fp=r['fp'],
                                                         tags=json.loads(r['tags']), parent=r['parent']))

    def query(self, n=None, below=None, above=None, status=None, source=None, tag=None, nongrid=False, limit=None,
              order='s', feasible=True):
        """Rows (sqlite Row, no coords) matching the filters, ordered by side."""
        w, a = [], []
        if n is not None:
            w.append('p.n=?'); a.append(n)
        if below is not None:
            w.append('p.s<?'); a.append(below)
        if above is not None:
            w.append('p.s>?'); a.append(above)
        if status:
            st = [status] if isinstance(status, str) else list(status)
            w.append(f'p.status IN ({",".join("?" * len(st))})'); a += st
        elif feasible:
            w.append("p.status!='infeasible'")
        if nongrid:
            w.append('p.grid=0')
        if tag:
            w.append('p.tags LIKE ?'); a.append(f'%"{tag}"%')
        if source:
            w.append('p.id IN (SELECT packing FROM sighting WHERE source LIKE ?)'); a.append(source)
        q = ('SELECT p.id,p.n,p.s,p.s_str,p.status,p.s_exact,p.ntilted,p.grid,p.pen,p.parent,p.tags,p.fp,'
             '(SELECT group_concat(DISTINCT source) FROM sighting WHERE packing=p.id) AS sources,'
             '(SELECT group_concat(DISTINCT finder) FROM sighting WHERE packing=p.id) AS finders FROM packing p')
        if w:
            q += ' WHERE ' + ' AND '.join(w)
        q += f' ORDER BY {order}'
        if limit:
            q += f' LIMIT {int(limit)}'
        return self.db.execute(q, a).fetchall()

    def best(self, n, **kw) -> Packing | None:
        r = self.query(n=n, limit=1, **kw)
        return self.get(r[0]['id']) if r else None

    def claims(self, n=None, pending=True):
        """{n: lowest claim row}: the `frontier` table is the CLAIMS table (10-17): one row per (n, source), from the
        register (pending=0) and issues (pending=1), with the packing that backs it when one matched."""
        w = ([] if pending else ['pending=0']) + ([] if n is None else ['n=?'])
        q = 'SELECT * FROM frontier' + (' WHERE ' + ' AND '.join(w) if w else '')
        out = {}
        for r in self.db.execute(q, (n,) if n is not None else ()):
            if r['n'] not in out or r['s'] < out[r['n']]['s']:
                out[r['n']] = r
        return out

    def public_best(self, n=None):
        """{n: (s, packing id, source, finder)}: the best feasible packing someone has published (sighting from a
        PUBLIC source), per n"""
        q = ("SELECT p.n, p.s, p.id, g.source, g.finder, g.kind, g.added FROM sighting g JOIN packing p ON p.id = g.packing "
             f"WHERE {public_sql('g')} AND p.status != 'infeasible' AND (p.pen IS NULL OR p.pen <= 1e-9)"
             + (' AND p.n = ?' if n is not None else ''))
        # credit: the sighting that claims it, then the register / a watched repo / the catalogue, a seed pool, and a
        # file that merely sat in some issue's linked dir ('context') last; earliest first
        rank = {'claim': 0, 'record': 0, 'watch': 1, 'catalogue-current': 1, 'seed-pool': 2, 'context': 3}
        out, key = {}, {}
        for r in self.db.execute(q, (n,) if n is not None else ()):
            k = (r['s'], rank.get(r['kind'], 2), r['added'] or '')
            if r['n'] not in out or k < key[r['n']]:
                out[r['n']], key[r['n']] = (r['s'], r['id'], r['source'], r['finder']), k
        return out

    def frontier(self, n=None, pending=True):
        """{n: row} : the best known side at n (10-17: a derived view).  A packing is its own proof, so the frontier is
        the best feasible PUBLISHED packing, or a claim no published packing backs (row 'verified' False: an assertion
        we can't check; it still counts, so a record check stays conservative).  Rows are dicts with n, s, s_str,
        finder, source, ref, pending, packing, verified.  pending=False: the register's own rows only (claims table)."""
        if not pending:
            return {k: dict(r, verified=r['packing'] is not None) for k, r in self.claims(n, pending=False).items()}
        best = self.public_best(n)
        C = collections.defaultdict(list)
        for r in self.db.execute('SELECT * FROM frontier' + (' WHERE n=?' if n is not None else ''),
                                 (n,) if n is not None else ()):
            C[r['n']].append(r)
        out = {}
        for k in set(best) | set(C):
            cands = []
            if k in best:
                s, pid, src, fnd = best[k]
                m = [r for r in C.get(k, ()) if r['packing'] == pid] or [r for r in C.get(k, ()) if abs(r['s'] - s) < 1e-6]
                c = min(m, key=lambda r: (r['pending'], r['s'])) if m else None
                cands.append(dict(n=k, s=s, s_str=None, finder=(c['finder'] if c else fnd), source=(c['source'] if c else src),
                                  ref=(c['ref'] if c else None), pending=(c['pending'] if c else int(src != 'register')),
                                  packing=pid, verified=True))
            for r in C.get(k, ()):
                if k not in best or best[k][0] > r['s'] + 1e-9:          # a claim nothing published backs
                    cands.append(dict(r, verified=False))
            out[k] = min(cands, key=lambda d: d['s'])
        return out

    def unverified(self, n=None):
        """claim rows no published packing backs (side above every public packing at n is not 'unverified': it is
        just superseded)"""
        best = self.public_best(n)
        q = 'SELECT * FROM frontier' + (' WHERE n=?' if n is not None else '')
        return [r for r in self.db.execute(q, (n,) if n is not None else ())
                if r['n'] not in best or best[r['n']][0] > r['s'] + 1e-9]

    def find_hash(self, p: Packing):
        """store id of a packing equal to p up to box symmetry (canonical hash), else None"""
        r = self.db.execute('SELECT id FROM packing WHERE hash=?', (p.canon_hash(),)).fetchone()
        return r['id'] if r else None

    def record_expansion(self, pid, n, mechanism, cpu_s, params=None, run=None, nprop=None, best_s=None, frontier_s=None,
                         commit=True):
        """ledger row: `cpu_s` measured CPU spent expanding store packing `pid` (its lineage in one run) by `mechanism`"""
        self.db.execute('INSERT INTO expansion(packing,n,mechanism,params,run,cpu_s,nprop,best_s,frontier_s,added) '
                        'VALUES (?,?,?,?,?,?,?,?,?,?)', (pid, n, mechanism, json.dumps(params) if params is not None else None,
                                                         run, float(cpu_s), nprop, best_s, frontier_s, now()))
        if commit:
            self.db.commit()

    def expansions(self, n=None, pid=None):
        """{packing id: (total cpu_s, rows, best_s)} from the ledger"""
        q, a = 'SELECT packing, sum(cpu_s) c, count(*) k, min(best_s) b FROM expansion', []
        w = []
        if n is not None:
            w.append('n=?'); a.append(n)
        if pid is not None:
            w.append('packing=?'); a.append(pid)
        q += (' WHERE ' + ' AND '.join(w) if w else '') + ' GROUP BY packing'
        return {r['packing']: (r['c'], r['k'], r['b']) for r in self.db.execute(q, a)}

    # 10-17: what each outside source was last synced to (register main sha, issue update cursor, watched repo HEADs):
    # record checks read the store plus a cheap probe of these (watch.probe), not the sources themselves
    def set_sync(self, key, value, commit=True):
        self.db.execute('INSERT INTO sync(key,value,at) VALUES (?,?,?) ON CONFLICT(key) DO UPDATE SET value=excluded.value, '
                        'at=excluded.at', (key, None if value is None else str(value), now()))
        if commit:
            self.db.commit()

    def get_sync(self, key=None):
        """{key: (value, at)} for keys starting with `key` (all if None)"""
        rows = self.db.execute('SELECT * FROM sync' + (' WHERE key LIKE ?' if key else ''), ((key + '%',) if key else ()))
        return {r['key']: (r['value'], r['at']) for r in rows}

    def sightings(self, pid):
        return self.db.execute('SELECT * FROM sighting WHERE packing=? ORDER BY added', (pid,)).fetchall()
