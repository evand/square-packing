"""What the store keeps, and how much of it (10-20; Evan: "some restrictions on what goes into the store ... available
as intake rules and pruning rules").  One classifier, two uses:

* intake (pk/sink, when a run's results are imported): drop = not stored at all (counted in the run's yield only);
* pruning (`pk.py prune`): an existing row is never deleted (its hash is dedupe memory, its parent link is lineage);
  it is kept `full` (coordinates) or turned into a `stub` (hash, n, side, fingerprint, parent, tags, sightings; no
  coordinates).

Tiers: full > stub > drop.  Tests, in order (the first that decides wins; `reasons` names it):
  protected  public sighting (store.PUBLIC), status certified / weak, an exact side, a near-miss queue row, the
             frontier's packing, tag in cfg.protect_tags                                         -> full
  infeasible                                                                                     -> stub
  n range    n > cfg.nmax and n not in cfg.extra_n                                               -> stub (drop at intake)
  grid       a full line / staggered chain of axis squares forces side >= ceil(sqrt n)           -> stub (drop at intake)
  distance   area excess a = s^2 - F^2 (F = min(frontier, our best); a in unit squares, comparable across n):
             screened (unpolished): a <= cfg.screened_full -> full, else stub (drop at intake beyond cfg.drop)
             polished and better:   a <= cfg.full          -> full, else stub (drop at intake beyond cfg.drop)
  duplicate  handled by the store itself (same canon hash -> a sighting on the existing row)

Targets (a campaign's n) get their own thresholds via cfg.targets {n: {full: ..., ...}}.
Config: DEFAULT below, overridden by pk/policy.yaml if present.
"""
from __future__ import annotations
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = dict(
    nmax=324,                 # beyond: only n listed in extra_n are hunted
    extra_n=[],
    full=1.0,                 # polished rows kept with coordinates up to this area excess (unit squares)
    screened_full=0.1,        # unpolished (screened) rows: tighter
    drop=8.0,                 # intake only: beyond this, not stored at all
    protect_tags=['kingbird', 'record', 'finish', 'candidate'],
    targets={},               # {n: {full: .., screened_full: .., drop: ..}}
)
PROTECT_STATUS = ('certified', 'weak')


def load(path=None):
    cfg = json.loads(json.dumps(DEFAULT))
    path = path or os.path.join(HERE, 'policy.yaml')
    if os.path.exists(path):
        import yaml
        cfg.update(yaml.safe_load(open(path)) or {})
    cfg['targets'] = {int(k): v for k, v in (cfg.get('targets') or {}).items()}
    cfg['extra_n'] = set(int(x) for x in cfg.get('extra_n') or ())
    return cfg


def area_excess(n, s, F):
    return s * s - F * F


def classify(row, F, cfg, protected=False, intake=False):
    """row: mapping with n, s, status, grid.  F: reference side at n (None = unknown -> keep full).
    Returns (tier, reason)."""
    if protected:
        return 'full', 'protected'
    st = row['status']
    if st == 'infeasible':
        return ('drop' if intake else 'stub'), 'infeasible'
    n = row['n']
    t = dict(cfg, **cfg['targets'].get(n, {}))
    if n > t['nmax'] and n not in t['extra_n'] and n not in cfg['targets']:
        return ('drop' if intake else 'stub'), 'n range'
    if (row['grid'] or 0) > 0:
        return ('drop' if intake else 'stub'), 'grid line'
    if F is None:
        return 'full', 'no reference'
    a = area_excess(n, row['s'], F)
    lim = t['screened_full'] if st in ('screened', 'raw') else t['full']
    if a <= lim:
        return 'full', 'near'
    if intake and a > t['drop']:
        return 'drop', 'far'
    return 'stub', 'far'


def references(st):
    """{n: F}: the better of the public frontier and our own best feasible row."""
    F = {n: r['s'] for n, r in st.frontier().items()}
    for n, s in st.db.execute("SELECT n, MIN(s) FROM packing WHERE status != 'infeasible' AND (pen IS NULL OR pen <= 1e-9)"
                              " GROUP BY n"):
        F[n] = min(F.get(n, s), s)
    return F


def protected_ids(st, cfg):
    """Row ids no rule may thin: public sightings, certified / weak, exact sides, near-miss queue, frontier packings,
    protected tags."""
    from .store import public_sql
    P = set()
    q = [f'SELECT DISTINCT g.packing FROM sighting g WHERE {public_sql("g")}',
         f"SELECT id FROM packing WHERE status IN ({','.join('?' * len(PROTECT_STATUS))}) OR s_exact IS NOT NULL",
         'SELECT packing FROM nearmiss', 'SELECT packing FROM frontier WHERE packing IS NOT NULL']
    for i, sql in enumerate(q):
        P.update(r[0] for r in st.db.execute(sql, PROTECT_STATUS if i == 1 else ()))
    for tag in cfg['protect_tags']:
        P.update(r[0] for r in st.db.execute('SELECT id FROM packing WHERE tags LIKE ?', (f'%"{tag}"%',)))
    return P
