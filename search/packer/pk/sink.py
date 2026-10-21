"""One sink (10-20; Evan: "If they're distinct packings, track them"): every evaluated result reaches the store, under
one intake policy (pk/policy.py + pk/policy.yaml), with detection on import.

What it reads (adapters; every run dir under runs/, ext/ excluded):
  explore    <run>/archive.jsonl (explore.py, dgen / hborder / chain sub-runs, time-travel worlds ...): every entry,
             parent links kept, exactsolve verdicts from exact.json.  Incremental: the state remembers how many archive
             lines were imported, so a live run is imported in slices (only lines older than SETTLE are taken).
  outputs    compose2 near.jsonl, results.jsonl / polish.jsonl rows (compose2, fresh, harden): nearmiss.read_outputs
             (thr = inf: every row).
Dirs with a manifest but no adapter are listed by `pk.py sink --report` (by manifest name), never guessed at.

Intake: policy.classify(intake=True): 'drop' rows are counted, not stored.  The store dedupes by canon hash (a repeat
is a sighting on the existing row) and follows pk/ident's dedup table.

Detection (on the rows this pass created or sighted):
  hit        s < public frontier(n) - 1e-9 and no public sighting -> tag 'candidate-hit', line in runs/sink/hits.jsonl
             (candidate.py is not run here: it is the checker, launched by hand / the finisher)
  near miss  nearmiss.scan(ids=...) -> the queue (pk.py finish drains it)

  pk.py sink [--dry-run] [--runs DIR ...] [--procs 12] [--report]
State: runs/sink/state.json {run: lines imported}; log runs/sink/log.jsonl (one line per pass).
"""
from __future__ import annotations
import glob, json, os, time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)
RUNS = os.path.join(PK, 'runs')
SDIR = os.path.join(RUNS, 'sink')
SETTLE = 600.0          # s: archive lines / output files younger than this wait for the next pass
HIT = 1e-9


def _rel(d):
    return os.path.relpath(os.path.abspath(d), PK)


def load_state():
    try:
        return json.load(open(os.path.join(SDIR, 'state.json')))
    except FileNotFoundError:
        return dict(explore={}, outputs={})


def save_state(S):
    os.makedirs(SDIR, exist_ok=True)
    tmp = os.path.join(SDIR, 'state.json.tmp')
    json.dump(S, open(tmp, 'w'), indent=0, sort_keys=True); os.replace(tmp, os.path.join(SDIR, 'state.json'))


def explore_runs(root=RUNS):
    return sorted(os.path.dirname(a) for a in glob.glob(os.path.join(root, '**', 'archive.jsonl'), recursive=True)
                  if '/ext/' not in a)


def output_dirs(root=RUNS):
    out = set()
    for pat in ('near.jsonl', 'results.jsonl', 'polish.jsonl'):
        for f in glob.glob(os.path.join(root, '**', pat), recursive=True):
            if '/ext/' not in f and not os.path.exists(os.path.join(os.path.dirname(f), 'archive.jsonl')):
                out.add(os.path.dirname(f))
    return sorted(out)


def _prep_explore(args):
    """worker: archive lines [start, ...) of one run, settled ones only -> (run, next_start, items)"""
    run, start, now = args
    from .sources import _prep_run
    try:
        if now - os.path.getmtime(os.path.join(run, 'archive.jsonl')) < SETTLE:
            lines = open(os.path.join(run, 'archive.jsonl')).read().splitlines()
            stop = len(lines) - 1                     # the newest line may be half written; take it next pass
        else:
            stop = sum(1 for _ in open(os.path.join(run, 'archive.jsonl')))
    except OSError:
        return run, start, []
    if stop <= start:
        return run, start, []
    _, items = _prep_run((run, None, start, stop))
    return run, stop, items


def _parent_lookup(st, src, i):
    r = st.db.execute('SELECT packing FROM sighting WHERE source=? AND ref=?', (src, str(i))).fetchone()
    return r[0] if r else None


def run(st=None, runs=None, dry=False, procs=12, log=print, outputs=True):
    from . import policy as P
    from .store import Store, public_sql
    st = st or Store()
    cfg = P.load()
    F = P.references(st)
    pubF = {n: r['s'] for n, r in st.frontier().items()}
    S = load_state()
    now = time.time()
    runs = explore_runs() if runs is None else [r for r in runs if os.path.exists(os.path.join(r, 'archive.jsonl'))]
    todo = [(r, S['explore'].get(_rel(r), 0), now) for r in runs]
    tot = dict(runs=0, entries=0, new=0, seen=0, drop=0)
    touched = set()
    with ProcessPoolExecutor(procs) as ex:
        for r, stop, items in ex.map(_prep_explore, todo, chunksize=4):
            if not items:
                if stop != S['explore'].get(_rel(r), 0):
                    S['explore'][_rel(r)] = stop
                continue
            src = f'run:{_rel(r)}'
            ids = {}
            for e, pk, pre, status, Sx in items:
                row = dict(n=pk.n, s=float(pk.s), status=status, grid=pre['grid'])
                tier, _ = P.classify(row, F.get(pk.n), cfg, intake=True)
                if tier == 'drop':
                    tot['drop'] += 1
                    continue
                if dry:
                    tot['entries'] += 1
                    continue
                par = e.get('parent')
                parent = ids.get(par) if par in ids else (_parent_lookup(st, src, par) if par is not None else None)
                pid, new = st.add(pk, src, 'ours', str(e['i']), e.get('kind'), status=status, parent=parent,
                                  check=False, commit=False, pre=pre, s_exact=Sx, fp=e.get('cfp'))
                ids[e['i']] = pid; touched.add(pid)
                tot['entries'] += 1; tot['new' if new else 'seen'] += 1
            if not dry:
                st.db.commit()
                S['explore'][_rel(r)] = stop
                save_state(S)
            tot['runs'] += 1
    if outputs:
        from . import nearmiss as NM
        for d in output_dirs() if runs is None or outputs == 'all' else []:
            key = _rel(d)
            mt = max(os.path.getmtime(os.path.join(d, f)) for f in ('near.jsonl', 'results.jsonl', 'polish.jsonl')
                     if os.path.exists(os.path.join(d, f)))
            if now - mt < SETTLE or S['outputs'].get(key) == mt:
                continue
            try:
                R = NM.read_outputs(d, thr=float('inf'), nmax=10 ** 6)
            except Exception as e:
                log(f'sink: {key}: {e!r}'); continue
            for kind, n, s, p in R:
                if dry:
                    tot['entries'] += 1; continue
                src = ('harden:' if kind == 'harden' else 'compose2:') + key
                pid, new = st.add(p, src, finder='ours', ref=key, kind=kind,
                                  status='polished' if kind != 'near' else 'raw', tol=1e-7 if kind == 'near' else 1e-9,
                                  commit=False)
                touched.add(pid); tot['entries'] += 1; tot['new' if new else 'seen'] += 1
            if not dry:
                st.db.commit(); S['outputs'][key] = mt; save_state(S)
    hits = []
    if not dry and touched:
        ids = sorted(touched)
        for k in range(0, len(ids), 900):
            ch = ids[k:k + 900]
            for r in st.db.execute(f"SELECT id, n, s FROM packing WHERE id IN ({','.join('?' * len(ch))}) "
                                   f"AND status != 'infeasible' AND NOT EXISTS (SELECT 1 FROM sighting g WHERE "
                                   f"g.packing=packing.id AND {public_sql('g')})", ch):
                f = pubF.get(r['n'])
                if f is not None and r['s'] < f - HIT:
                    hits.append(dict(id=r['id'], n=r['n'], s=r['s'], frontier=f, gap=r['s'] - f))
        if hits:
            os.makedirs(SDIR, exist_ok=True)
            with open(os.path.join(SDIR, 'hits.jsonl'), 'a') as fh:
                for h in hits:
                    st.tag(h['id'], 'candidate-hit')
                    fh.write(json.dumps(dict(h, t=time.strftime('%Y-%m-%dT%H:%M:%S'))) + '\n')
        from . import nearmiss as NM
        fr = lambda n: pubF.get(n)
        items = (NM.scan(st.db, fr, ids=ids) if len(ids) <= 20000 else
                 NM.scan(st.db, fr, ns={r[0] for k in range(0, len(ids), 900) for r in st.db.execute(
                     f"SELECT DISTINCT n FROM packing WHERE id IN ({','.join('?' * len(ids[k:k + 900]))})", ids[k:k + 900])}))
        tot['nearmiss'] = NM.enqueue(st, items, note='sink')
    tot['hits'] = len(hits)
    if not dry:
        os.makedirs(SDIR, exist_ok=True)
        with open(os.path.join(SDIR, 'log.jsonl'), 'a') as fh:
            fh.write(json.dumps(dict(t=time.strftime('%Y-%m-%dT%H:%M:%S'), **tot)) + '\n')
    for h in hits:
        log(f"sink: HIT n={h['n']} s={h['s']:.12f} ({h['gap']:+.3e} vs frontier) id {h['id']}")
    return tot


def report(root=RUNS):
    """manifest names of run dirs no adapter reads (counts + an example)"""
    import collections
    covered = set(explore_runs(root)) | set(output_dirs(root))
    pref = set()
    for c in covered:                                  # a launcher whose sub-runs are covered (chain, lcut, dgen) is too
        d = os.path.dirname(c)
        while d.startswith(root) and d not in pref:
            pref.add(d); d = os.path.dirname(d)
    C = collections.defaultdict(list)
    for m in glob.glob(os.path.join(root, '**', 'manifest.json'), recursive=True):
        d = os.path.dirname(m)
        if '/ext/' in m or d in covered or d in pref:
            continue
        try:
            name = json.load(open(m)).get('name')
        except Exception:
            name = '?'
        C[name].append(_rel(d))
    return sorted(((len(v), k, v[0]) for k, v in C.items()), reverse=True)
