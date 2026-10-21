"""Explorer archive coordinate files (runs/**/b#####.txt, one per archive entry) deduplicated against the store
(10-20; Evan: "dedupe the per-result coords against the store").

A file is removed only when the store row its import created (sighting source run:<dir>, ref = entry index)
regenerates it byte for byte (mcmin.write_deg's format: `n s!r` then `x!r y!r a!r` rows).  Each removal is written first to
<run>/coords_in_store.jsonl {i, file, packing, sha1}; `ensure(path)` (called by pk.packing.read and mcmin.load_deg on a
missing file) and `restore(run)` write the file back from the store and check the sha1.  Lossless as long as the store
row keeps its coordinates: pk/prune must treat these rows as backing files (`backing_ids`).

Only idle runs: archive.jsonl untouched for min_age_h (a live explorer appends to it and reads its own files)."""
from __future__ import annotations
import hashlib, json, os, sqlite3, time

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = 'coords_in_store.jsonl'


def _store_path():
    return os.environ.get('PK_STORE', os.path.join(PKG, 'runs', 'store.sqlite'))


def _ro(path=None):
    db = sqlite3.connect(f'file:{path or _store_path()}?mode=ro', uri=True, timeout=60)
    return db


def text_of(s, coords):
    """mcmin.write_deg's bytes for (s, coords)."""
    out = [f'{len(coords)} {s!r}\n']
    for x, y, a in coords:
        out.append(f'{x!r} {y!r} {a % 90!r}\n')
    return ''.join(out).encode()


def _row(db, pid):
    import numpy as np
    r = db.execute('SELECT s, coords FROM packing WHERE id=?', (pid,)).fetchone()
    if r is None or not r[1]:
        return None
    return float(r[0]), [tuple(map(float, q)) for q in np.frombuffer(r[1], '<f8').reshape(-1, 3)]


def _entry_file(run, e):
    p = e['path'] if os.path.isabs(e['path']) else os.path.join(PKG, e['path'])
    if os.path.dirname(os.path.abspath(p)) != os.path.abspath(run):
        return None                                   # only the run's own files (starts copied elsewhere are not ours)
    return p


def dedupe_run(run, apply=False, min_age_h=12.0, store=None):
    """-> dict(run, entries, removed, bytes, kept_not_in_store, kept_differs, skipped)."""
    out = dict(run=run, entries=0, removed=0, bytes=0, not_in_store=0, differs=0, missing=0, skipped=None)
    arch = os.path.join(run, 'archive.jsonl')
    try:
        if time.time() - os.path.getmtime(arch) < min_age_h * 3600:
            out['skipped'] = 'active'; return out
        E = [json.loads(l) for l in open(arch)]
    except Exception as ex:
        out['skipped'] = f'unreadable archive ({type(ex).__name__})'; return out
    db = _ro(store)
    src = 'run:' + os.path.relpath(os.path.abspath(run), PKG)
    ids = dict(db.execute('SELECT ref, packing FROM sighting WHERE source=?', (src,)).fetchall())
    man = open(os.path.join(run, MANIFEST), 'a') if apply else None
    try:
        for e in E:
            out['entries'] += 1
            f = _entry_file(run, e)
            if f is None or not os.path.exists(f):
                out['missing'] += 1; continue
            pid = ids.get(str(e['i']))
            if pid is None:
                out['not_in_store'] += 1; continue
            r = _row(db, pid)
            b = open(f, 'rb').read()
            if r is None or text_of(*r) != b:
                out['differs'] += 1; continue
            out['removed'] += 1; out['bytes'] += os.stat(f).st_blocks * 512      # disk use, not length
            if apply:
                man.write(json.dumps(dict(i=e['i'], file=os.path.basename(f), packing=pid,
                                          sha1=hashlib.sha1(b).hexdigest())) + '\n')
                man.flush(); os.fsync(man.fileno())
                os.remove(f)
    finally:
        if man:
            man.close()
        db.close()
    return out


def _manifest(run):
    p = os.path.join(run, MANIFEST)
    if not os.path.exists(p):
        return {}
    return {d['file']: d for d in map(json.loads, open(p))}


def ensure(path):
    """If `path` is an archive coordinate file removed by dedupe_run, write it back from the store.  Returns True if
    the file exists afterwards."""
    if os.path.exists(path):
        return True
    run, fn = os.path.split(os.path.abspath(path))
    d = _manifest(run).get(fn)
    if d is None:
        return False
    db = _ro()
    try:
        r = _row(db, d['packing'])
        rep = None
        if r is None:                                   # the row was merged by pk/ident's dedupe: its representative
            from .ident import resolve_old
            rep = resolve_old(db, d['packing'])
            r = _row(db, rep) if rep is not None else None
    finally:
        db.close()
    if r is None:
        raise FileNotFoundError(f'{path}: removed by archive dedupe, but store row {d["packing"]} has no coordinates')
    b = text_of(*r)
    if rep is not None:                                 # same minimum, not the same bytes: say so next to the file
        with open(os.path.join(run, 'restored_from_dedup.txt'), 'a') as f:
            f.write(f'{fn}: store row {d["packing"]} was merged into {rep} (pk/ident); written from {rep}\n')
    elif hashlib.sha1(b).hexdigest() != d['sha1']:
        raise ValueError(f'{path}: store row {d["packing"]} no longer regenerates the removed file')
    tmp = f'{path}.restore{os.getpid()}'
    with open(tmp, 'wb') as f:
        f.write(b)
    os.replace(tmp, path)
    return True


def restore(run):
    n = 0
    for fn in _manifest(run):
        p = os.path.join(run, fn)
        if not os.path.exists(p):
            ensure(p); n += 1
    return n


def backing_ids(root=None):  # noqa: E302
    """Store rows whose coordinates back a removed archive file."""
    import glob
    root = root or os.path.join(PKG, 'runs')
    S = set()
    for m in glob.glob(f'{root}/**/{MANIFEST}', recursive=True):
        S.update(json.loads(l)['packing'] for l in open(m))
    return S


def runs_under(root):
    import glob
    return sorted(os.path.dirname(p) for p in glob.glob(f'{root}/**/archive.jsonl', recursive=True)
                  if '/ext/' not in p)


def dedupe(roots, apply=False, min_age_h=12.0, procs=16, log=print):
    from concurrent.futures import ProcessPoolExecutor
    runs = [r for root in roots for r in (runs_under(root) if not os.path.exists(os.path.join(root, 'archive.jsonl'))
                                          else [root])]
    tot = dict(runs=0, entries=0, removed=0, bytes=0, not_in_store=0, differs=0, missing=0, skipped=0)
    with ProcessPoolExecutor(procs) as ex:
        for r in ex.map(dedupe_run, runs, [apply] * len(runs), [min_age_h] * len(runs), chunksize=4):
            tot['runs'] += 1
            if r['skipped']:
                tot['skipped'] += 1; continue
            for k in ('entries', 'removed', 'bytes', 'not_in_store', 'differs', 'missing'):
                tot[k] += r[k]
    log(f"{'removed' if apply else 'would remove'} {tot['removed']:,} files ({tot['bytes'] / 1e9:.2f} GB) of "
        f"{tot['entries']:,} entries in {tot['runs'] - tot['skipped']:,} runs ({tot['skipped']} skipped: active / unreadable); "
        f"kept: not in store {tot['not_in_store']:,}, store copy differs {tot['differs']:,}; already gone {tot['missing']:,}")
    return tot
