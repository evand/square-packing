"""Importers: files / dirs in any format, the jlevy register (witnesses + frontier table), explorer run archives.

10-17: issues and the repos they link are synced by pk/watch.py (parsing in pk/issues.py, every network read in
pk/github.py); their files come in through import_paths here (strict parse for untrusted files).
"""
from __future__ import annotations
import collections, glob, json, os, re, subprocess
import yaml
from .packing import Packing, read, side_val
from .store import Store

from .github import MAX_MEMBER              # one untrusted file's size cap (github.py has the rest)
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
FILE_RE = re.compile(r'\.(txt|yaml|yml|gz|cert|json|sqcert|sqpk)$')
SKIP_RE = re.compile(r'check|output|cert50|summary|SOURCES|README|archive\.jsonl|proposals|state\.json|exact\.json')


def import_paths(st: Store, paths, source, finder=None, ref=None, kind=None, status='raw', ns=None, tol=1e-9,
                 quiet=False, strict=False, collect=None, skip_path=None):
    """Every recognisable packing file under paths -> store.  Returns counts.  strict: untrusted files (pk.packing
    strict parse + validate, size cap MAX_MEMBER).  collect: list to append (path, n, s, id, new) to.  skip_path: regex;
    files whose path (below the given dir) matches are ignored."""
    c = dict(seen=0, new=0, dup=0, infeasible=0, unreadable=0)
    files = []
    for p in paths:
        if os.path.isfile(p):
            files.append(p)
        else:
            files += [f for f in sorted(glob.glob(f'{p}/**/*', recursive=True))
                      if not (skip_path and skip_path.search(os.path.relpath(f, p)))]
    for f in files:
        b = os.path.basename(f)
        if not os.path.isfile(f) or not FILE_RE.search(b) or SKIP_RE.search(b):
            continue
        try:
            pk = read(f, strict=True, max_bytes=MAX_MEMBER) if strict else read(f)
        except Exception as e:
            c['unreadable'] += 1
            if not quiet: print(f'  unreadable {f}: {e!r}'[:160])
            continue
        if pk is None or pk.n == 0 or (ns and pk.n not in ns):
            continue
        c['seen'] += 1
        pid, new = st.add(pk, source, finder, ref or f, kind, status=status, tol=tol, commit=False)
        r = st.db.execute('SELECT status FROM packing WHERE id=?', (pid,)).fetchone()
        if r['status'] == 'infeasible':
            c['infeasible'] += 1
            if not quiet and new: print(f'  infeasible n={pk.n} s={pk.s:.12f} {f}')
        c['new' if new else 'dup'] += 1
        if collect is not None:
            collect.append((f, pk.n, pk.s, pid, new))
    st.db.commit()
    return c


# ---------- jlevy register ----------
def _frontmatter(path):
    return yaml.safe_load(open(path).read().split('---')[1])['packing']


def import_register(st: Store, packing_dir, ref):
    """packing_dir = <squares checkout>/packing (frontier/n-*.md, witnesses/known-best/n-*.yaml)."""
    c = dict(rows=0, witnesses=0, new=0, failed=[])
    for f in sorted(glob.glob(f'{packing_dir}/frontier/n-*.md')):
        fm = _frontmatter(f); n = fm['n']; rep = fm.get('reported_upper_bound') or {}
        v = rep.get('value')
        if v is None:
            continue
        try:                                     # exact rational form when there is one (symbolic forms: keep the decimal)
            v = str(rep['exact_form']) if rep.get('exact_form') and side_val(rep['exact_form']) else v
        except (ValueError, ZeroDivisionError):
            pass
        w = f'{packing_dir}/witnesses/known-best/n-{n:03d}.yaml'
        pid = None
        finder = ', '.join(rep.get('found_by') or []) or None
        if os.path.exists(w):
            try:
                pk = read(w)
                pid, new = st.add(pk, 'register', finder, f'jlevy/squares@{ref}', 'record', commit=False)
                c['witnesses'] += 1; c['new'] += new
            except Exception as e:
                c['failed'].append((n, repr(e)[:80]))
        st.set_frontier(n, side_val(v), str(v), finder, 'register', f'jlevy/squares@{ref}', False, pid)
        c['rows'] += 1
    st.db.commit()
    return c


def fetch_register(dest):
    """Sparse shallow clone of jlevy/squares packing/{frontier,witnesses/known-best} into dest (hooks off, no checkout of
    anything else).  Returns the commit sha."""
    g = ['git', '-c', 'core.hooksPath=/dev/null', '-C', dest]
    if not os.path.isdir(f'{dest}/.git'):
        os.makedirs(dest, exist_ok=True)
        subprocess.run(['git', 'init', '-q', '--template=', dest], check=True)
        subprocess.run(g + ['remote', 'add', 'origin', 'https://github.com/jlevy/squares'], check=True)
        subprocess.run(g + ['sparse-checkout', 'set', '--no-cone', '/packing/frontier/', '/packing/witnesses/known-best/'],
                       check=True)
    subprocess.run(g + ['fetch', '-q', '--depth', '1', '--filter=blob:none', 'origin', 'main'], check=True)
    subprocess.run(g + ['checkout', '-q', '-f', 'FETCH_HEAD'], check=True)
    return subprocess.run(g + ['rev-parse', 'HEAD'], capture_output=True, text=True, check=True).stdout.strip()


# ---------- our explorer runs ----------
EXACT_STATUS = {'certified': 'certified', 'weak': 'weak', 'not-min': 'not-min', 'unresolved': 'unresolved'}


def _prep_run(args):
    """Worker: read one run's archive -> list of (entry, Packing, precomputed) (no db access)."""
    run, below = args
    out = []
    try:
        E = [json.loads(l) for l in open(f'{run}/archive.jsonl')]
    except Exception:
        return run, out
    X = {}
    if os.path.exists(f'{run}/exact.json'):
        try: X = json.load(open(f'{run}/exact.json'))
        except Exception: X = {}
    for e in E:
        if below is not None and e['s'] >= below:
            continue
        path = e['path'] if os.path.isabs(e['path']) else os.path.join(PKG, e['path'])
        if not os.path.exists(path):
            path = os.path.join(run, os.path.basename(e['path']))
        if not os.path.exists(path):
            continue
        try:
            pk = read(path)
        except Exception:
            continue
        if pk is None:
            continue
        x = X.get(str(e['i']), {})
        status = EXACT_STATUS.get(x.get('cls')) or ('screened' if e.get('unpolished') else 'polished')
        out.append((e, pk, dict(hash=pk.canon_hash(), grid=pk.grid_lines(), ntilted=pk.ntilted()), status, x.get('S')))
    return run, out


def import_runs(st: Store, runs, below=None, procs=12, below_k=False):
    """explore.py archives -> store (parallel parse, serial insert), parent links and exact results kept."""
    from concurrent.futures import ProcessPoolExecutor
    tot = dict(runs=0, entries=0, new=0)
    with ProcessPoolExecutor(procs) as ex:
        for run, items in ex.map(_prep_run, [(r, below) for r in runs], chunksize=1):
            src = f'run:{os.path.relpath(os.path.abspath(run), PKG)}'
            ids = {}
            for e, pk, pre, status, S in items:
                if below_k and pk.s >= pk.k:
                    continue
                pid, new = st.add(pk, src, 'ours', str(e['i']), e.get('kind'), status=status,
                                  parent=ids.get(e.get('parent')), check=False, commit=False, pre=pre, s_exact=S,
                                  fp=e.get('cfp'))       # fq contact fingerprint (archives from 10-14)
                ids[e['i']] = pid; tot['new'] += new
            st.db.commit()
            tot['runs'] += 1; tot['entries'] += len(ids)
    return tot


def import_run(st: Store, run, below=None):
    return import_runs(st, [run], below=below, procs=1)


def import_finals(st: Store, n, only_verdict=('certified', 'weak')):
    """Battery / corpus final packings (runs/battery/*/trials.jsonl 'fin', runs/corpus/*/fin_*.npz) never reached the
    `packing` table (review 10-11).  Import those whose side key has an exactsolve verdict in `basin` (default certified /
    weak: the representatives of minima), status = the verdict, source battery:<name> / corpus:<name>/<tag>."""
    import numpy as np
    from .identity import side_key
    B = dict(st.db.execute('SELECT key, cls FROM basin WHERE n=?', (n,)).fetchall())
    root = os.path.join(HERE, '..', 'runs')
    c = collections.Counter()
    for f in sorted(glob.glob(f'{root}/battery/*/trials.jsonl')):
        name = os.path.basename(os.path.dirname(f))
        for i, l in enumerate(open(f)):
            r = json.loads(l)
            if not r.get('fin') or r.get('s') is None or len(r['fin']) != n:
                continue
            cls = B.get(side_key(r['s']))
            if cls not in only_verdict:
                continue
            _, new = st.add(Packing(r['s'], np.array(r['fin'], float)), source=f'battery:{name}', finder='ours',
                            ref=f'{os.path.relpath(f, root)}:{i}', status=cls, commit=False)
            c['battery new' if new else 'battery seen'] += 1
    for f in sorted(glob.glob(f'{root}/corpus/*/fin_*.npz')):
        name = os.path.basename(os.path.dirname(f)); tag = os.path.basename(f)[4:-4]
        z = np.load(f)
        for i in z.files:
            x = z[i]
            if (len(x) - 3) // 3 != n:
                continue
            s = float(x[0]); cls = B.get(side_key(s))
            if cls not in only_verdict:
                continue
            _, new = st.add(Packing(s, x[3:].reshape(-1, 3)), source=f'corpus:{name}/{tag}', finder='ours',
                            ref=f'{os.path.relpath(f, root)}:{i}', status=cls, commit=False)
            c['corpus new' if new else 'corpus seen'] += 1
    st.db.commit()
    return dict(c)
