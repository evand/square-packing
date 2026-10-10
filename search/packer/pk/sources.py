"""Importers: files / dirs in any format, the jlevy register (witnesses + frontier table), explorer run archives, and
pending registration requests (cert repos at a pinned commit, listed in pk/pending.yaml).

Third-party repos are fetched as GitHub tarballs into runs/ext/<owner>-<repo>-<sha7>/ (data only, parsed, never executed;
files and dirs made read-only; runs/ext/PROVENANCE.txt line appended), following ~/math/_untrusted-third-party's rules
(that tree itself stays read-only).
"""
from __future__ import annotations
import glob, json, math, os, re, subprocess, tarfile, tempfile
import yaml
from .packing import Packing, read, side_val
from .store import Store

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EXT = os.path.join(PKG, 'runs', 'ext')
FILE_RE = re.compile(r'\.(txt|yaml|yml|gz|cert|json)$')
SKIP_RE = re.compile(r'check|output|cert50|summary|SOURCES|README|archive\.jsonl|proposals|state\.json|exact\.json')


def import_paths(st: Store, paths, source, finder=None, ref=None, kind=None, status='raw', ns=None, tol=1e-9,
                 quiet=False):
    """Every recognisable packing file under paths -> store.  Returns counts."""
    c = dict(seen=0, new=0, dup=0, infeasible=0, unreadable=0)
    files = []
    for p in paths:
        files += [p] if os.path.isfile(p) else sorted(glob.glob(f'{p}/**/*', recursive=True))
    for f in files:
        b = os.path.basename(f)
        if not os.path.isfile(f) or not FILE_RE.search(b) or SKIP_RE.search(b):
            continue
        try:
            pk = read(f)
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


# ---------- third-party repos (pending registrations) ----------
def fetch_repo(owner_repo, commit):
    """GitHub tarball of owner/repo@commit -> runs/ext/<owner>-<repo>-<sha7>/ (cached).  Returns the dir."""
    owner, repo = owner_repo.split('/')
    if len(commit) != 40:
        commit = subprocess.run(['gh', 'api', f'repos/{owner_repo}/commits/{commit}', '--jq', '.sha'],
                                capture_output=True, text=True, check=True).stdout.strip()
    os.makedirs(EXT, exist_ok=True)
    d = f'{EXT}/{owner}-{repo}-{commit[:7]}'
    if os.path.isdir(d):
        return d
    with tempfile.TemporaryDirectory(dir=EXT) as tmp:
        tgz = f'{tmp}/repo.tgz'
        with open(tgz, 'wb') as f:
            subprocess.run(['gh', 'api', f'repos/{owner_repo}/tarball/{commit}'], stdout=f, check=True)
        with tarfile.open(tgz) as t:
            t.extractall(tmp + '/x', filter='data')
        top = glob.glob(f'{tmp}/x/*')[0]
        os.rename(top, d)
    for root, dirs, files in os.walk(d):
        for fn in files:
            p = os.path.join(root, fn)
            if not os.path.islink(p):
                os.chmod(p, 0o444)
    for root, dirs, files in os.walk(d, topdown=False):
        os.chmod(root, 0o555)
    with open(f'{EXT}/PROVENANCE.txt', 'a') as f:
        f.write(f'\n{owner}-{repo}-{commit[:7]} — tarball {owner_repo}@{commit} fetched by pk/sources.py '
                f'(data only, read-only, nothing executed)\n')
    return d


def issue_claims(number, repo='jlevy/squares'):
    """(author, {n: claimed side}) from the markdown tables of an issue body (first decimal column after n)."""
    r = json.loads(subprocess.run(['gh', 'issue', 'view', str(number), '-R', repo, '--json', 'body,author,state'],
                                  capture_output=True, text=True, check=True).stdout)
    claims = {}
    for line in r['body'].split('\n'):
        cells = [x.strip(' *`$') for x in line.strip().strip('|').split('|')]
        if len(cells) < 2 or not re.fullmatch(r'\d{2,4}', cells[0]):
            continue
        m = re.match(r'(\d+\.\d{6,})', cells[1])
        if m and float(m.group(1)) > 4:
            n = int(cells[0]); claims[n] = min(claims.get(n, 9e9), float(m.group(1)))
    return r['author']['login'], r['state'], claims


def import_pending(st: Store, cfg_path=os.path.join(HERE, 'pending.yaml'), only=None):
    """pending.yaml: list of {issue, finder, repo, commit, paths: [globs], ns: optional}.  For each: claims from the issue
    table -> frontier rows (pending=1, source 'pending:#<issue>'); packings from the repo -> store; report each claim's
    best stored packing (claim not matched by a feasible file is flagged)."""
    cfg = yaml.safe_load(open(cfg_path))
    report = []
    for e in cfg:
        if only and e['issue'] not in only:
            continue
        author, state, claims = issue_claims(e['issue'])
        src = f"pending:#{e['issue']}"
        ref = f"{e['repo']}@{e['commit'][:7]}"
        d = fetch_repo(e['repo'], e['commit'])
        paths = [p for g in e.get('paths', ['.']) for p in glob.glob(os.path.join(d, g), recursive=True)]
        c = import_paths(st, paths, src, e.get('finder', author), ref, 'pending', quiet=True)
        for n, s in sorted(claims.items()):
            row = st.query(n=n, source=src, limit=1)
            pid = row[0]['id'] if row else None
            ok = bool(row) and row[0]['s'] <= s + 1e-9
            st.set_frontier(n, s, f'{s}', e.get('finder', author), src, ref, state == 'OPEN', pid)
            report.append((e['issue'], n, s, row[0]['s'] if row else None, ok))
        st.db.commit()
        print(f"#{e['issue']} {author} ({state}): {len(claims)} claims; files {c}")
    bad = [r for r in report if not r[4]]
    for r in bad:
        print(f'  claim not matched by a feasible stored packing: #{r[0]} n={r[1]} claim {r[2]} best file {r[3]}')
    return report


# ---------- our explorer runs ----------
EXACT_STATUS = {'certified': 'certified', 'not-min': 'not-min', 'unresolved': 'unresolved'}


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
                                  parent=ids.get(e.get('parent')), check=False, commit=False, pre=pre, s_exact=S)
                ids[e['i']] = pid; tot['new'] += new
            st.db.commit()
            tot['runs'] += 1; tot['entries'] += len(ids)
    return tot


def import_run(st: Store, run, below=None):
    return import_runs(st, [run], below=below, procs=1)
