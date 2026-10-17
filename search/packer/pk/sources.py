"""Importers: files / dirs in any format, the jlevy register (witnesses + frontier table), explorer run archives, and
pending registration requests (cert repos at a pinned commit, listed in pk/pending.yaml).

Third-party repos are fetched as GitHub tarballs into runs/ext/<owner>-<repo>-<sha7>/ (data only, parsed, never executed;
files and dirs made read-only; runs/ext/PROVENANCE.txt line appended), following ~/math/_untrusted-third-party's rules
(that tree itself stays read-only).
"""
from __future__ import annotations
import collections, glob, json, math, os, re, subprocess, tarfile, tempfile
import yaml
from .packing import Packing, read, side_val
from .store import Store

# Caps for untrusted downloads (pk/README.md, Security)
MAX_TARBALL = 200 << 20          # compressed tarball bytes
MAX_MEMBER = 20 << 20            # one extracted file (bigger members are skipped, not extracted)
MAX_MEMBERS = 20000              # tar entries
MAX_UNPACKED = 1 << 30           # total extracted bytes
MAX_ATTACH = 20 << 20            # one issue attachment
REPO_NAME = re.compile(r'[A-Za-z0-9](?:[A-Za-z0-9-]{0,38})/[A-Za-z0-9._-]{1,100}')
SHA = re.compile(r'[0-9a-f]{7,40}')
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EXT = os.path.join(PKG, 'runs', 'ext')
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


# ---------- third-party repos (pending registrations) ----------
def _gh(args, timeout=120):
    """gh api call (authenticated), args a list (no shell).  Returns stdout text; raises on failure."""
    return subprocess.run(['gh', 'api'] + args, capture_output=True, text=True, check=True, timeout=timeout).stdout


def resolve_sha(owner_repo, ref):
    """Full sha of a commit / branch on owner/repo (GitHub API), or None if it does not exist there."""
    if not REPO_NAME.fullmatch(owner_repo) or not re.fullmatch(r'[A-Za-z0-9._/-]{1,200}', ref) or '..' in ref:
        raise ValueError(f'bad repo / ref {owner_repo!r} {ref!r}')
    try:
        sha = _gh([f'repos/{owner_repo}/commits/{ref}', '--jq', '.sha']).strip()
    except subprocess.CalledProcessError:
        return None
    return sha if re.fullmatch(r'[0-9a-f]{40}', sha) else None


def _download_capped(cmd, dest, cap, timeout=600):
    """Run cmd (list) with stdout -> dest, killing it once more than cap bytes arrive.  Raises on failure / cap."""
    n = 0
    with open(dest, 'wb') as f, subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as pr:
        try:
            while True:
                b = pr.stdout.read(1 << 16)
                if not b:
                    break
                n += len(b)
                if n > cap:
                    raise ValueError(f'download over cap ({cap} bytes): {cmd[:3]}')
                f.write(b)
            if pr.wait(timeout=timeout) != 0:
                raise subprocess.CalledProcessError(pr.returncode, cmd)
        finally:
            if pr.poll() is None:
                pr.kill()
    return n


def _extract_capped(tgz, dest):
    """tarfile filter='data' (no links out, no devices, no absolute paths, no setuid), plus caps: member count, per-file
    size (bigger files skipped), total size.  Returns (extracted, skipped)."""
    skipped = []
    with tarfile.open(tgz) as t:
        members = []
        for m in t:                                         # streams headers; stops at the cap
            members.append(m)
            if len(members) > MAX_MEMBERS:
                raise ValueError(f'tarball has more than {MAX_MEMBERS} entries')
        total = sum(m.size for m in members if m.isfile() and m.size <= MAX_MEMBER)
        if total > MAX_UNPACKED:
            raise ValueError(f'tarball unpacks to {total} bytes > {MAX_UNPACKED}')

        def flt(m, path):
            m = tarfile.data_filter(m, path)
            if m is not None and m.isfile() and m.size > MAX_MEMBER:
                skipped.append(m.name)
                return None
            return m
        t.extractall(dest, members=members, filter=flt)
    return len(members) - len(skipped), skipped


def fetch_repo(owner_repo, commit):
    """GitHub tarball of owner/repo@commit -> runs/ext/<owner>-<repo>-<sha7>/ (cached).  Returns the dir.
    Untrusted data: size-capped download and extraction (MAX_*), tarfile filter 'data', everything made read-only,
    nothing in it is ever executed."""
    if not REPO_NAME.fullmatch(owner_repo) or not re.fullmatch(r'[0-9a-fA-F]{7,40}', commit):
        raise ValueError(f'bad repo / commit {owner_repo!r} {commit!r}')
    owner, repo = owner_repo.split('/')
    if len(commit) != 40:
        commit = resolve_sha(owner_repo, commit)
        if commit is None:
            raise ValueError(f'{owner_repo}: commit not found')
    commit = commit.lower()
    os.makedirs(EXT, exist_ok=True)
    d = f'{EXT}/{owner}-{repo}-{commit[:7]}'
    if os.path.isdir(d):
        return d
    with tempfile.TemporaryDirectory(dir=EXT) as tmp:
        tgz = f'{tmp}/repo.tgz'
        nbytes = _download_capped(['gh', 'api', f'repos/{owner_repo}/tarball/{commit}'], tgz, MAX_TARBALL)
        nfiles, skipped = _extract_capped(tgz, tmp + '/x')
        top = glob.glob(f'{tmp}/x/*')
        if len(top) != 1 or not os.path.isdir(top[0]) or os.path.islink(top[0]):
            raise ValueError(f'{owner_repo}@{commit[:7]}: unexpected tarball layout')
        os.rename(top[0], d)
    _read_only(d)
    with open(f'{EXT}/PROVENANCE.txt', 'a') as f:
        f.write(f'\n{owner}-{repo}-{commit[:7]} — tarball {owner_repo}@{commit} fetched by pk/sources.py '
                f'(data only, read-only, nothing executed; {nbytes} bytes, {nfiles} entries'
                f'{f", skipped over {MAX_MEMBER >> 20} MB: {skipped[:5]}" if skipped else ""})\n')
    return d


def _read_only(d):
    for root, dirs, files in os.walk(d):
        for fn in files:
            p = os.path.join(root, fn)
            if not os.path.islink(p):
                os.chmod(p, 0o444)
    for root, dirs, files in os.walk(d, topdown=False):
        os.chmod(root, 0o555)


# Issue attachments: only https://github.com/user-attachments/files/<id>/<name>; the redirect is followed by hand, once,
# and only to GitHub's own object hosts.
ATTACH_RE = re.compile(r'https://github\.com/user-attachments/files/(\d{1,12})/([A-Za-z0-9._-]{1,120})')
ATTACH_HOSTS = re.compile(r'(objects\.githubusercontent\.com|github-production-[a-z0-9-]+\.s3\.amazonaws\.com)')
ATTACH_DIR = os.path.join(EXT, 'attachments')


def fetch_attachment(url, issue):
    """One issue attachment -> runs/ext/attachments/issue-<n>/<id>-<name> (cached, read-only).  Returns (path, sha256)."""
    import hashlib, urllib.parse
    m = ATTACH_RE.fullmatch(url)
    if not m:
        raise ValueError(f'not an attachment url: {url[:80]}')
    d = f'{ATTACH_DIR}/issue-{int(issue)}'
    dest = f'{d}/{m.group(1)}-{m.group(2)}'
    if not os.path.exists(dest):
        os.makedirs(d, exist_ok=True)
        h = subprocess.run(['curl', '-sS', '--proto', '=https', '--max-redirs', '0', '-o', '/dev/null', '-D', '-', url],
                           capture_output=True, text=True, timeout=60).stdout
        loc = re.search(r'^location:\s*(\S+)', h, re.M | re.I)
        if not loc:
            raise ValueError(f'{url}: no redirect')
        u = urllib.parse.urlsplit(loc.group(1))
        if u.scheme != 'https' or not ATTACH_HOSTS.fullmatch(u.hostname or ''):
            raise ValueError(f'{url}: redirect to an unexpected host {u.hostname!r}')
        tmp = dest + '.part'
        try:
            _download_capped(['curl', '-sS', '-f', '--proto', '=https', '--max-redirs', '0',
                              '--max-filesize', str(MAX_ATTACH), loc.group(1)], tmp, MAX_ATTACH, timeout=120)
            os.replace(tmp, dest)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)
        os.chmod(dest, 0o444)
        with open(f'{EXT}/PROVENANCE.txt', 'a') as f:
            f.write(f'\nattachments/issue-{int(issue)}/{os.path.basename(dest)} — {url} (jlevy/squares#{issue}), fetched by '
                    f'pk/sources.py (data only, read-only, nothing executed)\n')
    return dest, hashlib.sha256(open(dest, 'rb').read()).hexdigest()


def body_claims(body):
    """{n: claimed side} from the markdown tables of an issue body.  10-17: the smallest cell in the row that can be a
    side for n (sqrt(n) <= v <= sqrt(n) + 3), not the first decimal column: #668's table reads "n | website-listed
    side | certified upper bound | reduction", so the first column was the OLD record.  Tables with or without the
    outer pipes."""
    claims = {}
    for line in (body or '').split('\n'):
        cells = [x.strip(' *`$') for x in line.strip().strip('|').split('|')]
        m0 = re.fullmatch(r'\[?(\d{2,4})\]?(?:\([^)]*\))?', cells[0]) if len(cells) >= 2 else None   # 237 or [237](url)
        if not m0:
            continue
        n = int(m0.group(1))
        v = [float(m.group(1)) for c in cells[1:] for m in [re.match(r'(\d+\.\d{6,})', c)] if m]
        v = [x for x in v if math.sqrt(n) - 1e-9 <= x <= math.sqrt(n) + 3]
        if v:
            claims[n] = min(claims.get(n, 9e9), min(v))
    return claims


TITLE_CLAIM = re.compile(r's\(\s*(\d{1,5})\s*\)\s*(?:≤|<=|\\le)\s*(\d+\.\d{6,})')


def issue_text(number, repo='jlevy/squares'):
    """(author, state, title, body, [the author's comment bodies]) of an issue (10-17)"""
    r = json.loads(subprocess.run(['gh', 'issue', 'view', str(int(number)), '-R', repo, '--json',
                                   'title,body,author,state,comments'],
                                  capture_output=True, text=True, check=True, timeout=120).stdout)
    who = r['author']['login']
    # claims: the author's own comments only (the maintainer's verdict tables restate other issues' values, #503);
    # links and attachments are read from every comment by watch.py
    return who, r['state'], r.get('title') or '', r.get('body') or '', \
        [c.get('body') or '' for c in r.get('comments') or [] if (c.get('author') or {}).get('login') == who]


def text_claims(title, body, comments=()):
    """{n: smallest claimed side}: markdown tables in the body and every comment, plus `s(n) <= x` in the title"""
    claims = {}
    for txt in [body, *comments]:
        for n, s in body_claims(txt).items():
            claims[n] = min(claims.get(n, 9e9), s)
    for n, s in TITLE_CLAIM.findall(title or ''):
        claims[int(n)] = min(claims.get(int(n), 9e9), float(s))
    return claims


def issue_claims(number, repo='jlevy/squares'):
    """(author, state, {n: claimed side}) from the issue's tables (body and comments) and title"""
    author, state, title, body, comments = issue_text(number, repo)
    return author, state, text_claims(title, body, comments)


CLAIM_TOL = 1e-6        # a file matches a claim when claim - CLAIM_TOL <= s <= claim + 1e-9 (claims are rounded up)


def reconcile(st, src, claims):
    """10-17: sort an issue's packings (sightings with source `src`) into 'claim' (its n is claimed and its side matches
    the claim) and 'context' (anything else the linked files held: other counts, older packings, parents).  Every
    packing stays in the store and counts for record checks; the kind only says what the issue claimed.
    Returns dict(matched={n: pid}, unmatched_claims={n: s}, context=k, below_claim=[(n, s, claim, pid)])."""
    rows = st.db.execute('SELECT p.id, p.n, p.s, p.status FROM packing p JOIN sighting g ON g.packing = p.id '
                         'WHERE g.source = ?', (src,)).fetchall()
    out = dict(matched={}, unmatched_claims={}, context=0, below_claim=[])
    for r in rows:
        c = claims.get(r['n'])
        ok = r['status'] != 'infeasible' and c is not None and c - CLAIM_TOL <= r['s'] <= c + 1e-9
        st.db.execute('UPDATE sighting SET kind=? WHERE packing=? AND source=?', ('claim' if ok else 'context', r['id'], src))
        if ok:
            out['matched'][r['n']] = r['id']
        else:
            out['context'] += 1
            if c is not None and r['status'] != 'infeasible' and r['s'] < c - CLAIM_TOL:
                out['below_claim'].append((r['n'], r['s'], c, r['id']))
    out['unmatched_claims'] = {n: s for n, s in claims.items() if n not in out['matched']}
    return out


def import_pending(st: Store, cfg_path=os.path.join(HERE, 'pending.yaml'), only=None):
    """pending.yaml: list of {issue, finder, repo, commit, paths: [globs], attachments: [urls], claims}; repo: null =
    no repo (claims and attachments only).  An issue may have several entries (10-17: one per linked commit, plus
    supplements).  Per ISSUE: every entry's files -> store (source 'pending:#<issue>'); claims (issue body, comments,
    title, and the entries' own) -> frontier rows (pending=1); then reconcile once: each packing 'claim' or 'context',
    and every claim without a matching file is reported."""
    cfg = yaml.safe_load(open(cfg_path))
    by_issue = collections.OrderedDict()
    for e in cfg:
        if not only or e['issue'] in only:
            by_issue.setdefault(e['issue'], []).append(e)
    report = []
    for issue, E in by_issue.items():
        author, state, claims = issue_claims(issue)
        src = f"pending:#{issue}"
        finder = E[0].get('finder', author)
        refs_, counts = [], []
        for e in E:
            for n_, s_ in (e.get('claims') or {}).items():  # claims stated outside a table (10-17: min, not overwrite)
                claims[int(n_)] = min(claims.get(int(n_), 9e9), float(s_))
            if e.get('repo'):
                ref = f"{e['repo']}@{e['commit'][:7]}"
                d = fetch_repo(e['repo'], e['commit'])
                paths = [p for g in e.get('paths', ['.']) for p in glob.glob(os.path.join(d, g), recursive=True)
                         if os.path.realpath(p).startswith(os.path.realpath(d))]
                c = import_paths(st, paths, src, e.get('finder', author), ref, 'pending', quiet=True, strict=True)
                refs_.append(ref); counts.append(c)
            for url in e.get('attachments') or ():     # files attached to the issue (github user-attachments only)
                try:
                    f, sha = fetch_attachment(url, issue)
                except Exception as ex:
                    print(f"  #{issue} attachment not fetched: {url} ({ex!r})"[:200]); continue
                aref = f"#{issue} attachment {os.path.basename(f)}"
                ca = import_paths(st, [f], src, e.get('finder', author), aref, 'pending', quiet=True, strict=True)
                refs_.append(aref); counts.append(ca)
                print(f"  #{issue} attachment {os.path.basename(f)} sha256 {sha[:16]}: {ca}")
        seen = sum(c['seen'] for c in counts)
        R = reconcile(st, src, claims)
        if R['context'] or R['below_claim']:
            print(f"  #{issue}: {len(R['matched'])} claims matched by a file, {R['context']} context packings"
                  + (f"; BELOW their claim: {[(n_, round(s_, 12), c_) for n_, s_, c_, _ in R['below_claim']]}"
                     if R['below_claim'] else ''))
        if claims and refs_ and not seen:
            print(f"  #{issue}: WARNING {len(claims)} claims but no packing files at {', '.join(refs_)}")
        ref = refs_[0] if refs_ else f"issue #{issue}"
        # 10-17: rows from an earlier parse that the issue no longer claims (e.g. maintainer tables in comments, #573)
        gone = [r['n'] for r in st.db.execute('SELECT n FROM frontier WHERE source=?', (src,)) if r['n'] not in claims]
        if gone:
            st.db.execute(f'DELETE FROM frontier WHERE source=? AND n IN ({",".join("?" * len(gone))})', (src, *gone))
            print(f"  #{issue}: dropped {len(gone)} frontier rows no longer claimed: {sorted(gone)}")
        for n, s in sorted(claims.items()):
            row = st.query(n=n, source=src, limit=1)
            pid = R['matched'].get(n) or (row[0]['id'] if row else None)
            st.set_frontier(n, s, f'{s}', finder, src, ref, state == 'OPEN', pid)
            report.append((issue, n, s, row[0]['s'] if row else None, n in R['matched']))
        st.db.commit()
        print(f"#{issue} {author} ({state}): {len(claims)} claims, {len(E)} entries; files seen {seen}"
              + ('' if refs_ else ' (no repo or attachment)'))
    bad = [r for r in report if not r[4]]
    for r in bad:
        print(f'  claim not matched by a stored file: #{r[0]} n={r[1]} claim {r[2]} best file {r[3]}')
    return report


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
