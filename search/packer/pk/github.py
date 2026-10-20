"""Every network read the store's ingestion makes (10-17: moved here from sources.py / watch.py): the GitHub API (gh),
commit resolution, repo files (git), issue attachments, issues and their comments, repo HEADs.

Untrusted data (pk/README.md, Security): repos live as bare partial git clones runs/ext/git/<owner>-<repo>.git (10-20;
before: a tarball snapshot per commit), read with git cat-file, size-capped, PROVENANCE.txt line per fetch, nothing
executed; attachments
only from github.com/user-attachments, redirected only to GitHub's object hosts, size-capped.  Issue text is returned
as data for pk/issues.py's regexes; nothing here interprets it.
"""
from __future__ import annotations
import contextlib, datetime as dt, glob, json, os, re, subprocess, tempfile

UPSTREAM = 'jlevy/squares'
# Caps for untrusted downloads (pk/README.md, Security)
MAX_TARBALL = 200 << 20          # bytes one repo fetch may add (git, MAX_FETCH)
MAX_MEMBER = 20 << 20            # one extracted file (bigger members are skipped, not extracted)
MAX_MEMBERS = 20000              # files read from one commit
MAX_UNPACKED = 1 << 30           # total bytes read from one commit
MAX_ATTACH = 20 << 20            # one issue attachment
REPO_NAME = re.compile(r'[A-Za-z0-9](?:[A-Za-z0-9-]{0,38})/[A-Za-z0-9._-]{1,100}')
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EXT = os.path.join(PKG, 'runs', 'ext')


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


# ---------- repos: one partial bare git clone per repo (10-20; Evan: "use git as it is meant to be used") ----------
# Replaces per-commit tarball snapshots (runs/ext/<owner>-<repo>-<sha7>/, ~300 MB per commit of a big repo).  Each
# repo is a bare clone runs/ext/git/<owner>-<repo>.git; each needed commit is fetched shallow (depth 1) trees-only, then
# just the blobs we read (selected files, by id; already-present blobs are not refetched), so a new commit costs its
# trees plus its changed packing files.  (Plain depth-1 fetches with blobs resent a whole snapshot per commit.)  Untrusted data: bare repo (no work tree, no
# checkout), created with an empty template (no hooks) and every git call runs with core.hooksPath=/dev/null and lazy
# fetching off (a missing blob is skipped, never fetched behind our back).  Files are read with `git cat-file` and
# written as plain read-only files into a temp dir for the parsers; nothing is ever executed.
GIT_DIR = os.path.join(EXT, 'git')
MAX_FETCH = MAX_TARBALL          # bytes one fetch may add to a repo's git dir (killed beyond)
_GIT_ENV = dict(os.environ, GIT_NO_LAZY_FETCH='1', GIT_TERMINAL_PROMPT='0', GIT_CONFIG_NOSYSTEM='1',
                GIT_ASKPASS='/bin/false')
_SHA = re.compile(r'[0-9a-f]{40}')


def _git(gd, args, timeout=600, text=True, check=True, input=None):
    cmd = ['git', '-c', 'core.hooksPath=/dev/null', '-c', 'protocol.allow=never', '-c', 'protocol.https.allow=always',
           '-c', 'core.fsmonitor=false'] + (['--git-dir', gd] if gd else []) + args
    return subprocess.run(cmd, capture_output=True, text=text, check=check, timeout=timeout, env=_GIT_ENV, input=input)


def _dir_bytes(d):
    return sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(d) for f in fs
               if not os.path.islink(os.path.join(r, f)))


def _mirror(owner_repo):
    if not REPO_NAME.fullmatch(owner_repo):
        raise ValueError(f'bad repo {owner_repo!r}')
    owner, repo = owner_repo.split('/')
    gd = f'{GIT_DIR}/{owner}-{repo}.git'
    if not os.path.isdir(gd):
        os.makedirs(GIT_DIR, exist_ok=True)
        tmp = tempfile.mkdtemp(dir=GIT_DIR, prefix='.init-')
        _git(None, ['init', '--bare', '--quiet', '--template=', tmp])
        for k, v in (('remote.origin.url', f'https://github.com/{owner_repo}.git'), ('remote.origin.promisor', 'true'),
                     ('remote.origin.partialclonefilter', 'blob:none'), ('gc.auto', '0')):
            _git(tmp, ['config', k, v])
        os.rename(tmp, gd)
    return gd


def _has_commit(gd, sha):
    return _git(gd, ['cat-file', '-e', f'{sha}^{{commit}}'], check=False).returncode == 0


def _fetch(owner_repo, sha):
    """Make commit sha (40 hex) present in the repo's bare clone (shallow, blobs <= MAX_MEMBER).  Returns the git dir."""
    gd = _mirror(owner_repo)
    if _has_commit(gd, sha):
        return gd
    before = _dir_bytes(gd)
    cmd = ['git', '-c', 'core.hooksPath=/dev/null', '-c', 'protocol.allow=never', '-c', 'protocol.https.allow=always',
           '--git-dir', gd, 'fetch', '--quiet', '--no-tags', '--no-write-fetch-head', '--depth', '1',
           '--filter=blob:none', 'origin', sha]
    _run_capped(cmd, gd, f'{owner_repo}@{sha[:7]}')
    if not _has_commit(gd, sha):
        raise ValueError(f'{owner_repo}@{sha[:7]}: commit not present after fetch')
    owner, repo = owner_repo.split('/')
    with open(f'{EXT}/PROVENANCE.txt', 'a') as f:
        f.write(f'\n{owner}-{repo}@{sha[:7]} — git fetch {owner_repo}@{sha} (depth 1, trees only) into '
                f'git/{owner}-{repo}.git by pk/github.py (data only, bare repo, nothing executed; '
                f'{_dir_bytes(gd) - before} bytes)\n')
    return gd


def _run_capped(cmd, gd, what, cap=None):
    """Run a git fetch, killing it once the git dir grows by more than cap (MAX_FETCH) bytes or after 900 s."""
    before, cap = _dir_bytes(gd), cap or MAX_FETCH
    with subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, env=_GIT_ENV) as pr:
        t0 = dt.datetime.now()
        while pr.poll() is None:
            try:
                pr.wait(timeout=2)
            except subprocess.TimeoutExpired:
                pass
            if _dir_bytes(gd) - before > cap or (dt.datetime.now() - t0).total_seconds() > 900:
                pr.kill(); pr.wait()
                raise ValueError(f'{what}: fetch over cap ({cap} bytes) or 900 s')
        if pr.returncode != 0:
            raise subprocess.CalledProcessError(pr.returncode, cmd[-6:], stderr=pr.stderr.read().decode()[:300])


def _want_blobs(owner_repo, gd, blobs):
    """Fetch the listed blobs that are not present yet (one request per 500), within MAX_FETCH per call.  Commits are
    fetched trees-only (blob:none), so a new commit costs its trees plus the blobs we read that changed."""
    if not blobs:
        return
    chk = _git(gd, ['cat-file', '--batch-check'], input=''.join(f'{b}\n' for b in blobs)).stdout.splitlines()
    miss = [b for b, line in zip(blobs, chk) if line.endswith(' missing')]
    before = _dir_bytes(gd)
    for i in range(0, len(miss), 500):
        _run_capped(['git', '-c', 'core.hooksPath=/dev/null', '-c', 'protocol.allow=never', '-c',
                     'protocol.https.allow=always', '--git-dir', gd, 'fetch', '--quiet', '--no-tags',
                     '--no-write-fetch-head', '--filter=blob:none', 'origin', *miss[i:i + 500]], gd,
                    f'{owner_repo} blobs', cap=MAX_FETCH - (_dir_bytes(gd) - before))
    if miss:
        owner, repo = owner_repo.split('/')
        with open(f'{EXT}/PROVENANCE.txt', 'a') as f:
            f.write(f'\n{owner}-{repo} — {len(miss)} blobs fetched into git/{owner}-{repo}.git by pk/github.py '
                    f'({_dir_bytes(gd) - before} bytes; data only, nothing executed)\n')


def _full_sha(owner_repo, commit):
    if not REPO_NAME.fullmatch(owner_repo) or not re.fullmatch(r'[0-9a-fA-F]{7,40}', commit or ''):
        raise ValueError(f'bad repo / commit {owner_repo!r} {commit!r}')
    commit = commit.lower()
    if len(commit) != 40:
        gd = _mirror(owner_repo)               # a short sha we already have needs no API call
        r = _git(gd, ['rev-parse', '--verify', '--quiet', f'{commit}^{{commit}}'], check=False)
        full = r.stdout.strip() if r.returncode == 0 else resolve_sha(owner_repo, commit)
        if not full:
            raise ValueError(f'{owner_repo}: commit {commit} not found')
        commit = full
    if not _SHA.fullmatch(commit):
        raise ValueError(f'bad commit {commit!r}')
    return commit


def _safe_path(p):
    parts = p.split('/')
    return (p and not p.startswith('/') and '\x00' not in p and
            all(x not in ('', '.', '..') for x in parts) and parts[0] != '.git')


def repo_files(owner_repo, commit, name_re=None, skip_re=None, changed_since=None):
    """[(path, blob sha, size)] of regular files at owner/repo@commit whose basename matches name_re and not skip_re;
    changed_since (a commit) -> only files added / modified since it (all files if it cannot be fetched).  Blobs over
    MAX_MEMBER or missing after the blob fetch) are left out; at most MAX_MEMBERS files, MAX_UNPACKED bytes."""
    sha = _full_sha(owner_repo, commit)
    gd = _fetch(owner_repo, sha)
    only = None
    if changed_since:
        try:
            base = _full_sha(owner_repo, changed_since)
            _fetch(owner_repo, base)
            out = _git(gd, ['diff-tree', '-r', '-z', '--no-renames', '--name-only', '--diff-filter=ACMRT', base, sha],
                       text=False).stdout
            only = {x.decode('utf-8', 'replace') for x in out.split(b'\x00') if x}
        except Exception:
            only = None                           # unknown base: import everything (as a first import)
    out = _git(gd, ['ls-tree', '-r', '-z', '--full-tree', sha], text=False).stdout
    cand = []
    for ent in out.split(b'\x00'):
        if not ent:
            continue
        meta, _, path = ent.partition(b'\t')
        mode, typ, blob = meta.decode().split()
        path = path.decode('utf-8', 'replace')
        if typ != 'blob' or mode not in ('100644', '100755') or not _safe_path(path):
            continue                              # symlinks, submodules, odd paths: never materialised
        b = path.rsplit('/', 1)[-1]
        if (name_re and not name_re.search(b)) or (skip_re and skip_re.search(b)) or (only is not None and path not in only):
            continue
        cand.append((path, blob))
    if len(cand) > MAX_MEMBERS:
        raise ValueError(f'{owner_repo}@{sha[:7]}: {len(cand)} files > {MAX_MEMBERS}')
    _want_blobs(owner_repo, gd, sorted({b for _, b in cand}))
    chk = _git(gd, ['cat-file', '--batch-check'], input=''.join(f'{b}\n' for _, b in cand)).stdout.splitlines()
    R, total = [], 0
    for (path, blob), line in zip(cand, chk):
        f = line.split()
        if len(f) != 3 or f[1] != 'blob' or int(f[2]) > MAX_MEMBER:
            continue                              # missing (over the fetch filter) or too big: skipped
        total += int(f[2])
        if total > MAX_UNPACKED:
            raise ValueError(f'{owner_repo}@{sha[:7]}: files over {MAX_UNPACKED} bytes')
        R.append((path, blob, int(f[2])))
    return R


def read_blob(owner_repo, blob):
    gd = _mirror(owner_repo)
    if not _SHA.fullmatch(blob):
        raise ValueError(f'bad blob {blob!r}')
    return _git(gd, ['cat-file', 'blob', blob], text=False).stdout


def read_file(owner_repo, commit, path):
    """bytes of one file at owner/repo@commit (fetched if needed; MAX_MEMBER cap)."""
    if not _safe_path(path):
        raise ValueError(f'bad path {path!r}')
    sha = _full_sha(owner_repo, commit)
    gd = _fetch(owner_repo, sha)
    blob = _git(gd, ['rev-parse', '--verify', f'{sha}:{path}']).stdout.strip()
    _want_blobs(owner_repo, gd, [blob])
    sz = int(_git(gd, ['cat-file', '-s', blob]).stdout)
    if sz > MAX_MEMBER:
        raise ValueError(f'{path}: {sz} bytes > {MAX_MEMBER}')
    return read_blob(owner_repo, blob)


@contextlib.contextmanager
def checkout(owner_repo, commit, name_re=None, skip_re=None, changed_since=None):
    """Temp dir holding just the selected files of owner/repo@commit (repo_files' selection) as plain read-only files,
    at their repo paths, so path-based importers (sources.import_paths, globs) work unchanged.  Removed on exit."""
    files = repo_files(owner_repo, commit, name_re, skip_re, changed_since)
    os.makedirs(EXT, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=EXT, prefix='.co-') as d:
        gd = _mirror(owner_repo)
        if files:                                 # one cat-file --batch stream for all blobs
            pr = subprocess.run(['git', '-c', 'core.hooksPath=/dev/null', '--git-dir', gd, 'cat-file', '--batch'],
                                input=''.join(f'{b}\n' for _, b, _ in files).encode(), capture_output=True, check=True,
                                env=_GIT_ENV, timeout=600).stdout
            i = 0
            for path, blob, size in files:
                j = pr.index(b'\n', i)
                hdr = pr[i:j].split()
                if hdr[0].decode() != blob or int(hdr[2]) != size:
                    raise ValueError(f'cat-file stream out of step at {path}')
                data = pr[j + 1:j + 1 + size]; i = j + 1 + size + 1
                dest = os.path.join(d, path)
                if not os.path.realpath(dest).startswith(os.path.realpath(d) + os.sep):
                    continue
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                with open(dest, 'wb') as f:
                    f.write(data)
                os.chmod(dest, 0o444)
        yield d


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


# ---------- issues ----------
def issue(number):
    """(issue json, [(login, body)] comments oldest first) for jlevy/squares#number"""
    r = json.loads(_gh([f'repos/{UPSTREAM}/issues/{int(number)}']))
    return r, comments(number) if r.get('comments') else []


def comments(number):
    """[(login, body)] of an issue's comments, oldest first"""
    out = []
    for line in _gh(['--paginate', f'repos/{UPSTREAM}/issues/{int(number)}/comments?per_page=100',
                     '--jq', '[.[] | [.user.login, .body]]']).splitlines():
        if line.strip():
            out += [(w, b or '') for w, b in json.loads(line)]
    return out


def updated_since(cursor, pages=10):
    """issues and PRs of jlevy/squares updated after `cursor` (ISO Z; None = the last day), newest update first"""
    since = cursor or (dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=1)).strftime('%Y-%m-%dT%H:%M:%SZ')
    out = []
    for page in range(1, pages + 1):
        rows = json.loads(_gh([f'repos/{UPSTREAM}/issues?state=all&sort=updated&direction=desc&per_page=100'
                               f'&page={page}&since={since}']))
        out += [r for r in rows if r['updated_at'] > since]
        if len(rows) < 100:
            break
    return out


def head(repo):
    """full sha of a repo's default-branch HEAD, or None"""
    return resolve_sha(repo, 'HEAD')
