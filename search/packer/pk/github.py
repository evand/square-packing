"""Every network read the store's ingestion makes (10-17: moved here from sources.py / watch.py): the GitHub API (gh),
commit resolution, repo tarballs, issue attachments, issues and their comments, repo HEADs.

Untrusted data (pk/README.md, Security): repo tarballs go to runs/ext/<owner>-<repo>-<sha7>/ (size-capped download
and extraction, tarfile filter 'data', made read-only, a PROVENANCE.txt line appended, nothing executed); attachments
only from github.com/user-attachments, redirected only to GitHub's object hosts, size-capped.  Issue text is returned
as data for pk/issues.py's regexes; nothing here interprets it.
"""
from __future__ import annotations
import datetime as dt, glob, json, os, re, subprocess, tarfile, tempfile

UPSTREAM = 'jlevy/squares'
# Caps for untrusted downloads (pk/README.md, Security)
MAX_TARBALL = 200 << 20          # compressed tarball bytes
MAX_MEMBER = 20 << 20            # one extracted file (bigger members are skipped, not extracted)
MAX_MEMBERS = 20000              # tar entries
MAX_UNPACKED = 1 << 30           # total extracted bytes
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
        f.write(f'\n{owner}-{repo}-{commit[:7]} — tarball {owner_repo}@{commit} fetched by pk/github.py '
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
