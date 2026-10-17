"""Keep the store current with jlevy/squares issues and the repos they link (10-14; rebuilt 10-17).

  pk.py watch                    dry run: updated issues and what each would import; watched repos' new commits
  pk.py watch --apply            sync them (fetch, import, reconcile), import new repo commits, record the sync state
  pk.py sync-issues [N ...] [--all]   re-sync given issues (or every issue the store has seen + pk/overrides.yaml)

Issues.  Every issue updated since the cursor is synced from its CURRENT text: github.issue -> issues.parse ->
issues.units (+ pk/overrides.yaml) -> each unit's files imported (source `pending:#N`, kind 'pending'), attachments
too -> reconcile (each of the issue's packings 'claim' or 'context'; claims with no matching file and files below their
claim reported) -> frontier rows for its claims (rows it no longer claims dropped).  Nothing parsed is kept by hand:
the `issue` table records each sync (units, claims, notes) for reading, and is rewritten on the next one.
Pins: a named commit resolves once and is cached (`resolved` table); "HEAD" and branch names resolve at an issue's
first sight and stay pinned for that issue.

Repos.  Every repo an issue unit names (except ours) is watched: new commits on its default branch since the last
import (first time: its HEAD), newest MAX_COMMITS, EVERY parseable packing file imported (source
`watch:<owner>/<repo>@<sha7>`, kind 'watch'; test / example dirs skipped).

State is all in the store (`sync` table): `watch.issues_updated` (cursor), `watch.poll`, `repo:<name>` (last imported
HEAD), `register.sha` (sync-register).  probe() compares them with GitHub (one call for the register, one for issues,
one per repo) for candidate.live.  One line per --apply in runs/watch/log.jsonl.

Security (pk/README.md): issue text and repo files are untrusted data, read only by pk/issues.py's regexes and
pk.packing.parse(strict=True); network reads only through pk/github.py (size caps, read-only dirs, nothing executed).
Nothing is posted.
"""
from __future__ import annotations
import collections, datetime as dt, glob, json, os, re, subprocess

import yaml

from . import github, issues
from .sources import import_paths

HERE = os.path.dirname(os.path.abspath(__file__))
OVERRIDES = os.path.join(HERE, 'overrides.yaml')
WDIR = os.path.join(HERE, '..', 'runs', 'watch')
MAX_COMMITS = 20                                    # commits fetched per repo per run (newest first)
TEST_DIRS = re.compile(r'(^|/)(tests?|test_?data|fixtures|examples?)/', re.I)   # checker test inputs: not packings
CLAIM_TOL = 1e-6        # a file matches a claim when claim - CLAIM_TOL <= s <= claim + 1e-9 (claims are rounded up)

SCHEMA = """
CREATE TABLE IF NOT EXISTS issue (
  number INTEGER PRIMARY KEY, author TEXT, state TEXT, title TEXT, updated TEXT, synced TEXT, units TEXT,
  claims TEXT, notes TEXT, created TEXT, credit TEXT);
CREATE TABLE IF NOT EXISTS resolved (
  repo TEXT NOT NULL, ref TEXT NOT NULL, issue INTEGER NOT NULL, sha TEXT NOT NULL, at TEXT, UNIQUE(repo, ref, issue));
"""


def now():
    return dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def _clean(t, k=90):
    """Untrusted text for the terminal: one line, printable, short."""
    return re.sub(r'[^\x20-\x7e]', '?', ' '.join(str(t).split()))[:k]


def overrides():
    """{issue: override dict} from pk/overrides.yaml"""
    return {int(k): v or {} for k, v in (yaml.safe_load(open(OVERRIDES)) or {}).items()}


def ensure_schema(st):
    st.db.executescript(SCHEMA)
    have = {r[1] for r in st.db.execute('PRAGMA table_info(issue)')}
    for col in ('created', 'credit'):                     # 10-17: added after the first issue syncs
        if col not in have:
            st.db.execute(f'ALTER TABLE issue ADD COLUMN {col} TEXT')


def resolver(st):
    """resolve(repo, ref, issue) -> full sha or None, cached in the store: a commit ref (hex) once for everyone; HEAD or
    a branch name once per issue (its first sight).  Misses are not cached (a commit can be pushed later, #583)."""
    ensure_schema(st)

    def resolve(repo, ref, issue):
        key = 0 if issues.HEX.fullmatch(ref) else int(issue)
        r = st.db.execute('SELECT sha FROM resolved WHERE repo=? AND ref=? AND issue=?', (repo, ref, key)).fetchone()
        if r:
            return r['sha']
        try:
            sha = github.resolve_sha(repo, ref)
        except (subprocess.CalledProcessError, ValueError):
            sha = None
        if sha:
            st.db.execute('INSERT OR IGNORE INTO resolved(repo, ref, issue, sha, at) VALUES (?,?,?,?,?)',
                          (repo, ref, key, sha, now()))
        return sha
    return resolve


def reconcile(st, src, claims):
    """Sort an issue's packings (sightings with source `src`) into 'claim' (its n is claimed and its side matches the
    claim) and 'context' (anything else the linked files held: other counts, older packings, parents).  Every packing
    stays in the store and counts for record checks; the kind only says what the issue claimed.
    -> dict(matched={n: pid}, unmatched_claims={n: s}, context=k, below_claim=[(n, s, claim, pid)])"""
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


def sync_issue(st, number, ov=None, data=None, apply=True, log=print):
    """Sync one issue from its current text (see the module doc).  data = (issue json, comments) if already fetched.
    apply=False: parse and plan only (resolution still goes through the cache).  -> summary dict"""
    r, cm = data or github.issue(number)
    P = issues.parse(r, cm)
    o = (ov if ov is not None else overrides()).get(number, {})
    U, att, claims, notes = issues.units(P, resolver(st), o)
    tag = f"#{number} {P['author']} ({P['state']}) {_clean(P['title'], 60)}"
    S = dict(issue=number, units=U, attachments=att, claims=claims, notes=notes)
    if not apply:
        log(f'{tag}: {len(U)} units {[(u["repo"], u["commit"][:7], len(u["paths"])) for u in U]}, '
            f'{len(att)} attachments, {len(claims)} claims' + (f'; {"; ".join(notes)}' if notes else ''))
        return S
    from .people import credit
    src, finder = f'pending:#{number}', o.get('finder') or credit(st, P['author'])
    seen, refs_ = 0, []
    for u in U:
        ref = f"{u['repo']}@{u['commit'][:7]}"
        try:
            d = github.fetch_repo(u['repo'], u['commit'])
        except Exception as ex:
            notes.append(f'{ref}: fetch failed {_clean(repr(ex), 80)}'); continue
        paths = [p for g in u['paths'] for p in glob.glob(os.path.join(d, g), recursive=True)
                 if os.path.realpath(p).startswith(os.path.realpath(d))]
        c = import_paths(st, paths, src, finder, ref, 'pending', quiet=True, strict=True)
        seen += c['seen']; refs_.append(ref); u['files'] = c['seen']
    for url in att:
        try:
            f, sha = github.fetch_attachment(url, number)
        except Exception as ex:
            notes.append(f'attachment not fetched: {_clean(url, 80)} ({_clean(repr(ex), 60)})'); continue
        aref = f"#{number} attachment {os.path.basename(f)}"
        seen += import_paths(st, [f], src, finder, aref, 'pending', quiet=True, strict=True)['seen']; refs_.append(aref)
    R = reconcile(st, src, claims)
    gone = [x['n'] for x in st.db.execute('SELECT n FROM frontier WHERE source=?', (src,)) if x['n'] not in claims]
    if gone:
        st.db.execute(f'DELETE FROM frontier WHERE source=? AND n IN ({",".join("?" * len(gone))})', (src, *gone))
    ref0 = refs_[0] if refs_ else f'issue #{number}'
    for n, s in sorted(claims.items()):
        st.set_frontier(n, s, f'{s}', finder, src, ref0, P['state'] == 'open', R['matched'].get(n))
    st.db.execute('INSERT OR REPLACE INTO issue(number, author, state, title, updated, synced, units, claims, notes, '
                  'created, credit) VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                  (number, P['author'], P['state'], P['title'], P['updated'], now(), json.dumps(U), json.dumps(claims),
                   json.dumps(notes), P['created'], finder))
    st.db.commit()
    S.update(seen=seen, matched=sorted(R['matched']), unmatched=R['unmatched_claims'], context=R['context'],
             below=R['below_claim'], dropped=gone)
    flags = []
    if R['unmatched_claims'] and refs_:
        flags.append(f"claims without a matching file: {sorted(R['unmatched_claims'])}")
    if R['below_claim']:
        flags.append(f"files BELOW their claim: {[(n, round(s, 12), c) for n, s, c, _ in R['below_claim']]}")
    if gone:
        flags.append(f'dropped frontier rows no longer claimed: {sorted(gone)}')
    log(f'{tag}: {len(U)} units, {len(att)} attachments; {seen} files, {len(claims)} claims '
        f'({len(R["matched"])} matched), {R["context"]} context' + ''.join(f'; {x}' for x in flags + notes))
    return S


# ---------- repos ----------
def poll_repo(repo, base):
    """(head sha, default branch, [new commit shas newest first] or None if base unknown / not comparable)."""
    br = json.loads(github._gh([f'repos/{repo}', '--jq', '{b: .default_branch}']))['b']
    head = github.resolve_sha(repo, br)
    if not base or not head:
        return head, br, None
    if head.startswith(base) or base.startswith(head):
        return head, br, []
    try:
        c = json.loads(github._gh([f'repos/{repo}/compare/{base}...{head}', '--jq', '{s: .status, c: [.commits[].sha]}']))
    except subprocess.CalledProcessError:
        return head, br, None
    return head, br, list(reversed(c['c'])) if c['s'] in ('ahead', 'diverged') else []


def watched(st):
    """{repo: finder} : every repo (except ours) that has given us a packing, through an issue unit or an earlier
    watch import.  10-17: repos an issue links that hold no packings (code, Lean proofs, the maintainer's other repos)
    are imported once with the issue and not watched."""
    ensure_schema(st)
    W = {}
    for r in st.db.execute('SELECT credit, author, units FROM issue'):
        for u in json.loads(r['units'] or '[]'):
            if u['repo'] not in issues.OWN and u.get('files'):
                W.setdefault(u['repo'], r['credit'] or r['author'])
    for r in st.db.execute("SELECT DISTINCT substr(source, 7, instr(source, '@') - 7) AS repo FROM sighting "
                           "WHERE source LIKE 'watch:%'"):
        if r['repo'] and r['repo'] not in issues.OWN:
            W.setdefault(r['repo'], None)
    return W


def sync_repos(st, apply, max_commits=MAX_COMMITS, log=print):
    Y = st.get_sync('repo:')
    C = collections.Counter()
    for repo, finder in sorted(watched(st).items()):
        prev = (Y.get(f'repo:{repo}') or (None,))[0]
        try:
            head, br, new = poll_repo(repo, prev)
        except subprocess.CalledProcessError as ex:
            log(f'repo {repo}: not reachable ({_clean(ex.stderr or ex, 80)})'); continue
        if head is None or (prev and head.startswith(prev[:7])):
            continue
        todo = ([head] if new is None or not new else new)[:max_commits]
        log(f'repo {repo}: {br} {head[:7]}; last imported {prev[:7] if prev else None}; import {len(todo)}: '
            f'{" ".join(c[:7] for c in todo)}')
        if not apply:
            continue
        for sha in todo:
            try:
                d = github.fetch_repo(repo, sha)
            except Exception as ex:
                log(f'  {repo}@{sha[:7]}: fetch failed ({_clean(repr(ex), 120)})'); C['fetch_failed'] += 1; continue
            c = import_paths(st, [d], f'watch:{repo}@{sha[:7]}', finder or repo.split('/')[0], f'{repo}@{sha[:7]}',
                             'watch', quiet=True, strict=True, skip_path=TEST_DIRS)
            log(f'  {repo}@{sha[:7]}: {c}')
            C['commits'] += 1
            for k, v in c.items():
                C[k] += v
        st.set_sync(f'repo:{repo}', head)
    return C


# ---------- driver ----------
def run(st=None, apply=False, only=None, everything=False, repos=True, max_commits=MAX_COMMITS, log=print):
    """Updated issues since the cursor (or `only` issue numbers, or with everything=True every issue the store has seen
    plus every issue in overrides.yaml), then watched repos."""
    from .store import Store
    st = st or Store()
    ensure_schema(st)
    started = now()
    started_local = dt.datetime.now().isoformat(timespec='seconds')     # the store's sighting.added clock
    ov = overrides()
    cursor = (st.get_sync('watch.issues_updated').get('watch.issues_updated') or (None,))[0]
    if only or everything:
        nums = set(only or ())
        if everything:
            nums |= {r['number'] for r in st.db.execute('SELECT number FROM issue')} | set(ov)
        upd, todo = [], [(n, None) for n in sorted(nums)]
    else:
        upd = github.updated_since(cursor)
        todo = [(r['number'], r) for r in sorted(upd, key=lambda r: r['number']) if not r.get('pull_request')]
        for r in upd:
            if r.get('pull_request') and r['user']['login'] != github.UPSTREAM.split('/')[0]:
                log(f"PR #{r['number']} {r['user']['login']} {_clean(r['title'], 70)}: not handled (PRs are not watched)")
    rows = []
    for n, r in todo:
        try:
            data = (r, github.comments(n) if r.get('comments') else []) if r else None
            rows.append(sync_issue(st, n, ov, data, apply=apply, log=log))
        except (subprocess.CalledProcessError, ValueError, KeyError) as ex:
            log(f'#{n}: error {_clean(repr(ex), 120)}')
    C = sync_repos(st, apply, max_commits, log) if repos else collections.Counter()
    if not apply:
        log(f'dry run: {len(todo)} issues, nothing fetched or written beyond the resolution cache; --apply to import')
        return rows
    if only or everything:                          # an explicit re-sync does not move the cursor
        pass
    elif upd:
        st.set_sync('watch.issues_updated', max(r['updated_at'] for r in upd), commit=False)
    elif cursor is None:
        st.set_sync('watch.issues_updated', started, commit=False)
    st.set_sync('watch.poll', started)
    os.makedirs(WDIR, exist_ok=True)
    with open(os.path.join(WDIR, 'log.jsonl'), 'a') as f:
        f.write(json.dumps(dict(start=started, start_local=started_local, issues=[x['issue'] for x in rows],
                                counts=dict(C), end=now())) + '\n')
    log(f'watch --apply: {len(rows)} issues synced; repos {dict(C)}')
    log(f'first sighting of this run: {started_local} local (python3 fresh.py new --since {started_local[:16]})')
    return rows


def probe(st=None):
    """Is the store current?  Cheap reads only (one call for the register, one for issues, one per watched repo):
    -> dict(stale=[...], as_of={...}).  Stale = register main moved since the last sync-register, an issue / PR updated
    since the last watch --apply, or a watched repo whose HEAD is not the commit last imported."""
    from .store import Store
    st = st or Store()
    Y = st.get_sync()
    stale, as_of = [], {}
    reg = (Y.get('register.sha') or (None,))[0]
    main = github.resolve_sha(github.UPSTREAM, 'main') or ''
    as_of['register'] = f'{(reg or "never")[:8]} (main {main[:8]})'
    if main != reg:
        stale.append(f'register main {main[:8]} != synced {(reg or "never")[:8]}')
    cur = (Y.get('watch.issues_updated') or (None,))[0]
    as_of['issues_updated'] = cur
    as_of['last_issue'] = (st.db.execute('SELECT max(number) m FROM issue').fetchone() or {'m': None})['m'] \
        if st.db.execute("SELECT 1 FROM sqlite_master WHERE name='issue'").fetchone() else None
    if cur is None:
        stale.append('issues: never watched into this store')
    else:
        for r in github.updated_since(cur, pages=1):
            kind = 'PR' if r.get('pull_request') else 'issue'
            stale.append(f"{kind} #{r['number']} updated {r['updated_at']}: {_clean(r['title'], 60)}")
    W = watched(st)
    for k, (v, at) in sorted(Y.items()):
        if k.startswith('repo:') and k[5:] in W:
            h = github.head(k[5:])
            if h and v and not h.startswith(v[:7]):
                stale.append(f'repo {k[5:]}: HEAD {h[:7]} != imported {v[:7]}')
    for repo in sorted(set(W) - {k[5:] for k in Y if k.startswith('repo:')}):
        stale.append(f'repo {repo}: watched, never imported')
    as_of['repos'] = sum(1 for k in Y if k.startswith('repo:'))
    return dict(stale=stale, as_of=as_of)
