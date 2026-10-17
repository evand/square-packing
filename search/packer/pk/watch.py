"""Watch jlevy/squares issues and other people's packing repos; import every packing they publish (10-14).

  pk.py watch            dry run: what is new and what --apply would do (GitHub API reads only, nothing fetched)
  pk.py watch --apply    fetch + import, append pending.yaml entries for new issues, save the state

(a) Issues: jlevy/squares issues (open, or closed within --recent days) numbered above the last one seen and not already
    in pending.yaml.  From each body: github.com/OWNER/REPO references (plus /commit|tree|blob/<sha>[/path] and
    `@sha` / `commit <sha>` mentions), claims (markdown table: sources.body_claims; `s(n) <= x` in the title), and
    github user-attachment file links.  Proposed pending.yaml entry: pin = the commit the issue names (checked to exist
    on that repo), else the default-branch HEAD now (the comment says which); paths = the blob / tree paths the issue
    links, else the whole repo.  Issues with no repo, attachment or claim are recorded as skipped.  --apply appends the
    entries and runs sources.import_pending on them (source pending:#N, as for hand-written entries).
(b) Repos: every repo in pending.yaml (except ours) plus repos found in new issues.  New commits since the last import
    (first time: since the newest pending.yaml pin) are fetched, newest MAX_COMMITS per repo, and EVERY parseable
    packing file is imported (strict parse), source watch:<owner>/<repo>@<sha7>, kind 'watch'.

State: runs/watch/state.json (last issue seen, per-repo last imported sha); one line per --apply in runs/watch/log.jsonl.
Security (pk/README.md): issue bodies and repo files are untrusted data.  Only regexes read issue bodies; the only links
followed are github.com/OWNER/REPO (via the GitHub API) and github user-attachment files; nothing fetched is executed;
downloads are size-capped (sources.MAX_*); files are parsed only by pk.packing.parse(strict=True).  Nothing is posted.
"""
from __future__ import annotations
import collections, datetime as dt, json, os, re, subprocess, urllib.parse

import yaml

from . import sources
from .sources import REPO_NAME, ATTACH_RE, body_claims, resolve_sha, _gh

HERE = os.path.dirname(os.path.abspath(__file__))
BACKTICK = re.compile(r'`([A-Za-z0-9._+/-]{1,200}\.[A-Za-z0-9]{1,10})`')   # `file.ext` or `dir/file.ext` in issue bodies
PENDING = os.path.join(HERE, 'pending.yaml')
WDIR = os.path.join(HERE, '..', 'runs', 'watch')
STATE = os.path.join(WDIR, 'state.json')
UPSTREAM = 'jlevy/squares'
OWN = {'evand/square-packing'}                      # ours: never watched (pending entries for our issues are by hand)
MAX_COMMITS = 20                                    # commits fetched per repo per run (newest first)
RESERVED = {'user-attachments', 'orgs', 'settings', 'issues', 'pulls', 'marketplace', 'sponsors', 'apps', 'features',
            'topics', 'notifications', 'login', 'about', 'pricing', 'enterprise', 'collections', 'search',
            'advisories', 'site', 'contact', 'explore', 'trending', 'new', 'organizations', 'join', 'security'}
# github.com/OWNER/REPO[/commit|tree|blob/REF[/path]]: the only link form read from issue bodies
GH_RE = re.compile(r'(?<![\w./-])(?:https?://)?(?:www\.)?github\.com/([A-Za-z0-9][A-Za-z0-9-]{0,38})/([A-Za-z0-9._-]{1,100})'
                   r'(?:/(commit|tree|blob)/([0-9A-Za-z._-]{1,100})((?:/[^\s/?#)\]>"\'`]{1,200}){0,20}))?'
                   r'(?:@([0-9a-f]{7,40})\b)?')
SHA_MENTION = re.compile(r'(?:commit|pinned at|pin|sha)[\s:`*]*([0-9a-f]{7,40})\b', re.I)
TEST_DIRS = re.compile(r'(^|/)(tests?|test_?data|fixtures|examples?)/', re.I)   # checker test inputs: not packings
TITLE_CLAIM = re.compile(r's\(\s*(\d{1,5})\s*\)\s*(?:≤|<=|\\le)\s*(\d+\.\d{6,})')


def now():
    return dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def _clean(t, k=90):
    """Untrusted text for a yaml comment / terminal: one line, printable, short."""
    return re.sub(r'[^\x20-\x7e]', '?', ' '.join(str(t).split()))[:k]


def load_state():
    try:
        return json.load(open(STATE))
    except FileNotFoundError:
        return dict(last_issue=None, repos={}, issues={})


def save_state(S):
    os.makedirs(WDIR, exist_ok=True)
    tmp = STATE + '.tmp'
    json.dump(S, open(tmp, 'w'), indent=1, sort_keys=True)
    os.replace(tmp, STATE)


def pending_cfg():
    return yaml.safe_load(open(PENDING)) or []


# ---------- (a) issues ----------
def list_issues(after, recent_days, pages=10):
    """jlevy/squares issues (not PRs) numbered > after: open, or closed within recent_days.  Newest first.
    Also returns PRs by people other than the owner (reported, not handled)."""
    cutoff = (dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=recent_days)).strftime('%Y-%m-%dT%H:%M:%SZ')
    out, prs = [], []
    for page in range(1, pages + 1):
        rows = json.loads(_gh([f'repos/{UPSTREAM}/issues?state=all&sort=created&direction=desc&per_page=100&page={page}']))
        for r in rows:
            if r['number'] <= after:
                continue
            if r.get('pull_request'):
                if r['user']['login'] != UPSTREAM.split('/')[0]:
                    prs.append(r)
                continue
            if r['state'] == 'open' or (r.get('closed_at') or '') >= cutoff:
                out.append(r)
        if not rows or rows[-1]['number'] <= after:
            break
    return out, prs


def refs(text):
    """{repo: {'shas': [...], 'paths': [...]}} for github.com/OWNER/REPO links in text (upstream, ours-as-upstream and
    reserved paths excluded), in order of appearance; plus standalone sha mentions."""
    R = collections.OrderedDict()
    for m in GH_RE.finditer(text or ''):
        owner, repo, kind, ref, path, at = m.groups()
        repo = re.sub(r'(\.git)?\.*$', '', repo)
        name = f'{owner}/{repo}'
        if owner.lower() in RESERVED or name.lower() == UPSTREAM or not repo or not REPO_NAME.fullmatch(name):
            continue
        e = R.setdefault(name, dict(shas=[], paths=[], trees=[]))
        if kind and re.fullmatch(r'[0-9a-f]{7,40}', ref):
            e['shas'].append(ref)
        if at:
            e['shas'].append(at)
        if kind in ('tree', 'blob') and path:
            p = urllib.parse.unquote(path.strip('/')).rstrip('.,;:')
            if p and '..' not in p.split('/') and re.fullmatch(r'[A-Za-z0-9._/ +-]{1,300}', p):
                e['paths'].append(p)
                if kind == 'tree':
                    e['trees'].append(p)
    loose = [s for s in SHA_MENTION.findall(text or '')]
    return R, loose


def propose(issue):
    """One issue -> (entry dict for pending.yaml or None, note)."""
    body, title = issue.get('body') or '', issue.get('title') or ''
    R, loose = refs(title + '\n' + body)
    R = collections.OrderedDict((k, v) for k, v in R.items() if k not in OWN or issue['user']['login'] == 'evand')
    claims = body_claims(body)
    for n, s in TITLE_CLAIM.findall(title):
        claims[int(n)] = min(claims.get(int(n), 9e9), float(s))
    att = sorted(set(m.group(0) for m in ATTACH_RE.finditer(body)))
    e = dict(issue=issue['number'], finder=issue['user']['login'])
    note = []
    if R:
        repo, info = next(iter(R.items()))              # first repo linked is the subject (others listed in the note)
        if len(R) > 1:
            note.append(f'also links {", ".join(list(R)[1:4])}')
        pin = None
        for c in info['shas'] + loose:
            pin = resolve_sha(repo, c)
            if pin:
                note.append(f'pin = commit {c[:12]} named in the issue')
                break
        if not pin:
            br = json.loads(_gh([f'repos/{repo}', '--jq', '{b: .default_branch}']))['b']
            pin = resolve_sha(repo, br)
            named = info['shas'] + loose
            note.append(f'pin = {br} HEAD at {now()} ' +
                        (f'(named commit {named[0][:12]} not found)' if named else '(issue names no commit)'))
        if not pin:
            return None, f'repo {repo} not reachable'
        # 10-17: a tree link (a directory) used to pin the whole dir, so every older packing in it became a sighting
        # of this issue (#659 / #664 tagged ~40 each).  When the body names files in backticks, narrow each linked dir
        # to those names (bare names under the dir, or full paths inside it); missing names just glob to nothing.
        named = [x.strip() for x in BACKTICK.findall(body)]
        paths = []
        for p in info['paths']:
            if p in info.get('trees', []):
                sub = [f'{p}/{x}' for x in named if '/' not in x] + [x for x in named if x.startswith(p + '/')]
                if sub:
                    paths += sub
                    note.append(f'tree {p} narrowed to {len(sub)} named files')
                    continue
            paths.append(p)
        e.update(repo=repo, commit=pin, paths=sorted(set(paths)) or ['.'])
    else:
        e.update(repo=None)
    if claims:
        e['claims'] = {int(n): float(s) for n, s in sorted(claims.items())}
    if att:
        e['attachments'] = att
    if not R and not claims and not att:
        return None, 'no repo, attachment or claim'
    return e, '; '.join(note)


def append_pending(entries):
    """Insert entries (newest first) above the existing ones in pending.yaml, one flow-style line each + comment."""
    lines = open(PENDING).read().split('\n')
    i = next((k for k, l in enumerate(lines) if l.startswith('- ')), len(lines))
    new = []
    for e, note, title in entries:
        y = '- ' + yaml.safe_dump(e, default_flow_style=True, sort_keys=False, width=100000).strip()
        new.append(f'{y}  # auto (pk.py watch {now()[:10]}): {_clean(note, 120)}; title: {_clean(title)}')
    open(PENDING, 'w').write('\n'.join(lines[:i] + new + lines[i:]))


# ---------- (b) repos ----------
def watch_list(S, cfg):
    """{repo: {'base': newest pending pin or None, 'finder': name}} for every watched repo."""
    W = {}
    for e in sorted(cfg, key=lambda e: e.get('issue') or 0):
        r = e.get('repo')
        if r and r not in OWN and REPO_NAME.fullmatch(r):
            W[r] = dict(base=str(e['commit']), finder=e.get('finder'))
    for r, v in S['repos'].items():
        W.setdefault(r, dict(base=None, finder=v.get('finder')))
    return W


def poll_repo(repo, base):
    """(head sha, default branch, [new commit shas newest first] or None if base unknown / not comparable)."""
    br = json.loads(_gh([f'repos/{repo}', '--jq', '{b: .default_branch}']))['b']
    head = resolve_sha(repo, br)
    if not base or not head:
        return head, br, None
    if head.startswith(base) or base.startswith(head):
        return head, br, []
    try:
        c = json.loads(_gh([f'repos/{repo}/compare/{base}...{head}', '--jq', '{s: .status, c: [.commits[].sha]}']))
    except subprocess.CalledProcessError:
        return head, br, None
    return head, br, list(reversed(c['c'])) if c['s'] in ('ahead', 'diverged') else []


# ---------- driver ----------
def run(st=None, apply=False, recent_days=7, max_commits=MAX_COMMITS):
    S = load_state()
    cfg = pending_cfg()
    known = {e['issue'] for e in cfg if e.get('issue')}
    if S['last_issue'] is None:                     # first run: everything up to the newest pending.yaml issue was seen
        S['last_issue'] = max(known)
    started = now()
    started_local = dt.datetime.now().isoformat(timespec='seconds')     # the store's sighting.added clock (local time)
    log = dict(start=started, start_local=started_local, apply=apply, issues=[], skipped_issues=[], prs=[], repos=[], counts={})

    # (a) issues
    issues, prs = list_issues(S['last_issue'], recent_days)
    props = []
    for r in sorted(issues, key=lambda r: r['number']):
        if r['number'] in known:
            continue
        try:
            e, note = propose(r)
        except (subprocess.CalledProcessError, ValueError, KeyError) as ex:
            e, note = None, f'error {_clean(repr(ex), 80)}'
        tag = f"#{r['number']} {r['user']['login']} ({r['state']}) {_clean(r['title'], 70)}"
        if e is None:
            print(f'issue {tag}: skip ({note})')
            log['skipped_issues'].append(dict(issue=r['number'], why=note))
        else:
            print(f'issue {tag}: propose {json.dumps(e)[:300]}' + (f'  [{note}]' if note else ''))
            props.append((e, note, r['title']))
            log['issues'].append(dict(entry=e, note=note))
    for r in prs:
        print(f"PR #{r['number']} {r['user']['login']} {_clean(r['title'], 70)}: not handled (PRs are not watched)")
        log['prs'].append(r['number'])

    # (b) repos
    W = watch_list(S, cfg)
    for e, _, _ in props:
        if e.get('repo') and e['repo'] not in OWN:
            W.setdefault(e['repo'], dict(base=e['commit'], finder=e['finder']))
    plan = []
    for repo, w in sorted(W.items()):
        prev = S['repos'].get(repo, {}).get('imported')
        base = prev or w['base']
        try:
            head, br, new = poll_repo(repo, base)
        except subprocess.CalledProcessError as ex:
            print(f'repo {repo}: not reachable ({_clean(ex.stderr or ex, 80)})'); continue
        if head is None:
            print(f'repo {repo}: no HEAD'); continue
        if prev and head.startswith(prev[:7]):
            print(f'repo {repo}: no new commits ({head[:7]})'); continue
        if new is None:
            todo = [head]
        else:
            todo = ([head] if not prev and not new else new)[:max_commits]
        print(f'repo {repo}: {br} {head[:7]}; base {base[:7] if base else None}'
              f'{"" if prev else " (newest pending pin; first watch import)"}; {len(new) if new is not None else "?"} '
              f'new commits; import {len(todo)}: {" ".join(c[:7] for c in todo)}')
        plan.append((repo, w, head, todo, len(new) if new is not None else None))
        log['repos'].append(dict(repo=repo, head=head, base=base, new_commits=len(new) if new is not None else None,
                                 fetch=todo))

    if not apply:
        print(f'dry run: {len(props)} new issue entries, {sum(len(p[3]) for p in plan)} commits to fetch in '
              f'{len(plan)} repos; --apply to fetch and import')
        return log

    from .store import Store
    st = st or Store()
    C = collections.Counter()
    coll = []
    for repo, w, head, todo, nnew in plan:
        for sha in todo:
            try:
                d = sources.fetch_repo(repo, sha)
            except Exception as ex:
                print(f'  {repo}@{sha[:7]}: fetch failed ({_clean(repr(ex), 120)})'); C['fetch_failed'] += 1; continue
            c = sources.import_paths(st, [d], f'watch:{repo}@{sha[:7]}', w.get('finder') or repo.split('/')[0],
                                     f'{repo}@{sha[:7]}', 'watch', quiet=True, strict=True, collect=coll,
                                     skip_path=TEST_DIRS)
            print(f'  {repo}@{sha[:7]}: {c}')
            C['commits'] += 1
            for k, v in c.items():
                C[k] += v
        S['repos'][repo] = dict(imported=head, finder=w.get('finder'), checked=now())
    if props:
        append_pending(props)
        nums = [e['issue'] for e, _, _ in props]
        print(f'pending.yaml: appended {nums}; importing them (sync-pending)')
        before = len(coll)
        n0 = st.db.execute('SELECT max(id) FROM packing').fetchone()[0] or 0
        sources.import_pending(st, only=nums)
        for r in st.db.execute('SELECT id, n, s, status FROM packing WHERE id > ?', (n0,)):
            coll.append((f'pending', r['n'], r['s'], r['id'], True))
        C['pending_new'] = len(coll) - before
    for e, _, _ in props:
        S['issues'][str(e['issue'])] = dict(at=now(), action='appended')
    for d in log['skipped_issues']:
        S['issues'][str(d['issue'])] = dict(at=now(), action='skipped: ' + d['why'])
    nums = [r['number'] for r in issues]
    if nums:
        S['last_issue'] = max(S['last_issue'], max(nums))
    S['last_poll'] = started
    save_state(S)
    ids = {pid for _, _, _, pid, new in coll if new}
    rows = [st.db.execute('SELECT n, status FROM packing WHERE id=?', (i,)).fetchone() for i in ids]
    byn = collections.Counter(r['n'] for r in rows if r['status'] != 'infeasible')
    C['new_distinct'] = sum(byn.values()); C['new_infeasible'] = sum(1 for r in rows if r['status'] == 'infeasible')
    log['counts'] = dict(C); log['new_by_n'] = dict(sorted(byn.items())); log['end'] = now()
    with open(os.path.join(WDIR, 'log.jsonl'), 'a') as f:
        f.write(json.dumps(log) + '\n')
    print(f'watch --apply: {dict(C)}')
    print('new distinct feasible packings by n:', ' '.join(f'{n}:{k}' for n, k in sorted(byn.items())) or 'none')
    print(f'first sighting of this run: {started_local} local (python3 fresh.py new --since {started_local[:16]})')
    return log
