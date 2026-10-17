"""Watch jlevy/squares issues and other people's packing repos; import every packing they publish (10-14).

  pk.py watch            dry run: what is new and what --apply would do (GitHub API reads only, nothing fetched)
  pk.py watch --apply    fetch + import, append pending.yaml entries for new issues, save the state

(a) Issues: every jlevy/squares issue updated since the cursor (10-17; see below) not already in pending.yaml gets an
    entry.  From its body and comments: github.com/OWNER/REPO references (plus /commit|tree|blob/<sha>[/path] and
    `@sha` / `commit <sha>` mentions), claims (markdown table: sources.body_claims; `s(n) <= x` in the title), and
    github user-attachment file links.  Proposed pending.yaml entry: pin = the commit the issue names (checked to exist
    on that repo), else the default-branch HEAD now (the comment says which); paths = the blob / tree paths the issue
    links, else the whole repo.  Issues with no repo, attachment or claim are recorded as skipped.  --apply appends the
    entries and runs sources.import_pending on them (source pending:#N, as for hand-written entries).
(b) Repos: every repo in pending.yaml (except ours) plus repos found in new issues.  New commits since the last import
    (first time: since the newest pending.yaml pin) are fetched, newest MAX_COMMITS per repo, and EVERY parseable
    packing file is imported (strict parse), source watch:<owner>/<repo>@<sha7>, kind 'watch'.

10-17 (Evan: import more packings, wherever they appear; record checks = the store + a check for anything newer):
(a) now reads every issue UPDATED since the cursor (state 'issues_updated'), body and comments together.  A new issue
    -> entry as before; a known issue -> a supplemental entry for attachments / lower claims that appeared later
    (supplement()).  Every repo linked anywhere in an updated issue is watched (whole-repo imports).  Tree links are not
    narrowed: sources.reconcile marks each issue's packings 'claim' (n claimed, side matches) or 'context'.
probe(): register main sha, issues / PRs updated since the cursor, watched repos' HEADs vs the store's sync table
    (store.set_sync, written by sync-register and watch --apply).  candidate.live() = probe (+ sync when stale) + store.
State: runs/watch/state.json (last issue seen, per-repo last imported sha); one line per --apply in runs/watch/log.jsonl.
Security (pk/README.md): issue bodies and repo files are untrusted data.  Only regexes read issue bodies; the only links
followed are github.com/OWNER/REPO (via the GitHub API) and github user-attachment files; nothing fetched is executed;
downloads are size-capped (sources.MAX_*); files are parsed only by pk.packing.parse(strict=True).  Nothing is posted.
"""
from __future__ import annotations
import collections, datetime as dt, json, os, re, subprocess, urllib.parse

import yaml

from . import sources
from .sources import REPO_NAME, ATTACH_RE, resolve_sha, _gh

HERE = os.path.dirname(os.path.abspath(__file__))
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
def updated_since(cursor, pages=10):
    """10-17: jlevy/squares issues and PRs updated after `cursor` (ISO Z; None = the last day), newest update first."""
    since = cursor or (dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=1)).strftime('%Y-%m-%dT%H:%M:%SZ')
    out = []
    for page in range(1, pages + 1):
        rows = json.loads(_gh([f'repos/{UPSTREAM}/issues?state=all&sort=updated&direction=desc&per_page=100'
                               f'&page={page}&since={since}']))
        out += [r for r in rows if r['updated_at'] > since]
        if len(rows) < 100:
            break
    return out


def comments(number):
    """[(login, body)] of an issue's comments, oldest first"""
    out = []
    for line in _gh(['--paginate', f'repos/{UPSTREAM}/issues/{int(number)}/comments?per_page=100',
                     '--jq', '[.[] | [.user.login, .body]]']).splitlines():
        if line.strip():
            out += [(w, b or '') for w, b in json.loads(line)]
    return out


def probe(st=None):
    """10-17: is the store current?  Cheap reads only (one call for the register, one for issues, one per watched repo):
    -> dict(stale=[...], as_of={...}).  Stale = register main moved since the last sync-register, an issue / PR updated
    since the last watch --apply, or a watched repo whose HEAD is not the commit last imported."""
    from .store import Store
    st = st or Store()
    Y = st.get_sync()
    S = load_state()
    stale, as_of = [], {}
    reg = (Y.get('register.sha') or (None,))[0]
    main = _gh([f'repos/{UPSTREAM}/commits/main', '--jq', '.sha']).strip()
    as_of['register'] = f'{(reg or "never")[:8]} (main {main[:8]})'
    if main != reg:
        stale.append(f'register main {main[:8]} != synced {(reg or "never")[:8]}')
    cur = (Y.get('watch.issues_updated') or (None,))[0]
    as_of['issues_updated'] = cur
    as_of['last_issue'] = S.get('last_issue')
    for r in updated_since(cur, pages=1) if cur else [None]:
        if r is None:
            stale.append('issues: never watched into this store'); break
        kind = 'PR' if r.get('pull_request') else 'issue'
        stale.append(f"{kind} #{r['number']} updated {r['updated_at']} by activity since the last watch: {_clean(r['title'], 60)}")
    for k, (v, at) in sorted(Y.items()):
        if not k.startswith('repo:'):
            continue
        repo = k[5:]
        try:
            head = resolve_sha(repo, 'HEAD')
        except (subprocess.CalledProcessError, ValueError):
            head = None
        if head and v and not head.startswith(v[:7]):
            stale.append(f'repo {repo}: HEAD {head[:7]} != imported {v[:7]}')
    as_of['repos'] = sum(1 for k in Y if k.startswith('repo:'))
    return dict(stale=stale, as_of=as_of)


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
        e = R.setdefault(name, dict(shas=[], paths=[], links=[]))
        if kind and re.fullmatch(r'[0-9a-f]{7,40}', ref):
            e['shas'].append(ref)
        if at:
            e['shas'].append(at)
        if kind in ('tree', 'blob') and path:
            p = urllib.parse.unquote(path.strip('/')).rstrip('.,;:')
            if p and '..' not in p.split('/') and re.fullmatch(r'[A-Za-z0-9._/ +-]{1,300}', p):
                e['paths'].append(p)
                e['links'].append((ref, p))          # 10-17: which commit each path was linked at (#607)
    loose = [s for s in SHA_MENTION.findall(text or '')]
    return R, loose


def propose(issue):
    """One issue -> (entry dict for pending.yaml or None, note)."""
    body, title = issue.get('body') or '', issue.get('title') or ''
    R, loose = refs(title + '\n' + body)
    R = collections.OrderedDict((k, v) for k, v in R.items() if k not in OWN or issue['user']['login'] == 'evand')
    claims = sources.text_claims(title, body)        # body = issue body + comments (run() joins them)
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
        # 10-17: no narrowing.  A tree link imports the whole dir (more parents); sources.reconcile marks which of the
        # packings the issue actually claims ('claim') and which just sat in the linked files ('context')
        # 10-17 (#607): paths linked at another commit of the same repo get their own entry at that commit (one pin
        # per issue put 23 certificate paths at a commit that lacks them)
        by_ref = collections.OrderedDict()
        for ref_, p_ in info['links']:
            by_ref.setdefault(ref_, []).append(p_)
        mine, extra = [], []
        for ref_, ps in by_ref.items():
            sha = pin if pin.startswith(ref_) or ref_ == pin else resolve_sha(repo, ref_)
            if sha is None or sha == pin:
                mine += ps
            else:
                extra.append(dict(issue=e['issue'], finder=e['finder'], repo=repo, commit=sha, paths=sorted(set(ps))))
        e.update(repo=repo, commit=pin, paths=sorted(set(mine)) or (['.'] if not extra else []))
        if extra:
            note.append(f'{len(extra)} more commit(s) linked: {", ".join(x["commit"][:7] for x in extra)}')
            e['_extra'] = extra
    else:
        e.update(repo=None)
    if claims:
        e['claims'] = {int(n): float(s) for n, s in sorted(claims.items())}
    if att:
        e['attachments'] = att
    if not R and not claims and not att:
        return None, 'no repo, attachment or claim'
    return e, '; '.join(note)


def supplement(issue, entries):
    """10-17: a known issue was updated.  -> (supplemental pending entry or None, note): attachments and claims in the
    issue (body + comments) that none of its pending.yaml entries has yet.  New repos are watched by run()."""
    body, title = issue.get('body') or '', issue.get('title') or ''
    have_att = {a for x in entries for a in (x.get('attachments') or [])}
    have_cl = {}
    for x in entries:
        for n_, s_ in (x.get('claims') or {}).items():
            have_cl[int(n_)] = min(have_cl.get(int(n_), 9e9), float(s_))
    att = sorted(set(m.group(0) for m in ATTACH_RE.finditer(body)) - have_att)
    from .sources import text_claims
    cl = {n_: s_ for n_, s_ in text_claims(title, body).items() if s_ < have_cl.get(n_, 9e9) - 1e-12}
    if not att and not cl:
        return None, 'no new attachments or claims'
    e = dict(issue=issue['number'], finder=issue['user']['login'], repo=None)
    if cl:
        e['claims'] = {int(n_): float(s_) for n_, s_ in sorted(cl.items())}
    if att:
        e['attachments'] = att
    return e, f'supplement: {len(att)} new attachments, {len(cl)} new / lower claims (issue updated)'


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
def run(st=None, apply=False, max_commits=MAX_COMMITS):
    S = load_state()
    cfg = pending_cfg()
    known = {e['issue'] for e in cfg if e.get('issue')}
    if S['last_issue'] is None:                     # first run: everything up to the newest pending.yaml issue was seen
        S['last_issue'] = max(known)
    started = now()
    started_local = dt.datetime.now().isoformat(timespec='seconds')     # the store's sighting.added clock (local time)
    log = dict(start=started, start_local=started_local, apply=apply, issues=[], skipped_issues=[], prs=[], repos=[], counts={})

    # (a) issues: every issue updated since the cursor (10-17: comments too, and known issues get supplemental entries
    #     for attachments / claims that appeared later; every repo any of them links is watched)
    cursor = S.get('issues_updated') or S.get('last_poll')
    upd = updated_since(cursor)
    prs = [r for r in upd if r.get('pull_request') and r['user']['login'] != UPSTREAM.split('/')[0]]
    issues = [r for r in upd if not r.get('pull_request')]
    props, extra_repos, touched = [], {}, []
    for r in sorted(issues, key=lambda r: r['number']):
        tag = f"#{r['number']} {r['user']['login']} ({r['state']}) {_clean(r['title'], 70)}"
        try:
            cm = comments(r['number']) if r.get('comments') else []
            own = [b for who, b in cm if who == r['user']['login']]       # claims: the author's comments only
            full = dict(r, body='\n\n'.join([r.get('body') or '', *own]))
            if r['number'] not in known:
                e, note = propose(full)
            else:
                e, note = supplement(full, [x for x in cfg if x.get('issue') == r['number']])
            Rr, _ = refs((r.get('title') or '') + '\n' + '\n\n'.join([r.get('body') or '', *(b for _, b in cm)]))
            for repo in Rr:
                if repo not in OWN:
                    extra_repos.setdefault(repo, r['user']['login'])
        except (subprocess.CalledProcessError, ValueError, KeyError) as ex:
            e, note = None, f'error {_clean(repr(ex), 80)}'
        if e is None:
            print(f'issue {tag}: {"known; nothing new" if r["number"] in known and not note.startswith("error") else "skip"} ({note})')
            if r['number'] not in known:
                log['skipped_issues'].append(dict(issue=r['number'], why=note))
        else:
            extra = e.pop('_extra', [])
            for x in [e, *extra]:
                print(f'issue {tag}: propose {json.dumps(x)[:300]}' + (f'  [{note}]' if note else ''))
                props.append((x, note, r['title']))
                log['issues'].append(dict(entry=x, note=note))
    for r in prs:
        print(f"PR #{r['number']} {r['user']['login']} {_clean(r['title'], 70)}: not handled (PRs are not watched)")
        log['prs'].append(r['number'])

    # (b) repos
    W = watch_list(S, cfg)
    for e, _, _ in props:
        if e.get('repo') and e['repo'] not in OWN:
            W.setdefault(e['repo'], dict(base=e['commit'], finder=e['finder']))
    for repo, who in extra_repos.items():                 # 10-17: every repo linked anywhere in an updated issue
        if REPO_NAME.fullmatch(repo):
            W.setdefault(repo, dict(base=None, finder=who))
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
    if props or touched:
        append_pending(props)
        nums = sorted({e['issue'] for e, _, _ in props} | set(touched))
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
    if upd:
        S['issues_updated'] = max(r['updated_at'] for r in upd)
    S['issues_updated'] = S.get('issues_updated') or started
    st.set_sync('watch.issues_updated', S['issues_updated'], commit=False)
    st.set_sync('watch.poll', started, commit=False)
    for repo, v in S['repos'].items():
        if v.get('imported'):
            st.set_sync(f'repo:{repo}', v['imported'], commit=False)
    st.db.commit()
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
