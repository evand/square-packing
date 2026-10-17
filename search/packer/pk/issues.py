"""Issue text -> what to import (10-17).  Pure functions: no network, no store, no files.  Issue text is untrusted
data; only these regexes read it.

  parse(issue, comments)          -> Parsed: links per repo and commit, attachments, the author's claims
  units(parsed, resolve, override) -> import units [(repo, commit sha, paths)] + notes

Re-run on every update of an issue (pk/watch.py), so nothing parsed is stored by hand and nothing goes stale.  What a
human knows that the text doesn't say lives in pk/overrides.yaml (see `units`).
"""
from __future__ import annotations
import collections, math, re, urllib.parse

UPSTREAM = 'jlevy/squares'
OWN = {'evand/square-packing'}                      # ours: imported only from our own issues
OWN_LOGINS = {'evand'}
REPO_NAME = re.compile(r'[A-Za-z0-9](?:[A-Za-z0-9-]{0,38})/[A-Za-z0-9._-]{1,100}')
RESERVED = {'user-attachments', 'orgs', 'settings', 'issues', 'pulls', 'marketplace', 'sponsors', 'apps', 'features',
            'topics', 'notifications', 'login', 'about', 'pricing', 'enterprise', 'collections', 'search',
            'advisories', 'site', 'contact', 'explore', 'trending', 'new', 'organizations', 'join', 'security'}
# github.com/OWNER/REPO[/commit|tree|blob/REF[/path]][@sha]: the only link form read from issue text
GH_RE = re.compile(r'(?<![\w./-])(?:https?://)?(?:www\.)?github\.com/([A-Za-z0-9][A-Za-z0-9-]{0,38})/([A-Za-z0-9._-]{1,100})'
                   r'(?:/(commit|tree|blob)/([0-9A-Za-z._-]{1,100})((?:/[^\s/?#)\]>"\'`]{1,200}){0,20}))?'
                   r'(?:@([0-9a-f]{7,40})\b)?')
SHA_MENTION = re.compile(r'(?:commit|pinned at|pin|sha)[\s:`*]*([0-9a-f]{7,40})\b', re.I)
TITLE_CLAIM = re.compile(r's\(\s*(\d{1,5})\s*\)\s*(?:≤|<=|\\le)\s*(\d+\.\d{6,})')
# Issue attachments: only https://github.com/user-attachments/files/<id>/<name> (fetched by github.fetch_attachment)
ATTACH_RE = re.compile(r'https://github\.com/user-attachments/files/(\d{1,12})/([A-Za-z0-9._-]{1,120})')
HEX = re.compile(r'[0-9a-f]{7,40}')
# 10-17: tables are read as upper-bound claims only in claim issues (the 669-issue backfill: the maintainer's working
# issues, lower-bound and proof issues carry tables of sides that are not claims).  Files are imported regardless.
CLAIM_TITLE = re.compile(r'registration request|s\(\s*\d+\s*\)\s*(?:≤|<=|\\le)|upper bound|below (?:the )?(?:record|register)|'
                         r'new (?:best|record|packing|upper|results?)|improved|certified s\s*=|certificate refinement', re.I)
NOT_CLAIM_TITLE = re.compile(r'lower bound|optimality|lean proof|exact value', re.I)
MAINTAINERS = {UPSTREAM.split('/')[0]}


def claim_issue(author, title):
    """does this issue claim upper bounds (so its tables are claims)?"""
    t = title or ''
    return (author not in MAINTAINERS and bool(CLAIM_TITLE.search(t))
            and not (NOT_CLAIM_TITLE.search(t) and not re.search('upper', t, re.I)))


def refs(text):
    """({repo: {'shas': [named commits], 'links': [(ref, path)]}}, [loose sha mentions]) for github.com/OWNER/REPO
    links in text, in order of appearance (upstream and reserved names excluded).  A link with no path (repo, commit,
    @sha) adds no (ref, path)."""
    R = collections.OrderedDict()
    for m in GH_RE.finditer(text or ''):
        owner, repo, kind, ref, path, at = m.groups()
        repo = re.sub(r'(\.git)?\.*$', '', repo)
        name = f'{owner}/{repo}'
        if owner.lower() in RESERVED or name.lower() == UPSTREAM or not repo or not REPO_NAME.fullmatch(name):
            continue
        e = R.setdefault(name, dict(shas=[], links=[]))
        if kind and HEX.fullmatch(ref):
            e['shas'].append(ref)
        if at:
            e['shas'].append(at)
        if kind in ('tree', 'blob') and path:
            p = urllib.parse.unquote(path.strip('/')).rstrip('.,;:')
            if p and '..' not in p.split('/') and re.fullmatch(r'[A-Za-z0-9._/ +-]{1,300}', p):
                e['links'].append((ref, p))
    return R, SHA_MENTION.findall(text or '')


def body_claims(body):
    """{n: claimed side} from markdown table rows (with or without outer pipes; the count may be a link, `[237](url)`).
    The claim is the smallest cell that can be a side for n (sqrt(n) <= v <= sqrt(n) + 3), not the first decimal
    column: #668's table reads "n | website-listed side | certified upper bound | reduction"."""
    claims = {}
    for line in (body or '').split('\n'):
        cells = [x.strip(' *`$') for x in line.strip().strip('|').split('|')]
        m0 = re.fullmatch(r'\[?(\d{2,4})\]?(?:\([^)]*\))?', cells[0]) if len(cells) >= 2 else None
        if not m0:
            continue
        n = int(m0.group(1))
        v = [float(m.group(1)) for c in cells[1:] for m in [re.match(r'(\d+\.\d{6,})', c)] if m]
        v = [x for x in v if math.sqrt(n) - 1e-9 <= x <= math.sqrt(n) + 3]
        if v:
            claims[n] = min(claims.get(n, 9e9), min(v))
    return claims


def text_claims(title, body, comments=()):
    """{n: smallest claimed side}: tables in the body and the given comments, plus `s(n) <= x` in the title"""
    claims = {}
    for txt in [body, *comments]:
        for n, s in body_claims(txt).items():
            claims[n] = min(claims.get(n, 9e9), s)
    for n, s in TITLE_CLAIM.findall(title or ''):
        claims[int(n)] = min(claims.get(int(n), 9e9), float(s))
    return claims


def parse(issue, comments=()):
    """issue = GitHub issue json; comments = [(login, body)].  Links and attachments from all text (anyone's comments
    may point at files); claims only from the author's own text (the maintainer's verdict tables restate other
    issues' values, #503 / #573)."""
    who = issue['user']['login']
    title, body = issue.get('title') or '', issue.get('body') or ''
    alltext = '\n\n'.join([title, body, *(b for _, b in comments)])
    R, loose = refs(alltext)
    return dict(number=issue['number'], author=who, state=issue['state'], title=title,
                updated=issue.get('updated_at'), created=issue.get('created_at'), pr=bool(issue.get('pull_request')),
                repos=R, loose=loose,
                attachments=sorted(set(m.group(0) for m in ATTACH_RE.finditer(alltext))),
                claim_issue=claim_issue(who, title),
                claims=text_claims(title, body, [b for w, b in comments if w == who]) if claim_issue(who, title) else {})


def units(P, resolve, override=None):
    """Parsed issue -> (units, attachments, claims, notes).  units = [dict(repo, commit, paths)]: one per (repo,
    commit) the issue links; paths linked at a commit that doesn't resolve go with the repo's pin.  A repo linked
    without paths is imported whole at its pin.  Pin = the first named commit (link or mention) that resolves, else
    the default-branch HEAD at first sight (resolve is expected to cache it per issue, so it stays put).  A repo linked
    bare (no path) is imported whole only from a claim issue (10-17: the maintainer's issues link big code repos).
    resolve(repo, ref, issue) -> full sha or None.
    override (pk/overrides.yaml, one issue): skip (import nothing), claims (merged by min: prose claims), units
    (added), only_override_units (drop the derived ones), finder (display name)."""
    o = override or {}
    notes = []
    claims = dict(P['claims'])
    for n, s in (o.get('claims') or {}).items():
        claims[int(n)] = min(claims.get(int(n), 9e9), float(s))
    if o.get('skip'):
        return [], [], claims, [f"skipped by override: {o.get('note', '')}"]
    U = []
    if not o.get('only_override_units'):
        for repo, info in P['repos'].items():
            if repo in OWN and P['author'] not in OWN_LOGINS:
                continue
            by, loose_paths = collections.OrderedDict(), []
            for ref, p in info['links']:
                sha = resolve(repo, ref, P['number'])
                if sha:
                    by.setdefault(sha, []).append(p)
                else:
                    loose_paths.append(p)
            if not info['links'] and not P.get('claim_issue', True):
                notes.append(f'{repo}: linked bare in a non-claim issue: not imported'); continue
            if loose_paths or not info['links']:          # paths at a dead ref, or the repo linked bare: need a pin
                pin = None
                for c in info['shas'] + P['loose']:
                    pin = resolve(repo, c, P['number'])
                    if pin:
                        break
                if pin is None:
                    pin = resolve(repo, 'HEAD', P['number'])
                    notes.append(f'{repo}: pin = HEAD at first sight ' +
                                 ('(named commit not found)' if info['shas'] or loose_paths else '(issue names no commit)'))
                if pin is None:
                    notes.append(f'{repo}: not reachable'); continue
                by.setdefault(pin, []).extend(loose_paths or ['.'])
            U += [dict(repo=repo, commit=sha, paths=sorted(set(ps))) for sha, ps in by.items()]
    for u in o.get('units') or ():
        U.append(dict(repo=u['repo'], commit=u['commit'], paths=list(u.get('paths') or ['.'])))
    return U, P['attachments'], claims, notes
