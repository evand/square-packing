#!/usr/bin/env python3
"""Registration-request draft from a candidate.py work dir (10-20; style guide and template: SUBMITTING.md).

Fills in everything that is a number or a fact: the claim table (live check), the certificate files at their commit
with SHA-256, checker results, local optimality, the relationship to the prior record (candidate.py's vs_prior) and to
every credited start (--start: polish the start, compare with ours: same minimum = refinement, else a new packing found
from it, with jlevy's matched-move numbers), and our open questions on earlier issues.  The lineage story is left as
TODO.  Nothing is posted (CLAUDE.md: Evan's OK per post).

  python3 submission.py NAME [--start STORE_ID ...] [--live [--sync]] [--questions] [--out PATH] [--force]

NAME = the candidate name (runs/candidates/NAME, candidates/NAME.*).  Output: outreach/draft-jlevy-n<N>-<date>.md.
"""
from __future__ import annotations
import argparse, datetime, hashlib, json, os, re, subprocess, sys
from decimal import Decimal, ROUND_CEILING, getcontext
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
CRED = 'Evan Daniel (GitHub @evand)'
OURS = re.compile(r'evand|Evan Daniel|\bours\b', re.I)
getcontext().prec = 90


def sh(*a, cwd=ROOT):
    r = subprocess.run(a, capture_output=True, text=True, cwd=cwd)
    return r.stdout.strip() if r.returncode == 0 else ''


def up12(x):
    return str(Decimal(str(x)).quantize(Decimal('1e-12'), rounding=ROUND_CEILING))


def dec(q):
    """a/b -> 30-place decimal string (the regularizer reports fractions)"""
    if '/' not in str(q):
        return q
    a, b = str(q).split('/')
    return str((Decimal(a) / Decimal(b)).quantize(Decimal('1e-30')))


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest() if os.path.exists(path) else None


def commit_of(paths):
    """newest commit touching the staged files, and whether origin has it"""
    rel = [os.path.relpath(p, ROOT) for p in paths if os.path.exists(p)]
    c = sh('git', 'log', '-1', '--format=%h', '--', *rel)
    dirty = sh('git', 'status', '--porcelain', '--', *rel)
    pushed = bool(c) and bool(sh('git', 'branch', '-r', '--contains', c))
    return c, bool(dirty), pushed


# ------------------------------------------------------------------ relationship to a credited start

def relation(work, pid):
    """polish store packing pid (fq POLISH, 300 s CPU) and compare with our KKT point.  Same minimum = the polished
    start lands on our side (1e-9) with no load-bearing square moving > 1e-2 under the best matching over the 8 box
    symmetries.  Also the raw start vs ours moves (jlevy's 0.1 rule is geometric)."""
    from pk.store import Store
    from pk.engine import FQ, POLISH
    from pk.packing import read
    from pk import distinct as D
    st = Store(); P = st.get(pid)
    who = sorted({f'{f} ({s})' for s, f in st.db.execute('SELECT source, finder FROM sighting WHERE packing=?', (pid,))
                  if not str(s).startswith(('run:', 'compose', 'harden', 'finish', 'candidate'))})   # public sightings
    O = read(os.path.join(work, 'cand.exact.txt')); Oq, sO = np.array(O.sq, float), float(O.s)
    K = json.load(open(os.path.join(work, 'cand.contacts.json')))
    free = set(K.get('free') or [])
    held = np.array([i not in free for i in range(len(Oq))])

    def mv(Q, sQ):
        md, k, perm, Pm = D.best_match(Oq, sO, np.array(Q, float), sQ)
        m, _ = D.moves(Oq, Pm[perm])
        return dict(max=float(m[held].max()), gt1e1=int((m[held] > 0.1).sum()), gt1e2=int((m[held] > 1e-2).sum()),
                    gt1e1_all=int((m > 0.1).sum()), rms=float(md))
    raw = mv(P.sq, float(P.s))
    q = FQ().quench(P, POLISH, timeout=1800, cpu_limit=300)
    out = dict(pid=pid, n=P.n, s=float(P.s), who=who, raw=raw, s_pol=None, pol=None)
    if q.p is not None:
        out['s_pol'] = float(q.s); out['pol'] = mv(q.p.sq, float(q.s))
    same = (out['s_pol'] is not None and P.n == len(Oq) and abs(out['s_pol'] - sO) < 1e-9 and out['pol']['gt1e2'] == 0)
    out['verdict'] = ('same minimum (a refinement of it)' if same else
                      'different local minimum (a new packing found from it)' if P.n == len(Oq) else
                      f'start at n = {P.n} (a new packing found from it)')
    return out


def rel_line(r):
    who = '; '.join(r['who'][:3]) or f'store {r["pid"]}'
    s = f"* **{who}** (store {r['pid']}, n = {r['n']}, side {r['s']:.12f}): **{r['verdict']}**."
    if r['n'] == r.get('_n'):
        g = r['raw']
        s += (f"  Matched to ours (best of the 8 box symmetries, load-bearing squares): max move {g['max']:.3g}, "
              f"{g['gt1e1']} move > 0.1, {g['gt1e2']} > 0.01")
        if r['pol']:
            s += (f"; after polishing it: side {r['s_pol']:.12f} ({r['s_pol'] - r['_s']:+.2e} vs ours), "
                  f"max move {r['pol']['max']:.3g}, {r['pol']['gt1e1']} > 0.1")
        s += '.'
        if 'different' in r['verdict'] and g['gt1e1'] == 0:
            s += ('  We call this a new packing although no square moves more than 0.1 (your refinement threshold): '
                  'it is a different local minimum, found by TODO(search).')
    return s


# ------------------------------------------------------------------ open questions on our earlier issues

def open_questions():
    """our jlevy/squares issues whose latest comment is theirs and contains a question"""
    js = sh('gh', 'issue', 'list', '-R', 'jlevy/squares', '--author', '@me', '--state', 'all', '--limit', '30',
            '--json', 'number,title,comments')
    out = []
    for i in json.loads(js or '[]'):
        c = i['comments']
        if c and c[-1]['author']['login'] != 'evand' and '?' in c[-1]['body']:
            qs = [s.strip() for s in re.split(r'(?<=[.?])\s+', c[-1]['body']) if s.strip().endswith('?')]
            out.append((i['number'], i['title'][:60], qs[:3]))
    return out


# ------------------------------------------------------------------ draft

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name'); ap.add_argument('--start', type=int, nargs='*', default=[])
    ap.add_argument('--live', action='store_true', help='re-run the live check now (read-only unless --sync)')
    ap.add_argument('--sync', action='store_true'); ap.add_argument('--questions', action='store_true')
    ap.add_argument('--out'); ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    work = os.path.join(HERE, 'runs', 'candidates', a.name)
    S = json.load(open(os.path.join(work, 'summary.json')))
    n = S['n']; Sx = S['S_exact']; Sc = S['S_cert']
    C = os.path.join(HERE, 'candidates', a.name)
    files = [(f'{a.name}.cert', 'the certificate'), (f'{a.name}.exact.txt', 'the 80-digit KKT point'),
             (f'{a.name}.input.txt', 'the f64 packing it was solved from'), (f'{a.name}.json', "the solver's report"),
             (f'{a.name}_regularized.cert', 'the same packing at the same side, regularized for display')]
    files = [(f, d) for f, d in files if os.path.exists(os.path.join(HERE, 'candidates', f))]
    commit, dirty, pushed = commit_of([os.path.join(HERE, 'candidates', f) for f, _ in files])
    todo = []
    if not files:
        todo.append('stage the files (candidate.py --stage)')
    if not commit or dirty:
        todo.append('commit the staged files')
    elif not pushed:
        todo.append(f'push {commit}')
    # live check
    if a.live:
        import candidate
        L = candidate.live(n, float(Sx[:30]), sync=a.sync)
    else:
        L = S['live']
    others = {k: v for k, v in L['sources'].items() if not OURS.search(k)}
    blabel, bval = min(others.items(), key=lambda kv: kv[1]) if others else ('TODO', float('nan'))
    asof = L.get('as_of', {})
    regj = os.path.join(C + '_regularized.json')
    by = json.load(open(regj)).get('by') if os.path.exists(regj) else None
    if by is not None and not any(OURS.search(b) for b in by):
        todo.append(f'regularized view says "by": {by}: re-run candidate.py (10-20 fix) or edit it')
    O = S['optimality']
    # relationships
    vp = C + '_vs_prior.json' if os.path.exists(C + '_vs_prior.json') else os.path.join(work, 'vs_prior.json')
    rel = []
    if os.path.exists(vp):
        from pk import distinct as D
        V = json.load(open(vp)); d = V['difference']
        term = {'distinct minimum': 'a different local minimum: a new packing, not a refinement',
                'same minimum': 'the same local minimum: a refinement', }.get(V['verdict'], V['verdict'] + ' (look closer before choosing a term)')
        rel.append(f"* **The previous record**: **{term}**.  Exact sides differ by {abs(d['side_gap']):.3e}; "
                   f"matched (regularized views, load-bearing squares): max move {d['max_move']:.3g}, "
                   f"{d.get('moved_gt_1e1', '?')} move > 0.1, {d['moved_gt_1e2']} > 0.01; tilt angles shared: "
                   f"{', '.join(d['shared_angles']) or 'none'}.")
        if V['verdict'] == 'distinct minimum' and d.get('moved_gt_1e1', 1) == 0:
            rel.append('  We call this a new packing although no square moves more than 0.1 (your refinement threshold): '
                       'it is a different local minimum, found by TODO(search).')
        rel.append('')
        rel.append(D.table_md(V, 'this packing', 'previous record'))
    else:
        todo.append('no vs_prior comparison (candidate.py without --no-prior)')
    for pid in a.start:
        r = relation(work, pid); r['_n'] = n; r['_s'] = float(Sx[:30])
        rel.append(rel_line(r))
    if not a.start:
        todo.append('name the credited start(s) with --start STORE_ID to get their relationship lines')
    # draft
    date = datetime.date.today().isoformat()
    ver = ', '.join(f"{k} {'VALID' if v['ok'] else 'FAILED'}" for k, v in S['checkers'].items())
    opt = (f"exact side S* = {Sx[:24]}…; {'a KKT point' if O.get('kkt') else 'TODO: not a KKT point?'}; "
           f"{O.get('status')}; corner–corner MILP descent: {'none' if (O.get('vv_milp_dS') or 0) >= 0 else O['vv_milp_dS']}; "
           f"{O.get('flat_motions')} exact flat motions, {O.get('corner_corner')} corner–corner contacts, {O.get('rattlers')} rattlers.")
    L_ = [f'# DRAFT {date} (submission.py; style guide search/packer/SUBMITTING.md).  Not posted: needs Evan\'s OK.',
          f'# Title: s({n}) ≤ {up12(Sc)}: registration request', '',
          *(['**Before posting:** ' + '; '.join(todo) + '.', ''] if todo else []),
          f'**Credit name.** {CRED}.', '',
          '**Claim.** One new upper bound.  The side is the certificate\'s exact rational side, rounded up at 12 places.', '',
          '| n | New upper bound | Best public before it | Improvement | That bound\'s source |',
          '|---:|---:|---:|---:|---|',
          f'| {n} | {up12(Sc)} | {up12(bval)} | {bval - float(Sc[:30]):.2e} | {blabel} |', '',
          f"Comparison = the lowest of the register (main @{asof.get('register', '?')}), open and closed issues through "
          f"#{L.get('latest_issue', '?')}, David Ellsworth's catalogue, and every value those give through the monotone "
          f"envelope.  Checked {date if a.live else 'TODO(date of the candidate run; re-run with --live)'}.", '',
          f"**Certificate.** In https://github.com/evand/square-packing at commit "
          f"[{commit or 'TODO'}](https://github.com/evand/square-packing/tree/{commit or 'TODO'}/search/packer/candidates):",
          ', '.join(f'`{f}` ({d})' for f, d in files) + '.', '',
          '| file | S\' (certificate side) | SHA-256 |', '|---|---|---|',
          *[f"| `{f}` | {Sc if not f.endswith('regularized.cert') else dec(S.get('regularized', {}).get('S_cert', '?'))} | "
            f"{sha(os.path.join(HERE, 'candidates', f))} |" for f, _ in files if f.endswith('.cert')], '',
          '**Credit and relationship.**  TODO(one or two sentences: "Please credit ... as start(s)"; "a new packing '
          '(a different local minimum) found by <search> from starts including X").', '',
          '<details><summary>Checking</summary>', '',
          f'**Checking.**  `search/exact/verify_cert.py` (separating axes) and `search/exact/verify_cert2.py` (polygon '
          f'intersection): {ver}; exact rational, standard library only: '
          f'`python3 search/exact/verify_cert.py search/packer/candidates/{a.name}.cert`.  They share an author (us) and '
          f'the certificate format; the geometry tests differ.', '',
          '</details>', '', '<details><summary>Local optimality</summary>', '',
          f'**Local optimality** (numerical evidence at 80 digits, not interval proofs): {opt}  TODO(explorer restarts '
          f'that found nothing lower).', '', '</details>', '',
          '<details><summary>Relationship to earlier packings: the numbers</summary>', '',
          '*Refinement* here means the same local minimum only; otherwise a new packing found from it.', '', *rel, '',
          '</details>', '', '<details><summary>How it was found</summary>', '',
          '**How it was found.**  TODO(short lineage: the search, and N moves from X; or hunt: method, generations from the credited start, '
          'intermediate non-record packings, CPU-h).  Please credit TODO(start packings, with source) as start(s).',
          '1. TODO', '2. TODO',
          'N. **Polish and certify:** TODO(restarts from it found nothing lower); `exactsolve.py` then certified it.', '',
          '</details>', '']
    if a.questions:
        Q = open_questions()
        if Q:
            L_ += ['<!-- Open questions from jlevy on our earlier issues (answer there, or here if related):']
            L_ += [f'  #{k} {t}: ' + ' | '.join(q) for k, t, q in Q] + ['-->', '']
    out = a.out or os.path.join(ROOT, 'outreach', f'draft-jlevy-n{n}-{date}.md')
    if os.path.exists(out) and not a.force:
        sys.exit(f'{out} exists (--force to overwrite)')
    open(out, 'w').write('\n'.join(L_))
    print(f'wrote {out}' + (f'  (TODO: {"; ".join(todo)})' if todo else ''))


if __name__ == '__main__':
    main()
