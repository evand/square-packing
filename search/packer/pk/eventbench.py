#!/usr/bin/env python3
"""Event benchmark harness (PLAN-2026-10-12 item 0; built 10-12).

A benchmark = (cutoff T, operator, budget).  Inputs = everything public at T.  For each n:
  frontier_at(n, T) = min(base(n), every event 'new' side dated < T), where base = the frozen 10-05 register
  (../exact/batch/inputs/n-<n>.txt) for n <= 324; for n > 324 the first event's 'old' (grid, or Couzo's 375 from 10-08),
  as in EVENTS-2026-10-12.md; for n > 324 without events the store's frontier (flagged 'store').
  frontier_now(n) = frontier_at(n, +infinity) (events through 10-13 00:42 UTC).
Success at n = a certified packing with side <= frontier_at(n, T) - 1e-8 (candidate.py).  Secondary = recovered fraction
(frontier(T) - ours) / (frontier(T) - frontier_now), capped at 1.  At an n with no event after T, success = a live record
candidate ('live'): check the live register (candidate.py does) before anyone says "record".

  eventbench.py build-events [--md EVENTS-2026-10-12.md]   parse the events table -> bench/events.json
  eventbench.py cutoffs                                     batch cutoffs, n moved in the next window / any time after
  eventbench.py frontier N [--cutoff T]
  eventbench.py starts --cutoff T --out DIR                 as-of-T start set (n-<n>.txt, chain.py --src), missing.json
  eventbench.py score --cutoff T RUNDIR... [--all] [--certify] [--json OUT]
      RUNDIR = a chain.py run (state.json) or a directory of packing files (n-<n>.txt, best_<n>.txt, <n>.txt, any
      format pk.packing.read knows).  Our best per n = min over all RUNDIRs.

Cutoff T: a name from `cutoffs` (reg1005, 1008, 1009, 1011, 1012, ...) or an ISO time (UTC), e.g. 2026-10-09T20:00.
The store is opened read-only.  Nothing here writes to it.
"""
import argparse, glob, json, math, os, re, sqlite3, subprocess, sys
from datetime import datetime, timezone
HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)                                   # search/packer
sys.path.insert(0, PK)
from pk.packing import Packing, read                          # noqa: E402

REG = os.path.join(PK, '..', 'exact', 'batch', 'inputs')
EVENTS_MD = os.path.join(PK, 'EVENTS-2026-10-12.md')
EVENTS_JSON = os.path.join(PK, 'bench', 'events.json')
STORE = os.environ.get('PK_STORE', os.path.join(PK, 'runs', 'store.sqlite'))
THR = 1e-8            # success threshold (gain)
MATCH = 1e-7          # store row matches an event's side
YEAR = 2026

# batch cutoffs from EVENTS "Batch cutoffs" (bold rows = PLAN's five); window end = the next window's end in that table
CUTOFFS = {
    'reg1005': ('2026-10-05T00:00', '2026-10-07T18:00', '10-05 register (frozen exact/batch/inputs)'),
    '1007':    ('2026-10-07T18:00', '2026-10-08T00:00', '10-07 18:00'),
    '1008':    ('2026-10-08T00:00', '2026-10-09T00:00', '10-08 00:00'),
    '1009':    ('2026-10-09T20:00', '2026-10-10T12:00', '10-09 20:00'),
    '1010':    ('2026-10-10T12:00', '2026-10-11T05:00', '10-10 12:00'),
    '1011':    ('2026-10-11T05:00', '2026-10-12T04:00', '10-11 05:00'),
    '1012':    ('2026-10-12T04:00', '2026-10-13T01:00', '10-12 04:00'),
}
MAIN_CUTOFFS = ('reg1005', '1008', '1009', '1011', '1012')


def ts(s):
    """'2026-10-09T20:00' (UTC, optional Z / seconds) -> aware datetime"""
    s = s.rstrip('Z')
    return datetime.fromisoformat(s).replace(tzinfo=timezone.utc)


def cutoff(T):
    return ts(CUTOFFS[T][0]) if T in CUTOFFS else ts(T)


# ------------------------------------------------------------------ events.json from the markdown table

def build_events(md=EVENTS_MD, out=EVENTS_JSON):
    L = open(md).read().splitlines()
    ev, prec, stale, sec = [], [], [], None
    for l in L:
        if l.startswith('## Events'):
            sec = 'ev'
        elif l.startswith('### Not counted'):
            sec = 'not'
        elif l.startswith('## ') and sec:
            sec = None
        if sec == 'ev' and l.startswith('| ') and re.match(r'\| \d+ \|', l):
            c = [x.strip() for x in l.strip('|').split('|')]
            n, date, author, issue, oldnew, gain, fam, par, pub, arr, ours = c
            m = re.match(r'(\d\d)-(\d\d) (\d\d):(\d\d)$', date)
            o, nw = (float(x) for x in oldnew.split('->'))
            ev.append(dict(n=int(n), date=f'{YEAR}-{m[1]}-{m[2]}T{m[3]}:{m[4]}Z', author=author, issue=issue, old=o, new=nw,
                           gain=float(gain), family=fam, parents=par, parent_public=pub, arr=arr, ours=ours == 'Y',
                           unposted='unposted' in issue))
        if sec == 'not' and l.startswith('- **'):
            (prec if 'Precision' in l else stale).append(l)
        elif sec == 'not' and l.startswith('  ') and (prec or stale):
            (stale if stale else prec)[-1] += ' ' + l.strip()
    clean = lambda s: re.sub(r'^- \*\*[^*]+\*\*\s*', '', s)
    D = dict(source=os.path.relpath(md, PK), built=datetime.now(timezone.utc).isoformat(timespec='seconds'),
             note='rebuilt from the EVENTS markdown table by pk/eventbench.py build-events; dates UTC (issue post time); '
                  'unposted=true marks our big1 375 (counted in the running frontier as the table does)',
             events=ev, precision=[clean(x) for x in prec], stale=[clean(x) for x in stale])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(D, open(out, 'w'), indent=1)
    return D


_EV = None
def events():
    global _EV
    if _EV is None:
        _EV = json.load(open(EVENTS_JSON))['events']
        for e in _EV:
            e['t'] = ts(e['date'])
    return _EV


# ------------------------------------------------------------------ store (read-only)

_DB = None
def db():
    global _DB
    if _DB is None:
        _DB = sqlite3.connect(f'file:{STORE}?mode=ro', uri=True, timeout=60)
        _DB.row_factory = sqlite3.Row
    return _DB


def store_match(n, s, tol=MATCH):
    """the feasible store packing at n whose side is nearest s, if within tol"""
    r = db().execute("SELECT id, s, coords, s_str FROM packing WHERE n=? AND s BETWEEN ? AND ? AND status!='infeasible' "
                     "ORDER BY abs(s-?) LIMIT 1", (n, s - tol, s + tol, s)).fetchone()
    if r is None:
        return None
    import numpy as np
    return Packing(r['s'], np.frombuffer(r['coords'], '<f8').reshape(-1, 3), r['s_str'], meta=dict(id=r['id']))


PUBLIC = False      # set_public(True): only what others can see (posted packings / claims); ours unposted left out
_PUB_SRC = "(g.source LIKE 'register%' OR g.source LIKE 'pending:%' OR g.source LIKE 'kingbird%' OR g.source LIKE 'pool:%')"


EXCLUDE = ()        # sighting sources left out in public mode (e.g. 'pending:#584' to replay a pass before we posted)


def set_public(on=True, exclude=()):
    """Public frontier mode: store rows count only with a public sighting, frontier rows only if not ours, events
    only if posted.  For live checks and for replaying a pass as if our unposted results did not exist."""
    global PUBLIC, EXCLUDE
    PUBLIC, EXCLUDE = on, tuple(exclude)
    _LC.clear(); _BASE.clear()


def store_frontier(n):
    """min of the frontier table and every feasible store packing at n (catalogue rows such as Goebel's s(331) are not
    in the frontier table).  Not as-of-T: rows added after T can only make success harder."""
    if PUBLIC:
        ex = ''.join(f" AND source != '{e}'" for e in EXCLUDE)
        gx = ''.join(f" AND g.source != '{e}'" for e in EXCLUDE)
        r = db().execute(f"SELECT min(s) FROM frontier WHERE n=? AND source NOT LIKE 'ours%'{ex}", (n,)).fetchone()
        q = db().execute("SELECT min(p.s) FROM packing p WHERE p.n=? AND p.status!='infeasible' AND (p.pen IS NULL OR "
                         f"p.pen<=1e-9) AND EXISTS (SELECT 1 FROM sighting g WHERE g.packing=p.id AND {_PUB_SRC}{gx})",
                         (n,)).fetchone()
    else:
        r = db().execute('SELECT min(s) FROM frontier WHERE n=?', (n,)).fetchone()
        q = db().execute("SELECT min(s) FROM packing WHERE n=? AND status!='infeasible' AND (pen IS NULL OR pen<=1e-9)", (n,)).fetchone()
    v = [x for x in (r[0], q[0]) if x is not None]
    return min(v) if v else None


# ------------------------------------------------------------------ frontier

_BASE = {}
def base(n):
    """(side, kind): 10-05 register for n <= 324; n > 324: first event's old, else the store frontier"""
    if n not in _BASE:
        p = f'{REG}/n-{n}.txt'
        if os.path.exists(p):
            _BASE[n] = (float(open(p).readline().split()[1]), 'reg1005')
        else:
            E = sorted((e for e in events() if e['n'] == n), key=lambda e: e['t'])
            if E:
                _BASE[n] = (E[0]['old'], 'first-event-old')
            else:
                _BASE[n] = (store_frontier(n) or float(math.ceil(math.sqrt(n) - 1e-12)), 'store')
    return _BASE[n]


ENVELOPE = 80        # s(n) <= s(m) for m > n (drop squares): the frontier is the envelope over m in [n, n + ENVELOPE]


def _point_at(n, T):
    t = cutoff(T) if isinstance(T, str) else T
    return min([base(n)[0]] + [e['new'] for e in events() if e['n'] == n and (t is None or e['t'] < t)
                               and not (PUBLIC and e.get('unposted'))])


def frontier_at(n, T):
    """Monotone envelope (10-12: without it, 325-330 'beat' the grid while Goebel's s(331) already gives 18.7071)."""
    return min(_point_at(m, T) for m in range(n, n + ENVELOPE + 1))


_LC = {}
def _l_closure(nmax=1200):
    """Known values closed under Friedman's L-step (DS7 sec. 2): n squares in side s (s not an integer) ->
    n + 2 floor(s) + 1 squares in side s + 1 (a 1-wide L along two edges).  Missing it made 26 of 36 of our 10-13
    n > 324 "records" stale (s(290) + L = 18.5355 at 325 beat our 130/7).  Inputs: events envelope + store."""
    if not _LC:
        v = {}
        for n in range(2, nmax + 1):
            x = [frontier_at(n, None) if n <= nmax - ENVELOPE else math.inf, store_frontier(n) or math.inf]
            v[n] = min(x + [float(math.ceil(math.sqrt(n) - 1e-12))])
        for n in range(2, nmax + 1):                     # ascending: chains propagate
            s = v[n]
            if abs(s - round(s)) > 1e-12:
                m = n + 2 * math.floor(s) + 1
                if m <= nmax and s + 1 < v[m]:
                    v[m] = s + 1
        for n in range(nmax - 1, 1, -1):                 # monotone envelope s(n) <= s(n + 1)
            v[n] = min(v[n], v[n + 1])
        _LC.update(v)
    return _LC


def frontier_now(n):
    """Live frontier: events + store (frontier table incl. pending claims, every feasible packing, posts newer than
    the events table), monotone envelope, closed under the L-step."""
    return _l_closure()[n]


def moved_ns(t0, t1=None):
    return sorted({e['n'] for e in events() if e['t'] >= t0 and (t1 is None or e['t'] < t1)})


# ------------------------------------------------------------------ starts

def grid_packing(n):
    k = math.ceil(math.sqrt(n) - 1e-12)
    return Packing(float(k), [(i % k + 0.5, i // k + 0.5, 0.0) for i in range(n)], str(k))


def starts(T, out):
    t = cutoff(T)
    os.makedirs(out, exist_ok=True)
    missing, src = {}, {}
    ns = sorted({int(re.search(r'n-(\d+)', p)[1]) for p in glob.glob(f'{REG}/n-*.txt')} | {e['n'] for e in events()})
    for n in ns:
        E = sorted((e for e in events() if e['n'] == n and e['t'] < t), key=lambda e: e['new'])
        got = None
        for i, e in enumerate(E):                          # best public side first, fall back to older events
            p = store_match(n, e['new'])
            if p is not None:
                got = (p, f'store id {p.meta["id"]} = {e["author"]} {e["issue"]} ({e["date"]})')
                if i:
                    missing[n] = dict(want=E[0]['new'], want_by=f'{E[0]["author"]} {E[0]["issue"]}', have=p.s,
                                      gap=p.s - E[0]['new'])
                break
        if got is None:
            if os.path.exists(f'{REG}/n-{n}.txt'):
                p = read(f'{REG}/n-{n}.txt'); got = (p, 'reg1005')
            else:
                b = base(n)[0]
                p = store_match(n, b) if b != round(b) else grid_packing(n)
                if p is None:
                    p = grid_packing(n)
                got = (p, 'grid' if p.s == round(p.s) else f'store id {p.meta["id"]} (base)')
            if E:
                missing[n] = dict(want=E[0]['new'], want_by=f'{E[0]["author"]} {E[0]["issue"]}', have=got[0].s,
                                  gap=got[0].s - E[0]['new'])
        p, why = got
        p.write(f'{out}/n-{n}.txt')
        src[n] = dict(s=p.s, frontier=frontier_at(n, t), src=why)
    json.dump(dict(cutoff=T, t=t.isoformat(), missing=missing), open(f'{out}/missing.json', 'w'), indent=1)
    json.dump(src, open(f'{out}/starts.json', 'w'), indent=1)
    return src, missing


# ------------------------------------------------------------------ score

def our_best(dirs):
    """{n: (side, path)} over chain run dirs (state.json) and dirs of packing files"""
    B = {}
    def put(n, s, p):
        if n not in B or s < B[n][0]:
            B[n] = (s, p)
    for d in dirs:
        if os.path.exists(f'{d}/state.json'):
            S = json.load(open(f'{d}/state.json'))['n']
            for n, v in S.items():
                p = v.get('path')
                p = p if p and os.path.isabs(p) else os.path.join(d, os.path.basename(p or f'best_{n}.txt'))
                put(int(n), v['best'], p)
        else:
            for p in glob.glob(f'{d}/*'):
                m = re.match(r'(?:n-|best_|n)?(\d+)\.(?:txt|cert|json)$', os.path.basename(p))
                if not m:
                    continue
                q = read(p)
                if q is not None and q.n == int(m[1]):
                    put(q.n, float(q.s), p)
    return B


def certify(n, path, tag):
    """candidate.py on one packing (no regularize); verdict + certified side"""
    name = f'eventbench_{tag}_n{n}'
    r = subprocess.run([sys.executable, os.path.join(PK, 'candidate.py'), path, '--name', name, '--no-regularize'],
                       capture_output=True, text=True, cwd=PK, timeout=7200)
    sm = os.path.join(PK, 'runs', 'candidates', name, 'summary.json')
    cj = os.path.join(PK, 'runs', 'candidates', name, 'cand.json')
    if r.returncode == 0 and os.path.exists(sm):
        S = json.load(open(sm))
        return dict(verdict='certified', S_cert=float(S['S_cert']), live_beaten=S['live'].get('beaten'), work=name)
    msg = (r.stderr.strip().splitlines() or r.stdout.strip().splitlines() or ['?'])[-1]
    v = dict(verdict='FAILED', why=msg[:200], work=name)
    if os.path.exists(cj):
        J = json.load(open(cj)); v['S_cert'] = float(J['S_cert_decimal']) if J.get('S_cert_decimal') else None
    return v


def score(T, dirs, all_n=False, do_cert=False, tag=None):
    t = cutoff(T)
    B = our_best(dirs)
    rows = []
    for n in sorted(B):
        s, p = B[n]
        fT, fN = frontier_at(n, t), frontier_now(n)
        g = fT - s
        live = fN >= fT - 1e-15
        later = fT - fN
        r = dict(n=n, frontier_T=fT, ours=s, gain=g, success=g >= THR, live=live, later_gain=later,
                 recovered=None if live else min(1.0, max(0.0, g / later)),
                 full=(not live) and s <= fN + 1e-9, below_now=s < fN - 1e-9, path=p, base=base(n)[1])
        if r['success'] and do_cert:
            r['cert'] = certify(n, p, tag or (T if T in CUTOFFS else 'T'))
            sc = r['cert'].get('S_cert')
            r['certified'] = r['cert']['verdict'] == 'certified' and sc is not None and sc <= fT - THR
        if all_n or g > 1e-12 or not live:
            rows.append(r)
    return rows


def print_score(rows, T):
    print(f'cutoff {T} ({cutoff(T).isoformat()}); success = gain >= {THR:g} (uncertified unless --certify); '
          f'live = nobody improved n after T')
    print(f'{"n":>4} {"frontier(T)":>15} {"ours":>15} {"gain":>9} {"succ":>11} {"recov":>5} {"later":>8}  note')
    for r in rows:
        su = ('cert' if r.get('certified') else 'CERT-FAIL' if 'cert' in r else 'uncert') if r['success'] else '-'
        note = []
        if r['live']:
            note.append('LIVE (no later event)' + (': check live register' if r['success'] else ''))
        if r['below_now']:
            note.append('below frontier_now')
        if r['full']:
            note.append('full')
        if 'cert' in r:
            note.append(f"cert {r['cert']['verdict']} S'={r['cert'].get('S_cert')} live_beaten={r['cert'].get('live_beaten')}")
        rec = '  n/a' if r['recovered'] is None else f'{r["recovered"]:5.2f}'
        print(f'{r["n"]:4d} {r["frontier_T"]:15.10f} {r["ours"]:15.10f} {r["gain"]:9.2e} {su:>11} {rec} {r["later_gain"]:8.1e}  '
              + '; '.join(note))
    ev = [r for r in rows if not r['live']]
    print(f'n with later events: {len(ev)}; success {sum(r["success"] for r in ev)} {[r["n"] for r in ev if r["success"]]}; '
          f'full recovery {sum(r["full"] for r in ev)}; mean recovered {sum(r["recovered"] for r in ev) / max(1, len(ev)):.3f}; '
          f'moved (gain > 0) {sum(r["gain"] > 1e-12 for r in ev)}')
    miss = [n for n in moved_ns(cutoff(T)) if n not in {r['n'] for r in rows}]
    print(f'n improved after T that these runs did not produce a packing for ({len(miss)}): {miss}')
    lv = [r for r in rows if r['live'] and r['success']]
    print(f'live-n successes (record candidates, check the live register): {[r["n"] for r in lv]}')


# ------------------------------------------------------------------ CLI

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('build-events'); p.add_argument('--md', default=EVENTS_MD); p.add_argument('--out', default=EVENTS_JSON)
    sp.add_parser('cutoffs')
    p = sp.add_parser('frontier'); p.add_argument('n', type=int); p.add_argument('--cutoff')
    p = sp.add_parser('starts'); p.add_argument('--cutoff', required=True); p.add_argument('--out', required=True)
    p = sp.add_parser('score'); p.add_argument('--cutoff', required=True); p.add_argument('dirs', nargs='+')
    p.add_argument('--all', action='store_true', help='also list n where ours did not move and nobody improved later')
    p.add_argument('--certify', action='store_true', help='candidate.py on each success (exactsolve + checkers + live check)')
    p.add_argument('--tag', help='candidate work-dir tag (runs/candidates/eventbench_<tag>_n<n>)')
    p.add_argument('--json')
    a = ap.parse_args()
    if a.cmd == 'build-events':
        D = build_events(a.md, a.out)
        print(f'{len(D["events"])} events, {len({e["n"] for e in D["events"]})} n -> {a.out}')
    elif a.cmd == 'cutoffs':
        print(f'{"name":>8} {"T (UTC)":>17} {"window end":>17} {"next window: n (<=324)":>23} {"after T: n (<=324)":>19}  label')
        for k, (t0, t1, lab) in CUTOFFS.items():
            w, aft = moved_ns(ts(t0), ts(t1)), moved_ns(ts(t0))
            print(f'{k + ("*" if k in MAIN_CUTOFFS else ""):>8} {t0:>17} {t1:>17} {len(w):>14} ({sum(n <= 324 for n in w):>3})'
                  f' {len(aft):>12} ({sum(n <= 324 for n in aft):>3})  {lab}')
        print('* = PLAN batch cutoffs')
    elif a.cmd == 'frontier':
        print(f'base {base(a.n)}; frontier_at({a.cutoff}) = {frontier_at(a.n, a.cutoff) if a.cutoff else "-"}; '
              f'frontier_now = {frontier_now(a.n)}; store frontier = {store_frontier(a.n)}')
    elif a.cmd == 'starts':
        src, miss = starts(a.cutoff, a.out)
        ch = {n: v for n, v in src.items() if v['src'] != 'reg1005'}
        print(f'{len(src)} start files -> {a.out}; {len(ch)} differ from the 10-05 register; {len(miss)} missing from the store '
              f'(start weaker than public): {sorted(miss)}')
        bad = [n for n, v in src.items() if v['s'] > v['frontier'] + MATCH]
        if bad:
            print(f'start above frontier_at(T) by > {MATCH:g}: {bad}')
    elif a.cmd == 'score':
        rows = score(a.cutoff, a.dirs, a.all, a.certify, a.tag)
        print_score(rows, a.cutoff)
        if a.json:
            json.dump(rows, open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    main()
