"""Frozen reference sets and holdouts for census-style experiments (10-11; CONFUSIONS X03, TODO census v2).

  snapshot  freeze the store's minima of n (table `minimum` below k, incl. weak) into refsets/<name>.json (tracked): per minimum
            key, S, cls, a representative packing id, and its copy side keys; sha256 of the content.  Reports that score
            against a snapshot cannot drift as the live tables grow.
  holdout   split a snapshot into visible / withheld minima:
              random  --frac F --seed K           withhold a random fraction
              side    --lo A --hi B               withhold minima with S - best in [A, B) (e.g. 0 .. 1e-3: the near-record set)
              list    --sides FILE                visible = minima matched by a json list of sides (e.g. a dated census list,
                                                  runs/known110_all.2026-10-08.json); the rest withheld
            and materialise it: <dir>/starts/*.txt (representatives of a visible sample or all), <dir>/novel_ref.json
            (sides of visible minima + their copies: explorer --novel-ref, so withheld minima count as novel).
  score     after a run's sides are certified (pk.py import-run + census_certify.py certify): the run's minima (side keys ->
            copy -> minimum) vs the frozen holdout: visible re-found, withheld rediscovered (+ first time), genuinely new
            (not in the snapshot), per measured CPU-h from the runs' manifests.

  python pk/refset.py snapshot --n 110 --name census110_1011
  python pk/refset.py holdout --ref census110_1011 --rule random --frac 0.5 --seed 1 --name h110_rand50 [--starts 40]
  python pk/refset.py score --holdout h110_rand50 RUN_DIR ...
"""
from __future__ import annotations
import argparse, glob, hashlib, json, math, os, random, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from pk.store import Store
from pk import identity as ID

ROOT = os.path.join(HERE, '..', 'refsets')     # tracked in git (frozen references are part of the experiment record)


def _sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def representative(st, n, m):
    """Lowest-pen feasible packing at S (within 1e-9), else at one of the minimum's copy side keys."""
    S = float(m['S'][:20])
    r = st.db.execute("SELECT id, s FROM packing WHERE n=? AND s BETWEEN ? AND ? AND status!='infeasible' ORDER BY pen LIMIT 1",
                      (n, S - 1e-9, S + 1e-9)).fetchone()
    if r:
        return r['id'], r['s']
    for (k,) in st.db.execute('SELECT side_key FROM copy WHERE n=? AND min_key=?', (n, m['key'])):
        r = st.db.execute("SELECT id, s FROM packing WHERE n=? AND s BETWEEN ? AND ? AND status!='infeasible' ORDER BY pen LIMIT 1",
                          (n, float(k) - 5e-10, float(k) + 5e-10)).fetchone()
        if r:
            return r['id'], r['s']
    return None, None


def snapshot(n, name, k=None):
    st = Store()
    k = k or math.ceil(math.sqrt(n) - 1e-12)
    M = []
    for m in st.minima(n):
        if float(m['S'][:20]) >= k:
            continue
        pid, s_rep = representative(st, n, m)
        sides = sorted(r[0] for r in st.db.execute('SELECT side_key FROM copy WHERE n=? AND min_key=?', (n, m['key'])))
        M.append(dict(key=m['key'], S=m['S'], cls=m['cls'], rep=pid, s_rep=s_rep, sides=sides))
    from pk.manifest import git_state
    d = dict(name=name, n=n, k=k, taken=time.strftime('%Y-%m-%dT%H:%M:%S'), git=git_state()['commit'], minima=M)
    d['sha256'] = _sha(M)
    os.makedirs(ROOT, exist_ok=True)
    json.dump(d, open(f'{ROOT}/{name}.json', 'w'), indent=0)
    nc = sum(m['cls'] == ID.CERTIFIED for m in M)
    print(f'{ROOT}/{name}.json: n={n} {len(M)} minima ({nc} certified, {len(M) - nc} weak), '
          f'{sum(m["rep"] is None for m in M)} without a representative, sha {d["sha256"][:12]}')
    return d


def load(name):
    return json.load(open(name if name.endswith('.json') else f'{ROOT}/{name}.json'))


def holdout(ref, rule, name, frac=0.5, seed=1, lo=0.0, hi=1e-3, sides=None, nstarts=None):
    R = load(ref)
    M = R['minima']
    best = min(float(m['S'][:20]) for m in M)
    if rule == 'random':
        rng = random.Random(seed)
        wh = {m['key'] for m in M if rng.random() < frac}
    elif rule == 'side':
        wh = {m['key'] for m in M if lo <= float(m['S'][:20]) - best < hi}
    elif rule == 'list':
        L = sorted(float(x) for x in json.load(open(sides)))
        import bisect
        def hit(m):
            for x in [float(m['S'][:20])] + [float(s) for s in m['sides']]:
                i = bisect.bisect_left(L, x - 1e-8)
                if i < len(L) and abs(L[i] - x) < 1e-8:
                    return True
            return False
        wh = {m['key'] for m in M if not hit(m)}
    else:
        raise SystemExit(f'rule {rule}?')
    vis = [m for m in M if m['key'] not in wh]
    d = dict(name=name, ref=R['name'], ref_sha256=R['sha256'], n=R['n'], k=R['k'], rule=rule,
             params=dict(frac=frac, seed=seed, lo=lo, hi=hi, sides=sides), visible=sorted(m['key'] for m in vis),
             withheld=sorted(wh), made=time.strftime('%Y-%m-%dT%H:%M:%S'))
    out = f'{ROOT}/{name}'
    os.makedirs(f'{out}/starts', exist_ok=True)
    json.dump(d, open(f'{out}/holdout.json', 'w'), indent=0)
    ref_sides = sorted({float(m['S'][:20]) for m in vis} | {float(s) for m in vis for s in m['sides']})
    json.dump(ref_sides, open(f'{out}/novel_ref.json', 'w'))
    st = Store()
    rng = random.Random(seed + 1)
    pool = [m for m in vis if m['rep']]
    pick = pool if not nstarts or nstarts >= len(pool) else rng.sample(pool, nstarts)
    for m in pick:
        st.get(m['rep']).write(f'{out}/starts/m{float(m["S"][:20]):.10f}.txt')
    wbest = min((float(m['S'][:20]) for m in M if m['key'] in wh), default=None)
    print(f'{out}: {rule} -> {len(vis)} visible, {len(wh)} withheld (best withheld {wbest}; record '
          f'{"withheld" if any(float(m["S"][:20]) == best and m["key"] in wh for m in M) else "visible"}); '
          f'{len(pick)} starts, novel_ref {len(ref_sides)} sides')
    return d


def score(hname, runs, n=None):
    H = json.load(open(f'{ROOT}/{hname}/holdout.json'))
    R = load(H['ref'])
    assert R['sha256'] == H['ref_sha256'], 'snapshot changed since the holdout was made'
    n = H['n']; k = H['k']
    st = Store()
    SM = st.side_minima(n)
    vis, wh, allk = set(H['visible']), set(H['withheld']), {m['key'] for m in R['minima']}
    found, first = {}, {}
    cpu = 0.0; unl = 0; notmin = 0
    for d in runs:
        mf = f'{d}/manifest.json'
        if os.path.exists(mf):
            cpu += json.load(open(mf)).get('cpu_s') or 0.0
        for l in open(f'{d}/archive.jsonl'):
            e = json.loads(l)
            if e['s'] >= k or e.get('parent', 0) < 0:
                continue
            sk = ID.side_key(e['s'])
            if sk not in SM:
                b = st.basin(n, sk)
                if b is None: unl += 1
                elif b['cls'] == 'not-min': notmin += 1
                continue
            mk, cls, S = SM[sk]
            if cls != ID.CERTIFIED:
                continue
            found.setdefault(mk, S)
            first[mk] = min(first.get(mk, 1e18), e.get('t', 0))
    f = set(found)
    out = dict(holdout=hname, runs=runs, cpu_h=cpu / 3600, visible_refound=len(f & vis), withheld_rediscovered=len(f & wh),
               withheld_total=len(wh), new=len(f - allk), unlabelled_sides=unl, notmin_sides=notmin,
               withheld_first_t=sorted(first[x] for x in f & wh), new_S=sorted(found[x] for x in f - allk))
    h = out['cpu_h'] or float('nan')
    print(f'{hname}: {len(runs)} runs, {out["cpu_h"]:.2f} measured CPU-h; withheld rediscovered {out["withheld_rediscovered"]} / '
          f'{len(wh)} ({out["withheld_rediscovered"] / h:.1f} / CPU-h); new {out["new"]} ({out["new"] / h:.2f} / CPU-h); '
          f'visible re-found {out["visible_refound"]}; uncertified sides {unl} (certify first), not-min {notmin}')
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('snapshot'); p.add_argument('--n', type=int, required=True); p.add_argument('--name', required=True)
    p = sp.add_parser('holdout'); p.add_argument('--ref', required=True); p.add_argument('--name', required=True)
    p.add_argument('--rule', choices=['random', 'side', 'list'], required=True); p.add_argument('--frac', type=float, default=0.5)
    p.add_argument('--seed', type=int, default=1); p.add_argument('--lo', type=float, default=0.0); p.add_argument('--hi', type=float, default=1e-3)
    p.add_argument('--sides'); p.add_argument('--starts', type=int)
    p = sp.add_parser('score'); p.add_argument('--holdout', required=True); p.add_argument('runs', nargs='+')
    p.add_argument('--json')
    a = ap.parse_args()
    if a.cmd == 'snapshot':
        snapshot(a.n, a.name)
    elif a.cmd == 'holdout':
        holdout(a.ref, a.rule, a.name, a.frac, a.seed, a.lo, a.hi, a.sides, a.starts)
    else:
        out = score(a.holdout, a.runs)
        if a.json:
            json.dump(out, open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    main()
