#!/usr/bin/env python3
"""TT3 feature-test scoring (pre-registered in PACKER.md "TT3", 10-15).

Arms = chain.py --schedule replays of TT2 arm S's 49 rounds on 39 n <= 240 (runs/tt3/schedule.json), plus TT2 S itself
restricted to those rounds.  Per arm and n: near-frontier new minima = distinct side keys (2e-9) of polished / certified,
grid-free, sub-k store rows sighted in the arm's rounds, within NEAR (1e-3) of frontier_at(n, T) and > 2e-9 from every
as-of-T reference side (as analysis/tt2_census.py).  CPU = sum of the rounds' explore manifests (cpu_s; the same for every
arm, so seed-job CPU is left out everywhere).  Ratio of rates with a bootstrap over n (B = 4000).  Also: proposals per
CPU-s, polished share, best gap per n.

  python3 analysis/tt3_score.py runs/tt3/G runs/tt3/A [--base runs/tt2/S] [--json runs/tt3/score.json]
"""
import argparse, collections, glob, json, math, os, random, re, sqlite3, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
sys.path.insert(0, PK); sys.path.insert(0, HERE)
import tt2_census as TC                                                      # noqa: E402
from pk import eventbench as EB, parents as PA                               # noqa: E402

NEAR = 1e-3


def arm_rounds(d, sched):
    want = {(r, n) for r, n in sched}
    R = {}
    for rd in glob.glob(f'{d}/r*_n*'):
        m = re.search(r'/r(\d+)_n(\d+)$', rd)
        if m and (int(m[1]), int(m[2])) in want:
            R[(int(m[1]), int(m[2]))] = rd
    return R


def score_arm(db, d, sched, ref, fT):
    R = arm_rounds(d, sched)
    per = collections.defaultdict(lambda: dict(near=0, cpu=0.0, props=0, pol=0, best=None))
    for (r, n), rd in R.items():
        try:
            per[n]['cpu'] += json.load(open(f'{rd}/manifest.json'))['cpu_s']
        except Exception:
            pass
        if os.path.exists(f'{rd}/proposals.jsonl'):
            for l in open(f'{rd}/proposals.jsonl'):
                p = json.loads(l); per[n]['props'] += 1
                per[n]['pol'] += p['st'] == 'new?' and not p.get('unp')
    rel = os.path.relpath(d, PK)
    for n in {n for _, n in sched}:
        k = math.ceil(math.sqrt(n) - 1e-12)
        got = sorted(s for (s,) in db.execute(
            "SELECT DISTINCT p.s FROM packing p JOIN sighting g ON g.packing=p.id WHERE p.n=? AND p.s<? AND p.s<? "
            "AND p.status IN ('polished','certified') AND (p.grid IS NULL OR p.grid=0) AND g.source LIKE ?",
            (n, k, fT[n] + NEAR, f'run:{rel}/r%_n{n}%')))
        keys = []
        for s in got:
            if (not keys or s - keys[-1] > 2e-9) and not TC.near_key(ref[n], s):
                keys.append(s)
        per[n]['near'] = len(keys)
        per[n]['best'] = (min(keys) - fT[n]) if keys else None
    return dict(per), len(R)


def boot(pa, pb, ns, B=4000, seed=1):
    ca, cb = sum(pa[n]['cpu'] for n in ns), sum(pb[n]['cpu'] for n in ns)
    rng, rs = random.Random(seed), []
    for _ in range(B):
        smp = [rng.choice(ns) for _ in ns]
        a, b = sum(pa[n]['near'] for n in smp), sum(pb[n]['near'] for n in smp)
        if a and b:
            rs.append((a / ca) / (b / cb))
    rs.sort()
    a, b = sum(pa[n]['near'] for n in ns), sum(pb[n]['near'] for n in ns)
    if not (a and b and rs):
        return None, None, None
    return (a / ca) / (b / cb), rs[int(.025 * len(rs))], rs[int(.975 * len(rs))]


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('arms', nargs='+'); ap.add_argument('--base', default='runs/tt2/S')
    ap.add_argument('--schedule', default='runs/tt3/schedule.json'); ap.add_argument('--json')
    a = ap.parse_args()
    os.chdir(PK)
    sched = json.load(open(a.schedule)); ns = sorted({n for _, n in sched})
    db = sqlite3.connect(f'file:{TC.DB}?mode=ro', uri=True)
    asof = PA.asof_local(TC.T)
    fT = {n: EB.frontier_at(n, TC.T) for n in ns}
    ref = {n: sorted(s for (s,) in db.execute("SELECT s FROM packing WHERE n=? AND s<? AND status!='infeasible' AND added<? "
                                               "AND (grid IS NULL OR grid=0)", (n, math.ceil(math.sqrt(n) - 1e-12), asof)))
           for n in ns}
    res = {}
    for d in [a.base] + a.arms:
        per, nr = score_arm(db, d, sched, ref, fT)
        cpu = sum(v['cpu'] for v in per.values()); props = sum(v['props'] for v in per.values())
        near = sum(v['near'] for v in per.values())
        res[d] = dict(per=per, rounds=nr, cpu_h=cpu / 3600, props=props, near=near)
        print(f'{d}: rounds {nr}/{len(sched)}  CPU-h {cpu / 3600:.2f}  proposals {props} ({props / max(cpu, 1):.2f} per CPU-s, '
              f'polished {sum(v["pol"] for v in per.values()) / max(props, 1):.0%})  near-frontier new {near} '
              f'({near / max(cpu / 3600, 1e-9):.1f} per CPU-h)  n with any {sum(1 for v in per.values() if v["near"])}')
    names = [a.base] + a.arms
    for i, x in enumerate(names):
        for y in names[:i]:
            r, lo, hi = boot(res[x]['per'], res[y]['per'], ns)
            print(f'  {x} / {y}: rate ratio {r:.2f}  95 % bootstrap over n [{lo:.2f}, {hi:.2f}]' if r else f'  {x} / {y}: no events')
    if a.json:
        json.dump(dict(T=TC.T, schedule=a.schedule, arms=res), open(a.json, 'w'), indent=1)
