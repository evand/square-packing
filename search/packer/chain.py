#!/usr/bin/env python3
"""Chaining across sizes (10-08; SQUISH-style lineage): a window of n, each starting from its frozen 10-05 register record
(../exact/batch/inputs: pre-SQUISH, so SQUISH's sides are ground truth for "how much does the chain recover").

Round for one n: an explorer run (P procs x T min) from best(n) plus removal seeds best(n + k) - k, k = 1, 2 (best of 4
random removals each, fq-quenched).  A new best at n marks n - 1 and n - 2 dirty (fresh seeds).  Scheduling: dirty n
first (most recently dirtied), otherwise the n with the best recent gain per round (ties: fewest rounds).  The top of the
window (n_hi + 1, n_hi + 2) stays at its register record (seed source only).

  chain.py --lo 123 --hi 131 --minutes 4 --procs-per 2 --slots 8 --hours 1 --out runs/ch1
  chain.py --report runs/ch1
"""
import argparse, glob, json, math, os, random, shutil, signal, subprocess, sys, tempfile, time
from jobpool import run_jobs
import mcmin, hop
from pk.engine import QOpt
from pk.moves import remove_k, insert_k

HERE = os.path.dirname(os.path.abspath(__file__))
TRIES, SEED_FIRST = 12, False                # removal tries per seed; v2 (10-08): non-shrinking removals rejected
SEED_Q = QOpt(loosen=1.0)                    # seed quench: full ALM + polish, no dilation


# a CANDIDATE must beat the register by more than re-polishing the same packing can (10-12: n=207 re-polished
# BuddenPing's #503 by 1.5e-9 at matching distance 0)
CAND_TOL = 1e-8


def asof_local(T):
    """UTC cutoff (eventbench name or ISO) -> the store's 'added' format (local wall time, no zone; TODO note: store
    times are local, not UTC)."""
    from pk import eventbench as EB
    return EB.cutoff(T).astimezone().strftime('%Y-%m-%dT%H:%M:%S')


_STORE_SIDES = {}
def store_sides(n, asof=None):
    """Sides of every feasible store packing at n (added before `asof`, store-local ISO): the global-novelty reference
    (10-08 decision: the bandit reward must be global novelty; ex3 showed run-local novelty piles onto kicksym)."""
    if n not in _STORE_SIDES:
        import sqlite3
        from pk.store import DEFAULT as DB
        db = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
        q, args = "SELECT s FROM packing WHERE n=? AND status != 'infeasible'", [n]
        if asof:
            q += ' AND added < ?'; args.append(asof)
        _STORE_SIDES[n] = [r[0] for r in db.execute(q, args)]
        db.close()
    return _STORE_SIDES[n]


def novel_ref(out, n, asof, path):
    """store sides at n + every archive side of this chain's earlier rounds at n -> sorted json list for --novel-ref"""
    ref = set(store_sides(n, asof))
    for f in glob.glob(f'{out}/r*_n{n}/archive.jsonl'):
        for l in open(f):
            try:
                ref.add(json.loads(l)['s'])
            except Exception:
                pass
    json.dump(sorted(ref), open(path, 'w'))
    return len(ref)


def ledger_round(rd, n, E, params, fr):
    """Expansion ledger (10-14): one store row per lineage root of this round: CPU spent on its descendants' proposals,
    proposal count, best side reached.  Roots are added to the store (dedup by hash) so the ledger keys on store ids."""
    from pk.store import Store
    from pk.packing import read
    if not E or not os.path.exists(f'{rd}/proposals.jsonl'):
        return 0
    root = {e['i']: e.get('root', e['i']) for e in E}
    agg = {}
    for l in open(f'{rd}/proposals.jsonl'):
        q = json.loads(l)
        r = root.get(q.get('par'))
        if r is None:
            continue
        c = q.get('cpu') if q.get('cpu') is not None else q.get('sec', 0.0)
        a = agg.setdefault(r, [0.0, 0]); a[0] += c or 0.0; a[1] += 1
    best = {}
    for e in E:
        r = root[e['i']]
        best[r] = min(best.get(r, math.inf), e['s'])
    st = Store()
    try:
        for e in E:
            if e['parent'] >= 0 or e['i'] not in agg:
                continue
            p = read(os.path.join(rd, os.path.basename(e['path'])))
            pid = st.find_hash(p) or st.add(p, source=f'run:{os.path.relpath(rd, HERE)}', kind=e['kind'], status='polished',
                                            commit=False)[0]
            cpu, k = agg[e['i']]
            st.record_expansion(pid, n, 'chain', cpu, params=dict(params, start=e['kind']), run=os.path.relpath(rd, HERE),
                                nprop=k, best_s=best.get(e['i']), frontier_s=fr, commit=False)
        st.db.commit()
    finally:
        st.close()
    return len(agg)


BANDS = (70, 130, 200)        # n bands for pooled move-kind statistics (analysis/bandit_timescale.md: kind effects vary mildly with n)


def bandit_file(a, n):
    if a.bandit == 'band':
        return f'{a.out}/bandit_band{sum(n >= b for b in BANDS)}.json'
    return f'{a.out}/bandit_{n}.json'


def frontier_s(n, asof_T):
    """best known side for the value reward: eventbench frontier (as of T for time travel, else live)"""
    from pk import eventbench as EB
    return EB.frontier_at(n, asof_T) if asof_T else EB.frontier_now(n)

def seed_job(j):
    """k > 0: best of TRIES random k-removals from best(n + k); k < 0: insert |k| squares into best(n + k) at pose holes
    after scaling the packing up by a random gap.  fq-quenched; a removal that does not shrink is a copy and is rejected."""
    src, k, n, out, seed = j
    s1, sq1 = mcmin.load_deg(src)
    rng = random.Random(seed)
    tmp = tempfile.mkdtemp()
    best = None
    for t in range(TRIES):
        if k > 0:
            _, sq2 = remove_k((s1, sq1), k, rng)
            r = hop.quench(s1, sq2, tmp, opt=SEED_Q)
            if not r or r[0] > s1 - 1e-7:      # did not shrink (removed a rattler / polished back): a copy, not a seed
                continue
        else:
            s2, sq2 = insert_k((s1, sq1), -k, rng, holes='pose')   # dilate, then any of the 12 best pose holes
            if len(sq2) != n:
                continue
            r = hop.quench(s2, sq2, tmp, opt=SEED_Q)
            if not r:
                continue
        if best is None or r[0] < best[0]:
            best = (r[0], r[1])
        if SEED_FIRST and best:
            break
    if best:
        mcmin.write_deg(out, *best)
        return best[0]
    return None


LIB_CAP = 300


def lib_load(out, n):
    """per-n library of distinct sub-k packings: [(s, path)] from <out>/lib/n<n>/index.json"""
    f = f'{out}/lib/n{n}/index.json'
    return [tuple(e) for e in json.load(open(f))] if os.path.exists(f) else []


def lib_add(out, n, items):
    """add (s, path) packings (copied) to n's library; distinct = side differs by >= identity.SAME_SIDE (2e-9; 10-11, was
    > 1e-9); keeps the LIB_CAP best"""
    d = f'{out}/lib/n{n}'; os.makedirs(d, exist_ok=True)
    L = lib_load(out, n)
    sides = sorted(x for x, _ in L)
    import bisect
    for s, p in sorted(items):
        from pk.identity import SAME_SIDE
        i = bisect.bisect_left(sides, s - SAME_SIDE)
        if i < len(sides) and abs(sides[i] - s) < SAME_SIDE or not os.path.exists(p):
            continue
        dst = f'{d}/{s:.12f}.txt'; shutil.copy(p, dst); L.append((s, dst)); bisect.insort(sides, s)
    L.sort()
    for s, p in L[LIB_CAP:]:
        if os.path.exists(p): os.remove(p)
    L = L[:LIB_CAP]
    json.dump(L, open(f'{d}/index.json.tmp', 'w')); os.replace(f'{d}/index.json.tmp', f'{d}/index.json')
    return len(L)


def archive_best(d):
    best = None
    if os.path.exists(f'{d}/archive.jsonl'):
        for l in open(f'{d}/archive.jsonl'):
            e = json.loads(l)
            if best is None or e['s'] < best[0]:
                best = (e['s'], e['path'])
    return best


def main(a):
    os.makedirs(a.out, exist_ok=True)
    from pk import manifest                  # explore.py rounds write their own manifests in their round dirs
    man = manifest.start(a.out, 'chain', args=vars(a), procs=a.slots * a.procs_per)
    stp = f'{a.out}/state.json'
    T = hop.squish_targets()
    if os.path.exists(stp):
        S = json.load(open(stp))
    else:
        S = {'n': {}, 'log': [], 't0': time.time(), 'round': 0}
        src = a.src or hop.BATCH
        ns = json.load(open(a.ns_file)) if a.ns_file else list(range(a.lo, a.hi + 1))
        ours = dict(x.split(':') for x in (a.ours or []))
        want = sorted(set(ns) | {m + k for m in ns for k in a.seed_ks})
        miss = [n for n in want if not os.path.exists(f'{src}/n-{n}.txt')]
        if a.src and miss and not a.allow_missing:     # 10-14: no silent fallback to the frozen 10-05 register
            sys.exit(f'--src {a.src} lacks n = {miss} (pk.py export them, or --allow-missing to skip those n)')
        for n in want:
            p = f'{a.out}/best_{n}.txt'
            sp = f'{src}/n-{n}.txt'
            if not os.path.exists(sp):
                print(f'skip n={n}: not in {src}', flush=True)
                continue
            shutil.copy(sp, p)
            reg, _ = mcmin.load_deg(p)
            if str(n) in ours:
                so, _ = mcmin.load_deg(ours[str(n)])
                if so < reg - 1e-12:
                    shutil.copy(ours[str(n)], p)
            s, _ = mcmin.load_deg(p)
            S['n'][str(n)] = dict(best=s, reg=reg, start=s, squish=None if a.src else T.get(n), path=p, rounds=0,
                                  dirty=(n in ns), gains=[], dirtied=0.0, active=(n in ns))
        if a.lib_init:
            for n_, ps in json.load(open(a.lib_init)).items():
                k_ = math.ceil(math.sqrt(int(n_)) - 1e-12)
                its = []
                for p in ps:
                    try:
                        its.append((mcmin.load_deg(p)[0], p))
                    except Exception:
                        pass
                lib_add(a.out, int(n_), [(s_, p) for s_, p in its if s_ < k_ - 1e-9])
        for q, x in enumerate(a.carry or []):       # extra lineage starts (e.g. our best non-grid packing above a grid bound)
            n, pth = x.split(':')
            dst = f'{a.out}/carry_{n}_x{q}.txt'; shutil.copy(pth, dst)
            S['n'][n].setdefault('carry', []).append(dst)
        if a.schedule:                              # TT3 replay (10-15): fixed [[round, n], ...] list, run in order
            S['sched'] = [[int(r), int(n)] for r, n in json.load(open(a.schedule))]
            bad = [n for _, n in S['sched'] if str(n) not in S['n']]
            if bad:
                sys.exit(f'--schedule names n = {sorted(set(bad))} that are not in the n set')
    live = {}
    t_end = time.time() + 3600 * a.hours
    def pick():
        cand = [int(n) for n, v in S['n'].items() if v.get('active', int(n) <= (a.hi or 0)) and int(n) not in live.values()]
        dirty = [n for n in cand if S['n'][str(n)]['dirty']]
        if dirty:
            return max(dirty, key=lambda n: (S['n'][str(n)]['dirtied'], -n))
        if not cand:
            return None
        if a.first_pass and any(S['n'][str(n)]['rounds'] == 0 for n in cand):
            cand = [n for n in cand if S['n'][str(n)]['rounds'] == 0]
        return max(cand, key=lambda n: (sum(S['n'][str(n)]['gains'][-2:]), -S['n'][str(n)]['rounds']))
    while (time.time() < t_end and (not a.schedule or S['sched'])) or live:
        while time.time() < t_end and len(live) < a.slots:
            if a.schedule:
                # replay: the first scheduled round whose n is not running (keeps each n's rounds in order); the round
                # number is the scheduled one, so seeds, removal seeds and library picks match the replayed run
                nxt = next((i for i, (_, m) in enumerate(S['sched']) if m not in live.values()), None)
                if nxt is None:
                    break
                rr, n = S['sched'].pop(nxt)
            else:
                n = pick()
                if n is None:
                    break
            v = S['n'][str(n)]
            v['dirty'] = False
            S['round'] = rr if a.schedule else S['round'] + 1
            rd = f'{a.out}/r{S["round"]:04d}_n{n}'
            if os.path.exists(rd):      # left by a run killed mid-round: set it aside (lineage ignores *.partial*)
                os.rename(rd, f'{rd}.partial{int(time.time())}')
            os.makedirs(rd)
            seeds = [(S['n'][str(n + k)]['path'], k, n, f'{rd}/seed_k{k}.txt', S['round'] * 10 + k) for k in a.seed_ks if str(n + k) in S['n']]
            res = run_jobs(seed_job, seeds, procs=a.procs_per, timeout=1800)
            carried = []
            for q, cp in enumerate(v.get('carry', [])):         # lineage carry from the previous round(s)
                dst = f'{rd}/seed_c{q}.txt'; shutil.copy(cp, dst); carried.append(dst)
            libs = []
            if a.lib == 'store':                                 # 10-14: store-wide, uniform over clusters, ledger-weighted
                try:
                    from pk import parents
                    got = parents.starts(n, a.lib_k, rd, policy=a.lib_policy, asof=a.as_of, seed=S['round'], prefix='seed_l')
                    libs = [g[0] for g in got]
                    S['log'].append(dict(ev='lib', n=n, round=S['round'], picks=[[g[1], g[2], g[3]] for g in got]))
                except Exception as ex:
                    print(f'lib n={n}: {ex!r}', flush=True)
            else:                                                # run-local sub-grid library (Evan 10-09)
                L = [e for e in lib_load(a.out, n) if e[0] > v['best'] + 1e-9]
                for q, (ls, lp) in enumerate(random.Random(S['round']).sample(L, min(a.lib_k, len(L)))):
                    dst = f'{rd}/seed_l{q}.txt'; shutil.copy(lp, dst); libs.append(dst)
            starts = [v['path']] + carried + libs + [sd[3] for sd, r in zip(seeds, res) if isinstance(r, float)]
            cmd = [sys.executable, os.path.join(HERE, 'explore.py'), '--n', str(n), '--starts', *starts, '--procs', str(a.procs_per),
                   '--minutes', str(a.minutes), '--seed', str(S['round']), '--out', rd, '--adapt', '--elite-share', '0.3',
                   '--frontier', f'{rd}/seed', '--frontier-share', str(a.frontier_share), '--bandit-state', bandit_file(a, n),
                   '--reward', a.reward] + (['--bandit-half', 'inf', '--bandit-carry', '1'] if a.bandit == 'band' else [])
            if a.novel == 'store':
                novel_ref(a.out, n, asof_local(a.as_of) if a.as_of else None, f'{rd}/novel_ref.json')
                cmd += ['--novel-ref', f'{rd}/novel_ref.json']
            if a.reward == 'value':
                cmd += ['--frontier-s', repr(frontier_s(n, a.as_of))]
            cmd += a.explore_extra.split() if a.explore_extra else []
            p = subprocess.Popen(cmd, stdout=open(f'{rd}/out', 'w'), stderr=subprocess.STDOUT, start_new_session=True)
            live[p] = n
            S['log'].append(dict(t=round(time.time() - S['t0']), ev='start', n=n, round=S['round'], seeds=[r for r in res],
                                 seeds_src={str(sd[1]): S['n'][str(n + sd[1])]['best'] for sd in seeds}))   # lineage (chain_lineage.py)
        time.sleep(3)
        for p, n in list(live.items()):
            if p.poll() is None:
                continue
            del live[p]
            v = S['n'][str(n)]
            rd = p.args[p.args.index('--out') + 1]
            b = archive_best(rd)
            v['rounds'] += 1
            gain = 0.0
            if b and b[0] < v['best'] - 1e-9:
                gain = v['best'] - b[0]
                shutil.copy(b[1], v['path'])
                v['best'] = b[0]
                if b[0] < min(v['reg'], v.get('start', v['reg'])) - CAND_TOL:   # beyond the live register and our own start
                    cp = f'{a.out}/cand_{n}_{b[0]:.10f}.txt'; shutil.copy(b[1], cp)
                    print(f'CANDIDATE n={n} {b[0]:.12f} (register {v["reg"]:.12f}, start {v.get("start", v["reg"]):.12f}) -> {cp}', flush=True)
                for m in {n - k for k in a.seed_ks}:
                    if str(m) in S['n'] and S['n'][str(m)].get('active', m >= (a.lo or 0)):
                        S['n'][str(m)]['dirty'] = True
                        S['n'][str(m)]['dirtied'] = time.time()
            # carry: best 3 distinct basins of the non-record lineages (roots = seed / carried starts), below the ceiling
            E = [json.loads(l) for l in open(f'{rd}/archive.jsonl')] if os.path.exists(f'{rd}/archive.jsonl') else []
            roots = {e['i'] for e in E if e['parent'] < 0 and e['kind'].startswith('start:seed')}
            k_n = math.ceil(math.sqrt(n) - 1e-12)
            lin = sorted((e for e in E if e.get('root', e['i']) in roots and e['s'] < k_n + 0.05), key=lambda e: e['s'])
            car, seen = [], []
            for e in lin:
                if all(abs(e['s'] - x) > 1e-7 for x in seen):
                    seen.append(e['s']); dst = f'{a.out}/carry_{n}_{len(car)}.txt'; shutil.copy(e['path'], dst + '.tmp'); car.append(dst)
                if len(car) == 3:
                    break
            for c in car: os.replace(c + '.tmp', c)
            v['carry'] = car
            lib_add(a.out, n, [(e['s'], e['path']) for e in E if e['s'] < k_n - 1e-9 and not e.get('unpolished')])
            try:                                     # 10-14: every round goes into the store (novelty ref, parents)
                from pk.store import Store
                from pk import sources
                _st = Store(); sources.import_run(_st, rd); _st.close()
            except Exception as ex:
                print(f'import n={n}: {ex!r}', flush=True)
            try:
                ledger_round(rd, n, E, dict(reward=a.reward, novel=a.novel, as_of=a.as_of, minutes=a.minutes, procs=a.procs_per),
                             frontier_s(n, a.as_of))
            except Exception as ex:                  # the ledger must never stop a hunt
                print(f'ledger n={n}: {ex!r}', flush=True)
            lb = seen[0] if seen else None
            lgain = (v.get('lin_best', float('inf')) - lb) if lb is not None and v.get('lin_best') is not None else 0.0
            if lb is not None and (v.get('lin_best') is None or lb < v['lin_best']): v['lin_best'] = lb
            gain = gain + max(0.0, lgain) * 0.1        # scheduling credit for lineage progress too (record gains dominate)
            v['gains'].append(gain)
            S['log'].append(dict(t=round(time.time() - S['t0']), ev='done', n=n, best=v['best'], gain=gain,
                                 vs_reg=v['best'] - v['reg'], vs_squish=(v['best'] - v['squish']) if v['squish'] else None))
            print(f'{time.time() - S["t0"]:7.0f}s n={n} lineage {v.get("lin_best") or 0:.10f} best {v["best"]:.10f} gain {gain:.2e} vs reg {v["best"] - v["reg"]:+.2e}'
                  + (f' vs SQUISH {v["best"] - v["squish"]:+.2e}' if v['squish'] else ''), flush=True)
            json.dump(S, open(stp + '.tmp', 'w')); os.replace(stp + '.tmp', stp)
    json.dump(S, open(stp + '.tmp', 'w')); os.replace(stp + '.tmp', stp)
    report(a.out)
    man.finish()


def report(out):
    S = json.load(open(f'{out}/state.json'))
    print(f'== {out}: {S["round"]} rounds')
    cands = [(n, v) for n, v in S['n'].items() if v['best'] < min(v['reg'], v.get('start', v['reg'])) - CAND_TOL]
    print(f'  candidates below register and start: {[(n, round(v["best"], 10), round(v["reg"] - v["best"], 10)) for n, v in cands]}')
    print(f'  {"n":>4} {"rounds":>6} {"register":>15} {"best":>15} {"gain":>9} {"SQUISH":>15} {"frac of SQUISH gain":>19} {"lineage best":>15}')
    for n, v in sorted(S['n'].items(), key=lambda q: int(q[0])):
        sq = v['squish']
        fr = (v['reg'] - v['best']) / (v['reg'] - sq) if sq and v['reg'] - sq > 1e-12 else None
        print(f'  {n:>4} {v["rounds"]:6d} {v["reg"]:15.10f} {v["best"]:15.10f} {v["reg"] - v["best"]:9.2e} '
              f'{(f"{sq:15.10f}") if sq else " " * 15} {(f"{fr:19.2f}") if fr is not None else " " * 19} {v.get("lin_best") or 0:15.10f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--lo', type=int); ap.add_argument('--hi', type=int); ap.add_argument('--minutes', type=float, default=4)
    ap.add_argument('--procs-per', type=int, default=2); ap.add_argument('--slots', type=int, default=8)
    ap.add_argument('--hours', type=float, default=1); ap.add_argument('--frontier-share', type=float, default=0.7);
    ap.add_argument('--explore-extra', default='', help='extra explore.py args (e.g. "--melt 2")'); ap.add_argument('--first-pass', action='store_true', help='every n once before any repeat'); ap.add_argument('--ns-file', help='json list of n to work on (instead of --lo/--hi)'); ap.add_argument('--src', help='start packings dir (n-<n>.txt), e.g. ../exact/batch/inputs_live')
    ap.add_argument('--seed-ks', type=lambda x: [int(v) for v in x.split(',')], default=[1, 2], help='seed offsets k: removal from best(n + k) (k > 0), insertion into it (k < 0)')
    ap.add_argument('--carry', nargs='*', help='n:path extra lineage start for n (carried like a seed lineage)'); ap.add_argument('--ours', nargs='*', help='n:path of our own better packings (used as start if better than src)'); ap.add_argument('--out'); ap.add_argument('--report')
    ap.add_argument('--lib-k', type=int, default=3, help='library entries (random distinct sub-k packings) per round as starts')
    ap.add_argument('--lib-init', help='json {n: [paths]}: sub-k packings to pre-load into the library')
    ap.add_argument('--seed-tries', type=int, default=TRIES, help='quench tries per seed job')
    ap.add_argument('--reward', default='value', choices=['legacy', 'value'],
                    help='explore bandit reward (10-14 default value: novelty graded by gap to the frontier; legacy = pre-10-14 chain)')
    ap.add_argument('--novel', default='store', choices=['store', 'run'],
                    help='novelty reference: store (+ this chain\'s earlier rounds; default) or run (pre-10-14: run-local)')
    ap.add_argument('--as-of', help='time travel: eventbench cutoff name or UTC ISO; novelty ref = store rows added before it, '
                    'value-reward frontier = frontier_at(n, T).  Starts come from --src (eventbench starts(T)).')
    ap.add_argument('--lib', default='store', choices=['store', 'run'],
                    help='library starts: store (pk/parents, default 10-14) or run (pre-10-14 run-local sub-k library)')
    ap.add_argument('--lib-policy', default='banded', help='pk/parents policy for --lib store (banded | cluster | spread | best)')
    ap.add_argument('--bandit', default='band', choices=['band', 'n'],
                    help='move-kind statistics: pooled per n band, no decay, full carry (default 10-14, A1) or per n (pre-10-14: 600 s wall half-life, x0.5 carry)')
    ap.add_argument('--schedule', help='json [[round, n], ...]: replay exactly these rounds in this order with these round '
                    'numbers (seeds, library picks), instead of the dirty-first scheduler (TT3, 10-15); ends when done')
    ap.add_argument('--allow-missing', action='store_true', help='skip n missing from --src instead of failing')
    a = ap.parse_args()
    TRIES = a.seed_tries
    if a.report:
        report(a.report)
    else:
        import scratch
        scratch.scope('chain')               # seed-job temp dirs under one run dir, removed at exit (explore.py scopes its own)
        main(a)
