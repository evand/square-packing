"""Screen-then-promote construction (10-25): relax as a slow move with a cheap low-fidelity version.

A chain arm's construction is `build stages | relax.X`.  The relax is ~all of its CPU (10-25 profile: the build stages
~0 s, relax.round.15c.p5 9-23 s median at n 549-978).  A cheap screen predicts the relax: fq's rounded compress quench
alone, a few outer ALM iterations (`--round-r0 0.15 --round-compress --round-steps 0 --no-alm --no-polish`), scored by
its side s_soft.  10-25 bench (one construction arm at large n, ~40 builds per n): 2 outer iterations, 0.1-0.25 s, Spearman
0.70-0.81 vs the full relax; the best full result ranked 0-4 of ~40 by the screen.

construct(): interleaves build + screen and full relaxes.  While screens < ratio x (full relaxes + 1), build and
screen one more; otherwise fully relax the best unrelaxed candidate by the screen score at `k` outer iterations, or
with probability `audit` a uniformly random unrelaxed one (keeps the screen's calibration measurable online).  Every
screen runs to `kmax` outer iterations with fq --traj, so the record holds the score at every k <= kmax and any (k,
ratio) policy can be backtested at equal CPU (backtest()).  ratio 1 and audit 1 fully relax everything: the dedicated
calibration dataset, same code.

Goal field `screen: {k: 2, kmax: 8, ratio: 5, audit: 0.1}` (None: the old build-and-relax-each loop).
"""
from __future__ import annotations
import math, os, random, tempfile
import numpy as np

DEFAULT = dict(k=2, kmax=8, ratio=16, audit=0.05)   # 10-26: lean, untested (CLAIMS C53)


def split_chain(spec):
    """'a | b | relax.X' -> ('a | b', 'relax.X'); None if the chain does not end in a relax stage"""
    parts = [p.strip() for p in str(spec).split('|')]
    if len(parts) < 2 or not parts[-1].startswith('relax.'):
        return None
    return ' | '.join(parts[:-1]), parts[-1]


def run_screen(p, kmax, cap=30.0):
    """the compress screen of Packing p -> dict(traj=[(outer, s, sec)], sec, failed)"""
    from .engine import QOpt, fq_quench
    opt = QOpt(round_r0=0.15, round_compress=True, round_steps=0, no_alm=True, no_polish=True, outer=int(kmax),
               traj=True, deadline=float(cap))
    with tempfile.TemporaryDirectory(prefix='scr_') as td:
        try:
            r = fq_quench(p, opt, tmp=td, mod90=True, timeout=cap + 60)
            d = r.stats(whole=True)
        except Exception as e:
            return dict(traj=[], sec=None, failed=repr(e)[:80])
    tr = [(int(q[0]), float(q[1]), float(q[4])) for q in d.get('traj') or []]
    failed = 'diverged' if d.get('diverged') else 'deadline' if d.get('deadline') else None
    return dict(traj=tr, sec=d.get('sec'), failed=failed)


def score_at(traj, k):
    """screen score at k outer iterations: the side after outer k (the last one reached if the quench stopped early)"""
    v = [s for o, s, _ in traj if o <= k]
    return v[-1] if v else None


def construct(chain, ctx, budget_cpu, cpu, rng, pick_start, keep, cfg=None, log=None):
    """Screen-then-promote construction.  chain: the arm's chain (must end in relax.X); ctx: stage ctx (n, fq_timeout
    set per relax by the caller's cap function ctx['cap']()); budget_cpu: CPU s for this phase; cpu(): CPU clock;
    pick_start() -> (pid, State or None); keep(st, lg) -> bool (the caller's filter for a finished build).
    Returns (kept [(s, pid, sq)], records)."""
    from . import stages as S
    cfg = dict(DEFAULT, **(cfg or {}))
    pre, last = split_chain(chain)
    c0 = cpu(); recs = []; pool = []; kept = []; fulls = 0; screens = 0
    while cpu() - c0 < budget_cpu:
        if not pool or screens < cfg['ratio'] * (fulls + 1):
            pid, st0 = pick_start()
            t = cpu()
            try:
                st, lg = S.run_chain(pre, dict(ctx), rng.randrange(1 << 30), st=st0, keep_best=False)
            except Exception as e:
                recs.append(dict(pid=pid, failed=repr(e)[:80])); continue
            if any(r['info'].get('failed') for r in lg) or st.p.n != int(ctx['n']):
                recs.append(dict(pid=pid, failed='build')); continue
            sc = run_screen(st.p, cfg['kmax'])
            screens += 1
            rec = dict(pid=pid if pid is not None else st.meta.get('parent'), raw_s=float(st.p.s), traj=sc['traj'],
                       scr_sec=sc['sec'], scr_cpu=cpu() - t, scr_failed=sc['failed'])
            recs.append(rec)
            if sc['failed'] == 'diverged' or score_at(sc['traj'], cfg['k']) is None:
                continue
            pool.append((rec, st))
            continue
        if rng.random() < cfg['audit']:
            i = rng.randrange(len(pool)); how = 'audit'
        else:
            i = min(range(len(pool)), key=lambda j: score_at(pool[j][0]['traj'], cfg['k'])); how = 'rank'
        rec, st = pool.pop(i)
        t = cpu()
        st2, lg = S.run_chain(last, dict(ctx, fq_timeout=ctx['cap']()), rng.randrange(1 << 30), st=st)
        fulls += 1
        info = lg[0]['info']; bad = S.failed(lg)                 # a keep.best row: the build itself was the best
        rec.update(promoted=how, full_s=float(st2.p.s) if not bad else None, full_cpu=cpu() - t,
                   full_failed=info.get('failed'), deadline=bool(info.get('deadline')),
                   kept_best=lg[-1]['info'].get('src') if lg[-1]['stage'] == 'keep.best' else None)
        if ctx.get('on_relax'):
            ctx['on_relax'](lg)
        if not bad and keep(st2, lg):
            kept.append((float(st2.p.s), rec['pid'], np.asarray(st2.p.sq, float)))
    if log is not None:
        log.update(screened=screens, full_relaxed=fulls, screen_cfg=cfg)
    return kept, recs


def construct_serial(chain, ctx, budget_cpu, cpu, rng, pick_start, keep, cfg=None, hist=(), log=None):
    """Serial screen-then-promote (10-26, for short units): build one, screen it, and fully relax it at once if its
    score at k is in the best 1 / ratio of the screens seen at this arm (`hist`, earlier units' scores, plus this
    run's); otherwise drop it.  An audit share (`audit` of the full relaxes) relaxes a build regardless of its score.
    No pool, so a unit can stop at any point.  Equal-CPU backtest on the audit-1 calibration set (30 arm x n cells,
    91 CPU-s budgets): same best side per unit as construct() (+-0.02 sd), ~10 % more top-5 % results per CPU-h.
    Until 1 / ratio scores exist (no history), only audits are relaxed.  Returns (kept, records) like construct()."""
    from . import stages as S
    cfg = dict(DEFAULT, **(cfg or {}))
    pre, last = split_chain(chain)
    q = 1.0 / cfg['ratio']; k = cfg['k']
    p_audit = q * cfg['audit'] / max(1e-9, 1.0 - cfg['audit'])     # audits are ~`audit` of all full relaxes
    seen = sorted(x for x in hist if x is not None)
    c0 = cpu(); recs = []; kept = []; fulls = screens = 0
    while cpu() - c0 < budget_cpu:
        pid, st0 = pick_start()
        t = cpu()
        try:
            st, lg = S.run_chain(pre, dict(ctx), rng.randrange(1 << 30), st=st0, keep_best=False)
        except Exception as e:
            recs.append(dict(pid=pid, failed=repr(e)[:80])); continue
        if any(r['info'].get('failed') for r in lg) or st.p.n != int(ctx['n']):
            recs.append(dict(pid=pid, failed='build')); continue
        sc = run_screen(st.p, cfg['kmax'])
        screens += 1
        rec = dict(pid=pid if pid is not None else st.meta.get('parent'), raw_s=float(st.p.s), traj=sc['traj'],
                   scr_sec=sc['sec'], scr_cpu=cpu() - t, scr_failed=sc['failed'])
        recs.append(rec)
        s = score_at(sc['traj'], k)
        if sc['failed'] == 'diverged' or s is None:
            continue
        import bisect
        bisect.insort(seen, s)
        cut = seen[max(0, int(q * len(seen)) - 1)] if len(seen) >= 1.0 / q else -math.inf
        if s <= cut:
            how = 'rank'
        elif rng.random() < p_audit:
            how = 'audit'
        else:
            continue
        t = cpu()
        st2, lg = S.run_chain(last, dict(ctx, fq_timeout=ctx['cap']()), rng.randrange(1 << 30), st=st)
        fulls += 1
        info = lg[0]['info']; bad = S.failed(lg)
        rec.update(promoted=how, full_s=float(st2.p.s) if not bad else None, full_cpu=cpu() - t,
                   full_failed=info.get('failed'), deadline=bool(info.get('deadline')), cut=cut,
                   kept_best=lg[-1]['info'].get('src') if lg[-1]['stage'] == 'keep.best' else None)
        if ctx.get('on_relax'):
            ctx['on_relax'](lg)
        if not bad and keep(st2, lg):
            kept.append((float(st2.p.s), rec['pid'], np.asarray(st2.p.sq, float)))
    if log is not None:
        log.update(screened=screens, full_relaxed=fulls, screen_cfg=dict(cfg, serial=True))
    return kept, recs


# ------------------------------------------------------------------ backtest

def _sim(R, k, ratio, budget, rng, full=False):
    """one simulated construction phase over records R (each with full_s, full_cpu, scr_cpu, traj): returns min full_s"""
    best = math.inf; spent = 0.0
    if full:                                               # no screen: build and fully relax every one
        while spent < budget:
            r = R[rng.randrange(len(R))]
            spent += r['full_cpu'] + r['build_cpu']; best = min(best, r['full_s'])
        return best
    pool = []; screens = fulls = 0
    while spent < budget:
        if not pool or screens < ratio * (fulls + 1):
            r = R[rng.randrange(len(R))]
            spent += r['build_cpu'] + r['scr_cpu_k'][k]; screens += 1
            if score_at(r['traj'], k) is not None:
                pool.append(r)
            continue
        i = min(range(len(pool)), key=lambda j: score_at(pool[j]['traj'], k))
        r = pool.pop(i); spent += r['full_cpu']; fulls += 1; best = min(best, r['full_s'])
    return best


def backtest(recs, budgets, ks=(1, 2, 3, 5, 8), ratios=(1, 2, 3, 5, 8, 12), reps=300, seed=0):
    """Equal-CPU comparison of screen policies on records whose builds were all fully relaxed (an audit-1 dataset of
    one arm at one n).  Screen CPU at k = the fq time to outer k (traj) + the fixed per-call overhead.  Returns rows
    (budget, policy, mean best side, P(best side <= the pool's best + 1e-9))."""
    R = [r for r in recs if r.get('full_s') is not None and r.get('traj')]
    if len(R) < 10:
        return []
    for r in R:
        over = max(0.0, (r['scr_cpu'] or 0) - (r['traj'][-1][2] if r['traj'] else 0))
        r['build_cpu'] = 0.0
        r['scr_cpu_k'] = {k: over + next((sec for o, s, sec in reversed(r['traj']) if o <= k), r['traj'][-1][2]) for k in ks}
    top = min(r['full_s'] for r in R)
    out = []
    for B in budgets:
        for pol in [('full', None, None)] + [('screen', k, q) for k in ks for q in ratios]:
            rng = random.Random(seed)
            v = [_sim(R, pol[1], pol[2], B, rng, full=pol[0] == 'full') for _ in range(reps)]
            out.append(dict(budget=B, policy='full' if pol[0] == 'full' else f'k{pol[1]} r{pol[2]}',
                            mean=float(np.mean(v)), p_top=float(np.mean([x <= top + 1e-9 for x in v]))))
    return out


# ------------------------------------------------------------------ report

def load(goal_dir):
    """every job's screen.json under a goal dir -> {(arm kind, n): [records]}"""
    import glob, json, re
    out = {}
    for f in sorted(glob.glob(os.path.join(goal_dir, 'j*', 'screen.json'))):
        arm = json.load(open(os.path.join(os.path.dirname(f), 'arm.json')))
        key = (re.match(r'[a-z_]*', arm['id'].split(':')[1]).group(0), arm['n'])
        out.setdefault(key, []).extend(json.load(open(f)))
    return out


def report(goal_dir, ks=(1, 2, 3, 5, 8), budget_fulls=(3, 6, 12)):
    """per (arm kind, n): screens, promotions, failures; Spearman of the score at k vs the full relax (rank-promoted
    and audited separately: audits are the unbiased sample); with >= 20 fully relaxed builds, the equal-CPU backtest
    at budgets of 3 / 6 / 12 median full relaxes (biased toward screen-promoted builds unless the run used audit 1)."""
    from scipy.stats import spearmanr
    L = []
    for (kind, n), R in sorted(load(goal_dir).items()):
        scr = [r for r in R if 'traj' in r]
        full = [r for r in scr if r.get('full_s') is not None]
        aud = [r for r in full if r.get('promoted') == 'audit']
        fails = {}
        for r in R:
            for f in (r.get('failed'), r.get('scr_failed'), r.get('full_failed')):
                if f: fails[f] = fails.get(f, 0) + 1
        sc = [r['scr_cpu'] for r in scr if r.get('scr_cpu')]
        fc = [r['full_cpu'] for r in full]
        L.append(f'{kind} n={n}: screened {len(scr)} (median {np.median(sc) if sc else 0:.2f} CPU-s), fully relaxed '
                 f'{len(full)} ({len(aud)} audit; median {np.median(fc) if fc else 0:.1f} CPU-s), failures {fails}')
        for name, S_ in (('all relaxed', full), ('audit only', aud)):
            if len(S_) >= 5:
                row = []
                for k in ks:
                    x = [score_at(r['traj'], k) for r in S_]
                    ok = [i for i, v in enumerate(x) if v is not None]
                    if len(ok) >= 5:
                        row.append(f'k{k} {spearmanr([x[i] for i in ok], [S_[i]["full_s"] for i in ok])[0]:+.2f}')
                L.append(f'    spearman ({name}, {len(S_)}): ' + '  '.join(row))
        if len(full) >= 20:
            med = float(np.median(fc))
            rows = backtest(full, [b * med for b in budget_fulls], ks=ks)
            for b in budget_fulls:
                rb = sorted([r for r in rows if abs(r['budget'] - b * med) < 1e-9], key=lambda r: r['mean'])
                base = next(r for r in rb if r['policy'] == 'full')
                L.append(f'    backtest at {b} full relaxes of CPU: full-only mean {base["mean"]:.6f} p_top {base["p_top"]:.2f}; '
                         f'best: ' + ', '.join(f'{r["policy"]} {r["mean"]:.6f} ({r["p_top"]:.2f})' for r in rb[:4]))
    return '\n'.join(L) if L else 'no screen records'
