#!/usr/bin/env python3
"""Move-kind bandit timescale study (2026-10-14): how fast do per-kind reward rates drift, and is per-n learning viable?

Reads existing explorer logs only (proposals.jsonl / archive.jsonl / manifest.json under runs/) and the store
(runs/store.sqlite, opened read-only).  Writes nothing except, optionally, a JSON of the tables (--json).
Findings: bandit_timescale.md next to this file.

Rewards recomputed per proposal (all three for every log, whatever the run actually used):
  L  legacy: a new archive entry (new to this explorer round) below k = ceil(sqrt n)          [what chain rewarded]
  G  global: L and the side is not within 2e-9 of any store packing at n added before the run's start, nor of any
     archive entry of an earlier round of the same chain run at n                              [--novel-ref done right]
  V  value: a new entry, globally novel as in G (no k cut), times min(1, 1e-3 / gap), gap = side - best known at n
     (store before the run + earlier rounds of the run), as explore.py --reward value without the robustness factor.
Pre-store logs (before 2026-10-09 20:15: hunt1, hunt1b, hunt2, hunt3/chain, lineage) have no usable reference: L only.

Analyses
  A  half-life scan: one-step-ahead Poisson predictive log-likelihood of each proposal's reward, kind model
     lam_k = (R_k + m) / (C_k + m / lam_all) (shrunk to the pooled rate, m = 2 pseudo-rewards) minus the pooled model
     lam_all = (1 + R) / (60 + C); counts decayed with half-life h in units of proposals / cost-seconds (`sec`, the
     bandit's cost) / wall seconds of round time, carry factor at chain round boundaries in {1, 0.5, 0}.  The gain
     (nats) is the information the kind statistics carry about the next reward; its argmax over h is the timescale.
  B  base-rate drift: reward per cost-second by quarter of each round (all kinds pooled).
  C  pooling: for chain round j >= 2 at n, predict the split of round j's rewards among kinds (multinomial given the
     round's own cost per kind) from (M1) earlier rounds at n (carry 0.5), (M2) other n of the run, earlier in time,
     indirectly standardised, (M3) M2 prior + per-n data with m pseudo-rewards; nats vs no kind information.
  D  kind ranking per reward (indirectly standardised rate ratio, pooled over n), and per-n information.
  E  value of kind choice: held-out reward rate of the policy 1/3 uniform + 2/3 on the estimated best kind, relative
     to uniform, using M1 / M2 / M3 estimates from earlier data and round j's empirical per-kind rates.

Usage (from search/packer):  python3 analysis/bandit_timescale.py [--procs 8] [--json out.json]
"""
from __future__ import annotations
import argparse, collections, datetime as dt, glob, json, math, os, re, sqlite3, sys
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)
RUNS = os.path.join(PK, 'runs')
TOL = 2e-9
G0 = 1e-3
STORE_T0 = dt.datetime(2026, 10, 9, 20, 15)
REW = ('L', 'G', 'V')

# (dataset label, tier, kind, glob of run dirs).  chain: run dir holds r####_n### round dirs; single: run dir is a round.
SOURCES = [
    ('night_nh', 1, 'chain', 'night1014/nh'),
    ('night_nb', 1, 'chain', 'night1014/nb'),
    ('hunt1012', 1, 'chain', 'hunt1012/h*'),
    ('timetravel', 1, 'chain', 'timetravel/*'),
    ('census', 1, 'single', 'census_v[234]/*'),
    ('hunt_old', 2, 'chain', 'hunt1'), ('hunt_old', 2, 'chain', 'hunt1b'), ('hunt_old', 2, 'chain', 'hunt2'),
    ('hunt_old', 2, 'chain', 'hunt3/chain'),
    ('lineage', 2, 'single', 'lineage/*/*'),
]


def parse_t(s):
    return dt.datetime.fromisoformat(s)


def load_round(arg):
    """One explorer round dir -> dict of per-proposal arrays + archive sides + metadata."""
    d, n_hint = arg
    M = {}
    try:
        M = json.load(open(os.path.join(d, 'manifest.json')))
    except Exception:
        pass
    argv = M.get('argv', [])
    n = int(argv[argv.index('--n') + 1]) if '--n' in argv else n_hint
    p = os.path.join(d, 'proposals.jsonl')
    t, kind, sec, cpu, new, s = [], [], [], [], [], []
    for line in open(p):
        try:
            r = json.loads(line)
        except Exception:
            continue
        t.append(r.get('t', 0.0)); kind.append(r['kind']); sec.append(r.get('sec') or 0.0)
        cpu.append(r.get('cpu') if r.get('cpu') is not None else np.nan)
        new.append(r.get('new') is not None); s.append(r['s'] if r.get('s') is not None else np.nan)
    arch = []
    try:
        for line in open(os.path.join(d, 'archive.jsonl')):
            try:
                arch.append(json.loads(line)['s'])
            except Exception:
                pass
    except FileNotFoundError:
        pass
    if M.get('start'):
        start = parse_t(M['start'])
    else:                                   # pre-manifest logs: file mtime minus the last log time
        start = dt.datetime.fromtimestamp(os.path.getmtime(p) - (t[-1] if t else 0.0))
    return dict(dir=d, n=n, start=start, t=np.array(t), kind=kind, sec=np.array(sec), cpu=np.array(cpu),
                new=np.array(new, bool), s=np.array(s), arch=np.array(arch), adapt=('--adapt' in argv) if argv else None)


def discover():
    jobs = []                                # (label, tier, run_id, round_dir, n_hint)
    for label, tier, typ, pat in SOURCES:
        for rd in sorted(glob.glob(os.path.join(RUNS, pat))):
            if not os.path.isdir(rd):
                continue
            if typ == 'chain':
                for d in sorted(glob.glob(os.path.join(rd, 'r[0-9]*_n[0-9]*'))):
                    if re.search(r'_n\d+$', d) and os.path.exists(os.path.join(d, 'proposals.jsonl')):
                        jobs.append((label, tier, rd, d, int(re.search(r'_n(\d+)$', d).group(1))))
            elif os.path.exists(os.path.join(rd, 'proposals.jsonl')):
                jobs.append((label, tier, rd, rd, 110))
    return jobs


def load_store(ns):
    con = sqlite3.connect(f'file:{os.path.join(RUNS, "store.sqlite")}?mode=ro', uri=True)
    out = {}
    q = ("select s, added from packing where n = ? and status not in ('infeasible', 'unresolved')")
    for n in ns:
        rows = con.execute(q, (n,)).fetchall()
        out[n] = (np.array([r[0] for r in rows]), np.array([parse_t(r[1]).timestamp() for r in rows]))
    con.close()
    return out


def near(ref, s):
    """bool array: s within TOL of some element of sorted ref."""
    if len(ref) == 0:
        return np.zeros(len(s), bool)
    i = np.searchsorted(ref, s - TOL)
    i2 = np.minimum(i, len(ref) - 1)
    return (i < len(ref)) & (np.abs(ref[i2] - s) < TOL)


def build_streams(jobs, rounds, store):
    """Group rounds into streams (run, n), order by start, attach rewards."""
    by = collections.defaultdict(list)
    for (label, tier, rid, d, _), R in zip(jobs, rounds):
        if R is None or len(R['t']) == 0:
            continue
        by[(label, tier, rid, R['n'])].append(R)
    run_start = {}
    for (label, tier, rid, n), Rs in by.items():
        run_start[rid] = min(run_start.get(rid, Rs[0]['start']), min(R['start'] for R in Rs))
    streams = []
    for (label, tier, rid, n), Rs in sorted(by.items(), key=lambda kv: str(kv[0])):
        Rs.sort(key=lambda R: R['start'])
        k = math.ceil(math.sqrt(n) - 1e-12)
        T0 = run_start[rid].timestamp()
        if tier == 1:
            ss, sa = store.get(n, (np.zeros(0), np.zeros(0)))
            ref = np.sort(ss[sa < T0])
        else:
            ref = np.zeros(0)
        best = ref.min() if len(ref) else np.inf
        for R in Rs:
            s = R['s']; isnew = R['new'] & np.isfinite(s)
            R['L'] = (isnew & (s < k)).astype(float)
            if tier == 1:
                nov = isnew & ~near(ref, np.nan_to_num(s, nan=-1.0))
                R['G'] = (nov & (s < k)).astype(float)
                gap = np.maximum(np.nan_to_num(s, nan=np.inf) - best, 1e-12)
                R['V'] = np.where(nov, np.minimum(1.0, G0 / gap), 0.0)
            else:
                R['G'] = R['V'] = np.full(len(s), np.nan)
            R['ref_size'] = len(ref)
            ref = np.sort(np.concatenate([ref, R['arch'][np.isfinite(R['arch'])]]))
            if len(R['arch']):
                best = min(best, np.nanmin(R['arch']))
        streams.append(dict(label=label, tier=tier, run=rid, n=n, rounds=Rs))
    return streams


# ---------------------------------------------------------------- A: half-life scan (worker globals via fork)
FLAT = []                                    # per stream: (kind idx list, sec, prop, cost, wall, rnd, Y[3][N], K)


BASE_H = 200.0                               # base (all-kind) rate: half-life in proposals, x0.5 at round starts, prior 5 rewards


def base_expect(sec, rnd, y, prior_rate):
    """expected reward of each proposal under the decayed all-kind rate (no kind information), computed online."""
    f = 0.5 ** (1.0 / BASE_H); R = 5.0; C = 5.0 / prior_rate; lr = rnd[0]; out = np.empty(len(sec))
    for i in range(len(sec)):
        if rnd[i] != lr:
            lr = rnd[i]; R *= 0.5; C *= 0.5
        out[i] = sec[i] * (R + 0.5) / (C + 0.5 / prior_rate)
        R = R * f + y[i]; C = C * f + sec[i]
    return out


def flatten(streams):
    out = []
    for S in streams:
        kinds = sorted({k for R in S['rounds'] for k in R['kind']})
        ki = {k: i for i, k in enumerate(kinds)}
        kk, sec, wall, rnd, Y = [], [], [], [], [[], [], []]
        off = 0.0
        for j, R in enumerate(S['rounds']):
            kk += [ki[k] for k in R['kind']]; sec.append(R['sec']); wall.append(off + R['t']); rnd.append(np.full(len(R['t']), j))
            off += (R['t'][-1] if len(R['t']) else 0.0) + 1.0
            for q, r in enumerate(REW):
                Y[q].append(R[r])
        sec = np.concatenate(sec); rnd = np.concatenate(rnd)
        Y = [np.nan_to_num(np.concatenate(y)) for y in Y]
        nq = 1 if S['tier'] == 2 else 3
        E = []
        for q in range(nq):
            pr = max(Y[q].sum(), 1.0) / max(sec.sum(), 1e-9)
            E.append(base_expect(sec.tolist(), rnd.tolist(), Y[q].tolist(), pr).tolist())
        out.append(dict(kind=kk, sec=sec.tolist(), prop=list(range(len(kk))), cost=np.cumsum(sec).tolist(),
                        wall=np.concatenate(wall).tolist(), rnd=rnd.tolist(), E=E,
                        Y=[y.tolist() for y in Y[:nq]], K=len(kinds), tier=S['tier']))
    return out


def scan_one(cfg):
    """kind-information gain: Poisson log-lik of y_i under rho_k * E_i minus under E_i, where E_i is the all-kind base
    expectation and rho_k = (R_k + m) / (E_k + m) from the kind's decayed rewards / base expectations (m = 2)."""
    unit, h, carry = cfg
    M_R = 2.0
    res = []                                  # per stream: [gain_q, base ll_q, nrew_q]
    for F in FLAT:
        K = F['K']; kind = F['kind']; sec = F['sec']; tv = F[unit]; rnd = F['rnd']; Y = F['Y']; E = F['E']
        nq = len(Y)
        R = [[0.0] * K for _ in range(nq)]; X = [[0.0] * K for _ in range(nq)]
        T0 = tv[0]; lr = rnd[0]
        inv = 0.0 if h == math.inf else 1.0 / h
        gain = [0.0] * nq; llb = [0.0] * nq; nr = [0.0] * nq
        for i in range(len(kind)):
            T = tv[i]
            fac = 1.0
            if rnd[i] != lr:                    # chain round boundary: carry
                lr = rnd[i]; fac = carry
            x = (T - T0) * inv
            if x > 600.0:                       # renormalise the inflated sums (stored value = true value * 2^x)
                fac *= 2.0 ** (-x); T0 = T; x = 0.0
            if fac != 1.0:
                for q in range(nq):
                    Rq = R[q]; Xq = X[q]
                    for j in range(K):
                        Rq[j] *= fac; Xq[j] *= fac
            g = 2.0 ** (-x)
            k = kind[i]
            if sec[i] > 0:
                w = 1.0 / g
                for q in range(nq):
                    e = E[q][i]; y = Y[q][i]
                    if e <= 0:
                        continue
                    rho = (R[q][k] * g + M_R) / (X[q][k] * g + M_R)
                    gain[q] += y * math.log(rho) - (rho - 1.0) * e
                    llb[q] += y * math.log(e) - e; nr[q] += y
                    R[q][k] += y * w; X[q][k] += e * w
        res.append((gain, llb, nr))
    return cfg, res


HGRID = {'prop': [30, 100, 300, 1000, 3000, 10000, 30000, 100000, math.inf],
         'cost': [10, 30, 100, 300, 1000, 3000, 10000, 30000, math.inf],
         'wall': [10, 30, 100, 300, 600, 1000, 3000, 10000, math.inf]}
CARRY = [1.0, 0.5, 0.0]


def nbucket(n):
    return '<70' if n < 70 else '70-129' if n < 130 else '130-199' if n < 200 else '>=200'


# ---------------------------------------------------------------- C/D/E: round-level kind effects
def round_tables(S, kinds):
    """per round: cost per kind (C), reward per kind (Y[q])."""
    ki = {k: i for i, k in enumerate(kinds)}
    out = []
    for R in S['rounds']:
        idx = np.array([ki[k] for k in R['kind']])
        C = np.bincount(idx, weights=R['sec'], minlength=len(kinds))
        Y = {q: np.bincount(idx, weights=np.nan_to_num(R[q]), minlength=len(kinds)) for q in REW}
        out.append(dict(start=R['start'], C=C, Y=Y, nprop=np.bincount(idx, minlength=len(kinds))))
    return out


def smr(pairs):
    """indirectly standardised kind rate ratio: sum R_k / sum E_k, E_k = C_k * (R / C) of the same unit."""
    num = 0.0; den = 0.0
    for C, Y in pairs:
        tot = C.sum()
        if tot <= 0 or Y.sum() <= 0:
            continue
        num = num + Y; den = den + C * (Y.sum() / tot)
    return num, den


def mlogl(Y, C, rho):
    """multinomial log-lik of the split Y given cost C and relative rates rho (kinds with C = 0 dropped)."""
    m = C > 0
    p = rho[m] * C[m]; p = p / p.sum()
    return float((Y[m] * np.log(np.maximum(p, 1e-300))).sum())


def pooling(streams, kinds, q, mgrid=(1, 3, 10, 30, 100)):
    """C + E for chain streams of tier 1 (q in G, V) or all chains (L)."""
    chains = [S for S in streams if S['label'] not in ('census', 'lineage') and (q == 'L' or S['tier'] == 1)]
    tabs = {id(S): round_tables(S, kinds) for S in chains}
    byrun = collections.defaultdict(list)
    for S in chains:
        byrun[S['run']].append(S)
    acc = collections.defaultdict(float)
    ev = collections.defaultdict(lambda: [0.0, 0.0])          # policy -> [held-out reward of policy, of uniform] (rate x cost)
    eq = collections.defaultdict(list)
    nrounds = 0; nrew = 0.0
    for S in chains:
        T = tabs[id(S)]
        for j in range(0, len(T)):
            Yj = T[j]['Y'][q]; Cj = T[j]['C']
            if Yj.sum() <= 0 or (Cj > 0).sum() < 2:
                continue
            # M1: earlier rounds at n, carry 0.5 per round
            num1 = np.zeros(len(kinds)); den1 = np.zeros(len(kinds))
            for a, Tr in enumerate(T[:j]):
                w = 0.5 ** (j - 1 - a)
                if Tr['C'].sum() > 0 and Tr['Y'][q].sum() > 0:
                    num1 += w * Tr['Y'][q]; den1 += w * Tr['C'] * Tr['Y'][q].sum() / Tr['C'].sum()
            # M2: other n of the same run, rounds started before this one
            pairs = [(Tr['C'], Tr['Y'][q]) for S2 in byrun[S['run']] if S2 is not S
                     for Tr in tabs[id(S2)] if Tr['start'] < T[j]['start']]
            num2, den2 = smr(pairs)
            if np.isscalar(num2):
                continue
            rho2 = (num2 + 1.0) / (den2 + 1.0)
            # M2 near: other n within x1.25 of n (same run, earlier in time)
            pn = [(Tr['C'], Tr['Y'][q]) for S2 in byrun[S['run']] if S2 is not S and abs(math.log(S2['n'] / S['n'])) < math.log(1.25)
                  for Tr in tabs[id(S2)] if Tr['start'] < T[j]['start']]
            numn, denn = smr(pn)
            rhon = rho2 if np.isscalar(numn) else (numn + 10 * rho2) / (denn + 10)
            # M4: other runs only (offline fixed weights from history)
            po = [(Tr['C'], Tr['Y'][q]) for S2 in chains if S2['run'] != S['run'] for Tr in tabs[id(S2)]]
            num4, den4 = smr(po)
            rho4 = rho2 if np.isscalar(num4) else (num4 + 1.0) / (den4 + 1.0)
            ests = {'M1 per-n': (num1 + 1.0) / (den1 + 1.0), 'M2 pooled': rho2, 'M2 near-n': rhon, 'M4 other runs': rho4}
            for m in mgrid:
                ests[f'M3 m={m}'] = (num1 + m * rho2) / (den1 + m)
            for m in (3, 10, 30):
                ests[f'M3 near m={m}'] = (num1 + m * rhon) / (den1 + m)
            base = mlogl(Yj, Cj, np.ones(len(kinds)))
            for name, rho in ests.items():
                acc[name] += mlogl(Yj, Cj, rho) - base
            # E: policy value on round j (kinds with cost in round j only)
            mk = Cj > 0
            lam = np.where(mk, Yj / np.where(mk, Cj, 1.0), 0.0)
            u = mk / mk.sum()
            for name, rho in list(ests.items()) + [('oracle (in-sample)', np.where(mk, lam, 0.0)), ('realized (logged)', None)]:
                if rho is None:                              # the allocation the logged bandit actually made in round j
                    w = Cj / Cj.sum()
                else:
                    r_ = np.where(mk, rho, -np.inf)
                    w = u / 3.0; w = w.copy(); w[int(np.argmax(r_))] += 2.0 / 3.0
                ev[name][0] += float((w * lam).sum()) * Cj.sum()
                ev[name][1] += float((u * lam).sum()) * Cj.sum()
                if (u * lam).sum() > 0:                      # per-round ratio, rounds weighted equally (n-balanced)
                    eq[name].append(float((w * lam).sum() / (u * lam).sum()))
            nrounds += 1; nrew += Yj.sum()
    return dict(nats={k: v for k, v in acc.items()}, rounds=nrounds, rewards=nrew,
                policy={k: v[0] / v[1] for k, v in ev.items() if v[1] > 0},
                policy_eq={k: (float(np.mean(v)), float(np.std(v) / math.sqrt(len(v)))) for k, v in eq.items()})


def dl_tau(groups):
    """DerSimonian-Laird between-unit SD from groups of (estimates, variances), each group centred on its own mean."""
    Q = df = den = 0.0
    for y, v in groups:
        if len(y) < 2:
            continue
        w = 1.0 / v; yb = (w * y).sum() / w.sum()
        Q += (w * (y - yb) ** 2).sum(); df += len(y) - 1; den += w.sum() - (w ** 2).sum() / w.sum()
    if den <= 0:
        return None
    return float(math.sqrt(max(0.0, (Q - df) / den)))


def per_n_info(streams, kinds, q, label='night_nh'):
    """per-n kind log rate ratios after all rounds: sampling SE vs between-n SD (method of moments) vs between-kind SD."""
    rows = []
    for S in streams:
        if S['label'] != label:
            continue
        T = round_tables(S, kinds)
        C = sum(t['C'] for t in T); Y = sum(t['Y'][q] for t in T)
        if Y.sum() <= 0:
            continue
        E = C * Y.sum() / C.sum()
        rows.append((S['n'], len(T), Y, E, C, T))
    if not rows:
        return {}
    Ys = np.array([r[2] for r in rows]); Es = np.array([r[3] for r in rows])
    num = Ys.sum(0); den = Es.sum(0)
    rho = num / np.maximum(den, 1e-12)
    out = {'n_count': len(rows), 'rounds_median': float(np.median([r[1] for r in rows])),
           'reward_per_n_median': float(np.median(Ys.sum(1))), 'kinds': {}}
    for j, k in enumerate(kinds):
        ok = (Es[:, j] > 0.5)
        if ok.sum() < 5 or num[j] <= 0:
            continue
        y = np.log((Ys[ok, j] + 0.5) / (Es[ok, j] + 0.5))
        v = 1.0 / (Ys[ok, j] + 0.5) + 1.0 / (Ys[ok].sum(1) + 0.5)    # Poisson var of the log ratio
        tau_n = dl_tau([(y, v)])                                     # between-n (DerSimonian-Laird)
        groups = []                                                  # between rounds at the same n (drift + extra noise)
        for r in rows:
            T = r[5]
            yy = []; vv = []
            for t in T:
                Yt = t['Y'][q]; Ct = t['C']
                if Ct[j] <= 0 or Yt.sum() <= 0:
                    continue
                Et = Ct[j] * Yt.sum() / Ct.sum()
                yy.append(math.log((Yt[j] + 0.5) / (Et + 0.5))); vv.append(1.0 / (Yt[j] + 0.5) + 1.0 / (Yt.sum() + 0.5))
            if len(yy) >= 2:
                groups.append((np.array(yy), np.array(vv)))
        tau_r = dl_tau(groups)
        se2 = 1.0 / (Ys[ok, j] + 0.5)
        out['kinds'][k] = dict(rho=float(rho[j]), med_rew_per_n=float(np.median(Ys[ok, j])), se_med=float(np.sqrt(np.median(se2))),
                               tau=tau_n, tau_round=tau_r, n=int(ok.sum()))
    ks = [k for k in out['kinds']]
    lrho = np.log([out['kinds'][k]['rho'] for k in ks])
    out['between_kind_sd'] = float(lrho.std())
    out['tau_median'] = float(np.median([out['kinds'][k]['tau'] for k in ks]))
    out['tau_round_median'] = float(np.median([out['kinds'][k]['tau_round'] for k in ks if out['kinds'][k]['tau_round'] is not None]))
    out['se_median'] = float(np.median([out['kinds'][k]['se_med'] for k in ks]))
    # n-trend: correlation of per-n log ratio with log n, per kind
    tr = {}
    for j, k in enumerate(kinds):
        ok = (Es[:, j] > 0.5)
        if ok.sum() < 10:
            continue
        lr = np.log((Ys[ok, j] + 0.5) / (Es[ok, j] + 0.5)); ln = np.log([r[0] for r, o in zip(rows, ok) if o])
        tr[k] = float(np.polyfit(ln, lr, 1)[0])
    out['slope_per_log_n'] = tr
    return out


def ranking(streams, kinds, q, labels):
    pairs = []; Ctot = np.zeros(len(kinds))
    for S in streams:
        if S['label'] not in labels or (q != 'L' and S['tier'] != 1):
            continue
        for t in round_tables(S, kinds):
            pairs.append((t['C'], t['Y'][q])); Ctot += t['C']
    num, den = smr(pairs)
    if np.isscalar(num):
        return {}
    return {k: dict(rho=float(num[j] / den[j]) if den[j] > 0 else None, rew=float(num[j]), cpu_share=float(Ctot[j] / Ctot.sum()))
            for j, k in enumerate(kinds) if Ctot[j] > 0}


def base_drift(streams, q):
    """reward per cost-s by quarter of round time (pooled over kinds), normalised to the first quarter; per dataset."""
    out = {}
    for lab in sorted({S['label'] for S in streams}):
        num = np.zeros(4); den = np.zeros(4)
        for S in streams:
            if S['label'] != lab or (q != 'L' and S['tier'] != 1):
                continue
            for R in S['rounds']:
                if len(R['t']) < 8:
                    continue
                qt = np.minimum((4 * R['t'] / (R['t'][-1] + 1e-9)).astype(int), 3)
                num += np.bincount(qt, weights=np.nan_to_num(R[q]), minlength=4)
                den += np.bincount(qt, weights=R['sec'], minlength=4)
        if den.sum() > 0 and num[0] > 0:
            r = num / np.maximum(den, 1e-12)
            out[lab] = [float(x / r[0]) for x in r] + [float(num.sum())]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--json', help='write all tables here')
    ap.add_argument('--no-scan', action='store_true', help='skip the half-life scan (A)')
    a = ap.parse_args()
    jobs = discover()
    print(f'{len(jobs)} round dirs', flush=True)
    with Pool(a.procs) as P:
        rounds = P.map(load_round, [(j[3], j[4]) for j in jobs], chunksize=4)
    store = load_store(sorted({R['n'] for R in rounds}))
    streams = build_streams(jobs, rounds, store)
    del rounds
    kinds = sorted({k for S in streams for R in S['rounds'] for k in R['kind']})
    out = dict(kinds=kinds)

    # inventory
    inv = collections.defaultdict(lambda: collections.Counter())
    for S in streams:
        c = inv[S['label']]
        c['streams'] += 1; c['rounds'] += len(S['rounds'])
        for R in S['rounds']:
            c['props'] += len(R['t']); c['sec'] += float(R['sec'].sum()); c['cpu'] += float(np.nansum(R['cpu']))
            c['cpu_logged'] += float(np.isfinite(R['cpu']).sum())
            for q in REW:
                c[q] += float(np.nansum(R[q]))
    out['inventory'] = {k: dict(v) for k, v in inv.items()}
    print('\n== inventory (rewards = events / value sum)')
    for k, v in out['inventory'].items():
        print(f"{k:11s} streams {v['streams']:4d} rounds {v['rounds']:4d} props {v['props']:8d} sec {v['sec'] / 3600:6.1f} h "
              f"cpu/sec {v['cpu'] / max(1e-9, v['sec']) if v['cpu_logged'] else float('nan'):.2f} "
              f"L {v['L']:8.0f} ({v['L'] / v['props']:.1%})  G {v['G']:7.0f} ({v['G'] / v['props']:.2%})  V {v['V']:7.1f}")
    # proposals per round vs n (night_nh)
    pr = [(S['n'], len(R['t']), float(R['sec'].mean())) for S in streams if S['label'] == 'night_nh' for R in S['rounds']]
    out['night_props_per_round'] = {b: [float(np.median([p for n, p, _ in pr if nbucket(n) == b])), float(np.median([s for n, _, s in pr if nbucket(n) == b]))]
                                    for b in ('<70', '70-129', '130-199', '>=200')}
    print('night_nh median proposals / round, sec / proposal by n bucket:', out['night_props_per_round'])

    out['base_drift'] = {q: base_drift(streams, q) for q in REW}
    print('\n== B: pooled reward rate by quarter of round time (rel. to Q1; last = events)')
    for q in REW:
        for lab, v in out['base_drift'][q].items():
            print(f'  {q} {lab:11s}', ' '.join(f'{x:5.2f}' for x in v[:4]), f'({v[4]:.0f})')

    out['ranking'] = {}
    print('\n== D: kind rate ratio (indirectly standardised; 1 = average kind at that n/round), CPU share')
    for labs, name in ((('night_nh', 'night_nb', 'hunt1012', 'timetravel'), 'chain'), (('census',), 'census'), (('hunt_old', 'lineage'), 'old')):
        out['ranking'][name] = {q: ranking(streams, kinds, q, labs) for q in REW}
        print(f' [{name}]')
        for k in kinds:
            row = [out['ranking'][name][q].get(k) for q in REW]
            if not row[0]:
                continue
            print(f'  {k:10s} cpu {row[0]["cpu_share"]:5.1%}  ' + '  '.join(
                f'{q} {r["rho"]:5.2f} ({r["rew"]:7.0f})' if r and r['rho'] is not None else f'{q}   -  ' for q, r in zip(REW, row)))
    try:
        from scipy.stats import spearmanr
        for name in out['ranking']:
            R_ = out['ranking'][name]
            ks = [k for k in kinds if all(R_[q].get(k) and R_[q][k]['rho'] for q in REW)]
            if len(ks) > 3:
                for q2 in ('G', 'V'):
                    rr = spearmanr([R_['L'][k]['rho'] for k in ks], [R_[q2][k]['rho'] for k in ks])[0]
                    print(f'  spearman L vs {q2} [{name}]: {rr:.2f}')
                    out.setdefault('spearman', {})[f'{name}:L-{q2}'] = float(rr)
    except ImportError:
        pass

    out['per_n_info'] = {q: per_n_info(streams, kinds, q) for q in REW}
    print('\n== D2: night_nh per-n information (all rounds at n)')
    for q in REW:
        I = out['per_n_info'][q]
        if not I:
            continue
        print(f"  {q}: n {I['n_count']}, rounds/n {I['rounds_median']:.0f}, rewards/n {I['reward_per_n_median']:.0f}; "
              f"between-kind SD of log ratio {I['between_kind_sd']:.2f}, between-n SD (tau, median kind) {I['tau_median']:.2f}, "
              f"between rounds at same n {I['tau_round_median']:.2f}, sampling SE per (n, kind) {I['se_median']:.2f}")
        print('     per-kind tau (between n / between rounds):', {k: (round(v['tau'], 2), v['tau_round'] and round(v['tau_round'], 2)) for k, v in I['kinds'].items()})
        print('     slope of log ratio per log n:', {k: round(v, 2) for k, v in I['slope_per_log_n'].items()})

    out['pooling'] = {q: pooling(streams, kinds, q) for q in REW}
    print('\n== C/E: predict round j >= 2 kind split (nats vs no kind info), policy value (1/3 uniform + 2/3 best) vs uniform')
    for q in REW:
        P_ = out['pooling'][q]
        print(f"  {q}: {P_['rounds']} rounds, {P_['rewards']:.0f} rewards")
        for k, v in P_['nats'].items():
            print(f"     {k:14s} {v:9.1f} nats ({1000 * v / max(P_['rewards'], 1):6.1f} per 1000 rewards)  policy x{P_['policy'].get(k, float('nan')):.3f} (rounds equal: x{P_['policy_eq'][k][0]:.3f} +- {P_['policy_eq'][k][1]:.3f})")
        print(f"     realized (logged bandit allocation) x{P_['policy']['realized (logged)']:.3f} (rounds equal x{P_['policy_eq']['realized (logged)'][0]:.3f})")
        print(f"     oracle (in-sample) policy x{P_['policy'].get('oracle (in-sample)', float('nan')):.3f} (rounds equal x{P_['policy_eq']['oracle (in-sample)'][0]:.3f})")

    if not a.no_scan:
        global FLAT
        FLAT = flatten(streams)
        meta = [(S['label'], S['n'], S['tier']) for S in streams]
        cfgs = [(u, h, c) for u in HGRID for h in HGRID[u] for c in CARRY]
        with Pool(a.procs) as P:
            res = P.map(scan_one, cfgs, chunksize=1)
        print('\n== A: half-life scan: kind-information gain (nats per 1000 rewards), by dataset group / n bucket')
        groups = {}
        for i, (lab, n, tier) in enumerate(meta):
            for g in ((lab if lab in ('census', 'lineage') else 'chain:' + nbucket(n)) if tier == 1 else 'old:' + (lab if lab == 'lineage' else nbucket(n)),):
                groups.setdefault(g, []).append(i)
        table = {}
        for cfg, rs in res:
            for g, idx in groups.items():
                for q in range(3):
                    vals = [rs[i] for i in idx if len(rs[i][0]) > q]
                    if not vals:
                        continue
                    G_ = sum(v[0][q] for v in vals); N_ = sum(v[2][q] for v in vals); LP = sum(v[1][q] for v in vals)
                    table[(g, REW[q], cfg)] = (G_, N_, LP)
        out['scan'] = {f'{g}|{q}|{c[0]}|{c[1]}|{c[2]}': v for (g, q, c), v in table.items()}
        for g in sorted(groups):
            for q in REW:
                keys = [(c, v) for (gg, qq, c), v in table.items() if gg == g and qq == q]
                if not keys or keys[0][1][1] <= 0:
                    continue
                N_ = keys[0][1][1]
                cur = table.get((g, q, ('wall', 600, 0.5)))
                inf_ = table.get((g, q, ('wall', math.inf, 1.0)))
                line = f'  {g:16s} {q} rewards {N_:8.0f} | current(wall600,x.5) {1000 * cur[0] / N_:6.1f}  no-decay {1000 * inf_[0] / N_:6.1f} |'
                for u in HGRID:
                    for c in CARRY:
                        best = max(((cc, v) for cc, v in keys if cc[0] == u and cc[2] == c), key=lambda z: z[1][0])
                        line += f' {u[0]}{c:g}: h*={best[0][1]:g} {1000 * best[1][0] / N_:5.1f}'
                print(line)
            # curve for the main unit
        print('\n  gain curves (nats / 1000 rewards), carry 0.5:')
        for g in sorted(groups):
            for q in REW:
                for u in HGRID:
                    vals = [table.get((g, q, (u, h, 0.5))) for h in HGRID[u]]
                    if not vals[0] or vals[0][1] <= 0:
                        continue
                    print(f'   {g:16s} {q} {u:4s} ' + ' '.join(f'{h:g}:{1000 * v[0] / v[1]:.1f}' for h, v in zip(HGRID[u], vals)))
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
