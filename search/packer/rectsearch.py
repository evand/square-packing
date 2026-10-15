#!/usr/bin/env python3
"""From-scratch search for m unit squares in an (h + d) x h rectangle (PLAN-2026-10-12 item 2; 10-12 night).

Why: the composition gap is piece quality.  Ostashev's decisive piece (38 squares, (h + 2) x h, h = 5.76417431748659,
#528) came from random starts in the rectangle; our best re-optimised crop for that cell is 5.8712 (his best crop
5.8203).  This is a multistart + monotonic basin hopping search with fq's rectangle mode (`--rect-dx d`).

Starts: positions uniform in a loose box (area ~ m x loose), angles: half the squares upright, the rest from a tilt
alphabet (45 deg, atan 1/2, atan 1/3, atan 2/3, 24.4, 36.87, 28.07 deg; one angle per start, +- mirror) or all random.
Each start: fq ALM quench + polish in the rectangle.  Every distinct minimum (h to 1e-7) is logged to minima_<seed>.jsonl
and imported into the store's piece table at the end (pk.pieces; compose2 reads from there).  Then MBH on the elite: gaussian kick (sigma drawn from
{0.02, 0.05, 0.15}) + quench, accept if lower.  Per-cell results: runs/rectsearch/<tag>/m<m>_d<d>/best.txt + log.jsonl.

  python3 rectsearch.py run --cells 38:2 --minutes 10 --procs 8 --out runs/rectsearch/t1
  python3 rectsearch.py run --cells-file cells.json ...       (list of [m, d])
  python3 rectsearch.py report runs/rectsearch/t1
"""
from __future__ import annotations
import argparse, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RELOC = 0.0                      # share of elite moves that relocate 1-3 squares (run --reloc; set before the pool forks)
sys.path.insert(0, HERE)
ALPHABET = [45.0, math.degrees(math.atan(0.5)), math.degrees(math.atan(1 / 3)), math.degrees(math.atan(2 / 3)),
            24.4, 36.87, 28.07]
KICKS = (0.02, 0.05, 0.15)


def h_guess(m, d, loose):
    """h with (h + d) h = m * loose."""
    return (-d + math.sqrt(d * d + 4 * m * loose)) / 2


def random_start(m, d, rng, loose):
    h = h_guess(m, d, loose)
    W = h + d
    x = rng.uniform(0.5, W - 0.5, m); y = rng.uniform(0.5, h - 0.5, m)
    mode = rng.integers(3)
    if mode == 0:
        a = rng.uniform(0, 90, m)
    else:
        t = ALPHABET[rng.integers(len(ALPHABET))]
        t = t if rng.random() < 0.5 else 90 - t
        a = np.where(rng.random(m) < 0.5, 0.0, t)
        if mode == 2:                       # tilted squares in a diagonal band, uprights elsewhere
            c = rng.uniform(-0.3, 0.3) * W
            band = np.abs((x - W / 2) - (y - h / 2) * W / h - c) < rng.uniform(0.15, 0.4) * W
            a = np.where(band, t, 0.0)
    return h, np.c_[x, y, a]


def worker(args):
    m, d, seconds, seed, out, loose, start_frac, a_soft, adaptive, *rest = args
    elite = [(h, np.array(sq, float)) for h, sq in (rest[0] if rest else [])]   # continue from these (bandit episodes)
    sig_ad, h_par = 0.05, math.inf
    from pk.engine import FQ, QOpt
    from pk.packing import Packing
    rng = np.random.default_rng(seed)
    fq = FQ()
    Q = QOpt(loosen=1.0, rect_dx=float(d), flip_top=16)
    t0 = time.time()
    log = open(os.path.join(out, f'log_{seed}.jsonl'), 'w')
    mins, seen = open(os.path.join(out, f'minima_{seed}.jsonl'), 'w'), set()   # every distinct local minimum
    nq = 0
    while time.time() - t0 < seconds:
        if not elite or (time.time() - t0) < start_frac * seconds or rng.random() < 0.2:
            h0, sq = random_start(m, d, rng, loose)
            kind = 'start'
        else:
            k = min(len(elite) - 1, int(rng.exponential(1.5)))
            h0, sq = elite[k]
            h_par = h0                               # the parent's side, before the dilation below
            if RELOC and rng.random() < RELOC:            # relocate 1-3 squares (mostly uprights) to random spots
                sq = sq.copy(); W = h0 + d
                up = np.abs(((sq[:, 2] + 45) % 90) - 45) < 0.5
                kk = int(rng.integers(1, 4))
                pool = np.nonzero(up)[0] if up.any() and rng.random() < 0.7 else np.arange(m)
                idx = rng.choice(pool, size=min(kk, len(pool)), replace=False)
                tilt = sq[~up, 2]
                for i in idx:
                    sq[i, 0] = rng.uniform(0.5, W - 0.5); sq[i, 1] = rng.uniform(0.5, h0 - 0.5)
                    sq[i, 2] = 0.0 if (rng.random() < 0.5 or not len(tilt)) else float(np.median(tilt))
                sq[:, :2] += rng.normal(0, 0.01, (m, 2))
                h0 = h0 * 1.01
                kind = f'reloc{len(idx)}'
            else:
                sig = sig_ad if adaptive else KICKS[rng.integers(len(KICKS))]
                sq = sq.copy(); sq[:, :2] += rng.normal(0, sig, (m, 2)); sq[:, 2] += rng.normal(0, 20 * sig, m)
                h0 = h0 * (1 + sig / 10)
                kind = f'kick{sig:.3g}'
        q = Q
        if kind == 'start' and rng.random() < a_soft:
            q = Q.but(round_r0=float(rng.choice([0.1, 0.2, 0.3, 0.45])), round_compress=True, round_steps=20)
            kind = 'soft'
        r = fq.quench(Packing(float(h0), sq), q, timeout=300)
        nq += 1
        if not r.ok:
            continue
        h, s2 = r.s, np.array(r.p.sq, float)
        log.write(json.dumps(dict(t=round(time.time() - t0, 2), kind=kind, h=h)) + '\n')
        sig = signature(h, s2)
        if sig not in seen:
            seen.add(sig); mins.write(json.dumps(dict(h=h, sq=s2.tolist())) + '\n'); mins.flush()
        if adaptive and kind.startswith('kick'):      # Gensane-Ryckelynck: double on success, halve on failure
            sig_ad = min(0.5, sig_ad * 2) if h < h_par - 1e-9 else max(0.003, sig_ad * 0.5)
        if all(abs(h - e[0]) > 1e-9 for e in elite):
            elite.append((h, s2))
        elite.sort(key=lambda e: e[0]); del elite[12:]
    log.close(); mins.close()
    return dict(seed=seed, nq=nq, best=[(e[0], e[1].tolist()) for e in elite[:12]])


def signature(h, sq):
    """Identity of a local minimum for de-duplication: h to 1e-7.  Mirror / transposed copies and rattler moves share
    it; angle multisets don't work as a key (free rattlers take arbitrary angles: 241 "distinct" minima at m = 20 were
    ~25 values of h).  Distinct minima at the same h to 1e-7 are rare enough to lose."""
    return round(h, 7)


def run(a):
    from pk import manifest
    import scratch
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'rectsearch', args=vars(a), procs=a.procs)
    scratch.scope('rectsearch')
    global RELOC
    RELOC = a.reloc
    cells = [tuple(map(int, c.split(':'))) for c in a.cells.split(',')] if a.cells else \
        [tuple(c) for c in json.load(open(a.cells_file))]
    per = max(1, a.procs // len(cells)) if len(cells) <= a.procs else 1
    init = {}
    for idir in (a.init or []):                          # continue: elite = the best distinct minima of earlier runs
        import glob as _g
        for f in _g.glob(os.path.join(idir, 'm*_d*', 'minima_*.jsonl')):
            mm, dd = map(int, os.path.basename(os.path.dirname(f))[1:].split('_d'))
            for l in open(f):
                e = json.loads(l); init.setdefault((mm, dd), []).append((e['h'], e['sq']))
    for k in init:
        L, out_ = sorted(init[k], key=lambda e: e[0]), []
        for e in L:
            if all(abs(e[0] - x[0]) > 1e-9 for x in out_):
                out_.append(e)
            if len(out_) >= 12:
                break
        init[k] = out_
    tasks = []
    for m, d in cells:
        od = os.path.join(a.out, f'm{m}_d{d}'); os.makedirs(od, exist_ok=True)
        for k in range(per):
            el = init.get((m, d), [])
            if el and a.init_split:                      # worker k starts from the k-th best (diverse chains)
                el = [el[k % len(el)]]
            tasks.append((m, d, a.minutes * 60, 1000 * m + 10 * d + k + a.seed, od, a.loose, a.start_frac, a.soft,
                          a.adaptive, el))
    with ProcessPoolExecutor(a.procs) as ex:
        res = list(ex.map(worker, tasks))
    summ = {}
    for (m, d, *_), r in zip(tasks, res):
        key = f'{m}:{d}'
        b = summ.setdefault(key, dict(m=m, d=d, nq=0, h=math.inf))
        b['nq'] += r['nq']
        if r['best'] and r['best'][0][0] < b['h']:
            b['h'] = r['best'][0][0]
            with open(os.path.join(a.out, f'm{m}_d{d}', 'best.txt'), 'w') as f:
                f.write(f'{m} {b["h"]!r}\n' + ''.join(f'{x!r} {y!r} {t!r}\n' for x, y, t in r['best'][0][1]))
    from pk.pieces import import_paths                  # every distinct minimum -> the store's piece table
    summ['_store'] = dict(zip(('seen', 'new'), import_paths([a.out])))
    json.dump(summ, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    for v in (v for k, v in summ.items() if not k.startswith('_')):
        print(f"m={v['m']} d={v['d']}: best h {v['h']:.10f}  ({v['nq']} quenches)")
    man.finish(summary=summ)


def bandit(a):
    """Cells chosen by Thompson sampling (pk.bandit.CellBandit), one short episode at a time (PLAN P4, 10-13).

    Cells file: JSON list of {m, d, target} (target = the h at which the cell's piece gives a record somewhere; from
    `compose2.py rank` / the probe), or [m, d].  Each cell's elite (12) persists across episodes and is seeded from the
    store's best pieces in that cell (crops and earlier searches: mixed seeds).  Reward of an episode = progress of the
    cell's best h toward its target, (h_before - h_after) / (target gap at the start), capped at 1; cost = CPU seconds.
    A cell that reaches its target retires (and is printed: compose at once).  Stops at --cpu-hours."""
    from pk import manifest
    from pk.bandit import CellBandit
    from pk.store import Store
    import random, scratch
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'rectsearch-bandit', args=vars(a), procs=a.procs)
    scratch.scope('rectsearch')
    C = {}
    S = Store()
    for c in json.load(open(a.cells_file)):
        m, d, tgt = (c['m'], int(c['d']), c.get('target', 0.0)) if isinstance(c, dict) else (c[0], int(c[1]), 0.0)
        el = [(p.s, p.sq.tolist()) for p in S.pieces(m=m, d=float(d))[:a.seed_pieces]]
        od = os.path.join(a.out, f'm{m}_d{d}'); os.makedirs(od, exist_ok=True)
        C[(m, d)] = dict(m=m, d=d, target=tgt, elite=el, h0=el[0][0] if el else math.inf, eps=0, od=od)
    B = CellBandit(half=1e9, uniform=0.15, b0=a.episode)
    rng = random.Random(a.seed)
    budget, used, k = a.cpu_hours * 3600, 0.0, 0
    active = lambda: [c for c, v in C.items() if not (v['elite'] and v['elite'][0][0] <= v['target'] + 1e-9)]
    log = open(os.path.join(a.out, 'episodes.jsonl'), 'a')
    from concurrent.futures import FIRST_COMPLETED, wait
    with ProcessPoolExecutor(a.procs) as ex:
        fut = {}
        def submit():
            nonlocal k
            cells = active()
            if not cells or used + len(fut) * a.episode >= budget:
                return False
            c = B.choose(sorted(cells), rng); v = C[c]; k += 1
            fut[ex.submit(worker, (v['m'], v['d'], a.episode, a.seed * 100000 + k, v['od'], a.loose,
                                   0.1 if v['elite'] else 0.4, a.soft, True, v['elite']))] = (c, time.time())
            return True
        for _ in range(a.procs):
            submit()
        while fut:
            done, _ = wait(fut, return_when=FIRST_COMPLETED)
            for f in done:
                c, t1 = fut.pop(f); v = C[c]; r = f.result(); sec = time.time() - t1; used += sec
                hb = v['elite'][0][0] if v['elite'] else math.inf
                pool = {round(h, 9): (h, sq) for h, sq in v['elite'] + [tuple(e) for e in r['best']]}
                v['elite'] = sorted(pool.values(), key=lambda e: e[0])[:12]; v['eps'] += 1
                ha = v['elite'][0][0] if v['elite'] else math.inf
                gap0 = max(v['h0'] - v['target'], 1e-6) if math.isfinite(v['h0']) else 1.0
                rew = min(1.0, max(0.0, (hb - ha) / gap0)) if math.isfinite(hb) else (1.0 if math.isfinite(ha) else 0.0)
                B.update(c, sec, rew)
                met = ha <= v['target'] + 1e-9
                log.write(json.dumps(dict(cell=c, ep=v['eps'], h_before=hb, h_after=ha, target=v['target'], rew=rew,
                                          sec=round(sec, 1), nq=r['nq'], met=met)) + '\n'); log.flush()
                print(f"  ({c[0]},{c[1]}) ep {v['eps']:3d} h {ha:.9f} target {v['target']:.9f}"
                      f"{'  MET' if met else ''}  used {used / 3600:.2f} CPU-h", flush=True)
                if ha < hb - 1e-12:
                    with open(os.path.join(v['od'], 'best.txt'), 'w') as g:
                        g.write(f"{v['m']} {ha!r}\n" + ''.join(f'{x!r} {y!r} {t!r}\n' for x, y, t in v['elite'][0][1]))
                submit()
    from pk.pieces import import_paths
    seen, new = import_paths([a.out])
    summ = {f'{c[0]}:{c[1]}': dict(h0=v['h0'], h=v['elite'][0][0] if v['elite'] else None, target=v['target'],
                                     eps=v['eps']) for c, v in C.items()}
    json.dump(summ, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    man.finish(summary=dict(cells=len(C), store_new=new, met=[k for k, v in summ.items() if v['h'] and v['h'] <= v['target']]))


def export(a):
    """Every distinct minimum of a run -> <dir>/pieces.jsonl (rectlib's format, for `rectlib.py atlas`).  The store
    already has them (run() imports them); this is only for the file-based tools."""
    from pk.pieces import from_rectsearch
    n = 0
    with open(os.path.join(a.dir, 'pieces.jsonl'), 'w') as f:
        for i, (p, src, ref) in enumerate(from_rectsearch(a.dir)):
            f.write(json.dumps(dict(key=dict(src=f'rs:{p.n}:{p.dx:g}', piece=i, src_n=0), m=p.n, d=p.dx, W0=p.s + p.dx,
                                    H0=p.s, s0=p.s, h=p.s, ok=True, how='search', sq=p.sq.tolist())) + '\n')
            n += 1
    print(f'{n} pieces -> {a.dir}/pieces.jsonl')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run')
    r.add_argument('--cells'); r.add_argument('--cells-file'); r.add_argument('--out', required=True)
    r.add_argument('--minutes', type=float, default=10); r.add_argument('--procs', type=int, default=8)
    r.add_argument('--seed', type=int, default=0); r.add_argument('--loose', type=float, default=1.25)
    r.add_argument('--soft', type=float, default=0.5, help='share of fresh starts that go through soft-to-rigid')
    r.add_argument('--adaptive', action='store_true', help='adaptive kick size (x2 on success, /2 on failure)')
    r.add_argument('--start-frac', type=float, default=0.4, help='share of the time spent on fresh starts only')
    r.add_argument('--init', nargs='*', help='earlier rectsearch out dirs: their best minima seed the elite (continue)')
    r.add_argument('--reloc', type=float, default=0.0, help='share of elite moves that relocate 1-3 squares')
    r.add_argument('--init-split', action='store_true', help='each worker gets one of the init minima (k-th best)')
    b = sp.add_parser('bandit'); b.add_argument('--cells-file', required=True); b.add_argument('--out', required=True)
    b.add_argument('--episode', type=float, default=120); b.add_argument('--cpu-hours', type=float, default=1.0)
    b.add_argument('--procs', type=int, default=14); b.add_argument('--seed', type=int, default=1)
    b.add_argument('--loose', type=float, default=1.25); b.add_argument('--soft', type=float, default=0.5)
    b.add_argument('--seed-pieces', type=int, default=6, help='store pieces per cell seeding the elite')
    e = sp.add_parser('export'); e.add_argument('dir')
    a = ap.parse_args()
    {'run': run, 'bandit': bandit, 'export': export}[a.cmd](a)
