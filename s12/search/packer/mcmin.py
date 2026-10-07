#!/usr/bin/env python3
"""Metropolis over local minima with replica exchange ("PT on minima", 2026-10-06).

State = a quenched packing (soft penalty squeeze, then slp2).  Energy E = side s (+ line_pen if the minimum has a full line of k
axis squares).  Each replica runs a Metropolis chain at its own temperature T (side units); a proposal applies 1 + Poisson(lam)
random moves, quenches, and accepts with min(1, exp(-(E' - E) / T)).  Between epochs (each replica proposes for `--epoch-sec` seconds in its
own process), adjacent replicas swap states with the usual PT rule.  The quench is a deterministic map, so this is basin hopping
with Metropolis acceptance, not an exact MH chain on configurations: occupancies are protocol-relative.

Moves (all in degrees; content-agnostic unless noted):
  kick       Gaussian sigma on every coordinate (angles 20 sigma deg), sigma log-uniform in [0.003, 0.05]
  lkick      Gaussian sigma in [0.05, 0.25] on the squares within r in [1, 3] of a random point
  reinsert   remove a random square; re-add it at a hole chosen with weight clearance^alpha (holes = local maxima of the
             clearance field: distance to the nearest square or wall), angle 0 / random / nearest neighbour's
  crot       rotate the squares within r in [1, 2.5] of a random point about it by phi (+-2..45 deg, or 90)
  aswap      swap the angles of two neighbouring squares whose angles differ by > 3 deg
  band       movegen.move (structure-aware band moves; s(110)-style angle windows); falls back to kick if no bands

  mcmin.py start.txt [start2 ...] --n 110 --k 11 --out runs/pt1 --replicas 12 --tmin 3e-5 --tmax 1e-2 --hours 2
"""
import argparse, collections, json, math, os, random, shutil, subprocess, sys, tempfile, time
import numpy as np
from jobpool import run_jobs

PACKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'target/release/packer')
H = 0.5
KINDS = dict(kick=2, lkick=2, reinsert=3, crot=2, aswap=1, band=1)
SLP_BUDGET = 60


# ---------- io (file angles in degrees) ----------
def load_deg(path):
    L = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    n, s = int(L[0][0]), float(L[0][1])
    return s, [(float(a), float(b), float(c)) for a, b, c, *_ in L[1:n + 1]]


def write_deg(path, s, sq):
    with open(path, 'w') as f:
        f.write(f'{len(sq)} {s!r}\n')
        for x, y, a in sq:
            f.write(f'{x!r} {y!r} {a % 90!r}\n')


# ---------- holes ----------
def clearance(s, sq, h=0.04):
    """Clearance field on a grid: Euclidean distance to the nearest square (negative inside) or wall."""
    g = np.arange(h / 2, s, h)
    X, Y = np.meshgrid(g, g, indexing='ij')
    P = np.stack([X.ravel(), Y.ravel()], 1)
    r = np.minimum.reduce([P[:, 0], P[:, 1], s - P[:, 0], s - P[:, 1]])
    for x, y, a in sq:
        t = math.radians(a)
        dx, dy = P[:, 0] - x, P[:, 1] - y
        m = (np.abs(dx) < 1.3) & (np.abs(dy) < 1.3)
        if not m.any():
            continue
        u = np.abs(dx[m] * math.cos(t) + dy[m] * math.sin(t)) - H
        v = np.abs(-dx[m] * math.sin(t) + dy[m] * math.cos(t)) - H
        d = np.where((u > 0) | (v > 0), np.hypot(np.maximum(u, 0), np.maximum(v, 0)), np.maximum(u, v))
        r[m] = np.minimum(r[m], d)
    return g, r.reshape(len(g), len(g))


def holes(s, sq, rmin=0.01):
    """Local maxima of the clearance field: [(x, y, r)]."""
    g, R = clearance(s, sq)
    P = np.pad(R, 1, constant_values=-1)
    c = P[1:-1, 1:-1]
    m = c > rmin
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                m &= c >= P[1 + di:P.shape[0] - 1 + di, 1 + dj:P.shape[1] - 1 + dj]
    I, J = np.nonzero(m)
    out = []
    for i, j in sorted(zip(I, J), key=lambda t: -c[t]):
        x, y = g[i], g[j]
        if all(math.hypot(x - a, y - b) > 0.1 for a, b, _ in out):       # merge plateau duplicates
            out.append((float(x), float(y), float(c[i, j])))
    return out


# ---------- holes in pose space ----------
def pose_field(s, sq, h=0.05, dth=3.0):
    """E(x, y, th) = overlap energy of one unit square at that pose against sq: sum of squared SAT penetrations (pair) and
    squared wall protrusions, on a grid (x, y step h over [H, s - H]; th step dth over [0, 90)).  Returns (gx, gth, E)."""
    gx = np.arange(H, s - H + 1e-9, h)
    gth = np.arange(0.0, 90.0, dth)
    nx_, nt = len(gx), len(gth)
    E = np.zeros((nx_, nx_, nt))
    th = np.radians(gth)
    c, sn = np.cos(th), np.sin(th)
    w = H * (np.abs(c) + np.abs(sn))                           # half extent of the rotated square
    wx = np.maximum(w[None, :] - gx[:, None], 0) ** 2 + np.maximum(gx[:, None] + w[None, :] - s, 0) ** 2   # (nx, nt)
    E += wx[:, None, :] + wx[None, :, :]
    for xj, yj, aj in sq:
        i0, i1 = np.searchsorted(gx, [xj - 1.42, xj + 1.42])
        j0, j1 = np.searchsorted(gx, [yj - 1.42, yj + 1.42])
        if i0 >= i1 or j0 >= j1:
            continue
        dx = (xj - gx[i0:i1])[:, None, None]                    # d = B - A
        dy = (yj - gx[j0:j1])[None, :, None]
        tj = math.radians(aj)
        cj, sj = math.cos(tj), math.sin(tj)
        rel = th - tj
        cc = H * (1 + np.abs(np.cos(rel)) + np.abs(np.sin(rel)))      # (nt,)
        p = np.maximum(np.maximum(np.abs(dx * c + dy * sn), np.abs(-dx * sn + dy * c)),
                       np.maximum(np.abs(dx * cj + dy * sj), np.abs(-dx * sj + dy * cj)))
        pen = np.maximum(cc - p, 0.0)
        E[i0:i1, j0:j1, :] += pen * pen
    return gx, gth, E


def pose_holes(s, sq, exclude=None):
    """Local minima of pose_field (th periodic): [(E, x, y, deg)], best first.  exclude = (x, y): drop poses within 0.5."""
    gx, gth, E = pose_field(s, sq)
    m = np.ones(E.shape, bool)
    P = np.pad(E, ((1, 1), (1, 1), (0, 0)), constant_values=np.inf)
    P = np.concatenate([P[:, :, -1:], P, P[:, :, :1]], axis=2)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                if di or dj or dk:
                    m &= E <= P[1 + di:P.shape[0] - 1 + di, 1 + dj:P.shape[1] - 1 + dj, 1 + dk:P.shape[2] - 1 + dk]
    I, J, K = np.nonzero(m)
    out = sorted((float(E[i, j, k]), float(gx[i]), float(gx[j]), float(gth[k])) for i, j, k in zip(I, J, K))
    if exclude is not None:
        out = [o for o in out if math.hypot(o[1] - exclude[0], o[2] - exclude[1]) > 0.5]
    return out


def contact_counts(s, sq):
    """Load-bearing proxy per square: number of tight incidences (g < 1e-6) and corner-corner touches it is in."""
    import inc
    rows, disj = inc.model_fast(s, [(x, y, math.radians(a)) for x, y, a in sq], 1e-6)
    G, I, J, V = rows
    n = len(sq)
    cnt = np.zeros(n)
    sqs = collections.defaultdict(set)
    for r, j in zip(I.tolist(), J.tolist()):
        if j < 3 * n:
            sqs[r].add(j // 3)
    for r, S in sqs.items():
        for i in S:
            cnt[i] += 1
    for (i, j), _ in disj:
        cnt[i] += 1; cnt[j] += 1
    return cnt


def pick_removal(s, sq, rng, a):
    if getattr(a, 'remove', 'uniform') == 'loose':
        c = contact_counts(s, sq)
        return rng.choices(range(len(sq)), weights=[1.0 / (1 + x) ** 2 for x in c])[0]
    return rng.randrange(len(sq))


# ---------- moves ----------
def near(sq, p, r):
    return [i for i, q in enumerate(sq) if math.hypot(q[0] - p[0], q[1] - p[1]) < r]


def mv_kick(s, sq, rng, a):
    sig = math.exp(rng.uniform(math.log(0.003), math.log(0.05)))
    return [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), t + rng.gauss(0, 20 * sig)) for x, y, t in sq], f'kick {sig:.3f}'


def mv_lkick(s, sq, rng, a):
    p = (rng.uniform(0, s), rng.uniform(0, s)); r = rng.uniform(1, 3); sig = rng.uniform(0.05, 0.25)
    idx = set(near(sq, p, r))
    out = [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), t + rng.gauss(0, 20 * sig)) if i in idx else (x, y, t)
           for i, (x, y, t) in enumerate(sq)]
    return out, f'lkick {len(idx)} {sig:.2f}'


def mv_reinsert(s, sq, rng, a):
    if getattr(a, 'reinsert', 'clear') == 'pose':
        return mv_reinsert_pose(s, sq, rng, a)
    hs = holes(s, sq)
    if not hs:
        return mv_kick(s, sq, rng, a)
    w = [h[2] ** a.alpha for h in hs]
    x, y, r = rng.choices(hs, weights=w)[0]
    i = pick_removal(s, sq, rng, a)
    rest = sq[:i] + sq[i + 1:]
    mode = rng.choice(['axis', 'rand', 'nbr'])
    ang = 0.0 if mode == 'axis' else rng.uniform(0, 90) if mode == 'rand' else \
        min(rest, key=lambda q: math.hypot(q[0] - x, q[1] - y))[2]
    return rest + [(x, y, ang)], f'reinsert r={r:.3f} {mode} (of {len(hs)})'


def mv_reinsert_pose(s, sq, rng, a):
    """Remove a square (uniform or load-weighted), re-add it at a local minimum of the pose-space overlap field E(x, y, th)
    of the rest, chosen with weight exp(-(E - E_min) / tau) (tau = 0: the best pose)."""
    i = pick_removal(s, sq, rng, a)
    rest = sq[:i] + sq[i + 1:]
    hs = pose_holes(s, rest, exclude=sq[i][:2])
    if not hs:
        return mv_kick(s, sq, rng, a)
    tau = getattr(a, 'tau', 0.0)
    if tau <= 0:
        E, x, y, ang = hs[0]
    else:
        e0 = hs[0][0]
        E, x, y, ang = rng.choices(hs, weights=[math.exp(-(h[0] - e0) / tau) for h in hs])[0]
    nc = '' if getattr(a, 'remove', 'uniform') != 'loose' else ' loose'
    return rest + [(x, y, ang)], f'reinsert_pose E={E:.4f} th={ang:.0f} (of {len(hs)}){nc}'


def mv_kicksym(s, sq, rng, a):
    """Kick with equal corner displacement from translation and rotation: sigma_th = sigma / (H sqrt 2) rad."""
    sig = math.exp(rng.uniform(math.log(0.003), math.log(0.05)))
    st = math.degrees(sig / (H * math.sqrt(2)))
    return [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), t + rng.gauss(0, st)) for x, y, t in sq], f'kicksym {sig:.3f}'


def mv_crot(s, sq, rng, a):
    p = (rng.uniform(0, s), rng.uniform(0, s)); r = rng.uniform(1, 2.5)
    ph = 90.0 if rng.random() < 0.15 else rng.choice([-1, 1]) * rng.uniform(2, 45)
    c, sn = math.cos(math.radians(ph)), math.sin(math.radians(ph))
    idx = set(near(sq, p, r))
    out = [(p[0] + c * (x - p[0]) - sn * (y - p[1]), p[1] + sn * (x - p[0]) + c * (y - p[1]), t + ph) if i in idx else (x, y, t)
           for i, (x, y, t) in enumerate(sq)]
    return out, f'crot {len(idx)} {ph:+.0f}'


def mv_aswap(s, sq, rng, a):
    d = lambda u, v: min((u - v) % 90, (v - u) % 90)
    for _ in range(50):
        i = rng.randrange(len(sq))
        js = [j for j in near(sq, sq[i], 1.6) if d(sq[j][2], sq[i][2]) > 3]
        if js:
            j = rng.choice(js)
            out = list(sq)
            out[i] = (sq[i][0], sq[i][1], sq[j][2]); out[j] = (sq[j][0], sq[j][1], sq[i][2])
            return out, f'aswap {d(sq[i][2], sq[j][2]):.0f}'
    return mv_kick(s, sq, rng, a)


def mv_band(s, sq, rng, a):
    import movegen
    G = movegen.split(sq)
    L, B = G['L'], G['B']
    if len(L) < 3 or len(B) < 3:
        return mv_kick(s, sq, rng, a)
    p, q = min(((p, q) for p in L for q in B), key=lambda t: math.hypot(t[0][0] - t[1][0], t[0][1] - t[1][1]))
    G = {'L': list(L), 'B': list(B), 'J': list(G['J']) + list(G['O'])}
    desc = movegen.move(G, s, rng, ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2))
    T = G['L'] + G['B'] + G['J']
    out = T + list(movegen.split(sq)['A'])
    # fix the count against the axis squares (as movegen.evaluate_keep)
    from gen import pen
    A = list(movegen.split(sq)['A'])
    over = lambda q, S: sum(max(0.0, pen(q, o)) for o in S)
    while len(T) + len(A) > len(sq):
        A.remove(max(A, key=lambda q: over(q, T)))
    if len(T) + len(A) < len(sq):
        hs = holes(s, T + A) or [(s / 2, s / 2, 0)]
        while len(T) + len(A) < len(sq):
            x, y, _ = hs.pop(0) if hs else (rng.uniform(.5, s - .5), rng.uniform(.5, s - .5), 0)
            A.append((x, y, 0.0))
    return T + A, 'band ' + desc


MOVES = dict(kicksym=mv_kicksym, kick=mv_kick, lkick=mv_lkick, reinsert=mv_reinsert, crot=mv_crot, aswap=mv_aswap, band=mv_band)


# ---------- quench ----------
def quench(s, sq, k, tmp, loosen=1.02, mu0='1e3'):
    """Soft penalty squeeze from a (possibly overlapping) state, then slp2.  Returns (s, sq_deg, info)."""
    from rigid import load as load_rad
    from slp2 import slp2
    from slp import save
    from layout import full_lines
    s1 = s * loosen
    P = [(min(max(x * loosen, H), s1 - H), min(max(y * loosen, H), s1 - H), t) for x, y, t in sq]
    p, q, r = (os.path.join(tmp, f) for f in ('p.txt', 'q.txt', 'r.txt'))
    write_deg(p, s1, P)
    subprocess.run([PACKER, 'relax', '--in', p, '--squeeze-pen', '--mu0', mu0, '--out', q], capture_output=True, timeout=300)
    s2, sqr = load_rad(q)
    if not math.isfinite(s2):
        return None
    info = {}
    s3, sq3, _ = slp2(s2, sqr, R=1e-3, rmin=1e-8, budget=SLP_BUDGET, info=info)
    save(r, s3, sq3)
    info['lines'] = full_lines(r, k)
    return s3, [(x, y, math.degrees(t) % 90) for x, y, t in sq3], info


def energy(s, lines, a):
    return s + (a.line_pen if lines else 0.0)


# ---------- one replica epoch (runs in a worker process) ----------
def epoch(job):
    rid, walker, T, path, E, steps, seed, a = job
    rng = random.Random(seed)
    s, sq = load_deg(path)
    tmp = tempfile.mkdtemp()
    log = []
    kinds, wts = zip(*[(kk, w) for kk, w in a.kinds.items() if w > 0])
    te = time.time() + a.epoch_sec
    st = 0
    while st == 0 or time.time() < te:
        t0 = time.time()
        m = 1 + np.random.default_rng(seed * 1000 + st).poisson(a.lam)
        new, desc = sq, []
        for _ in range(m):
            kd = rng.choices(kinds, weights=wts)[0]
            try:
                new, d = MOVES[kd](s, new, rng, a)
            except Exception as e:                            # a move that cannot apply: skip it
                d = f'{kd} failed {e!r}'[:80]
            desc.append(d)
        res = None
        try:
            res = quench(s, new, a.k, tmp)
        except Exception as e:
            desc.append(f'quench error {e!r}'[:80])
        rec = dict(rid=rid, w=walker, T=T, step=st, moves=desc, dt=round(time.time() - t0, 2))
        if res is None:
            rec.update(acc=False, status='quench failed')
            log.append(rec); st += 1; continue
        s2, sq2, info = res
        E2 = energy(s2, info['lines'], a)
        acc = E2 <= E or (T > 0 and rng.random() < math.exp(-(E2 - E) / T))
        rec.update(s=s2, lines=info['lines'], jam=info['jammed'], timeout=info['timeout'], acc=acc, E_old=E)
        if s2 < a.save_below and not info['lines']:
            fn = os.path.join(a.out, 'min', f'{s2:.10f}.txt')
            if not os.path.exists(fn):
                write_deg(fn, s2, sq2)
            rec['saved'] = fn
        if acc:
            s, sq, E = s2, sq2, E2
        log.append(rec)
        st += 1
    out = os.path.join(a.out, 'state', f'w{walker}.txt')
    write_deg(out + '.tmp', s, sq); os.replace(out + '.tmp', out)
    shutil.rmtree(tmp, ignore_errors=True)
    return dict(path=out, E=E, s=s, log=log)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('start', nargs='+'); ap.add_argument('--n', type=int, required=True); ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--replicas', type=int, default=12)
    ap.add_argument('--tmin', type=float, default=3e-5); ap.add_argument('--tmax', type=float, default=1e-2)
    ap.add_argument('--tinf', action='store_true', help='top replica at T = inf (accept everything)')
    ap.add_argument('--epoch-sec', type=float, default=60, help='each replica proposes until this much wall time has passed')
    ap.add_argument('--sweeps', type=int, default=None, help='swap sweeps per epoch (default: replicas)')
    ap.add_argument('--lam', type=float, default=0.25, help='moves per proposal = 1 + Poisson(lam)')
    ap.add_argument('--alpha', type=float, default=2.0, help='reinsert hole weight = clearance^alpha')
    ap.add_argument('--reinsert', choices=['clear', 'pose'], default='clear', help='holes: clearance maxima / pose-space minima')
    ap.add_argument('--tau', type=float, default=0.0, help='pose reinsert: weight exp(-(E - E_min) / tau); 0 = best pose')
    ap.add_argument('--remove', choices=['uniform', 'loose'], default='uniform', help='reinsert: which square to remove')
    ap.add_argument('--kinds', default=','.join(f'{k}={v}' for k, v in KINDS.items()))
    ap.add_argument('--line-pen', type=float, default=0.0); ap.add_argument('--save-below', type=float, default=None)
    ap.add_argument('--hours', type=float, default=1.0); ap.add_argument('--stop-below', type=float, default=-1)
    ap.add_argument('--procs', type=int, default=None); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--slp-budget', type=float, default=SLP_BUDGET, help='slp2 wall budget per quench (s)')
    a = ap.parse_args()
    globals()['SLP_BUDGET'] = a.slp_budget
    a.kinds = {kv.split('=')[0]: float(kv.split('=')[1]) for kv in a.kinds.split(',')}
    a.save_below = a.k if a.save_below is None else a.save_below
    for d in ('state', 'min'):
        os.makedirs(os.path.join(a.out, d), exist_ok=True)
    R = a.replicas
    Ts = [a.tmin * (a.tmax / a.tmin) ** (i / max(R - 1, 1)) for i in range(R)]
    if a.tinf:
        Ts[-1] = float('inf')
    # walker w starts in replica w; slot[r] = walker at temperature r
    walkers = []
    for w in range(R):
        st = a.start[w % len(a.start)]
        p = os.path.join(a.out, 'state', f'w{w}.txt'); shutil.copy(st, p)
        s, _ = load_deg(p)
        walkers.append(dict(path=p, E=s, s=s, start=st))      # starts are minima without full lines (assumed)
    slot = list(range(R))
    best = min(w['s'] for w in walkers)
    meta = dict(args={k: v for k, v in vars(a).items()}, Ts=Ts, t0=time.time())
    flog = open(os.path.join(a.out, 'log.jsonl'), 'a')
    flog.write(json.dumps(dict(meta=meta)) + '\n')
    lastdir = [+1] * R            # round-trip bookkeeping: direction label per walker (+1 = last visited bottom)
    trips = [0] * R
    t_end = time.time() + 3600 * a.hours
    ep = 0
    while time.time() < t_end:
        jobs = [(r, slot[r], Ts[r], walkers[slot[r]]['path'], walkers[slot[r]]['E'], 0,
                 a.seed * 10 ** 7 + ep * 1000 + r, a) for r in range(R)]
        res = run_jobs(epoch, jobs, procs=a.procs or R, timeout=a.epoch_sec + 4 * SLP_BUDGET + 300,
                       on_timeout=lambda j: None, on_error=lambda j, e: None)
        for r, x in enumerate(res):
            if x is None:
                continue
            w = slot[r]
            walkers[w].update(path=x['path'], E=x['E'], s=x['s'])
            for rec in x['log']:
                flog.write(json.dumps(dict(ep=ep, **rec)) + '\n')
                if 's' in rec and rec['s'] < best - 1e-10 and not rec['lines']:
                    best = rec['s']
                    print(f'  ep {ep} new best {best:.10f} (T={Ts[r]:.1e}, w{w}, {rec["moves"]})', flush=True)
        # replica exchange, alternating even / odd pairs
        sw = []
        for r in (r for sweep in range(a.sweeps or R) for r in range((ep + sweep) % 2, R - 1, 2)):
            wi, wj = slot[r], slot[r + 1]
            Ei, Ej = walkers[wi]['E'], walkers[wj]['E']
            bi, bj = 1 / Ts[r], (0.0 if math.isinf(Ts[r + 1]) else 1 / Ts[r + 1])
            x = (bi - bj) * (Ei - Ej)
            ok = x >= 0 or random.random() < math.exp(x)
            if ok:
                slot[r], slot[r + 1] = wj, wi
            sw.append((r, ok))
            for rr in (0, R - 1):
                w = slot[rr]
                d = +1 if rr == 0 else -1
                if lastdir[w] != d:
                    if d == +1:
                        trips[w] += 1
                    lastdir[w] = d
        flog.write(json.dumps(dict(ep=ep, swaps=sw, slot=slot, E=[walkers[slot[r]]['E'] for r in range(R)], trips=trips)) + '\n')
        flog.flush()
        Es = ' '.join(f'{walkers[slot[r]]["E"]:.5f}' for r in range(R))
        print(f'ep {ep} t={(time.time() - meta["t0"]) / 60:.1f}m best {best:.8f} E: {Es} swaps {sum(o for _, o in sw)}/{len(sw)} props {sum(len(x["log"]) for x in res if x)} '
              f'trips {sum(trips)}', flush=True)
        ep += 1
        if best < a.stop_below:
            print(f'stop: best {best} < {a.stop_below}'); break


if __name__ == '__main__':
    main()
