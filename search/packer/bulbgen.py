#!/usr/bin/env python3
"""Band / bulb start generator (template, 10-19; notes: BULB.md).  Evan: tilted-band constructions for c*(k) (Ryu's L,
our R1') and the two-bulb records near k^2 (s(308)) are one family: wall bands of b-wide tilted rows at about
2 arctan(1/b), with the ends tapered and the far corner handled as an L junction, a V bend or a plain grid block.

Every parameter is a distribution (`sample_params`), not a constant: one call = one random start.  Pipeline per start:
  sample params -> build arms in the (k - delta) box -> upright fill + top-up placer (pk/fill.py, 10-20) to a target
  count n + extra (extra in -2..2), or until full -> optional fq quench.  The output count is a distribution: each
  start file is named with the count it reached (s<seed>_n<count>.txt) and is a start for that n.
Output is a distribution of starts for refinement (CLAUDE.md rule 10-18), never "one packing".

Status: TEMPLATE.  Geometry and fill run; distributions are first guesses; no evaluation yet (BULB.md, "Next").

  python3 bulbgen.py sample --k 18 --n 308 --count 20 --out runs/bulb/t308 [--seed 0] [--quench shrink]
"""
from __future__ import annotations
import argparse, json, math, os, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen  # noqa: E402  (pen, wall_viol, axis_fill, write)

TOPOLOGIES = ('L', 'V', 'bulb')


# ---------------------------------------------------------------- parameters (all random)

PRIORS = ('broad', 'record')


def sample_params(rng, k, n, prior='broad'):
    """One draw of every structural parameter.  Edit the distributions here; keep each one explicit.
    prior 'record' (10-21, 274 campaign): centred on the two-bulb records next to the target (272 / 273 / 307 / 308:
    6-row bands at ~18.5-19 deg = theta0(6) 18.9, balanced, mostly the bulb layout, side just below k)."""
    if prior == 'record':
        b1 = int(rng.choice([5, 6, 7], p=[0.2, 0.6, 0.2]))
        b2 = int(b1 if rng.random() < 0.5 else np.clip(b1 + rng.choice([-1, 1]), 4, 8))
        shared = rng.random() < 0.7
        sd = 0.7
        topo_p = [0.1, 0.1, 0.8]
        delta = float(10 ** rng.uniform(-4, -2.5))
    else:
        b1 = int(rng.choice([b for b in range(3, 9)], p=None))      # right-arm width (rows of b1 tilted squares)
        b2 = int(np.clip(b1 + rng.choice([-1, 0, 0, 1]), 3, 9))    # bottom arm: same width or +-1 (211: 5 and 4)
        shared = rng.random() < 0.6                                 # one angle for both arms (211) or each its own
        sd = 1.5
        topo_p = [0.25, 0.25, 0.5]
        delta = float(10 ** rng.uniform(-4, -1.7))                  # side s = k - delta
    th1 = 2 * math.degrees(math.atan(1 / b1)) + rng.normal(0, sd)
    th2 = th1 if shared else 2 * math.degrees(math.atan(1 / b2)) + rng.normal(0, sd)
    P = dict(
        k=k, n=n, prior=prior,
        delta=delta,
        topology=str(rng.choice(TOPOLOGIES, p=topo_p)),
        arms=[dict(b=b1, theta=th1), dict(b=b2, theta=th2)],
        # bulb: how far each arm stops short of the far wall / the corner block (in units; 308: ~1 and ~b)
        end_gap=[float(rng.uniform(0, 2)) for _ in range(4)],
        corner_gap=float(rng.uniform(-1, 2)),
        # taper per end: rows over which the width ramps from w_tip to b; side the rows lose squares from
        taper_rows=[int(rng.integers(0, 5)) for _ in range(4)],
        taper_tip=[int(rng.integers(1, 4)) for _ in range(4)],
        taper_side=[str(rng.choice(['inner', 'outer', 'both'])) for _ in range(4)],
        wall_gap=[float(abs(rng.normal(0, 0.02))) for _ in range(2)],  # arm against its wall, or slightly off
        row_shear=float(abs(rng.normal(0, 0.05))),                      # per-row slide noise (slid rows / shear)
        # fill (pk/fill.py, 10-20): uprights around the bands, then a top-up placer; the count is a distribution
        # (Evan 10-20: end at n +- 1 or 2; each start is a start for the n it reached)
        base=str(rng.choice(['exact', 'greedy'])),                      # fillfast uprights around the tilted squares
        base_tau=float(rng.choice([0.0, 0.01, 0.03, 0.06])),            # greedy only: tolerated overlap
        extra=int(rng.integers(-2, 3)),                                  # target count n + extra
        stop=str(rng.choice(['count', 'full'])),                         # top up to the target, or until full (<= target + 2)
        topup=str(rng.choice(['bestofK', 'leastbad', 'farthest'], p=[0.45, 0.45, 0.1])),
        topup_K=int(rng.choice([1, 4, 16, 64])),
        drop=str(rng.choice(['random', 'overlap'])),                     # base fill above the target: which uprights go
    )
    return _record(rng, P) if prior == 'record' else P


def _record(rng, P):
    """10-21: overrides fitted to the two-bulb records 273 / 306 / 307 / 308 (PACKER "Overnight 10-20/21", band
    profiles): 9-12 rows, widths ramp 3 -> 6 -> 3 over ~3 rows at each end, corner block ~6.3-7.4 from the corner
    (= arm width ~6.0 + corner_gap 0.3-1.4), far ends ~0.5-1 from the far wall.  Still distributions, wider than the
    four records."""
    P.update(
        corner_gap=float(rng.uniform(0.2, 1.4)),
        end_gap=[float(rng.uniform(0.0, 0.9)) for _ in range(4)],
        rows_mode='round',
        taper_rows=[int(rng.choice([1, 2, 3, 4], p=[0.1, 0.25, 0.45, 0.2])) for _ in range(4)],
        taper_tip=[int(rng.choice([2, 3, 4], p=[0.2, 0.6, 0.2])) for _ in range(4)],
        row_shear=float(abs(rng.normal(0, 0.02))),
    )
    return P


# ---------------------------------------------------------------- arms

def arm(b, theta, length, taper, shear_sd, rng, rows_mode='floor'):
    """A vertical band in local coords: rows of up to b squares at angle theta (deg), stacked along +y with spacing
    sec theta (Ryu Lemma 3.3: in the rotated frame square (i, j) is [i + j tan, i + 1 + j tan] x [j, j + 1]).
    taper = ((rows, tip, side) at the low end, (rows, tip, side) at the high end).  Returns (x, y, deg) list; the
    caller translates it against a wall."""
    t = math.radians(theta)
    c, s = math.cos(t), math.sin(t)
    u, v = np.array([c, s]), np.array([-s, c])
    nrows = max(1, int(length * c) if rows_mode == 'floor' else int(round(length * c)))   # 'round': record prior (10-21)
    out = []
    for j in range(nrows):
        w, lo = b, 0
        for (rows, tip, side), d in ((taper[0], j), (taper[1], nrows - 1 - j)):
            if rows and d < rows:
                wj = min(b, tip + round((b - tip) * d / rows))
                cut = b - wj
                lo_j = 0 if side == 'inner' else cut if side == 'outer' else cut // 2
                if wj < w:
                    w, lo = wj, lo_j
        sh = j * math.tan(t) + rng.normal(0, shear_sd)
        for i in range(lo, lo + w):
            p = (i + 0.5 + sh) * u + (j + 0.5) * v
            out.append((float(p[0]), float(p[1]), float(theta)))
    return out


def _extent(sq):
    C = np.array([pt for q in sq for pt in gen.corners(*q)])
    return C.min(0), C.max(0)


def place_right(sq, s, wall_gap, y0):
    """Translate so the arm's right extent is at s - wall_gap and its bottom extent at y0.  "Inner" = toward -x."""
    lo, hi = _extent(sq)
    dx, dy = s - wall_gap - hi[0], y0 - lo[1]
    return [(float(x + dx), float(y + dy), a) for x, y, a in sq]


def reflect(sq, s):
    """Mirror in the anti-diagonal, (x, y) -> (s - y, s - x): a right-wall arm becomes a bottom-wall arm and the angle
    changes sign (308: -19 / +18 deg).  Arm y-range [y0, y1] becomes x-range [s - y1, s - y0]."""
    return [(s - y, s - x, (-a) % 90) for x, y, a in sq]


def build(P, rng):
    """Tilted squares for one parameter draw, inside the box [0, s]^2.  Arm A on the right wall, arm B on the bottom
    (308's layout; symmetry is free, the quench does not care)."""
    if P['topology'] == 'multi':
        return build_multi(P, rng)
    s = P['k'] - P['delta']
    A, B = P['arms']
    tp = [(P['taper_rows'][i], P['taper_tip'][i], P['taper_side'][i]) for i in range(4)]
    wA = A['b'] * math.cos(math.radians(A['theta'])) + abs(math.sin(math.radians(A['theta'])))   # f(theta)
    wB = B['b'] * math.cos(math.radians(B['theta'])) + abs(math.sin(math.radians(B['theta'])))
    topo = P['topology']
    if topo == 'L':          # A runs the full wall, B stops at A
        a_lo, a_hi = 0.0, s
        b_lo, b_hi = 0.0, s - wA
    elif topo == 'V':        # both run into the corner; B squares that hit A are dropped below
        a_lo, a_hi = 0.0, s
        b_lo, b_hi = 0.0, s
    else:                    # bulb: neither reaches the corner; grid block in [s - wA, s] x [0, wB]
        a_lo, a_hi = wB + P['corner_gap'], s - P['end_gap'][1]
        b_lo, b_hi = P['end_gap'][2], s - wA - P['corner_gap']
    rm = P.get('rows_mode', 'floor')
    TA = arm(A['b'], A['theta'], a_hi - a_lo, (tp[0], tp[1]), P['row_shear'], rng, rm) if a_hi - a_lo > 1.5 else []
    TB = arm(B['b'], B['theta'], b_hi - b_lo, (tp[2], tp[3]), P['row_shear'], rng, rm) if b_hi - b_lo > 1.5 else []
    TA = place_right(TA, s, P['wall_gap'][0], a_lo) if TA else []
    TB = reflect(place_right(TB, s, P['wall_gap'][1], s - b_hi), s) if TB else []   # built on the right, mirrored
    T = []
    for q in TA + TB:        # drop out-of-box and overlapping (V's elbow; any taper collision)
        if gen.wall_viol(q, s) <= 1e-9 and all(gen.pen(q, p) <= 1e-9 for p in T):
            T.append(q)
    return s, T


# ---------------------------------------------------------------- multi-band layouts (10-21)
# Evan 10-21: "Systematically trying different counts / widths that have been undersampled ... trying same-tilt more ...
# Same-tilt also gives a very different set of junction patterns."  A design = one band per chosen wall (0 right,
# 1 bottom, 2 left, 3 top), each with its own width b and signed tilt; every arm is built on the right wall and rotated
# by w quarter turns, (x, y) -> (y, s - x), which keeps the tilt sign (mod 90), so tilt signs are free per arm.
# Corners (0 BR, 1 BL, 2 TL, 3 TR): arm on wall w has its low end at corner w and its high end at corner w - 1.
# A corner shared by two arms draws a junction (bulb: both stop short by the other's width + gap, leaving a corner
# block; L: one runs through, the other stops at it; V: both run in, overlaps dropped); a free end stops end_gap short.

def _rot(sq, s, w):
    for _ in range(w % 4):
        sq = [(y, s - x, a) for x, y, a in sq]
    return sq


def sample_design(rng, k, n, walls, signs, bs, sd=0.7, junction_p=(0.6, 0.2, 0.2)):
    """params for a fixed design (walls, tilt signs, widths); geometry details stay random (record-style tapers/gaps)"""
    P = sample_params(rng, k, n, 'record')
    arms = []
    for w, sg, b in zip(walls, signs, bs):
        th = 2 * math.degrees(math.atan(1 / b)) + rng.normal(0, sd)
        arms.append(dict(wall=int(w), b=int(b), theta=float(sg * th)))
    J = {c: str(rng.choice(['bulb', 'L', 'V'], p=list(junction_p))) for c in range(4)}
    P.update(arms=arms, topology='multi', junction=J, runner={c: int(rng.integers(0, 2)) for c in range(4)},
             corner_gaps=[float(rng.uniform(0.2, 1.4)) for _ in range(4)],
             end_gaps=[float(rng.uniform(0.2, 1.4)) for _ in range(8)],
             tapers=[[(int(rng.choice([1, 2, 3, 4], p=[0.1, 0.25, 0.45, 0.2])), int(rng.choice([2, 3, 4], p=[0.2, 0.6, 0.2])),
                       str(rng.choice(['inner', 'outer', 'both']))) for _ in range(2)] for _ in arms])
    return P


def build_multi(P, rng):
    s = P['k'] - P['delta']
    arms = P['arms']
    W = {a['wall']: a['b'] * math.cos(math.radians(abs(a['theta']))) + abs(math.sin(math.radians(a['theta']))) for a in arms}
    T = []
    for i, a in enumerate(arms):
        w = a['wall']
        res = []
        for end, c in ((0, w % 4), (1, (w - 1) % 4)):
            other = (c + 1) % 4 if end == 0 else c         # the other wall at corner c
            if other in W:
                J = P['junction'][c]
                if J == 'bulb':
                    res.append(W[other] + P['corner_gaps'][c])
                elif J == 'L':
                    # runner index 0 -> the lower-numbered wall of the pair runs through
                    runs = (min(w, other) == w) == (P['runner'][c] == 0)
                    res.append(0.0 if runs else W[other])
                else:
                    res.append(0.0)
            else:
                res.append(P['end_gaps'][2 * i + end])
        lo, hi = res[0], s - res[1]
        if hi - lo <= 1.5:
            continue
        A = arm(a['b'], a['theta'], hi - lo, tuple(P['tapers'][i]), P['row_shear'], rng, 'round')
        A = place_right(A, s, P['wall_gap'][i % 2], lo)
        T += _rot(A, s, w)
    out = []
    for q in T:
        if gen.wall_viol(q, s) <= 1e-9 and all(gen.pen(q, p) <= 1e-9 for p in out):
            out.append(q)
    return s, out


# ---------------------------------------------------------------- fill and count (pk/fill.py)

def overlap_score(sq):
    X = np.array(sq)
    sc = np.zeros(len(sq))
    for i in range(len(sq)):
        near = np.nonzero(np.hypot(*(X[:, :2] - X[i, :2]).T) < 1.42)[0]
        sc[i] = sum(max(0.0, gen.pen(sq[i], sq[j])) for j in near if j != i)
    return sc


def fill_to_count(s, T, P, prng):
    """Uprights around the bands (fillfast exact / greedy), then: above the target, drop uprights; below, top up
    with a pk.fill placer, to the target ('count') or until it reports full, capped at target + 2 ('full').
    The final count is whatever this reaches: no forced insertion.  -> (sq, fill info)"""
    from pk import fill as F
    target = P['n'] + P['extra']
    kw = dict(tau=P['base_tau']) if P['base'] == 'greedy' else {}
    X, i0 = F.fill(s, T, prng, placer=P['base'], m=None, **kw)
    base = len(X)
    ups = [tuple(q) for q in X[len(T):]]
    while len(T) + len(ups) > target and ups:
        if P['drop'] == 'overlap':
            sc = overlap_score(list(T) + ups)[len(T):]
            ups.pop(int(np.argmax(sc)))
        else:
            ups.pop(prng.randrange(len(ups)))
    sq = np.array(list(T) + ups, float).reshape(-1, 3)
    info = dict(base_count=base, topup=None)
    need = target - len(sq)
    if need > 0:
        tp = P['topup']
        kw = dict(K=P['topup_K']) if tp in ('bestofK', 'leastbad') else {}
        if P['stop'] == 'full' and tp != 'farthest':
            kw['m_max'] = need + 2
            X, i1 = F.fill(s, sq, prng, placer=tp, m=None, **kw)
        else:
            X, i1 = F.fill(s, sq, prng, placer=tp, m=need, **kw)
        sq = X
        info['topup'] = {k: i1[k] for k in ('placer', 'params', 'added', 'pen')}
    info.update(count=len(sq), target=target, pen=F.packing_pen(s, sq))
    return [tuple(map(float, q)) for q in sq], info


# ---------------------------------------------------------------- driver

def one(k, n, seed, quench=None, prior='broad', side=None, design=None):
    import random
    rng = np.random.default_rng(seed)
    prng = random.Random(seed)                       # pk.fill draws (random.Random, as pk.moves)
    P = sample_design(rng, k, n, *design) if design else sample_params(rng, k, n, prior)
    if side:                                         # 10-21: build in a larger box (side k + U(lo, hi)), quench compresses
        P['delta'] = -float(rng.uniform(*side))
    s, T = build(P, rng)
    sq, finfo = fill_to_count(s, T, P, prng)
    rec = dict(seed=seed, params=P, s0=s, tilted=len(T), n_target=n, n_out=len(sq), fill=finfo, sq=sq)
    if quench:
        from pk.engine import fq_quench, read_cfg
        from relax import OPS
        with tempfile.TemporaryDirectory() as tmp:
            r = fq_quench((s, sq), OPS[quench], tmp=tmp)
            if r.out:
                s1, sq1 = read_cfg(r.out)
                rec.update(s1=float(s1), sq=[tuple(map(float, q)) for q in sq1])
        from pk.packing import Packing
        q = Packing(rec.get('s1', s), np.array(rec['sq'], float))
        rec.update(pen1=float(q.max_pen()), lines=int(q.grid_lines()))
    return rec


def _one(t):
    return one(*t)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('sample')
    a.add_argument('--k', type=int, required=True)
    a.add_argument('--n', type=int, required=True)
    a.add_argument('--count', type=int, default=20)
    a.add_argument('--seed', type=int, default=0)
    a.add_argument('--out', required=True)
    a.add_argument('--quench', default=None, help='relax.OPS key, e.g. shrink, soft05, loose')
    a.add_argument('--prior', default='broad', choices=PRIORS)
    a.add_argument('--procs', type=int, default=1)
    a.add_argument('--side', default=None, help='lo,hi: build at side k + U(lo, hi) instead of k - delta (then quench)')
    e = sub.add_parser('explore', help='explorer runs from the line-free starts of a sample dir (10-21: composers are '
                       'sources of starts; each explorer run is a walk from them)')
    e.add_argument('--dir', required=True); e.add_argument('--n', type=int, required=True)
    e.add_argument('--runs', type=int, default=6); e.add_argument('--starts-per-run', type=int, default=6)
    e.add_argument('--cpu-h', type=float, default=6.0, help='total, split over the runs')
    e.add_argument('--smax', type=float, default=None); e.add_argument('--seed', type=int, default=0)
    g = ap.parse_args()
    if g.cmd == 'explore':
        return explore(g)
    os.makedirs(g.out, exist_ok=True)
    from concurrent.futures import ProcessPoolExecutor
    from pk import manifest
    man = manifest.start(g.out, 'bulbgen', seed=g.seed, args=vars(g), procs=g.procs)
    with open(os.path.join(g.out, 'starts.jsonl'), 'a') as f, ProcessPoolExecutor(g.procs) as ex:
        for r in ex.map(_one, [(g.k, g.n, g.seed + i, g.quench, g.prior,
                                     tuple(map(float, g.side.split(','))) if getattr(g, 'side', None) else None)
                                    for i in range(g.count)], chunksize=1):
            gen.write(os.path.join(g.out, f's{r["seed"]:05d}_n{r["n_out"]}.txt'), r.get('s1', r['s0']), r['sq'])
            sq = r.pop('sq')
            f.write(json.dumps(r) + '\n')
            P = r['params']
            print(f"seed {r['seed']}: {P['topology']:4s} b={P['arms'][0]['b']},{P['arms'][1]['b']} "
                  f"th={P['arms'][0]['theta']:.1f},{P['arms'][1]['theta']:.1f} tilted={r['tilted']} "
                  f"base={r['fill']['base_count']} n_out={r['n_out']} (target {r['fill']['target']}) pen={r['fill']['pen']:.3f} s0={r['s0']:.5f}"
                  + (f" s1={r['s1']:.6f} pen1={r['pen1']:.1e} lines={r['lines']}" if 's1' in r else ''), flush=True)
    man.finish()


def explore(g):
    """Starts = line-free, feasible, count n outputs of `sample` (any side), spread over the runs round-robin in side
    order so every run gets a range; explorer as SW (--adapt --adapt-cells --sparse-share 0.5), CPU-budgeted."""
    import subprocess
    from pk import manifest
    R = [json.loads(l) for l in open(os.path.join(g.dir, 'starts.jsonl'))]
    ok = sorted((r for r in R if r['n_out'] == g.n and r.get('lines') == 0 and r.get('pen1', 1) <= 1e-9 and 's1' in r),
                key=lambda r: r['s1'])
    seen, S = set(), []
    for r in ok:
        if round(r['s1'], 7) not in seen:
            seen.add(round(r['s1'], 7)); S.append(r)
    S = S[:g.runs * g.starts_per_run]
    print(f'{len(ok)} line-free n={g.n} starts ({len(S)} used; sides {S[0]["s1"]:.5f}..{S[-1]["s1"]:.5f})' if S else 'no starts',
          flush=True)
    if not S:
        return
    smax = g.smax or max(r['s1'] for r in S) + 0.01
    man = manifest.start(os.path.join(g.dir, 'explore'), 'bulbgen-explore', seed=g.seed, args=vars(g), procs=g.runs)
    P = []
    for i in range(g.runs):
        mine = S[i::g.runs]
        if not mine:
            continue
        out = os.path.join(g.dir, 'explore', f'x{i}'); os.makedirs(out, exist_ok=True)
        paths = [os.path.join(g.dir, f's{r["seed"]:05d}_n{g.n}.txt') for r in mine]
        cmd = [sys.executable, os.path.join(HERE, 'explore.py'), '--n', str(g.n), '--starts', *paths,
               '--cpu-h', str(g.cpu_h / g.runs), '--procs', '1', '--seed', str(g.seed + 900 + i), '--adapt',
               '--adapt-cells', '--sparse-share', '0.5', '--smax', str(smax), '--frontier-s', str(float(math.ceil(math.sqrt(g.n) - 1e-12))),
               '--final-polish', '8', '--out', out]
        P.append(subprocess.Popen(cmd, stdout=open(os.path.join(out, 'log.txt'), 'w'), stderr=subprocess.STDOUT))
    for p in P:
        p.wait()
    man.finish()


if __name__ == '__main__':
    main()
