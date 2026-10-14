#!/usr/bin/env python3
"""Rectangle library (PLAN-2026-10-12 item 2; 10-12).

Pieces = corner crops of known packings (the Mondrian atlas, `mondrian.py atlas`), re-optimised as packings of m unit
squares in an (h + d) x h rectangle (d integer) with fq's rectangle mode (`--rect-dx d`, h minimised).  This is the step
Ostashev's library adds over crops (#528 / #562): a piece relaxes in its own rectangle, so it is no longer pinned to
the source's side.  Output = an atlas in the Mondrian schema (W = h + d, H = h, m, sq = the re-optimised piece as its
own source), so `mondrian.py compose --atlas OUT` (exact guillotine upright counts: every win is a construction) runs
on it unchanged.

  python3 rectlib.py build --atlas runs/mondrian/atlas_v2.npz --out runs/rectlib/lib1 [--mmin 8 --mmax 250]
                           [--per-cell 1] [--procs 12] [--limit N]
  python3 rectlib.py table runs/rectlib/lib1          (best h per (m, d) vs the crop's H; gain histogram)
  python3 rectlib.py atlas runs/rectlib/lib1 --out runs/rectlib/lib1/atlas.npz [--with-crops ATLAS]
"""
from __future__ import annotations
import argparse, json, math, os, subprocess, sys, tempfile, time
from concurrent.futures import ProcessPoolExecutor

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FQ = os.path.join(HERE, 'target', 'release', 'fq')
FT = 1e-9


def crop_squares(sq, s, o, W, H):
    """Squares of source (sq, s) in orientation o (mondrian.orient) lying wholly in [0, W] x [0, H]."""
    import mondrian
    q = mondrian.orient(np.asarray(sq, float), s, o & 3)
    if o & 4:
        q = np.c_[q[:, 1], q[:, 0], (-q[:, 2]) % 90.0]
    h = mondrian.half_ext(q[:, 2])
    ok = (q[:, 0] - h >= -1e-7) & (q[:, 1] - h >= -1e-7) & (q[:, 0] + h <= W + 1e-7) & (q[:, 1] + h <= H + 1e-7)
    return q[ok]


def _fq(inp, outp, d, extra):
    """fq rectangle quench ([0, s + d] x [0, s]); extra = fq flags on top (e.g. ['--no-alm'])."""
    from pk.engine import QOpt, fq_quench
    r = fq_quench(inp, QOpt(loosen=1.0, rect_dx=float(d), flip_top=32).with_args(extra), out=outp, binary=FQ, timeout=900)
    try:
        js = json.loads(r.stdout.strip().splitlines()[-1])
        L = open(outp).read().split('\n')
        sq = np.array([[float(t) for t in l.split()[:3]] for l in L[1:] if l.strip()])
    except Exception:
        return None, (r.stderr or r.stdout)[-300:]
    if not (js['min_gap'] > -1e-9 and js['min_wall'] > -1e-9):
        return None, 'infeasible output'
    return (float(js['s']), js['sec'], sq), None


def reopt(args):
    """fq rectangle quench of one crop (full ALM quench + polish; if that ends above the crop's own start side, polish
    only from the feasible crop).  Returns dict(h, d, m, sq, how) or an error."""
    key, P, W, H, d = args
    m = len(P)
    s0 = max(H, W - d)
    with tempfile.TemporaryDirectory() as td:
        inp, outp = os.path.join(td, 'in.txt'), os.path.join(td, 'out.txt')
        with open(inp, 'w') as f:
            f.write(f'{m} {float(s0)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(a)!r}\n' for x, y, a in P))
        best, err = _fq(inp, outp, d, [])
        how, sec = 'alm', best[1] if best else 0.0
        if best is None or best[0] > s0 + 1e-9:
            pol, err2 = _fq(inp, outp, d, ['--no-alm'])
            if pol is not None and (best is None or pol[0] < best[0]):
                best, how = pol, 'polish'
            sec += pol[1] if pol else 0.0
            err = err or err2
        if best is None:
            return dict(key=key, err=err)
    h, _, sq = best
    return dict(key=key, m=m, d=d, W0=float(W), H0=float(H), s0=float(s0), h=h, ok=True, how=how, sec=sec,
                sq=sq.tolist())


def select(A, mmin, mmax, per_cell, hmax):
    """Atlas pieces to re-optimise: m in range, H <= hmax, best per_cell per (m, floor W, floor H) by area."""
    W, H, m = A['W'], A['H'], A['m']
    idx = np.nonzero((m >= mmin) & (m <= mmax) & (H <= hmax))[0]
    cells = {}
    for i in idx:
        k = (int(m[i]), int(math.floor(W[i] + FT)), int(math.floor(H[i] + FT)))
        cells.setdefault(k, []).append(i)
    out = []
    for k, L in cells.items():
        L.sort(key=lambda i: W[i] * H[i])
        out += L[:per_cell]
    return sorted(out)


def cmd_build(a):
    import mondrian
    from pk import manifest
    import scratch
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'rectlib', args=vars(a), procs=a.procs)
    scratch.scope('rectlib')
    A = mondrian.load_atlas(a.atlas)
    sel = select(A, a.mmin, a.mmax, a.per_cell, a.hmax)
    if a.limit:
        rng = np.random.default_rng(a.seed)
        sel = sorted(rng.choice(sel, size=min(a.limit, len(sel)), replace=False))
    tasks = []
    for i in sel:
        sid = int(A['sid'][i]); sq = A['sq'][sid]; s = float(A['src_s'][sid])
        W, H = float(A['W'][i]), float(A['H'][i])
        P = crop_squares(sq, s, int(A['o'][i]), W, H)
        if len(P) != int(A['m'][i]):
            continue                                   # orientation convention mismatch: skip, counted below
        for d in sorted({int(math.floor(W - H + FT)), int(math.ceil(W - H - FT))} & set(range(a.dmax + 1))):
            tasks.append((dict(piece=int(i), src=str(A['labels'][sid]), src_n=int(A['src_n'][sid])), P, W, H, d))
    print(f'{len(sel)} pieces selected, {len(tasks)} re-optimisations', flush=True)
    res_path = os.path.join(a.out, 'pieces.jsonl')
    t0 = time.time(); nd = 0
    with ProcessPoolExecutor(a.procs) as ex, open(res_path, 'w') as f:
        for r in ex.map(reopt, tasks, chunksize=4):
            f.write(json.dumps(r) + '\n'); nd += 1
            if nd % 500 == 0:
                print(f'  {nd}/{len(tasks)} ({time.time() - t0:.0f} s)', flush=True)
    man.finish(summary=dict(selected=len(sel), tasks=len(tasks)))


def load(dirp):
    return [json.loads(l) for l in open(os.path.join(dirp, 'pieces.jsonl'))]


def cmd_table(a):
    R = [r for r in load(a.dir) if 'h' in r and r['ok']]
    gain = np.array([r['s0'] - r['h'] for r in R])
    print(f'{len(R)} re-optimised pieces; h below the crop start s0 by: median {np.median(gain):.2e}, '
          f'p90 {np.percentile(gain, 90):.2e}, max {gain.max():.3f}; > 1e-3: {(gain > 1e-3).sum()}, > 1e-2: {(gain > 1e-2).sum()}')
    best = {}
    for r in R:
        k = (r['m'], r['d'])
        if k not in best or r['h'] < best[k]['h']:
            best[k] = r
    print(f'{len(best)} (m, d) cells')
    if a.show:
        for k in sorted(best)[:a.show]:
            r = best[k]
            print(f"  m={k[0]:4d} d={k[1]}  h {r['h']:.9f}  (crop {r['H0']:.6f}x{r['W0']:.6f} from {r['key']['src']})")


def cmd_atlas(a):
    """Library -> Mondrian atlas schema; each re-optimised piece is its own source (orientation 0, at the origin)."""
    R = [r for r in load(a.dir) if 'h' in r and r['ok']]
    W = np.array([r['h'] + r['d'] for r in R]); H = np.array([r['h'] for r in R]); m = np.array([r['m'] for r in R])
    keep = m > np.floor(W + FT) * np.floor(H + FT)               # beats the grid in its own box
    R = [r for r, k in zip(R, keep) if k]
    W, H, m = W[keep], H[keep], m[keep]
    import mondrian
    W, H, m, idx = mondrian.pareto_prune(W, H, m, np.arange(len(R)))
    R = [R[i] for i in idx]
    # mondrian convention: W >= H, source box = W x W with the piece at its lower-left corner (o = 0)
    sq = mondrian._objarr([np.array(r['sq']) for r in R])
    Z = dict(W=W, H=H, m=m, o=np.zeros(len(R), int), sid=np.arange(len(R)),
             labels=np.array([f"rect:{r['key']['src']}:{r['key']['piece']}:d{r['d']}" for r in R]),
             src_n=np.array([r['key']['src_n'] for r in R]), src_s=W.copy(), sq=sq)
    print(f'{len(R)} library pieces')
    if a.with_crops:                              # merged atlas: crop pieces + library pieces (sources concatenated)
        C = mondrian.load_atlas(a.with_crops)
        off = len(C['sq'])
        Z = dict(W=np.r_[C['W'], Z['W']], H=np.r_[C['H'], Z['H']], m=np.r_[C['m'], Z['m']], o=np.r_[C['o'], Z['o']],
                 sid=np.r_[C['sid'], Z['sid'] + off], labels=np.r_[C['labels'].astype(str), Z['labels']],
                 src_n=np.r_[C['src_n'], Z['src_n']], src_s=np.r_[C['src_s'], Z['src_s']],
                 sq=mondrian._objarr(list(C['sq']) + list(Z['sq'])))
        print(f'merged with {len(C["W"])} crop pieces')
    np.savez_compressed(a.out, **Z)


def _column(piece, c, s, top):
    """Squares of one composer column [0, c] x [0, s]: the piece (lower-left frame, w x h) at the bottom-left corner
    (top-left mirrored to the top if top) plus the better of the two upright variants (as mondrian.Composer.column)."""
    fl = lambda v: math.floor(v + FT)
    out = []
    if piece is None:
        for i in range(fl(c)):
            for j in range(fl(s)):
                out.append((i + 0.5, j + 0.5, 0.0))
        return out
    P, w, h = piece
    v1 = fl(s - h) * fl(c) + fl(c - w) * fl(h)
    v2 = fl(c - w) * fl(s) + fl(s - h) * fl(w)
    out += [tuple(q) for q in P]
    if v1 >= v2:
        out += [(i + 0.5, h + j + 0.5, 0.0) for i in range(fl(c)) for j in range(fl(s - h))]
        out += [(w + i + 0.5, j + 0.5, 0.0) for i in range(fl(c - w)) for j in range(fl(h))]
    else:
        out += [(w + i + 0.5, j + 0.5, 0.0) for i in range(fl(c - w)) for j in range(fl(s))]
        out += [(i + 0.5, h + j + 0.5, 0.0) for i in range(fl(w)) for j in range(fl(s - h))]
    if top:
        out = [(x, s - y, (-a) % 90.0) for x, y, a in out]
    return out


def build_layout(A, r):
    """Coordinates of one `mondrian.py compose` result r (dict with s, c, left, right) on atlas A."""
    s, c = r['s'], r['c']

    def piece(info):
        if not info:
            return None
        i = info['piece']; sid = int(A['sid'][i])
        P = crop_squares(A['sq'][sid], float(A['src_s'][sid]), int(A['o'][i]), float(A['W'][i]), float(A['H'][i]))
        if info['tr']:
            P = np.c_[P[:, 1], P[:, 0], (-P[:, 2]) % 90.0]
        return P, info['w'], info['h']
    L = _column(piece(r['left']), c, s, top=False)
    R = [(c + x, y, a) for x, y, a in _column(piece(r['right']), s - c, s, top=True)]
    return np.array(L + R, float)


def cmd_construct(a2):
    """compose JSON -> packings (trimmed to n, side = s), fq polish, f64 check; writes OUT/n-<n>.txt + summary."""
    import mondrian
    A = mondrian.load_atlas(a2.atlas)
    rep = json.load(open(a2.compose))
    os.makedirs(a2.out, exist_ok=True)
    summ = []
    for r in rep:
        if r['gap'] > a2.margin:
            continue
        X = build_layout(A, r)
        n = r['n']
        if len(X) < n:
            summ.append(dict(n=n, err=f'built {len(X)} < n')); continue
        X = mondrian.trim(X, n)
        raw = os.path.join(a2.out, f'raw-{n}.txt'); outp = os.path.join(a2.out, f'n-{n}.txt')
        with open(raw, 'w') as f:
            f.write(f'{n} {float(r["s"])!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in X))
        chk = subprocess.run([FQ, 'check', '--in', raw], capture_output=True, text=True).stdout.split()
        from pk.engine import fq_quench, POLISH_DEEP
        q = fq_quench(raw, POLISH_DEEP, out=outp, binary=FQ)
        try:
            js = json.loads(q.stdout.strip().splitlines()[-1])
        except Exception:
            summ.append(dict(n=n, err=q.stderr[-200:])); continue
        e = dict(n=n, s_comp=r['s'], s_pol=js['s'], raw_gap=float(chk[2]), raw_wall=float(chk[4]),
                 front_store=r['front'], left=(r['left'] or {}).get('src'), right=(r['right'] or {}).get('src'))
        summ.append(e)
        print(f"{n:4d} composed {r['s']:.10f} (raw gap {e['raw_gap']:.1e} wall {e['raw_wall']:.1e}) polished "
              f"{js['s']:.10f}  store front {r['front']:.10f}", flush=True)
    json.dump(summ, open(os.path.join(a2.out, 'summary.json'), 'w'), indent=1)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    b = sp.add_parser('build')
    b.add_argument('--atlas', default=os.path.join(HERE, 'runs', 'mondrian', 'atlas_v2.npz'))
    b.add_argument('--out', required=True); b.add_argument('--mmin', type=int, default=8)
    b.add_argument('--mmax', type=int, default=250); b.add_argument('--hmax', type=float, default=18.5)
    b.add_argument('--per-cell', type=int, default=1); b.add_argument('--procs', type=int, default=12)
    b.add_argument('--limit', type=int, default=0); b.add_argument('--dmax', type=int, default=8); b.add_argument('--seed', type=int, default=1)
    t = sp.add_parser('table'); t.add_argument('dir'); t.add_argument('--show', type=int, default=0)
    t2 = sp.add_parser('atlas'); t2.add_argument('dir'); t2.add_argument('--out', required=True)
    t2.add_argument('--with-crops')
    t3 = sp.add_parser('construct'); t3.add_argument('--atlas', required=True); t3.add_argument('--compose', required=True)
    t3.add_argument('--out', required=True); t3.add_argument('--margin', type=float, default=1e-3)
    a = ap.parse_args()
    {'build': cmd_build, 'table': cmd_table, 'atlas': cmd_atlas, 'construct': cmd_construct}[a.cmd](a)
