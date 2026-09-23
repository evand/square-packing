#!/usr/bin/env python3
"""Near-tight pose set of a certificate scaled slightly beyond its critical container.

Let a certificate have weights w with total W (numerators over the file's weight denominator
WD), and at container side t let m(t) be the EXACT minimum captured weight of a closed unit
square over all placements (the verifier's `min`, in the same numerator units).  Above the
critical side m < WD, and w/m is a cover of total W/m = 12 + E(t), E(t) = W/m - 12.

Any packing of 12 unit squares into a container of side t' < t, scaled up to side t, gives 12
closed unit squares with pairwise disjoint closures; square i captures c_i * m with c_i >= 1 and
sum_i c_i m <= W (no point is counted twice), so  sum_i (c_i - 1) <= W/m - 12 = E.  Every term is
>= 0, hence EVERY square satisfies  captured weight <= m (1 + E) = m + (W - 12 m) = W - 11 m.  (*)
So each of the 12 squares lies in the near-tight set T_E = {poses Q : weight(Q) <= W - 11 m}.
In numerator units the threshold W - 11 m is an exact integer -- no rounding is needed.

The verifier examines, per angle bin [theta_k, theta_{k+1}], the sigma_k-shrunk square at
theta_k; a true unit square at any angle in the bin contains that shrunk square and so captures
at least what it captures.  Hence the set of verifier cells with sum <= threshold is a SUPERSET
of the true T_E: the safe direction for a case analysis (nothing a packing could use is left out).

For each t in a list this script
  1. rescales the certificate to container t (same integer coordinates times 1000, new
     denominator -- the arithmetic of search/nu_f.py::rescale_cert),
  2. runs the verifier to get m(t) exactly,
  3. runs the verifier again with TIGHT_DUMP / TIGHT_THRESH = W - 11 m (see verify/src/main.rs):
     per angle bin it records the number of arrangement cells with captured weight <= threshold,
     their area (clipped to the rotated admissible centre box), a G x G area histogram of their
     midpoints over the container, and a uniform reservoir sample of the cells themselves
     (TIGHT_MAX in total) with exact integer coordinates and sums,
  4. analyses the dump: pose-space measure of T_E = sum_k area_k x (theta_{k+1} - theta_k) as
     a fraction of the admissible pose space sum_k (s - w_k)^2 x width_k over [0, 45 deg],
     angle histogram, spatial distribution (folded to a quadrant by the C4 rotations, the
     pose-wise symmetry of the [0, 45 deg] sweep), connected components of the set of
     (x, y, theta) voxels of the histogram (26-connectivity), and, on the sampled cells, an
     independent re-clipping of the areas as a consistency check,
  5. plots the sampled cell midpoints coloured by angle and the measure map into
     runs/tight_<t>.png,
and writes a table to search/TIGHTSET.md.

The pseudo-value  t = crit  uses the shipped certificate itself with threshold = m (the exactly
tight set at the critical size; there E < 0 and (*) is vacuous -- this is the reference set).

Usage:
    python3 search/tightset.py [--ts crit,3.9687,3.9690,3.9695,3.970,3.972,3.975,3.980,3.990]
                               [--N 6000] [--threads 32] [--cert certificates/s12_lower_3.9686.txt]
                               [--outdir runs] [--grid 48] [--max-cells 300000] [--dtheta 0.5]
                               [--no-plot] [--md search/TIGHTSET.md] [--min-only] [--stop-after-dump]
Intermediate files (runs/tight_<t>_cert.txt, _min.txt, _cells.txt, _stats.json) are reused
if present; delete them to recompute.
"""
import argparse, json, math, os, subprocess, sys, time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VERIFY = os.path.join(ROOT, 'verify', 'target', 'release', 'verify')


# ----------------------------------------------------------------------------- rescale
def read_cert(path):
    t = open(path).read().split(); it = iter(t)
    sn, sd, D, WD, m = [int(next(it)) for _ in range(5)]
    rows = np.array([[int(next(it)), int(next(it)), int(next(it))] for _ in range(m)], dtype=np.int64)
    return sn, sd, D, WD, rows


def rescale_cert(path, s_target, out):
    """Identical to search/nu_f.py::rescale_cert (kept local so this script needs no scipy.optimize):
    integer coordinates times 1000, denominator D' = round(K*1000/s_target) where K = s*D."""
    sn, sd, D, WD, rows = read_cert(path)
    K = Fraction(sn, sd) * D; assert K.denominator == 1; K = int(K)
    MUL = 1000; Dp = int(round(K * MUL / s_target)); sp_ = Fraction(K * MUL, Dp)
    with open(out, 'w') as f:
        f.write(f"{sp_.numerator} {sp_.denominator}\n{Dp}\n{WD}\n{len(rows)}\n")
        for a, b, c in rows: f.write(f"{a * MUL} {b * MUL} {c}\n")
    return rows[:, 2].sum() / WD, float(sp_)


# ----------------------------------------------------------------------------- verifier
def run_verifier(path, N, threads, n=12, dump=None, thresh=None, max_cells=None, grid=None):
    env = dict(os.environ)
    if dump is not None:
        env['TIGHT_DUMP'] = dump; env['TIGHT_THRESH'] = str(thresh)
        if max_cells is not None: env['TIGHT_MAX'] = str(max_cells)
        if grid is not None: env['TIGHT_GRID'] = str(grid)
    r = subprocess.run([VERIFY, path, str(n), str(N), str(threads), '0'], capture_output=True, text=True, env=env)
    if r.returncode != 0: raise RuntimeError(r.stdout[-800:] + r.stderr[-800:])
    out = r.stdout
    tot = None; mn = None; ncell = None
    for l in out.split('\n'):
        if l.startswith('atoms='):
            fr = l.split('total weight =')[1].split('=')[0].strip(); a, b = fr.split('/'); tot = (int(a), int(b))
        if l.startswith('min covered'):
            fr = l.split('=')[1].strip().split()[0]; a, b = fr.split('/'); mn = (int(a), int(b))
        if l.startswith('TIGHT_DUMP:'): ncell = int(l.split()[1])
    if tot is None or mn is None: raise RuntimeError(out[-800:])
    return tot, mn, ncell, out


# ----------------------------------------------------------------------------- dump parsing
def parse_dump(path):
    """Returns (N, G, bins, hist, cells).
    bins: dict k -> dict(deg, DEN, cn, sn, g, L, ntight, area, area_box, nwritten)
    hist: dict k -> list of (i, j, area)
    cells: arrays with float positions (container units) of the sampled cells."""
    bins = {}; hist = {}
    ks = []; a_ = []; b_ = []; c0_ = []; c1_ = []; sm = []; cx = []; cy = []; ar = []
    N = None; G = None
    with open(path) as f:
        for line in f:
            ch = line[0]
            if ch == 'c':
                p = line.split()
                k = int(p[1]); DEN = bins[k]['DEN']
                ks.append(k); a_.append(int(p[2]) / DEN); b_.append(int(p[3]) / DEN)
                c0_.append(int(p[4]) / DEN); c1_.append(int(p[5]) / DEN)
                sm.append(int(p[6])); cx.append(float(p[7])); cy.append(float(p[8])); ar.append(float(p[9]))
            elif ch == 'h':
                p = line.split(); hist[int(p[1])].append((int(p[2]), int(p[3]), float(p[4])))
            elif ch == 'b':
                p = line.split()
                k = int(p[1]); N = int(p[2]); DEN = int(p[4])
                bins[k] = dict(deg=float(p[3]), DEN=DEN, cn=int(p[5]), sn=int(p[6]), g=int(p[7]),
                               L=int(p[11]) / (2 * int(p[12])), ntight=int(p[13]), area=float(p[14]),
                               area_box=float(p[15]), nwritten=int(p[16]))
                hist[k] = []
            elif ch == '#':
                if 'grid=' in line: G = int(line.split('grid=')[1].split()[0])
    cells = dict(k=np.array(ks, dtype=int), a=np.array(a_), b=np.array(b_), c0=np.array(c0_), c1=np.array(c1_),
                 sum=np.array(sm, dtype=float), cx=np.array(cx), cy=np.array(cy), area=np.array(ar))
    return N, G, bins, hist, cells


def clip_area(rect, poly):
    """Sutherland-Hodgman: area of the rectangle (x0,x1,y0,y1) clipped to the convex polygon
    poly (list of (x,y), counter-clockwise)."""
    x0, x1, y0, y1 = rect
    out = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    n = len(poly)
    for i in range(n):
        px, py = poly[i]; qx, qy = poly[(i + 1) % n]
        ex, ey = qx - px, qy - py
        inp = out; out = []
        if not inp: break
        def side(pt): return ex * (pt[1] - py) - ey * (pt[0] - px)   # >= 0 : inside (left of edge)
        s_prev = inp[-1]; d_prev = side(s_prev)
        for cur in inp:
            d_cur = side(cur)
            if d_cur >= 0:
                if d_prev < 0:
                    tt = d_prev / (d_prev - d_cur); out.append((s_prev[0] + tt * (cur[0] - s_prev[0]), s_prev[1] + tt * (cur[1] - s_prev[1])))
                out.append(cur)
            elif d_prev >= 0:
                tt = d_prev / (d_prev - d_cur); out.append((s_prev[0] + tt * (cur[0] - s_prev[0]), s_prev[1] + tt * (cur[1] - s_prev[1])))
            s_prev, d_prev = cur, d_cur
    if len(out) < 3: return 0.0
    A = 0.0
    for i in range(len(out)):
        x, y = out[i]; xn, yn = out[(i + 1) % len(out)]; A += x * yn - xn * y
    return abs(A) / 2


def runs_of(sorted_ints):
    out = []
    if not sorted_ints: return out
    st = pv = sorted_ints[0]
    for k in sorted_ints[1:]:
        if k != pv + 1: out.append((st, pv)); st = k
        pv = k
    out.append((st, pv)); return out


def analyse(N, G, bins, hist, cells, s, dtheta_deg=0.5):
    from scipy import ndimage
    kmax = max(bins)
    width = {k: 2 * (math.atan((k + 1) / N) - math.atan(k / N)) for k in bins}
    meas_k = {k: bins[k]['area'] * width[k] for k in bins}
    measure = sum(meas_k.values())
    total_pose = sum(bins[k]['area_box'] * width[k] for k in bins)
    ntight = sum(b['ntight'] for b in bins.values()); nwritten = sum(b['nwritten'] for b in bins.values())
    present = sorted(k for k in bins if bins[k]['ntight'] > 0)
    def deg_hi(k): return 2 * math.degrees(math.atan((k + 1) / N))
    intervals = [(bins[a]['deg'], deg_hi(b), sum(meas_k[k] for k in range(a, b + 1)) / measure if measure else 0.0)
                 for a, b in runs_of(present)]
    edges = np.arange(0, 46, 1.0)
    hist_deg = np.zeros(len(edges) - 1)
    for k in bins:
        i = min(int(bins[k]['deg']), len(hist_deg) - 1); hist_deg[i] += meas_k[k]
    # measure map over the container: H[i, j] (box size s/G), plus a 3-d voxel set for clusters
    H = np.zeros((G, G)); nth = int(math.ceil(45.0 / dtheta_deg)) + 1
    V = np.zeros((G, G, nth), dtype=bool)
    for k, lst in hist.items():
        if not lst: continue
        w = width[k]; ti = int(bins[k]['deg'] / dtheta_deg)
        for i, j, a in lst:
            H[i, j] += a * w; V[i, j, ti] = True
    _, ncl = ndimage.label(V, structure=np.ones((3, 3, 3), dtype=bool))
    _, ncl_xy = ndimage.label(V.any(axis=2), structure=np.ones((3, 3), dtype=bool))
    lab, _ = ndimage.label(V, structure=np.ones((3, 3, 3), dtype=bool))
    # measure carried by each component
    comp = {}
    for k, lst in hist.items():
        ti = int(bins[k]['deg'] / dtheta_deg); w = width[k]
        for i, j, a in lst:
            c = int(lab[i, j, ti]); comp[c] = comp.get(c, 0.0) + a * w
    comp_sizes = sorted(comp.values(), reverse=True)
    # fold by C4: (x, y) -> (s - y, x) maps box (i, j) -> (G-1-j, i)
    Hf = H + np.rot90(H, 1) + np.rot90(H, 2) + np.rot90(H, 3)
    h = G // 2
    Q = Hf[:h, :h] / max(measure, 1e-300)                  # quadrant [0, s/2]^2, fractions of the measure
    # radial-ish descriptors: distance of the box centre from the nearest wall, in bands of s/G
    box = s / G
    ii = (np.arange(G) + 0.5) * box
    dwall = np.minimum(np.minimum(ii[:, None], ii[None, :]), np.minimum(s - ii[:, None], s - ii[None, :]))
    bands = np.arange(0.5, s / 2 + box, box)
    wall_hist, _ = np.histogram(dwall.ravel(), bins=bands, weights=H.ravel() / max(measure, 1e-300))
    cxs = (np.arange(G) + 0.5) * box
    # marginal of the folded centre coordinate (both coordinates pooled), in boxes
    marg = (Hf.sum(axis=1)[:h] + Hf.sum(axis=0)[:h]) / (2 * max(measure, 1e-300))
    # corner / wall / centre shares
    near_corner = float(Q[cxs[:h] < 0.75][:, cxs[:h] < 0.75].sum())
    near_wall = float(Q[cxs[:h] < 0.75].sum() + Q[:, cxs[:h] < 0.75].sum() - near_corner)
    near_centre = float(Q[cxs[:h] > s / 2 - 0.75][:, cxs[:h] > s / 2 - 0.75].sum())
    # consistency check: the Rust f64 polygon clip vs an independent Python clip, on the sampled cells
    ratio = None
    if len(cells['k']):
        polys = {}
        for k in set(cells['k'].tolist()):
            bn = bins[k]; L = bn['L']; U = s - L; cn, sn = bn['cn'] / bn['g'], bn['sn'] / bn['g']
            polys[k] = [(cn * x + sn * y, -sn * x + cn * y) for x, y in [(L, L), (U, L), (U, U), (L, U)]]
        ex = 0.0; raw = 0.0
        for i in range(len(cells['k'])):
            k = cells['k'][i]; bn = bins[k]; L = bn['L']; U = s - L
            r = (cells['a'][i], cells['b'][i], cells['c0'][i], cells['c1'][i])
            ex += clip_area(r, polys[k]); raw += cells['area'][i]
        ratio = ex / raw if raw > 0 else None
    return dict(ncells=ntight, nwritten=nwritten, measure=float(measure), total_pose=float(total_pose),
                fraction=float(measure / total_pose), nbins_present=len(present), nbins=len(bins),
                angle_intervals=intervals, hist_deg=hist_deg.tolist(), hist_edges=edges.tolist(),
                nclusters=int(ncl), nclusters_xy=int(ncl_xy), cluster_measure=comp_sizes[:24],
                n_voxels=int(V.sum()), n_voxels_total=int(sum(b['area_box'] > 0 for b in bins.values())),
                grid=G, dtheta_deg=dtheta_deg, H=H.tolist(), Q=Q.tolist(), wall_bands=bands.tolist(), wall_hist=wall_hist.tolist(),
                marg=marg.tolist(), near_corner=near_corner, near_wall=near_wall, near_centre=near_centre,
                py_over_rust_area=ratio, mean_deg=float(sum(meas_k[k] * bins[k]['deg'] for k in bins) / measure) if measure else None)


def plot(tstr, s, cells, bins, st, out, thresh_str):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(15, 7.4))
    ax = axs[0]
    degs = np.array([bins[k]['deg'] for k in cells['k']])
    area = cells['area']
    def frame(ax):
        ax.add_patch(plt.Rectangle((0, 0), s, s, fill=False, lw=1.2, color='k'))
        ax.add_patch(plt.Rectangle((0.5, 0.5), s - 1, s - 1, fill=False, lw=0.6, color='gray', ls='--'))
        ax.add_patch(plt.Rectangle((0, 0), s / 2, s / 2, fill=False, lw=0.9, color='tab:red', ls='-.'))
        ax.add_patch(plt.Polygon([(0, 0), (s / 2, 0), (s / 2, s / 2)], closed=True, fill=True, alpha=0.08, color='tab:red', lw=0))
        ax.plot([0, s / 2], [0, s / 2], color='tab:red', lw=0.9, ls='-.')
        ax.set_xlim(-0.05, s + 0.05); ax.set_ylim(-0.05, s + 0.05); ax.set_aspect('equal'); ax.set_xlabel('cx'); ax.set_ylabel('cy')
    frame(ax)
    if len(degs):
        order = np.argsort(area)
        sz = 1.5 + 30 * np.sqrt(area / area.max()) if area.max() > 0 else 3
        sc = ax.scatter(cells['cx'][order], cells['cy'][order], c=degs[order], cmap='viridis', vmin=0, vmax=45,
                        s=np.asarray(sz)[order] if np.ndim(sz) else sz, alpha=0.55, linewidths=0)
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03); cb.set_label('bin angle theta_k (deg); sweep covers [0, 45]')
    ax.set_title(f"t = {tstr}: {len(degs)} sampled cells of {st['ncells']} with weight <= {thresh_str}\n"
                 f"E = {st['E']:.4f}; T_E = {100 * st['fraction']:.2f}% of pose space; {st['nclusters']} clusters "
                 f"(grid {s / st['grid']:.3f}, {st['dtheta_deg']} deg)", fontsize=9)
    ax = axs[1]; frame(ax)
    H = np.array(st['H']); G = st['grid']
    im = ax.imshow(H.T, origin='lower', extent=(0, s, 0, s), cmap='magma', interpolation='nearest')
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.set_label('measure of T_E per box (all angles)')
    ax.set_title("measure map (cap-independent); red: C4 quadrant (pose symmetry of the sweep)\nand D4 fundamental triangle of the container", fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def fmt_intervals(iv, maxn=5):
    iv2 = sorted(iv, key=lambda x: -x[2])
    parts = [f"{x:.2f}-{y:.2f} ({100 * z:.0f}%)" for x, y, z in iv2[:maxn]]
    return ', '.join(parts) + (f", +{len(iv) - maxn} more" if len(iv) > maxn else '')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ts', default='crit,3.9687,3.9690,3.9695,3.970,3.972,3.975,3.980,3.990')
    ap.add_argument('--N', type=int, default=6000)
    ap.add_argument('--threads', type=int, default=os.cpu_count() or 16)
    ap.add_argument('--cert', default=os.path.join(ROOT, 'certificates', 's12_lower_3.9686.txt'))
    ap.add_argument('--outdir', default=os.path.join(ROOT, 'runs'))
    ap.add_argument('--md', default=os.path.join(HERE, 'TIGHTSET.md'))
    ap.add_argument('--grid', type=int, default=48, help='TIGHT_GRID: G x G area histogram')
    ap.add_argument('--max-cells', type=int, default=300000, help='TIGHT_MAX: cells written per run (reservoir sample)')
    ap.add_argument('--dtheta', type=float, default=0.5, help='angle slab (deg) for the voxel clustering')
    ap.add_argument('--no-plot', action='store_true')
    ap.add_argument('--stop-after-dump', action='store_true')
    ap.add_argument('--min-only', action='store_true', help='only rescale and compute m(t), E(t)')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    ts = [x for x in args.ts.split(',') if x]
    rows = []
    for tstr in ts:
        tag = os.path.join(args.outdir, f'tight_{tstr}')
        cert = tag + '_cert.txt'; minf = tag + '_min.txt'; dump = tag + '_cells.txt'; statf = tag + '_stats.json'
        crit = tstr == 'crit'
        if crit: cert = args.cert
        elif not os.path.exists(cert): rescale_cert(args.cert, float(tstr), cert)
        sn, sd, D, WD, _ = read_cert(cert); s = sn / sd
        if os.path.exists(minf):
            j = json.load(open(minf)); tot = tuple(j['total']); mn = tuple(j['min'])
        else:
            t0 = time.time(); tot, mn, _, out = run_verifier(cert, args.N, args.threads)
            json.dump(dict(total=tot, min=mn, N=args.N, cert=cert, s=[sn, sd], seconds=time.time() - t0, stdout=out), open(minf, 'w'), indent=1)
        Wn, m = tot[0], mn[0]
        assert tot[1] == WD and mn[1] == WD
        E = Fraction(Wn, m) - 12
        thresh = m if crit else Wn - 11 * m           # exact integer numerator: m(1+E) = W - 11 m
        print(f"[{tstr}] s = {sn}/{sd} = {s:.6f}  m = {m}/{WD} = {m / WD:.7f}  E = {float(E):.6f}  thresh = {thresh}/{WD} = {thresh / WD:.7f}", flush=True)
        if args.min_only: continue
        if not os.path.exists(dump):
            t0 = time.time()
            try:
                tot2, mn2, ncell, out = run_verifier(cert, args.N, args.threads, dump=dump, thresh=thresh, max_cells=args.max_cells, grid=args.grid)
            except Exception:
                if os.path.exists(dump): os.remove(dump)
                raise
            assert mn2 == mn and tot2 == tot, (mn2, mn)
            print(f"[{tstr}] {ncell} tight cells ({time.time() - t0:.0f}s, {os.path.getsize(dump) / 1e6:.0f} MB)", flush=True)
        if args.stop_after_dump: continue
        if os.path.exists(statf):
            st = json.load(open(statf))
        else:
            t0 = time.time()
            N, G, bins, hist, cells = parse_dump(dump)
            st = analyse(N, G, bins, hist, cells, s, dtheta_deg=args.dtheta)
            st.update(t=tstr, E=float(E), m=m, W=Wn, WD=WD, thresh=thresh, s=s, s_frac=[sn, sd], N=N)
            print(f"[{tstr}] analysed in {time.time() - t0:.0f}s: fraction {st['fraction']:.4e}, clusters {st['nclusters']}", flush=True)
            if not args.no_plot:
                plot(tstr, s, cells, bins, st, tag + '.png', f"{thresh}/{WD}")
            json.dump(st, open(statf, 'w'), indent=1)
        rows.append(st)
    if args.stop_after_dump or args.min_only or not rows: return
    write_md(rows, args)


def write_md(rows, args):
    L = []
    W, WD = rows[0]['W'], rows[0]['WD']
    L.append("# Near-tight pose set of the shipped certificate beyond its critical container\n")
    L.append(f"Generated by `search/tightset.py` (its docstring has the argument).  Certificate `{os.path.relpath(args.cert, ROOT)}`, "
             f"W = {W}/{WD} = {W / WD:.7f}, verifier net N = {args.N} (D4-symmetric set, so the sweep covers angles [0, 45 deg], "
             f"{rows[0]['nbins']} bins).  Threshold = W - 11 m exactly (integer numerators); for `crit` the threshold is m itself.\n")
    L.append("| t | container s | m(t) | E(t) = W/m - 12 | threshold W - 11m | tight cells | T_E share of pose space | angle bins present | clusters (x,y,theta) | clusters (x,y) | dominant angles deg (share) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        L.append(f"| {r['t']} | {r['s_frac'][0]}/{r['s_frac'][1]} = {r['s']:.6f} | {r['m'] / WD:.7f} | {r['E']:+.6f} | {r['thresh'] / WD:.7f} | {r['ncells']:,} | {100 * r['fraction']:.3f}% | {r['nbins_present']}/{r['nbins']} | {r['nclusters']} | {r['nclusters_xy']} | {fmt_intervals(r['angle_intervals'])} |")
    L.append("")
    L.append("Columns.  `m(t)`: the verifier's exact minimum captured weight over all placements at container t (shrunk-square net: a lower bound on the true minimum).  "
             "`tight cells`: arrangement cells (over all angle bins) with captured weight <= threshold.  `T_E share`: sum over cells of (area clipped to the admissible centre box) x (bin width), "
             "divided by the admissible pose space sum_k (s - w_k)^2 x width_k over [0, 45 deg] (the cell areas are clipped to the rotated centre box in f64; the Python/Rust column below re-clips the sampled cells independently).  `clusters`: 26-connected components of the (x, y, theta) voxels (box %.3f, %.1f deg) that carry T_E measure; components come in C4 orbits "
             "(rotating the container by 90 deg preserves the bin angle), so the number of distinct clusters up to symmetry is roughly a quarter.  `dominant angles`: maximal runs of consecutive bins present, with their share of the measure.\n"
             % (rows[0]['s'] / rows[0]['grid'], rows[0]['dtheta_deg']))
    L.append("## Spatial distribution\n")
    L.append("Shares of the T_E measure by centre position (folded into the quadrant [0, s/2]^2 by the C4 rotations):\n")
    L.append("| t | centre within 0.75 of two walls (corner) | within 0.75 of one wall only | both folded coords > s/2 - 0.75 (middle) | largest cluster | 4 largest | Python/Rust clipped area on the sample | mean angle (deg) |")
    L.append("|---|---|---|---|---|---|---|---|")
    for r in rows:
        cm = r['cluster_measure']; tot = r['measure'] or 1.0
        L.append(f"| {r['t']} | {100 * r['near_corner']:.1f}% | {100 * r['near_wall']:.1f}% | {100 * r['near_centre']:.1f}% | {100 * cm[0] / tot:.1f}% | {100 * sum(cm[:4]) / tot:.1f}% | "
                 f"{r['py_over_rust_area'] if r['py_over_rust_area'] is None else f'{r['py_over_rust_area']:.6f}'} | {r['mean_deg']:.2f} |")
    L.append("")
    L.append("Angle histogram of the T_E measure (1-degree buckets, percent; `.` = below 0.05%):\n")
    L.append('| t | ' + ' | '.join(f"{int(e)}" for e in rows[0]['hist_edges'][:-1]) + ' |')
    L.append('|' + '---|' * len(rows[0]['hist_edges']))
    for r in rows:
        tot = r['measure'] or 1.0
        L.append(f"| {r['t']} | " + ' | '.join(f"{100 * h / tot:.1f}" if h / tot >= 0.0005 else '.' for h in r['hist_deg']) + ' |')
    L.append("")
    L.append("Distance of the centre from the nearest wall (bands of width s/G starting at 0.5, percent of the measure):\n")
    wb = rows[0]['wall_bands']
    L.append('| t | ' + ' | '.join(f"{wb[i]:.2f}" for i in range(len(wb) - 1)) + ' |'); L.append('|' + '---|' * len(wb))
    for r in rows:
        L.append(f"| {r['t']} | " + ' | '.join(f"{100 * h:.1f}" if h >= 0.0005 else '.' for h in r['wall_hist']) + ' |')
    L.append("")
    L.append("Marginal of the folded centre coordinate (x and y pooled; boxes of s/G from 0 to s/2, percent of the measure):\n")
    G = rows[0]['grid']; s = rows[0]['s']
    L.append('| t | ' + ' | '.join(f"{(i + 0.5) * s / G:.2f}" for i in range(G // 2)) + ' |'); L.append('|' + '---|' * (G // 2 + 1))
    for r in rows:
        L.append(f"| {r['t']} | " + ' | '.join(f"{100 * h:.1f}" if h >= 0.0005 else '.' for h in r['marg']) + ' |')
    L.append("")
    L.append("Plots: `runs/tight_<t>.png` -- left: sampled cell midpoints coloured by bin angle (marker size ~ sqrt of the cell area), right: the cap-independent measure map; red: C4 quadrant and D4 fundamental triangle.\n")
    # keep a hand-written reading of the results, if any, below a '## Reading' heading
    tail = ''
    if os.path.exists(args.md):
        old = open(args.md).read()
        if '\n## Reading' in old: tail = old[old.index('\n## Reading'):]
    open(args.md, 'w').write('\n'.join(L) + '\n' + tail)
    print('wrote', args.md)


if __name__ == '__main__':
    main()
