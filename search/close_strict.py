#!/usr/bin/env python3
"""close_strict.py -- closing loop whose separation oracle IS the strict protocol (task s45-close, 2026-09-25;
search/S45_COVER.md sec 9).  HEURISTIC throughout (float LP over sampled rows, float scans).

rung2_close.py separates with `closed4.stress` (pitch 0.003, integer + near-tile angles) and prunes rows, so a hole
it closed can reopen, and the strict minimum (`s21_cover_eval.py hardscan`, esp. at pitch 0.002) is only measured
after the fact.  Here every round runs the strict protocol itself on the LP weights and keeps what it finds:

  loop CERT TAG [--cols-from CERT2,...] [--rows-from POSES,...] [--lazy-seed POOLS] [--hs-pitch 0.002] ...
      columns = the points of CERT (+ --cols-from), fixed (no pricing), D4 orbits;
      initial rows = 0.01 lattice (closed4.angle_list(0.5)) at the CERT weights, thr 1 + --init-thr, plus
      (--init-hs PITCH) the hardscan-angle lattice poses with CERT weight < 1 + --init-hs-thr, thinned, plus
      --rows-from poses (PERMANENT) plus lazy-pool violations;
      each round: IPM LP solve -> export runs/TAG_rR.txt -> strict hardscan at --hs-pitch on the LP weights
      (exactly the s21_cover_eval.hardscan protocol: 0, 0.02..3 step 0.02, 3..45 step 0.25 deg; 30 worst per angle,
      polish the 1400 worst with the same chunking and seeds) -> item-3 family at pitch 0.0005.
      PERMANENT rows (never pruned): the polished dip poses, the lattice seeds, the worst --nvisit polish-visited
      poses (the polished neighbourhood), the worst --nfam family poses; all appended to runs/TAG_dips.txt.
      Prunable rows: per-angle thinned lattice violations, lazy-pool violations, 0.01-lattice re-separation.
      Strict min of round R = min(hardscan polished, hardscan lattice, family) on the LP weights (the exported
      file is rounded UP, so its strict min is >= this: the logged cost total/min is conservative).
      Optionally (--eval-004) each checkpoint also gets the standard `s21_cover_eval.py hardscan` (pitch 0.004)
      in a background child process on --eval-nproc workers (overlaps the next LP solve).
      Stops when the strict min of the last --stable rounds spans < --stable-tol (relative), at --time, or --rounds.
  harvest OUT CERT... [--pitch 0.004]
      run the same strict hardscan on each CERT and write every dip pose (polished, seeds, worst visited) to OUT
      (`cx cy theta_rad value` rows): a --rows-from file, e.g. from the checkpoints of an earlier closing loop.
  confirm CERT DIPS [--pitch 0.001] [--box 0.1] [--dtilt 0.3] [--dnear 0.05] [--full-angles K]
      denser confirmation scan around the dips (DIPS: `cx cy theta_rad [value]` rows): a pitch-0.001 lattice in a
      +-box c-space box around each dip at angles dip +- dtilt (step 0.01 deg) / +- dnear for theta < 3 deg
      (step 0.002 deg), plus (--full-angles 'deg,lo:hi:step,...') a whole-container pitch-0.001 scan (1 x 1 tiles)
      at the given angles; then polish of the 1400 worst.  Estimates the scan-vs-true gap.
"""
import sys, os, math, time, json, argparse, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import multiprocessing as mp
import closed4 as C
import family_rows as F

HERE = os.path.dirname(os.path.abspath(__file__))


def hs_degs():
    return [0.0] + list(np.arange(0.02, 3.0, 0.02)) + list(np.arange(3.0, 45.0 + 1e-9, 0.25))


def _hs_scan(args):
    """one angle: the 30 worst lattice poses (as s21_cover_eval._scan) + thinned violated poses (rows)."""
    P, w, s, th, pitch, perang, block = args
    v, p = C.scan_angle(P, w, s, th, pitch)
    if len(v) == 0: return np.zeros(0), np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0)
    o = np.argsort(v)[:30]
    bad = np.nonzero(v < 1 - 1e-9)[0]
    bad = bad[np.argsort(v[bad], kind='stable')]
    key = np.floor(p[bad, :2] / block).astype(np.int64); kk = key[:, 0] * 100000 + key[:, 1]
    _, first = np.unique(kk, return_index=True)
    sel = bad[np.sort(first)][:perang]
    return v[o], p[o], p[sel], v[sel]


def _hs_pol(args):
    P, w, s, seeds, seed = args
    fp, fv, cur, curv = C.polish(P, w, s, seeds, rounds=9, nper=64, rng=np.random.default_rng(seed))
    return fp, fv, cur, curv


def _cap(args):
    P, w, Q = args
    return C.captured(P, w, Q)


def strict_scan(P, w, s, pool, nproc, pitch=0.002, perang=20, block=0.05, nvisit=1500):
    """the s21_cover_eval.hardscan protocol on weights w, returning the rows it found.  Returns dict."""
    t0 = time.time(); degs = hs_degs()
    outs = pool.map(_hs_scan, [(P, w, s, math.radians(d), pitch, perang, block) for d in degs])
    V = np.concatenate([o[0] for o in outs]); Q = np.concatenate([o[1] for o in outs])
    angle_min = np.array([o[0].min() if len(o[0]) else 9e9 for o in outs])
    lat = float(V.min()); order = np.argsort(V)[:1400]; seeds = Q[order]; seedv = V[order]
    res = pool.map(_hs_pol, [(P, w, s, ch, i) for i, ch in enumerate(np.array_split(seeds, 4 * nproc))])
    fp = np.concatenate([r[0] for r in res]); fv = np.concatenate([r[1] for r in res])
    cur = np.concatenate([r[2] for r in res]); curv = np.concatenate([r[3] for r in res])
    j = int(curv.argmin())
    vo = np.argsort(fv)[:nvisit]
    latrows = np.concatenate([o[2] for o in outs]); latv = np.concatenate([o[3] for o in outs])
    return dict(lat=lat, pol=float(curv[j]), pose=cur[j].copy(), min=min(lat, float(curv[j])),
                cur=cur, curv=curv, seeds=seeds, seedv=seedv, visit=fp[vo], visitv=fv[vo],
                latrows=latrows, latv=latv, degs=np.array(degs), angle_min=angle_min, t=time.time() - t0)


def family_scan(P, w, s, pool, nproc, pitch=0.0005, nfam=300):
    angles = [0.0] + [10 ** e for e in (-7, -6, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5)]
    Q = np.array(F.item3(s, pitch, angles))
    v = np.concatenate(pool.map(_cap, [(P, w, ch) for ch in np.array_split(Q, 8 * nproc)]))
    j = int(v.argmin()); bad = np.nonzero(v < 1 - 1e-9)[0]; bad = bad[np.argsort(v[bad])[:nfam]]
    return float(v[j]), Q[j], Q[bad], v[bad], int((v < 1 - 1e-9).sum())


def dip_rows(hs, fam=None):
    """permanent rows of one strict scan: polished finals < 1, seeds < 1, worst visited, family worst."""
    parts = [(hs['cur'][hs['curv'] < 1 - 1e-9], hs['curv'][hs['curv'] < 1 - 1e-9]),
             (hs['seeds'][hs['seedv'] < 1 - 1e-9], hs['seedv'][hs['seedv'] < 1 - 1e-9]),
             (hs['visit'], hs['visitv'])]
    if fam is not None: parts.append((fam[2], fam[3]))
    Pq = np.concatenate([p for p, _ in parts]); Vq = np.concatenate([v for _, v in parts])
    return Pq.reshape(-1, 3), Vq


def write_poses(path, Pq, Vq, mode='a', note=None):
    with open(path, mode) as f:
        if note: f.write(f"# {note}\n")
        for (cx, cy, th), v in zip(Pq, Vq): f.write(f"{cx:.9f} {cy:.9f} {th:.12f} {v:.8f}\n")


def fmt_pose(p):
    return f"({p[0]:.5f}, {p[1]:.5f}, {math.degrees(p[2]):.4f} deg)"


# ----------------------------------------------------------------------------- loop
def loop(a):
    os.makedirs('runs', exist_ok=True)
    lf = open(f"runs/{a.tag}.log", 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"close_strict.py loop {' '.join(sys.argv[2:])}   pid {os.getpid()}  {time.strftime('%F %T')}")
    C.LPSOLVER = 'ipm'
    s, P0, w0c = C.load_points(a.cert)
    m = C.Model(s, D=1000, sym=True)
    ncol = sum(m.add_point(x, y, 'cert') for x, y in P0)
    for f in (a.cols_from or '').split(','):
        if f.strip():
            s2, P2, _ = C.load_points(f.strip()); assert abs(s2 - s) < 1e-12
            ncol += sum(m.add_point(x, y, 'cols') for x, y in P2)
    wmap = {}
    for (x, y), wt in zip(P0, w0c): wmap[(int(round(x * m.D)), int(round(y * m.D)))] = wt
    w0 = np.array([wmap.get((int(round(px * m.D)), int(round(py * m.D))), 0.0) for px, py in m.P])
    log(f"[{a.tag}] s={s} columns: {ncol} orbits / {len(m.P)} atoms; CERT total {w0.sum():.6f}")
    pool = mp.get_context('fork').Pool(a.nproc)
    dipf = f"runs/{a.tag}_dips.txt"
    prot = []
    def add(poses, protect):
        poses = np.asarray(poses, dtype=float).reshape(-1, 3)
        if protect and len(prot):            # a dip that is already a (prunable) row becomes permanent
            idx = {m.canon(*rw): i for i, rw in enumerate(m.rows)}
            for q in poses:
                i = idx.get(m.canon(*q))
                if i is not None: prot[i] = True
        n = m.add_rows(poses); prot.extend([protect] * n); return n
    # permanent rows from files
    for f in (a.rows_from or '').split(','):
        if f.strip():
            R = C.read_poses(f.strip(), s); n = add(R, True)
            log(f"[{a.tag}] --rows-from {f.strip()}: {len(R)} poses, {n} new permanent rows")
    thetas = [math.radians(d) for d in C.angle_list(0.5)]
    gmin, nviol, poses, _ = C.separate(m.P, w0, s, thetas, 0.01, thr=1.0 + a.init_thr, perang=800, block=0.03, pool=pool)
    n = add(poses, False)
    log(f"[{a.tag}] initial 0.01-lattice rows at the CERT weights (thr {1 + a.init_thr}): {n}, lattice min {gmin:.6f}")
    if a.init_hs:
        # rows where the CERT weights are tight on the strict-protocol lattice (the hardscan angle set at pitch
        # --init-hs): a warm start that keeps the previous loop's closed holes closed without its rows
        gmin, nviol, poses, _ = C.separate(m.P, w0, s, [math.radians(d) for d in hs_degs()], a.init_hs,
                                           thr=1.0 + a.init_hs_thr, perang=a.init_hs_perang, block=a.init_hs_block, pool=pool)
        n = add(poses, False)
        log(f"[{a.tag}] initial hardscan-lattice rows (pitch {a.init_hs}, thr {1 + a.init_hs_thr}, {a.init_hs_perang}/angle, "
            f"block {a.init_hs_block}): {n}, lattice min {gmin:.6f}")
    lazy = np.zeros((0, 3))
    for f in (a.lazy_seed or '').split(','):
        if f.strip(): lazy = np.r_[lazy, C.read_poses(f.strip(), s)]
    if len(lazy):
        Sv, smin = C.pool_violations(m.P, w0, lazy, pool, a.nproc, cap=3000); n = add(Sv, False)
        log(f"[{a.tag}] lazy pool {len(lazy)} poses: {n} rows (min {smin:.6f} at the CERT weights)")
    log(f"[{a.tag}] rows {len(m.rows)} ({sum(prot)} permanent)")

    t0 = time.time(); hist = []; kids = []
    for r in range(a.rounds):
        tl = time.time(); out = m.solve()
        if out is None: log("LP failed"); break
        val, x, y, Ax = out; wt = x[m.own]; tlp = time.time() - tl
        ck = f"runs/{a.tag}_r{r}.txt"
        tot, npts = C.export(m, x, ck, WD=10 ** 7, up=True)
        if a.eval_004:
            lg = open(f"runs/{a.tag}_eval_r{r}.log", 'w')
            kids.append(subprocess.Popen([sys.executable, os.path.join(HERE, 's21_cover_eval.py'), 'hardscan', ck,
                                          '--nproc', str(a.eval_nproc)], stdout=lg, stderr=subprocess.STDOUT))
        hs = strict_scan(m.P, wt, s, pool, a.nproc, pitch=a.hs_pitch, perang=a.perang, block=a.block, nvisit=a.nvisit)
        tf = time.time(); fam = family_scan(m.P, wt, s, pool, a.nproc, nfam=a.nfam); tfam = time.time() - tf
        smin = min(hs['min'], fam[0])
        worst3 = np.argsort(hs['curv'])[:5]
        rec = dict(round=r, t=round(time.time() - t0), LP=val, total=tot, points=npts, rows=len(m.rows), perm=int(sum(prot)),
                   hs_lat=hs['lat'], hs_pol=hs['pol'], hs_pose=[float(hs['pose'][0]), float(hs['pose'][1]), math.degrees(hs['pose'][2])],
                   fam=fam[0], fam_pose=[float(fam[1][0]), float(fam[1][1]), math.degrees(fam[1][2])], fam_nviol=fam[4],
                   strict=smin, cost=tot / smin, t_lp=round(tlp), t_hs=round(hs['t']), t_fam=round(tfam),
                   worst=[(float(hs['curv'][k]), float(hs['cur'][k, 0]), float(hs['cur'][k, 1]), math.degrees(hs['cur'][k, 2])) for k in worst3])
        hist.append(rec)
        log(f"[{a.tag}] round {r} LP={val:.6f} total={tot:.6f} rows={len(m.rows)} ({rec['perm']} perm) | hardscan p{a.hs_pitch} "
            f"lat {hs['lat']:.7f} pol {hs['pol']:.7f} at {fmt_pose(hs['pose'])} | family {fam[0]:.7f} ({fam[4]} < 1) | "
            f"strict {smin:.7f} cost {tot / smin:.4f} | LP {tlp:.0f}s hs {hs['t']:.0f}s fam {tfam:.0f}s t={time.time() - t0:.0f}s")
        for k in worst3[1:]: log(f"      {hs['curv'][k]:.7f} at {fmt_pose(hs['cur'][k])}")
        json.dump(dict(s=s, args=vars(a), hist=hist), open(f"runs/{a.tag}.json", 'w'), indent=1)
        # rows
        Dp, Dv = dip_rows(hs, fam)
        write_poses(dipf, Dp, Dv, note=f"round {r} strict {smin:.7f}")
        n1 = add(Dp, True)
        n2 = add(hs['latrows'], False)
        n3 = 0
        if len(lazy):
            Sv, _ = C.pool_violations(m.P, wt, lazy, pool, a.nproc, cap=3000); n3 = add(Sv, False)
        g2, nv2, poses, _ = C.separate(m.P, wt, s, thetas, 0.01, thr=1.0 - 1e-7, perang=100, block=0.08, pool=pool)
        n4 = add(poses, False)
        log(f"   +rows {n1} permanent (dips {len(Dp)}) {n2} hardscan-lattice {n3} lazy {n4} 0.01-lattice (min {g2:.6f}); rows now {len(m.rows)}")
        if a.prune_at and len(m.rows) > a.prune_at:
            keep = np.r_[(y > 1e-12) | (Ax <= 1 + a.prune_slack), np.zeros(len(m.rows) - len(y), dtype=bool)]
            keep[-a.prune_keep:] = True; keep |= np.array(prot, dtype=bool)
            m.prune(keep); prot = [p_ for p_, k in zip(prot, keep) if k]
            log(f"   pruned to {len(m.rows)} rows ({sum(prot)} permanent)")
        if len(hist) >= a.stable:
            last = [h['strict'] for h in hist[-a.stable:]]
            if (max(last) - min(last)) / max(last) < a.stable_tol:
                log(f"[{a.tag}] STABLE: strict min of the last {a.stable} rounds {', '.join(f'{v:.6f}' for v in last)}"); break
        if smin >= 1 - 1e-7: log(f"[{a.tag}] strict min >= 1"); break
        if time.time() - t0 > a.time: log(f"[{a.tag}] time limit"); break
    pool.close(); pool.join()
    for k in kids: k.wait()
    b = min(hist, key=lambda h: h['cost'])
    log(f"[{a.tag}] DONE {len(hist)} rounds; best cost {b['cost']:.4f} at round {b['round']} (strict {b['strict']:.7f}, total {b['total']:.6f})")


# ----------------------------------------------------------------------------- harvest
def harvest(a):
    pool = mp.get_context('fork').Pool(a.nproc)
    open(a.out, 'w').close()
    for c in a.certs:
        s, P, w = C.load_points(c)
        hs = strict_scan(P, w, s, pool, a.nproc, pitch=a.pitch, perang=a.perang, block=a.block, nvisit=a.nvisit)
        Dp, Dv = dip_rows(hs)
        write_poses(a.out, Dp, Dv, note=f"{c} hardscan p{a.pitch} min {hs['min']:.7f}")
        print(f"{c}: total {w.sum():.6f} hardscan p{a.pitch} lat {hs['lat']:.7f} pol {hs['pol']:.7f} at {fmt_pose(hs['pose'])}; "
              f"{len(Dp)} dip poses; {hs['t']:.0f}s", flush=True)
    pool.close(); pool.join()


# ----------------------------------------------------------------------------- confirm
def scan_box(P, w, s, th, pitch, x0, x1, y0, y1, tol=None):
    """closed4.scan_angle restricted to centres in [x0,x1] x [y0,y1] (and the open admissible box, offset WALL);
    rotated-frame lattice anchored at 0 (so different boxes share one lattice); plus the wall-band rows in the box."""
    if tol is None: tol = C.TOL
    ct, st = math.cos(th), math.sin(th)
    lo, hi = C.admissible_box(s, th); lo += C.WALL; hi -= C.WALL
    bx0, bx1, by0, by1 = max(x0, lo), min(x1, hi), max(y0, lo), min(y1, hi)
    if bx1 < bx0 or by1 < by0: return np.zeros(0), np.zeros((0, 3))
    near = (P[:, 0] > bx0 - 0.75) & (P[:, 0] < bx1 + 0.75) & (P[:, 1] > by0 - 0.75) & (P[:, 1] < by1 + 0.75)
    P = P[near]; w = w[near]
    q0 = P[:, 0] * ct + P[:, 1] * st; q1 = -P[:, 0] * st + P[:, 1] * ct
    cr = np.array([[bx0, by0], [bx1, by0], [bx1, by1], [bx0, by1]])
    cu0 = cr[:, 0] * ct + cr[:, 1] * st; cu1 = -cr[:, 0] * st + cr[:, 1] * ct
    g0 = np.arange(math.floor(cu0.min() / pitch), math.ceil(cu0.max() / pitch) + 1) * pitch
    g1 = np.arange(math.floor(cu1.min() / pitch), math.ceil(cu1.max() / pitch) + 1) * pitch
    G0, G1 = np.meshgrid(g0, g1, indexing='ij')
    cx = G0 * ct - G1 * st; cy = G0 * st + G1 * ct
    ok = (cx >= bx0) & (cx <= bx1) & (cy >= by0) & (cy <= by1)
    h = 0.5 + tol
    oy = np.argsort(q1); q1s = q1[oy]; wy = w[oy]; q0y = q0[oy]
    lo1 = np.searchsorted(q1s, g1 - h, 'left'); hi1 = np.searchsorted(q1s, g1 + h, 'right')
    vals = np.full((len(g0), len(g1)), 9e9)
    for i in range(len(g0)):
        if not ok[i].any(): continue
        mask = np.abs(q0y - g0[i]) <= h
        cw = np.r_[0.0, np.cumsum(np.where(mask, wy, 0.0))]
        vals[i] = cw[hi1] - cw[lo1]
    vals = vals[ok]; poses = np.c_[cx[ok], cy[ok], np.full(int(ok.sum()), th)]
    bp = []
    for fixed, lo_, hi_, axis in ((lo, by0, by1, 0), (hi, by0, by1, 0), (lo, bx0, bx1, 1), (hi, bx0, bx1, 1)):
        if (axis == 0 and bx0 <= fixed <= bx1) or (axis == 1 and by0 <= fixed <= by1):
            t = np.arange(lo_, hi_ + 1e-12, pitch)
            bp.append(np.c_[np.full(len(t), fixed), t] if axis == 0 else np.c_[t, np.full(len(t), fixed)])
    if bp:
        bp = np.concatenate(bp); bp = np.c_[bp, np.full(len(bp), th)]
        vals = np.r_[vals, C.captured(P, w, bp, tol)]; poses = np.r_[poses, bp]
    return vals, poses


def _box_worker(args):
    P, w, s, th, pitch, box = args
    v, p = scan_box(P, w, s, th, pitch, *box)
    if len(v) == 0: return np.zeros(0), np.zeros((0, 3))
    o = np.argsort(v)[:30]
    return v[o], p[o]


def _full_worker(args):
    P, w, s, th, pitch = args
    v, p = C.scan_angle(P, w, s, th, pitch)
    if len(v) == 0: return np.zeros(0), np.zeros((0, 3))
    o = np.argsort(v)[:30]
    return v[o], p[o]


def confirm(a):
    s, P, w = C.load_points(a.cert); t0 = time.time()
    D = np.loadtxt(a.dips, comments='#').reshape(-1, 4 if a.has_values else 3)
    if a.has_values:
        D = D[np.argsort(D[:, 3])]
    # distinct dip neighbourhoods (fold to the fundamental angle range [0,45] is already the case for hardscan dips)
    keep = []
    for d in D:
        if all(abs(d[0] - k[0]) > a.box / 2 or abs(d[1] - k[1]) > a.box / 2 or abs(math.degrees(d[2] - k[2])) > a.dtilt / 2 for k in keep):
            keep.append(d)
        if len(keep) >= a.ndips: break
    jobs = []
    for d in keep:
        deg = math.degrees(d[2]); box = (d[0] - a.box, d[0] + a.box, d[1] - a.box, d[1] + a.box)
        if deg < 3.0: angs = np.arange(max(0.0, deg - a.dnear), deg + a.dnear + 1e-12, 0.002)
        else: angs = np.arange(deg - a.dtilt, min(45.0, deg + a.dtilt) + 1e-12, 0.01)
        if deg < 3.0 and deg - a.dnear > 0: angs = np.r_[0.0, angs]
        jobs += [(P, w, s, math.radians(t), a.pitch, box) for t in angs]
    print(f"confirm {a.cert}: total {w.sum():.6f}; {len(keep)} dip neighbourhoods (box +-{a.box}), {len(jobs)} box-angle scans "
          f"at pitch {a.pitch}", flush=True)
    with mp.get_context('fork').Pool(a.nproc) as pool:
        outs = pool.map(_box_worker, jobs, chunksize=4)
        V = np.concatenate([o[0] for o in outs]); Q = np.concatenate([o[1] for o in outs])
        print(f"   dip boxes: lattice min {V.min():.7f} at {fmt_pose(Q[V.argmin()])}  ({time.time() - t0:.0f}s)", flush=True)
        if a.full_angles:
            fa = []
            for tok in a.full_angles.split(','):     # degrees, or lo:hi:step ranges
                if ':' in tok:
                    lo_, hi_, st_ = map(float, tok.split(':')); fa += list(np.arange(lo_, hi_ + 1e-9, st_))
                else: fa.append(float(tok))
            fa = sorted(set(round(min(max(t, 0.0), 45.0), 6) for t in fa))
            # whole container, tiled into 1 x 1 boxes (the pitch-0.001 full lattice at 45 deg would need ~5 GB per worker)
            tiles = [(i, i + 1.0, j, j + 1.0) for i in range(int(s)) for j in range(int(s))]
            outs2 = pool.map(_box_worker, [(P, w, s, math.radians(t), a.pitch, b) for t in fa for b in tiles], chunksize=8)
            V2 = np.concatenate([o[0] for o in outs2]); Q2 = np.concatenate([o[1] for o in outs2])
            print(f"   full container at {len(fa)} angles: lattice min {V2.min():.7f} at {fmt_pose(Q2[V2.argmin()])} "
                  f"({time.time() - t0:.0f}s)", flush=True)
            V = np.r_[V, V2]; Q = np.r_[Q, Q2]
        seeds = Q[np.argsort(V)[:1400]]
        res = pool.map(_hs_pol, [(P, w, s, ch, i) for i, ch in enumerate(np.array_split(seeds, 4 * a.nproc))])
    cur = np.concatenate([r[2] for r in res]); curv = np.concatenate([r[3] for r in res])
    j = int(curv.argmin()); mn = min(float(V.min()), float(curv[j]))
    print(f"{a.cert}: confirm pitch {a.pitch}: lattice min {V.min():.7f}; polished min {curv[j]:.7f} at {fmt_pose(cur[j])}; "
          f"cost total/min = {w.sum() / mn:.6f}; {time.time() - t0:.0f}s")
    for k in np.argsort(curv)[:8]: print(f"   {curv[k]:.7f} at {fmt_pose(cur[k])}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='mode', required=True)
    L = sub.add_parser('loop'); L.add_argument('cert'); L.add_argument('tag')
    L.add_argument('--cols-from', default=None); L.add_argument('--rows-from', default=None); L.add_argument('--lazy-seed', default=None)
    L.add_argument('--hs-pitch', type=float, default=0.002); L.add_argument('--perang', type=int, default=20)
    L.add_argument('--block', type=float, default=0.05); L.add_argument('--nvisit', type=int, default=1500)
    L.add_argument('--nfam', type=int, default=300); L.add_argument('--init-thr', type=float, default=0.02)
    L.add_argument('--prune-at', type=int, default=0); L.add_argument('--prune-keep', type=int, default=20000)
    L.add_argument('--prune-slack', type=float, default=0.003)
    L.add_argument('--rounds', type=int, default=40); L.add_argument('--time', type=float, default=36000)
    L.add_argument('--stable', type=int, default=3); L.add_argument('--stable-tol', type=float, default=0.001)
    L.add_argument('--init-hs', type=float, default=0.0, help='also initial rows from the hardscan angle set at this pitch')
    L.add_argument('--init-hs-thr', type=float, default=0.005); L.add_argument('--init-hs-perang', type=int, default=300)
    L.add_argument('--init-hs-block', type=float, default=0.02)
    L.add_argument('--eval-004', action='store_true'); L.add_argument('--eval-nproc', type=int, default=6)
    L.add_argument('--nproc', type=int, default=7)
    H = sub.add_parser('harvest'); H.add_argument('out'); H.add_argument('certs', nargs='+')
    H.add_argument('--pitch', type=float, default=0.004); H.add_argument('--perang', type=int, default=20)
    H.add_argument('--block', type=float, default=0.05); H.add_argument('--nvisit', type=int, default=1500)
    H.add_argument('--nproc', type=int, default=7)
    K = sub.add_parser('confirm'); K.add_argument('cert'); K.add_argument('dips')
    K.add_argument('--pitch', type=float, default=0.001); K.add_argument('--box', type=float, default=0.1)
    K.add_argument('--dtilt', type=float, default=0.3); K.add_argument('--dnear', type=float, default=0.05)
    K.add_argument('--ndips', type=int, default=60); K.add_argument('--full-angles', default=None)
    K.add_argument('--no-values', dest='has_values', action='store_false'); K.add_argument('--nproc', type=int, default=7)
    a = ap.parse_args()
    {'loop': loop, 'harvest': harvest, 'confirm': confirm}[a.mode](a)


if __name__ == '__main__':
    main()
