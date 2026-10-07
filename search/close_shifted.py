#!/usr/bin/env python3
"""close_shifted.py -- closing loop measured by the strict protocol, separated on SHIFTED lattices
(task s32-close, 2026-09-25; search/S32_COVER.md sec 9).  HEURISTIC throughout (float LP over sampled rows, float scans).

Why: `rung2_close.py` separates with `closed4.stress` and prunes rows, so a closed hole can reopen, and the strict
minimum (`s21_cover_eval.py hardscan`, pitch 0.002) is only measured after the fact.  Here every round
  1. MEASURES the exported round cover with exactly the strict protocol: `s21_cover_eval.hardscan` at --hs-pitch
     (angles 0, 0.02..3 step 0.02, 3..45 step 0.25 deg; 30 worst per angle; polish of the 1400 worst in 28 chunks
     with seeds 0..27, i.e. the chunking of the default --nproc 7) and `s21_cover_eval.family` (pitch 0.0005);
     strict min = min of the three; cost = file total / strict min;
  2. SEPARATES on a lattice of the same pitch but a random sub-pitch offset and randomly jittered angles (a fresh
     draw every round), so the rows do not sit on the measurement lattice: per angle the violated poses thinned to the
     worst per --block cell (at most --perang), permanent when below --protect;
  3. keeps as PERMANENT rows (never pruned) every dip below --protect (0.995) among: the worst --npol polished finals
     of the measurement (each fixes a whole cell of the piecewise-constant capture function), the worst --nvisit
     polish-visited poses, the worst --nfam family poses, the --protect-perang worst shifted-lattice poses per angle;
     plus a stencil around the --nstencil deepest distinct polished dips.  The permanent dips are appended to
     runs/TAG_dips.txt (a --rows-from file for a restart).  Prunable: the same kinds of poses in [--protect, 1),
     the other shifted-lattice violations, lazy-pool violations, the 0.01-lattice re-separation;
  4. re-solves the LP (fixed columns: CERT + --cols-from; HiGHS IPM by default, --simplex for warm dual simplex).
Round -1 is the CERT itself (its measurement must reproduce `s21_cover_eval.py hardscan --pitch P` + `family`).
Stops when the strict min changed by < --stable-tol (relative) over --stable consecutive rounds, at --time, or --rounds.

  loop CERT TAG [--cols-from C2,..] [--rows-from POSES,..] [--lazy-seed POOLS] [--hs-pitch 0.002] [--eval-004] ...
  confirm CERT [--dips runs/TAG_dips.txt] [--pitch 0.001] [--full] [--interleaved] [--tiles] [--box 0.08] [--ndips 80]
      denser confirmation scan: (--full) the whole container at pitch 0.001 at every hardscan angle, (--interleaved) the
      same at the angles half-way between them, (--tiles) every tile pose +-0.015 at pitch 0.0005 and 0..1.5 deg step
      0.01, plus pitch-0.001
      boxes (+-box) around the deepest distinct dips of DIPS and of the CERT's own hardscan-0.002 at fine angles
      (+-0.3 deg step 0.01 above 3 deg; +-0.06 deg step 0.002 below), then polish of the 1400 worst.
Outputs: runs/TAG.log, runs/TAG.json, runs/TAG_r{R}.txt (every round, certificate format, weights rounded UP),
runs/TAG_dips.txt, runs/TAG_eval004_r{R}.log (--eval-004: s21_cover_eval.py hardscan at 0.004, background child).
"""
import sys, os, math, time, json, argparse, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import multiprocessing as mp
import closed4 as C
import family_rows as F

HERE = os.path.dirname(os.path.abspath(__file__))
NEAR = np.arange(0.02, 3.0, 0.02)
TILT = np.arange(3.0, 45.0 + 1e-9, 0.25)


def hs_degs():
    """the s21_cover_eval.hardscan angle list (degrees)."""
    return [0.0] + list(NEAR) + list(TILT)


def scan_shift(P, w, s, th, pitch, sh=(0.0, 0.0), box=None, tol=None):
    """closed4.scan_angle with the rotated-frame lattice moved by sh (in units of the pitch, [0,1)^2) and optionally
    restricted to centres in box = (x0, x1, y0, y1); plus the four wall bands (also moved by sh[0]).  With sh = 0 and
    box = None this is exactly closed4.scan_angle."""
    if tol is None: tol = C.TOL
    ct, st = math.cos(th), math.sin(th)
    lo, hi = C.admissible_box(s, th); lo += C.WALL; hi -= C.WALL
    if hi <= lo: return np.zeros(0), np.zeros((0, 3))
    bx0, bx1, by0, by1 = (lo, hi, lo, hi) if box is None else (max(box[0], lo), min(box[1], hi), max(box[2], lo), min(box[3], hi))
    if bx1 < bx0 or by1 < by0: return np.zeros(0), np.zeros((0, 3))
    if box is not None:
        near = (P[:, 0] > bx0 - 0.75) & (P[:, 0] < bx1 + 0.75) & (P[:, 1] > by0 - 0.75) & (P[:, 1] < by1 + 0.75)
        P = P[near]; w = w[near]
    q0 = P[:, 0] * ct + P[:, 1] * st; q1 = -P[:, 0] * st + P[:, 1] * ct
    corners = np.array([[bx0, by0], [bx1, by0], [bx1, by1], [bx0, by1]])
    cu0 = corners[:, 0] * ct + corners[:, 1] * st; cu1 = -corners[:, 0] * st + corners[:, 1] * ct
    if box is None:            # anchored at the box corner, as closed4.scan_angle
        g0 = np.arange(cu0.min(), cu0.max() + 1e-12, pitch) + sh[0] * pitch
        g1 = np.arange(cu1.min(), cu1.max() + 1e-12, pitch) + sh[1] * pitch
    else:                      # anchored at 0, so that overlapping boxes share one lattice
        g0 = (np.arange(math.floor(cu0.min() / pitch), math.ceil(cu0.max() / pitch) + 1) + sh[0]) * pitch
        g1 = (np.arange(math.floor(cu1.min() / pitch), math.ceil(cu1.max() / pitch) + 1) + sh[1]) * pitch
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
    for fixed, a0, a1, ax in ((lo, by0, by1, 0), (hi, by0, by1, 0), (lo, bx0, bx1, 1), (hi, bx0, bx1, 1)):
        if (ax == 0 and bx0 <= fixed <= bx1) or (ax == 1 and by0 <= fixed <= by1):
            if box is None:
                t = np.arange(lo, hi + 1e-12, pitch)
                if t[-1] < hi - 1e-9: t = np.r_[t, hi]
                if sh[0]: t = np.clip(t + sh[0] * pitch, lo, hi)
            else:
                t = np.arange(a0, a1 + 1e-12, pitch)
            bp.append(np.c_[np.full(len(t), fixed), t] if ax == 0 else np.c_[t, np.full(len(t), fixed)])
    if bp:
        bp = np.concatenate(bp); bp = np.c_[bp, np.full(len(bp), th)]
        vals = np.r_[vals, C.captured(P, w, bp, tol)]; poses = np.r_[poses, bp]
    return vals, poses


def _worst30(args):
    P, w, s, th, pitch, sh, box = args
    v, p = scan_shift(P, w, s, th, pitch, sh, box)
    if len(v) == 0: return np.zeros(0), np.zeros((0, 3))
    o = np.argsort(v)[:30]
    return v[o], p[o]


def _sep(args):
    """shifted-lattice separation at one angle: poses with value < thr, worst per `block` cell, at most `perang`,
    worst first."""
    P, w, s, th, pitch, sh, perang, block, thr = args
    v, p = scan_shift(P, w, s, th, pitch, sh)
    bad = np.nonzero(v < thr)[0]
    if len(bad) == 0: return np.zeros((0, 3)), np.zeros(0), (float(v.min()) if len(v) else 9e9)
    bad = bad[np.argsort(v[bad], kind='stable')]
    key = np.floor(p[bad, :2] / block).astype(np.int64); kk = key[:, 0] * 100000 + key[:, 1]
    _, first = np.unique(kk, return_index=True)
    sel = bad[np.sort(first)][:perang]
    return p[sel], v[sel], float(v.min())


def _pol(args):
    P, w, s, seeds, seed = args
    return C.polish(P, w, s, seeds, rounds=9, nper=64, rng=np.random.default_rng(seed))


def _cap(args):
    P, w, Q = args
    return C.captured(P, w, Q)


def measure(P, w, s, pool, pitch=0.002, nvisit=1500, nfam=300, family=True):
    """the strict protocol (s21_cover_eval.hardscan at `pitch` + family), keeping what it found."""
    t0 = time.time()
    outs = pool.map(_worst30, [(P, w, s, math.radians(d), pitch, (0.0, 0.0), None) for d in hs_degs()])
    V = np.concatenate([o[0] for o in outs]); Q = np.concatenate([o[1] for o in outs])
    lat = float(V.min()); seeds = Q[np.argsort(V)[:1400]]
    res = pool.map(_pol, [(P, w, s, ch, i) for i, ch in enumerate(np.array_split(seeds, 28))])
    fp = np.concatenate([r[0] for r in res]); fv = np.concatenate([r[1] for r in res])
    cur = np.concatenate([r[2] for r in res]); curv = np.concatenate([r[3] for r in res])
    j = int(curv.argmin()); vo = np.argsort(fv)[:nvisit]
    out = dict(lat=lat, pol=float(curv[j]), pose=cur[j].copy(), cur=cur, curv=curv, visit=fp[vo], visitv=fv[vo],
               t_hs=time.time() - t0)
    fmin, fpose, fbad, fbadv, fn = 9e9, np.zeros(3), np.zeros((0, 3)), np.zeros(0), 0
    if family:
        t1 = time.time()
        angles = [0.0] + [10 ** e for e in (-7, -6, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5)]
        FQ = np.array(F.item3(s, 0.0005, angles))
        v = np.concatenate(pool.map(_cap, [(P, w, ch) for ch in np.array_split(FQ, 56)]))
        k = int(v.argmin()); bad = np.nonzero(v < 1 - 1e-9)[0]; fn = len(bad); bad = bad[np.argsort(v[bad])[:nfam]]
        fmin, fpose, fbad, fbadv = float(v[k]), FQ[k], FQ[bad], v[bad]
        out['t_fam'] = time.time() - t1
    out.update(fam=fmin, fam_pose=fpose, fam_rows=fbad, fam_rowsv=fbadv, fam_n=fn,
               hs=min(lat, out['pol']), strict=min(lat, out['pol'], fmin))
    return out


def stencil(poses, vals, k, h, a):
    """(3x3 in c-space at +-h) x {0} + centre x {+-a}: 10 extra poses around each of the k deepest distinct dips."""
    o = np.argsort(vals); pick = []
    for i in o:
        q = poses[i]
        if all(abs(q[0] - p[0]) > 2 * h or abs(q[1] - p[1]) > 2 * h or abs(q[2] - p[2]) > 2 * a for p in pick): pick.append(q)
        if len(pick) >= k: break
    out = []
    for q in pick:
        for dx in (-h, 0, h):
            for dy in (-h, 0, h):
                if dx or dy: out.append((q[0] + dx, q[1] + dy, q[2]))
        out += [(q[0], q[1], q[2] + a), (q[0], q[1], max(0.0, q[2] - a))]
    return np.array(out).reshape(-1, 3)


def write_poses(path, Pq, Vq, note=None):
    with open(path, 'a') as f:
        if note: f.write(f"# {note}\n")
        for (cx, cy, th), v in zip(Pq, Vq): f.write(f"{cx:.9f} {cy:.9f} {th:.12f} {v:.8f}\n")


def fmt(p):
    return f"({p[0]:.5f}, {p[1]:.5f}, {math.degrees(p[2]):.4f} deg)"


def clip_poses(Q, s):
    Q = np.array(Q, dtype=float).reshape(-1, 3)
    for i in range(len(Q)):
        lo, hi = C.admissible_box(s, Q[i, 2]); lo += C.WALL; hi -= C.WALL
        Q[i, 0] = min(max(Q[i, 0], lo), hi); Q[i, 1] = min(max(Q[i, 1], lo), hi)
    return Q


# ----------------------------------------------------------------------------- loop
def loop(a):
    os.makedirs('runs', exist_ok=True)
    lf = open(f"runs/{a.tag}.log", 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"close_shifted.py loop {' '.join(sys.argv[2:])}   pid {os.getpid()}  {time.strftime('%F %T')}")
    if not a.simplex: C.LPSOLVER = 'ipm'
    rng = np.random.default_rng(a.seed)
    s, P0, w0c = C.load_points(a.cert)
    m = C.Model(s, D=1000, sym=True)
    ncol = sum(m.add_point(x, y, 'cert') for x, y in P0)
    for f in (a.cols_from or '').split(','):
        if f.strip():
            s2, P2, _ = C.load_points(f.strip()); assert abs(s2 - s) < 1e-12
            ncol += sum(m.add_point(x, y, 'cols') for x, y in P2)
    wmap = {(int(round(x * m.D)), int(round(y * m.D))): wt for (x, y), wt in zip(P0, w0c)}
    wt = np.array([wmap.get((int(round(px * m.D)), int(round(py * m.D))), 0.0) for px, py in m.P])
    log(f"[{a.tag}] s={s} columns: {ncol} orbits / {len(m.P)} atoms; CERT total {wt.sum():.6f}; LP {C.LPSOLVER}")
    pool = mp.get_context('fork').Pool(a.nproc)
    dipf = f"runs/{a.tag}_dips.txt"; prot = []
    def add(poses, protect):
        poses = np.asarray(poses, dtype=float).reshape(-1, 3)
        if protect and len(prot):            # a pose that is already a prunable row becomes permanent
            idx = {m.canon(*rw): i for i, rw in enumerate(m.rows)}
            for q in poses:
                i = idx.get(m.canon(*q))
                if i is not None: prot[i] = True
        n = m.add_rows(poses); prot.extend([protect] * n); return n
    for f in (a.rows_from or '').split(','):
        if f.strip():
            R = C.read_poses(f.strip(), s); n = add(R, True)
            log(f"[{a.tag}] --rows-from {f.strip()}: {len(R)} poses, {n} new permanent rows")
    thetas = [math.radians(d) for d in C.angle_list(0.5)]
    gmin, _, poses, _ = C.separate(m.P, wt, s, thetas, 0.01, thr=1.0 + a.init_thr, perang=800, block=0.03, pool=pool)
    n = add(poses, False)
    log(f"[{a.tag}] initial 0.01-lattice rows at the CERT weights (thr {1 + a.init_thr}): {n}, lattice min {gmin:.6f}")
    lazy = np.zeros((0, 3))
    for f in (a.lazy_seed or '').split(','):
        if f.strip(): lazy = np.r_[lazy, C.read_poses(f.strip(), s)]

    t0 = time.time(); hist = []; kids = []; r = -1; cert = a.cert; y = Ax = None; val = float('nan'); tlp = 0.0
    while True:
        # ---- 1. measure the exported cover of this round (round -1: the CERT itself)
        sP, sw = C.load_points(cert)[1:]
        if a.eval_004 and r >= 0:
            lg = open(f"runs/{a.tag}_eval004_r{r}.log", 'w')
            kids.append(subprocess.Popen(['taskset', '-c', a.eval_cpus, sys.executable, os.path.join(HERE, 's21_cover_eval.py'),
                                          'hardscan', cert, '--nproc', str(a.eval_nproc)], stdout=lg, stderr=subprocess.STDOUT))
        ms = measure(sP, sw, s, pool, pitch=a.hs_pitch, nvisit=a.nvisit, nfam=a.nfam)
        tot = float(sw.sum()); smin = ms['strict']
        worst = np.argsort(ms['curv'])
        rec = dict(round=r, t=round(time.time() - t0), cert=cert, LP=val, total=tot, points=len(sw), rows=len(m.rows),
                   perm=int(sum(prot)), hs_lat=ms['lat'], hs_pol=ms['pol'], hs_pose=[float(ms['pose'][0]), float(ms['pose'][1]),
                   math.degrees(ms['pose'][2])], fam=ms['fam'], fam_n=ms['fam_n'],
                   fam_pose=[float(ms['fam_pose'][0]), float(ms['fam_pose'][1]), math.degrees(ms['fam_pose'][2])],
                   strict=smin, cost=tot / smin, t_lp=round(tlp), t_hs=round(ms['t_hs']), t_fam=round(ms.get('t_fam', 0)),
                   worst=[(float(ms['curv'][k]), float(ms['cur'][k, 0]), float(ms['cur'][k, 1]), math.degrees(ms['cur'][k, 2]))
                          for k in worst[:8]])
        hist.append(rec)
        log(f"[{a.tag}] round {r} LP={val:.6f} total={tot:.6f} rows={len(m.rows)} ({rec['perm']} perm) | hardscan p{a.hs_pitch} "
            f"lat {ms['lat']:.7f} pol {ms['pol']:.7f} at {fmt(ms['pose'])} | family {ms['fam']:.7f} ({ms['fam_n']} < 1) | "
            f"strict {smin:.7f} cost {tot / smin:.4f} | LP {tlp:.0f}s hs {ms['t_hs']:.0f}s fam {ms.get('t_fam', 0):.0f}s "
            f"t={time.time() - t0:.0f}s")
        seen = [ms['pose']]
        for k in worst[1:]:
            q = ms['cur'][k]
            if all(np.abs(q - p_).max() > 0.01 for p_ in seen):
                seen.append(q); log(f"      {ms['curv'][k]:.7f} at {fmt(q)}")
            if len(seen) >= 6: break
        json.dump(dict(s=s, args=vars(a), hist=hist), open(f"runs/{a.tag}.json", 'w'), indent=1)
        # ---- stop?
        if len(hist) > a.stable:
            last = [h_['strict'] for h_ in hist[-(a.stable + 1):]]
            if hist[-(a.stable + 1)]['round'] >= 0 and all(abs(u - v) / v < a.stable_tol for u, v in zip(last[1:], last[:-1])):
                log(f"[{a.tag}] STABLE: strict min over the last {a.stable + 1} rounds {', '.join(f'{v:.6f}' for v in last)}"); break
        if smin >= 1 - 1e-7: log(f"[{a.tag}] strict min >= 1"); break
        if time.time() - t0 > a.time: log(f"[{a.tag}] time limit"); break
        if r + 1 >= a.rounds: log(f"[{a.tag}] round limit"); break
        # ---- 2./3. rows
        po = np.argsort(ms['curv'])[:a.npol]; po = po[ms['curv'][po] < 1 - 1e-9]
        Dp = [ms['cur'][po], ms['visit'], ms['fam_rows']]; Dv = [ms['curv'][po], ms['visitv'], ms['fam_rowsv']]
        sm = po[ms['curv'][po] < a.protect]
        st = clip_poses(stencil(ms['cur'][sm], ms['curv'][sm], a.nstencil, a.hs_pitch / 2, math.radians(a.stencil_deg)), s)
        ts = time.time()
        sh = rng.random(2); jit = [0.0] + list(NEAR + 0.02 * rng.random(len(NEAR))) + list(TILT[:-1] + 0.25 * rng.random(len(TILT) - 1)) + [45.0]
        thr = 1 + a.init_thr if r < 0 else 1 - 1e-9     # round -1: the near-tight poses of the CERT (keeps the LP near it)
        outs = pool.map(_sep, [(m.P, wt, s, math.radians(d), a.hs_pitch, tuple(sh), a.perang * (3 if r < 0 else 1), a.block, thr)
                               for d in jit])
        LP_ = np.concatenate([o[0] for o in outs]); LV = np.concatenate([o[1] for o in outs]); smin_sh = min(o[2] for o in outs)
        # permanent: per angle the --protect-perang worst below --protect (the lists are worst first)
        lowm = np.concatenate([(np.arange(len(o[1])) < a.protect_perang) & (o[1] < a.protect) for o in outs])
        Dp.append(LP_[lowm]); Dv.append(LV[lowm])
        Dp = np.concatenate(Dp).reshape(-1, 3); Dv = np.concatenate(Dv); pm = Dv < a.protect
        write_poses(dipf, Dp[pm], Dv[pm], note=f"round {r} strict {smin:.7f}")
        n1 = add(Dp[pm], True); n1s = add(st, True)
        n2 = add(Dp[~pm], False) + add(LP_[~lowm], False)
        n3 = 0
        if len(lazy):
            Sv, _ = C.pool_violations(m.P, wt, lazy, pool, a.nproc, cap=3000); n3 = add(Sv, False)
        g2, _, poses, _ = C.separate(m.P, wt, s, thetas, 0.01, thr=1.0 - 1e-7, perang=100, block=0.08, pool=pool)
        n4 = add(poses, False)
        log(f"   +rows {n1} permanent dips (< {a.protect} of {len(Dp)}: {len(po)} polished, {len(ms['visit'])} visited, {len(ms['fam_rows'])} family, "
            f"{int(lowm.sum())} shifted) + {n1s} stencil | {n2} other dips + shifted-lattice < {thr:.4f} (shift {sh[0]:.3f},{sh[1]:.3f}; "
            f"min {smin_sh:.7f}) {n3} lazy {n4} 0.01-lattice (min {g2:.6f}) | rows {len(m.rows)} ({sum(prot)} perm) {time.time() - ts:.0f}s")
        if a.prune_at and len(m.rows) > a.prune_at and y is not None:
            keep = np.r_[(y > 1e-12) | (Ax <= 1 + a.prune_slack), np.zeros(len(m.rows) - len(y), dtype=bool)]
            keep[-a.prune_keep:] = True; keep |= np.array(prot, dtype=bool)
            m.prune(keep); prot = [p_ for p_, k in zip(prot, keep) if k]
            log(f"   pruned to {len(m.rows)} rows ({sum(prot)} permanent)")
        # ---- 4. LP
        r += 1; tl = time.time()
        out = m.solve_warm() if a.simplex else m.solve()
        if out is None: log("LP failed"); break
        val, x, y, Ax = out; wt = x[m.own]; tlp = time.time() - tl
        cert = f"runs/{a.tag}_r{r}.txt"; C.export(m, x, cert, WD=10 ** 7, up=True)
    pool.close(); pool.join()
    for k in kids: k.wait()
    b = min(hist, key=lambda h_: h_['cost'])
    log(f"[{a.tag}] DONE {len(hist)} measurements; best cost {b['cost']:.4f} at round {b['round']} "
        f"(strict {b['strict']:.7f}, total {b['total']:.6f}); last {hist[-1]['cost']:.4f} {time.strftime('%F %T')}")


# ----------------------------------------------------------------------------- confirm
def confirm(a):
    s, P, w = C.load_points(a.cert); t0 = time.time(); tot = w.sum()
    pool = mp.get_context('fork').Pool(a.nproc)
    print(f"confirm {a.cert}: total {tot:.6f}  {time.strftime('%F %T')}", flush=True)
    V = np.zeros(0); Q = np.zeros((0, 3))
    # the cover's own hardscan-0.002 dips, plus DIPS
    ms = measure(P, w, s, pool, pitch=0.002, family=False)
    print(f"   hardscan p0.002 lat {ms['lat']:.7f} pol {ms['pol']:.7f} at {fmt(ms['pose'])} ({time.time() - t0:.0f}s)", flush=True)
    D = [np.c_[ms['cur'], ms['curv']]]
    if a.dips:
        dd = np.loadtxt(a.dips, comments='#').reshape(-1, 4); D.append(dd)
    D = np.concatenate(D); D = D[np.argsort(D[:, 3])]
    keep = []
    for d in D:
        if all(abs(d[0] - k[0]) > a.box or abs(d[1] - k[1]) > a.box or abs(math.degrees(d[2] - k[2])) > 0.3 for k in keep): keep.append(d)
        if len(keep) >= a.ndips: break
    jobs = []
    for d in keep:
        deg = math.degrees(d[2]); box = (d[0] - a.box, d[0] + a.box, d[1] - a.box, d[1] + a.box)
        angs = np.arange(max(0.0, deg - 0.06), deg + 0.06 + 1e-12, 0.002) if deg < 3.0 else \
            np.arange(deg - 0.3, min(45.0, deg + 0.3) + 1e-12, 0.01)
        if deg < 3.0: angs = np.r_[0.0, angs]
        jobs += [(P, w, s, math.radians(t), a.pitch, (0.0, 0.0), box) for t in angs]
    print(f"   {len(keep)} dip neighbourhoods (box +-{a.box}), {len(jobs)} box scans at pitch {a.pitch}", flush=True)
    outs = pool.map(_worst30, jobs, chunksize=4)
    Vb = np.concatenate([o[0] for o in outs]); Qb = np.concatenate([o[1] for o in outs])
    print(f"   dip boxes: lattice min {Vb.min():.7f} at {fmt(Qb[Vb.argmin()])} ({time.time() - t0:.0f}s)", flush=True)
    V = np.r_[V, Vb]; Q = np.r_[Q, Qb]
    if a.tiles:                # every tile pose (i + 1/2, j + 1/2), centres within +-0.015 on a fine lattice, 0..1.5 deg step 0.01
        n = int(round(s)); tj = []
        for i in range(n):
            for j in range(n):
                box = (i + 0.5 - 0.015, i + 0.5 + 0.015, j + 0.5 - 0.015, j + 0.5 + 0.015)
                tj += [(P, w, s, math.radians(t), a.tile_pitch, (0.0, 0.0), box) for t in np.arange(0.0, 1.5 + 1e-9, 0.01)]
        outs = pool.map(_worst30, tj, chunksize=8)
        Vt = np.concatenate([o[0] for o in outs]); Qt = np.concatenate([o[1] for o in outs])
        print(f"   near-tile: {n * n} tiles x 151 angles 0..1.5 deg, box +-0.015 at pitch {a.tile_pitch}: lattice min {Vt.min():.7f} "
              f"at {fmt(Qt[Vt.argmin()])} ({time.time() - t0:.0f}s)", flush=True)
        V = np.r_[V, Vt]; Q = np.r_[Q, Qt]
    for name, degs in (('hardscan', hs_degs() if a.full else []),
                       ('interleaved', (list(np.arange(0.01, 3.0, 0.02)) + list(np.arange(3.125, 45.0, 0.25))) if a.interleaved else [])):
        if not degs: continue
        outs = pool.map(_worst30, [(P, w, s, math.radians(d), a.pitch, (0.0, 0.0), None) for d in degs])
        Vf = np.concatenate([o[0] for o in outs]); Qf = np.concatenate([o[1] for o in outs])
        print(f"   full container, {len(degs)} {name} angles, pitch {a.pitch}: lattice min {Vf.min():.7f} at "
              f"{fmt(Qf[Vf.argmin()])} ({time.time() - t0:.0f}s)", flush=True)
        for k in np.argsort(Vf)[:4]: print(f"      {Vf[k]:.7f} at {fmt(Qf[k])}", flush=True)
        V = np.r_[V, Vf]; Q = np.r_[Q, Qf]
    seeds = Q[np.argsort(V)[:1400]]
    res = pool.map(_pol, [(P, w, s, ch, i) for i, ch in enumerate(np.array_split(seeds, 28))])
    pool.close(); pool.join()
    cur = np.concatenate([r[2] for r in res]); curv = np.concatenate([r[3] for r in res])
    j = int(curv.argmin()); mn = min(float(V.min()), float(curv[j]), ms['hs'])
    print(f"{a.cert}: confirm pitch {a.pitch}: lattice min {V.min():.7f}; polished min {curv[j]:.7f} at {fmt(cur[j])}; "
          f"overall min (incl. hardscan 0.002) {mn:.7f}; cost total/min = {tot / mn:.6f}; {time.time() - t0:.0f}s")
    for k in np.argsort(curv)[:8]: print(f"   {curv[k]:.7f} at {fmt(cur[k])}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='mode', required=True)
    L = sub.add_parser('loop'); L.add_argument('cert'); L.add_argument('tag')
    L.add_argument('--cols-from', default=None); L.add_argument('--rows-from', default=None); L.add_argument('--lazy-seed', default=None)
    L.add_argument('--hs-pitch', type=float, default=0.002)
    L.add_argument('--perang', type=int, default=40); L.add_argument('--block', type=float, default=0.02)
    L.add_argument('--protect', type=float, default=0.995, help='shifted-lattice violations below this are permanent rows ...')
    L.add_argument('--protect-perang', type=int, default=4, help='... at most this many (the worst) per angle')
    L.add_argument('--nvisit', type=int, default=1000); L.add_argument('--npol', type=int, default=600); L.add_argument('--nfam', type=int, default=300)
    L.add_argument('--nstencil', type=int, default=40); L.add_argument('--stencil-deg', type=float, default=0.01)
    L.add_argument('--init-thr', type=float, default=0.005)
    L.add_argument('--prune-at', type=int, default=0); L.add_argument('--prune-keep', type=int, default=25000)
    L.add_argument('--prune-slack', type=float, default=0.003)
    L.add_argument('--rounds', type=int, default=60); L.add_argument('--time', type=float, default=36000)
    L.add_argument('--stable', type=int, default=3); L.add_argument('--stable-tol', type=float, default=0.001)
    L.add_argument('--eval-004', action='store_true'); L.add_argument('--eval-nproc', type=int, default=6)
    L.add_argument('--eval-cpus', default='0-6,16-22')
    L.add_argument('--simplex', action='store_true', help='warm dual simplex (closed4.Model.solve_warm) instead of cold IPM')
    L.add_argument('--seed', type=int, default=1); L.add_argument('--nproc', type=int, default=7)
    K = sub.add_parser('confirm'); K.add_argument('cert'); K.add_argument('--dips', default=None)
    K.add_argument('--pitch', type=float, default=0.001); K.add_argument('--box', type=float, default=0.08)
    K.add_argument('--ndips', type=int, default=80); K.add_argument('--full', action='store_true')
    K.add_argument('--interleaved', action='store_true', help='full container at the angles between the hardscan ones '
                   '(0.01..2.99 step 0.02, 3.125..44.875 step 0.25)')
    K.add_argument('--tiles', action='store_true', help='fine near-tile scan: every tile pose +-0.015, 0..1.5 deg step 0.01')
    K.add_argument('--tile-pitch', type=float, default=0.0005)
    K.add_argument('--nproc', type=int, default=7)
    a = ap.parse_args()
    {'loop': loop, 'confirm': confirm}[a.mode](a)


if __name__ == '__main__':
    main()
