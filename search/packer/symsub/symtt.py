#!/usr/bin/env python3
"""Time travel for the symmetric subspace (TT-S, 10-18).  Start = an earlier public packing at an n whose current
record is Md-symmetric (kingbird revisions / older issues, 0.003-0.05 above today's record).  No record information:
the start's own best diagonal mirror (Md or Ma; Ma rotated to Md) is used.

Arms (same hop budget, same seeds):
  sym : symmetrise the start (Hungarian pairing with its mirror images, average), then monotone basin hopping: kick
        (symmetrised, sigma from {0.05, 0.1, 0.2}), harden reharden --sym-md 1 (r0 0.1, box free while soft), free fq
        polish (= release).  Current state = the symmetric harden output; accept when the released side improves.
  free: the same hopping from the unsymmetrised start, free kicks, free harden, free polish.
Outcome: best released side per run vs today's record (hit: within 1e-7 or below).

  symtt.py run --out runs/symsub/tt1 [--hops 80] [--reps 4] [--procs 15]
"""
import sys, os, json, argparse, subprocess, tempfile, time, collections
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from symscan import invariant
from symharden import write_cfg, HARDEN
from symfunnel import sym_kick, partners

SIGS = (0.05, 0.1, 0.2)


def mirror(sq, s, el):
    x, y, t = sq[:, 0], sq[:, 1], sq[:, 2]
    if el == 'Md':
        return np.c_[y, x, (-t) % 90]
    return np.c_[s - y, s - x, (-t) % 90]                       # Ma


def best_diag_mirror(sq, s):
    """(element, mean matched distance) by Hungarian between squares and their images"""
    out = []
    for el in ('Md', 'Ma'):
        im = mirror(sq, s, el)
        C = np.sqrt(((sq[:, None, :2] - im[None, :, :2]) ** 2).sum(-1))
        r, c = linear_sum_assignment(C)
        out.append((C[r, c].mean(), el))
    d, el = min(out)
    return el, d


def to_md(sq, s, el):
    """rotate the box 90 degrees when the start's mirror is Ma: (x, y) -> (s - y, x) maps Ma-invariant to Md-invariant"""
    if el == 'Md':
        return sq.copy()
    return np.c_[s - sq[:, 1], sq[:, 0], (sq[:, 2] + 90) % 90]


def symmetrize(sq, s):
    im = mirror(sq, s, 'Md')
    C = np.sqrt(((sq[:, None, :2] - im[None, :, :2]) ** 2).sum(-1))
    r, c = linear_sum_assignment(C)
    out = sq.copy()
    for i, j in zip(r, c):                                       # i's mirror-matched partner image is im[j]
        out[i, :2] = 0.5 * (sq[i, :2] + im[j, :2])
        d = ((im[j, 2] - sq[i, 2] + 45) % 90) - 45
        out[i, 2] = (sq[i, 2] + 0.5 * d) % 90
    # squares matched to themselves: put on the diagonal at 0 / 45
    k = partners(out)
    for i in range(len(out)):
        if k[i] == i:
            m = 0.5 * (out[i, 0] + out[i, 1]); out[i, :2] = m
            a = out[i, 2] % 90; out[i, 2] = 45.0 if 22.5 <= a < 67.5 else 0.0
    k = partners(out)
    ok = bool((k[k] == np.arange(len(out))).all())
    return out, ok


def duplicate_half(sq, s, n, half, rng):
    """Md-symmetric packing from one half: keep the squares on one side of the diagonal (half 'lo': y < x, 'hi':
    y > x) and mirror them; squares whose centre is within 0.36 of the diagonal (|x - y| / sqrt 2 < 0.36: they cross it)
    go onto the diagonal (snapped, 0 / 45 degrees) if from the kept side.  Then fix the count to n: too many -> drop
    diagonal squares, then mirror pairs, nearest the diagonal first; too few -> add mirror pairs at random poses (the
    soft stage resolves the overlap)."""
    d = (sq[:, 0] - sq[:, 1]) / np.sqrt(2)                      # signed distance from the diagonal (lo side > 0)
    sgn = 1 if half == 'lo' else -1
    keep = sq[sgn * d >= 0.36]
    near = sq[(np.abs(d) < 0.36) & (sgn * d >= 0)]
    diag = []
    for x, y, a in near:
        m = 0.5 * (x + y); a = a % 90
        diag.append((m, m, 45.0 if 22.5 <= a < 67.5 else 0.0))
    diag = np.array(diag).reshape(-1, 3)
    img = np.c_[keep[:, 1], keep[:, 0], (-keep[:, 2]) % 90]
    cnt = len(diag) + 2 * len(keep)
    while cnt > n:
        if len(diag) and (cnt - n) % 2 == 1:
            diag = diag[:-1]; cnt -= 1
        elif len(keep):
            k = np.argmin(np.abs(keep[:, 0] - keep[:, 1])); keep = np.delete(keep, k, 0); img = np.delete(img, k, 0); cnt -= 2
        else:
            diag = diag[:-1]; cnt -= 1
    extra = []
    while cnt < n:
        if (n - cnt) % 2 == 1:
            m = rng.uniform(0.5, s - 0.5); extra.append((m, m, 0.0)); cnt += 1
        else:
            x, y = sorted(rng.uniform(0.5, s - 0.5, 2))[::-1]; a = rng.uniform(0, 90)
            extra += [(x, y, a), (y, x, (-a) % 90)]; cnt += 2
    out = np.r_[diag, keep, img, np.array(extra).reshape(-1, 3)]
    return out


def relax_once(arm, s, sq, seed, tmp, r0=0.1):
    """one symmetric (or free) soft-to-rigid relaxation from sq, no kick, then free polish -> (harden s, polished s)"""
    from pk.engine import fq_quench, POLISH
    from pk.packing import read
    a, b = os.path.join(tmp, 'a'), os.path.join(tmp, 'b')
    write_cfg(a, s * 1.03, sq + np.r_[[[0.015 * s, 0.015 * s, 0]]] * 1)
    args = [HARDEN, 'reharden', '--in', a, '--out', b, '--seed', str(seed), '--hold-move', '1', '--r0', str(r0)]
    if arm == 'sym':
        args += ['--sym-md', '1']
    subprocess.run(args, capture_output=True, text=True, timeout=900)
    if not os.path.exists(b):
        return None
    h = read(b)
    p = read(fq_quench(b, POLISH, tmp=tmp, timeout=900).out)
    return h.s, p.s, (p.sq if p.max_pen() <= 1e-9 else None)


def hop_once(arm, s, sq, sig, seed, tmp):
    from pk.engine import fq_quench, POLISH
    from pk.packing import read
    rng = np.random.default_rng(seed)
    L = 1.03; c = s / 2
    x = sq.copy(); x[:, :2] = c * L + L * (x[:, :2] - c)
    x = sym_kick(x, sig, rng) if arm == 'sym' else x + np.c_[rng.normal(0, sig, (len(x), 2)), rng.normal(0, 20 * sig, len(x))]
    a, b = os.path.join(tmp, 'a'), os.path.join(tmp, 'b')
    write_cfg(a, s * L, x)
    args = [HARDEN, 'reharden', '--in', a, '--out', b, '--seed', str(seed), '--hold-move', '1', '--r0', '0.1']
    if arm == 'sym':
        args += ['--sym-md', '1']
    subprocess.run(args, capture_output=True, text=True, timeout=900)
    if not os.path.exists(b):
        return None
    h = read(b)
    run = fq_quench(b, POLISH, tmp=tmp, timeout=900)
    p = read(run.out) if run.out else None
    if p is None or p.max_pen() > 1e-9:
        return None
    return h.s, h.sq, p.s, p.sq


def run_one(t):
    arm, n, rec, start_id, s0, sq0, hops, rep = t
    sq = np.array(sq0); s = s0
    best_rel, best_sq = None, None
    cur_h, cur_hs = sq, s
    log = []
    t0 = time.process_time()
    with tempfile.TemporaryDirectory() as tmp:
        for h in range(hops):
            sig = SIGS[h % len(SIGS)]
            r = hop_once(arm, cur_hs, cur_h, sig, 1000 * rep + h + 7 * n, tmp)
            if r is None:
                continue
            hs, hsq, ps, psq = r
            if best_rel is None or ps < best_rel - 1e-10:
                best_rel, best_sq = ps, psq
                cur_h, cur_hs = hsq, hs                          # monotone: move to the state that gave the best release
                log.append((h, sig, ps - rec))
    return dict(arm=arm, n=n, rec=rec, start=start_id, s_start=s0, rep=rep, best=best_rel,
                gap=None if best_rel is None else best_rel - rec, hits=bool(best_rel is not None and best_rel < rec + 1e-7),
                log=log, cpu=time.process_time() - t0, sq=None if best_sq is None else np.asarray(best_sq).tolist())


def starts():
    from pk.store import Store, public_sql
    st = Store()
    sc = json.load(open('runs/symsub/symscan.json'))
    out = []
    for n, r in sc.items():
        n = int(n)
        if n < 30 or 'Md' not in r['inv'] or not r['ntilt'] or abs(r['s'] - round(r['s'])) < 1e-9:
            continue
        q = ("SELECT DISTINCT p.id, p.s FROM sighting g JOIN packing p ON p.id = g.packing WHERE p.n = ? AND "
             + public_sql('g') + " AND p.status != 'infeasible' AND p.s > ? AND p.s < ? ORDER BY p.s")
        seen = set()
        for row in st.db.execute(q, (n, r['s'] + 1e-7, r['s'] + 0.1)):
            k = round(row['s'], 7)
            if k in seen:
                continue
            seen.add(k)
            p = st.get(row['id'])
            out.append((n, r['s'], row['id'], p.s, p.sq))
            if len(seen) >= 2:
                break
    return out


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd')
    r = sub.add_parser('run'); r.add_argument('--out', required=True); r.add_argument('--hops', type=int, default=80)
    r.add_argument('--reps', type=int, default=4); r.add_argument('--procs', type=int, default=15)
    sub.add_parser('starts')
    q = sub.add_parser('report'); q.add_argument('out')
    a = ap.parse_args()
    if a.cmd == 'starts':
        for n, rec, pid, s, sq in starts():
            el, d = best_diag_mirror(sq, s)
            sy, ok = symmetrize(to_md(sq, s, el), s)
            print(n, pid, f'gap {s-rec:.4f}', el, f'mismatch {d:.3f}', 'Md-inv' if invariant(sq, s, "Md", 1e-6)[0] else '', 'pairing ok' if ok else 'PAIRING BAD')
        return
    if a.cmd == 'report':
        return report(a.out)
    from multiprocessing import Pool
    from pk import manifest
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'symsub-tt', args=vars(a), procs=a.procs)
    T = []
    for n, rec, pid, s, sq in starts():
        el, d = best_diag_mirror(sq, s)
        sy, ok = symmetrize(to_md(sq, s, el), s)
        if not ok:
            print('skip', n, pid, 'pairing'); continue
        for rep in range(a.reps):
            T.append(('sym', n, rec, pid, s, sy.tolist(), a.hops, rep))
            T.append(('free', n, rec, pid, s, sq.tolist(), a.hops, rep))
    print(len(T), 'runs', flush=True)
    with Pool(a.procs) as pool, open(os.path.join(a.out, 'runs.jsonl'), 'w') as fo:
        for r in pool.imap_unordered(run_one, T):
            fo.write(json.dumps(r) + '\n'); fo.flush()
    man.finish()
    report(a.out)


def report(out):
    R = [json.loads(l) for l in open(os.path.join(out, 'runs.jsonl'))]
    g = collections.defaultdict(list)
    for r in R:
        g[(r['n'], r['start'], r['arm'])].append(r)
    print('n start s_start-rec | free: hits / best gap | sym: hits / best gap')
    for (n, sid) in sorted({(k[0], k[1]) for k in g}):
        f, s_ = g[(n, sid, 'free')], g[(n, sid, 'sym')]
        if not f or not s_:
            continue
        bg = lambda v: min((r['gap'] for r in v if r['gap'] is not None), default=float('nan'))
        print(n, sid, f"{f[0]['s_start'] - f[0]['rec']:+.4f} | {sum(r['hits'] for r in f)}/{len(f)} {bg(f):+.2e} | "
              f"{sum(r['hits'] for r in s_)}/{len(s_)} {bg(s_):+.2e}")
    for arm in ('free', 'sym'):
        v = [r for r in R if r['arm'] == arm]
        print('TOTAL', arm, 'runs', len(v), 'hits', sum(r['hits'] for r in v),
              'starts improved', len({(r['n'], r['start']) for r in v if r['gap'] is not None and r['gap'] < r['s_start'] - r['rec'] - 1e-9}),
              'cpu-h', round(sum(r['cpu'] for r in v) / 3600, 2))


if __name__ == '__main__':
    main()
