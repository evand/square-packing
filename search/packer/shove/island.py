#!/usr/bin/env python3
"""Island shove (IS1, 10-17; Evan: "shove just the grid squares, or also shove the central non-grid section").

The corner shove (harness.shove) moves a corner's grid block with its corner and everything else with the opposite
corner, so the strip opens along the block's staircase: at the 171 / 198 family the tilted islands never reach it
(CS2, L198: every return is an exact-side re-routing of the upright frame).  Here each tilted island (a connected
cluster of tilted squares) takes a fraction t of the displacement instead:

    t = 0    the island moves with the corner block (only the grid squares on the far side move: the gap opens between
             the island and the rest of the jumble)
    t = 1/2  the island moves halfway (slack on both of its faces)
    t = 1    the island moves with the rest (= the plain corner shove)

Per square: offset = block offset + t * (rest offset - block offset), t = 0 for the block, 1 for other non-island
squares.  meta carries 'offsets' (n x 2) so fill_resolve.unshove (and the wallL control) can undo it.

  python3 shove/island.py scan --src FILE                       # islands of a packing
  python3 shove/island.py build --out runs/shove/is1 --from-bench runs/shove/l198 [--levels ...] [--fracs 0,0.5,1]
"""
from __future__ import annotations
import argparse, itertools, json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
sys.path.insert(0, PK); sys.path.insert(0, HERE)
import harness as H

TILT_MIN = 0.5        # degrees: a square counts as tilted above this
LINK = 1.6            # centre distance that joins two tilted squares into one island


def tilt(a):
    a = np.asarray(a, float) % 90.0
    return np.minimum(a, 90.0 - a)


def islands(s, sq, tilt_min=TILT_MIN, link=LINK, min_size=3):
    """list of index arrays, largest first: connected clusters of tilted squares (size >= min_size)"""
    sq = np.asarray(sq, float)
    T = np.nonzero(tilt(sq[:, 2]) > tilt_min)[0]
    if len(T) == 0:
        return []
    P = sq[T, :2]
    D = np.hypot(P[:, None, 0] - P[None, :, 0], P[:, None, 1] - P[None, :, 1])
    lab = -np.ones(len(T), int); c = 0
    for i in range(len(T)):
        if lab[i] >= 0:
            continue
        st = [i]; lab[i] = c
        while st:
            j = st.pop()
            for k in np.nonzero((D[j] < link) & (lab < 0))[0]:
                lab[k] = c; st.append(k)
        c += 1
    out = [T[lab == k] for k in range(c)]
    return sorted([o for o in out if len(o) >= min_size], key=len, reverse=True)


def ishove(s, sq, corner, d, frac, block=None, isl=None, fix=None):
    """frac: list of t per island (islands(s, sq) order unless isl is given).  -> (S, sq2, meta)"""
    sq = np.asarray(sq, float).copy()
    if block is None:
        block = H.corner_blocks(s, sq)[corner]
    if isl is None:
        isl = islands(s, sq)
    bx = d if corner[1] == 'R' else 0.0
    by = d if corner[0] == 'U' else 0.0
    ox, oy = d - bx, d - by
    t = np.ones(len(sq))
    for I, f in zip(isl, frac):
        t[np.asarray(I, int)] = float(f)
    t[list(block)] = 0.0
    nfix = 0
    if (H.shove_fix_on() if fix is None else fix):          # 10-18 collision fix (harness.fix_levels)
        t, nfix = H.fix_levels(sq, t, (bx, by), (ox, oy), groups=isl)
        block = np.nonzero(t < 1e-12)[0]
    off = np.c_[bx + t * (ox - bx), by + t * (oy - by)]
    sq[:, :2] += off
    S = s + d
    blk = np.zeros(len(sq), bool); blk[list(block)] = True
    meta = dict(corner=corner, d=d, block=sorted(int(i) for i in block), s0=float(s), n0=len(sq),
                k0=int(math.floor(s + 1e-9)), strip=H.strip_rects(S, sq, blk, corner, d), offsets=off.tolist(),
                frac=list(map(float, frac)), islands=[sorted(map(int, I)) for I in isl], fixed=nfix, t=t.tolist())
    meta['cells'] = H.clean_cells(S, sq, meta)
    return S, sq, meta


def from_task(s, sq, task, d):
    """bench hook: task['ishove'] = {'frac': [...]} (islands recomputed on the source, same order)"""
    return ishove(s, sq, task['corner'], d, task['ishove']['frac'])


def scan(a):
    import mcmin
    s, sq = mcmin.load_deg(a.src); sq = np.asarray(sq, float)
    for I in islands(s, sq):
        P = sq[I]
        print(len(I), 'squares, centre (%.2f, %.2f)' % tuple(P[:, :2].mean(0)), 'angles', sorted(set(np.round(P[:, 2] % 90, 2))))


def build(a):
    """tasks: every (source task of --from-bench, island fraction vector except all-1) ; sources filtered by --ns"""
    T = json.load(open(os.path.join(a.from_bench, 'tasks.json')))
    ns = set(map(int, a.ns.split(','))) if a.ns else None
    fr = [float(x) for x in a.fracs.split(',')]
    import mcmin, shutil
    from pk import eventbench as EB
    os.makedirs(os.path.join(a.out, 'starts'), exist_ok=True)
    out, seen = [], {}
    for t in T['tasks']:
        if ns and t['n'] not in ns:
            continue
        src = os.path.join(a.from_bench, t['src'])
        dst = os.path.join(a.out, 'starts', os.path.basename(t['src']))
        if not os.path.exists(dst):
            shutil.copy(src, dst)
        if src not in seen:
            s, sq = mcmin.load_deg(src); seen[src] = islands(s, np.asarray(sq, float))
        k = min(len(seen[src]), a.max_islands)
        if k == 0:
            continue
        for v in itertools.product(fr, repeat=k):
            if all(abs(x - 1) < 1e-12 for x in v):
                continue
            if a.fracs_only and ','.join(f'{x:g}' for x in v) not in a.fracs_only.split(';'):
                continue
            frac = list(v) + [1.0] * (len(seen[src]) - k)
            nt = dict(t, src='starts/' + os.path.basename(t['src']), ishove=dict(frac=frac))
            if a.ms_shift or a.ms_plus:                     # 10-18: d = 2 tasks (counts shifted by ~ one more L)
                base = t.get('ms') or [2 * t['k'], 2 * t['k'] + 1, 2 * t['k'] + 2]
                ms = sorted(set([m + a.ms_shift for m in base] + [max(base) + a.ms_shift + j for j in range(1, a.ms_plus + 1)]))
                nt['ms'] = ms
                nt['targets'] = {str(m): dict(n2=t['n'] + m, fT=EB.frontier_now(t['n'] + m), fNow=EB.frontier_now(t['n'] + m))
                                 for m in ms}
            nt['tags'] = dict(t.get('tags') or {}, frac=','.join(f'{x:g}' for x in frac), nisl=len(seen[src]))
            nt['id'] = len(out)
            out.append(nt)
    json.dump(dict(cutoff=T.get('cutoff', 'now'), t=T.get('t'), tasks=out), open(os.path.join(a.out, 'tasks.json'), 'w'),
              indent=0, default=str)
    print(len(out), 'tasks from', len(seen), 'sources')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    b = sp.add_parser('scan'); b.add_argument('--src', required=True)
    b = sp.add_parser('build'); b.add_argument('--out', required=True); b.add_argument('--from-bench', required=True)
    b.add_argument('--ns', help='comma list of source n to keep'); b.add_argument('--fracs', default='0,0.5,1')
    b.add_argument('--max-islands', type=int, default=2)
    b.add_argument('--ms-shift', type=int, default=0); b.add_argument('--ms-plus', type=int, default=0)
    b.add_argument('--fracs-only', help='keep only these fraction vectors, ";"-separated, e.g. "0.5,0.5;0,0"')
    a = ap.parse_args()
    {'scan': scan, 'build': build}[a.cmd](a)
