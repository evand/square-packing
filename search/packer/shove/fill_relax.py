"""Fill B (10-16): infeasible fill + expand / relax (Evan: "fill with infeasible packings (least-bad among random
options?) and then expand/relax from there").

1. Least-bad infeasible placement.  The m new squares go in one at a time; for each, draw K random poses and keep
   the one with the least TOTAL penetration (sum of SAT depths against every square already present, old or new,
   plus wall excess).  K is the "best of k random options" knob (Evan, 10-16): K = 1 is pure random, large K is near
   greedy.  Positions: uncovered raster points, mostly in / near the opened strip (meta['strip'] dilated by 0.6),
   some anywhere uncovered; angles: 0, a nearby square's tilt (+ small jitter) or uniform.
2. Expand / relax.  fq quench from the overlapping state with a pre-stage: rounded squares (round_r0, compress with
   the side free, then r -> 0 at the held box), or a held box slightly above S with an SLP polish (no ALM), or plain
   ALM at loosens 1.1 / 1.2 (fq's loosen dilates positions and box before squeezing).  Dev tests (NOTES_relax.md):
   plain ALM at any loosen, --shrink, soft squares and one-at-a-time insertion with quenches between all collapse to
   the k + 1 grid; the rounded pre-stages are the only ones that keep non-grid structure.
3. Return the relaxed feasible states (best internal side first); the bench's evaluate quenches them again.

Knobs (env vars, read per call):
  SHOVE_K        best-of-K poses per new square (default 16)
  SHOVE_RELAX    comma list of relax schedules (keys of RELAX; default 'round.15c,round.25c,holdP1.01')
  SHOVE_CELLS    1 = put the clean upright cells (meta['cells']) in first, least-bad poses for the rest (default 1)
  SHOVE_NPLACE   independent least-bad placements per call, each relaxed by every schedule (default 1)
  SHOVE_RAW      1 = also return the raw infeasible state (no relax) as a candidate (default 0)
  SHOVE_LOG      path: append one json line per call with the internal sides (for notes)
"""
import json, math, os, tempfile
import numpy as np
import harness as H


def _k():
    return int(os.environ.get('SHOVE_K', '64'))


# 10-20: penetration + least-bad placement moved to pk/fill.py (placer 'leastbad'); same draws and output.
from pk.fill import pen_vec, total_pen, _sample_pts, _angle, leastbad_place  # noqa: F401


def place(S, X, strip, m, rng, K, resample_every=4):
    """add m squares to X one at a time, each the least-total-penetration of K random poses -> (X', total pen)"""
    X2, tot, _ = leastbad_place(S, X, strip, m, rng, K, resample_every)
    return X2, tot


# relax schedules: pk/schedules.py (10-20; this file's RELAX table moved there, same QOpts)
from pk.schedules import RELAX_NAMES, SCHEDULES
RELAX = {k: SCHEDULES[k] for k in RELAX_NAMES}


def _quench(S, X, tmp, kw):
    import hop
    r = hop.quench(S, [tuple(map(float, q)) for q in X], tmp, opt=kw)
    return None if r is None else (float(r[0]), np.asarray(r[1], float))


def fill(S, sq, meta, m, rng):
    """least-bad placement (clean cells first if SHOVE_CELLS=1), then each relax schedule in SHOVE_RELAX"""
    K = _k()
    names = [v for v in os.environ.get('SHOVE_RELAX', 'round.15c,round.25c,holdP1.01').split(',') if v]
    cells = os.environ.get('SHOVE_CELLS', '1') == '1'
    raw = os.environ.get('SHOVE_RAW', '0') == '1'
    sq = np.asarray(sq, float)
    strip = meta['strip']
    base = sq
    if cells and meta['cells']:
        base = np.vstack([sq, np.array(meta['cells'][:m])])
    nplace = int(os.environ.get('SHOVE_NPLACE', '1'))
    log = dict(n0=meta['n0'], s0=meta['s0'], corner=meta['corner'], d=meta['d'], m=m, K=K, cells=cells,
               ncells=len(base) - len(sq), pen=[])
    out = []
    with tempfile.TemporaryDirectory(prefix='relax_') as tmp:
        for ip in range(nplace):
            X0, pen0 = place(S, base, strip, m - (len(base) - len(sq)), rng, K)
            log['pen'].append(round(pen0, 4))
            if raw:
                out.append((S + 1e3, 'raw', X0))
            for nm in names:
                r = _quench(S, X0, tmp, RELAX[nm])
                log[f'{nm}#{ip}' if nplace > 1 else nm] = None if r is None else r[0]
                if r:
                    out.append((r[0], f'{nm}#{ip}', r[1]))
    out.sort(key=lambda t: t[0])
    log['order'] = [t[1] for t in out]
    if os.environ.get('SHOVE_LOG'):
        with open(os.environ['SHOVE_LOG'], 'a') as f:
            f.write(json.dumps(log) + '\n')
    return [t[2] for t in out]
