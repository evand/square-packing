"""Example fill (interface check, not a contender): the clean upright cells (meta['cells']), topped up to m with the
least-penetrating of 20 random poses at random uncovered points; extras dropped.  One candidate."""
import numpy as np
import harness as H


def fill(S, sq, meta, m, rng):
    sq = np.asarray(sq, float)
    add = list(meta['cells'])[:m]
    xs, ys, free = H.free_mask(S, np.vstack([sq] + ([np.array(add)] if add else [])), h=0.1)
    pts = np.argwhere(free)
    while len(add) < m and len(pts):
        cur = np.vstack([sq] + ([np.array(add)] if add else []))
        best = None
        for _ in range(20):
            j, i = pts[rng.randrange(len(pts))]
            pose = (float(xs[i]), float(ys[j]), rng.uniform(0, 90))
            p = H.pen_one(S, cur, pose)
            if best is None or p < best[0]:
                best = (p, pose)
        add.append(best[1])
    return [np.vstack([sq, np.array(add)])] if add else [sq]
