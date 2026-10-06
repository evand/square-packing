"""Many restarts in chosen classes, keeping every optimum (no early stop), plus basin hopping from the best.

  python3 intensify.py NWORK RESTARTS_PER_TASK TASKS_PER_CLASS h class [class ...]
      class written as 9 digits k0k1k2k3 e0e1e2e3 z (e.g. 111111114)
Appends to runs/intensify.jsonl one line per optimum (margin, poses, max tilt, kind).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys, json, zlib
import numpy as np
from multiprocessing import Pool
from opt import Problem
from occ import regions, assignment, label
from census import random_seed, tiling_seed, pinwheel_seed, library_seed, load_library, RUNS


def task(args):
    h, c, R, t, lib = args
    rng = np.random.default_rng(zlib.crc32(f"int{h}{c}{t}".encode()))
    boxes = [regions(h)[r] for r in assignment(c)]
    pr = Problem(boxes)
    out, best = [], (-np.inf, None)
    for r in range(R):
        kind = r % 6
        if kind == 5 and best[1] is not None:          # hop: perturb best
            P0 = best[1].copy()
            k = rng.integers(1, 4)
            idx = rng.choice(12, k, replace=False)
            P0[idx] = random_seed([boxes[i] for i in idx], rng)
            P0[:, :2] += rng.normal(0, 0.03, (12, 2))
        elif kind == 0:
            P0 = random_seed(boxes, rng)
        elif kind == 1:
            P0 = tiling_seed(boxes, rng)
        elif kind == 2:
            P0 = pinwheel_seed(boxes, rng, r % 2)
        elif kind == 3 and lib:
            P0 = library_seed(lib, boxes, rng)
        else:
            P0 = random_seed(boxes, rng); P0[:, 2] = rng.normal(0, 0.15, 12) % (np.pi / 2)
        m, P = pr.solve(P0)
        if m > best[0]:
            best = (m, P)
        tilt = np.degrees(np.minimum(P[:, 2], np.pi / 2 - P[:, 2]))
        out.append({"h": h, "class": list(c), "kind": kind, "margin": m, "maxtilt": float(tilt.max()),
                    "ntilted": int((tilt > 0.01).sum()), "poses": P.tolist()})
    return out


if __name__ == "__main__":
    nw, R, T, h = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
    cls = [tuple(int(ch) for ch in s) for s in sys.argv[5:]]
    lib = load_library()
    jobs = [(h, c, R, t, lib) for c in cls for t in range(T)]
    with Pool(nw) as p, open(os.path.join(RUNS, "intensify.jsonl"), "a") as f:
        for res in p.imap_unordered(task, jobs):
            for r in res:
                f.write(json.dumps(r) + "\n")
            f.flush()
            b = max(res, key=lambda r: r["margin"])
            print(label(tuple(b["class"])), f"h={h} best {b['margin']:+.6f} of {len(res)}", flush=True)
