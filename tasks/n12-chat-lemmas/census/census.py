"""Occupancy census at side 4.

  python3 census.py free  NJOBS RESTARTS      # stage A: no class constraint; library of near-0 optima
  python3 census.py class NWORK               # stage B: every (h, class), adaptive restarts

Outputs under runs/: free.jsonl (stage A), classes.jsonl (stage B, one line per (h, class)),
witnesses/<tag>.json (best poses per (h, class)).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys, json, time, itertools, zlib
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linear_sum_assignment
from opt import Problem, sat_margin, pair_sat, wall_gap, L
from occ import enumerate_classes, regions, assignment, classify, label, strips, canon

HS = [1.1213, 1.15, 1.18, 1.2071]
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
WIT = os.path.join(RUNS, "witnesses")
ZERO = -1e-9

PINWHEEL = None   # 26 deg pinwheel of 5 (filled lazily from runs/pin5.json if present)


def d4(P, g):
    """Apply D4 element g in 0..7 to poses (theta taken mod pi/2)."""
    P = P.copy()
    if g & 4:
        P[:, 0] = L - P[:, 0]; P[:, 2] = -P[:, 2]
    for _ in range(g & 3):
        x, y = P[:, 0].copy(), P[:, 1].copy()
        P[:, 0], P[:, 1] = L - y, x
    P[:, 2] = np.mod(P[:, 2], np.pi / 2)
    return P


def fit_to_class(P, boxes):
    """Assign the poses of P to the box slots (Hungarian on distance to box) and clip."""
    b = np.array(boxes)
    cost = np.zeros((len(P), len(b)))
    for i, (x, y, _) in enumerate(P):
        dx = np.maximum(0, np.maximum(b[:, 0] - x, x - b[:, 1]))
        dy = np.maximum(0, np.maximum(b[:, 2] - y, y - b[:, 3]))
        cost[i] = dx ** 2 + dy ** 2
    r, c = linear_sum_assignment(cost)
    Q = np.zeros((len(b), 3))
    Q[c] = P[r]
    return Q


GRID = [(i, j) for i in range(4) for j in range(4)]


def tiling_seed(boxes, rng):
    """A 12-subset of the 4x4 tiling grid fitted to the boxes, jittered."""
    sub = rng.choice(16, 12, replace=False)
    P = np.array([[0.5 + GRID[s][0], 0.5 + GRID[s][1], 0.0] for s in sub])
    P = fit_to_class(P, boxes)
    P[:, :2] += rng.normal(0, 0.03, (len(P), 2))
    P[:, 2] = np.mod(rng.normal(0, 0.05, len(P)), np.pi / 2)
    return P


def random_seed(boxes, rng):
    b = np.array(boxes)
    return np.column_stack([rng.uniform(b[:, 0], b[:, 1]), rng.uniform(b[:, 2], b[:, 3]),
                            rng.uniform(0, np.pi / 2, len(b))])


def pinwheel_seed(boxes, rng, kind):
    """5 central squares as a 45 deg quincunx or a ~26 deg pinwheel; others random in their boxes."""
    P = random_seed(boxes, rng)
    zi = [i for i, bb in enumerate(boxes) if bb[0] > 0.5 + 1e-9 and bb[1] < L - 0.5 - 1e-9 and bb[2] > 0.5 + 1e-9 and bb[3] < L - 0.5 - 1e-9]
    if not zi:
        zi = list(range(5))
    c = L / 2
    if kind == 0:
        a = np.pi / 4; s = 0.78 + rng.uniform(0, 0.05)
        pts = [(c, c)] + [(c + s * dx, c + s * dy) for dx, dy in ((1, 1), (1, -1), (-1, 1), (-1, -1))]
    else:
        a = np.radians(26.57 + rng.normal(0, 2)); s = 1.0
        # pinwheel: centre square + 4 around at offsets rotated by a, spacing ~1/cos-ish
        r = 1.0 + rng.uniform(0, 0.1)
        pts = [(c, c)] + [(c + r * np.cos(a + q * np.pi / 2 + np.arctan(0.5)), c + r * np.sin(a + q * np.pi / 2 + np.arctan(0.5))) for q in range(4)]
    rng.shuffle(pts)
    for i, p in zip(zi, pts):
        P[i] = (p[0], p[1], a % (np.pi / 2))
    # also non-Z squares near walls: axis-parallel, pushed toward their walls
    for i, bb in enumerate(boxes):
        if i in zi:
            continue
        P[i, 2] = rng.normal(0, 0.05) % (np.pi / 2)
    return P


def load_library():
    lib = []
    f = os.path.join(RUNS, "free.jsonl")
    if os.path.exists(f):
        for line in open(f):
            r = json.loads(line)
            if r["margin"] >= -0.02:
                lib.append(np.array(r["poses"]))
    return lib


def library_seed(lib, boxes, rng):
    P = lib[rng.integers(len(lib))]
    P = d4(P, int(rng.integers(8)))
    P = fit_to_class(P, boxes)
    P[:, :2] += rng.normal(0, 0.02, (len(P), 2))
    return P


def run_free(args):
    seed, R = args
    rng = np.random.default_rng(10_000 + seed)
    pr = Problem([(0.5, L - 0.5, 0.5, L - 0.5)] * 12)
    out = []
    for r in range(R):
        kind = r % 3
        if kind == 0:
            P0 = random_seed(pr.boxes, rng)
        elif kind == 1:
            P0 = pinwheel_seed(list(map(tuple, pr.boxes)), rng, r % 2)
            # in the free problem every box is the full square; put 5 at the centre, others random
        else:
            P0 = tiling_seed(list(map(tuple, pr.boxes)), rng)
            P0[:, 2] = rng.uniform(0, np.pi / 2, 12) * (rng.random(12) < 0.3)
        m, P = pr.solve(P0)
        out.append({"seed": seed, "r": r, "kind": kind, "margin": m, "poses": P.tolist(),
                    "classes": {str(h): list(classify(P, h)) for h in HS}})
    return out


def run_class(args):
    h, c, rmax, lib = args
    t0 = time.time()
    tag = f"h{h}_" + "".join(map(str, c))
    rng = np.random.default_rng(zlib.crc32(tag.encode()))
    reg = regions(h)
    boxes = [reg[r] for r in assignment(c)]
    pr = Problem(boxes)
    best = (-np.inf, None, None)
    hist = []
    nr = 0
    for r in range(rmax):
        kind = r % 5
        if kind == 0 or (kind == 4 and not lib):
            P0 = random_seed(boxes, rng)
        elif kind == 1:
            P0 = tiling_seed(boxes, rng)
        elif kind == 2:
            P0 = pinwheel_seed(boxes, rng, r % 2)
        elif kind == 3:
            P0 = random_seed(boxes, rng); P0[:, 2] = rng.normal(0, 0.1, 12) % (np.pi / 2)
        else:
            P0 = library_seed(lib, boxes, rng)
        m, P = pr.solve(P0)
        nr += 1
        hist.append(round(m, 6))
        if m > best[0]:
            best = (m, P, kind)
        if best[0] >= ZERO and nr >= 6:
            break
    m, P, kind = best
    # independent recheck
    G = pair_sat(P)
    ok_class = tuple(classify(P, h)) == tuple(c) or canon(classify(P, h)) == canon(c)
    os.makedirs(WIT, exist_ok=True)
    with open(os.path.join(WIT, tag + ".json"), "w") as f:
        json.dump({"h": h, "class": list(c), "label": label(c), "margin": m, "boxes": boxes,
                   "poses_x_y_thetarad": P.tolist(), "tilts_deg": np.degrees(P[:, 2]).tolist(),
                   "pair_min": float(G[np.isfinite(G)].min()), "wall_min": float(wall_gap(P).min()),
                   "class_recheck": ok_class}, f, indent=1)
    tilt = np.degrees(np.minimum(P[:, 2], np.pi / 2 - P[:, 2]))
    return {"h": h, "class": list(c), "label": label(c), "k": sum(c[:4]), "e": sum(c[4:8]), "z": c[8],
            "strips": list(strips(c)), "margin": m, "restarts": nr, "best_kind": kind, "maxtilt": float(tilt.max()),
            "hist": hist, "secs": time.time() - t0, "witness": f"runs/witnesses/{tag}.json"}


if __name__ == "__main__":
    os.makedirs(RUNS, exist_ok=True)
    mode = sys.argv[1]
    if mode == "free":
        nj, R = int(sys.argv[2]), int(sys.argv[3])
        with Pool(nj) as p, open(os.path.join(RUNS, "free.jsonl"), "a") as f:
            for res in p.imap_unordered(run_free, [(s, R) for s in range(nj * 4)]):
                for r in res:
                    f.write(json.dumps(r) + "\n")
                f.flush()
    elif mode == "class":
        nw = int(sys.argv[2])
        hs = [float(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else HS
        lib = load_library()
        done = set()
        fo = os.path.join(RUNS, "classes.jsonl")
        if os.path.exists(fo):
            for line in open(fo):
                r = json.loads(line); done.add((r["h"], tuple(r["class"])))
        jobs = []
        for h in hs:
            for c in enumerate_classes():
                if (h, c) in done:
                    continue
                rmax = 200 if c[8] <= 7 else 60
                jobs.append((h, c, rmax, lib))
        # hardest (z small, many restarts) first for load balance
        jobs.sort(key=lambda j: j[1][8])
        with Pool(nw) as p, open(fo, "a") as f:
            for res in p.imap_unordered(run_class, jobs):
                f.write(json.dumps(res) + "\n"); f.flush()
                print(f"h={res['h']} {res['label']:40s} margin={res['margin']:+.6f} R={res['restarts']} {res['secs']:.0f}s", flush=True)
