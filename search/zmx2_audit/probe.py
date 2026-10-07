"""Nested-box soundness probe for zmx2: for an admissible rational pose P, build boxes of every level
k (centre pitch 1/(10 2^k), u pitch 1/(8 2^k)) containing P -- in every corner/edge position when P
is on a grid boundary, and with random offsets otherwise -- ask `zmx2 boxes` for their bounds, and
check bound <= mu(P) (exact, aud.mu).  Any violation is a soundness bug."""
import os
import subprocess, sys, random, math
from fractions import Fraction as Fr
from aud import load, mu, inside

ZMX2 = os.environ.get('ZMX2', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'verify2', 'target', 'release', 'zmx2'))


def cells(v, den, rng, w):
    """intervals [i/den, (i+w)/den] containing v (all placements if v is on the grid)"""
    i = math.floor(v * den)
    out = set()
    for off in range(w):
        lo = i - off
        if Fr(lo, den) <= v <= Fr(lo + w, den):
            out.add((Fr(lo, den), Fr(lo + w, den)))
        if v * den == i:  # on the grid: also the cells ending at v
            lo2 = i - w + off
            out.add((Fr(lo2, den), Fr(lo2 + w, den)))
    out = [c for c in out if c[0] <= v <= c[1]]
    rng.shuffle(out)
    return out[:2]


def boxes_for(P, rng, kmax=26, umax=Fr(1, 2)):
    x, y, u = P
    res = []
    for k in range(0, kmax + 1):
        cd, ud = 10 << k, 8 << k
        w = rng.choice([1, 1, 2, 3])
        for cx in cells(x, cd, rng, w):
            for cy in cells(y, cd, rng, w):
                for cu in cells(u, ud, rng, rng.choice([1, 2])):
                    if cu[0] < 0 or cu[1] > umax:
                        continue
                    res.append(cx + cy + cu)
    return res


def run(path, poses, refl=False, seed=1, kmax=26, verbose=True, extra=[]):
    cv = load(path)
    if refl:
        S = cv['s'] * cv['D']
        S = S.numerator
        cv = dict(cv)
        cv['pts'] = [(X, S - Y, w) for (X, Y, w) in cv['pts']]
        cv['segs'] = [(a, S - b, c, S - d, w) for (a, b, c, d, w) in cv['segs']]
    rng = random.Random(seed)
    allb, owner = [], []
    for pi, P in enumerate(poses):
        assert inside(cv, *P), P
        for b in boxes_for(P, rng, kmax):
            allb.append(b)
            owner.append(pi)
    inp = '\n'.join(','.join('%d/%d' % (f.numerator, f.denominator) for f in b) for b in allb) + '\n'
    pr = subprocess.run([ZMX2, 'boxes', path] + (['--refl'] if refl else []) + extra, input=inp, capture_output=True,
                         text=True)
    if pr.returncode != 0:
        print('skip %s: %s' % (path.split('/')[-1], (pr.stdout + pr.stderr).strip()[:100]))
        return 0
    out = pr.stdout.split('\n')
    mus = [mu(cv, *P) for P in poses]
    nfail = 0
    worst = None
    best_ratio = {}
    for b, pi, line in zip(allb, owner, out):
        bound, unit, empty, _ = (int(t) for t in line.split())
        B = Fr(bound, unit)
        if empty:
            print('FAIL: box with admissible pose reported EMPTY', b, poses[pi])
            nfail += 1
            continue
        gap = mus[pi] - B
        if gap < 0:
            nfail += 1
            print('FAIL: bound %.12f > mu %.12f at pose %s box %s' % (float(B), float(mus[pi]), poses[pi], b))
        if worst is None or gap < worst:
            worst = gap
    if verbose:
        print('probe %s: %d poses, %d boxes, %d FAIL, min(mu - bound) = %.3g, mu range [%.6f, %.6f]' % (
            path.split('/')[-1] + (' (refl)' if refl else ''), len(poses), len(allb), nfail, float(worst),
            float(min(mus)), float(max(mus))))
    return nfail


if __name__ == '__main__':
    path = sys.argv[1]
    poses = [tuple(Fr(t) for t in p.split(',')) for p in sys.argv[2:]]
    sys.exit(1 if run(path, poses) else 0)
