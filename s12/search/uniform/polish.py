#!/usr/bin/env python3
"""Polish a uniform certificate: move points to maximise the critical container.

The point multiset is kept D4-symmetric (orbit representatives are moved, images follow; a
representative on a symmetry axis stays on it).  The objective is the critical denominator
D_crit -- the smallest D at which the exact verifier (N=2000) still passes with the integer
coordinates fixed, i.e. the largest container K/D -- minimised lexicographically with the
number of failing angle bins one step beyond it.  Random small moves, accepted when the score
does not get worse (simulated-annealing style hill climb; every accepted state is a verified
certificate).

    python3 search/uniform/polish.py CERT OUT [--iters 3000] [--Dw 4000] [--N 2000]
        [--threads 4] [--seed 0] [--time 1800]
"""
import argparse, math, os, random, subprocess, sys, tempfile, time
from fractions import Fraction
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
VERIFY = REPO / "verify" / "target" / "release" / "verify"


def read_cert(path):
    t = Path(path).read_text().split()
    sn, sd, D, W, m = map(int, t[:5])
    rows = [(int(t[5+3*i]), int(t[6+3*i]), int(t[7+3*i])) for i in range(m)]
    return Fraction(sn, sd), D, W, rows


SYM = 'd4'


def images(x, y, K):
    if SYM == 'none': return {(x, y)}
    return {(x, y), (K-x, y), (x, K-y), (K-x, K-y), (y, x), (K-y, x), (y, K-x), (K-y, K-x)}


class Polisher:
    def __init__(self, cert, Dw, N, threads, log):
        s, D, W, rows = read_cert(cert)
        self.k = W; self.N = N; self.threads = threads; self.log = log
        K0 = s * D
        if K0.denominator == 1 and D >= 1000:
            # exact: keep the certificate's own grid (doubled if the container is odd in grid units)
            f = max(1, -(-Dw // D))
            if (K0.numerator * f) % 2: f *= 2
            self.Dw = D * f; self.K = K0.numerator * f
            pts = [(x * f, y * f) for x, y, w in rows for _ in range(w)]
        else:
            self.Dw = Dw; self.K = 2 * int(round(float(s) * Dw / 2))
            pts = [(int(round(x * Dw / D)), int(round(y * Dw / D))) for x, y, w in rows for _ in range(w)]
        # group into D4 orbits
        self.orbits = []   # list of representative (x,y)
        seen = set()
        for p in pts:
            if p in seen: continue
            img = images(p[0], p[1], self.K)
            seen |= img
            self.orbits.append(min(img))
        m = sum(len(images(x, y, self.K)) for x, y in self.orbits)
        if m != len(pts): log(f"  WARNING: input not D4-symmetric ({len(pts)} points, orbits give {m}); symmetrised")
        self.m = m
        self.td = tempfile.mkdtemp(); self.tmp = os.path.join(self.td, "c.txt"); self.wit = os.path.join(self.td, "w.txt")
        self.calls = 0

    def points(self, orbits):
        P = set()
        for x, y in orbits: P |= images(x, y, self.K)
        return sorted(P)

    def write(self, orbits, D, path):
        P = self.points(orbits)
        with open(path, "w") as f:
            f.write(f"{self.K} {D}\n{D}\n{self.k}\n{len(P)}\n")
            for x, y in P: f.write(f"{x} {y} 1\n")

    def check(self, orbits, D):
        """returns (verified, nfail) at denominator D"""
        self.write(orbits, D, self.tmp); self.calls += 1
        r = subprocess.run([str(VERIFY), self.tmp, "12", str(self.N), str(self.threads), "1", self.wit],
                           capture_output=True, text=True)
        ok = r.returncode == 0 and "VERIFIED" in r.stdout and "NOT VERIFIED" not in r.stdout
        nf = 0
        for l in r.stdout.split("\n"):
            if l.startswith("wrote"): nf = int(l.split()[1])
        return ok, nf

    def score(self, orbits, Dhint):
        """(D_crit, nfail at D_crit-1); Dhint = a D known/expected to pass"""
        D = Dhint
        ok, nf = self.check(orbits, D)
        if not ok:
            step = 1
            while not ok:
                D += step; step *= 2; ok, nf = self.check(orbits, D)
                if D > 4 * self.Dw: return None
            lo = D - step // 2; hi = D
        else:
            step = 1
            while ok:
                lo = D; D -= step; step *= 2; ok, nf = self.check(orbits, D)
                if D < self.Dw // 2: return None
            hi = lo; lo = D
        # invariant: hi passes, lo fails
        while hi - lo > 1:
            mid = (lo + hi) // 2
            ok, nf = self.check(orbits, mid)
            if ok: hi = mid
            else: lo = mid
        _, nf = self.check(orbits, hi - 1)
        return hi, nf

    def worst_witness(self, orbits, D):
        """centre (grid units at denominator D) of the worst placement at denominator D"""
        self.check(orbits, D)
        for line in open(self.wit):
            q = line.split()
            if len(q) == 4:
                return float(q[2]) * D, float(q[3]) * D
        return None

    def fill(self, orbits, mmax, sc, log):
        """add orbits (size 8, then 4, then the centre) at the worst placements until m = mmax"""
        o = list(orbits); K = self.K
        while True:
            m = sum(len(images(x, y, K)) for x, y in o)
            room = mmax - m
            if room < 0:
                # too many points: drop the orbit whose removal costs least
                best = None
                for i, (x, y) in enumerate(o):
                    if len(images(x, y, K)) > -room: continue
                    o2 = o[:i] + o[i+1:]
                    sc2 = self.score(o2, sc[0])
                    if sc2 is not None and (best is None or sc2 < best[0]): best = (sc2, i)
                if best is None: break
                sc, i = best; x, y = o.pop(i)
                log(f"  drop: -{len(images(x, y, K))} at {(x, y)} -> m={m - len(images(x, y, K))} D_crit={sc[0]} s={K/sc[0]:.6f} nfail={sc[1]}")
                continue
            if room <= 0: break
            w = self.worst_witness(o, sc[0] - 1)
            if w is None: break
            x, y = int(round(w[0])), int(round(w[1]))
            x = min(max(x, 0), K); y = min(max(y, 0), K)
            if SYM == 'none':
                cand = (x, y)
            elif room >= 8:
                if x == y: x += 1
                if 2 * x == K: x += 1
                if 2 * y == K: y += 1
                cand = (x, y)
            elif room >= 4:
                # project onto the nearest symmetry axis (size-4 orbit)
                opts = [(x, K // 2), (K // 2, y), ((x + y) // 2, (x + y) // 2), ((x + K - y) // 2, (K - x + y) // 2)]
                opts = [p for p in opts if len(images(p[0], p[1], K)) == 4]
                cand = min(opts, key=lambda p: abs(p[0] - x) + abs(p[1] - y))
            elif room >= 1 and (K // 2, K // 2) not in self.points(o):
                cand = (K // 2, K // 2)
            else:
                break
            rep = min(images(cand[0], cand[1], K))
            if rep in o: break
            o.append(rep)
            sc2 = self.score(o, sc[0])
            if sc2 is None: o.pop(); break
            sc = sc2
            log(f"  fill: +{len(images(*cand, K))} at {cand} -> m={m + len(images(*cand, K))} D_crit={sc[0]} s={K/sc[0]:.6f} nfail={sc[1]}")
        self.m = sum(len(images(x, y, K)) for x, y in o)
        return o, sc

    def kick(self, orbits, rng, D):
        """relocate one orbit (same size) to the worst placement of the remaining set at denominator D"""
        o = list(orbits); K = self.K
        i = rng.randrange(len(o)); x, y = o.pop(i)
        size = len(images(x, y, K))
        w = self.worst_witness(o, D)
        if w is None: return None
        x, y = int(round(w[0])), int(round(w[1]))
        x = min(max(x, 0), K); y = min(max(y, 0), K)
        if SYM == 'none':
            cand = (x, y)
        elif size == 8:
            if x == y: x += 1
            if 2 * x == K: x += 1
            if 2 * y == K: y += 1
            cand = (x, y)
        elif size == 4:
            opts = [(x, K // 2), (K // 2, y), ((x + y) // 2, (x + y) // 2), ((x + K - y) // 2, (K - x + y) // 2)]
            opts = [p for p in opts if len(images(p[0], p[1], K)) == 4]
            cand = min(opts, key=lambda p: abs(p[0] - x) + abs(p[1] - y))
        else:
            return None
        rep = min(images(cand[0], cand[1], K))
        if rep in o or len(images(*cand, K)) != size: return None
        o.append(rep)
        return o

    def move(self, orbits, rng, big):
        o = list(orbits); i = rng.randrange(len(o)); x, y = o[i]
        K = self.K
        mag = rng.choice([1, 2, 4, 8, 16, 32] + ([64, 128] if big else []))
        dx = rng.randint(-mag, mag); dy = rng.randint(-mag, mag)
        if SYM == 'd4':
            if x == y:                     # on the diagonal: stay on it
                dy = dx
            if 2 * x == K: dx = 0
            if 2 * y == K: dy = 0
            if 2 * x == K and 2 * y == K: return None
        nx, ny = x + dx, y + dy
        if not (0 <= nx <= K and 0 <= ny <= K) or (dx == 0 and dy == 0): return None
        if len(images(nx, ny, K)) != len(images(x, y, K)): return None
        o[i] = min(images(nx, ny, K))
        if len(set(o)) < len(o): return None
        return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert"); ap.add_argument("out")
    ap.add_argument("--iters", type=int, default=100000); ap.add_argument("--Dw", type=int, default=4000)
    ap.add_argument("--N", type=int, default=2000); ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0); ap.add_argument("--time", type=float, default=1800)
    ap.add_argument("--ties", type=float, default=0.5, help="probability of accepting an equal-score move")
    ap.add_argument("--sym", default="d4", choices=["d4", "none"])
    ap.add_argument("--mmax", type=int, default=None, help="fill up to this many points before polishing")
    ap.add_argument("--pkick", type=float, default=0.03, help="probability of a kick move (relocate an orbit to the worst placement)")
    a = ap.parse_args()
    global SYM; SYM = a.sym
    t0 = time.time()
    def log(msg): print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)
    P = Polisher(a.cert, a.Dw, a.N, a.threads, log)
    rng = random.Random(a.seed)
    cur = list(P.orbits)
    sc = P.score(cur, P.Dw)
    if sc is None: sys.exit("cannot score the input")
    log(f"start: k={P.k} m={P.m} orbits={len(cur)} K={P.K} D_crit={sc[0]} s={P.K/sc[0]:.6f} nfail={sc[1]}")
    if a.mmax is not None and a.mmax != P.m:
        cur, sc = P.fill(cur, a.mmax, sc, log)
        log(f"after fill: m={P.m} orbits={len(cur)} D_crit={sc[0]} s={P.K/sc[0]:.6f} nfail={sc[1]}")
    best = sc; acc = 0
    for it in range(a.iters):
        if time.time() - t0 > a.time: break
        if rng.random() < a.pkick:
            cand = P.kick(cur, rng, sc[0] - 1)
        else:
            cand = P.move(cur, rng, big=(it % 5 == 0))
        if cand is None: continue
        ok, nf = P.check(cand, sc[0] - 1)
        if not ok and nf > sc[1]:
            continue                       # worse: same or larger D_crit and more failing bins
        if not ok:
            ok2, _ = P.check(cand, sc[0])
            if not ok2: continue
            new = (sc[0], nf)
        else:
            new = P.score(cand, sc[0] - 1)
            if new is None: continue
        if new < sc or (new == sc and rng.random() < a.ties):
            cur = cand; acc += 1
            if new < sc:
                sc = new
                P.write(cur, sc[0], a.out)
                log(f"it{it}: D_crit={sc[0]} s={P.K/sc[0]:.6f} nfail={sc[1]} acc={acc} calls={P.calls}")
            else:
                sc = new
    P.write(cur, sc[0], a.out)
    log(f"final: D_crit={sc[0]} s={P.K/sc[0]:.6f} = {Fraction(P.K, sc[0])} nfail={sc[1]} acc={acc} calls={P.calls} -> {a.out}")
    print(f"POLISH k={P.k} m={P.m} s={Fraction(P.K, sc[0])} = {P.K/sc[0]:.6f} file={a.out}")


if __name__ == "__main__":
    main()
