#!/usr/bin/env python3
"""Rigorous-in-float lower bound on the best closed point cover of the strip-leaf region.

Given a pose set S' (feasible closed unit squares in the region), solve the fractional packing LP
    max sum mu_Q  s.t.  sum_{Q ∋ v} mu_Q <= 1  for EVERY vertex v of the arrangement of S'
(square corners + pairwise edge crossings).  Max closed coverage of a finite family is attained at an
arrangement vertex, so mu is a feasible fractional packing and its mass is <= the weight of ANY closed
point cover (weak duality).  Float arithmetic only (tolerance 1e-9 on containment and crossings).

Usage: python3 search/strip/region_dual.py TAG [--cov 1.05] [--max 12000]
  reads runs/strip_region/TAG.npz (P, w, S, mu from region_lp.py), takes poses with primal coverage <= cov.
"""
import argparse, ctypes, os, subprocess, sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import region_lp as RL

CSRC = r"""
#include <math.h>
#include <stdlib.h>
/* squares k: centre, angle.  Writes crossing points of edges of squares i<j with |ci-cj|_inf<1.42 */
long crossings(int m, const double*cx, const double*cy, const double*th, double*out, long cap){
  long n=0;
  for(int i=0;i<m;i++){ double ci=cos(th[i]), si=sin(th[i]); double Vi[4][2];
    double a[4]={1,-1,-1,1}, b[4]={1,1,-1,-1};
    for(int q=0;q<4;q++){Vi[q][0]=cx[i]+0.5*(a[q]*ci-b[q]*si); Vi[q][1]=cy[i]+0.5*(a[q]*si+b[q]*ci);}
    for(int j=i+1;j<m;j++){ if(fabs(cx[i]-cx[j])>1.42||fabs(cy[i]-cy[j])>1.42) continue;
      double cj=cos(th[j]), sj=sin(th[j]); double Vj[4][2];
      for(int q=0;q<4;q++){Vj[q][0]=cx[j]+0.5*(a[q]*cj-b[q]*sj); Vj[q][1]=cy[j]+0.5*(a[q]*sj+b[q]*cj);}
      for(int p=0;p<4;p++){ double p0x=Vi[p][0],p0y=Vi[p][1],d1x=Vi[(p+1)%4][0]-p0x,d1y=Vi[(p+1)%4][1]-p0y;
        for(int q=0;q<4;q++){ double q0x=Vj[q][0],q0y=Vj[q][1],d2x=Vj[(q+1)%4][0]-q0x,d2y=Vj[(q+1)%4][1]-q0y;
          double den=d1x*d2y-d1y*d2x; if(fabs(den)<1e-12) continue;
          double rx=q0x-p0x, ry=q0y-p0y;
          double t=(rx*d2y-ry*d2x)/den, u=(rx*d1y-ry*d1x)/den;
          if(t<-1e-12||t>1+1e-12||u<-1e-12||u>1+1e-12) continue;
          if(n<cap){out[2*n]=p0x+t*d1x; out[2*n+1]=p0y+t*d1y;} n++; } } } }
  return n; }
/* containment pairs: vertices sorted by x */
long pairs(int nv, const double*vx, const double*vy, int m, const double*cx, const double*cy, const double*th,
           double hs, int*pi, int*pk, long cap){
  long n=0;
  for(int k=0;k<m;k++){ double c=cos(th[k]), s=sin(th[k]); double lo=cx[k]-0.7072, hi=cx[k]+0.7072;
    int a=0,b=nv; while(a<b){int mid=(a+b)/2; if(vx[mid]<lo) a=mid+1; else b=mid;}
    for(int j=a;j<nv && vx[j]<=hi;j++){ double dx=vx[j]-cx[k], dy=vy[j]-cy[k];
      if(fabs(dy)>0.7072) continue;
      double u=dx*c+dy*s, v=-dx*s+dy*c;
      if(fabs(u)<=hs&&fabs(v)<=hs){ if(n<cap){pi[n]=j;pk[n]=k;} n++; } } }
  return n; }
"""


def lib():
    so = os.path.join(RL.OUT, "_region_dual.so"); src = os.path.join(RL.OUT, "_region_dual.c")
    if not os.path.exists(so):
        open(src, "w").write(CSRC)
        subprocess.check_call(["gcc", "-O3", "-march=native", "-shared", "-fPIC", src, "-o", so, "-lm"])
    L = ctypes.CDLL(so)
    D = ctypes.POINTER(ctypes.c_double); I = ctypes.POINTER(ctypes.c_int)
    L.crossings.argtypes = [ctypes.c_int, D, D, D, D, ctypes.c_long]; L.crossings.restype = ctypes.c_long
    L.pairs.argtypes = [ctypes.c_int, D, D, ctypes.c_int, D, D, D, ctypes.c_double, I, I, ctypes.c_long]
    L.pairs.restype = ctypes.c_long
    return L


def P_(a, t=ctypes.c_double):
    return a.ctypes.data_as(ctypes.POINTER(t))


def support_vertices(L, S):
    S = np.ascontiguousarray(S, float)
    cx, cy, th = [np.ascontiguousarray(S[:, i]) for i in range(3)]
    cap = 20_000_000
    buf = np.zeros(2 * cap)
    n = L.crossings(len(S), P_(cx), P_(cy), P_(th), P_(buf), cap)
    assert n <= cap, n
    c, s = np.cos(th), np.sin(th)
    corners = [np.c_[cx + 0.5 * (a * c - b * s), cy + 0.5 * (a * s + b * c)] for a, b in ((1, 1), (-1, 1), (-1, -1), (1, -1))]
    return np.unique(np.round(np.vstack([buf[:2 * n].reshape(-1, 2)] + corners), 12), axis=0)


def incidence(V, S, hs):
    rows, cols = [], []
    for k, (cx, cy, th) in enumerate(S):
        c, s = np.cos(th), np.sin(th)
        dx, dy = V[:, 0] - cx, V[:, 1] - cy
        j = np.nonzero((np.abs(dx * c + dy * s) <= hs) & (np.abs(-dx * s + dy * c) <= hs))[0]
        rows.append(j); cols.append(np.full(len(j), k))
    r = np.concatenate(rows); cc = np.concatenate(cols)
    return sp.csr_matrix((np.ones(len(r)), (r, cc)), shape=(len(V), len(S)))


def dual_bound(S, V0, log=print, hs=0.5 + 1e-9, iters=60, region=None, price_every=0):
    """max sum mu over poses S, coverage <= 1 at rows; rows grown by arrangement vertices of supp(mu)."""
    L = lib()
    V = V0.copy()
    best = 0.0
    for it in range(iters):
        A = incidence(V, S, hs)
        t = time.time()
        res = linprog(-np.ones(len(S)), A_ub=A, b_ub=np.ones(A.shape[0]), bounds=(0, 1), method="highs-ipm")
        assert res.status == 0, res.message
        mu = res.x; sup = mu > 1e-10
        Ss, ms = S[sup], mu[sup]
        Vs = support_vertices(L, Ss)
        cov = np.asarray(incidence(Vs, Ss, hs) @ ms).ravel()
        M = max(cov.max(), ms.max())
        Lb = ms.sum() / M
        best = max(best, Lb)
        bad = Vs[cov > 1 + 1e-9]
        log(f"  dual it {it}: rows {len(V)}  LP {-res.fun:.5f}  supp {sup.sum()}  support vertices {len(Vs)}  "
            f"max cov {M:.5f}  -> L {Lb:.5f}  (best {best:.5f})  +rows {len(bad)}  [{time.time()-t:.0f}s]")
        newS = np.zeros((0, 3))
        if region is not None and (len(bad) == 0 or (price_every and it % price_every == price_every - 1)):
            y = -res.ineqlin.marginals
            rng = np.random.default_rng(it)
            newS, cmin, nb = RL.separate(region, V, y, hs, 0.02, 0.5, rng)
            log(f"    pricing: min pose cost {cmin:.4f}, {nb} improving poses, +{len(newS)}")
        if len(bad) == 0 and len(newS) == 0:
            return Lb, mu, S
        if len(bad):
            V = np.vstack([V, bad])
        if len(newS):
            S = np.vstack([S, newS])
    return best, mu, S


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--cov", type=float, default=1.05)
    ap.add_argument("--max", type=int, default=12000)
    ap.add_argument("--price", action="store_true", help="also price poses (column generation) on the lattice")
    ap.add_argument("--iters", type=int, default=60)
    a = ap.parse_args()
    RL.LIB = RL.kernel(); RL.NT = 2
    fn = os.path.join(RL.OUT, a.tag + ".npz")
    if os.path.exists(fn):
        d = np.load(fn); P, w, S, mu = d["P"], d["w"], d["S"], d["mu"]
    else:  # mid-run state: re-solve the primal LP
        d = np.load(os.path.join(RL.OUT, a.tag + "_state.npz")); P, S = d["P"], d["S"]
        RL.METHOD = "highs-ipm"
        val, w, mu = RL.solve_lp(RL.build_A(P, S, 0.5 + 1e-9), len(P))
        print(f"re-solved state LP: {val:.5f}")
    sup = w > 1e-12
    cv = RL.cov_poses(P[sup], w[sup], S, 0.5 + 1e-9)
    idx = np.argsort(cv)
    idx = idx[cv[idx] <= a.cov][:a.max]
    idx = np.union1d(idx, np.nonzero(mu > 1e-12)[0])
    import json
    reg = None
    if a.price:
        ar = json.load(open(os.path.join(RL.OUT, a.tag + ".json")))["args"]
        reg = RL.Region(ar["x0"], ar["xn"], ar["yn"], ar["T"], notch=not ar["nonotch"])
    Lb, mu2, S2 = dual_bound(S[idx], P, region=reg, price_every=3 if a.price else 0, iters=a.iters)
    np.savez(os.path.join(RL.OUT, a.tag + f"_dual_{a.cov}.npz"), S=S2, mu=mu2)
    print(f"RESULT {a.tag} L={Lb:.5f}")


if __name__ == "__main__":
    main()
