"""Rigorous certificates via cell cuts in the ROTATED frame (tightest erosion).

For an angle bin [tm-dt/2, tm+dt/2] work in u = R_{-tm} c coordinates.  A cell is a square
of side eta in u-space.  An atom p is guaranteed to lie in the INTERIOR of the unit square
Q(c,t) for every (c,t) in the cell provided
        | R_{-tm}(p) - u_mid |_inf  <  hE := 1/2 - eta/2 - 0.35359*dt.
   (centre shift: |u - u_mid|_inf <= eta/2 ;  angle shift: |p-c| * dt/2 <= 0.7072*dt/2.)
Cells whose u-square meets the admissible centre region cover every placement, so if every
such cell has guaranteed weight >= 1 then every unit square in C covers weight >= 1.

Command line
------------
    python3 search/lp_search.py S FINE ETA DT TLIMIT TAG [options]

positional (unchanged from the original script; the README command still works):
    S        container side (e.g. 3.92)
    FINE     atom grid spacing
    ETA      cell side in u-space
    DT       angle bin width
    TLIMIT   wall-clock limit in seconds (ignored when --iters is given)
    TAG      output files are runs/cert_TAG.txt, runs/snap_TAG.txt, runs/xsep_TAG.txt

options (all new; every one is optional and the defaults reproduce the old behaviour):
    --iters N            run exactly N cutting-plane iterations (N LP solves) and stop.
                         Replaces the wall-clock limit, and removes the timeout on the
                         verifier subprocess, so the run is a pure function of its
                         arguments (and of the software versions -- see REPRODUCIBILITY.md).
                         The one remaining early exit, "certificate with weight < 11.999",
                         is itself deterministic.
    --seed N             seed python's `random` and numpy's global RNG, and pass N to HiGHS as
                         `random_seed`.  NOTE: the search itself draws no random numbers; the
                         seed only pins the LP solver's internal perturbations (HiGHS default
                         random_seed is 0, which is what you get without --seed).
    --dump-lp PATH       write the LP whose solution is exported (objective, constraint
                         matrix, cell and orbit definitions, solution and duals) to PATH in
                         the exact text format described in search/REPRODUCIBILITY.md.
                         `--resolve PATH` reads such a dump back, re-solves it, and compares.
    --verifier PATH      exact verifier binary (default: <repo>/verify/target/release/verify).
                         If it is missing the exact check is skipped, which changes the
                         search (the verifier's worst placements are fed back as cuts), so
                         the log records whether it was found.
    --verify-threads N   threads for the verifier (default 4).  Its output is canonically
                         sorted before use, so N does not affect the result.
    --no-verify          never call the verifier (useful for quick tests).
    --runs DIR           output directory (default runs/, created if needed).
    --highs-threads N    HiGHS `threads` option (default 1).

    python3 search/lp_search.py --resolve DUMP    re-solve an LP dump and print the objective.

Nondeterminism audit (what --iters/--seed fix, and what they do not):
    * wall-clock termination (`time.time()-t0 > tlimit`)      -> fixed by --iters
    * verifier subprocess timeout=600 s                       -> removed under --iters
    * verifier multithreading: witnesses are merged in thread completion order and only
      stably sorted by value, so ties are in arbitrary order  -> canonical full-tuple sort
    * stale runs/xsep_TAG.txt from an earlier run with the same tag was consumed on
      iteration 0                                             -> deleted at start-up
    * HiGHS: scipy forces the serial dual simplex; we additionally pin threads=1 and
      random_seed; the LP is dumped so it can be re-solved by anything
    * no python/numpy RNG is used anywhere; dict/set are only used for membership tests,
      never iterated, so hash order is irrelevant
    * floating-point results depend on numpy/scipy/HiGHS versions and on the CPU's SIMD
      code paths (np.cumsum, argsort tie-breaking, HiGHS pivoting).  Same machine + same
      wheels => byte-identical; across machines the certificate may differ (but every
      certificate produced is independently checkable).
"""
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
import math, sys, time, os, warnings, random, argparse, hashlib
from fractions import Fraction

HIGHS_OPTS = {'threads': 1, 'random_seed': 0}   # overridden from the command line

def _linprog(c, A, b):
    with warnings.catch_warnings():
        # scipy warns that `threads`/`random_seed` are not among its named options; they are
        # passed to HiGHS verbatim, which is exactly what we want.
        warnings.simplefilter("ignore")
        return linprog(c=c, A_ub=-A, b_ub=-b, bounds=(0, None), method='highs',
                       options=dict(HIGHS_OPTS))

class M:
    def __init__(s_,s,fine,eta,dt):
        s_.s=s; s_.fine=fine; s_.eta=eta; s_.dt=dt
        s_.hE=0.5-eta/2-0.35359*dt
        n=int(round(s/fine)); s_.n=n; s_.step=s/n
        s_.gx=(np.arange(n)+0.5)*s_.step
        s_.orbits=[]; s_.okey={}; s_._ca=None
        s_.cells=[]; s_.ckey=set(); s_.R=[];s_.C=[];s_.V=[]
    def add_orbit(s_,i,j):
        n=s_.n
        cand=sorted(set([(i,j),(n-1-i,j),(i,n-1-j),(n-1-i,n-1-j),(j,i),(n-1-j,i),(j,n-1-i),(n-1-j,n-1-i)]))
        key=cand[0]
        if key in s_.okey: return s_.okey[key],False
        k=len(s_.orbits); s_.okey[key]=k
        s_.orbits.append(np.array([[s_.gx[a],s_.gx[b]] for a,b in cand])); s_._ca=None
        if s_.cells:
            Q=np.array(s_.cells)   # (u0,u1,tm)
            ct=np.cos(Q[:,2]); st=np.sin(Q[:,2]); tot=np.zeros(len(Q),dtype=np.int64)
            for (px,py) in s_.orbits[k]:
                q0=px*ct+py*st; q1=-px*st+py*ct
                tot+=(np.maximum(np.abs(q0-Q[:,0]),np.abs(q1-Q[:,1]))<s_.hE)
            hit=np.nonzero(tot)[0]
            if hit.size:
                s_.R.append(hit.astype(np.int32)); s_.C.append(np.full(hit.size,k,dtype=np.int32)); s_.V.append(tot[hit].astype(float))
        return k,True
    def atoms(s_):
        if s_._ca is None or s_._ca[2]!=len(s_.orbits):
            P=np.concatenate(s_.orbits); own=np.concatenate([np.full(len(Q),k) for k,Q in enumerate(s_.orbits)])
            s_._ca=(P,own,len(s_.orbits))
        return s_._ca[0],s_._ca[1]
    def add_cells(s_,u0,u1,tm):
        """add many cells at once for one angle"""
        P,own=s_.atoms(); ct,st=math.cos(tm),math.sin(tm)
        q0=P[:,0]*ct+P[:,1]*st; q1=-P[:,0]*st+P[:,1]*ct
        newr=[];
        for a,b in zip(u0,u1):
            key=(round(a/s_.eta),round(b/s_.eta),round(tm/s_.dt))
            if key in s_.ckey: continue
            s_.ckey.add(key); r=len(s_.cells); s_.cells.append((a,b,tm))
            ok=np.maximum(np.abs(q0-a),np.abs(q1-b))<s_.hE
            cnt=np.bincount(own[ok],minlength=len(s_.orbits)); idx=np.nonzero(cnt)[0]
            if idx.size:
                s_.R.append(np.full(idx.size,r,dtype=np.int32)); s_.C.append(idx.astype(np.int32)); s_.V.append(cnt[idx].astype(float))
            newr.append(r)
        return len(newr)
    def add_points(s_,pts,hcl):
        """cuts at exact placements (cx,cy,theta) with the CLOSED square of half-side hcl"""
        P,own=s_.atoms(); added=0
        for cx,cy,tm in pts:
            key=(round(cx*1e6),round(cy*1e6),round(tm*1e7))
            if key in s_.ckey: continue
            s_.ckey.add(key); r=len(s_.cells); s_.cells.append((cx*math.cos(tm)+cy*math.sin(tm),
                                                                -cx*math.sin(tm)+cy*math.cos(tm), tm))
            ct,st=math.cos(tm),math.sin(tm)
            q0=P[:,0]*ct+P[:,1]*st; q1=-P[:,0]*st+P[:,1]*ct
            u0=cx*ct+cy*st; u1=-cx*st+cy*ct
            ok=np.maximum(np.abs(q0-u0),np.abs(q1-u1))<=hcl
            cnt=np.bincount(own[ok],minlength=len(s_.orbits)); idx=np.nonzero(cnt)[0]
            if idx.size:
                s_.R.append(np.full(idx.size,r,dtype=np.int32)); s_.C.append(idx.astype(np.int32)); s_.V.append(cnt[idx].astype(float))
            added+=1
        return added
    def prune(s_,keep):
        idx=np.full(len(s_.cells),-1,dtype=np.int64); idx[np.nonzero(keep)[0]]=np.arange(int(keep.sum()))
        R=np.concatenate(s_.R);C=np.concatenate(s_.C);V=np.concatenate(s_.V)
        nr=idx[R]; sel=nr>=0
        s_.R=[nr[sel].astype(np.int32)];s_.C=[C[sel]];s_.V=[V[sel]]
        s_.cells=[c for c,k in zip(s_.cells,keep) if k]
        s_.ckey=set((round(a/s_.eta),round(b/s_.eta),round(t/s_.dt)) for a,b,t in s_.cells)
    def solve(s_):
        R=np.concatenate(s_.R);C=np.concatenate(s_.C);V=np.concatenate(s_.V)
        A=sp.coo_matrix((V,(R,C)),shape=(len(s_.cells),len(s_.orbits))).tocsr()
        sizes=np.array([len(P) for P in s_.orbits],dtype=float)
        s_.lastA=A; s_.lastc=sizes          # kept for --dump-lp
        res=_linprog(sizes,A,np.ones(len(s_.cells)))
        if not res.success: return None
        return res.fun,res.x,-res.ineqlin.marginals

def scan_angle(P,w,s,tm,hE,eta):
    """guaranteed weight for every cell centre on the eta-grid in u-space; returns vals,U0,U1"""
    ct,st=math.cos(tm),math.sin(tm)
    q0=P[:,0]*ct+P[:,1]*st; q1=-P[:,0]*st+P[:,1]*ct
    wid=ct+abs(st); lo,hi=wid/2, s-wid/2           # admissible centre box in c-space
    if hi<=lo: return None
    corners=np.array([[lo,lo],[hi,lo],[hi,hi],[lo,hi]])
    cu0=corners[:,0]*ct+corners[:,1]*st; cu1=-corners[:,0]*st+corners[:,1]*ct
    g0=np.arange(cu0.min()-eta, cu0.max()+eta+1e-12, eta)
    g1=np.arange(cu1.min()-eta, cu1.max()+eta+1e-12, eta)
    o=np.argsort(q0); q0s=q0[o]; q1o=q1[o]; wo=w[o]
    r1=np.argsort(np.argsort(q1o))                  # rank in q1 order
    ordy=np.argsort(q1o); q1s=q1o[ordy]; wy=wo[ordy]; rank_in_q0=np.argsort(q0[o][ordy]*0+np.arange(len(o))[ordy])
    idx_in_q0=np.arange(len(o))[ordy]               # position in q0-sorted order, for q1-sorted atoms
    a=np.searchsorted(q0s,g0-hE,'right'); b=np.searchsorted(q0s,g0+hE,'left')
    vals=np.empty((len(g0),len(g1)))
    lo1=np.searchsorted(q1s,g1-hE,'right'); hi1=np.searchsorted(q1s,g1+hE,'left')
    for i in range(len(g0)):
        if b[i]<=a[i]: vals[i,:]=0.0; continue
        mask=(idx_in_q0>=a[i])&(idx_in_q0<b[i])
        cw=np.r_[0.0,np.cumsum(np.where(mask,wy,0.0))]
        vals[i,:]=cw[hi1]-cw[lo1]
    return vals,g0,g1,(lo,hi,ct,st)

def admissible(U0,U1,box,tol=0.0):
    """cell (u-square of side eta centred at (U0,U1)) may contain an admissible centre.
    Conservative: half-extent of the rotated cell in c-space is <= 0.70711*eta."""
    lo,hi,ct,st=box
    cx=U0*ct-U1*st; cy=U0*st+U1*ct
    return (cx>=lo-tol)&(cx<=hi+tol)&(cy>=lo-tol)&(cy<=hi+tol)

REPO=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_VERIFIER=os.path.join(REPO,"verify","target","release","verify")

def read_xsep(path):
    """parse the verifier's witness file and return the violated placements as (cx,cy,theta)
    in a CANONICAL order.  The verifier merges per-thread lists in completion order and then
    sorts stably by value only, so ties would otherwise come back in a thread-dependent order."""
    pts=[]
    for line in open(path):
        q=line.split()
        if len(q)==4 and float(q[0])<1.0: pts.append((float(q[0]),float(q[1]),float(q[2]),float(q[3])))
    pts.sort()
    return [(cx,cy,th) for v,th,cx,cy in pts]

def run(s,fine=0.005,eta=0.005,dt=0.005,tlimit=7200,log=print,perang=260,tag='x',
        iters=None,verifier=DEFAULT_VERIFIER,verify_threads=4,no_verify=False,runs='runs'):
    t0=time.time(); m=M(s,fine,eta,dt)
    os.makedirs(runs,exist_ok=True)
    snapf=os.path.join(runs,f"snap_{tag}.txt"); xf=os.path.join(runs,f"xsep_{tag}.txt"); winf=os.path.join(runs,f"WIN_{tag}.txt")
    for f_ in (snapf,xf):                       # stale state from an earlier run with the same tag
        if os.path.exists(f_): os.remove(f_)
    use_verifier=(not no_verify) and os.path.isfile(verifier) and os.access(verifier,os.X_OK)
    log(f"  verifier: {'using '+verifier if use_verifier else 'NOT USED (missing or --no-verify); exact cuts disabled'}")
    log(f"  mode: {'deterministic, iters='+str(iters) if iters is not None else 'wall-clock limit '+str(tlimit)+'s'}  highs={HIGHS_OPTS}")
    thetas=np.arange(dt/2, math.pi/4+dt/2+1e-12, dt)
    k=max(1,int(round(0.125/fine)))
    for i in range(k//2,m.n,k):
        for j in range(k//2,m.n,k): m.add_orbit(i,j)
    P,own=m.atoms()
    for tm in thetas[::max(1,len(thetas)//6)]:
        r=scan_angle(P,np.ones(len(P)),s,tm,m.hE,eta)
        if r is None: continue
        vals,g0,g1,box=r
        G0,G1=np.meshgrid(g0,g1,indexing='ij')
        sel=admissible(G0.ravel(),G1.ravel(),box,0.70711*eta)
        u0=G0.ravel()[sel][::40]; u1=G1.ravel()[sel][::40]
        m.add_cells(u0,u1,tm)
    best=None; last=None
    for it in range(400 if iters is None else iters):
        out=m.solve()
        if out is None: log("  LP failed"); break
        val,x,y=out
        P,own=m.atoms(); w=x[own]
        last=snapshot(m,val,x,y)
        tot_viol=0; gmin=9e9; added=0; pend=[]
        for tm in thetas:
            r=scan_angle(P,w,s,tm,m.hE,eta)
            if r is None: continue
            vals,g0,g1,box=r
            G0,G1=np.meshgrid(g0,g1,indexing='ij')
            ok=admissible(G0.ravel(),G1.ravel(),box,0.70711*eta)
            v=vals.ravel(); v=np.where(ok,v,9e9)
            gmin=min(gmin,v.min())
            bad=np.nonzero(v<1-1e-9)[0]
            tot_viol+=len(bad)
            if len(bad):
                # spread the cuts: take the worst cell in each B x B block of the grid
                B=12
                VV=v.reshape(vals.shape)
                n0,n1=VV.shape
                p0=(-n0)%B; p1=(-n1)%B
                VP=np.pad(VV,((0,p0),(0,p1)),constant_values=9e9)
                bl=VP.reshape(VP.shape[0]//B,B,VP.shape[1]//B,B).transpose(0,2,1,3).reshape(-1,B*B)
                mn=bl.min(axis=1); am=bl.argmin(axis=1)
                sb=np.nonzero(mn<1-1e-9)[0]
                if len(sb)>perang: sb=sb[np.argsort(mn[sb])[:perang]]
                nb1=VP.shape[1]//B
                bi,bj=np.divmod(sb,nb1); ii,jj=np.divmod(am[sb],B)
                r0=bi*B+ii; r1=bj*B+jj
                okk=(r0<n0)&(r1<n1); r0=r0[okk]; r1=r1[okk]
                sel=r0*n1+r1
                pend.append((G0.ravel()[sel].copy(),G1.ravel()[sel].copy(),tm))
        exact=None
        try:
            if val<12.6 and use_verifier:
                tw=write_cert(m,x,s,snapf)
                import subprocess
                r=subprocess.run([verifier,snapf,"12","2000",str(verify_threads),"6",xf],
                                 capture_output=True,text=True,timeout=None if iters is not None else 600)
                ln=[l for l in r.stdout.split(chr(10)) if l.startswith("min covered")]
                exact=float(ln[0].split("=")[-2].split()[0]) if ln else None
                if exact is not None and exact>=1.0 and tw<12.0:
                    log(f"  *** EXACT-VERIFIED CERTIFICATE  s={s}  total={tw:.6f} < 12  (mincov={exact}) ***")
                    import shutil; shutil.copy(snapf, winf)
        except Exception as e: log(f"   [exact check failed: {e}]")
        log(f"  it{it} LP={val:.5f} exact={exact} cells={len(m.cells)} orb={len(m.orbits)} minguar={gmin:.4f} viol={tot_viol} +{added} t={time.time()-t0:.0f}s")
        if tot_viol==0:
            if best is None or val<best.val:
                best=last
                log(f"  *** RIGOROUS CERTIFICATE at s={s}: TOTAL WEIGHT = {val:.6f} ***")
            if val<11.999: break
        if iters is None and time.time()-t0>tlimit: break
        act=np.nonzero(y>1e-12)[0]
        if len(act):
            nn=m.n; Pr=np.zeros((nn,nn)); GX,GY=np.meshgrid(m.gx,m.gx,indexing='ij')
            Q=np.array(m.cells)
            for r_ in act[:4000]:
                u0,u1,tm=Q[r_]; ct,st=math.cos(tm),math.sin(tm)
                cx=u0*ct-u1*st; cy=u0*st+u1*ct
                i0=max(0,int((cx-0.8)/s*nn)); i1=min(nn,int((cx+0.8)/s*nn)+2)
                j0=max(0,int((cy-0.8)/s*nn)); j1=min(nn,int((cy+0.8)/s*nn)+2)
                if i1<=i0 or j1<=j0: continue
                dx=GX[i0:i1,j0:j1]; dy=GY[i0:i1,j0:j1]
                qq0=dx*ct+dy*st; qq1=-dx*st+dy*ct
                Pr[i0:i1,j0:j1]+=y[r_]*(np.maximum(np.abs(qq0-u0),np.abs(qq1-u1))<m.hE)
            if Pr.max()>1+1e-7:
                for f_ in np.argsort(Pr.ravel())[::-1][:200]:
                    i,j=divmod(int(f_),nn)
                    if Pr[i,j]<=1+1e-7: break
                    m.add_orbit(i,j)
        # drop non-binding cells (validity is certified by the full scan, not by the LP rows)
        if len(m.cells)>26000:
            keep=(y>1e-12); keep[-9000:]=True; m.prune(keep)
        for u0,u1,tm in pend: added+=m.add_cells(u0,u1,tm)
        # extra cuts at the exact worst placements found by the integer verifier
        try:
            if os.path.exists(xf):
                pts=read_xsep(xf)
                if pts: added+=m.add_points(pts[:4000], 0.49975)
                os.remove(xf)
        except Exception as e: pass
    return best if best is not None else last

class Snap:
    """a frozen LP iterate: everything needed to write the certificate and to dump the LP.
    (`M` keeps growing after the iterate is taken, so the raw (val,x,m) triple is not enough.)"""
    __slots__=("val","x","y","P","own","step","cells","orbits","A","c")
    def __init__(s_,**kw):
        for k_,v_ in kw.items(): setattr(s_,k_,v_)

def snapshot(m,val,x,y):
    P,own=m.atoms()
    return Snap(val=val,x=x.copy(),y=y.copy(),P=P.copy(),own=own.copy(),step=m.step,
                cells=list(m.cells),orbits=list(m.orbits),A=m.lastA,c=m.lastc)

def write_cert(m,x,s,path,WD=10**7):
    """m may be an `M` (live model) or a `Snap`."""
    if isinstance(m,Snap): P,own,step=m.P,m.own,m.step
    else: P,own=m.atoms(); step=m.step
    w=x[own]
    keep=w>1e-10; P=P[keep]; w=w[keep]
    fr=Fraction(s).limit_denominator(10**6); D=int(round(2/step))
    X=np.round(P[:,0]*D).astype(np.int64); Y=np.round(P[:,1]*D).astype(np.int64)
    assert np.max(np.abs(P[:,0]*D-X))<1e-7 and (fr*D).denominator==1
    W=np.ceil(w*WD).astype(np.int64)
    with open(path,'w') as f:
        f.write(f"{fr.numerator} {fr.denominator}\n{D}\n{WD}\n{len(X)}\n")
        for a,b,c in zip(X,Y,W): f.write(f"{a} {b} {c}\n")
    return W.sum()/WD

# ----------------------------------------------------------------------------- LP dump
# Text format (see search/REPRODUCIBILITY.md).  All floats are written as C99 hex literals
# (float.hex()) so that reading them back reproduces the exact double.
def write_lp_dump(snap,path,meta):
    A=snap.A.tocoo(); A.sum_duplicates()
    with open(path,'w') as f:
        f.write("LPDUMP 1\n")
        for k_,v_ in meta.items(): f.write(f"META {k_} {v_}\n")
        f.write("PROBLEM minimize c.x  subject to  A x >= 1,  x >= 0\n")
        f.write(f"SIZE {A.shape[0]} {A.shape[1]} {A.nnz}\n")
        f.write(f"OBJ {A.shape[1]}\n")
        for v_ in snap.c: f.write(f"{int(v_)}\n")
        f.write(f"ROWS {len(snap.cells)}\n")
        for u0,u1,tm in snap.cells: f.write(f"{float(u0).hex()} {float(u1).hex()} {float(tm).hex()}\n")
        f.write(f"COLS {len(snap.orbits)}\n")
        for Q in snap.orbits:
            f.write(str(len(Q))+"".join(f" {float(a).hex()} {float(b).hex()}" for a,b in Q)+"\n")
        f.write(f"ENTRIES {A.nnz}\n")
        for r_,c_,v_ in zip(A.row,A.col,A.data): f.write(f"{r_} {c_} {int(v_)}\n")
        f.write(f"SOLUTION {float(snap.val).hex()}\n")
        f.write(f"X {len(snap.x)}\n")
        for v_ in snap.x: f.write(f"{float(v_).hex()}\n")
        f.write(f"DUAL {len(snap.y)}\n")
        for v_ in snap.y: f.write(f"{float(v_).hex()}\n")
        f.write("END\n")

def read_lp_dump(path):
    """returns dict with keys meta, c, A (csr), cells, orbits, val, x, y"""
    out={'meta':{}}
    with open(path) as f:
        L=[l.rstrip("\n") for l in f]
    i=0
    assert L[i]=="LPDUMP 1"; i+=1
    while L[i].startswith("META "):
        _,k_,v_=L[i].split(" ",2); out['meta'][k_]=v_; i+=1
    assert L[i].startswith("PROBLEM"); i+=1
    _,nr,nc,nnz=L[i].split(); nr,nc,nnz=int(nr),int(nc),int(nnz); i+=1
    assert L[i]==f"OBJ {nc}"; i+=1
    out['c']=np.array([float(L[i+j]) for j in range(nc)]); i+=nc
    assert L[i]==f"ROWS {nr}"; i+=1
    out['cells']=[tuple(float.fromhex(t) for t in L[i+j].split()) for j in range(nr)]; i+=nr
    assert L[i]==f"COLS {nc}"; i+=1
    orbs=[]
    for j in range(nc):
        t=L[i+j].split(); n_=int(t[0]); orbs.append(np.array([[float.fromhex(t[1+2*q]),float.fromhex(t[2+2*q])] for q in range(n_)]))
    out['orbits']=orbs; i+=nc
    assert L[i]==f"ENTRIES {nnz}"; i+=1
    R=np.empty(nnz,dtype=np.int64);C=np.empty(nnz,dtype=np.int64);V=np.empty(nnz)
    for j in range(nnz):
        a,b,c_=L[i+j].split(); R[j]=int(a);C[j]=int(b);V[j]=float(c_)
    i+=nnz
    out['A']=sp.coo_matrix((V,(R,C)),shape=(nr,nc)).tocsr()
    assert L[i].startswith("SOLUTION "); out['val']=float.fromhex(L[i].split()[1]); i+=1
    assert L[i]==f"X {nc}"; i+=1
    out['x']=np.array([float.fromhex(L[i+j]) for j in range(nc)]); i+=nc
    assert L[i]==f"DUAL {nr}"; i+=1
    out['y']=np.array([float.fromhex(L[i+j]) for j in range(nr)]); i+=nr
    assert L[i]=="END"
    return out

def resolve_dump(path):
    d=read_lp_dump(path)
    res=_linprog(d['c'],d['A'],np.ones(d['A'].shape[0]))
    print(f"dump: {path}")
    for k_,v_ in d['meta'].items(): print(f"  {k_}: {v_}")
    print(f"  rows={d['A'].shape[0]} cols={d['A'].shape[1]} nnz={d['A'].nnz}")
    print(f"  recorded objective  = {d['val']!r}")
    print(f"  re-solved objective = {res.fun!r}  (status {res.status})")
    print(f"  max |x - x_recorded| = {np.max(np.abs(res.x-d['x'])):.3e}")
    print(f"  x identical: {np.array_equal(res.x,d['x'])}")
    return res

def _sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f: h.update(f.read())
    return h.hexdigest()

def main(argv=None):
    ap=argparse.ArgumentParser(description="LP cutting-plane certificate search (see module docstring)")
    ap.add_argument("s",type=float,nargs='?'); ap.add_argument("fine",type=float,nargs='?')
    ap.add_argument("eta",type=float,nargs='?'); ap.add_argument("dt",type=float,nargs='?')
    ap.add_argument("tlimit",type=float,nargs='?'); ap.add_argument("tag",nargs='?')
    ap.add_argument("--iters",type=int,default=None); ap.add_argument("--seed",type=int,default=None)
    ap.add_argument("--dump-lp",default=None); ap.add_argument("--resolve",default=None)
    ap.add_argument("--verifier",default=DEFAULT_VERIFIER); ap.add_argument("--verify-threads",type=int,default=4)
    ap.add_argument("--no-verify",action="store_true"); ap.add_argument("--runs",default="runs")
    ap.add_argument("--highs-threads",type=int,default=1)
    a=ap.parse_args(argv)
    HIGHS_OPTS['threads']=a.highs_threads
    if a.seed is not None:
        random.seed(a.seed); np.random.seed(a.seed); HIGHS_OPTS['random_seed']=a.seed
    if a.resolve is not None:
        resolve_dump(a.resolve); return
    if None in (a.s,a.fine,a.eta,a.dt,a.tlimit,a.tag): ap.error("need S FINE ETA DT TLIMIT TAG (or --resolve DUMP)")
    import scipy
    print(f"lp_search: python {sys.version.split()[0]} numpy {np.__version__} scipy {scipy.__version__}  args={' '.join(sys.argv[1:] if argv is None else argv)}")
    best=run(a.s,fine=a.fine,eta=a.eta,dt=a.dt,tlimit=a.tlimit,tag=a.tag,iters=a.iters,
             verifier=a.verifier,verify_threads=a.verify_threads,no_verify=a.no_verify,runs=a.runs)
    cert=os.path.join(a.runs,f"cert_{a.tag}.txt")
    tot=write_cert(best,best.x,a.s,cert)
    print(f"RESULT s={a.s} eta={a.eta} dt={a.dt} TOTAL={best.val:.5f} exported={tot:.5f}")
    print(f"wrote {cert}  sha256={_sha256(cert)}")
    if a.dump_lp:
        meta={'s':repr(a.s),'fine':repr(a.fine),'eta':repr(a.eta),'dt':repr(a.dt),'tag':a.tag,
              'iters':a.iters,'seed':a.seed,'highs':HIGHS_OPTS,'numpy':np.__version__,'scipy':scipy.__version__,
              'certificate':cert,'certificate_sha256':_sha256(cert)}
        write_lp_dump(best,a.dump_lp,meta)
        print(f"wrote {a.dump_lp}  sha256={_sha256(a.dump_lp)}")

if __name__=="__main__":
    main()
