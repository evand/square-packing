"""Rigorous certificates via cell cuts in the ROTATED frame (tightest erosion).

For an angle bin [tm-dt/2, tm+dt/2] work in u = R_{-tm} c coordinates.  A cell is a square
of side eta in u-space.  An atom p is guaranteed to lie in the INTERIOR of the unit square
Q(c,t) for every (c,t) in the cell provided
        | R_{-tm}(p) - u_mid |_inf  <  hE := 1/2 - eta/2 - 0.35359*dt.
   (centre shift: |u - u_mid|_inf <= eta/2 ;  angle shift: |p-c| * dt/2 <= 0.7072*dt/2.)
Cells whose u-square meets the admissible centre region cover every placement, so if every
such cell has guaranteed weight >= 1 then every unit square in C covers weight >= 1.
"""
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
import math, sys, time, os
from fractions import Fraction

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
        res=linprog(c=sizes,A_ub=-A,b_ub=-np.ones(len(s_.cells)),bounds=(0,None),method='highs')
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

def run(s,fine=0.005,eta=0.005,dt=0.005,tlimit=7200,log=print,perang=260,tag='x'):
    t0=time.time(); m=M(s,fine,eta,dt)
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
    best=None
    for it in range(400):
        out=m.solve()
        if out is None: log("  LP failed"); break
        val,x,y=out
        P,own=m.atoms(); w=x[own]
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
            if val<12.6:
                tw=write_cert(m,x,s,f"runs/snap_{tag}.txt")
                import subprocess
                r=subprocess.run(["./verify/target/release/verify",f"runs/snap_{tag}.txt","12","2000","4","6",f"runs/xsep_{tag}.txt"],
                                 capture_output=True,text=True,timeout=600)
                ln=[l for l in r.stdout.split(chr(10)) if l.startswith("min covered")]
                exact=float(ln[0].split("=")[-2].split()[0]) if ln else None
                if exact is not None and exact>=1.0 and tw<12.0:
                    log(f"  *** EXACT-VERIFIED CERTIFICATE  s={s}  total={tw:.6f} < 12  (mincov={exact}) ***")
                    import shutil; shutil.copy(f"runs/snap_{tag}.txt", f"runs/WIN_{tag}.txt")
        except Exception as e: log(f"   [exact check failed: {e}]")
        log(f"  it{it} LP={val:.5f} exact={exact} cells={len(m.cells)} orb={len(m.orbits)} minguar={gmin:.4f} viol={tot_viol} +{added} t={time.time()-t0:.0f}s")
        if tot_viol==0:
            if best is None or val<best[0]:
                best=(val,x,m)
                log(f"  *** RIGOROUS CERTIFICATE at s={s}: TOTAL WEIGHT = {val:.6f} ***")
            if val<11.999: break
        if time.time()-t0>tlimit: break
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
            xf=f"runs/xsep_{tag}.txt"
            if os.path.exists(xf):
                pts=[]
                for line in open(xf):
                    q=line.split()
                    if len(q)==4 and float(q[0])<1.0: pts.append((float(q[2]),float(q[3]),float(q[1])))
                if pts: added+=m.add_points(pts[:4000], 0.49975)
                os.remove(xf)
        except Exception as e: pass
    return best if best is not None else (val,x,m)

def write_cert(m,x,s,path,WD=10**7):
    P,own=m.atoms(); w=x[own]
    keep=w>1e-10; P=P[keep]; w=w[keep]
    fr=Fraction(s).limit_denominator(10**6); D=int(round(2/m.step))
    X=np.round(P[:,0]*D).astype(np.int64); Y=np.round(P[:,1]*D).astype(np.int64)
    assert np.max(np.abs(P[:,0]*D-X))<1e-7 and (fr*D).denominator==1
    W=np.ceil(w*WD).astype(np.int64)
    with open(path,'w') as f:
        f.write(f"{fr.numerator} {fr.denominator}\n{D}\n{WD}\n{len(X)}\n")
        for a,b,c in zip(X,Y,W): f.write(f"{a} {b} {c}\n")
    return W.sum()/WD

if __name__=="__main__":
    s=float(sys.argv[1]); fine=float(sys.argv[2]); eta=float(sys.argv[3]); dt=float(sys.argv[4]); tl=float(sys.argv[5]); tag=sys.argv[6]
    best=run(s,fine=fine,eta=eta,dt=dt,tlimit=tl,tag=tag)
    val,x,m=best
    tot=write_cert(m,x,s,f"runs/cert_{tag}.txt")
    print(f"RESULT s={s} eta={eta} dt={dt} TOTAL={val:.5f} exported={tot:.5f}")
