"""Excess-budget check for band LP solutions (FRIEDMAN.md section 8).  Usage: seam_budget.py <support.txt> <w>"""
import sys,ast,numpy as np
f,w=sys.argv[1],float(sys.argv[2])
H=1e-3; ys=np.arange(0,w+2,H); dens=np.zeros_like(ys)   # per-period mass density in y (atoms smeared over one bin)
seam=np.zeros_like(ys); tot=0
def orb(phi0,phi1=None):
    # mirror orbit size: self-mirror iff phi0 in {0, .5} for points/vsegs; for intervals iff symmetric
    return 1
for line in open(f):
    k,v=line.rsplit('\t',1); t=ast.literal_eval(k); v=float(v)
    if t[0]=='Pp':
        _,phi,y=t; n=1 if abs(phi)<1e-9 or abs(phi-.5)<1e-9 else 2
        i=int(round(y/H)); dens[i]+=v*n/H; tot+=v*n
        if abs(phi)<1e-9: seam[i]+=v/H
    elif t[0]=='Pv':
        _,phi,y0,y1=t; n=1 if abs(phi)<1e-9 or abs(phi-.5)<1e-9 else 2
        m=(ys>=y0)&(ys<y1); dens[m]+=v*n/(y1-y0); tot+=v*n
        if abs(phi)<1e-9: seam[m]+=v/(y1-y0)
    elif t[0]=='Ph':
        _,y,x0,x1=t; n=1 if abs(x0+x1)<1e-9 or abs(x0+x1-1)<1e-9 else 2
        i=int(round(y/H)); dens[i]+=v*n/H; tot+=v*n
    elif t[0]=='Pc':
        _,x0,x1,y0,y1=t; n=1 if abs(x0+x1)<1e-9 or abs(x0+x1-1)<1e-9 else 2
        m=(ys>=y0)&(ys<y1); dens[m]+=v*n/(y1-y0); tot+=v*n
    else: print('??',t)
print('per-period strip mass',round(tot,4),'(should be w =',w,')')
dens[ys>=w]+=1.0
cum=np.concatenate([[0],np.cumsum(dens)*H])
def win(a,b): return cum[min(int(round(b/H))+1,len(cum)-1)]-cum[int(round(a/H))]
cs=np.arange(0.5,w+0.5,0.01); E=np.array([win(c-.5,c+.5)-1 for c in cs])
print('min window excess %.4f  budget int E dc = %.4f (lemma: <= 0.5)'%(E.min(),E.sum()*0.01))
for j in range(int(w)):
    sel=(cs>=j+.5)&(cs<j+1.5)
    print(' row %d: seam %.3f   excess integral over c in [%d.5,%d.5): %.3f'%(j,seam[(ys>=j)&(ys<j+1)].sum()*H,j,j+1,E[sel].sum()*0.01))
