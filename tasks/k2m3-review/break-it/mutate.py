"""write D4-symmetric mutations of the box-7 cover.  usage in code."""
from fractions import Fraction as F
from math import lcm
import ev
SRC='/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt'
M=F(7)
maps=[lambda x,y:(x,y),lambda x,y:(M-x,y),lambda x,y:(x,M-y),lambda x,y:(M-x,M-y),
      lambda x,y:(y,x),lambda x,y:(M-y,x),lambda x,y:(y,M-x),lambda x,y:(M-y,M-x)]
def key(x0,y0,x1,y1): return frozenset([(x0,y0),(x1,y1)])
def orbit(seg):
    x0,y0,x1,y1=seg[:4]
    return {key(*g(x0,y0),*g(x1,y1)) for g in maps}
def base():
    d=ev.load(SRC); return d
def find(segs,x0,y0,x1,y1):
    k=key(F(x0),F(y0),F(x1),F(y1))
    for s in segs:
        if key(*s[:4])==k: return s
    raise KeyError
def write(path,segs,polys,comment):
    fr=[v for s in segs for v in s[:4]]+[v for w,V in polys for p in V for v in p]
    D=lcm(*[f.denominator for f in fr])
    ws=[s[4] for s in segs]+[w for w,V in polys]
    W=lcm(*[f.denominator for f in ws])
    with open(path,'w') as f:
        f.write(f"mixed 1\n# {comment}\n7 1\n{D}\n{W}\n0\n{len(segs)}\n")
        for s in segs: f.write(" ".join(str(int(v*D)) for v in s[:4])+f" {int(s[4]*W)}\n")
        f.write(f"{len(polys)}\n")
        for w,V in polys: f.write(f"{len(V)} {int(w*W)} "+" ".join(f"{int(x*D)} {int(y*D)}" for x,y in V)+"\n")
def lighten(segs,seg,delta):
    ob=orbit(seg); out=[]
    for s in segs:
        if key(*s[:4]) in ob: s=s[:4]+(s[4]-delta,)
        out.append(s)
    return out
def shift_break(segs,sA,sB,old,new):
    """sA, sB adjacent collinear segments sharing endpoint `old` (a point); move the shared point to `new`,
    keeping each segment's density (so masses change); applied on the whole D4 orbit."""
    out=[]; obA=orbit(sA); obB=orbit(sB)
    for g in maps:
        pass
    # do it per map image
    todo={}
    for g in maps:
        for s in (sA,sB):
            a=g(*s[:2]); b=g(*s[2:4]); o=g(*old); n=g(*new)
            todo[key(*a,*b)]=(o,n)
    for s in segs:
        k=key(*s[:4])
        if k in todo:
            o,n=todo[k]; p0=(s[0],s[1]); p1=(s[2],s[3])
            L=abs(p1[0]-p0[0])+abs(p1[1]-p0[1])
            if p0==o: p0=n
            elif p1==o: p1=n
            else: raise ValueError
            L2=abs(p1[0]-p0[0])+abs(p1[1]-p0[1])
            s=(p0[0],p0[1],p1[0],p1[1],s[4]*L2/L)
        out.append(s)
    return out
if __name__=='__main__':
    d=base(); segs=d['segs']; polys=d['polys']
    # M1: wall family: lighten y=1, x in [2.8,3.0] by 1e-6
    s=find(segs,F(14,5),1,3,1); write('mut1.txt',lighten(segs,s,F(1,10**6)),polys,'M1 lighten y=1 x[2.8,3] orbit by 1e-6')
    # M2: germ (3/2,3/2): lighten y=2 x[1.2,1.4] by 1e-6
    s=find(segs,F(6,5),2,F(7,5),2); write('mut2.txt',lighten(segs,s,F(1,10**6)),polys,'M2 lighten y=2 x[1.2,1.4] orbit by 1e-6')
    # M3: breakpoint y=1 between [0.8,1.0] (dens 1.9) and [1.0,1.2] (dens .12): move 1.0 -> 0.999
    a=find(segs,F(4,5),1,1,1); b=find(segs,1,1,F(6,5),1)
    write('mut3.txt',shift_break(segs,a,b,(F(1),F(1)),(F(999,1000),F(1))),polys,'M3 breakpoint (1,1) on y=1 -> x=0.999')
    # M4: same, move to 1.001 (gains mass)
    write('mut4.txt',shift_break(segs,a,b,(F(1),F(1)),(F(1001,1000),F(1))),polys,'M4 breakpoint (1,1) on y=1 -> x=1.001')
    for i in (1,2,3,4):
        c=ev.Cover(f'mut{i}.txt'); print(i, float(c.total-ev.Cover(SRC).total))
    s=find(segs,F(8,5),F(8,5),F(9,5),F(8,5)); write('mut5.txt',lighten(segs,s,F(1,10**4)),polys,'M5 lighten y=1.6 x[1.6,1.8] orbit by 1e-4 (interior of [1,2]^2)')
    write('mut6.txt',lighten(segs,s,F(1,10**7)),polys,'M6 lighten y=1.6 x[1.6,1.8] orbit by 1e-7')
