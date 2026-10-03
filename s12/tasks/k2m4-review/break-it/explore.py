from cover import *
from collections import Counter
print('S',S,'nseg',len(SEGS),'poly',PX0,PX1,PY0,PY1,'dens',PDENS,'total',TOTAL,float(TOTAL))
ori=Counter(); lens=Counter()
for x0,y0,x1,y1,w in SEGS:
    ori['H' if y0==y1 else 'V' if x0==x1 else 'D']+=1; lens[abs(x1-x0)+abs(y1-y0)]+=1
print(ori,lens)
# duplicates
key=lambda s:(min((s[0],s[1]),(s[2],s[3])),max((s[0],s[1]),(s[2],s[3])))
cnt=Counter(key(s) for s in SEGS); d=[k for k,v in cnt.items() if v>1]
print('dups',len(d)); print([(tuple(map(str,k[0]+k[1])),cnt[k]) for k in d][:20])
# D4 symmetry check of measure
from collections import defaultdict
meas=defaultdict(F)
for s in SEGS: meas[key(s)]+=s[4]
def tr(p,g):
    x,y=p; S9=S
    ops=[(x,y),(S9-x,y),(x,S9-y),(S9-x,S9-y),(y,x),(S9-y,x),(y,S9-x),(S9-y,S9-x)]
    return ops[g]
for g in range(8):
    ok=all(meas.get(tuple(sorted([tr(k[0],g),tr(k[1],g)])),None)==v for k,v in meas.items())
    print('sym',g,ok)
# lines carrying mass: per horizontal line y, total mass
hl=defaultdict(F);vl=defaultdict(F)
for x0,y0,x1,y1,w in SEGS:
    if y0==y1: hl[y0]+=w
    else: vl[x0]+=w
print('H lines',sorted((str(k),float(v)) for k,v in hl.items()))
print('V lines',sorted((str(k),float(v)) for k,v in vl.items()))
print('min w',float(min(s[4] for s in SEGS)),'max',float(max(s[4] for s in SEGS)))
