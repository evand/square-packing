#!/usr/bin/env python3
"""Incidence-based richness at n = 110 (10-12): which certified minima each census v2-v4 replicate reached (non-start
archive entries, side key -> minimum via the store copy table); Chao2 / ICE by gap band.  Run from search/packer."""
import sqlite3, json, glob, collections, math
db=sqlite3.connect('runs/store.sqlite')
copy=dict(db.execute("select side_key,min_key from copy where n=110"))
cls=dict(db.execute("select key,cls from minimum where n=110"))
REC=10.99678339663159164
runs=sorted(glob.glob('runs/census_v[234]/*/archive.jsonl'))
inc=collections.Counter(); ab=collections.Counter(); per=[]
for f in runs:
    hit=set()
    for l in open(f):
        e=json.loads(l)
        if e['kind'].startswith('start:') or e['s']>=11: continue
        k=None
        for d in (0,-1e-9,1e-9):
            k=copy.get(f"{e['s']+d:.9f}")
            if k: break
        if k and cls.get(k)=='certified':
            hit.add(k); ab[k]+=1
    per.append(hit)
    for k in hit: inc[k]+=1
T=len(runs); print('runs',T,'distinct certified minima hit (non-start)',len(inc),'of',sum(1 for v in cls.values() if v=='certified'))
def chao2(c,T):
    Q=collections.Counter(c.values()); q1,q2=Q[1],Q[2]; S=len(c)
    est=S+((T-1)/T)*(q1*q1/(2*q2) if q2 else q1*(q1-1)/2)
    return S,q1,q2,est
def ice(c,T,cut=10):
    Q=collections.Counter(c.values()); Sfreq=sum(1 for v in c.values() if v>cut)
    inf={k:v for k,v in c.items() if v<=cut}; Sinf=len(inf); Ninf=sum(inf.values()); q1=Q[1]
    if Ninf==0: return len(c)
    Cice=1-q1/Ninf; Tinf=len([p for p in per if any(k in inf for k in p)])
    g=max(Sinf/Cice*Tinf/(Tinf-1)*sum(j*(j-1)*Q[j] for j in range(1,cut+1))/Ninf**2-1,0) if Cice>0 else float('nan')
    return Sfreq+Sinf/Cice+q1/Cice*g if Cice>0 else float('inf')
bands=[(0,1e-5),(0,1e-4),(1e-4,1e-3),(1e-3,3e-3),(3e-3,1e-2),(0,1)]
for lo,hi in bands:
    c={k:v for k,v in inc.items() if lo<=float(k)-REC<hi}
    if not c: print(lo,hi,'none'); continue
    S,q1,q2,est=chao2(c,T)
    known=sum(1 for k,v in cls.items() if v=='certified' and lo<=float(k)-REC<hi)
    print(f'gap [{lo:g},{hi:g}): known {known}, hit by runs {S}, Q1 {q1} Q2 {q2}, Chao2 {est:.0f}, ICE {ice(c,T):.0f}, incidences {sorted(c.values(),reverse=True)[:12]}')
json.dump({'T':T,'inc':inc},open('runs/census/n110/incidence.json','w'))
