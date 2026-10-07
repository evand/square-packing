import json,glob,os,collections
from mpmath import mp,mpf,identify
mp.dps=60
d={r['n']:r for r in json.load(open('results.json'))}
rows=[]
for n in range(1,325):
    f=f'work/n-{n}/witness.exact.txt'
    if not os.path.exists(f): rows.append((n,None));continue
    L=open(f).read().split('\n')
    S=mpf(L[0].split()[1])
    th=[float(l.split()[2])%90 for l in L[1:] if l.strip()]
    tilt=[t for t in th if 1e-9<t<90-1e-9]
    ang=set()
    for t in tilt:
        a=min(t,90-t); ang.add(round(a,6))
    rows.append((n,S,len(tilt),len(ang),d[n]['status']))
json.dump([(r[0],str(r[1]) if r[1] is not None else None,*r[2:]) for r in rows],open('/tmp/claude-1000/-home-evand-math/e6f049c2-2e76-4498-80cb-37969d16b554/scratchpad/survey.json','w'))
c=collections.Counter()
for r in rows:
    if r[1] is None: continue
    c[(r[2]>0, min(r[3],5))]+=1
print(c)
# fractional parts clustering
fr=collections.defaultdict(list)
for r in rows:
    if r[1] is None or r[2]==0: continue
    fr[mp.nstr(r[1]-mp.floor(r[1]),25)].append(r[0])
dup=[(k,v) for k,v in fr.items() if len(v)>1]
print(len(fr),'distinct frac parts among tilted;',len(dup),'shared')
for k,v in sorted(dup,key=lambda x:-len(x[1]))[:30]: print(k,v)
