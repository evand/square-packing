import sys, random, time
sys.path.insert(0,'/home/evand/math/square-packing/public/s12/search')
from fractions import Fraction as F
import ev, zm_mixed as ZM, mixed_cover as MC
P='/home/evand/math/square-packing/public/s12/search/qx2_data/'
for fn in ('L4_k02_box7.txt','REJECTED_L2_k02_box7.txt'):
    my=ev.Cover(P+fn); zc=ZM.Cover(MC.load(P+fn))
    print(fn, 'total', my.total, float(my.total))
    random.seed(1); bad=0; t=time.time()
    poses=[(F(3,2)+F(1,10**12),F(3,2)+F(2,10**13),F(1,10**12)), (F(1,2),F(1,2),F(0)),(F(3),F(1,2),F(0)),(F(3),F(1,2)+F(1,10**9),F(1,10**9))]
    for _ in range(40):
        u=F(random.randint(-500,500),1000); c,s=ev.cs(u); e=(abs(c)+abs(s))/2
        cx=e+F(random.randint(0,10**6),10**6)*(7-2*e); cy=e+F(random.randint(0,10**6),10**6)*(7-2*e)
        poses.append((cx,cy,u))
    for cx,cy,u in poses:
        a=my.mass(cx,cy,u); 
        try: b=ZM.exact_mass(zc,cx,cy,u)
        except Exception as ex: b=repr(ex)
        if a!=b: bad+=1; print('DIFF',float(cx),float(cy),float(u),a,b)
    print('checked',len(poses),'diffs',bad, 'time',time.time()-t)
    print('germ', float(my.mass(*poses[0])))
