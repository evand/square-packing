# for box side b: max over 5 unit squares (centres in [0,b]^2, exact via sigmoid) of min pairwise SAT gap
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from ce import verts, gap
def mg(z,b):
    P=z.reshape(5,3).copy(); P[:,:2]=b/(1+np.exp(-P[:,:2])); m=9
    for i in range(5):
        for j in range(i+1,5):
            m=min(m,gap(verts(*P[i]),verts(*P[j]),P[i,2],P[j,2]))
    return m
def run(args):
    b,seed=args; rng=np.random.default_rng(seed); best=-9
    for k in range(8):
        z=np.column_stack([rng.normal(0,1.5,5),rng.normal(0,1.5,5),rng.uniform(0,np.pi/2,5)]).ravel()
        for _ in range(3):
            r=minimize(lambda z:-mg(z,b),z,method='Nelder-Mead',options={'maxiter':6000}); z=r.x
        best=max(best,mg(z,b))
    return best
if __name__=='__main__':
    bs=[1.60,1.62,1.64,1.66,1.68,1.70]
    with Pool(30) as p:
        res=p.map(run,[(b,s) for b in bs for s in range(10)])
    for i,b in enumerate(bs):
        print(f"b={b:.3f}  h=(4-b)/2={(4-b)/2:.4f}  best min-gap {max(res[i*10:(i+1)*10]):+.5f}",flush=True)
