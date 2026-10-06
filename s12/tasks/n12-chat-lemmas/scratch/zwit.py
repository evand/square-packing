import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from zthr import mg
from ce import verts, gap
def run(args):
    b,seed=args; rng=np.random.default_rng(seed); best=(-9,None)
    for k in range(8):
        z=np.column_stack([rng.normal(0,1.5,5),rng.normal(0,1.5,5),rng.uniform(0,np.pi/2,5)]).ravel()
        for _ in range(3):
            z=minimize(lambda z:-mg(z,b),z,method='Nelder-Mead',options={'maxiter':6000}).x
        v=mg(z,b)
        if v>best[0]: best=(v,z)
    return best
if __name__=='__main__':
    bs=[3-np.sqrt(2),1.5,1.42,1.35]
    with Pool(32) as p: res=p.map(run,[(b,s) for b in bs for s in range(8)])
    for i,b in enumerate(bs):
        v,z=max(res[i*8:(i+1)*8],key=lambda t:t[0])
        P=z.reshape(5,3).copy(); P[:,:2]=b/(1+np.exp(-P[:,:2]))
        # independent recheck: SAT on all pairs with explicit polygons
        m=min(gap(verts(*P[i_]),verts(*P[j]),P[i_,2],P[j],) if False else gap(verts(*P[i_]),verts(*P[j]),P[i_,2],P[j,2]) for i_ in range(5) for j in range(i_+1,5))
        print(f"b={b:.4f} best gap {v:+.5f} recheck {m:+.5f}")
        if b==bs[0]:
            for row in P: print("   centre (%.4f, %.4f) tilt %.2f deg"%(row[0],row[1],np.degrees(row[2])%90))
