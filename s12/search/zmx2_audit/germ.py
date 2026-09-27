from aud import *
import numpy as np, sys
def germscan(fc, gx, gy, cyr=0.25, ncy=101, nt=241, th=2e-6, side=1):
    """float min of mu near the tile germ (gx,gy) at theta=th: c_x = gx + t*th/2 (zeta sweeps), c_y in gy+-cyr"""
    best=(9,None)
    for cy in np.linspace(gy-cyr, gy+cyr, ncy):
        for t in np.linspace(-1.2,1.2,nt):
            cx = gx + t*th/2*side
            if not fc.adm(cx,cy,th): continue
            m=fc.mu(cx,cy,th)
            if m<best[0]: best=(m,(cx,cy,th))
    return best
if __name__=='__main__':
    cv=load(sys.argv[1]); fc=FCover(cv)
    for g in [(1.5,1.5),(2.5,2.5),(1.5,2.5),(0.5,1.5),(0.5,2.5),(0.5,0.5)]:
        print(g, germscan(fc,*g,ncy=41,nt=61))
