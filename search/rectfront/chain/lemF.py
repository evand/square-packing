# Two-square lemma test: x_b >= k_a + (h_a(phi)+h_b(phi) - |sin phi| D)/cos phi for the separating normal phi,
# D = H - k_a - k_b.  Want every normal term >= m(H) = min(1.5, 2+sqrt2-H).
import numpy as np
r2=np.sqrt(2)
def hs(t,phi): d=phi-t; return (np.abs(np.cos(d))+np.abs(np.sin(d)))/2
N=721
al,be=np.meshgrid(np.radians(np.linspace(0,90,N,endpoint=False)),np.radians(np.linspace(0,90,N,endpoint=False)),indexing='ij')
ka,kb=hs(al,0),hs(be,0)
worst=1e9
for H in np.linspace(1,1.9999,400):
    ok=(2*ka<=H)&(2*kb<=H); D=H-ka-kb; m=min(1.5,2+r2-H)
    for phi in (al,al-np.pi/2,be,be-np.pi/2):
        c=np.cos(phi); good=ok&(c>1e-9)
        t=ka+(hs(al,phi)+hs(be,phi)-np.abs(np.sin(phi))*D)/np.where(good,c,1)
        g=np.where(good,t-m,np.inf); i=np.unravel_index(np.argmin(g),g.shape)
        if g[i]<worst+1e-12 and g[i]<1e-3: pass
        if g[i]<worst: worst=g[i]; arg=(H,np.degrees(al[i]),np.degrees(be[i]),np.degrees(phi[i]))
print('min slack',worst,'at H,alpha,beta,phi',arg)
