import numpy as np
r2=np.sqrt(2)
def hs(t,phi): d=phi-t; return (np.abs(np.cos(d))+np.abs(np.sin(d)))/2
def Parg(beta,H,N=1441):
    al=np.radians(np.linspace(0,90,N,endpoint=False)); ka=hs(al,0); kb=hs(beta,0)
    best=(np.inf,)
    for yname,ya in (('low',ka),('high',H-ka)):
        dy=kb-ya
        for pname,phi in (('a0',al),('a-90',al-np.pi/2),('b0',np.full_like(al,beta)),('b-90',np.full_like(al,beta-np.pi/2)),('b+90',np.full_like(al,beta+np.pi/2))):
            c=np.cos(phi); good=(c>1e-9)&(2*ka<=H)
            t=np.where(good,ka+(hs(al,phi)+hs(beta,phi)-np.sin(phi)*dy)/np.where(good,c,1),np.inf)
            i=np.argmin(t)
            if t[i]<best[0]: best=(round(t[i],4),'alpha=%.2f'%np.degrees(al[i]),yname,pname)
    return best
for H in (1.2,1.5,1.7,1.85,1.9142,1.95):
    print('H',H)
    for bd in (-45,-30,-20,-10,-5,0,5,10,20,30,45):
        b=np.radians(bd)
        if 2*hs(b,0)>H: continue
        print('  beta %4d'%bd,Parg(b,H))
