import numpy as np
from scipy.optimize import minimize
rng=np.random.default_rng(1)
def verts(x,y,t):
    c,s=np.cos(t),np.sin(t); u=np.array([c,s]); v=np.array([-s,c])
    return np.array([[x,y]])+0.5*np.array([u+v,u-v,-u-v,-u+v])
def gap(P,Q,a,b):
    best=-9
    for t in (a,a+np.pi/2,b,b+np.pi/2):
        n=np.array([np.cos(t),np.sin(t)])
        p=P@n;q=Q@n
        best=max(best,q.min()-p.max(),p.min()-q.max())
    return best
def f(z,pen=1e4):
    ax,ay,a,ex,ey,e=z
    e=np.clip(e,-np.pi/4,np.pi/4)
    c,s=np.cos(e),np.sin(e)
    obj=(c-s)*ex+(c+s)*ey-(3.91*c-1.91)
    P=verts(ax,ay,a);Q=verts(ex,ey,e)
    viol=max(0,-gap(P,Q,a,e))+max(0,-P.min())+max(0,-Q.min())+max(0,ax-ex)+max(0,ay-1.21)+max(0,ey-1.21)
    return obj+pen*viol, obj, viol
