# min ∫rho s.t. every unit square with centre height yc<=h (any tilt, inside y>=0) has ∫rho(y)*chord(y) dy >= 1.
# Disjoint squares have disjoint chords on each line, so N <= L*∫rho; ∫rho <= 1 and L<4 => N <= 3.
import numpy as np
from scipy.optimize import linprog
def chord(th,yc,y):
    c,s=np.cos(th),np.sin(th)
    d=np.abs(y-yc); hw=(c+s)/2
    if th<1e-12: return np.where(d<=0.5,1.0,0.0)
    # width of rotated unit square at vertical offset d: min(1/c... ) piecewise
    lo=abs(c-s)/2
    w=np.where(d<=lo, 1/max(c,s), np.where(d<=hw, (hw-d)/(c*s), 0.0))
    return w
def solve(h,ny=400,nth=46,nyc=60):
    Y=np.linspace(0,h+0.75,ny); dy=Y[1]-Y[0]
    rows=[]
    for th in np.linspace(0,np.pi/4,nth):
        hw=(np.cos(th)+np.sin(th))/2
        for yc in np.linspace(hw,h,nyc):
            rows.append(chord(th,yc,Y)*dy)
    A=-np.array(rows); b=-np.ones(len(rows))
    r=linprog(np.full(ny,dy),A_ub=A,b_ub=b,bounds=(0,None),method='highs')
    return r.fun, Y, r.x
if __name__=="__main__":
 for h in (1.0,1.05,1.1,1.15,1.175,1.19,1.2,1.207):
    v,Y,x=solve(h); sup=Y[x>1e-9]
    print(f"h={h:.3f}  min ∫rho = {v:.5f}   support {sup.min():.3f}..{sup.max():.3f}",flush=True)
