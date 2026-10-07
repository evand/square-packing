import numpy as np
from rho import solve
lo,hi=1.10,1.15
while hi-lo>5e-4:
    m=(lo+hi)/2; v=solve(m,ny=600,nth=91,nyc=80)[0]
    (lo,hi)=(m,hi) if v<=1+1e-7 else (lo,m)
print("lines-only threshold h ~",lo,hi)
# single line: best y0 and max h with |yc-y0|<=B(th) for all poses
th=np.linspace(0,np.pi/4,2001); c,s=np.cos(th),np.sin(th); B=(c+s)/2-c*s; hw=(c+s)/2
# need y0 - hw <= B (lowest centre) and h - y0 <= B for all th
y0=np.min(hw+B); h=y0+np.min(B); print("single line y0=",y0,"max h=",h)
