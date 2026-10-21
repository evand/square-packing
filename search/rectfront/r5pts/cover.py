import numpy as np, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from escape import escape, rect_pts
r2 = np.sqrt(2)
print("task 3: rectangular 4-point set, esc<=0 means unavoidable (closed)")
print("a       b        a+b-(4+r2)  esc        x  y  deg")
rows = []
for a in np.linspace(2.45, 2 + 1/r2, 13):
    bb = 4 + r2 - a
    for db in [-0.05, -1e-3, 0.0, 1e-3, 0.01]:
        b = bb + db
        if b < a - 1e-12: continue
        e, x = escape(rect_pts(a, b), a, b)
        print(f"{a:.5f} {b:.5f} {db:+.4f}  {e:+.3e}  {x[0]:.4f} {x[1]:.4f} {np.degrees(x[2])%90:.2f}", flush=True)
