"""Random checks of the paper proofs in R5PTS.md (witness reductions, case inequalities, conclusions)."""
import numpy as np
rng = np.random.default_rng(7); r2 = np.sqrt(2); N = 400000
def uv(d, c, s): return d[0]*c + d[1]*s, -d[0]*s + d[1]*c
bad = dict(crect_concl=0, crect_wit=0, crect_case=0, side_concl=0, side_wit=0, side_case=0, side_g=0, corner=0)
for _ in range(N):
    # --- centre-rectangle lemma (closed), theta in [0,45] for the witness proof, any theta for the conclusion
    X = rng.uniform(0, 1); Y = rng.uniform(0, min(1, r2 - X))
    if rng.random() < .3: Y = min(1, r2 - X)   # on the boundary
    A = rng.uniform(0, X); B = rng.uniform(0, Y); A_, B_ = X - A, Y - B
    th = rng.uniform(0, np.pi/2); c, s = np.cos(th), np.sin(th)
    cor = [(A, B), (-A_, B), (-A_, -B_), (A, -B_)]
    g = [max(map(abs, uv(d, c, s))) for d in cor]
    if min(g) > 0.5 + 1e-12: bad['crect_concl'] += 1
    th = rng.uniform(0, np.pi/4); c, s = np.cos(th), np.sin(th)
    U = [uv(d, c, s) for d in cor]
    # allowed witnesses: corner1 {u,v}, corner2 {v,-u}, corner3 {-u,-v}, corner4 {u,-v}
    allowed = [lambda u, v: (u, v), lambda u, v: (v, -u), lambda u, v: (-u, -v), lambda u, v: (u, -v)]
    for k, (u, v) in enumerate(U):
        if max(abs(u), abs(v)) > 0.5 and max(allowed[k](u, v)) <= 0.5: bad['crect_wit'] += 1
    # case inequalities
    if not (X*c <= 1 + 1e-12 and Y*c <= 1 + 1e-12 and (X + Y)*(c + s) <= 2 + 1e-12 and (X + Y)*(c - s) <= 2): bad['crect_case'] += 1
    # --- side strip (DS7 orientation: wall y=0, P=(0,y0), Q=(x0,y0)), closed
    y0 = rng.uniform(0.5, 1); x0 = rng.uniform(0, min(1, 2*r2 - 2*y0))
    if rng.random() < .3: x0 = min(1, 2*r2 - 2*y0)
    th = rng.uniform(0, np.pi/2); c, s = np.cos(th), np.sin(th); k = (c + s)/2
    if k > y0: continue
    m1 = rng.uniform(0, x0); m2 = rng.uniform(k, y0)
    gP = max(map(abs, uv((-m1, y0 - m2), c, s))); gQ = max(map(abs, uv((x0 - m1, y0 - m2), c, s)))
    if min(gP, gQ) > 0.5 + 1e-12: bad['side_concl'] += 1
    th = rng.uniform(0, np.pi/4); c, s = np.cos(th), np.sin(th); k = (c + s)/2
    if k <= y0:
        m2 = rng.uniform(k, y0); A, A_, B = x0 - m1, m1, y0 - m2
        uQ, vQ = uv((A, B), c, s); uP, vP = uv((-A_, B), c, s)
        if max(abs(uQ), abs(vQ)) > .5 and max(uQ, vQ) <= .5: bad['side_wit'] += 1
        if max(abs(uP), abs(vP)) > .5 and max(vP, -uP) <= .5: bad['side_wit'] += 1
        if vQ > .5 + 1e-12: bad['side_case'] += 1           # Q above top edge is impossible
        if uQ > .5 and vP > .5 and s*uQ + c*vP <= (s + c)/2 + 1e-12: bad['side_case'] += 1
        p = c + s
        if x0*(p*p - 1)/2 + y0 - p > 1e-12: bad['side_g'] += 1
    # --- corner: walls x=0,y=0, point (1,1), centre in [k,1]^2
    th = rng.uniform(0, np.pi/2); c, s = np.cos(th), np.sin(th); k = (c + s)/2
    m = rng.uniform(k, 1, 2)
    if max(map(abs, uv((1 - m[0], 1 - m[1]), c, s))) > .5 + 1e-12: bad['corner'] += 1
print(N, "samples; failures:", bad)
