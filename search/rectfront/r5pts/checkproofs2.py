"""Random checks (>= 1e6 samples each) of every step in the paper proofs of §10 (lemmas In_side2, In_far, row, m' bound,
and the cell decomposition of the quadrant)."""
import numpy as np
R2 = np.sqrt(2); C0 = 0.91; X10 = 1.1
rng = np.random.default_rng(11); N = 2_000_000
def frame(n):
    th = rng.uniform(0, np.pi/2, n); return np.cos(th), np.sin(th)
def In(dx, dy, c, s, tol=1e-12):
    return (abs(dx*c+dy*s) <= .5+tol) & (abs(-dx*s+dy*c) <= .5+tol)
fails = {}
def chk(name, cond):
    bad = (~cond).sum(); fails[name] = bad; print(f"{name:60s} violations: {bad}")
# ---- In_side2 (hypotheses as stated, general y0) ----
c, s = frame(N); p = c+s; k = p/2
y0 = rng.uniform((3-R2)/2, 1, N); X = rng.uniform(0, 1, N); X = np.minimum(X, 2*(R2-y0))
A = rng.uniform(0, 1, N)*X; Ap = X-A
B = rng.uniform(-(R2-1)/2, 1, N); ok = (B <= y0-k)
chk('In_side2 conclusion (sampled hypotheses)', (In(A, B, c, s) | In(-Ap, B, c, s))[ok])
# key inequality B>=0 case: X(p^2-1)/2 + y0 - p <= 0
chk('In_side2 key g(p)<=0', (X*(p*p-1)/2 + y0 - p <= 1e-12))
chk('In_side2 linear majorant L(p)<=0', (X*(R2+1)*(p-1)/2 + y0 - p <= 1e-12))
chk('In_side2 B<0 (vQ-,uP-): X(p^2-1)/2+beta<=p/2 for beta<=(sqrt2-1)/2',
    (X*(p*p-1)/2 + (R2-1)/2 <= p/2 + 1e-12))
# ---- In_far (s or s^): beta = -B in [0, 1/2], X in [sqrt2-1, 1], X + 2 beta <= sqrt2 ----
c, s = frame(N); p = c+s
X = rng.uniform(R2-1, 1, N); beta = rng.uniform(0, 1, N)*(R2-X)/2
A = rng.uniform(0, 1, N)*X; Ap = X-A
chk('In_far conclusion', In(A, -beta, c, s) | In(-Ap, -beta, c, s))
chk('In_far (vQ-,uP-): X(p^2-1)/2 + beta <= p/2',
    X*(p*p-1)/2 + beta <= p/2 + 1e-12)
chk('In_far factor form (sqrt2-p)(1-X(sqrt2+p))/2 <= 0', (R2-p)*(1-X*(R2+p))/2 <= 1e-12)
# ---- row lemma: |dx|<=39/100, |dy|<=1-k => In ----
c, s = frame(N); k = (c+s)/2
dx = rng.uniform(-.39, .39, N); dy = rng.uniform(-1, 1, N)*(1-k)
chk('row lemma (delta=0.39)', In(dx, dy, c, s))
sv = rng.uniform(0, .78, N)
chk('row poly P(s)=(1-s)^3-(1+s)(39/50-s)^2 >= 0 on [0,0.78]', (1-sv)**3-(1+sv)*(.78-sv)**2 >= -1e-12)
# ---- m' bound in cell 3: endpoints ----
c, s = frame(N)
for al, ga, nm in [(1.9-R2, 2-R2, 'b/2=sqrt2'), (R2-1.1, .5, 'b/2=3/2')]:
    chk(f"m' endpoint {nm}: c(c-s)/2 - al c + ga s >= 0", c*(c-s)/2 - al*c + ga*s >= -1e-12)
# ---- cell decomposition on the line (closed): quadrant [k,a/2] x [k,b/2] ----
a = rng.uniform(1+R2, 4-R2, N); b = 4+R2-a; c, s = frame(N); k = (c+s)/2
x = rng.uniform(0, 1, N)*(a/2-k)+k; y = rng.uniform(0, 1, N)*(b/2-k)+k
P = lambda px, py: In(px-x, py-y, c, s)
r, S, m, rh, sh, mp = P(C0, 1), P(1, 1), P(1, b/2), P(C0, b-1), P(1, b-1), P(a-1, b/2)
cell1 = y <= 2-k; cell2 = (y >= 2-k) & (x <= X10); cell3 = (y >= 2-k) & (x >= X10)
chk('cell1 (y<=2-k): r and s', (r & S)[cell1])
chk('cell2 (y>=2-k, x<=1.1): m and (r or rh)', (m & (r | rh))[cell2])
chk('cell3 (y>=2-k, x>=1.1): m and mp and (s or sh)', (m & mp & (S | sh))[cell3])
print('cell sizes', cell1.sum(), cell2.sum(), cell3.sum())
