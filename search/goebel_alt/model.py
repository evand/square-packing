# Discrete model of Göbel-type packings at side S = K + r, r = sqrt(2)/2.
# Coordinates live on the index line t -> p(t) = t//2 + r*(t%2).
# Unit axis squares: indices 0..2K-1 per axis (p <= K-1+r).
# s(5) blocks (5 squares in a box of side 2+r): corner indices with p+2+r <= S.
import highspy, numpy as np, math, sys, itertools
r = math.sqrt(2)/2
def p(t): return t//2 + r*(t%2)
def overlap(a0,a1,b0,b1,eps=1e-9): return a0 < b1-eps and b0 < a1-eps

def build(K, fixed_blocks=None):
    S = K + r
    T = 2*K                         # square indices 0..T-1
    B = [t for t in range(2*K+2) if p(t)+2+r <= S+1e-9]
    h = highspy.Highs(); h.setOptionValue('output_flag', False)
    sq = {}; bl = {}
    for u in range(T):
        for v in range(T):
            sq[u,v] = len(sq)
    for a in B:
        for b in B:
            bl[a,b] = len(sq)+len(bl)
    n = len(sq)+len(bl)
    cost = np.array([-1.0]*len(sq) + [-5.0]*len(bl))
    h.addVars(n, np.zeros(n), np.ones(n))
    h.changeColsCost(n, np.arange(n), cost)
    h.changeColsIntegrality(n, np.arange(n), np.array([highspy.HighsVarType.kInteger]*n))
    def row(idx, ub):
        idx = sorted(set(idx)); h.addRow(-highspy.kHighsInf, ub, len(idx), np.array(idx), np.ones(len(idx)))
    # king-graph cliques among squares
    for u in range(T-1):
        for v in range(T-1):
            row([sq[u,v],sq[u+1,v],sq[u,v+1],sq[u+1,v+1]], 1)
    # block footprints
    X = {a: [u for u in range(T) if overlap(p(u),p(u)+1,p(a),p(a)+2+r)] for a in B}
    for (a,b),j in bl.items():
        for u in X[a]:
            for v in X[b]:
                row([j, sq[u,v]], 1)
    for (a,b),(c,d) in itertools.combinations(bl,2):
        if overlap(p(a),p(a)+2+r,p(c),p(c)+2+r) and overlap(p(b),p(b)+2+r,p(d),p(d)+2+r):
            row([bl[a,b],bl[c,d]],1)
    if fixed_blocks is not None:
        for k,j in bl.items():
            v = 1.0 if k in fixed_blocks else 0.0
            h.changeColBounds(j, v, v)
    return h, sq, bl

def solve(h):
    h.run(); sol = h.getSolution().col_value
    return -h.getInfo().objective_function_value, sol
