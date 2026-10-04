import sys, itertools, json, numpy as np
from model import *
def seqs():
    return [tuple(int(v) for v in np.cumsum((0,)+st)*2) for st in itertools.product((2,3),repeat=4) if sum(st)==10]
X = seqs()
shard, nsh = int(sys.argv[1]), int(sys.argv[2])
cands = [(xs,ys,pm) for xs in X for ys in X for pm in itertools.permutations(range(5))]
out = open(f'shard{shard}.jsonl','w')
for i,(xs,ys,pm) in enumerate(cands):
    if i % nsh != shard: continue
    L = {(int(xs[k]), int(ys[pm[k]])) for k in range(5)}
    pos=sorted(L); bad=any(abs(a-c)<5 and abs(b-d)<5 for (a,b),(c,d) in itertools.combinations(pos,2))
    if bad: continue          # two blocks overlap (corner gap < 2.5 in both axes)
    h, sq, bl = build(12, L); val, sol = solve(h)
    rec = {'L': sorted(L), 'val': round(val) if np.isfinite(val) else -1}
    if rec['val'] == 149:
        S = [j for j in sq.values() if sol[j] > .5]
        h.addRow(-highspy.kHighsInf, len(S)-1, len(S), np.array(S), np.ones(len(S)))
        rec['unique_fill'] = solve(h)[0] < 148.5
    out.write(json.dumps(rec)+'\n'); out.flush()
