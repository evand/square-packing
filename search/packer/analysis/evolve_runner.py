#!/usr/bin/env python3
"""Sandboxed runner for one LLM-written constructor (analysis/evolve.py).  Run as
  unshare -rn python3 -I evolve_runner.py PROG.py N SEED MAXC
in an empty temp dir (no network; -I: no user site / cwd imports), with CPU and memory limits set by the parent.
Prints one JSON line: {"cands": [[s, [[x, y, deg], ...]], ...]} or {"error": "..."}."""
import json, math, resource, sys

resource.setrlimit(resource.RLIMIT_AS, (3 << 30, 3 << 30))
resource.setrlimit(resource.RLIMIT_CPU, (40, 45))
prog, n, seed, maxc = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
partial = len(sys.argv) > 5 and sys.argv[5] == '1'     # fill mode: fewer than n squares allowed (uprights added)
try:
    import numpy as np
    g = {'__name__': 'constructor'}
    exec(compile(open(prog).read(), 'constructor.py', 'exec'), g)
    out = g['construct'](n, np.random.default_rng(seed))
    if isinstance(out, tuple) and len(out) == 2 and not isinstance(out[0], (tuple, list)):
        out = [out]
    cands = []
    for s, sq in list(out)[:maxc]:
        a = np.asarray(sq, float).reshape(-1, 3)
        if (len(a) != n and not (partial and 1 <= len(a) <= n)) or not np.isfinite(a).all() or not (0 < float(s) < 2 * math.sqrt(n) + 4):
            continue
        cands.append([float(s), a.tolist()])
    print(json.dumps({'cands': cands} if cands else {'error': f'no valid candidate (need {n} rows of x, y, deg)'}))
except BaseException as e:                              # noqa: B036  (report anything, incl. SystemExit)
    import traceback
    tb = traceback.format_exc().strip().splitlines()
    print(json.dumps({'error': (type(e).__name__ + ': ' + str(e))[:300] + ' | ' + ' / '.join(tb[-3:])[:400]}))
