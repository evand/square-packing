#!/usr/bin/env python3
"""What do we know that jlevy's register doesn't?  (10-23; Evan: notice these reliably; policy in SUBMITTING.md
"Beyond the register".)  For every register n, compares the register side with the store (every packing we hold, ours
or public) and with side-only public claims (Kaggle per-puzzle leaderboard: sides, no packings), and sorts each n where
something beats the register into:

  PUBLIC   a public packing (Kaggle CSV, kingbird, watched repo, a non-register pool) below the register: post it to
           jlevy crediting its finder / source (ours too if ours is the same packing).
  FIRST    our packing below the register, at a side that only a side-only claim reached first: we hold the first
           publicly known packing; claim finder credit, naming the side-only claim and its date.
  NEW      our packing below every public packing and side we know: a record candidate (candidate.py, SUBMITTING.md).
  PENDING  already claimed in an open jlevy issue (pending:#N) at or below our side: nothing to do but watch.

Our side "matches" a side-only claim within --tol (default 1e-8: float noise between copies of one packing reaches ~2e-9).  Every candidate still goes through candidate.py
(exact check, live check, relationship) before a draft.  Run `pk.py sync-register` and `pk.py watch --apply` first;
side-only lists: runs/kaggle/lb_perpuzzle_<date>.json (Kaggle staff per-puzzle leaderboard notebook output; newest
file is used; refresh by downloading that notebook's output again).

  regcheck.py [--lo 1 --hi 324] [--tol 1e-8] [--sides FILE]
"""
import argparse, glob, json, os, sqlite3
HERE = os.path.dirname(os.path.abspath(__file__))
STORE = os.environ.get('PK_STORE', os.path.join(HERE, 'runs/store.sqlite'))
FEAS = "p.status!='infeasible' AND (p.pen IS NULL OR p.pen<=1e-9)"
PUBLIC = ('register', 'pending', 'watch', 'kingbird', 'kaggle')   # + pool:<someone else>


def public(src, finder):
    k = src.split(':', 1)[0]
    return k in PUBLIC or (k == 'pool' and finder not in (None, 'ours'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lo', type=int, default=1); ap.add_argument('--hi', type=int, default=324)
    ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--sides', help='side-only claims json [{n, L, team, date}] (default: newest lb_perpuzzle)')
    a = ap.parse_args()
    fs = a.sides or max(glob.glob(os.path.join(HERE, 'runs/kaggle/lb_perpuzzle_*.json')), default=None)
    side_only = {int(r['n']): (float(r['L']), r.get('team'), r.get('date')) for r in json.load(open(fs))} if fs else {}
    db = sqlite3.connect(f'file:{STORE}?mode=ro', uri=True)
    reg = dict(db.execute(f"SELECT p.n, min(p.s) FROM packing p JOIN sighting g ON g.packing=p.id WHERE g.source='register' "
                          f"AND {FEAS} AND p.n BETWEEN ? AND ? GROUP BY p.n", (a.lo, a.hi)))
    print(f'side-only list: {fs}\n{"n":>4} {"class":8} {"register":>15} {"best":>15} {"d":>9}  details')
    for n in sorted(reg):
        rows = db.execute(f"SELECT p.id, p.s, g.source, g.finder FROM packing p JOIN sighting g ON g.packing=p.id "
                          f"WHERE p.n=? AND {FEAS} AND p.s < ? - ? ORDER BY p.s", (n, reg[n], a.tol)).fetchall()
        so = side_only.get(n)
        if not rows and not (so and so[0] < reg[n] - a.tol):
            continue
        pub = [r for r in rows if public(r[2], r[3]) and not r[2].startswith('pending')]
        pend = [r for r in rows if r[2].startswith('pending')]
        ours = [r for r in rows if not public(r[2], r[3])]
        best_pub = pub[0][1] if pub else None
        best_ours = ours[0][1] if ours else None
        cls, det = [], []
        if pend:
            cls.append('PENDING'); det.append(f'{pend[0][2]} {pend[0][1]:.10f}')
        if pub and not (pend and pend[0][1] <= pub[0][1] + a.tol):
            cls.append('PUBLIC'); det.append(f'{pub[0][2]} ({pub[0][3]}) {best_pub:.10f} id {pub[0][0]}')
        if best_ours is not None and (best_pub is None or best_ours < best_pub - a.tol) and \
                not (pend and pend[0][1] <= best_ours + a.tol):
            if so and so[0] <= best_ours + a.tol:
                cls.append('FIRST'); det.append(f'ours {best_ours:.10f} id {ours[0][0]} ({ours[0][2]}); side-only {so[0]:.10f} {so[1]} {so[2]}')
            else:
                cls.append('NEW'); det.append(f'ours {best_ours:.10f} id {ours[0][0]} ({ours[0][2]})'
                                              + (f'; side-only best {so[0]:.10f} {so[1]}' if so else ''))
        if so and so[0] < reg[n] - a.tol and (best_ours is None or so[0] < best_ours - a.tol) and \
                (best_pub is None or so[0] < best_pub - a.tol):
            cls.append('SIDEONLY'); det.append(f'side-only {so[0]:.10f} {so[1]} {so[2]} below every packing we hold')
        b = min(x for x in (best_ours, best_pub, pend[0][1] if pend else None) if x is not None) if rows else so[0]
        print(f'{n:>4} {"+".join(cls) or "-":8} {reg[n]:15.10f} {b:15.10f} {b - reg[n]:9.2e}  ' + ' | '.join(det))


if __name__ == '__main__':
    main()
