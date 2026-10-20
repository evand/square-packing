"""Is a candidate a different local minimum from the record it beats? (10-14; s(146) / jlevy#606 was done by hand)

Inputs are two candidate.py work dirs (ours and the prior record's), each with exactsolve's cand.exact.txt,
cand.json, cand.contacts.json and the regularized view regularized.json (search/regularize, conservative policy).

  prior_packing(n, exclude_s=None)      the best public packing at n in the store (register / pending / kingbird / watch /
                                        pool sightings by a finder other than us) -> (store id, Packing, sources)
  facts(work)                           exact side, optimality, contact counts, tilt counts on the regularized view
  compare(ours_work, prior_work)        facts of both + the matched diff of the regularized views + verdict
  table_md(C, ours='this packing', prior='prior')   the "different local minimum" table for an issue draft

Matched diff: optimal assignment (scipy linear_sum_assignment) of the prior's regularized squares onto ours over the 8 box
symmetries, cost = centre distance^2 + (corner radius x angle change)^2 (pk.cluster.md); the regularizer may reorient
(orientation_g) and slide squares along flat motions, so rattlers are mapped raw -> regularized by the same matching.
A square "moves" by max(centre distance, corner radius x angle change); moves are counted over squares load-bearing in
both (a rattler's place and angle are the regularizer's choice; squares that rattle in one and not the other are counted
as swaps; exactsolve's contact and rattler sets can differ between two inputs of one minimum).  The per-angle tilt
counts are over load-bearing squares only, for the same reason.
Verdict (10-19, Evan: "be cautious about magic number thresholds.  Changes in contact graphs can produce small local
minimum differences from small moves"): decided by the exact sides alone.  Both inputs are exactsolve KKT points at
80 digits; sides that differ by more than SIDE_RES (1e-30, far above the solver's precision and far below any real gap)
are different KKT points, however little moved; equal sides with different contact graphs -> 'same side, different
contact graph' (a symmetric / degenerate twin, or one minimum reached twice with different weakly active contacts: look
closer); equal sides and equal graphs -> 'same minimum'.  The move counts (> MOVE, > 1e-2) are descriptive only, and the
contact-graph diff (exactsolve's active sets, pair level: square-square pairs and square-wall sides, matched in the raw
frames) says how different the two are.  Neither says how far apart the basins are: for that, an expansion barrier
(landscape/morph.walk_carrot; analysis/cmp154.py).
"""
from __future__ import annotations
import collections, json, math, os, re, sqlite3
import numpy as np

SAME_SIDE = 2e-9                # legacy (descriptive): the old float-side rule
SIDE_RES = 1e-30                # exact sides (80-digit KKT points) closer than this are equal
MOVE = 1e-4
TILT = 1e-6                     # degrees: below = axis-parallel
ROUND = 4                       # tilt angles grouped at 1e-4 degrees
R = 0.5 * math.sqrt(2)          # corner radius of a unit square
from .store import public_sql
PUBLIC = public_sql('g')                          # 10-17: the store's one definition
OURS = "(g.finder = 'ours' OR g.finder LIKE '%(evand)%' OR g.finder = 'evand')"


def tilt(a):
    """|angle| folded to [0, 45] degrees"""
    return abs((float(a) + 45) % 90 - 45)


# ------------------------------------------------------------------ prior record

def prior_packing(n, path=None, exclude_s=None):
    """Best public packing at n by someone else: (store id, Packing, [sighting sources]) or None.
    exclude_s: skip packings within SAME_SIDE of this side (e.g. when ours was already posted and synced)."""
    from .store import Store, DEFAULT
    db = sqlite3.connect(f'file:{path or DEFAULT}?mode=ro', uri=True)
    q = (f"SELECT p.id, p.s FROM packing p WHERE p.n=? AND p.status != 'infeasible' AND (p.pen IS NULL OR p.pen <= 1e-9) "
         f"AND EXISTS (SELECT 1 FROM sighting g WHERE g.packing = p.id AND ({PUBLIC}) AND NOT {OURS}) ORDER BY p.s LIMIT 20")
    rows = db.execute(q, (n,)).fetchall()
    rows = [r for r in rows if exclude_s is None or abs(r[1] - exclude_s) > SAME_SIDE]
    if not rows:
        return None
    same = [r for r in rows if abs(r[1] - rows[0][1]) <= SAME_SIDE]   # copies of the best minimum (side key)

    def first_posted(pid):        # credit the first public poster of that minimum: lowest issue number, then earliest
        iss, add = [], []         # public sighting
        for src, a in db.execute(f'SELECT g.source, g.added FROM sighting g WHERE g.packing=? AND ({PUBLIC}) AND NOT {OURS}',
                                 (pid,)):
            m = re.match(r'pending:#(\d+)', src or '')
            if m:
                iss.append(int(m[1]))
            add.append(a or '')
        return (min(iss) if iss else 10 ** 9, min(add) if add else '')
    pid = min((r[0] for r in same), key=first_posted)
    src = [f'{s} ({f})' for s, f in db.execute('SELECT source, finder FROM sighting WHERE packing=?', (pid,))]
    others = [r[0] for r in same if r[0] != pid]
    db.close()
    st = Store(path) if path else Store()
    try:
        return pid, st.get(pid), src + [f'later copies of the same minimum: store ids {others}'] * bool(others)
    finally:
        st.close()


def write_float(p, path):
    """ours format with plain floats (np.float64 reprs would not parse)"""
    with open(path, 'w') as f:
        f.write(f'{p.n} {float(p.s)!r}\n')
        for x, y, a in p.sq:
            f.write(f'{float(x)!r} {float(y)!r} {float(a)!r}\n')


# ------------------------------------------------------------------ matching

def _sym(sq, s, k):
    from .cluster import _sym as S
    return S(sq, s, k)


def best_match(A, sA, B, sB):
    """B onto A over the 8 box symmetries: (rms cost, symmetry k, perm, Bk scaled to A's side) with A[i] ~ Bk[perm[i]]"""
    from scipy.optimize import linear_sum_assignment
    A = np.asarray(A, float); B = np.asarray(B, float)
    best = None
    for k in range(8):
        Bk = _sym(B, sB, k); Bk[:, :2] *= sA / sB
        dx = A[:, None, 0] - Bk[None, :, 0]; dy = A[:, None, 1] - Bk[None, :, 1]
        da = np.mod(A[:, None, 2] - Bk[None, :, 2] + 45, 90) - 45
        C = dx * dx + dy * dy + (R * np.radians(da)) ** 2
        r, c = linear_sum_assignment(C)
        v = float(math.sqrt(C[r, c].mean()))
        if best is None or v < best[0]:
            best = (v, k, c, Bk)
    return best


def moves(A, Bm):
    """per-square move between matched rows: max(centre distance, corner radius x angle change); and the angle change"""
    d = np.hypot(A[:, 0] - Bm[:, 0], A[:, 1] - Bm[:, 1])
    da = np.abs(np.mod(A[:, 2] - Bm[:, 2] + 45, 90) - 45)
    return np.maximum(d, R * np.radians(da)), da


# ------------------------------------------------------------------ one packing's facts

def _regularized(work):
    j = json.load(open(os.path.join(work, 'regularized.json')))
    P = j['packing']
    if isinstance(P, str):
        P = json.loads(P)
    return float(j['S_exact'][:40]), np.array([[float(v) for v in r[:3]] for r in P['sq']]), j


def facts(work):
    """exact side, optimality, contacts, rattlers and tilt counts (regularized view) of a candidate.py work dir"""
    from .packing import read
    J = json.load(open(os.path.join(work, 'cand.json')))
    K = json.load(open(os.path.join(work, 'cand.contacts.json')))
    raw = read(os.path.join(work, 'cand.exact.txt'))
    s_reg, reg, RJ = _regularized(work)
    rawsq = np.array(raw.sq, float)
    free_raw = sorted(K.get('free') or [])
    _, _, perm, _ = best_match(reg, s_reg, rawsq, float(raw.s))       # raw rattlers -> regularized indices
    inv = {int(perm[i]): i for i in range(len(perm))}
    free = sorted(inv[f] for f in free_raw)
    t = np.array([tilt(a) for a in reg[:, 2]])
    held = np.array([i not in set(free) for i in range(len(t))])
    ang = collections.Counter(round(float(x), ROUND) for x in t[(t > TILT) & held])
    so = J.get('second_order', {})
    return dict(
        S_exact=J['S_exact'], S_cert=J.get('S_cert_decimal'),
        optimality=dict(status=so.get('status'), flat_motions=so.get('n_zero'), negative_eigs=so.get('n_neg'),
                        corner_corner=J.get('n_vv_exact', J.get('n_vv')), rattlers=len(free_raw),
                        vv_milp_dS=J.get('vv_milp_dS'), min_lambda=J.get('lambda_A_maxmin')),
        contacts=dict(collections.Counter(e[0] for e in K['equations'])),
        regularized=dict(tilted=int((t > TILT).sum()), load_bearing_tilted=int(((t > TILT) & held).sum()),
                         angles_load_bearing={f'{a:.{ROUND}f}': c for a, c in sorted(ang.items())},
                         angle_groups=RJ.get('groups_after'), orientation_g=RJ.get('orientation_g')),
        _reg=(s_reg, reg, free))


def _wall(x, y, s):
    """the wall a square touching the box is nearest to (geometry, so no symmetry-label bookkeeping)"""
    return min((x, 'L'), (s - x, 'R'), (y, 'B'), (s - y, 'T'))[1]


def contact_diff(ours_work, prior_work):
    """contact graphs of two candidate.py work dirs (exactsolve's active sets, cand.contacts.json, indices of
    cand.exact.txt), pair level, prior mapped onto ours by best_match on the raw packings"""
    from .packing import read
    G = []
    for w in (ours_work, prior_work):
        K = json.load(open(os.path.join(w, 'cand.contacts.json')))
        P = read(os.path.join(w, 'cand.exact.txt'))
        G.append((K, np.array(P.sq, float), float(P.s)))
    (Ko, Oq, sO), (Kp, Pq, sP) = G
    _, k, perm, Pm = best_match(Oq, sO, Pq, sP)          # ours[i] ~ Pm[perm[i]]
    inv = {int(perm[i]): i for i in range(len(perm))}
    def edges(K, sq, s, m):
        E = set()
        for e in K['equations']:
            if e[0] == 'W':
                i = m(e[1]); E.add(('W', i, _wall(sq[e[1], 0], sq[e[1], 1], s)))
            else:
                E.add(('C',) + tuple(sorted((m(e[1]), m(e[3])))))
        return E
    Eo = edges(Ko, Oq, sO, lambda i: i)
    Ep = edges(Kp, Pm, sO, lambda j: inv[j])
    fo, fp = set(Ko.get('free') or []), {inv[j] for j in (Kp.get('free') or [])}
    return dict(symmetry=int(k), edges_ours=len(Eo), edges_prior=len(Ep), shared=len(Eo & Ep),
                ours_only=len(Eo - Ep), prior_only=len(Ep - Eo), rattlers_ours=len(fo), rattlers_prior=len(fp),
                rattlers_shared=len(fo & fp))


def _side_eq(a, b):
    from decimal import Decimal, getcontext
    getcontext().prec = 90
    return abs(Decimal(a) - Decimal(b)) <= Decimal(SIDE_RES)


def compare(ours_work, prior_work):
    O, P = facts(ours_work), facts(prior_work)
    sO, Oq, fo = O.pop('_reg'); sP, Pq, fp = P.pop('_reg')
    md, k, perm, Pm = best_match(Oq, sO, Pq, sP)
    mv, da = moves(Oq, Pm[perm])
    fps = set(fp)
    held = [i for i in range(len(Oq)) if i not in set(fo) and int(perm[i]) not in fps]     # load-bearing in both
    swaps = [i for i in range(len(Oq)) if (i in set(fo)) != (int(perm[i]) in fps)]
    ao, ap = set(O['regularized']['angles_load_bearing']), set(P['regularized']['angles_load_bearing'])
    gap = float(P['S_exact'][:40]) - float(O['S_exact'][:40])
    n1, n2 = int((mv[held] > MOVE).sum()), int((mv[held] > 1e-2).sum())
    n3, n3all = int((mv[held] > 0.1).sum()), int((mv > 0.1).sum())     # 10-20: jlevy's refinement threshold (0.1)
    try:
        cd = contact_diff(ours_work, prior_work)
    except Exception as e:                                # older work dirs without contacts
        cd = dict(error=repr(e))
    if not _side_eq(O['S_exact'], P['S_exact']):
        verdict = 'distinct minimum'
    elif cd.get('ours_only') == 0 and cd.get('prior_only') == 0:
        verdict = 'same minimum'
    else:
        verdict = 'same side, different contact graph'
    moved = sorted(held, key=lambda i: -mv[i])[:n1]
    return dict(
        ours=O, prior=P,
        difference=dict(side_gap=gap, matching_distance_regularized=md, symmetry=int(k),
                        moved_gt_1e4=n1, moved_gt_1e2=n2, moved_gt_1e1=n3, moved_gt_1e1_all=n3all,
                        max_move=float(mv[held].max()) if held else 0.0, max_move_all=float(mv.max()),
                        max_angle_change=float(da[held].max()) if held else 0.0,
                        rattler_status_swaps=len(swaps), shared_angles=sorted(ao & ap),
                        ours_only_angles=sorted(ao - ap), prior_only_angles=sorted(ap - ao),
                        moved_squares=[dict(i=int(i), x=round(float(Oq[i, 0]), 4), y=round(float(Oq[i, 1]), 4),
                                            move=float(mv[i]), dangle=float(da[i])) for i in moved],
                        contacts=cd),
        verdict=verdict)


def summary_lines(C):
    d, o, p = C['difference'], C['ours'], C['prior']
    c = d.get('contacts') or {}
    cg = (f"contact pairs ours {c['edges_ours']} / prior {c['edges_prior']}: {c['shared']} shared, {c['ours_only']} ours "
          f"only, {c['prior_only']} prior only; rattlers {c['rattlers_shared']} of {c['rattlers_ours']} / "
          f"{c['rattlers_prior']} shared") if 'shared' in c else f"contact diff: {c.get('error')}"
    L = [f"verdict: {C['verdict']}  (exact side gap {d['side_gap']:+.3e}; {cg})",
         f"  descriptive: matching distance {d['matching_distance_regularized']:.3g}; max move {d['max_move']:.3g}; "
         f"{d['moved_gt_1e4']} / {d['moved_gt_1e2']} held squares move > 1e-4 / > 1e-2; "
         f"{d['rattler_status_swaps']} rattler / held swaps"]
    for nm, F in (('ours ', o), ('prior', p)):
        r, q = F['regularized'], F['optimality']
        L.append(f"  {nm} S* {F['S_exact'][:22]}  {q['status'] and q['status'][:40]}; flat {q['flat_motions']}, "
                 f"corner-corner {q['corner_corner']}, rattlers {q['rattlers']}, contacts {F['contacts']}; tilted {r['tilted']} "
                 f"({r['load_bearing_tilted']} held) {r['angles_load_bearing']}; groups {r['angle_groups']}")
    L.append(f"  angles shared {d['shared_angles']}, ours only {d['ours_only_angles']}, prior only {d['prior_only_angles']}")
    return L


def _angles(F):
    return '; '.join(f'{a}° ×{c}' for a, c in F['regularized']['angles_load_bearing'].items()) or '-'


def table_md(C, ours='this packing', prior='prior'):
    o, p, d = C['ours'], C['prior'], C['difference']
    c = d.get('contacts') or {}
    row = lambda k, a, b: f'| {k} | {a} | {b} |'
    def st(F):
        q = F['optimality']
        return (f"KKT, {q['status'].replace('local minimum of S', 'local minimum').split(' (')[0]}" if q['status'] else '?')
    lines = [f'| | {ours} | {prior} |', '|---|---|---|',
             row('exact side S*', o['S_exact'][:24] + '…', p['S_exact'][:24] + '…'),
             row('verdict', st(o), st(p)),
             row('corner–corner MILP descent', 'none' if (o['optimality']['vv_milp_dS'] or 0) >= 0 else o['optimality']['vv_milp_dS'],
                 'none' if (p['optimality']['vv_milp_dS'] or 0) >= 0 else p['optimality']['vv_milp_dS']),
             row('contact equations (wall + square–square)', f"{o['contacts'].get('W', 0)} + {o['contacts'].get('C', 0)}",
                 f"{p['contacts'].get('W', 0)} + {p['contacts'].get('C', 0)}"),
             row('corner–corner contacts', o['optimality']['corner_corner'], p['optimality']['corner_corner']),
             row('rattlers', o['optimality']['rattlers'], p['optimality']['rattlers']),
             row('tilted squares (regularized; load-bearing)', f"{o['regularized']['tilted']} ({o['regularized']['load_bearing_tilted']})",
                 f"{p['regularized']['tilted']} ({p['regularized']['load_bearing_tilted']})"),
             row('load-bearing tilt angles', _angles(o), _angles(p)), '',
             f"Verdict: {C['verdict']} (exact sides {'equal' if C['verdict'] != 'distinct minimum' else 'differ by %.3e' % abs(d['side_gap'])}).  "
             + (f"Contact graphs (active sets, square pairs and walls, matched over the 8 box symmetries): "
                f"{c['shared']} pairs shared, {c['ours_only']} only in {ours}, {c['prior_only']} only in {prior}; "
                f"{c['rattlers_shared']} of {c['rattlers_ours']} rattlers shared.  " if 'shared' in c else '')
             + f"Geometry: matching distance {d['matching_distance_regularized']:.3g} (regularized views), largest square "
             f"move {d['max_move']:.3g}; tilt angles shared: {', '.join(d['shared_angles']) or 'none'}."]
    return '\n'.join(lines) + '\n'
