#!/usr/bin/env python3
"""Graft operator + level-1 benchmark on the graft events (10-15; tasks/composition-baseline/events.yaml subtype graft).

Graft (SQUISH #401 / #422: "Couzo s(102) placed into the Kingbird s(123) record, then basin-hopping search and polish"):
a whole donor packing D (m squares, side d) replaces the host's squares in a corner region of side d (+ margin): for each
host corner, each of D's 8 orientations and each margin, host squares whose centre lies in [0, d + margin]^2 at that
corner are removed; where exactly m were removed (count preserved), D goes in, and the overlapping state is quenched
(fq full ALM + polish at the host's side: the squeeze).  The best distinct results then get an explorer run.
Scheduling (Evan 10-15): every quench of every event in one 16-process queue (jobpool), then every explorer at once.

  python3 graft.py run --out runs/graft/g1 [--margins 0,0.3,0.6] [--keep 4] [--minutes 10]
  python3 graft.py score runs/graft/g1
"""
import argparse, glob, json, math, os, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np                                                     # noqa: E402
import mcmin                                                           # noqa: E402

EVY = os.path.join(HERE, '..', '..', 'tasks', 'composition-baseline', 'events.yaml')


def events():
    import yaml
    from pk import eventbench as EB
    E = yaml.safe_load(open(EVY)); E = E['events'] if isinstance(E, dict) else E
    out = []
    for e in E:
        if e.get('subtype') not in ('graft',):
            continue
        host = next((p for p in e['parents'] if 'host' in str(p.get('role', '')) or 'host' in str(p.get('description', ''))), None)
        donor = next((p for p in e['parents'] if p is not host), None)
        if not host or not donor or not host.get('store_ids') or not donor.get('store_ids'):
            out.append(dict(n=e['n'], issue=e['issue'], skip='parent missing in store')); continue
        T = e['date'].replace('Z', '')
        out.append(dict(n=e['n'], issue=str(e['issue']), T=T, new=e['new'], ante=EB.frontier_at(e['n'], T),
                        host=host['store_ids'][0], donor=donor['store_ids'][0], m=donor['n']))
    return out


def sym(sq, d, k):
    """the 8 symmetries of a packing in [0, d]^2 (x, y, deg)"""
    x, y, a = sq[:, 0], sq[:, 1], sq[:, 2]
    if k & 1: x = d - x; a = -a
    if k & 2: y = d - y; a = -a
    if k & 4: x, y, a = y, x, 90 - a
    return np.stack([x, y, np.mod(a, 90)], 1)


def corner_place(sq, d, S, c):
    """put a [0, d]^2 packing at corner c of [0, S]^2 (c: 0 LL, 1 LR, 2 UL, 3 UR)"""
    out = sq.copy()
    if c & 1: out[:, 0] += S - d
    if c & 2: out[:, 1] += S - d
    return out


def one(job):
    ev, Hs, Hsq, Ds, Dsq, c, k, mg, tmpd = job
    import hop
    from chain import SEED_Q
    Hsq = np.asarray(Hsq); D = corner_place(sym(np.asarray(Dsq), Ds, k), Ds, Hs, c)
    lo_x, hi_x = (0, Ds + mg) if not c & 1 else (Hs - Ds - mg, Hs)
    lo_y, hi_y = (0, Ds + mg) if not c & 2 else (Hs - Ds - mg, Hs)
    inside = (Hsq[:, 0] >= lo_x) & (Hsq[:, 0] <= hi_x) & (Hsq[:, 1] >= lo_y) & (Hsq[:, 1] <= hi_y)
    if inside.sum() != len(D):
        return dict(c=c, k=k, mg=mg, removed=int(inside.sum()), status='count')
    new = np.concatenate([Hsq[~inside], D])
    with tempfile.TemporaryDirectory(dir=tmpd) as tmp:
        r = hop.quench(Hs, [tuple(q) for q in new], tmp, opt=SEED_Q)
    if not r:
        return dict(c=c, k=k, mg=mg, removed=int(inside.sum()), status='infeasible')
    f = f'{tmpd}/g_c{c}_k{k}_m{mg}.txt'; mcmin.write_deg(f, r[0], r[1])
    return dict(c=c, k=k, mg=mg, removed=int(inside.sum()), status='ok', s=r[0], path=f)


def run(a):
    from pk.store import Store
    from pk import manifest
    from jobpool import run_jobs
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'graft', args=vars(a), procs=a.procs)
    st = Store()
    evs = [e for e in events() if 'skip' not in e]
    jobs, meta = [], []
    for ev in evs:
        d = f'{a.out}/n{ev["n"]}_{ev["issue"].replace(" ", "")}'; os.makedirs(d, exist_ok=True)
        H, D = st.get(ev['host']), st.get(ev['donor'])
        ev['host_s'] = float(H.s); ev['donor_s'] = float(D.s); ev['dir'] = d
        for c in range(4):
            for k in range(8):
                for mg in [float(x) for x in a.margins.split(',')]:
                    jobs.append((ev, float(H.s), H.sq.tolist(), float(D.s), D.sq.tolist(), c, k, mg, d)); meta.append(ev)
    st.close()
    print(f'{len(evs)} events, {len(jobs)} graft quenches', flush=True)
    res = run_jobs(one, jobs, procs=a.procs, timeout=900)
    by = {}
    for ev, r in zip(meta, res):
        if isinstance(r, dict):
            by.setdefault(ev['dir'], (ev, []))[1].append(r)
    # explorer stage: all events at once, spare cores to the events in order (one more replicate each)
    runs = []
    for d, (ev, R) in by.items():
        ok = sorted((r for r in R if r['status'] == 'ok'), key=lambda r: r['s'])
        keep, seen = [], []
        for r in ok:
            if any(abs(r['s'] - x) < 2e-9 for x in seen):
                continue
            seen.append(r['s']); keep.append(r['path'])
            if len(keep) >= a.keep:
                break
        ev['graft_best'] = ok[0]['s'] if ok else None
        ev['counts'] = {s: sum(1 for r in R if r['status'] == s) for s in ('ok', 'count', 'infeasible')}
        json.dump(dict(ev=ev, results=R), open(f'{d}/grafts.json', 'w'), indent=0)
        if keep and a.minutes > 0:
            runs.append((ev, keep, 0))
    extra, r = [], 1
    while len(runs) + len(extra) < a.procs and runs:
        extra += [(ev, keep, r) for ev, keep, _ in runs][:a.procs - len(runs) - len(extra)]; r += 1
    P = []
    for ev, keep, rep in runs + extra:
        cmd = [sys.executable, f'{HERE}/explore.py', '--n', str(ev['n']), '--adapt', '--reward', 'value', '--elite-share', '0.3',
               '--polish-gate', '2e-3', '--gate-ref', 'run', '--frontier-s', repr(ev['ante']), '--starts', *keep,
               '--minutes', str(a.minutes), '--procs', '1', '--seed', str(ev['n'] * 10 + rep), '--out', f'{ev["dir"]}/explore{rep}']
        P.append(subprocess.Popen(cmd, stdout=open(f'{ev["dir"]}/explore{rep}.log', 'w'), stderr=subprocess.STDOUT, cwd=HERE))
    for p in P:
        p.wait()
    man.finish()
    score(a)


def score(a):
    rows = []
    for f in sorted(glob.glob(f'{a.out}/*/grafts.json')):
        G = json.load(open(f)); ev = G['ev']; d = os.path.dirname(f)
        best = ev.get('graft_best')
        for ad in glob.glob(f'{d}/explore*/archive.jsonl'):
            for l in open(ad):
                e = json.loads(l)
                if not e.get('unpolished') and (best is None or e['s'] < best):
                    best = e['s']
        tier = 'REPLICATED' if best is not None and best <= ev['new'] + 1e-9 else \
            'BEAT-ANTE' if best is not None and best < ev['ante'] - 1e-9 else 'no'
        rows.append(dict(n=ev['n'], issue=ev['issue'], host_s=ev['host_s'], graft_best=ev.get('graft_best'), best=best,
                         ante=ev['ante'], new=ev['new'], tier=tier, counts=ev.get('counts')))
        print(f"n={ev['n']:3d} #{ev['issue']:8s} host {ev['host_s']:.8f} graft {ev.get('graft_best') or float('nan'):.8f} "
              f"final {best or float('nan'):.8f} ante {ev['ante']:.8f} posted {ev['new']:.8f} {tier} {ev.get('counts')}")
    json.dump(rows, open(f'{a.out}/score.json', 'w'), indent=1)


def live_candidates(lo=1e-7, hi=0.03, nmin=20, nmax=324):
    """live sweep: host = frontier_now(n) packing, donor = frontier_now(m) packing; count-preserving corner regions whose
    host block extent b exceeds the donor side d by lo < b - d < hi (the donor is a tighter version of a block the host
    already has).  -> [(gap, n, m, corner, margin, host_id, donor_id)]"""
    import sqlite3
    import mondrian as M
    from pk import eventbench as EB
    from pk.store import Store
    st = Store(); db = sqlite3.connect(f'file:{os.path.join(HERE, "runs", "store.sqlite")}?mode=ro', uri=True)
    F = {}
    for n in range(5, nmax + 1):
        f = EB.frontier_now(n)
        r = db.execute("SELECT id FROM packing WHERE n=? AND abs(s-?)<2e-9 AND status!='infeasible' LIMIT 1", (n, f)).fetchone()
        if r:
            F[n] = (f, np.array(st.get(r[0]).sq), r[0])
    out = []
    for n in range(nmin, nmax + 1):
        if n not in F:
            continue
        S, H, hid = F[n]; h = M.half_ext(H[:, 2])
        for m in range(5, n):
            if m not in F or F[m][0] > S - 0.9:
                continue
            d, _, did = F[m]
            for c in range(4):
                seen = set()
                for mg in (0.0, 0.3, 0.6):
                    lo_x, hi_x = (0, d + mg) if not c & 1 else (S - d - mg, S)
                    lo_y, hi_y = (0, d + mg) if not c & 2 else (S - d - mg, S)
                    ins = (H[:, 0] >= lo_x) & (H[:, 0] <= hi_x) & (H[:, 1] >= lo_y) & (H[:, 1] <= hi_y)
                    if ins.sum() != m or ins.tobytes() in seen:
                        continue
                    seen.add(ins.tobytes())
                    X, Y, hh = H[ins, 0], H[ins, 1], h[ins]
                    ex = (X + hh).max() if not c & 1 else S - (X - hh).min()
                    ey = (Y + hh).max() if not c & 2 else S - (Y - hh).min()
                    b = max(ex, ey)
                    if lo < b - d < hi:
                        out.append((float(b - d), n, m, c, mg, hid, did))
    return out


def live(a):
    from pk.store import Store
    from pk import manifest, eventbench as EB
    from jobpool import run_jobs
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'graft_live', args=vars(a), procs=a.procs)
    C = live_candidates(hi=a.hi)
    st = Store(); jobs = []
    for gap, n, m, c, mg, hid, did in C:
        H, D = st.get(hid), st.get(did)
        d = f'{a.out}/n{n}'; os.makedirs(d, exist_ok=True)
        ev = dict(n=n, m=m, gap=gap)
        for k in range(8):
            jobs.append((ev, float(H.s), H.sq.tolist(), float(D.s), D.sq.tolist(), c, k, mg, d))
    print(f'{len(C)} candidates, {len(jobs)} quenches', flush=True)
    res = run_jobs(one, jobs, procs=a.procs, timeout=900)
    hits = []
    for j, r in zip(jobs, res):
        if isinstance(r, dict) and r.get('status') == 'ok':
            n = j[0]['n']; f = EB.frontier_now(n)
            if r['s'] < f - 1e-9:
                hits.append(dict(n=n, m=j[0]['m'], s=r['s'], frontier=f, gain=f - r['s'], path=r['path']))
    json.dump(hits, open(f'{a.out}/hits.json', 'w'), indent=1)
    print(f'{len(hits)} below frontier_now:', sorted({(h["n"], round(h["gain"], 9)) for h in hits}), flush=True)
    man.finish()


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--out', required=True); r.add_argument('--margins', default='0,0.3,0.6')
    r.add_argument('--keep', type=int, default=4); r.add_argument('--minutes', type=float, default=10)
    r.add_argument('--procs', type=int, default=16)
    s = sp.add_parser('score'); s.add_argument('out')
    sp.add_parser('list')
    lv = sp.add_parser('live'); lv.add_argument('--out', required=True); lv.add_argument('--hi', type=float, default=0.03)
    lv.add_argument('--procs', type=int, default=16)
    a = ap.parse_args()
    if a.cmd == 'list':
        for e in events(): print(e)
    else:
        {'run': run, 'score': score, 'live': live}[a.cmd](a)
