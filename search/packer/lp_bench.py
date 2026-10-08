#!/usr/bin/env python3
"""Offline LP timing for fq's SLP polish (10-08): Clarabel (as recorded by `fq quench --dump-lp DIR`) vs HiGHS dual
simplex cold, HiGHS IPM, and HiGHS dual simplex warm-started from the previous LP's basis (statuses carried by row /
column key; new rows basic, new columns nonbasic at 0's nearest bound; basic count repaired).

LP (as fq's slp_lp): min w_s + kappa tau  s.t.  a_r . w - tau <= b_r,  -1 <= w <= 1,  tau >= 0.

  lp_bench.py DIR/lp-*.jsonl [--max 200]        (pin to one core for timing: taskset -c K)
"""
import argparse, collections, json, statistics as st, time
import numpy as np
import highspy

INF = highspy.kHighsInf
BS = highspy.HighsBasisStatus


def build(d):
    nx = d['nx']; m = len(d['rows'])
    cols = d['cols'] + ['tau']
    lo = np.full(nx, -1.0); hi = np.full(nx, 1.0); lo[-1] = 0.0; hi[-1] = INF
    c = np.zeros(nx); c[cols.index('s')] = 1.0; c[-1] = d['kappa']
    # CSR rows -> HiGHS row-wise
    start, idx, val, ub = [0], [], [], []
    for _, b, nz in d['rows']:
        acc = {}
        for j, a in nz:
            acc[j] = acc.get(j, 0.0) + a
        acc[nx - 1] = acc.get(nx - 1, 0.0) - 1.0
        for j in sorted(acc):
            idx.append(j); val.append(acc[j])
        start.append(len(idx)); ub.append(b)
    lp = highspy.HighsLp()
    lp.num_col_ = nx; lp.num_row_ = m
    lp.col_cost_ = c; lp.col_lower_ = lo; lp.col_upper_ = hi
    lp.row_lower_ = np.full(m, -INF); lp.row_upper_ = np.array(ub)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_ = np.array(start, dtype=np.int32); lp.a_matrix_.index_ = np.array(idx, dtype=np.int32)
    lp.a_matrix_.value_ = np.array(val)
    return lp, cols, [r[0] for r in d['rows']]


def solver(opts):
    h = highspy.Highs()
    h.setOptionValue('output_flag', False); h.setOptionValue('threads', 1)
    for k, v in opts.items():
        h.setOptionValue(k, v)
    return h


def solve(lp, opts, basis=None):
    h = solver(opts)
    t0 = time.perf_counter()
    h.passModel(lp)
    t1 = time.perf_counter()
    ok = None
    if basis is not None:
        ok = h.setBasis(basis)
    t2 = time.perf_counter()
    h.run()
    t3 = time.perf_counter()
    info = h.getInfo()
    return dict(obj=info.objective_function_value, it=info.simplex_iteration_count, ipm_it=info.ipm_iteration_count,
                sec=t3 - t2, pass_sec=t1 - t0, status=h.modelStatusToString(h.getModelStatus()), h=h, setb=ok)


def map_basis(prev, cols, rows, lp):
    """Basis for the new LP from the previous solve's statuses by key; repair the basic count to num_row."""
    pc, pr = prev
    b = highspy.HighsBasis()
    cs = [pc.get(k, BS.kZero) for k in cols]
    # nonbasic columns must sit at a finite bound: w in [-1, 1] -> kLower / kUpper; tau -> kLower
    cs = [BS.kLower if s == BS.kZero else s for s in cs]
    rs = [pr.get(k, BS.kBasic) for k in rows]
    m = len(rows)
    nb = sum(s == BS.kBasic for s in cs) + sum(s == BS.kBasic for s in rs)
    if nb > m:          # too many basics: demote new / basic rows first (slack nonbasic = row active at its upper bound)
        for i in range(m):
            if nb == m: break
            if rs[i] == BS.kBasic and rows[i] not in pr:
                rs[i] = BS.kUpper; nb -= 1
        for i in range(m):
            if nb == m: break
            if rs[i] == BS.kBasic:
                rs[i] = BS.kUpper; nb -= 1
    elif nb < m:        # too few: promote nonbasic rows to basic
        for i in range(m):
            if nb == m: break
            if rs[i] != BS.kBasic:
                rs[i] = BS.kBasic; nb += 1
    b.col_status = cs; b.row_status = rs; b.valid = True
    return b


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('files', nargs='+'); ap.add_argument('--max', type=int, default=10**9)
    a = ap.parse_args()
    for f in a.files:
        D = [json.loads(l) for l in open(f)][:a.max]
        R = dict(clarabel=[], dual=[], ipm=[], warm=[], pass_=[], keep_r=[], keep_c=[], it_cold=[], it_warm=[], dobj=[], setb_fail=0)
        prev = None; prev_rows = prev_cols = None; prev_step = None
        tags = {}
        T = collections.defaultdict(float)          # (method, tag) -> total s
        for d in D:
            lp, cols, rows = build(d)
            tg = d['tag']; tags[tg] = tags.get(tg, 0) + 1
            R['clarabel'].append(d['clarabel']['solve_time']); T['clarabel', tg] += d['clarabel']['solve_time']
            cd = solve(lp, dict(solver='simplex', simplex_strategy=1))
            ci = solve(lp, dict(solver='ipm', run_crossover='off'))
            R['dual'].append(cd['sec']); R['ipm'].append(ci['sec']); R['pass_'].append(cd['pass_sec']); R['it_cold'].append(cd['it'])
            T['dual', tg] += cd['sec']; T['ipm', tg] += ci['sec']
            ref = d['clarabel']['obj']
            R['dobj'].append(abs(cd['obj'] - ref) / max(1.0, abs(ref)))
            if prev is not None:
                R['keep_r'].append(len(set(rows) & prev_rows) / len(rows)); R['keep_c'].append(len(set(cols) & prev_cols) / len(cols))
                for meth, src, strat in (('warm', prev, 1), ('warm_p', prev, 4), ('warm_step', prev_step if tg == 0 else prev, 1)):
                    w = solve(lp, dict(solver='simplex', simplex_strategy=strat), basis=map_basis(src, cols, rows, lp))
                    if str(w['setb']) != 'HighsStatus.kOk': R['setb_fail'] += 1
                    if w['status'] != 'Optimal' or abs(w['obj'] - cd['obj']) > 1e-9 * max(1, abs(cd['obj'])):
                        R['bad_' + meth] = R.get('bad_' + meth, 0) + 1
                    R.setdefault(meth, []).append(w['sec']); T[meth, tg] += w['sec']
                    if meth == 'warm':
                        R['it_warm'].append(w['it']); R.setdefault('wst', []).append((w['status'], tg, round(w['sec'] * 1e3, 1), w['it']))
            bs = cd['h'].getBasis()
            prev = (dict(zip(cols, bs.col_status)), dict(zip(rows, bs.row_status)))
            if tg == 0:
                prev_step = prev
            prev_rows, prev_cols = set(rows), set(cols)
        m = [len(d['rows']) for d in D]; nx = [d['nx'] for d in D]
        md = lambda x: st.median(x) if x else float('nan')
        print(f'== {f}: n {D[0]["n"]}, {len(D)} LPs (tags {tags}), rows median {md(m)}, cols median {md(nx)}')
        for k in ('clarabel', 'dual', 'ipm', 'warm', 'warm_p', 'warm_step'):
            print(f'  {k:9s} median {md(R.get(k, [])) * 1e3:7.1f} ms  total {sum(R.get(k, [])):7.2f} s  by tag ' +
                  ' '.join(f'{t}:{T[k, t]:.2f}' for t in sorted(tags)) + f'  not-optimal/mismatch {R.get("bad_" + k, 0)}')
        print(f'  passModel median {md(R["pass_"]) * 1e3:.1f} ms; simplex iters cold {md(R["it_cold"])} warm {md(R["it_warm"])}; '
              f'setBasis failures {R["setb_fail"]}')
        q = lambda x, p: sorted(x)[int(p * (len(x) - 1))] if x else float('nan')
        print('  p90 ms: ' + ' '.join(f'{k} {q(R[k], 0.9) * 1e3:.1f}' for k in ('clarabel', 'dual', 'ipm', 'warm')))
        print('  warm status/tag:', dict(collections.Counter((x[0], x[1]) for x in R.get('wst', []))))
        print('  slowest warm (status, tag, ms, iters):', sorted(R.get('wst', []), key=lambda x: -x[2])[:5])
        print(f'  kept rows {md(R["keep_r"]):.3f} (min {min(R["keep_r"] or [0]):.3f}), kept cols {md(R["keep_c"]):.3f}; '
              f'|obj - clarabel| max {max(R["dobj"]):.1e}')


if __name__ == '__main__':
    main()
