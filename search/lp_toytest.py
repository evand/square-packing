#!/usr/bin/env python3
"""zero-change warm-start sanity test for lp_bench's warm/inc paths: same instance, no new rows or
columns, carried basis -> must re-solve in 0-1 simplex iterations.  Also checks basis lengths and
(for a dump pair) the row-key / column alignment that inc relies on.

    python3 search/lp_toytest.py runs/lptest_it1                 # toy
    python3 search/lp_toytest.py runs/lp1110_it3 --prev runs/lp1110_it2 --align-only
"""
import sys, os, argparse, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import lp_bench as LB
import highspy


def zero_change(d, prev):
    LB.PREV_SOLVER = prev
    h = LB.hs_new(1, prev); h.passModel(LB.hs_model(d)); h.run(); it0 = int(h.getInfo().simplex_iteration_count)
    h.setOptionValue('solver', 'simplex'); h.setOptionValue('simplex_strategy', 1); h.run()
    print(f"[{prev}] same object re-run: extra simplex iterations = {int(h.getInfo().simplex_iteration_count) - it0}, status {str(h.getModelStatus()).split('.')[-1]}")
    b = h.getBasis(); print(f"[{prev}] basis lengths: col_status {len(b.col_status)} (model cols {d['M'].shape[1]}), row_status {len(b.row_status)} (rows {d['M'].shape[0]}), valid={b.valid}")
    r = LB.run_inc(d, d, 1); print(f"[{prev}] run_inc zero-change: it={r['it']} t={r['t']:.2f}s obj={r['obj']:.9f} prev_obj={r['prev_obj']:.9f} new_rows={r['new_rows']} new_cols={r['new_cols']}")
    r = LB.run_warm(d, d, 1); print(f"[{prev}] run_warm zero-change: it={r['it']} t={r['t']:.2f}s obj={r['obj']:.9f} basic_in_basis={r['basic_in_basis']} nr={r['nr']}")


def align(d, dprev):
    kp, kd = dprev['keys'], d['keys']
    print(f"prev rows {len(kp)}, rows {len(kd)}; prev keys are a prefix of keys: {kd[:len(kp)] == kp}; prev keys all present: {len(set(kp) - set(kd)) == 0}; duplicates in keys: {len(kd) - len(set(kd))}")
    n0, n1 = dprev['M'].shape[1], d['M'].shape[1]; nl = int(d['nl'])
    # column identity: compare the old-row block of the first n0-nl columns
    pos = {k: i for i, k in enumerate(kd)}; old = np.array([pos[k] for k in kp])
    A0 = dprev['M'][:, :n0 - nl]; A1 = d['M'][old][:, :n0 - nl]
    diff = (A0 != A1).nnz; print(f"old-rows x old-cols block identical: {diff == 0} (differing entries {diff})")
    L0 = dprev['M'][:, n0 - nl:]; L1 = d['M'][old][:, n1 - nl:]; print(f"lambda block identical: {(L0 != L1).nnz == 0}")
    print(f"c/lo/hi equal on old point cols: {np.array_equal(dprev['c'][:n0-nl], d['c'][:n0-nl])} {np.array_equal(dprev['lo'][:n0-nl], d['lo'][:n0-nl])} {np.array_equal(dprev['hi'][:n0-nl], d['hi'][:n0-nl])}; lambda bounds equal: {np.array_equal(dprev['lo'][n0-nl:], d['lo'][n1-nl:])} {np.array_equal(dprev['hi'][n0-nl:], d['hi'][n1-nl:])}")
    print(f"b equal on old rows: {np.array_equal(dprev['b'], d['b'][old])}" + (f"; rflag equal: {np.array_equal(dprev['rflag'], d['rflag'][old])}" if 'rflag' in d else ""))
    # reduced costs of the new columns under prev's recorded dual (branch.py's y): how many are dual infeasible
    if len(dprev['y']):
        y = dprev['y']; newc = np.arange(n0 - nl, n1 - nl)
        rc = d['c'][newc] - (d['M'][old][:, newc].T @ y)
        print(f"new columns {len(newc)}: reduced cost < -1e-9 under prev's y: {(rc < -1e-9).sum()} (min {rc.min():.4f})")
        rc_old = dprev['c'][:n0-nl] - (A0.T @ y); print(f"old columns: reduced cost < -1e-7 under prev's y: {(rc_old < -1e-7).sum()} (min {rc_old.min():.2e}) -- y from the dump")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('dump'); ap.add_argument('--prev', default=None); ap.add_argument('--align-only', action='store_true')
    a = ap.parse_args()
    d = LB.load(a.dump)
    if a.prev: align(d, LB.load(a.prev))
    if not a.align_only:
        for prev in ('ipm', 'simplex'): zero_change(d, prev)


if __name__ == '__main__':
    main()
