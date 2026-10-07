#!/usr/bin/env python3
"""Process pool with a hard per-job wall-clock limit (a hung job is killed, not waited on).

multiprocessing.Pool cannot kill a worker stuck in C code (e.g. a HiGHS call that never returns: one s(130) census
trial ran 9 h this way), so each job runs in its own forked process; overdue processes are terminated and their result
is `on_timeout(job)`.  Results come back in job order.  BLAS/OpenMP threads are pinned to 1 per job (the
parallelism is across jobs).

  from jobpool import run_jobs
  res = run_jobs(fn, jobs, procs=14, timeout=600)
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import multiprocessing as mp
import sys
import time


def _run(fn, job, q, k):
    try:
        q.put((k, 'ok', fn(job)))
    except BaseException as e:            # report, don't hang the parent
        q.put((k, 'error', repr(e)))


def run_jobs(fn, jobs, procs=14, timeout=600.0, on_timeout=None, on_error=None, log=sys.stderr):
    ctx = mp.get_context('fork')
    q = ctx.Queue()
    jobs = list(jobs)
    out = [None] * len(jobs)
    pending = list(range(len(jobs)))[::-1]
    live = {}                              # k -> (process, start time)
    done = 0
    while pending or live:
        while pending and len(live) < procs:
            k = pending.pop()
            p = ctx.Process(target=_run, args=(fn, jobs[k], q, k), daemon=True)
            p.start()
            live[k] = (p, time.time())
        try:
            k, status, val = q.get(timeout=0.5)
            out[k] = val if status == 'ok' else (on_error(jobs[k], val) if on_error else None)
            if status != 'ok':
                print(f'jobpool: job {k} failed: {val}', file=log, flush=True)
            p, _ = live.pop(k)
            p.join()
            done += 1
        except Exception:                  # queue.Empty
            pass
        now = time.time()
        for k, (p, t0) in list(live.items()):
            if now - t0 > timeout:
                p.kill(); p.join()
                live.pop(k)
                out[k] = on_timeout(jobs[k]) if on_timeout else None
                print(f'jobpool: job {k} killed after {timeout:.0f}s: {jobs[k]!r}', file=log, flush=True)
            elif not p.is_alive() and p.exitcode not in (0, None):
                live.pop(k)                # died without reporting (segfault, OOM kill)
                out[k] = on_error(jobs[k], f'exit {p.exitcode}') if on_error else None
                print(f'jobpool: job {k} died (exit {p.exitcode})', file=log, flush=True)
    return out


if __name__ == '__main__':
    def f(x):
        if x == 3:
            time.sleep(100)
        if x == 4:
            raise ValueError('boom')
        return x * x
    t = time.time()
    print(run_jobs(f, range(8), procs=4, timeout=2, on_timeout=lambda j: 'TIMEOUT'), f'{time.time() - t:.1f}s')
