#!/usr/bin/env python3
"""Jobs with a hard per-job wall-clock limit (a hung job is killed, not waited on), as units of the compute pool.

A job stuck in C code (e.g. a HiGHS call that never returns: one s(130) census trial ran 9 h) is killed after `timeout`
s and its result is `on_timeout(job)`; a failed job's is `on_error(job, repr(exception))`; both default to None.
Results come back in job order.  Since 11-02 each job is a pk/tasks.Executor unit (no private pool; fn must be
importable, i.e. defined at module level, as for every Executor call).

  from jobpool import run_jobs
  res = run_jobs(fn, jobs, timeout=600)
"""
import sys
import time
from pk import tasks


def run_jobs(fn, jobs, timeout=600.0, on_timeout=None, on_error=None, log=sys.stderr, campaign=None, prio='normal'):
    out = []
    with tasks.Executor(campaign, prio, timeout=timeout) as ex:
        for k, (j, f) in enumerate([(j, ex.submit(fn, j)) for j in jobs]):
            try:
                out.append(f.result())
            except TimeoutError:
                print(f'jobpool: job {k} killed after {timeout:.0f}s: {j!r}', file=log, flush=True)
                out.append(on_timeout(j) if on_timeout else None)
            except Exception as e:
                print(f'jobpool: job {k} failed: {e!r}', file=log, flush=True)
                out.append(on_error(j, repr(e)) if on_error else None)
    return out


def _demo(x):
    if x == 3:
        time.sleep(100)
    if x == 4:
        raise ValueError('boom')
    return x * x


if __name__ == '__main__':
    t = time.time()
    print(run_jobs(_demo, range(8), timeout=4, on_timeout=lambda j: 'TIMEOUT', campaign='t-jobpool', prio='high'),
          f'{time.time() - t:.1f}s')
