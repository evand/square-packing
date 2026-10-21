"""CPU budgets (10-20; allocator prerequisite, tasks/allocator/DESIGN.md: cost is measured CPU, never wall x procs).

  tree_cpu(pid)   CPU seconds used so far by pid and every descendant, live or finished: for each live process in the
                  tree, utime + stime + cutime + cstime (/proc/<p>/stat).  A finished process that was waited for is in
                  its parent's cutime / cstime; a zombie still shows its own times.  (A finished child whose parent never
                  waits for it is lost: none of our launchers do that.)
  Budget(cpu_s)   .spent() / .left() / .over(), polled at most every `every` s (a /proc walk costs ~ms).
"""
from __future__ import annotations
import os, time

TCK = os.sysconf('SC_CLK_TCK')


def _stat(pid):
    try:
        with open(f'/proc/{pid}/stat') as f:
            v = f.read().rsplit(')', 1)[1].split()
        return int(v[1]), sum(int(x) for x in v[11:15]) / TCK      # ppid, utime + stime + cutime + cstime
    except (OSError, ValueError, IndexError):
        return None


def tree_cpu(pid=None):
    pid = pid or os.getpid()
    S = {}
    for d in os.listdir('/proc'):
        if d.isdigit():
            r = _stat(int(d))
            if r:
                S[int(d)] = r
    kids = {}
    for p, (pp, _) in S.items():
        kids.setdefault(pp, []).append(p)
    tot, todo = 0.0, [pid]
    while todo:
        p = todo.pop()
        if p in S:
            tot += S[p][1]
        todo += kids.get(p, [])
    return tot


class Budget:
    def __init__(self, cpu_s, pid=None, every=5.0):
        self.cpu_s, self.pid, self.every = cpu_s, pid or os.getpid(), every
        self.c0 = tree_cpu(self.pid); self._t = 0.0; self._last = 0.0

    def spent(self):
        if time.time() - self._t >= self.every:
            self._last = tree_cpu(self.pid) - self.c0; self._t = time.time()
        return self._last

    def left(self):
        return self.cpu_s - self.spent() if self.cpu_s else float('inf')

    def over(self):
        return self.cpu_s is not None and self.spent() >= self.cpu_s
