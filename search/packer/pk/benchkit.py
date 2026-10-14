"""Race benchmarks for search mechanisms (10-13): replicates restart on success, rates per CPU, sequential stopping.

A benchmark = an event a search run can reach (e.g. "a side below the best visible minimum" = entering a withheld
basin) + arms (ways of running the search).  Each slot runs replicates back to back: a replicate ends at its first event
or at a CPU cap, and a fresh one (new seed, new start sample) takes the slot.  Information comes from events, so this
spends no CPU after an event (E0 at n = 110 spent ~80 % of its 16 CPU-h after entry).

Statistics: time-to-event treated as exponential, censored at the cap.  Rate = events / CPU (the MLE).  Arm vs control:
rate ratio with a 95 % interval from the log-ratio's standard error sqrt(1/e1 + 1/e0) (+0.5 continuity per arm);
interim checks use 99.5 % (Z_INTERIM) against the inflation of repeated looks.  The budget counts running replicates' CPU.
Sequential rule (pre-registered by the caller): check whenever every arm has spent another `check_cpu_h`; stop when
every non-control arm's interval excludes 1 (a difference) or lies inside [1/margin, margin] (no difference that
matters), or at the CPU budget.

Usage (see bench110.py):
    arms = [Arm('control', args, starts_fn), Arm('nomirror', args + ['--weights', 'mirror=0'], starts_fn)]
    race(arms, launch_fn, event_fn, out_dir, slots=16, cap_cpu_s=1800, budget_cpu_h=16, check_cpu_h=2, margin=1.5)
    summary(read_results(out_dir), control='control')
launch_fn(arm, rep_dir, seed, starts) -> argv; event_fn(rep_dir) -> event dict or None (polled).
"""
from __future__ import annotations
import dataclasses, json, math, os, signal, subprocess, time
from typing import Callable, Optional


@dataclasses.dataclass
class Arm:
    name: str
    args: list                                  # extra arguments for the search command
    starts: Callable[[int, str], list]          # (seed, rep_dir) -> start files for one replicate
    env: dict = dataclasses.field(default_factory=dict)   # extra environment for the search command (e.g. FQ_* knobs)


@dataclasses.dataclass
class Rep:
    arm: Arm
    k: int
    seed: int
    dir: str
    proc: subprocess.Popen
    t0: float


def cpu_by_pgrp():
    """{process group id: CPU seconds} over live processes: utime + stime + cutime + cstime (children they reaped).
    One /proc scan serves every slot of a poll."""
    tck = os.sysconf('SC_CLK_TCK')
    out = {}
    for p in os.listdir('/proc'):
        if not p.isdigit():
            continue
        try:
            f = open(f'/proc/{p}/stat').read().rsplit(')', 1)[1].split()
        except OSError:
            continue
        g = int(f[2])                                   # pgrp
        out[g] = out.get(g, 0.0) + (int(f[11]) + int(f[12]) + int(f[13]) + int(f[14])) / tck
    return out


def race(arms, launch_fn, event_fn, out, slots=16, cap_cpu_s=1800.0, budget_cpu_h=16.0, check_cpu_h=2.0,
         margin=1.5, control=None, poll_s=2.0, log=print):
    """Run replicates of `arms` round-robin in `slots` concurrent slots until the sequential rule or the budget stops.
    Appends one JSON line per finished replicate to out/results.jsonl."""
    os.makedirs(out, exist_ok=True)
    res_path = os.path.join(out, 'results.jsonl')
    done = read_results(out)
    nrep = {a.name: sum(1 for r in done if r['arm'] == a.name) for a in arms}
    spent = lambda: {a.name: sum(r['cpu_s'] for r in read_results(out) if r['arm'] == a.name) for a in arms}
    running: list[Rep] = []
    next_check = check_cpu_h
    stop = False

    def launch():
        # the arm that has spent the least CPU gets the slot (keeps arms balanced under unequal event times)
        sp = spent()
        C = cpu_by_pgrp()
        live = {a.name: sum(C.get(r.proc.pid, 0.0) for r in running if r.arm.name == a.name) for a in arms}
        a = min(arms, key=lambda x: (sp[x.name] + live[x.name], arms.index(x)))
        k = nrep[a.name]; nrep[a.name] += 1
        seed = 100003 * (arms.index(a) + 1) + k
        d = os.path.join(out, f'{a.name}_{k:03d}')
        os.makedirs(d, exist_ok=True)
        argv = launch_fn(a, d, seed, a.starts(seed, d))
        p = subprocess.Popen(argv, stdout=open(os.path.join(d, 'log'), 'w'), stderr=subprocess.STDOUT,
                             start_new_session=True, env={**os.environ, **a.env})
        running.append(Rep(a, k, seed, d, p, time.time()))

    last_cpu = {}                                       # replicate pid -> CPU at the last poll (exited groups vanish)

    def finish(r, ev, why):
        cpu = max(cpu_by_pgrp().get(r.proc.pid, 0.0), last_cpu.get(r.proc.pid, 0.0))
        if r.proc.poll() is None:
            try:
                os.killpg(r.proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            for _ in range(50):
                if r.proc.poll() is not None:
                    break
                time.sleep(0.1)
            if r.proc.poll() is None:
                try:
                    os.killpg(r.proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        rec = dict(arm=r.arm.name, k=r.k, seed=r.seed, dir=r.dir, cpu_s=round(cpu, 1), wall_s=round(time.time() - r.t0, 1),
                   event=ev is not None, why=why, **({'ev': ev} if ev else {}))
        with open(res_path, 'a') as f:
            f.write(json.dumps(rec) + '\n')
        running.remove(r)
        log(f"{r.arm.name:12s} #{r.k:03d} {why:7s} cpu {cpu:7.0f} s" + (f"  {ev.get('kind', '')} {ev.get('s', '')}" if ev else ''))

    while True:
        while not stop and len(running) < slots:
            launch()
        time.sleep(poll_s)
        C = cpu_by_pgrp()
        for r in running:
            if r.proc.pid in C:
                last_cpu[r.proc.pid] = C[r.proc.pid]
        for r in list(running):
            ev = event_fn(r.dir)
            if ev is not None:
                finish(r, ev, 'event')
            elif r.proc.poll() is not None:
                finish(r, None, 'exited')
            elif C.get(r.proc.pid, 0.0) >= cap_cpu_s:
                finish(r, None, 'cap')
        sp = spent()
        if not stop and min(sp.values()) / 3600 >= next_check:
            next_check += check_cpu_h
            S = summary(read_results(out), control or arms[0].name, margin=margin, quiet=True, z=Z_INTERIM)
            verdicts = [v['verdict'] for k, v in S['vs'].items()]
            log(f'check at {min(sp.values()) / 3600:.1f} CPU-h per arm: ' +
                ', '.join(f"{k} {v['ratio']:.2f} [{v['lo']:.2f}, {v['hi']:.2f}] {v['verdict']}" for k, v in S['vs'].items()))
            if verdicts and all(v != 'open' for v in verdicts):
                log('sequential rule: every arm resolved; stopping'); stop = True
        live_total = sum(C.get(r.proc.pid, 0.0) for r in running)
        if not stop and (sum(sp.values()) + live_total) / 3600 >= budget_cpu_h:
            log('budget reached; stopping'); stop = True
        if stop:
            for r in list(running):
                ev = event_fn(r.dir)
                finish(r, ev, 'event' if ev else 'stopped')
            break
    return summary(read_results(out), control or arms[0].name, margin=margin)


def read_results(out):
    p = os.path.join(out, 'results.jsonl')
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


Z_FINAL = 1.96         # 95 % at the end
Z_INTERIM = 2.81       # 99.5 % at interim checks: repeated looks at 95 % gave false "faster" verdicts in R1 (10-13: at the
                       # 1 CPU-h check two arms were "faster" at 3.3-3.5x; both ended near 1.2x, open)


def summary(R, control, margin=1.5, quiet=False, z=Z_FINAL):
    """Per arm: replicates, events, CPU-h, rate per CPU-h; vs control: rate ratio + 95 % interval + verdict
    (faster / slower: interval excludes 1; same: interval inside [1/margin, margin]; open otherwise)."""
    arms = sorted({r['arm'] for r in R}, key=lambda a: (a != control, a))
    A = {}
    for a in arms:
        X = [r for r in R if r['arm'] == a]
        e = sum(r['event'] for r in X); h = sum(r['cpu_s'] for r in X) / 3600
        A[a] = dict(reps=len(X), events=e, cpu_h=h, rate=e / h if h else float('nan'),
                    median_cpu_s_to_event=_median([r['cpu_s'] for r in X if r['event']]))
    vs = {}
    if control in A:
        c = A[control]
        for a in arms:
            if a == control:
                continue
            x = A[a]
            ratio = ((x['events'] + 0.5) / x['cpu_h']) / ((c['events'] + 0.5) / c['cpu_h']) if x['cpu_h'] and c['cpu_h'] else float('nan')
            se = math.sqrt(1 / (x['events'] + 0.5) + 1 / (c['events'] + 0.5))
            lo, hi = ratio * math.exp(-z * se), ratio * math.exp(z * se)
            verdict = 'faster' if lo > 1 else 'slower' if hi < 1 else 'same' if (lo >= 1 / margin and hi <= margin) else 'open'
            vs[a] = dict(ratio=ratio, lo=lo, hi=hi, verdict=verdict)
    if not quiet:
        print(f'{"arm":14s} {"reps":>5s} {"events":>6s} {"CPU-h":>6s} {"rate/CPU-h":>10s} {"median CPU s":>12s}')
        for a, x in A.items():
            print(f"{a:14s} {x['reps']:5d} {x['events']:6d} {x['cpu_h']:6.2f} {x['rate']:10.2f} {x['median_cpu_s_to_event'] or '-':>12}")
        for a, v in vs.items():
            print(f"  {a} / {control}: {v['ratio']:.2f} [{v['lo']:.2f}, {v['hi']:.2f}]  {v['verdict']}")
    return dict(arms=A, vs=vs)


def _median(x):
    x = sorted(x)
    return round(x[len(x) // 2]) if x else None
