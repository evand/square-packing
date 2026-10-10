"""Proposal pipeline: screen -> classify -> (staged) polish -> grid test.  One function, explicit policy, full record.

Statuses (as explore.py): fail | screen-discard (above smax) | screen-grid | return (screened side within `same` of a known
basin) | new? (polished, or screened-only above pol_above) | late-return (staged polish reached a known basin early).
"""
from __future__ import annotations
import bisect, math, random, time
from dataclasses import dataclass, field, asdict
from .engine import QOpt, SCREEN, POLISH, fq
from .packing import Packing


@dataclass(frozen=True)
class Policy:
    loosen: tuple = (1.0, 1.02, 1.05)        # screen dilation drawn per proposal (pressure release)
    screen: QOpt = SCREEN
    polish: QOpt = POLISH
    smax_above_k: float = 0.05               # stepping-stone ceiling: k + this
    same: float = 1e-6                       # screened side within this of a known basin = return (no polish)
    pol_above: float | None = None           # polish only if screened side < best + this (None: always)
    stage_pit: int = 0                       # > 0: polish in chunks of this many SLP iterations, checking known basins
    late_same: float = 1e-7                  # staged polish: within this of a known basin (and stalled) = late-return
    late_slope: float = 1e-8                 # ... stalled = last chunk improved by less than this
    grid_check: bool = True

    def as_dict(self):
        d = asdict(self); d['screen'] = self.screen.args(); d['polish'] = self.polish.args(); return d


class Known:
    """Sorted list of known basin sides (one n) for nearest lookups."""
    def __init__(self, sides=()):
        self.v = sorted(sides)

    def near(self, s, tol):
        i = bisect.bisect_left(self.v, s - tol)
        return self.v[i] if i < len(self.v) and self.v[i] <= s + tol else None

    def add(self, s):
        bisect.insort(self.v, s)


def evaluate(prop: Packing, policy: Policy, known: Known, best: float, seed: int, k: int | None = None) -> dict:
    """Run one proposal through the pipeline.  Returns a record (no coordinates; the final packing under 'p')."""
    rng = random.Random(seed)
    k = k or prop.k
    smax = k + policy.smax_above_k
    loosen = rng.choice(policy.loosen)
    W = fq()
    rec = dict(loosen=loosen)
    t0 = time.time()
    r1 = W.quench(prop, policy.screen.but(loosen=loosen))
    rec['t_screen'] = time.time() - t0
    if not r1.ok:
        return dict(rec, st='fail', sec=time.time() - t0, err=r1.stats.get('error'))
    s1 = r1.s; rec['s_screen'] = s1
    if s1 > smax + 1e-3:
        return dict(rec, st='screen-discard', sec=time.time() - t0)
    if policy.grid_check and r1.p.grid_lines(k):
        return dict(rec, st='screen-grid', sec=time.time() - t0)
    hit = known.near(s1, policy.same)
    if hit is not None:
        return dict(rec, st='return', s=hit, sec=time.time() - t0)
    if policy.pol_above is not None and s1 > best + policy.pol_above:
        return dict(rec, st='new?', s=s1, unpolished=True, p=r1.p, sec=time.time() - t0)
    t1 = time.time()
    cur, its, chunks = r1.p, 0, 0
    if policy.stage_pit > 0:
        prev = s1
        while True:
            r2 = W.quench(cur, policy.polish.but(pit=policy.stage_pit, flip_top=0))
            if not r2.ok:
                return dict(rec, st='fail', sec=time.time() - t0, err='polish')
            cur = r2.p; chunks += 1; its += r2.stats.get('slp_it', 0)
            hit = known.near(cur.s, policy.late_same)
            if hit is not None and prev - cur.s < policy.late_slope:      # stalled next to a known basin
                rec.update(t_polish=time.time() - t1, chunks=chunks, slp_it=its)
                return dict(rec, st='late-return', s=hit, s_stage=cur.s, sec=time.time() - t0)
            if r2.stats.get('slp_it', 0) < policy.stage_pit or chunks > 60:
                break
            prev = cur.s
        r2 = W.quench(cur, policy.polish)          # finish: default polish incl. flip search
    else:
        r2 = W.quench(cur, policy.polish)
    if not r2.ok:
        return dict(rec, st='fail', sec=time.time() - t0, err='polish')
    its += r2.stats.get('slp_it', 0)
    rec.update(t_polish=time.time() - t1, chunks=chunks, slp_it=its, fp=r2.fp)
    lines = bool(r2.p.grid_lines(k)) if policy.grid_check else False
    return dict(rec, st='new?', s=r2.s, lines=lines, p=r2.p, sec=time.time() - t0)
