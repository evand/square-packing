"""Named fq relax schedules: one table (10-20; merged from relax.OPS and shove/fill_relax.RELAX).

A schedule is a QOpt (pk.engine): the fq quench run on a packing, usually an overlapping start from a fill or a
construction.  Refer to schedules by name, so runs log what they used and results compare across tools.

Two bases, kept as they were in their source tables (same names give byte-identical fq flags):
  * relax.OPS (10-13, compositions / band merges): built on QOpt() (fq defaults: full ALM + polish).
  * fill_relax.RELAX (10-16, CS1 agent B): built on engine.SCREEN (pit 8, no corner-flip search); the bench's
    evaluate polishes afterwards.

What we know (NOTES_relax.md, CS1 dev tests 10-16): from a heavily overlapping fill, plain ALM at any loosen
(1.0-1.2), `shrink` and soft squares collapse to the k + 1 grid; the rounded-square compress pre-stages
(round.15c / round.25c) are the only ones that kept non-grid structure.  The relax step mattered more than the
placement.  soft15c (relax.OPS) is round.15c on the QOpt() base with 40 round steps and a 0.02 kick.

  from pk.schedules import SCHEDULES, get
  opt = get('round.15c').but(seed=3)
"""
from __future__ import annotations
from .engine import QOpt, SCREEN

# relax.OPS (10-13): QOpt() base.  relax.py's CLI default --ops is exactly these, in this order.
OPS_NAMES = ('polish', 'shrink', 'kshrink', 'soft05', 'soft15', 'soft15c', 'loose')
_OPS = {
    'polish': QOpt(loosen=1.0, no_alm=True),                                     # rigid baseline
    'shrink': QOpt(loosen=1.0, shrink=True),                                     # fixed-side overlap relax + shrink hops
    'kshrink': QOpt(loosen=1.0, kick=0.02, shrink=True),
    'soft05': QOpt(loosen=1.0, round_r0=0.05, hold_f=1.0, round_steps=40, round_kick=0.03),
    'soft15': QOpt(loosen=1.0, round_r0=0.15, hold_f=1.0, round_steps=40, round_kick=0.03),
    'soft15c': QOpt(loosen=1.0, round_r0=0.15, round_compress=True, round_steps=40, round_kick=0.02),
    'loose': QOpt(loosen=1.02, kick=0.03),                                       # jitter + ALM from a dilated box
}

# shove/fill_relax.RELAX (10-16): SCREEN base.  fill_relax's default SHOVE_RELAX is round.15c,round.25c,holdP1.01.
RELAX_NAMES = ('L1.1', 'L1.2', 'round.15c', 'round.25c', 'round.15cL1.1', 'holdP1.01', 'shrinkhop')
_RELAX = {
    'L1.1':      dict(loosen=1.1),
    'L1.2':      dict(loosen=1.2),
    'round.15c': dict(loosen=1.0, round_r0=0.15, round_compress=True),
    'round.25c': dict(loosen=1.0, round_r0=0.25, round_compress=True, round_steps=40),
    'round.15cL1.1': dict(loosen=1.1, round_r0=0.15, round_compress=True),
    'holdP1.01': dict(loosen=1.0, round_r0=0.15, hold_f=1.01, no_alm=True, pit=300, flip_top=8),
    'shrinkhop': dict(loosen=1.1, shrink=True, hops=3),
    # 10-24 relax-schedule bench (PACKER "RS1"): round.15c variants.  s5: 5 hardening steps; nc: no compress at r0;
    # h1003: box held 0.3 % above the compressed side; g6: loose inner solves (gtol 1e-6, <= 300 L-BFGS its per step);
    # p: pressure (side free at every step, Nakajima-style); pk: p + kick noise peaking at r = 0.06 (angle lock-in).
    'round.15c.s5':    dict(loosen=1.0, round_r0=0.15, round_compress=True, round_steps=5),
    'round.15c.nc':    dict(loosen=1.0, round_r0=0.15),
    'round.15c.h1003': dict(loosen=1.0, round_r0=0.15, round_compress=True, hold_f=1.003),
    'round.15c.g6':    dict(loosen=1.0, round_r0=0.15, round_compress=True, round_gtol=1e-6, round_maxit=300),
    'round.15c.g6s5':  dict(loosen=1.0, round_r0=0.15, round_compress=True, round_gtol=1e-6, round_maxit=300, round_steps=5),
    'round.15c.p':     dict(loosen=1.0, round_r0=0.15, round_compress=True, round_press=True),
    'round.15c.p5':    dict(loosen=1.0, round_r0=0.15, round_compress=True, round_press=True, round_steps=5),
    'round.15c.pk':    dict(loosen=1.0, round_r0=0.15, round_compress=True, round_press=True, round_kick=0.01,
                            round_kick_peak=0.06, round_kick_w=0.03),
    # 10-24 RS2: press steps capped at K outer ALM iterations (fq --round-press-outer; uncapped they always run 60)
    'round.15c.p5o8':  dict(loosen=1.0, round_r0=0.15, round_compress=True, round_press=True, round_steps=5,
                            round_press_outer=8),
    'round.15c.p10o4': dict(loosen=1.0, round_r0=0.15, round_compress=True, round_press=True, round_steps=10,
                            round_press_outer=4),
    # 10-24 random starts: disks (r0 = 0.5) compressed, then hardened under pressure to squares (Nakajima's path in fq)
    'disk.p20':  dict(loosen=1.0, round_r0=0.5, round_compress=True, round_press=True, round_steps=20),
    'disk.p40':  dict(loosen=1.0, round_r0=0.5, round_compress=True, round_press=True, round_steps=40),
    'disk.p20k': dict(loosen=1.0, round_r0=0.5, round_compress=True, round_press=True, round_steps=20, round_kick=0.01,
                      round_kick_peak=0.06, round_kick_w=0.03),
    # plain screens at one loosen each (the explorer / shove harness evaluate takes the best of 1.0 / 1.02 / 1.05)
    'screen1.0':  dict(loosen=1.0),
    'screen1.02': dict(loosen=1.02),
    'screen1.05': dict(loosen=1.05),
}

SCHEDULES: dict[str, QOpt] = dict(_OPS)
SCHEDULES.update({k: SCREEN.but(**v) for k, v in _RELAX.items()})
SOURCE = {**{k: 'relax.OPS (QOpt base)' for k in _OPS}, **{k: 'fill_relax.RELAX (SCREEN base)' for k in _RELAX}}


def get(name: str) -> QOpt:
    try:
        return SCHEDULES[name]
    except KeyError:
        raise KeyError(f'unknown relax schedule {name!r}; known: {", ".join(SCHEDULES)}') from None


if __name__ == '__main__':
    for k, v in SCHEDULES.items():
        print(f'{k:14s} {SOURCE[k]:32s} {" ".join(v.args())}')
