#!/usr/bin/env python3
"""Kernel-cost model for the ZMTree / ZMTreeM Lean builds, its fit to the measured builds, and the stratified
extrapolation of a Lean-tree census sample (search/golf/lean_census.py) to the whole D4 region.

  python3 search/golf/costmodel.py fit
  python3 search/golf/costmodel.py estimate search/golf/data/lc_current.jsonl [more.jsonl ...]

Model (kernel CPU of the part files; `Pts`/`Cov`/top file are minutes and ignored):

  CPU = c_Z * (Z leaves) + c_adm * (ADM witness entries) + c_ch * (chain entries) + c_L * (L-lines)
        + c_T * (T-groups) + c_S * (S-blocks)

E leaves, C nodes and splits are ~free (a few Nat comparisons each).  c_adm, c_ch come from the s(13) heaviest-chunk
profile (LADDER.md: ADM entry ~0.5 ms, chain entry ~2.5 ms), c_L, c_T, c_S from the one-leaf T1 profile
(notes/lean-segments.md sec. 4: S-only leaf 23 ms, T-group leaf 27 ms, L-block leaf 125 ms ~ 40 ms per line);
c_Z (decode, walls, candidate walk, per-leaf context, region counts, file import amortised) is fitted on the three
full builds.  Memory: RSS ~ 6.6 GB (Mathlib) + m_d * (base-2^20 digits in the file); m_d fitted on S13 / S32Z.
"""
import json
import sys
from collections import defaultdict

# ---- the measured builds (LADDER.md, notes/lean-segments.md) ----
BUILDS = [
    # name, Z leaves, E leaves, ADM entries, chain entries, L-lines, T-groups, S-blocks, kernel CPU s, digits, files, peak RSS GB
    dict(name='S13 (s(13)=4)', Z=4079, E=1738, adm=1122165, ch=203367, L=0, T=0, S=0, cpu=1995,
         digits=1.80e6, files=4, rss=13.2),
    dict(name='S32Z (s(32)=6)', Z=89473, E=5127, adm=43.8e6 - 3.32e6 * 0, ch=3.32e6, L=0, T=0, S=0, cpu=49766,
         digits=55.1e6, files=96, rss=15.2),
    # T1 grid cover (s16_ge_4): 9,451 Z leaves, 8,256 with an L-block (~3 lines each: the one-leaf profile's average),
    # ~1,742 with T-groups; 6 files x ~420 s.  No points.
    dict(name='S16/T1 (segments only)', Z=9451, E=0, adm=0, ch=0, L=8256 * 3, T=1742, S=9451 - 8256, cpu=2520,
         digits=None, files=6, rss=11.0),
]
# In LADDER.md s(32)'s 47.1 M claimed entries = 43.8 M ADM witnesses + 3.32 M chain entries (+ 358 k pivots).

MICRO = dict(adm=0.5e-3, ch=2.5e-3, L=40e-3, T=4e-3, S=0.0)

# three parameter sets (low / central / high), see GOLF_PILOT.md sec. 1
MODELS = {
    'low':     dict(Z=0.15, adm=0.5e-3, ch=2.5e-3, L=40e-3, T=4e-3, S=0.0),
    'central': dict(Z=0.22, adm=0.5e-3, ch=2.5e-3, L=40e-3, T=4e-3, S=0.0),
    # calibrated on the s(21) pilot one-leaf chunks (0.69 s per T+230-point leaf, 0.89 s per L+points leaf):
    'high':    dict(Z=0.46, adm=1.0e-3, ch=2.5e-3, L=60e-3, T=4e-3, S=0.0),
}
RSS_BASE = 6.6
RSS_PER_DIGIT = 15.0e-6   # GB per digit in one file (fit below)
RSS_CAP = 15.0            # target peak per process (GB) -> digits per file


def cpu(m, f):
    return (m['Z'] * f['Z'] + m['adm'] * f['adm'] + m['ch'] * f['ch'] + m['L'] * f['L'] + m['T'] * f['T']
            + m['S'] * f['S'])


def fit():
    print("c_Z fitted per build (with the micro-profile entry costs fixed):")
    for b in BUILDS:
        rest = MICRO['adm'] * b['adm'] + MICRO['ch'] * b['ch'] + MICRO['L'] * b['L'] + MICRO['T'] * b['T']
        cz = (b['cpu'] - rest) / b['Z']
        print(f"  {b['name']:28s} CPU {b['cpu']:7.0f} s, entries+lines {rest:7.0f} s -> c_Z = {cz:.3f} s/leaf")
    for k, m in MODELS.items():
        print(f"model {k}: " + ", ".join(f"{b['name'].split()[0]} pred {cpu(m, b):.0f} / meas {b['cpu']} s"
                                         f" ({cpu(m, b) / b['cpu']:.2f}x)" for b in BUILDS))
    print("memory: RSS per file = 6.6 GB + m_d * digits/file")
    for b in BUILDS:
        if b['digits']:
            dpf = b['digits'] / b['files']
            print(f"  {b['name']:28s} {dpf / 1e3:.0f} k digits/file, {b['rss']} GB -> m_d = "
                  f"{(b['rss'] - RSS_BASE) / dpf * 1e6:.1f} KB/digit")


def feats(r):
    return dict(Z=r.get('Z', 0), E=r.get('E', 0), adm=r.get('adm_ents', 0), ch=r.get('chain_ents', 0),
                L=r.get('Llines', 0), T=r.get('Tgrp', 0), S=r.get('Sblk', 0), digits=r.get('digits'),
                boxes=r.get('boxes', 0), cpu=r.get('cpu', 0), UNCERT=r.get('UNCERT', 0),
                ptclaims=r.get('ptclaims_sum', 0), pairs=r.get('pairs', 0), segclaims=r.get('segclaims', 0))


if __name__ == '__main__':
    fit()   # the extrapolation of a census sample is search/golf/estimate.py
