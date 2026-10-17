#!/usr/bin/env python3
"""Waste model (10-16, Evan's idea): does a smooth model of record waste w = s^2 - n vs side s point at the n whose
records later improve?  Time-travel backtest on the event history (pk.eventbench).  Pre-registered: PACKER.md WM1.

  python3 analysis/wastemodel.py backtest [--now ISO]
  python3 analysis/wastemodel.py rank [--T ISO] [--top 40]      current ranking (targets)"""
import argparse, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path[:0] = [PK]
import numpy as np
from pk import eventbench as EB

CUTS = ['reg1005', '1008', '1009', '1011', '1012', '2026-10-14T00:00']
NOW = '2026-10-16T21:00'


def frontier(T, nmax=324):
    F = {}
    for n in range(5, nmax + 1):
        s = EB.frontier_at(n, T)
        if abs(s - round(s)) > 1e-9:
            F[n] = s
    return F


def scores(F, window=1.0):
    n = np.array(sorted(F)); s = np.array([F[k] for k in n]); w = s * s - n
    b, c = np.polyfit(np.log(s), np.log(w), 1)
    M0 = np.log(w) - (b * np.log(s) + c)
    ws = w / s
    M1 = np.array([ws[i] - np.median(ws[(np.abs(s - s[i]) <= window) & (np.arange(len(n)) != i)])
                   if ((np.abs(s - s[i]) <= window).sum() > 3) else 0.0 for i in range(len(n))])
    frac = s % 1
    Fall = {k: F.get(k, None) for k in range(5, max(F) + 2)}
    flat = np.array([1.0 if Fall.get(k + 1) is not None and abs(F[k] - Fall[k + 1]) < 1e-9 else 0.0 for k in n])
    L = []
    for k in n:
        m = k - 2 * math.floor(F[k] - 1) - 1
        L.append(1.0 if m in F and abs(F[m] + 1 - F[k]) < 1e-9 else 0.0)
    return n, dict(M0=M0, M1=M1, frac=frac, flat=flat, L=np.array(L)), (b, c)


def auc(score, y):
    y = np.asarray(y, bool)
    if y.all() or (~y).all():
        return float('nan')
    from scipy.stats import rankdata
    r = rankdata(score)
    return (r[y].sum() - y.sum() * (y.sum() + 1) / 2) / (y.sum() * (~y).sum())


def prec(score, y, k):
    o = np.argsort(-score, kind='stable')[:k]
    return float(np.mean(np.asarray(y)[o]))


def backtest(a):
    ev = EB.events()
    print(f'labels to {a.now}')
    agg = {}
    for T in CUTS:
        F = frontier(T); n, S, (b, c) = scores(F)
        Fn = {k: EB.frontier_at(k, a.now) for k in n}
        t = EB.cutoff(T)
        y_imp = np.array([Fn[k] < F[k] - 1e-8 for k in n])
        y_own = np.array([any(e['n'] == k and e['t'] >= t and e['new'] < F[k] - 1e-8 for e in ev) for k in n])
        for lab, y in (('improved', y_imp), ('own-event', y_own)):
            line = f'{T:17s} {lab:9s} n={len(n)} base {y.mean():.2f} | '
            for name, sc in S.items():
                A = auc(sc, y); agg.setdefault((lab, name), []).append(A)
                line += f'{name} AUC {A:.2f} p15 {prec(sc, y, 15):.2f}  '
            print(line)
        for lab, sl in (('n<150', n < 150), ('n>=150', n >= 150)):
            print(f'{"":17s} {"improved":9s} {lab:6s} base {y_imp[sl].mean():.2f} | M1 AUC {auc(S["M1"][sl], y_imp[sl]):.2f}  '
                  f'M0 AUC {auc(S["M0"][sl], y_imp[sl]):.2f}  L AUC {auc(S["L"][sl], y_imp[sl]):.2f}')
    print('\nmean AUC over cutoffs:')
    for (lab, name), v in sorted(agg.items()):
        v = np.array(v, float); print(f'  {lab:9s} {name:5s} mean {np.nanmean(v):.3f}  >= 0.6 at {(v >= 0.6).sum()}/{len(v)}  {np.round(v, 2)}')


def rank(a):
    F = frontier(a.T); n, S, (b, c) = scores(F)
    print(f'log w = {b:.3f} log s + {c:.3f}')
    o = np.argsort(-S['M1'])
    for i in o[:a.top]:
        k = n[i]; s = F[k]
        print(f'n={k:3d} s={s:.6f} frac={s % 1:.3f}  M1 {S["M1"][i]:+.3f} (dS ~ {S["M1"][i] / 2:+.4f})  M0 {S["M0"][i]:+.3f}'
              f'  {"L" if S["L"][i] else " "}{"F" if S["flat"][i] else " "}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    b = sp.add_parser('backtest'); b.add_argument('--now', default=NOW)
    r = sp.add_parser('rank'); r.add_argument('--T', default=NOW); r.add_argument('--top', type=int, default=40)
    a = ap.parse_args()
    {'backtest': backtest, 'rank': rank}[a.cmd](a)
