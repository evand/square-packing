#!/usr/bin/env python3
"""Plot a certificate's points (and list the distinct coordinate values).

    python3 search/uniform/plot_cert.py CERT [CERT ...] [-o OUT.png]
"""
import argparse, sys
from fractions import Fraction
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_cert(path):
    t = Path(path).read_text().split()
    sn, sd, D, W, m = map(int, t[:5])
    rows = [(int(t[5+3*i]), int(t[6+3*i]), int(t[7+3*i])) for i in range(m)]
    return Fraction(sn, sd), D, W, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("certs", nargs="+"); ap.add_argument("-o", default=None)
    a = ap.parse_args()
    n = len(a.certs)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.2), squeeze=False)
    for ax, c in zip(axes[0], a.certs):
        s, D, W, rows = read_cert(c)
        P = np.array([[x / D, y / D] for x, y, _ in rows]); w = np.array([r[2] for r in rows])
        sf = float(s)
        ax.add_patch(plt.Rectangle((0, 0), sf, sf, fill=False, lw=1.2))
        ax.add_patch(plt.Rectangle((0.5 * sf - 0.5, 0.5 * sf - 0.5), 1, 1, fill=False, ls=":", lw=0.8, color="gray"))
        ax.scatter(P[:, 0], P[:, 1], s=14 + 10 * (w - 1), c="C3", zorder=3)
        ax.set_xlim(-0.1, sf + 0.1); ax.set_ylim(-0.1, sf + 0.1); ax.set_aspect("equal")
        ax.set_title(f"{Path(c).name}\n{len(rows)} pts, k={W}, s={sf:.4f}", fontsize=9)
        ax.set_xticks(np.arange(0, sf + 0.01, 0.5)); ax.set_yticks(np.arange(0, sf + 0.01, 0.5)); ax.grid(alpha=0.25)
        xs = sorted(set(P[:, 0].round(4))); ys = sorted(set(P[:, 1].round(4)))
        print(f"{c}: m={len(rows)} k={W} s={sf:.6f}")
        print("  distinct x:", " ".join(f"{v:.4f}" for v in xs))
        print("  distinct y:", " ".join(f"{v:.4f}" for v in ys))
        cx = sf / 2
        r = np.sort(np.hypot(P[:, 0] - cx, P[:, 1] - cx).round(3))
        print("  radii from centre:", " ".join(f"{v:.3f}" for v in r))
    out = a.o or (Path(a.certs[0]).stem + ".png")
    fig.tight_layout(); fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
