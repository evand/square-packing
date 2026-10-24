"""Atlas triangle of the store's best packings for n <= NMAX (Kaggle range, 10-22), in jlevy's layout and colours
(reuses site/timelapse/render.py: row k = ceil(sqrt n) right-aligned, n = k^2 on the right edge; hue = angle class,
shade = full-side contacts).  Grid cells (best side >= k, or no store packing) are drawn as the k x k grid minus
squares.  Writes OUT (full size) and OUT with '-small' (downscaled to --small px wide).

  python3 triangle1000.py [--nmax 1000] [--cell 96] [--out runs/kaggle/triangle1000.png] [--procs 15]
"""
import argparse, math, os, sqlite3, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'site', 'timelapse'))
import render                                                           # noqa: E402  (tile, font, colours)

STORE = os.path.join(HERE, 'runs/store.sqlite')
FEAS = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9)"
K = lambda n: math.isqrt(n - 1) + 1


def grid_view(n):
    k = K(n)
    cells = [(i + .5, j + .5, 0.0) for j in range(k) for i in range(k)][:n]
    return dict(s=float(k), squares=cells)


def view(n):
    db = sqlite3.connect(f'file:{STORE}?mode=ro', uri=True)
    r = db.execute(f"SELECT s, coords FROM packing WHERE n=? AND {FEAS} ORDER BY s LIMIT 1", (n,)).fetchone()
    if r is None or r[0] >= K(n) - 1e-9:
        return n, grid_view(n), True
    sq = np.frombuffer(r[1], '<f8').reshape(-1, 3)
    return n, dict(s=r[0], squares=[tuple(map(float, q)) for q in sq]), False


def tile_job(a):
    n, cell = a
    n, v, g = view(n)
    return n, render.tile(v, cell).tobytes(), g, v['s']


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--nmax', type=int, default=1000)
    ap.add_argument('--cell', type=int, default=96); ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--out', default=os.path.join(HERE, 'runs/kaggle/triangle1000.png'))
    ap.add_argument('--small', type=int, default=2400)
    A = ap.parse_args()
    rows = K(A.nmax); cell = A.cell; pitch = int(cell * 1.12); rpitch = int(cell * 1.22)
    gapf = 0.45; mr, top, lab = 30, 30, max(9, cell // 7)
    W = int(mr * 2 + (2 * rows - 1 + gapf) * pitch); H = top + rows * rpitch + 20
    res = {}
    with ProcessPoolExecutor(A.procs) as ex:
        for n, b, g, s in ex.map(tile_job, [(n, cell) for n in range(1, A.nmax + 1)], chunksize=4):
            res[n] = (Image.frombytes('RGB', (cell, cell), b), g, s)
    im = Image.new('RGB', (W, H), render.BG); d = ImageDraw.Draw(im); f = render.font(lab)
    for k in range(1, rows + 1):
        ns = [n for n in range((k - 1) ** 2 + 1, k * k + 1) if n <= A.nmax]
        g0 = next((m for m in ns if all(res[q][1] for q in ns if q >= m)), k * k + 1)
        for n in ns:
            c = k * k - n
            x = int(W - mr - cell - c * pitch - (gapf * pitch if n < g0 else 0))
            y = top + (k - 1) * rpitch
            im.paste(res[n][0], (x, y))
            d.text((x, y + cell + 1), f'{n}', fill=render.MUTED, font=f)
    im.save(A.out)
    sm = im.resize((A.small, int(H * A.small / W)), Image.LANCZOS)
    sm.save(A.out.replace('.png', '-small.png'))
    print(A.out, im.size)


if __name__ == '__main__':
    main()
