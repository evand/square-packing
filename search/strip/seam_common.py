import numpy as np, json, glob
from rectpack import half
def load_top(sigwant, H):
    rows = [r for f in glob.glob(f'census4*_H{H}_s*.json') for r in json.load(open(f)) if r['sig'] == sigwant]
    r = min(rows, key=lambda r: r['W']); z = np.array(r['z']); W = r['W']
    x, y, t = z[:4], z[4:8], z[8:]
    for fx in (0, 1):
        for fy in (0, 1):
            xx = W - x if fx else x; yy = H - y if fy else y; tt = -t if fx ^ fy else t
            k = np.argmax(yy)
            if W - xx[k] - half(tt[k]) < 2e-3 and H - yy[k] - half(tt[k]) < 2e-3:
                return xx[k], yy[k], tt[k], W
    raise SystemExit("no orientation")
def load_strip(sigwant, H):
    rows = [r for f in glob.glob(f'census4*_H{H}_s*.json') for r in json.load(open(f)) if r['sig'] == sigwant]
    r = min(rows, key=lambda r: r['W']); z = np.array(r['z']); W = r['W']
    x, y, t = z[:4], z[4:8], z[8:]
    for fx in (0, 1):
        for fy in (0, 1):
            xx = W - x if fx else x; yy = H - y if fy else y; tt = -t if fx ^ fy else t
            k = np.argmax(yy)
            if W - xx[k] - half(tt[k]) < 2e-3 and H - yy[k] - half(tt[k]) < 2e-3:
                return xx, yy, tt, W
def rect_pen(x, y, t, cx, cy, hx, hy):
    dx, dy = x - cx, y - cy; best = np.full(len(x), np.inf)
    for ph in (np.zeros(len(x)), np.full(len(x), np.pi/2), t, t + np.pi/2):
        u = np.abs(np.cos(ph)*dx + np.sin(ph)*dy)
        r = half(t - ph) + hx*np.abs(np.cos(ph)) + hy*np.abs(np.sin(ph))
        best = np.minimum(best, r - u)
    return best
