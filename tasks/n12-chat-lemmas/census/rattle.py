"""Tilt of non-rattler squares: a square is 'active' if its own min gap (pairs + walls) is <= tol."""
import numpy as np
from opt import pair_sat, wall_gap

def active_tilts(P, tol=1e-7):
    P = np.asarray(P)
    G = pair_sat(P); w = wall_gap(P)
    own = np.minimum(np.where(np.isfinite(G), G, np.inf).min(1), w)
    tilt = np.degrees(np.minimum(P[:, 2], np.pi / 2 - P[:, 2]))
    act = own <= tol
    return tilt, act, own
