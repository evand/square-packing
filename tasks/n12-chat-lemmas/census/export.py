"""Write witness JSONs for (h, class) bests that came from runs/intensify.jsonl, and the equality-class plateau
points with the largest touching tilt.  Every witness is re-verified here (SAT margin, class membership)."""
import json, os
import numpy as np
from opt import sat_margin
from occ import classify, canon, label, regions, assignment
from rattle import active_tilts

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.join(HERE, "runs/witnesses")


def in_boxes(P, h, c):
    boxes = [regions(h)[r] for r in assignment(c)]
    return all(b[0] - 1e-12 <= x <= b[1] + 1e-12 and b[2] - 1e-12 <= y <= b[3] + 1e-12 for (x, y, _), b in zip(P, boxes))


top, plateau = {}, {}
for line in open(os.path.join(HERE, "runs/intensify.jsonl")):
    r = json.loads(line)
    key = (r["h"], tuple(r["class"]))
    if key not in top or r["margin"] > top[key]["margin"]:
        top[key] = r
    if r["class"] == [1, 1, 1, 1, 1, 1, 1, 1, 4] and r["margin"] >= -1e-9:
        t, a, _ = active_tilts(r["poses"])
        r["active_maxtilt"] = float(t[a].max())
        if r["h"] not in plateau or r["active_maxtilt"] > plateau[r["h"]]["active_maxtilt"]:
            plateau[r["h"]] = r

for (h, c), r in top.items():
    P = np.array(r["poses"])
    tag = f"h{h}_" + "".join(map(str, c)) + "_int"
    json.dump({"h": h, "class": list(c), "label": label(c), "margin_recheck": sat_margin(P), "in_boxes": in_boxes(P, h, c),
               "poses_x_y_thetarad": P.tolist()}, open(os.path.join(W, tag + ".json"), "w"), indent=1)
for h, r in plateau.items():
    P = np.array(r["poses"])
    tag = f"h{h}_111111114_plateau_maxtilt"
    json.dump({"h": h, "class": r["class"], "margin_recheck": sat_margin(P), "in_boxes": in_boxes(P, h, tuple(r["class"])),
               "touching_maxtilt_deg": r["active_maxtilt"], "poses_x_y_thetarad": P.tolist()},
              open(os.path.join(W, tag + ".json"), "w"), indent=1)
    print(h, "plateau touching max tilt %.2f" % r["active_maxtilt"], "margin %.2e" % sat_margin(P))

# recheck all witnesses
bad = 0
for f in os.listdir(W):
    d = json.load(open(os.path.join(W, f)))
    P = np.array(d["poses_x_y_thetarad"])
    m = sat_margin(P)
    ref = d.get("margin", d.get("margin_recheck"))
    if abs(m - ref) > 1e-9 or not in_boxes(P, d["h"], tuple(d["class"])):
        bad += 1; print("MISMATCH", f, m, ref)
print(len(os.listdir(W)), "witnesses, mismatches:", bad)
