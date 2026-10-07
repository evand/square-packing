"""Tables for RESULTS.md from runs/classes.jsonl (+ runs/intensify.jsonl if present: best margin per (h, class) is the max)."""
import json, os, sys
from collections import defaultdict
from occ import label, strips, canon

HERE = os.path.dirname(os.path.abspath(__file__))
ZERO = -1e-9


def load():
    best = {}
    for line in open(os.path.join(HERE, "runs/classes.jsonl")):
        r = json.loads(line)
        key = (r["h"], tuple(r["class"]))
        best[key] = dict(r, extra=0)
    f = os.path.join(HERE, "runs/intensify.jsonl")
    if os.path.exists(f):
        cnt = defaultdict(int); top = {}
        for line in open(f):
            r = json.loads(line)
            key = (r["h"], tuple(r["class"]))
            cnt[key] += 1
            if key not in top or r["margin"] > top[key]["margin"]:
                top[key] = r
        for key, r in top.items():
            b = best[key]
            b["extra"] = cnt[key]
            if r["margin"] > b["margin"]:
                b["margin"] = r["margin"]; b["witness"] = f"runs/witnesses/h{key[0]}_" + "".join(map(str, key[1])) + "_int.json"; b["maxtilt"] = r["maxtilt"]
    return best


if __name__ == "__main__":
    B = load()
    hs = sorted({h for h, _ in B})
    mode = sys.argv[1] if len(sys.argv) > 1 else "kz"
    if mode == "kz":
        print("| k | z | #classes | " + " | ".join(f"h={h}: best (#at 0)" for h in hs) + " |")
        print("|---|---|---|" + "---|" * len(hs))
        groups = sorted({(sum(c[:4]), c[8]) for _, c in B}, key=lambda t: (-t[1] + t[0], t))
        for k, z in sorted(groups, key=lambda t: (t[1] - t[0], t[0]), reverse=True):
            row = []
            ncl = None
            for h in hs:
                rs = [r for (hh, c), r in B.items() if hh == h and sum(c[:4]) == k and c[8] == z]
                ncl = len(rs)
                m = max(r["margin"] for r in rs)
                n0 = sum(r["margin"] >= ZERO for r in rs)
                row.append(f"{m:+.4f} ({n0})")
            print(f"| {k} | {z} | {ncl} | " + " | ".join(row) + " |")
    elif mode == "classes":
        zmin = int(sys.argv[2]); zmax = int(sys.argv[3])
        cls = sorted({c for _, c in B if zmin <= c[8] <= zmax}, key=lambda c: (-c[8], -sum(c[:4]), c))
        print("| class | " + " | ".join(f"h={h}" for h in hs) + " | restarts (h=1.15) | witness (h=1.15) |")
        print("|---|" + "---|" * len(hs) + "---|---|")
        for c in cls:
            row = []
            for h in hs:
                r = B[(h, c)]
                row.append(f"{r['margin']:+.4f}" + (f" ({r['maxtilt']:.1f}°)" if r["margin"] >= ZERO else ""))
            r = B[(1.15, c)]
            print(f"| `{label(c)}` | " + " | ".join(row) + f" | {r['restarts'] + r['extra']} | `{r['witness']}` |")
