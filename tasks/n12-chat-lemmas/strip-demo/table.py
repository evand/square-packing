"""h vs certificate size: search, then check exactly.  Writes certs/h<h>_<menu>.json and prints rows.

    python3 table.py MENU H [H ...]        MENU in {three, fine}
"""
import sys, time
from fractions import Fraction as Fr
import search, check

MENUS = {'three': [Fr(3, 4), Fr(19, 20)], 'fine': None}

if __name__ == '__main__':
    menu_name = sys.argv[1]
    for hs in sys.argv[2:]:
        h = Fr(hs)
        t = time.time()
        S, boxes, Z, v, ok, nraw, dt = search.run(h, Fr(457, 500), Fr(1, 100), 700, 1, True, False, MENUS[menu_name])
        if not ok:
            print(f"{menu_name} h={hs}: search FAILED at {len(boxes)} boxes (float value {v:.5f})", flush=True)
            continue
        path = f"certs/h{float(h):.3f}_{menu_name}.json"
        search.export(S, boxes, Z, path)
        out = check.check(path, quiet=True)
        print(f"{menu_name} h={hs} ({float(h):.4f}): boxes={out['boxes']} (before merge {nraw}), "
              f"H-boxes={out['H_boxes']}, leaf inequalities={out['leaf_inequalities']}, menu size={len(out['menu'])}, "
              f"exact min chain={float(out['min_chain']):.6f}  [{time.time() - t:.0f}s]  -> {path}", flush=True)
