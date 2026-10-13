"""Rectangle pieces (m squares in (h + d) x h) <-> the store's `piece` table.

One place pieces live, whatever made them: rectlib crops, rectsearch minima, composer by-products.  Readers here turn the
run-directory formats into store rows; compose2 reads only the store.

  python3 -m pk.pieces import runs/rectlib/lib1/pieces.jsonl runs/rectsearch/sw1 ...   (jsonl files or rectsearch dirs)
  python3 -m pk.pieces stats
"""
from __future__ import annotations
import glob, json, os, sys
from .packing import Packing
from .store import Store


def from_jsonl(path):
    """rectlib / rectsearch-export rows -> (Packing, source, ref)."""
    for l in open(path):
        r = json.loads(l)
        if 'h' not in r or not r.get('ok', True):
            continue
        yield Packing(r['h'], r['sq'], dx=float(r['d'])), r['key']['src'], f"{path}:{r['key'].get('piece', 0)}"


def from_rectsearch(d):
    """Every distinct minimum a rectsearch run logged (minima_*.jsonl; best.txt for runs that predate them)."""
    for cell in sorted(glob.glob(os.path.join(d, '**', 'm*_d*'), recursive=True)):
        dd = float(os.path.basename(cell).split('_d')[1])
        files = glob.glob(os.path.join(cell, 'minima_*.jsonl'))
        for mf in files:
            for r in map(json.loads, open(mf)):
                yield Packing(r['h'], r['sq'], dx=dd), 'rectsearch', mf
        if not files and os.path.exists(os.path.join(cell, 'best.txt')):
            L = [l.split() for l in open(os.path.join(cell, 'best.txt')) if l.strip()]
            sq = [[float(t) for t in x[:3]] for x in L[1:1 + int(L[0][0])]]
            yield Packing(float(L[0][1]), sq, dx=dd), 'rectsearch', os.path.join(cell, 'best.txt')


def import_paths(paths, store=None, tags=()):
    S = store or Store()
    new = seen = 0
    for p in paths:
        it = from_jsonl(p) if p.endswith('.jsonl') else from_rectsearch(p)
        for pk, src, ref in it:
            _, c = S.add_piece(pk, src, ref=ref, tags=tags, commit=False)
            new += c; seen += 1
    S.db.commit()
    return seen, new


if __name__ == '__main__':
    if sys.argv[1] == 'import':
        print('pieces seen %d, new %d' % import_paths(sys.argv[2:]))
    else:
        S = Store()
        q = ("SELECT substr(source, 1, instr(source || ':', ':') - 1) AS src, count(*), count(DISTINCT m || ':' || d), "
             "sum(status='infeasible') FROM piece GROUP BY src")
        print('source pieces cells infeasible')
        for r in S.db.execute(q):
            print(*r)
