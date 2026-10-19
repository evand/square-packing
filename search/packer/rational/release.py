#!/usr/bin/env python3
"""Release the fixed angles of rational designs (10-18): fq deep polish (free angles, SLP + corner-corner branch
search) from each kept design; report s_rel vs the rational s* and the frontier; import released packings with
parent = the design's store id, tags rational-release + the design's tri tag.

  release.py DIR [--procs 15] [--gap 0.03]       reads DIR/keep.jsonl, writes DIR/release.jsonl
"""
import sys, os, json, argparse, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def job(r):
    from pk.engine import fq_quench, POLISH_DEEP
    from pk.packing import read
    with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as tmp:
        run = fq_quench((r['s'], np.array(r['sq'])), POLISH_DEEP, tmp=tmp, timeout=900)
        if run.out is None:
            return r, None
        p = read(run.out)
        return r, (float(p.s), p.sq.tolist(), float(p.max_pen()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('dir'); ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--gap', type=float, default=0.03)
    a = ap.parse_args()
    from multiprocessing import Pool
    from pk.store import Store
    from pk.packing import Packing
    import ratgen
    fr = ratgen.frontier()
    rows = [json.loads(l) for l in open(os.path.join(a.dir, 'keep.jsonl'))]
    rows = [r for r in rows if r['s'] < fr.get(r['n'], 1e9) + a.gap]
    print(len(rows), 'designs', flush=True)
    st = Store()
    with Pool(a.procs) as pool, open(os.path.join(a.dir, 'release.jsonl'), 'w') as fo:
        for r, out in pool.imap_unordered(job, rows):
            if out is None:
                continue
            s, sq, pen = out
            f = fr.get(r['n'])
            rec = dict(n=r['n'], s_rat=r['s'], s_rel=s, front=f, pen=pen, gain=r['s'] - s, d=s - f,
                       tri=r.get('tri'), C=r['C'], L=r['L'], shear=r.get('shear'))
            if pen <= 1e-9:
                pd, _ = st.add(Packing(s=r['s'], sq=np.array(r['sq'])), source='ours:rational', finder='ours',
                               kind='rational-design', tags=['rational'], commit=False)
                tri = r.get('tri'); tags = ['rational-release'] + ([f'tri:{tri[0]}-{tri[1]}-{tri[2]}'] if tri else [])
                pid, _ = st.add(Packing(s=s, sq=np.array(sq)), source='ours:rational-release', finder='ours',
                                kind='rational-release', status='polished', parent=pd, tags=tags, commit=False)
                rec['pid'] = pid
            fo.write(json.dumps(rec) + '\n'); fo.flush()
            flag = ' <<< BELOW' if rec['d'] < -1e-9 else ''
            print(f"{r['n']} rat {r['s']:.10f} rel {s:.10f} front {f:.10f} d {rec['d']:+.2e}{flag}", flush=True)
    st.db.commit()


if __name__ == '__main__':
    main()
