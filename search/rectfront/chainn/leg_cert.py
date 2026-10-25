"""Certificates for the leg_gen sub-cases (same LP as band_cert: Handelman products + circle multiples).
Cache: runs/rectfront/band/leg_certs.pkl (key: tag)."""
import os, sys, pickle, time
import sympy as sp
import band_cert as BC
import leg_gen as LG
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../runs/rectfront/band/leg_certs.pkl')
def gens_for(sc):
    gens = BC.base_gens()
    for n, (kind, e) in enumerate(sc['hyps']): gens.append((f'h{n}', sp.expand(e)))
    return gens
if __name__ == '__main__':
    degs = [int(d) for d in sys.argv[1:]] or [2, 3, 4]
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    fails = []
    for sc in LG.subcases():
        if sc['tag'] in cache: continue
        if LG.numcheck(sc, nth=181, nH=3) is None: cache[sc['tag']] = ('empty',); continue
        t0 = time.time(); got = None
        for d in degs:
            gens = gens_for(sc); res = BC.certify(sc['P'], gens, d)
            if res and res[2] < 1e-7:
                got = (d, [(tuple(gens[k][0] for k in lab), w) for lab, w in res[0]], res[1], [gn for gn, _ in gens]); break
        if got: cache[sc['tag']] = got; print(f"  {sc['tag']}: deg {got[0]}, {len(got[1])} prods ({time.time()-t0:.1f}s)", flush=True)
        else: fails.append(sc['tag']); print(f"  {sc['tag']}: FAIL ({time.time()-t0:.1f}s)", flush=True)
        pickle.dump(cache, open(CACHE, 'wb'))
    print('failed:', fails)

def run_splits(degs=(4, 5, 6)):
    cache = pickle.load(open(CACHE, 'rb'))
    TIGHT = {"Sa_t_t'_Lz_A1", "Sb_t_t'_A1_Lz", "C0_t_t'_Lz_Lz", "C1_t_t'_A1_A1"}
    for sc in LG.subcases():
        if sc['tag'] in cache or sc['tag'] in TIGHT: continue
        t0 = time.time(); hy = [e for _, e in sc['hyps']]
        got = BC.certify_split(sc['P'], [(None, e) for e in hy], degs)
        print(f"  {sc['tag']}: {'split ' + str(got[0]) if got else 'FAIL'} ({time.time()-t0:.1f}s)", flush=True)
        if got: cache[sc['tag']] = ('split', got); pickle.dump(cache, open(CACHE, 'wb'))
