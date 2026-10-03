"""Re-run the checker's own primitive on sampled dumped leaves (record consistency, not independence):
EXACT* -> Exact.certify(box); PIECE -> piece_bound(box, bin_data) >= 1.  Also a power test: certify(box, tau = sampled
exact min + 1e-9) must FAIL for the EXACT leaves whose sampled min is attained (the bound is not vacuous)."""
import sys, json, random, time
sys.path.insert(0, '/home/evand/math/square-packing/public/s12/search')
from fractions import Fraction as F
import qx2_zm as QZ
import zm_mixed as ZM
import mixed_cover as MC
import zeromargin as zm

BOX = '/home/evand/math/square-packing/public/s12/runs/qx2_k4x_k008/sol_exact_box9.txt'
cov = ZM.Cover(MC.load(BOX))
chk = QZ.QXChecker(cov, exact_umax=F(1, 2), exact_from=3, use_exact=True, max_depth=18, use_chain=False,
                   cert_mode=False, dump=False)
rows = [json.loads(l) for l in open(sys.argv[1])]
per = int(sys.argv[2]); tcap = float(sys.argv[3])
rng = random.Random(11)
for k in (sys.argv[4].split(',') if len(sys.argv) > 4 else ('EXACT', 'EXACT0', 'EXACT45', 'PIECE')):
    R = [r for r in rows if r['kind'] == k]
    slow = [r for r in R if r['cpu'] >= 100]
    pick = rng.sample(slow, min(len(slow), per // 2)); pick += rng.sample([r for r in R if r not in pick], min(len(R), per - len(pick)))
    ok = 0; bad = []; t0 = time.time(); power = [0, 0]
    for r in pick:
        if time.time() - t0 > tcap: break
        b = tuple(F(v) for v in r['box'])
        if k == 'PIECE':
            L, thr = chk.piece_bound(b, zm.bin_data(b[4], b[5]))
            good = L >= 1
            if good and L > F(r['min']): bad.append(('PIECE bound above sampled mass!', r['box'], str(L), r['min']))
        else:
            good, why = chk.exact.certify(b)
            # power: tau slightly above the sampled minimum (only meaningful when the minimiser has u > u0 or u0 > 0)
            smin = F(r['min'])
            at_u = F(r['at'][2])
            if good and power[0] + power[1] < 40 and at_u > 0 and at_u * at_u + 2 * at_u - 1 < 0:
                g2, _ = chk.exact.certify(b, tau=smin + F(1, 10 ** 9))
                power[0 if not g2 else 1] += 1
                if g2: bad.append(('certified ABOVE a sampled exact mass', r['box'], r['min'], r['at']))
        ok += good
        if not good: bad.append(('not re-certified', r['box']))
    print(f'{k}: re-checked {ok}/{len(pick)} in {time.time()-t0:.0f}s; power test (tau = sampled min + 1e-9 fails / '
          f'certifies): {power}; problems {len(bad)}', flush=True)
    for x in bad[:8]: print('   ', x)
