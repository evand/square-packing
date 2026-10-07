import sys
from fractions import Fraction as F
import test_poly as TP
from test_poly import ZM, zm, MC, S
import os
cv = MC.load(os.path.join(S,'qx2_data/L4_k02_box7.txt')); cov = ZM.Cover(cv)
chk = ZM.MixedChecker(cov, cert_mode=False)
chk0 = ZM.MixedChecker(ZM.Cover(dict(cv, polygons=[])), cert_mode=True)
a,b = F(9,5), F(26,5)
for box in [tuple(F(x) for x in ['4153/2500', '4153/2500', '67/50', '13437/10000', '783/2500', '1261/4000']),
            tuple(F(x) for x in ['3263/2500', '6751/5000', '2199/1250', '8991/5000', '4239/20000', '4819/20000'])]:
    B = zm.bin_data(box[4], box[5])
    L,_ = chk.piece_bound(box,B); L0,_ = chk0.piece_bound(box,B)
    K,_ = TP.region_K(chk, box, B)
    print('polyL', L-L0, float(L-L0), 'K', [(float(x),float(y)) for x,y in K])
    for u in [box[4], (box[4]+box[5])/2, box[5]]:
        for cx in (box[0],box[1]):
            for cy in (box[2],box[3]):
                print('  ', float(u), float(cx), float(cy), float(TP.area_QU(cx,cy,u,a,b)), [ (float(x),float(y)) for x,y in TP.clip_rect(TP.square(cx,cy,u),a,b,a,b)])
