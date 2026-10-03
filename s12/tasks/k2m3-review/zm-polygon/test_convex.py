import os, sys
S = os.path.expanduser('~/math/square-packing/public/s12/search'); sys.path.insert(0, S)
import zm_mixed as ZM
cases = {
 'ccw square': [(0,0),(2,0),(2,2),(0,2)],
 'cw square': [(0,0),(0,2),(2,2),(2,0)],
 'ccw square + collinear': [(0,0),(1,0),(2,0),(2,2),(0,2)],
 'double wound pentagon': [(2,0),(1,2),(0,1),(2,1),(0,0)][::1],
 'pentagram': [(10,0),(3,9),(-8,6),(-8,-6),(3,-9)],
 'pentagram order': [(10,0),(-8,6),(3,-9),(3,9),(-8,-6)],
 'repeated vertex': [(0,0),(2,0),(2,2),(2,2),(0,2)],
 'collinear only': [(0,0),(1,0),(2,0)],
 'spike back': [(0,0),(2,0),(1,0),(1,1)],
}
for k, v in cases.items():
    try: ZM.check_simple_convex(v); r = 'accepted'
    except ValueError as e: r = 'rejected: ' + str(e)
    print(f'{k:28s} {r}')
