from wz_sub import *
from wz_sub3 import Z
run('Z: m', Z, lambda g: g['m'])
run('Z: r or mp', Z, lambda g: m_(g['r'], g['mp']))
run('Z: s or mp', Z, lambda g: m_(g['s'], g['mp']))
