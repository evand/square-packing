"""pk: packing search toolkit (store, engine, moves, replay bench).  See ../DESIGN.md."""
import os, sys
_PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PKG not in sys.path:                     # legacy modules (layout, mcmin, ...) live next to the package
    sys.path.insert(0, _PKG)
from .packing import Packing, read, parse          # noqa: E402
from .store import Store                           # noqa: E402

_ST = None


def record(parent, move, params, child, cpu_s, run=None, **kw):
    """one move-log row (+ the child, stored with the parent link) from a script that is not a stage: Store.record on a
    store opened once per process (11-02)"""
    global _ST
    _ST = _ST or Store()
    return _ST.record(parent, move, params, child, cpu_s, run=run, **kw)
