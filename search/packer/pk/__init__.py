"""pk: packing search toolkit (store, engine, moves, replay bench).  See ../DESIGN.md."""
import os, sys
_PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PKG not in sys.path:                     # legacy modules (layout, mcmin, ...) live next to the package
    sys.path.insert(0, _PKG)
from .packing import Packing, read, parse          # noqa: E402
from .store import Store                           # noqa: E402
