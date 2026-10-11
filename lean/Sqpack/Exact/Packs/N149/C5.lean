import Sqpack.Exact.Packs.N149.Data

namespace UnitSquarePacking.EC.N149

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 149, 125 ≤ i.val → i.val < 149 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N149
