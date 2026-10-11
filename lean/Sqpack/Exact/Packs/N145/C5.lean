import Sqpack.Exact.Packs.N145.Data

namespace UnitSquarePacking.EC.N145

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 145, 125 ≤ i.val → i.val < 145 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N145
