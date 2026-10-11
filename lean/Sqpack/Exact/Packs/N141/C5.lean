import Sqpack.Exact.Packs.N141.Data

namespace UnitSquarePacking.EC.N141

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 141, 125 ≤ i.val → i.val < 141 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N141
