import Sqpack.Exact.Packs.N202.Data

namespace UnitSquarePacking.EC.N202

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 202, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N202
