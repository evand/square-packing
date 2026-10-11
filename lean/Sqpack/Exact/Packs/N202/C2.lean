import Sqpack.Exact.Packs.N202.Data

namespace UnitSquarePacking.EC.N202

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 202, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N202
