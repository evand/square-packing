import Sqpack.Exact.Packs.N140.Data

namespace UnitSquarePacking.EC.N140

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 140, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N140
