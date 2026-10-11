import Sqpack.Exact.Packs.N185.Data

namespace UnitSquarePacking.EC.N185

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 185, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N185
