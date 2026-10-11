import Sqpack.Exact.Packs.N283.Data

namespace UnitSquarePacking.EC.N283

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 283, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N283
