import Sqpack.Exact.Packs.N283.Data

namespace UnitSquarePacking.EC.N283

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 283, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N283
