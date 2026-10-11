import Sqpack.Exact.Packs.N316.Data

namespace UnitSquarePacking.EC.N316

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 316, 25 ≤ i.val → i.val < 50 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N316
