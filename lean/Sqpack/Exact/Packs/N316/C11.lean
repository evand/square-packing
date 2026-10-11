import Sqpack.Exact.Packs.N316.Data

namespace UnitSquarePacking.EC.N316

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 316, 275 ≤ i.val → i.val < 300 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N316
