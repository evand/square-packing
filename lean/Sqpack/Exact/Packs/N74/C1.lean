import Sqpack.Exact.Packs.N74.Data

namespace UnitSquarePacking.EC.N74

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 74, 25 ≤ i.val → i.val < 50 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N74
