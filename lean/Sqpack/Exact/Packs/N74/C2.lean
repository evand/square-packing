import Sqpack.Exact.Packs.N74.Data

namespace UnitSquarePacking.EC.N74

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 74, 50 ≤ i.val → i.val < 74 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N74
