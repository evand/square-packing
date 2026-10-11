import Sqpack.Exact.Packs.N45.Data

namespace UnitSquarePacking.EC.N45

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 45, 25 ≤ i.val → i.val < 45 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N45
