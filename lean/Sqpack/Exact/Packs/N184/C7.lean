import Sqpack.Exact.Packs.N184.Data

namespace UnitSquarePacking.EC.N184

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 184, 175 ≤ i.val → i.val < 184 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N184
