import Sqpack.Exact.Packs.N171.Data

namespace UnitSquarePacking.EC.N171

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 171, 150 ≤ i.val → i.val < 171 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N171
