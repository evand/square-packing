import Sqpack.Exact.Packs.N188.Data

namespace UnitSquarePacking.EC.N188

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 188, 150 ≤ i.val → i.val < 175 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N188
