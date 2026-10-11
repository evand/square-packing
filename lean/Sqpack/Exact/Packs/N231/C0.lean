import Sqpack.Exact.Packs.N231.Data

namespace UnitSquarePacking.EC.N231

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 231, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N231
