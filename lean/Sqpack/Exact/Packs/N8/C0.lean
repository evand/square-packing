import Sqpack.Exact.Packs.N8.Data

namespace UnitSquarePacking.EC.N8

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 8, 0 ≤ i.val → i.val < 8 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N8
