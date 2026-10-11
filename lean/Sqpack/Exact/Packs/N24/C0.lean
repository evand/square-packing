import Sqpack.Exact.Packs.N24.Data

namespace UnitSquarePacking.EC.N24

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 24, 0 ≤ i.val → i.val < 24 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N24
