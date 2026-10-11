import Sqpack.Exact.Packs.N129.Data

namespace UnitSquarePacking.EC.N129

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 129, 80 ≤ i.val → i.val < 90 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N129
