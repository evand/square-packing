import Sqpack.Exact.Packs.N129.Data

namespace UnitSquarePacking.EC.N129

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 129, 0 ≤ i.val → i.val < 10 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N129
