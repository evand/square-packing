import Sqpack.Exact.Packs.N4.Data

namespace UnitSquarePacking.EC.N4

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 4, 0 ≤ i.val → i.val < 4 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N4
