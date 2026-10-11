import Sqpack.Exact.Packs.N13.Data

namespace UnitSquarePacking.EC.N13

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 13, 0 ≤ i.val → i.val < 13 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N13
