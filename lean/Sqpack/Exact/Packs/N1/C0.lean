import Sqpack.Exact.Packs.N1.Data

namespace UnitSquarePacking.EC.N1

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 1, 0 ≤ i.val → i.val < 1 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N1
