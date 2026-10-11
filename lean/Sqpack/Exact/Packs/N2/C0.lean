import Sqpack.Exact.Packs.N2.Data

namespace UnitSquarePacking.EC.N2

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 2, 0 ≤ i.val → i.val < 2 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N2
