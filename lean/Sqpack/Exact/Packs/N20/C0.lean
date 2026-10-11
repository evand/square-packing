import Sqpack.Exact.Packs.N20.Data

namespace UnitSquarePacking.EC.N20

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 20, 0 ≤ i.val → i.val < 20 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N20
