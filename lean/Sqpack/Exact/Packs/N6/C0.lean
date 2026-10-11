import Sqpack.Exact.Packs.N6.Data

namespace UnitSquarePacking.EC.N6

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 6, 0 ≤ i.val → i.val < 6 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N6
