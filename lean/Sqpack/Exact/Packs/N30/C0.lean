import Sqpack.Exact.Packs.N30.Data

namespace UnitSquarePacking.EC.N30

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 30, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N30
