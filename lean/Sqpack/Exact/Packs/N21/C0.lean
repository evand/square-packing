import Sqpack.Exact.Packs.N21.Data

namespace UnitSquarePacking.EC.N21

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 21, 0 ≤ i.val → i.val < 21 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N21
