import Sqpack.Exact.Packs.N7.Data

namespace UnitSquarePacking.EC.N7

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 7, 0 ≤ i.val → i.val < 7 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N7
