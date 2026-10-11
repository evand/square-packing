import Sqpack.Exact.Packs.N79.Data

namespace UnitSquarePacking.EC.N79

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 79, 75 ≤ i.val → i.val < 79 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N79
