import Sqpack.Exact.Packs.N312.Data

namespace UnitSquarePacking.EC.N312

set_option maxHeartbeats 0

theorem chunk_12 : ∀ i : Fin 312, 300 ≤ i.val → i.val < 312 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N312
