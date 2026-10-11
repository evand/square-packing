import Sqpack.Exact.Packs.N23.Data

namespace UnitSquarePacking.EC.N23

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 23, 0 ≤ i.val → i.val < 23 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N23
