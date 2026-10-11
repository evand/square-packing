import Sqpack.Exact.Packs.N3.Data

namespace UnitSquarePacking.EC.N3

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 3, 0 ≤ i.val → i.val < 3 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N3
