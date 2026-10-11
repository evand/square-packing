import Sqpack.Exact.Packs.N18.Data

namespace UnitSquarePacking.EC.N18

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 18, 0 ≤ i.val → i.val < 18 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N18
