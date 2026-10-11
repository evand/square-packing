import Sqpack.Exact.Packs.N40.Data

namespace UnitSquarePacking.EC.N40

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 40, 25 ≤ i.val → i.val < 40 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N40
