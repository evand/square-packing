import Sqpack.Exact.Packs.N94.Data

namespace UnitSquarePacking.EC.N94

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 94, 25 ≤ i.val → i.val < 50 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N94
