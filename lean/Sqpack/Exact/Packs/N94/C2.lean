import Sqpack.Exact.Packs.N94.Data

namespace UnitSquarePacking.EC.N94

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 94, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N94
