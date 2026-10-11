import Sqpack.Exact.Packs.N67.Data

namespace UnitSquarePacking.EC.N67

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 67, 50 ≤ i.val → i.val < 67 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N67
