import Sqpack.Exact.Packs.N174.Data

namespace UnitSquarePacking.EC.N174

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 174, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N174
