import Sqpack.Exact.Packs.N113.Data

namespace UnitSquarePacking.EC.N113

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 113, 100 ≤ i.val → i.val < 113 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N113
