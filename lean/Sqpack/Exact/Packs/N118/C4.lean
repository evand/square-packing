import Sqpack.Exact.Packs.N118.Data

namespace UnitSquarePacking.EC.N118

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 118, 100 ≤ i.val → i.val < 118 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N118
