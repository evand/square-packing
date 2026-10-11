import Sqpack.Exact.Packs.N158.Data

namespace UnitSquarePacking.EC.N158

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 158, 150 ≤ i.val → i.val < 158 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N158
