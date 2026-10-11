import Sqpack.Exact.Packs.N109.Data

namespace UnitSquarePacking.EC.N109

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 109, 100 ≤ i.val → i.val < 109 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N109
