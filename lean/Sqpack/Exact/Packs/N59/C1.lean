import Sqpack.Exact.Packs.N59.Data

namespace UnitSquarePacking.EC.N59

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 59, 25 ≤ i.val → i.val < 50 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N59
