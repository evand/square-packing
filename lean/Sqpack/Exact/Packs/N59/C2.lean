import Sqpack.Exact.Packs.N59.Data

namespace UnitSquarePacking.EC.N59

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 59, 50 ≤ i.val → i.val < 59 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N59
