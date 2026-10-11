import Sqpack.Exact.Packs.N115.Data

namespace UnitSquarePacking.EC.N115

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 115, 100 ≤ i.val → i.val < 115 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N115
