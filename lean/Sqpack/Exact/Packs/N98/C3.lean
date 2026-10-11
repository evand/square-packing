import Sqpack.Exact.Packs.N98.Data

namespace UnitSquarePacking.EC.N98

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 98, 75 ≤ i.val → i.val < 98 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N98
