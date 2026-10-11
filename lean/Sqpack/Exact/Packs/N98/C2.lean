import Sqpack.Exact.Packs.N98.Data

namespace UnitSquarePacking.EC.N98

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 98, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N98
