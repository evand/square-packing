import Sqpack.Exact.Packs.N86.Data

namespace UnitSquarePacking.EC.N86

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 86, 75 ≤ i.val → i.val < 86 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N86
