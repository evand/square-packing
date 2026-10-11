import Sqpack.Exact.Packs.N89.Data

namespace UnitSquarePacking.EC.N89

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 89, 75 ≤ i.val → i.val < 89 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N89
