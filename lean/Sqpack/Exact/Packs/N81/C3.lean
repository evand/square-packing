import Sqpack.Exact.Packs.N81.Data

namespace UnitSquarePacking.EC.N81

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 81, 75 ≤ i.val → i.val < 81 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N81
