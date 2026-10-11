import Sqpack.Exact.Packs.N218.Data

namespace UnitSquarePacking.EC.N218

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 218, 200 ≤ i.val → i.val < 218 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N218
