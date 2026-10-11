import Sqpack.Exact.Packs.N293.Data

namespace UnitSquarePacking.EC.N293

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 293, 200 ≤ i.val → i.val < 225 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N293
