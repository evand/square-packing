import Sqpack.Exact.Packs.N315.Data

namespace UnitSquarePacking.EC.N315

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 315, 275 ≤ i.val → i.val < 300 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N315
