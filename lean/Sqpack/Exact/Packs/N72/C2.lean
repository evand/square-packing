import Sqpack.Exact.Packs.N72.Data

namespace UnitSquarePacking.EC.N72

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 72, 50 ≤ i.val → i.val < 72 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N72
