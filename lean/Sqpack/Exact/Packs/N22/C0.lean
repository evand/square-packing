import Sqpack.Exact.Packs.N22.Data

namespace UnitSquarePacking.EC.N22

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 22, 0 ≤ i.val → i.val < 22 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N22
