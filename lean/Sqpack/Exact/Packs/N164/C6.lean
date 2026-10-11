import Sqpack.Exact.Packs.N164.Data

namespace UnitSquarePacking.EC.N164

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 164, 150 ≤ i.val → i.val < 164 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N164
