import Sqpack.Exact.Packs.N319.Data

namespace UnitSquarePacking.EC.N319

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 319, 175 ≤ i.val → i.val < 200 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N319
