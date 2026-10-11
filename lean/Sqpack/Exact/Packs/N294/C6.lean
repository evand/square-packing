import Sqpack.Exact.Packs.N294.Data

namespace UnitSquarePacking.EC.N294

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 294, 150 ≤ i.val → i.val < 175 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N294
