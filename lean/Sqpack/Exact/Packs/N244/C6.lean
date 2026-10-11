import Sqpack.Exact.Packs.N244.Data

namespace UnitSquarePacking.EC.N244

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 244, 150 ≤ i.val → i.val < 175 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N244
