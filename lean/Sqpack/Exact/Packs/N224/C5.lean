import Sqpack.Exact.Packs.N224.Data

namespace UnitSquarePacking.EC.N224

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 224, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N224
