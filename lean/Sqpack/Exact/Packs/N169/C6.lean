import Sqpack.Exact.Packs.N169.Data

namespace UnitSquarePacking.EC.N169

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 169, 150 ≤ i.val → i.val < 169 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N169
