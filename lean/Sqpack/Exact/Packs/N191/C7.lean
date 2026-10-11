import Sqpack.Exact.Packs.N191.Data

namespace UnitSquarePacking.EC.N191

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 191, 175 ≤ i.val → i.val < 191 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N191
