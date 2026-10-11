import Sqpack.Exact.Packs.N203.Data

namespace UnitSquarePacking.EC.N203

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 203, 175 ≤ i.val → i.val < 200 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N203
