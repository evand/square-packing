import Sqpack.Exact.Packs.N38.Data

namespace UnitSquarePacking.EC.N38

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 38, 25 ≤ i.val → i.val < 38 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N38
