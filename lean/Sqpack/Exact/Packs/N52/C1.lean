import Sqpack.Exact.Packs.N52.Data

namespace UnitSquarePacking.EC.N52

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 52, 25 ≤ i.val → i.val < 50 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N52
