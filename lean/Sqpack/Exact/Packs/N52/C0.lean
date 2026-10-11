import Sqpack.Exact.Packs.N52.Data

namespace UnitSquarePacking.EC.N52

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 52, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N52
