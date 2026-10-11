import Sqpack.Exact.Packs.N52.Data

namespace UnitSquarePacking.EC.N52

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 52, 50 ≤ i.val → i.val < 52 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N52
