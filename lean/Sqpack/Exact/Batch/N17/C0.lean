import Sqpack.Exact.Batch.N17.Data

namespace UnitSquarePacking.EC.N17

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 17, 0 ≤ i.val → i.val < 17 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N17
