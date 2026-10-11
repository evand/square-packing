import Sqpack.Exact.Packs.N320.Data

namespace UnitSquarePacking.EC.N320

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 320, 275 ≤ i.val → i.val < 300 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N320
