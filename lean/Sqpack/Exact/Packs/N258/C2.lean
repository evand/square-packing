import Sqpack.Exact.Packs.N258.Data

namespace UnitSquarePacking.EC.N258

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 258, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N258
