import Sqpack.Exact.Packs.N257.Data

namespace UnitSquarePacking.EC.N257

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 257, 250 ≤ i.val → i.val < 257 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N257
