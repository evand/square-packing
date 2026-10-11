import Sqpack.Exact.Packs.N265.Data

namespace UnitSquarePacking.EC.N265

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 265, 250 ≤ i.val → i.val < 265 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N265
