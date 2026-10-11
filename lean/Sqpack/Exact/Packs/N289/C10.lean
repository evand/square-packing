import Sqpack.Exact.Packs.N289.Data

namespace UnitSquarePacking.EC.N289

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 289, 250 ≤ i.val → i.val < 275 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N289
