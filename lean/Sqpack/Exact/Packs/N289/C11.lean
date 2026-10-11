import Sqpack.Exact.Packs.N289.Data

namespace UnitSquarePacking.EC.N289

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 289, 275 ≤ i.val → i.val < 289 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N289
