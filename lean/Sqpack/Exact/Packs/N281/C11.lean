import Sqpack.Exact.Packs.N281.Data

namespace UnitSquarePacking.EC.N281

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 281, 275 ≤ i.val → i.val < 281 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N281
