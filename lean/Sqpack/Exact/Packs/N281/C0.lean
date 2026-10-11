import Sqpack.Exact.Packs.N281.Data

namespace UnitSquarePacking.EC.N281

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 281, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N281
