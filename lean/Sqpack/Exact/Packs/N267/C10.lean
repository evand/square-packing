import Sqpack.Exact.Packs.N267.Data

namespace UnitSquarePacking.EC.N267

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 267, 250 ≤ i.val → i.val < 267 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N267
