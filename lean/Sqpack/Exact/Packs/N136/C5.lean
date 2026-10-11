import Sqpack.Exact.Packs.N136.Data

namespace UnitSquarePacking.EC.N136

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 136, 125 ≤ i.val → i.val < 136 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N136
