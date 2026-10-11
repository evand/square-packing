import Sqpack.Exact.Packs.N186.Data

namespace UnitSquarePacking.EC.N186

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 186, 175 ≤ i.val → i.val < 186 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N186
