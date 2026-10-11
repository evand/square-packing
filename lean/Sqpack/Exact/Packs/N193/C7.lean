import Sqpack.Exact.Packs.N193.Data

namespace UnitSquarePacking.EC.N193

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 193, 175 ≤ i.val → i.val < 193 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N193
