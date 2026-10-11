import Sqpack.Exact.Packs.N168.Data

namespace UnitSquarePacking.EC.N168

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 168, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N168
