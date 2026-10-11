import Sqpack.Exact.Packs.N120.Data

namespace UnitSquarePacking.EC.N120

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 120, 100 ≤ i.val → i.val < 120 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N120
