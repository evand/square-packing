import Sqpack.Exact.Packs.N222.Data

namespace UnitSquarePacking.EC.N222

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 222, 200 ≤ i.val → i.val < 222 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N222
