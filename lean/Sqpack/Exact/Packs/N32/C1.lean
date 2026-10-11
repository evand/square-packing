import Sqpack.Exact.Packs.N32.Data

namespace UnitSquarePacking.EC.N32

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 32, 25 ≤ i.val → i.val < 32 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N32
