import Sqpack.Exact.Packs.N12.Data

namespace UnitSquarePacking.EC.N12

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 12, 0 ≤ i.val → i.val < 12 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N12
