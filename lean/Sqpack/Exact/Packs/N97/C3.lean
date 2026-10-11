import Sqpack.Exact.Packs.N97.Data

namespace UnitSquarePacking.EC.N97

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 97, 75 ≤ i.val → i.val < 97 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N97
