import Sqpack.Exact.Packs.N279.Data

namespace UnitSquarePacking.EC.N279

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 279, 200 ≤ i.val → i.val < 225 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N279
