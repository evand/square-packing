import Sqpack.Exact.Packs.N279.Data

namespace UnitSquarePacking.EC.N279

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 279, 150 ≤ i.val → i.val < 175 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N279
