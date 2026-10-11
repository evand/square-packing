import Sqpack.Exact.Packs.N27.Data

namespace UnitSquarePacking.EC.N27

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 27, 25 ≤ i.val → i.val < 27 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N27
