import Sqpack.Exact.Packs.N60.Data

namespace UnitSquarePacking.EC.N60

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 60, 50 ≤ i.val → i.val < 60 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N60
