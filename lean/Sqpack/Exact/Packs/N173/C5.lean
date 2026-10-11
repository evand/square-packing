import Sqpack.Exact.Packs.N173.Data

namespace UnitSquarePacking.EC.N173

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 173, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N173
