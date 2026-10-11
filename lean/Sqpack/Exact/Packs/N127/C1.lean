import Sqpack.Exact.Packs.N127.Data

namespace UnitSquarePacking.EC.N127

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 127, 25 ≤ i.val → i.val < 50 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N127
