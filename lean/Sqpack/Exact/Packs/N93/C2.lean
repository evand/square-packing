import Sqpack.Exact.Packs.N93.Data

namespace UnitSquarePacking.EC.N93

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 93, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N93
