import Sqpack.Exact.Packs.N232R3.Data

namespace UnitSquarePacking.EC.N232R3

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 232, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N232R3
