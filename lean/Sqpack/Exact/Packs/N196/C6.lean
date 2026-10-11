import Sqpack.Exact.Packs.N196.Data

namespace UnitSquarePacking.EC.N196

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 196, 150 ≤ i.val → i.val < 175 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N196
