import Sqpack.Exact.Packs.N96.Data

namespace UnitSquarePacking.EC.N96

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 96, 75 ≤ i.val → i.val < 96 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N96
