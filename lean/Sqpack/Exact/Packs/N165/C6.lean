import Sqpack.Exact.Packs.N165.Data

namespace UnitSquarePacking.EC.N165

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 165, 150 ≤ i.val → i.val < 165 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N165
