import Sqpack.Exact.Packs.N278.Data

namespace UnitSquarePacking.EC.N278

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 278, 175 ≤ i.val → i.val < 200 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N278
