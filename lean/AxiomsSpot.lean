-- CI spot check (the full list is Axioms.lean; build everything locally with `lake build`).
-- These modules stay under ~7 GB per Lean process; Sqpack.Bentz (23.5 GB), S32 (11.8 GB) and their
-- dependents do not fit a 16 GB hosted runner, so CI does not build them.
import Sqpack.Spec
import Sqpack.SpecFC
import Sqpack.ZeroMargin
import Sqpack.ChordLine
import Sqpack.Cover
import Sqpack.Exact.N5c
import Sqpack.Exact.N11c
import Sqpack.LocalMinCheck
import Sqpack.ChainLocalMin
import Sqpack.Exact.N8Chain
import Sqpack.ConjecturesProofs
#print axioms SquarePacking.lemmaA
#print axioms UnitSquarePacking.setOf_packs_eq
#print axioms UnitSquarePacking.EC.N5c.packs
#print axioms UnitSquarePacking.EC.N5c.packs_closed
#print axioms UnitSquarePacking.EC.N5c.minSide_le_closed
#print axioms UnitSquarePacking.EC.N11c.packs
#print axioms UnitSquarePacking.LMC.isLocalMin_of_lcert
#print axioms UnitSquarePacking.isLocalMin_of_band
#print axioms UnitSquarePacking.Chain.N8.localMin
#print axioms UnitSquarePacking.Conjectures.tilingBound
#print axioms UnitSquarePacking.Conjectures.cStarStepsOfTwo
#print axioms UnitSquarePacking.Conjectures.tilingsEventuallyNotOptimal_of
