# bandcut-cost: the accounting of a cut — why a tilted band pays at `T >= 12` and what it would have to pay at `T = 4`  (2026-09-20)

**Why.**  `search/ARCH_TLEDGER.md` §2: for a wall-to-wall chain of `T` squares at a common tilt `t` with rise `R`,
`delta <= [(T-u) cos t + R sin t]/(T-1) - 1`, `u = cos t + sin t` — the `T`-dependence is in the rise.  Known counterexamples
(`T >= 12`, Arslanov–Mustafin–Shangitbayev EJC 28(4) 2021 P4.22; Cleemann `T = 17`) cut every row and column with a band of
four-chains at `2 arctan(1/4)`, where `4 cos t + sin t = 4`.  A proof of `s(12) = 4` must contain a statement that is false
at `T = 12`; the candidate is "a cut costs more than the container can afford".  Find that statement.

**Do.**  Pen-and-paper first, small numeric checks only.  (1) Generalise the chain inequality to **mixed tilts**: a
wall-to-wall path in the separation DAG through squares of arbitrary angles `theta_i`, links with widths
`1/2 + W(theta_j - theta_i)/2` along an edge normal of either square (`search/S6_SKELETON.md` §3.1); what is the exact
bound, what plays the role of the rise, and when can a path through a band be *shorter* than `T`?  (2) Read the Arslanov
et al. construction and do its accounting: what is gained per band chain, what is lost at the band's two interfaces, at the
elbow and at the walls; express the net as a function of `T` and find where it changes sign.  Does the accounting explain
`T = 12` (and `11`)?  (3) Turn it round: state the strongest lemma of the form "any configuration with no wall-to-wall
chain of `T` near-axis squares has `delta* <= -f(T) < 0` for `T <= T_0`" that the accounting *suggests*, with the constants,
and list exactly what would have to be proved.  Run the filter: it must be false at `T = 12`, and must not prove
`s(5) = 3` (`n = T^2 - 4`, `T = 3`).  Be explicit about what is proved, what is heuristic, and what is a guess.
**Deliverable.**  `notes/bandcut-cost.md` (+ `search/bandcut_cost*.py` for checks).

**Addendum (same night, after `search/BANDCUT_SCAN.md`).**  Read its §1, §4, §5 first: the published scheme is two integer
blocks + two *squeezable* rectangles with `waste(A) + waste(B) = T` exactly; smallest known squeezable rectangle is `(4,8)`
with 26 squares, waste `6`, hence `T >= 12`; the paper's mechanism and tables are already transcribed there, so the PDF need
not be re-read in full.  `T = 12` and `T = 11` are verified positive in exact arithmetic.  So a second candidate form of the
`T`-dependent statement is a **waste lower bound for squeezable regions** — assess it next to the chain/rise accounting: same
statement or not, provable in any regime (cf. Roth–Vaughan), filter.  At `T = 4` a wall-to-wall `(4,1)` stack at `28.59°`
beside eight axis-parallel squares reaches margin exactly `0` but still contains an axis-parallel chain of four.
**Working style:** the first attempt at this task died by exceeding the output limit in a single response.  Work in small
steps; write `notes/bandcut-cost.md` section by section (append), never in one giant message; keep derivations on disk.
