# packer: searching for square packings

Started 2026-10-04.  Goal: better packings of n unit squares in a square, first target s(90) < 10 (WISHLIST §P1); method goal:
a general solver (any n), since all known approaches fail in ways we now partly understand.  `PACKER.md` is the chronological lab
log (what was tried, in order, with numbers); this file is the map.  Everything here is **screening** (f64) unless stated.

## Status (10-06; 10-04 text below kept for the numbers)
* 10-06: s(110) calibration answered "where is the needle" (entering the funnel); s(90) first campaign 10.0095668 (no sub-10); now samplers that accept uphill steps (`mcmin.py`) and a neighbour record sweep (`sweep.py`).  Plan below.
* **Pass test not passed.**  Gate before spending compute on s(90): reproduce an s(110) < 11 packing without s(110)-specific
  hints.  Best from family structure alone: 11.00765 (jammed, no full lines).  The record's funnel was never entered from outside.
* **s(110) record = Couzo's 2026-09-27 packing (`couzo110`; older notes and run labels call it "cand" / "candidate"), which we rediscovered on 10-04** (jlevy register: 10.996783396634359).  Exact optimum of the
  basin certified at 10.996783396631591639612746267957 (`../exact/results/couzo110.cert`; `../exact/verify_cert.py`, `verify_cert2.py`),
  2.8e-12 below Couzo's printed bound.  Not a new record.
* **Clean census (10-05, `runs/cen7`, 1000 trials, exact dedupe):** 19 distinct sub-11 minima, top four hold 91 %, best =
  Couzo's (2nd most common); unseen mass ≈ 0.005 pooled, ≈ 0.09 at σ = 0.1; σ50 ≈ 0.11.  Entering the funnel is the hard part.
* **Landscape near s(110) (10-04, partly superseded by cen7):** many near-degenerate jammed minima (≥ 10 junction classes, ~2e-3 spread), separated by small ridges
  (Ellsworth's → Couzo's ≤ 6.7e-4 on one path); escape σ50 ≈ 0.1–0.15 under soft-squeeze + LP finish.

## Tools
| file | what |
|---|---|
| `src/main.rs` (`packer`) | Rust engine.  `relax [--squeeze \| --squeeze-pen]`: L-BFGS on SAT overlap energy at fixed side; bisection squeeze; penalty-continuation squeeze over (config, s) (μ0 = stiffness: low = soft/exploring, high = precise).  `hop`: basin hopping, descending target.  `hhop`: hard-square basin hopping with event chains (slide until contact, pass remainder on), `--forbid-lines k`, `--soft p` (soften/re-harden). |
| `slp2.py` | **The local minimiser.**  Trust-region SLP on min s s.t. linearised separation, incidence contact model (`inc.py`), second-order correction, repair by uniform scaling; final jam test exhaustive over corner–corner branches (MILP).  20–50 iterations, ~3–5 s per descent at n = 110; outputs go straight into `../exact/exactsolve.py`.  `budget=` wall cap. |
| `inc.py` | Incidence contact model (f64 twin of `../exact/find_contacts`): corner-on-side rows; corner–corner pairs as disjunctions (relax-and-round LP; `shrink_milp` exhaustive). |
| `jobpool.py` | Process pool with a hard per-job kill (a HiGHS call can hang in C; `multiprocessing.Pool` cannot recover).  All census drivers use it. |
| `rigid.py` | Old one-axis-per-pair jam LP (false jams: kept for `slp2 --axis` comparison) and per-coordinate rigidity.  Pins BLAS threads to 1. |
| `gen.py` | Geometry helpers; `mis_fill2`: exact max axis fill (corner lattices + staircase candidates; HiGHS MILP).  Recovers record counts of 54, 71, 87, 88, 110, 132 from their tilted squares. |
| `layout.py` | Channel-row constructions (w, θ, set-back) + exact fill; `hyb` mode: genome MC scored by squeezed side, rejecting full lines. |
| `ftmc.py` | Free-tilted MC: perturb tilted squares, re-solve axis fill exactly, squeeze, Metropolis on side. |
| `census_exact.py` | Exact dedupe of census3 runs (every sub-k output through `../exact/exactsolve.py`, classes by S to 1e-20) + coupon statistics (f1/N, Chao1) per run, σ, pooled. |
| `census_roles.py` | Square roles per census trial (bands / axis / other) vs the start, by index. |
| `movegen.py` | Structural band moves (add / del / transfer / row / shift / rotate / drop_J; `--focus` adds uniform / rebuild) with keep-axis evaluation → slp2; several starts per run. |
| `beam.py` | Count-pair beam: best packing per role-count key from movegen/crop runs → movegen from them (`--n --k --procs --focus`). |
| `crop.py` | Side k+1 → k seeds: remove a row and a column strip, shift, fix the count, squeeze, slp2 (132 → 110 pass test, 110 → 90). |
| `mcmin.py` | Minima PT: Metropolis over slp2-quenched minima with replica exchange; moves kick / lkick / reinsert (hole-weighted) / crot / aswap / band, 1 + Poisson(λ) per proposal; time-based epochs, swap sweeps; `log.jsonl` per proposal. |
| `sweep.py` | Neighbour sweep over register records (`../exact/batch/inputs`): n+1 minus a square, kicks of n, quench, compare to register. |
| `census3.py` | Perturb (σ) → soft squeeze → slp2 → side, jam, full lines, fingerprint.  `basin.py`, `census2.py`: earlier variants. |
| `contacts.py`, `stress.py` | Feature-level contact graphs (corner-on-side, side-side, corner-on-wall); force networks (jam-LP duals). |
| `homotopy.py` | Rotate chosen squares between two configurations with everything else re-minimised; prints s(t) (barrier estimates). |
| `channel.py`, `surgery.py`, `gadget_mc.py`, `mkseed.py`, `slp.py` | Earlier experiments (superseded; kept for the record). |
| `lit-*.md` | Literature: optimizer techniques, stat-mech / basin counting, provenance of records n = 50–135. |

Formats: text `n s` + `x y deg` per line (centres).  `site2txt.py` / `txt2site.py` convert the site JSON.  Build: `cargo build --release`.

## Lessons (the short version; evidence in PACKER.md)
0. **Compare against jlevy's register, not our site mirror**: the mirror (Ellsworth's page, 09-30) is stale at the frontier
   (Couzo 09-27: s(110), s(130), s(132)).
1. **Grid attractor.**  Any state with a full line of k axis squares has side ≥ k; these states are a huge attractor at every
   level (random search, hopping, event chains, layout genomes all drift to them).  Reject full lines explicitly.
2. **s < k needs a transversal**: the tilted part must cross every row and column.  Local moves don't create one; structure does.
3. **Descent map matters more than it looks.**  L-BFGS on the SAT penalty stalls at kinks (face-to-face contacts): its "minima"
   are not first-order jammed (LP finds shrink directions), which made early censuses and search comparisons noisy at ~1e-3.
   slp2 fixes this.  Soft penalty squeeze (low μ0) is a useful *coarse* descent (funnel-level compaction); slp2 is the finisher.
4. **Integer count is the wrong objective** (needle landscape); squeezed side with full-line rejection works.
5. **Records' channels are not rigid lattice rows** (staggered, mixed angles, resting on axis squares): row genomes can't express
   them; they find relatives in other funnels (our 11.0076).
6. Records near n = 50–135 mostly come from constructions + seeded annealing; Ellsworth: plain SA "almost always gets stuck just
   above the trivial size due to stacked rows and/or columns".

7. **(10-06)** Grid obstructions include *staggered* axis chains (k axis squares, consecutive ones overlapping in the other
   projection), not just straight lines (`layout.axis_chain`).  Under our moves the grid is absorbing: hot replicas sink into it.
8. **(10-06)** Temperature alone doesn't connect funnels (minima PT runs A, B at 110 missed); reinsert-into-hole moves are
   destructive in every variant tried (a jammed packing has no hole: best insert costs ~0.4 penetration).
9. **(10-06)** Kick + good descent is short-range: reliable to matching distance ~0.02, rare beyond 0.05 (rediscovery, 51 record
   pairs).  It still finds records: s(266), s(270), s(272) (jlevy#399).  Kicking *older* packings finds basins the current
   record's kicks don't (272).  At n ≈ 270 give slp2 ≥ 180 s (60 s left 13 % unconverged, and that changed detection).

## Plan (10-06, with Evan)
State: the s(110) calibration answered its main question (PACKER.md cen7, rev1, jx1, crop pass test): below 11 a smooth funnel
(≥ ~47 exact minima, best basins common, descent to the bottom easy once inside); above 11 a rough plateau of small, densely
connected basins with a few gateways.  Per-coordinate kicks can't move squares between bands; the greedy beam can't go uphill (s(90)
stalled at 10.0096, the 132 → 110 crop reached the rim only).  The 10-05 steps 1–2 are done; 3 is replaced by the items below; 4
(needle estimate for s(90)) stays as the success criterion for the samplers; 5 (case study) waits for the sampler results.

1. **Minima PT** (`mcmin.py`): Metropolis over quenched minima + replica exchange; content-agnostic moves (kicks, hole-weighted
   reinsert, cluster rotation, angle swap) beside band moves; 1 + Poisson(λ) moves per proposal.  Pre-registered calibration at 110
   (PACKER "Minima PT"): A rim → bottom, B 11.0076 funnel → sub-11.  Then: score move kinds by where they lead (not acceptance),
   re-space the ladder offline from per-pair acceptance, excursion N_eff = N / (1 + CV²).
2. **Barrier type.**  If the plateau → funnel transition is entropic (first-order), no temperature ladder opens it (a lesson from an earlier
   MH/PT project): bias in side (multicanonical / Wang–Landau), or gap-spanning proposals (constructors as independence proposals).
3. **Record attempt** (`sweep.py`): every open non-integer n, record(n+1) minus a square and kicks of record(n), quenched; then
   minima PT on whatever n looks soft.  Register check before any claim.
4. **Landscape measure** (content-agnostic): nested sampling / multilevel splitting on side from random starts, log X(s) at 110;
   comparable across n ("how hard is this n"), unlike Good–Turing, which is protocol-relative and can't see narrow deep basins.
5. **s(90)** after 1–2 give a detection rate.  Constructors for k = 10 other than cropping.

**Side bet: rectangle containers** (unchanged): slp2 with an (a−δ) × b container; Arslanov's 26 in (4−δ) × 8 control, then
jlevy H-049 (20 in (4−δ) × 6).  **Hygiene:** slp2 speed (contact-model pair loop ~1/3 of a descent); register before claims.

## Open threads (see also `../../TODO.md`)
* Exact solver (`../exact/`, done): KKT Newton + rational certificates; known gap: corner–corner contacts first-order only.
* Census: counting contact graphs honestly needs the force-bearing (max-support dual) network (parked).
* Funnel-level moves (junction flips; structure classes).
* Other targets for the same instrument: s(147) and the rest of WISHLIST §P; other "not yet analytically optimized" records.
