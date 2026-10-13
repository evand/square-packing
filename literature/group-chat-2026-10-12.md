# Group-chat reading list (2026-10-12)

Transcribed from three screenshots of a Zotero library that someone posted in the group chat.  The screenshots show only
title and first author, so this list has no years, venues or DOIs.  **Nothing below has been read or verified.**  The
"relevance" notes are my guesses from the titles [I].  "Covered" means our existing notes already cite the work or its
method.

## Likely relevant to the packer (continuous rotations, non-overlap formulations, containers as variables)

| Title | First author | Relevance [I] / covered |
|---|---|---|
| Out-of-the-Box Global Optimization for Packing Problems: New Models and Improved Solutions | Berthold et al. | **Covered**: arXiv:2605.04850; the lifted separator in `fq` comes from it (`search/packer/lit-optimizers.md`, REVIEW-2026-10-11). |
| Mathematical model and method of searching for a local extremum for the non-convex oriented polygons allocation problem | Stoyan et al. | Phi-function local search.  The phi-function family is covered in `lit-optimizers.md` §3 via Chernov–Stoyan–Romanova 2010; this specific paper is not. |
| An Optimization of Phi-function for Convex Polygons | Narendran Rajagopalan et al. | Phi-functions for convex polygons, which includes rotated unit squares. |
| Packing stretched convex polygons in an optimized rectangle | Bennell et al. | Phi-functions with a variable container. |
| Packing Soft Convex Polygons in an Optimized Convex Container | Litvinchev et al. | Variable container and soft (deformable) items; possibly a soft-to-rigid continuation (cf. `lit-2026-10-12.md` §6).  Appears twice in the screenshots. |
| Packing convex polygons in minimum-perimeter convex hulls | Kallrath et al. | MINLP with the container as a variable. |
| Cutting circles and polygons from area-minimizing rectangles | Kallrath | MINLP with the container as a variable. |
| Mixed integer quadratically-constrained programming model to solve the irregular strip packing problem with continuous rotations | Cherri et al. | MIQCP with continuous rotation; compare Berthold. |
| Packing Circles and Irregular Polygons using Separation Lines | Peralta et al. | Explicit separating lines: the same idea as SAT / Berthold's separator. |
| Method of sentinels for packing items within arbitrary convex regions | Birgin et al. | Sentinel-point non-overlap with smooth NLP (Birgin's group also did the classic circle-packing NLP work). |

## Squares (other problems)

| Title | First author | Relevance [I] |
|---|---|---|
| Packing, Hitting, and Colouring Squares | Caoduro and Sebő | Probably packing/hitting numbers of families of axis-parallel squares, not s(n).  Check before citing. |

## Surveys

| Title | First author | Relevance [I] |
|---|---|---|
| Two-dimensional irregular packing problems: A review | Guo et al. | Survey; a likely way into the nesting literature. |

## Heuristics, approximation, bin packing (2D)

| Title | First author | Relevance [I] |
|---|---|---|
| The optimal packing of arbitrarily-shaped polygons using simulated annealing and polynomial-time cooling schedules | Theodoracatos and Grimsley | Early SA nesting. |
| Optimizing Two-Dimensional Irregular Packing: A Hybrid Approach of Genetic Algorithm and Linear Programming | Liu et al. | GA + LP compaction. |
| Improved Approximations for Translational Packing of Convex Polygons | Kurpisz and Suter | Approximation algorithms, translations only. |
| A new approach for bin packing problem using knowledge reuse and improved heuristic | Fang et al. | Bin packing. |

## 3D (probably not relevant)

| Title | First author | Relevance [I] / covered |
|---|---|---|
| Dense Packings of Polyhedra: Platonic and Archimedean Solids | Torquato and Jiao | **Covered**: the ASC method (`lit-optimizers.md` §3, `lit-statmech.md`). |
| Packing of concave polyhedra with continuous rotations using nonlinear optimisation | Romanova et al. | 3D phi-functions. |
| Packing Oblique 3D Objects | Pankratov et al. | 3D phi-functions. |
| 3D Irregular Packing in an Optimized Cuboid Container | Litvinchev et al. | 3D, variable container. |
| Packing Polyhedra within Convex Sets | F. Callisaya and R. V. Kari | 3D. |
| HAPE3D—a new constructive algorithm for the 3D irregular packing problem | Liu et al. | 3D constructive. |
| Greedy Heuristic for Optimized Placement of Orthogonal Polyhedra in Irregular Cutting and Packing Problems | Chekanin and Chekanin | 3D, orthogonal. |
| Minimum Convex Container of Two Convex Polytopes under Translations | Ahn et al. | Computational geometry, two objects. |

## Follow-up

Read and compared with the packer in [`irregular-packing-2026-10-12.md`](irregular-packing-2026-10-12.md).
