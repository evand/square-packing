# TODO

One line per task; details in the linked notes.  Status log: `notes/review-2026-09-13.md` (latest), `notes/review-2026-09-12.md`, `notes/review-2026-09-11.md`, `notes/review-2026-09-08.md`, `notes/review-2026-09-07.md`,
`notes/review-2026-08-29.md`, `notes/TODO-archive-2026-08-28.md`.  Results: `search/*.md`.

Rules of thumb: **explore on the packing side, certify on the cover side**; nothing below `t = 4`
transfers to a proof of `s(12) = 4`, so measure every idea at `t = 4` (closed semantics) first;
launch anything longer than 10 minutes detached (`setsid nohup`) so it survives a session; never
`pkill`.  16 threads is the machine's useful width (8 → 16 gives 1.8×, 16 → 32 nothing).

## Where the proof stands (2026-09-13)

Certified: `s(12) >= 3.968616`.  **`s(13) = 4` case-free** (`certificates/rung2/`, `search/RUNG2.md`)
has two exhaustive checkers sharing nothing (`search/zeromargin.py`; `verify2/zmcheck`, full domain,
0 uncertified, 17 min), 23 rejection tests, Lean for every soundness lemma, CI green — written up
(`notes/s13-casefree.md`, README, `docs/`).

For `n = 12` at the container (`t = 4`, closed semantics), **every degree-1 certifiable family is
`>= 12`**, and the 2026-09-13 measurements (`search/BENTZ.md`, `notes/review-2026-09-13.md`) say what
they all relax:

| what | value | status |
|---|---|---|
| pure, points only | `>= 12.2688` exact | `COVER4.md` |
| pure, points + polygons | `12.173398` certified | `pgonly_pure_exact.txt` |
| corner `k = 3` (`1,1,1,0`), points + polygons + regions + chord | **`12.027273` certified** (`M < 1`, corner masses off by a `9e-4` boundary tie-break; LP `12.038` still rising) | `BENTZ.md` §7.2 |
| corner `k = 4`, same family | **`11.999999926` certified**, `<= 12` | `pgonly_corner_exact.txt`, `LEAF_CEILING.md` §5.2 |
| corner `k = 4`, **fully pinned Bentz pattern leaf** (four `{a_i,b_i}` corners, eight singleton `{c_j}`/`{d_j}`) | **`11.999999926 <= value <= 12` exactly** — the `PGCORN` measure lies inside the leaf | `BENTZ.md` §0 |
| corner `k = 4`, **axis-parallel poses only** | **exactly `9`** (integral witness; `{1,2,3}²` is a closed cover) | `BENTZ.md` §0, C1 |
| corner `k = 4`, `\|θ mod 90°\| <= 1°` | **`11.599390` certified** (`M = 1`) | `BENTZ.md` §7.2 |
| sound cliques of any positive-volume rule | honest cost `>= 18` | `ALLMEET.md` |
| level-2 slot leaves under `k = 4`, no cliques | `11.92`, rising, pose-set-starved, `M > 1` — the one count level not yet shown `>= 12` | `T4LEAF.md` |

Read together: the obstruction is **rotation, and infinitesimal rotation** — axis-parallel squares
are trivially `<= 9`, at least `2.6` of the missing `3` is bought by tilts under one degree (grid
squares nudged `1e-3` and tilted a fraction of a degree, dodging every grid-line point at `O(θ)`
bounding-box cost), and the LP averages over that continuum.  The integrality gap is **exactly one
square and is visible in integers**: the optimum's own 162-pose support admits only 11 pairwise-disjoint
squares (`α = 11`, no violated rank-family inequality); the four corner pairs, four `{c_j}` and four
`{d_j}` are each realisable alone, and **the eight non-corner singletons together admit only seven**.
The theorem a proof needs is that rank-8 statement — with four corner squares holding `{a_i, b_i}`,
at most seven further squares each contain one of `c_j, d_j` — true integrally (one wall square per
wall, then `s(4) = 2` on the interior) and inexpressible by any single-square inequality.  Count and
pattern branching cannot reach it: a count row cannot tell one square from a cloud of near-squares.

## Critical path

1. **Publish `s(13) = 4`** — the write-up is done; remaining are the doc follow-ups under Next and the
   CI decision (dispatch-only vs split) before the public repo.
2. **Write the `n = 12` negative result as a self-contained note** (`notes/n12-gap.md`): what every
   degree-1 family does at `t = 4` (the table above), the axis-parallel `9`, the sub-degree-tilt
   `>= 2.6`, the `α = 11` localisation and the rank-8 statement.  It is the precise statement of why
   the problem is hard and it is new; it is also the spec for anything that comes after.
3. **The rank-8 statement, as mathematics** — the only route left that needs no new relaxation:
   prove that with four corner squares holding `{a_i, b_i}`, at most seven squares can each contain
   one of the eight `c_j, d_j`.  Every square in it is near its tile and tilted by under a degree, so
   the tool is a **linearised rigidity argument around the tiling** (first-order disjointness along
   rows and columns; Nagamochi / Bentz Corollary-7 `"> 1"` lemmas for what closed semantics loses),
   with the machine checking leaves *outside* the neighbourhood where margins are positive.  Start by
   writing the first-order constraint system for the 12 perturbed tiles and testing its
   infeasibility numerically; then ask what neighbourhood it covers.
4. **Degree-2 (research option).**  Sherali–Adams level-2 on the twelve *labelled* squares of the
   fully pinned leaf — pair marginals supported on disjoint poses — measured on the 162-pose support
   first (finite, cheap: does it read 11?), then the continuum-crediting design question.  Not before
   item 3 has been tried on paper.
5. **`s(11)` as the hedge**: the same rotation-carried gap (`3.8143` certified vs Trump `3.877083`,
   LP crossing `~3.816`, `search/N11.md`), all of it `t < 4` work with the sweep verifier and the branch
   trailer.  First question: the anatomy of the mass-11 fractional packing at `t ≈ 3.82–3.85`
   (`dual_exact.py`-style measure extraction) — near-Trump smear, or the 3×3 grid in disguise?
6. V1 (disjunctive `t = 4` verifier) — **only after a certificate-shaped object below 12 exists.**
   None does; no count or pattern tree at `t = 4` will produce one.

## Running

- [x] `tasks/bentz-incidence` (Opus, 2026-09-13, merged at `493a530`; `search/bentz.py`, `search/BENTZ.md`,
      logs in `runs/bentz-2026-09-13/`): fully pinned Bentz leaf `[11.999999926, 12]` exactly; tree
      93 patterns / 86,403 leaves / 10,945 `D4` classes; axis-parallel exactly `9`; `|θ| <= 1°`
      `>= 11.599` certified; corner `k = 3` `>= 12.027` certified; `α = 11` on the 162-pose support.
      Nothing left running.  `cliquelever.py` gained an inert `--pose-filter` hook.
- [x] `PGPURE`, `PGCORN`: the no-clique measurement (`runs/launch_2026-09-12_pgonly.sh`), done;
      certified values in the table above, measures in `search/pgonly_*_exact.txt`.
- [x] `B40KM`, `PUREM`: the clique-only E1 instruments — no longer running as of 2026-09-13 (no `t4leaf` processes alive); their last values stand as recorded below and
      **not to be quoted as bounds on anything certifiable** (their clique rows are non-Helly).
      2026-09-12 20:00: `B40KM` `11.853128` (flat 20 iters), `PUREM` `11.852667` (creeping);
      `A0101M` **silent since 05:50** inside one call (1 thread at 100 %, 62 idle; its previous bad
      solve was `lp 2909s it296364`).  py-spy: `t4leaf.py:142` = `Hi.run` → HiGHS ipm+crossover
      with `time_limit 1e30` and no fallback under `--no-warm`, a crossover crawl.  **Killed
      2026-09-12 20:40** (SIGTERM, pid only; last checkpoint 05:37, LP `11.488236`): its number is a
      non-Helly-clique relaxation, not on a path to a proof.  `B40KM`, `PUREM` left running as cheap
      optionality for `helly-allmeet`.  If either is ever relaunched, cap the ipm `time_limit` in
      `Hi.run` and retry with `run_crossover off` on `kTimeLimit`.
- [x] `tasks/s13-writeup` (Opus, 2026-09-12 20:20–21:25, merged): README section, `VERIFICATION.md`
      entry, `verify.sh` runs `zmcheck` (1064 s / 8 thr, 2015 s / 4 thr) + `tests/rung2`, CI on every
      push, `docs/index.html`, `notes/s13-casefree.md`; full `verify.sh` exit 0, 25 VERIFIED;
      `search/S13_WRITEUP.md` §4 lists 9 doc inconsistencies — open ones below.
      CI (2026-09-13, dispatch-only on the research repo): green in 4 h 42 min, the sweep 13229 s on
      the runner (6.5× the local 4-thread time; `VERIFICATION.md`).  Decision pending: leave it
      dispatch-only, or split into fast-on-push / full-on-schedule before the public repo.
- [x] `tasks/helly-allmeet` (Opus, 2026-09-12 20:20–22:00, merged; `search/ALLMEET.md`, `allmeet.py`):
      **answered, no.**  Every dual-carrying clique row (`317/317`, `638/638`, 100 % of the clique dual)
      has an exactly certified positive-volume box `B_K` with `K ∪ B_K` pairwise meeting (median
      `vol 0.13–0.16`), and it recovers nothing: honest cost under `K ∪ B_K` is `1543` / `145`.  The
      clean statement: any sound rule's credited set is a clique containing `K`, hence `⊆ A(K)`, hence
      dominated pointwise by `meet` — **`honest ≥ 20.162 / 18.109` for every sound rule** on these
      duals.  Correction to the 09-12 narrative: the heavy rows are *robustly* pairwise-overlapping
      non-Helly families (least pair margin `+0.0004…+0.0008`, median `+0.7–0.95`, 0 grazing pairs at
      `1e-6`), not fans of grazing tangencies — the fan species exists (`E2Pg` rows 5/220/588) but
      carries little.  A robust non-Helly triple exists at `t = 4` (Farkas certificate), so the literal
      Helly statement in item 2 is false; the conjecture to state instead (`ALLMEET.md` §5): for any
      dual `w` of the `t = 4` packing LP on any finite pose set, with cliques credited by any rule whose
      set contains `K` and is pairwise meeting, `Θ(w) / min_S capture_w(S) ≥ 12`.
- [x] `zmcheck` reproduction in the main checkout: `runs/zmcheck_main_2026-09-12.log`;
      rejection suite `runs/rung2_rejection_2026-09-12.log`.

## Next

- [x] The all-meet box measurement — done, `search/ALLMEET.md`.
- [x] Bentz-incidence branching at `t = 4` — done, `search/BENTZ.md`; answer no.
- [x] Publication checklist for `s(13) = 4` — done 2026-09-12 (`search/S13_WRITEUP.md`).
- [ ] Doc follow-ups from `S13_WRITEUP.md` §4: `docs/index.html` §06 still says the LP ceiling is
      "about 3.95–3.96" (README: `s* ∈ [3.968616, 3.99)`); `ZEROMARGIN.md` §2 says four leaf tests
      (seven) and states the `u ∈ [0,½]` reduction without its symmetry precondition;
      `certificates/rung2/README.md` predates `zmcheck`/the Lean; zmcheck time quoted three ways
      (1017 / 1065 / "~20 min"); `lean-zeromargin.md` §6 item 2 not marked resolved.
- [ ] `../site/www/proofs.html` ("how Bentz ruled out thirteen"): add the case-free proof + link.
- [ ] The level-2 slot leaves at `t = 4` to a certified value (`PGCORN` configuration, `01010101`),
      to complete the table above — low priority; expected `>= 12` like its parents.
- [ ] `zeromargin.py` on Bentz's 16 unit-weight points: 3,516 boxes uncertified with no counterexample
      (`BENTZ.md` §3) — a checker-limits note for `ZEROMARGIN.md` (Theorem 1's edge), not a bug report yet.

## Compute efficiency

Done or in flight: sweep verifier (28×, `search/VERIFYSPEED.md`, a `t < 4` tool); the zero-margin
checker's per-box cost (vectorised, `ZEROMARGIN.md` §7); `cliquelever` restricted master
(`tasks/cl-master`, running).  Measured and not worth it (`LPSPEED.md`): warm dual simplex in any
form, PDLP, row sifting, LP threads (IPX is single-threaded).  The exact anchor scan (3.5 h) is
superseded by the QSTAB numbers.

## The consolation prize (`s(12) >= 3.98`) — off the critical path

- [ ] Finalize `1110` from `runs/branch_L1110f_*` and `k = 4` from `runs/branch_J16i_*` without
      `--max-iters`, 16 threads, detached; `--lam-hi 2` on `k = 4` (λ sits at the cap).  ~1 h/round
      each since the verifier speedups.  Then `verify_branch.sh`, README.

## Small

- [x] `lc_union2b` (corner-branch ceiling at `t = 4`, points only): LP mass pinned at 12.000000
      through polish 50, exact `M` never below 1.0007 — inconclusive, not a certified ceiling;
      `LEAF_CEILING.md` §5.2 stands as written.
- [ ] `leaf_ceiling.py exhaust` is complete-with-a-budget on real measures (`LEAF_CEILING.md` §7).
- [ ] Exhaustiveness of the leaf enumeration is Python, not Lean (`notes/branch-semantics.md` §5).
- [ ] `level2_capacity.py` is wrong at `t = 4` (nine squares fit in `[1,3]^2`); fix or retire.
- [ ] Exact re-checks of `certificates/branch/s12_t3.98_corner_k{1,2}.txt` (`xcheck.py`, ~8 h each).
- [ ] Prune the 29 agent worktrees (20 GB) once their `runs/` are no longer needed; the inputs the
      2026-09-11 tasks use are consolidated read-only in `runs/inputs-2026-09-11/`.

## Not worth it (decided)

- Interior (quadrant) branching on symmetric leaves: drop 0 by C4-averaging (`search/DEPTH.md`).
- Squeezing 3.975–3.978 from the corner level alone; more plus-region packing; SDP / Lovász-ϑ.
- Pure 3.99 to convergence (pure cannot go below 12 there; the clique gain is measured).
- Line-chord count cuts and the chord inequality as levers at `t = 4` (worth 0, tilt is immune).
- Any further count or incidence-pattern branch tree at `t = 4` (`BENTZ.md`: the fully pinned leaf is
  exactly 12; corner `k = 3` is `>= 12.03`; a count row cannot separate a square from a cloud of
  near-squares).
- Testing a new inequality on the axis-parallel case first: it is already `3` below 12 there.
- Level-2 branching by wall total or per wall (integral already).
