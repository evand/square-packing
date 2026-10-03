# rung2-xcheck: an independent re-implementation of the zero-margin checker (2026-09-12)

**Why.**  `certificates/rung2/s13_closed_cover_4.txt` (3,621 weighted points of `[0,4]²`, total
`2591194431/200000000 = 12.955972 < 13`) is a case-free proof of `s(13) = 4` — *if*
`search/zeromargin.py` is correct.  It is a single 1,015-line Python implementation with no second
check of its two load-bearing primitives (`ADM`, `CHAIN`) and no Lean; the `s(12)` bound, by
contrast, has Rust + `xcheck.py` + Lean.  Before `s(13) = 4` is stated anywhere public it needs an
independent checker, and the same code becomes the second implementation of V1.

**Semantics** (do not change): closed unit squares, closed containment, a point on `∂Q` counts;
a pose is admissible iff its closed square lies in `[0,4]²` (`search/ZEROMARGIN.md` §1).  The
claim to check: **every admissible pose captures total weight `≥ 1`.**

**Rules of independence.**  Work from the *definitions and lemmas* in `search/RUNG2.md`
(§1 notation, §3 Lemmas A–C for `ADM`, §4.6 `clip_bin`, §6 Lemmas E–H and the down-set rule for
`CHAIN`, §7 what is exact) and `search/ZEROMARGIN.md` §1–2, not from `zeromargin.py`.  Read the
Python only to resolve an ambiguity in the note, and record every such place in your write-up
(those are documentation bugs to fix).  You may design a *different* subdivision, splitting rule or
primitive set as long as every certified box is certified by a lemma you state and prove; you do
not have to reproduce the 16,872-box census.  Rust preferred (`verify/` shows the house style:
exact integer / rational arithmetic, no floats on the load-bearing path; `i128` or a bignum crate —
the certificate has `D = 1000`, `W = 10⁹`, and the polynomials of Lemma B are degree ≤ 4 in `u`
with rational coefficients, so plan the bit widths or use `num-rational`/`rug`).  Exact Python
`Fraction` is acceptable if Rust would not finish in the budget, but say why.

**Definition of done.**
1. Your checker reports **0 uncertified boxes** on the shipped certificate over the full admissible
   domain, deterministically, with its own leaf census and wall time.
2. It **refuses** mutated certificates: (a) one point's weight reduced by 1 %; (b) one point deleted;
   (c) the whole set scaled by `0.995`; (d) a point moved by `0.01`; (e) the `closed4_best_x103.txt`
   cover, which is *invalid* at `(3/2, 1461/2000, u = 1/40000)` by `0.9703` (`RUNG2.md` §4.3) — your
   checker must either report that box uncertified or, better, exhibit the violating pose.  Put
   these under `tests/rung2/` in the style of `tests/`.
3. A note listing every lemma you rely on, with a proof or a pointer to the proof in `RUNG2.md`,
   and every place the note was ambiguous or wrong.

**Inputs.**  `certificates/rung2/s13_closed_cover_4.txt` (format: line 1 `m sym`, line 2 `D`,
line 3 `W`, line 4 count, then `X Y w` with coordinates `X/D`, weights `w/W` — confirm against
`certificates/FORMAT.md` and the file header), and read-only
`/home/evand/math/square-packing/s12/runs/inputs-2026-09-11/closed4_best_x103.txt`.

**Budget.**  Up to 8 processes; launch anything over 10 minutes detached (`setsid nohup … &`,
`< /dev/null`); never `pkill`; pids 1109693–1109695 are other runs — do not touch.  Report within
~4 h even if not done, with the exact state (which boxes remain, and why).

**Deliverables.**  `verify2/` (or `verify/rung2/`, a separate crate; do not modify `verify/`'s
existing binaries or `verify.sh`), `tests/rung2/`, `search/RUNG2_XCHECK.md`, commit on your
worktree branch.  Do not edit `TODO.md`, `README.md`, `lean/`, `certificates/`, or other tasks'
files.
