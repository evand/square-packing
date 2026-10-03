# t3-chain: rehearse the chain architecture on `s(6) = 3` — does it close where the answer is known?  (2026-09-21)

**Why.**  The candidate theorem for `s(12) = 4` (`notes/proof-architecture.md` §0a item 14, `search/BANDCUT_K.md` §6) is
"every configuration has a tight wall-to-wall chain of `T` on which MT gives `<= 0`".  As stated it is not a certificate:
in the chain inequality of `BANDCUT_K.md` §1.2 the one positive term is `sum_s sin(tau_s) Dy_s`, so the chain kills only
when its **rise is bounded** (`R <= R*(T,t)`, about `1.01–1.13` at `T = 3`), and the rise is bounded by *transverse*
chains — which is what every dual we have read looks like (`BANDCUT_K.md` §5: a chain of `T` at weight `~1/(T+1)` per row
plus a transverse chain at a tenth of it; `notes/status.md`: the pinwheel dual is a staircase chain plus transverse links at
`t/3`, `t^2/3`).  So the lemma to prove is **H-shaped**: *a wall-to-wall chain of `T` plus the transverse chains that pin
its rise, with total value `<= 0`.*  And it says nothing yet about the far field (all squares well tilted, generic normals),
where `sin(tau) Dy` is large and nothing in the repo is a mechanism (`search/ARCH_FARFIELD.md`: the capped LP fails even at `T = 3`).
Nothing in this architecture has ever closed *any* case.  `T = 3` is a theorem (Kearney–Shiu 2002), six squares, chains of
three, hole `k in {3, 4}`: if the H-lemma cannot be proved there it will not be at `T = 4`, and if it can, the proof shows
how the far field is handled.

**Read first.**  `notes/proof-architecture.md` §0a; `search/BANDCUT_K.md` §1.2, §2, §5, §6; `notes/bandcut-cost.md` §1;
`search/S6_LOCAL.md` §3, §5; `search/S6_SKELETON.md` §0, §2, §3.1; `notes/proof-anatomy.md` §5.1–5.2 (Stromquist, Kearney–Shiu).
Instruments that exist: `search/s6skel.py` (exact `Decider`), `search/s6local.py`, `search/s6exact.py`, `search/s6cube.py`,
`search/bandcut_k.py` (chain detector, `cert` dual reader), `search/chains.py`.  Use the wall-carrying `delta` of `BANDCUT_K.md` §1.1.

**Do.**  `T = 3`, `n = 6` only.
1. **State the H-lemma exactly**, as a restricted dual: the Farkas multipliers of `S6_SKELETON.md` §3.1 with support confined
   to one wall-to-wall chain of three (links on any edge normals, staircases allowed) plus transverse wall-to-wall chains
   through its members.  Give the inequality it yields in closed form (chain row + how each transverse chain bounds a `Dy_s`).
2. **Test it as a statement over the whole angle space**, not just the small-angle cube: near `Z`, the hole, the far field
   (`k <= 2`, including all-`45°`, uniform large tilt, random generic tilts, mixed signs).  At each sampled optimum compare
   `delta*` with the best H-restricted dual bound.  Where is H `<= 0`?  Where it is not (or no chain exists), read the true
   dual: what shape is it?  That shape is the far-field mechanism we are missing — classify it, do not just report failure.
   Everything sampled is a measurement; report hit rates and under-explored cells.
3. **Translate Kearney–Shiu into this language** (`proof-anatomy.md` §5.2; the paper is Electron. J. Combin. 9 (2002) R14 if
   you need the text).  Their two 7-point lattices + "covers X ⇒ cannot avoid Y" steps: which step is the chain, which
   bounds the rise, and how do they dispose of all-tilted configurations without any chain?  Does any step have an obvious
   `T = 4` analogue or an obvious reason to fail at `T = 4`?
4. **Attempt the proof.**  A written proof of `s(6) = 3` through S0/S1 + H-lemma + whatever far-field statement step 2–3
   suggests.  Pen and paper, exact; a finite case split is fine if each case closes with an explicit dual whose weights are
   trig-rational and valid on a whole semi-algebraic piece (the closed forms of `S6_LOCAL.md` §3 are the model — they are
   global in `t`, not asymptotic).  If it does not close, say precisely which region of angle space is left and what
   statement would close it.  A proof with one clearly-stated open sub-lemma is a good outcome; a vague "mostly works" is not.
5. **Filter.**  State the lemmas `T`-generically.  They must not prove `s(5) = 3` (five squares in side 3; `s(5) = 2.7071`),
   and the set must contain one that is false at `T = 11`: say which, and where `3` (resp. `T`) enters as a number.

Label everything **[proved]** / **[measured]** / **[heuristic]** / **[guess]**.  An adversarial attitude to your own lemmas
is wanted: last session four of the coordinator's claims and one agent's were wrong (`notes/review-2026-09-20c.md`, "Errors caught").

**Deliverable.**  `notes/t3-chain.md` (verdict up front), scripts `search/t3_chain*.py`, runs `runs/t3_chain_*`.
Do not edit `TODO.md`, `notes/status.md`, `notes/proof-architecture.md`; do not commit.

**Working style.**  Work in small steps; write the note section by section (append), never in one giant message; keep
derivations on disk.  Anything over 10 minutes: launch detached (`setsid nohup`), 16 threads is the useful width, never
`pkill`, leave no `until … sleep` / `tail -f` watchers behind (a `pgrep -f <name>` wait loop matches itself and never exits).
