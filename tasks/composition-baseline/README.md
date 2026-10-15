# Composition baseline: family A / B record events at n <= 324 (built 2026-10-14)

`events.yaml`: one entry per record event in `search/packer/bench/events.json` with n <= 324, family A (rectangle-library
composition / graft) or B (soft-to-rigid re-hardening), gain >= 1e-8: **45 events** (41 A, 4 B; total gain 0.186).
Each entry has the cutoff T (event date, UTC, plus the local-time equivalent), old/new, what the post says was done (with
a quote), parents with store lookups, piece details, our past status with a file/section reference, and the component our
pipeline would lack.

Sources (read-only): jlevy/squares #401, #422, #481, #500, #528, #562, #575 (`gh issue view`, bodies + comments);
`search/packer/runs/ext/itsnaka-squish-certs-{d45669b,aeacc6c}/` and `dostashev-square-packing-b909e1b/` READMEs;
`MONDRIAN.md`, `EVENTS-2026-10-12.md`, `EVENTS-2026-10-15.md`, `REVIEW-2026-10-11.md` §4, `PACKER.md`; the store
(`runs/store.sqlite`, opened read-only).  Build script: session scratchpad, not kept.  Re-run its lookups by hand from
the `store_ids` / `kingbird_ref` fields.

Store lookup rule: a parent matches store rows with |s - side| <= 1e-9, or a Kingbird row by its SVG name.
`in_store_before_T` compares `packing.added` (local time, UTC-7) with T.  The store was created 2026-10-09 20:15 local
(10-10 03:15Z) and Kingbird revisions were imported 10-10 23:35 local.  So "not in store before T" is often an artifact
of when the store was built.  The data-gap list below separates the cases.

## Counts

**By method sub-type**

| sub-type | events |
|---|---|
| Mondrian v1 (two pieces, opposite corners, uprights, squeeze) #481 | 12 |
| Mondrian v1, same-source double pieces (153, 270, 307) | 3 |
| Mondrian v2, narrowed pieces #500 (incl. 305 same-source double) | 8 |
| Mondrian v2, interlock #500 (176) | 1 |
| graft (#401 update 123, 126, 239; #422 88, 199, 207, 236, 302) | 8 |
| carve (263 #401) and nearby search after carve (263 #422) | 2 |
| rect library: new searched block(s) + aligned (175, 176 repo, 204, 303) | 4 |
| rect library: new block + crop (232, 264 #528) | 2 |
| rect library: cut-out (129 #562) | 1 |
| soft-to-rigid V merge: two-band packing (211, 241) / composition (176) | 3 |
| soft-to-rigid relax of a composition, no V (304 #575 edit) | 1 |

**By our_status**

| our_status | events | gain | which |
|---|---|---|---|
| REPLICATED (rebuilt at or below their side) | 7 | 0.044 | 129 #562, 175 #528, 259 / 263 / 292 #481, 176 / 211 #575 |
| BEAT-ANTE (below frontier(T), above theirs) | 2 | 0.012 | 232 #481 (15.76872), 241 #575 (15.98798; uncertified) |
| no (tried, did not reach frontier(T)) | 22 | 0.075 | 11 #481 Mondrian v1 (MONDRIAN.md trimmed pass); 304 #575 (harden / fq relax); 10 graft / carve / nearby events with generic search only (hop0, TT1); the graft operator was never tried |
| not attempted | 14 | 0.055 | all 9 #500 events (REVIEW §4: "Not replicated by us"); 176 repo, 204, 303 #562; 232, 264 #528 |

**By missing component** (primary component per event; `events.yaml` lists a secondary for #500, seam repair)

| missing component | events | gain |
|---|---|---|
| squeeze / seam search after composing (our compose + fq / short explorer lands 1-5 squares short or in other basins) | 12 (+8 as secondary) | 0.057 |
| graft / carve into a host record, then nearby search | 9 | 0.029 |
| narrowing a piece in its own rectangle before composing (Mondrian v2) | 8 | 0.017 |
| interlock (piece minus a 2 x 3 corner block, other piece reaching in, + add-a-square search) | 1 | 0.010 |
| from-scratch rectangle block search for cell (62, 6) [303 #562] and (55, 8) / (68, 7) [304 #575 start layout] | 2 | 0.003 |
| soft-to-rigid V merge at 241 (works at 211 and 176, not 241) | 1 | 0.001 |
| nearby search from a public parent (263 #422) | 1 | 3e-5 |
| none known: replicated | 7 | 0.044 |
| none known (guess): every component exists but the event was never rebuilt (176 repo, 204, 232 #528, 264) | 4 | 0.025 |

Reading (guess): the squeeze after composing blocks the most events: all the #481 misses, and it is the second step of
every #500 event.  Grafting is the largest untried operator by event count.  Narrowing comes next and is the defining
step of Mondrian v2.  Rectangle-block search is not the bottleneck at these events: the (38, 2) and (44, 3) blocks
were reached by our rectsearch (the (44, 3) one only in the 10-13 bandit smoke test, store pieces 2928-2994).  The one
cell that fell short is (62, 6), our best 5.919-5.928 vs 5.90996.

## Data gaps: parents NOT in the store before T

Events with every parent in the store before T: **129 #562, 176 #500, 211 #575, 241 #575** (4 of 45).

1. **Not public before T** (no snapshot could have had them):
   - SQUISH unfiled sources released with the issue: 131 and 209 (s(155); 209 also s(240)), 305 #481 (s(306)); 181
     (s(154) a + b), 238 (s(305) 17.950859 + s(154) b), 239 #500 (s(240) b), 271 (s(305)), 305 #500 (s(272), same post).
   - Same-post parent: 207, 236 and 302 #422 (SQUISH's new s(88), posted in #422).
   - New pieces first published in the event: 38-square (h+2) x h block (175, 232, 264 #528), 44-square (h+3) x h block
     (204, 176 repo), 62-square (h+6) x h block (303 #562).  304 #575: the composer layout (17.939518) and its
     (55, 8) / (68, 7) pieces were never published.
2. **Public before T, not in our local mirror either**: Kingbird former records square-152_r2 (153, 236 and 302 #481),
   square-202_r3 (263 #481), square-241_r3 (307).  These are absent from the 2026-09-30 snapshot and were fetched 10-10
   (MONDRIAN.md, Replication harness).
3. **Public before T, store too late** (an artifact of when the store was built; the rows exist now):
   - Every #401 / #422 / #481 parent: register, Couzo / Xu repos, and Kingbird files that were in the 09-30 mirror
     (square-5, -11_r1a, -37, -41, -50_r5, -88, -123, -126, -239, -263).
   - 270, 272 and 304 #500: Kingbird square-19_r1 / -182_r1 / -11_r1a were imported 10-11 06:35Z, about 1 h after #500.
   - 176 repo (Ostashev s(175) #528, synced 10-13 07:19Z) and 176 #575 (Ostashev's 13.808663 composition, synced
     10-13 07:19Z, 6.6 h after #575).
4. **Host packing side absent from the store**: these hosts are guessed, since #422 does not name them.  199 host (the
   10-05 value 14.618988956907) and 302 host (17.885993892531): no store row within 1e-9.  The 207 host is only in our
   10-14 runs, the 236 host only in our 10-09 runs.

## Caveats

- "our_status" is taken from what MONDRIAN.md / PACKER.md record.  Runs not written up there (e.g. rectsearch `bsmoke`)
  are noted only where the store's piece table shows them.
- REPLICATED for 175 and 129 used searches run after the event: the piece search for 175 had #528's h as its success
  target, and the 129 rebuild used crops of the #481 s(302) at the 10-11 05:51 cutoff.  211's T1 hit (14.9936458) is
  below #575's posted value but above the later frontier.
- Graft hosts for #422 (199, 207, 236, 302) are not named in the issue.  The host in `events.yaml` is a guess: the
  record at T.
