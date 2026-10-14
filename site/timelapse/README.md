# Records time-lapse (n = 1..324)

A video of the best known packings of n unit squares in a square, from the grids before 1979 to the
register of October 2026, on jlevy's Atlas triangle with jlevy's colours.  Each record that falls flashes
its cell; the caption keeps the latest records.

Outputs (gitignored, in `out/`): `timelapse.mp4` (3840x2160, 30 fps, 101 s, 31 MB), `timelapse-1080p.mp4` (11 MB), `final.png`
(the last frame), `history.json`, `views/` (one regularised drawing per state), `regview.log`.

## Rebuild

```sh
cd site/timelapse
python3 history.py [--fetch-issue-dates]   # dated states per n -> out/history.json (seconds)
python3 regview.py --procs 3               # regularised views, cached by content -> out/views/, out/views.json
python3 render.py                          # -> out/timelapse.mp4 + out/final.png  (--still-only: final.png only)
ffmpeg -i out/timelapse.mp4 -vf scale=1920:1080:flags=lanczos -c:v libx264 -crf 20 -pix_fmt yuv420p out/timelapse-1080p.mp4
python3 render.py --strip 155 237          # debug: every state of n, source (top) and view (bottom)
```

Needs the packing store (`search/packer/runs/store.sqlite`), the register sparse clone
(`search/packer/runs/jlevy`), `site/data/credits.json`, `site/www/data/{hist.json,p/}`, `search/exact/batch/`,
`search/packer/bench/events.json`; imports `search/regularize/{reg,jstyle}.py` (not modified).

## Data and dating (`history.py`)

Per n, a list of states with strictly falling side:

1. **Grid** from the start (before 1979 every cell is a grid; n = k^2 stays one).
2. **Catalogue era** (to 2026-10-05): the upper-bound holders of `site/data/credits.json` (Ellsworth's catalogue prose
   continuing Friedman's survey, plus `site/data/records2026.json` for the September 2026 postings), with the
   drawings `site/www/data/hist.json` assigns (catalogue SVG parses or superseded records from the store).  Dates are
   what the prose gives: (year, month) -> the 15th; year only -> 1 July; "early/late" -> February/November; a day
   where the posting has one (de Winter 09-16, Casson 09-23, Couzo 09-27 and 10-03, ...).
3. **10-05 register** where it is below that chain (22 n: Couzo's September packings that credits does not list
   separately).  Dated by the register case file's priority-note claim at that side when the finder's own claim is
   there (155, 207), else 2026-09-27 (Couzo's repo report; approximate, flagged `approx_date`).  Drawing: the frozen
   batch input `search/exact/batch/inputs/n-N.txt`.
4. **Record events** 10-07 .. 10-13 00:42 UTC: `search/packer/bench/events.json` (136 events at n <= 324, issue
   post times), drawing = the store packing at the event's side (all 136 found; 134 used: at 228 and 306 the
   catalogue chain's Couzo 10-03 packing is already below the event's side).
5. **Later public postings** in the store (register / pending-issue sightings) below all of that: dated by the issue's
   `createdAt` (`issue_dates.json`, fetched read-only with `gh issue view`).  As of 10-14: 154 (#586), 207 (#573),
   211, 241, 302 (#583), 237 (#593), 304 (#575).
6. **Envelope** s(n) <= s(n+1): where n+1's best known packing (own or itself derived) is below n's own record, n's
   cell shows n+1's packing minus one square (fewest contacts first, then axis-parallel, then nearest the top-right
   corner), dated with n+1's record.  112 derived states in 66 n (e.g. 102, 103 from Göbel's 104 in 1979).  They
   change the picture but are not listed as records in the caption.

Left out: precision-only claims (< 1e-8), our unposted results, n > 324.

Known gaps:
- Catalogue holders without a side whose drawing is a later holder's (the catalogue redrew the entry in place, e.g.
  Cleemann's 272 of 1991-98 under Ellsworth's 2026 drawing) take the largest-side Kingbird revision in the store
  between their neighbours' sides; if there is none they are skipped (only 68: Allen Chang 2026-09).  Holders whose
  side equals the previous one (refinements credited as records) are skipped too.  `history.json` `notes` lists all.
- Year-only dates (1997, 2002, ...) are placed mid-year; month-only dates on the 15th.
- Most of the 22 register-only September records are dated 09-27 by assumption (Couzo's repo commits run 09-23..10-03).
- The envelope's removed square is a heuristic choice, not anyone's drawing.

## Regularisation (`regview.py`)

A display view only; the side never changes.  Per state:
1. angles equal to 1e-6 deg merged to one exact value (float witnesses); near-axis squares (|theta| < 1 deg)
   straightened to 0 together if the translation LP (`reg.place`) stays feasible at S + 1e-9;
2. orientation: the last state of each n keeps its source orientation (jlevy's atlas = source drawings); walking back,
   each earlier state takes the dihedral image that best matches the next one (Hungarian assignment on centres in
   units of the side, plus an angle term).  Grids stay canonical (holes top right);
3. placement at S + 1e-9: tilted clusters pulled together, straight tilted chains made exact where feasible, axis squares
   under gravity to the bottom left (REGULARIZE.md's "pretty" rules 2-3 via `reg.pretty_objective` / `reg.place`,
   without the wall pull, which split grid blocks, and without imposed symmetry); fallback plain gravity; accepted if
   every gap is >= -1e-7, else the angle-cleaned source is drawn;
4. a final state whose side matches `search/regularize/lists/best.json` is drawn from that certified regularised list.

## Rendering (`render.py`)

Triangle: row k = ceil(sqrt n) holds n = (k-1)^2+1 .. k^2, right-aligned (column c = k^2 - n from the right edge), with a
small gap before each row's grid suffix (as of the last frame), as in jlevy's poster (`square-packings-324-20261008.pdf`).
Colours: `search/regularize/jstyle.py` (port of jlevy's `render/color.py`: hue = tilt class, 0 deg and 45 deg fixed,
the others by class size within the packing; shade = full-side contacts, walls included).  Because hue order for the
other tilts is per packing (jlevy's rule), a band's colours can change when its record changes.

Clock: piecewise linear, 1979 -> 2024 in 13 s, 2024 -> 2026-09 in 20 s, Sept 2026 in 16 s, 10-05 -> 10-15 in 48 s;
steps (states sharing a timestamp) at least 0.27 s apart; 3 s intro, 6 s hold at the end.  A changed cell cross-fades
(0.3 s) and carries a crimson frame fading over 2 s.

## Counts (first build, 2026-10-14)

881 states: 324 grids, 282 catalogue holders, 22 register-only (10-05), 134 events, 7 later postings, 112 derived;
about 560 timestamps (steps).  Regularisation: ~9 CPU-min for all views (3 processes); 47 placements hit the LP
iteration cap ("not converged": still a checked, overlap-free placement), none fell back to the source.  Render ~7
CPU-min (Python compositing + libx264 at 3 threads).
