# Ceiling tests and lower-bound targets, n = 17–20

Lower bounds (site data): s(17) ≥ 4.66044 (R071), s(18) ≥ 4.679, s(19) ≥ 4.8, s(20) ≥ 4.85.
Best packings: s(17) ≤ 4.67553 (Bidwell), s(18) ≤ (7+√7)/2 ≈ 4.82288 (Hämäläinen), s(19) ≤ 3+4√2/3 ≈ 4.88562
(Wainwright), s(20) ≤ 5.  Milestone: the strict chain s(17) < s(18) < s(19) < s(20) needs
**s(19) > (7+√7)/2** and **s(20) > 3+4√2/3** (s(18) > s(17) already holds).

## Questions (per (n, t))
- **Cover side** (`search/line_cover.py loop --s t`, LINE_COVER.md; points + grid-line densities, closed semantics):
  the LP value of the best closed cover of `[0,t]²`.  Value < n ⇒ a certificate of s(n) > t is plausible
  (exactness overhead was 0.2 % at k = 5 after closing loops, 1.5 % in phase A, S60_COVER §4).
- **Ceiling side** (`search/cover4_cg.py cg ... --t`, exact restricted master via `dual_exact.py`): an exact
  fractional packing of mass ≥ n at side t proves **no** pure closed cover gives s(n) > t.
Targets: (17, 4.6604) [the open ceiling test in TODO], (17, 4.6755); (19, 4.8229), (19, 4.8856);
(20, 4.8856), (20, 5) [the last is known infeasible for covers, L(5) ≈ 4.25: sanity check only]; (18, 4.8229) if
cheap.  Report the bracket [ceiling, cover LP] for each.

## Then
If the cover LP at a milestone side (19 at 4.8229, 20 at 4.8856) is below n with room (≳ 0.5 % of n), build the
certificate: exact closing loop, verify with zm_mixed (and zmx2 if time), budget below.  Use the near-tight-row
lever (`notes/jlevy-s17-techniques.md` §2 C: verifier near-tight cells as exact LP rows).  If not, stop and report.

## Rules
CPU: no pinning; keep total busy processes ≲ 16; budget ~4 CPU-h for the LPs, ask before going past ~15 CPU-h total
(certificates included).  Write-up `search/CEILINGS_17_20.md`; runs `runs/ceil_*`; commit only your own files to
public main, don't push.  Don't post externally.
