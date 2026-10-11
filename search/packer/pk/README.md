# pk: packing toolkit (10-09; design and rationale: ../DESIGN.md)

| module | what |
|---|---|
| `packing.py` | `Packing` (n, s, sq = x, y, deg); `read` / `parse` any format (ours, exact `.cert`, Couzo, Ellsworth, SQUISH `cert.json`, register Witness v2: center-angle / center-basis / corners); writers; `max_pen` (f64 SAT); `grid_lines`; D4 + relabelling canonical hash |
| `store.py` | sqlite store `runs/store.sqlite`: `packing` (dedupe by canonical hash, status raw < screened < polished < certified / not-min / infeasible), `sighting` (provenance, many per packing), `frontier` (best known per n and source; pending flag), `basin` (certification keyed by side to 1e-9) |
| `sources.py` | importers: files, register (sparse clone), pending issues (`pending.yaml`: repo @ commit, tarball into `runs/ext/`, read-only, data only), explorer archives (parent links + exact.json statuses) |
| `engine.py` | `fq serve` client (one persistent child per process); `QOpt` = every fq knob; `SCREEN`, `POLISH` presets |
| `moves.py` | move registry with declared parameter spaces; params recorded per proposal; `ARMS` = 10-08 explorer kinds as presets |
| `pipeline.py` | `evaluate(proposal, Policy, Known, best, seed)`: screen -> grid / discard / return -> (staged) polish |
| `corpus.py` | replay bench: `build` (frozen proposals from store parents x arms), `replay VARIANT`, `certify` (exactsolve, shared cache), `report [--by arm|param:X|loosen|parent_gap]` |
| `kingbird.py` | Ellsworth's catalogue: `pk.py sync-kingbird` crawls all 10 catalogue pages (incl. `compared2`/`compared3`, the only pages linking most `_r*` revisions), dated snapshot `runs/kingbird/<date>/` (pages, svg, `manifest.json` with sha256, finder / date / full header comment, flags), site parser (`site/tools/parse_svg.py`), store import (source `kingbird`, ref `<file>@<date>`, tags `kingbird:current|revision|variant`); revisions reused from older snapshots; `latest(name)` |
| `jobs.py` | background jobs without command-line matching: `start NAME -- CMD` (own process group, pid recorded, forced log + meta with git commit), `status` (group CPU, last log line, load vs ~16), `stop NAME` (group kill), `log NAME`; records in `runs/jobs/<NAME>/` |

CLI: `../pk.py` (sync-register, sync-pending, import, import-run, frontier, ls, get, info, tag, export).
Tests: `python pk/tests/test_packing.py`.
