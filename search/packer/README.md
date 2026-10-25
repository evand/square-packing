# packer: searching for square packings

Packings of n unit squares in the smallest square, by search; results go through the store and, for records, through
`candidate.py` and `SUBMITTING.md`.  How it works: `DESIGN.md`.  What we have learned, with evidence: `CLAIMS.md`.
History (the lab log to 10-24): `archive/PACKER-2026-10.md`.

## Run something

| want | command |
|---|---|
| explore a set of n | write `<plugin>/goals/NAME.yaml` (format: `goals/example.yaml`), then `python3 pk/jobs.py start NAME --procs P -- python3 pk.py goal run NAME`; `pk.py goal report NAME` |
| see arms and donor pools first | `pk.py goal plan NAME` |
| add an idea | a stage (`DESIGN.md` "Adding an idea"); the planner or a plugin's `ARM_TEMPLATES` turns it into arms |
| packings at n | `pk.py ls N`, `pk.py get ID -o FILE`, `pk.py frontier` |
| keep the store current | `pk.py ingest --apply` (register, watched repos, fresh public packings) |
| which generators feed records | `pk.py yield` |
| a record candidate | `candidate.py FILE --name NAME`, then `submission.py NAME` (`SUBMITTING.md`) |
| check a packing independently | `satcheck_f64.py ID ...`; exact: `dilate_cert.py`, `../exact/` |

## Files

- `pk/`: the library (store, sink, converge, engine, stages, moves, fill, schedules, goal, parents, jobs, ingestion).
- `explore.py`: the within-n local search the `walk` arms run.  `compose2.py`, `mondrian.py`, `rectlib.py`: piece
  composition.  `dgen.py`, `bandflip.py`, `harden.py`, `lcut.py`, `ladd.py`, `crop.py`: generators and moves still to be
  ported as stages.  `mcmin.py`, `hop.py`, `layout.py`, `gen.py`: helpers (file I/O, quench wrapper, full-line test,
  exact upright fill).
- `src/`: Rust (`fq` the quench / polish engine, `anneal`, `harden`).  Build: `cargo build --release`.
- `goals/example.yaml`: the goal format.  Real goal files (each is its run's pre-registration and, after the run, its
  result) are research results: they live in a plugin's `goals/`, found by name.
- Formats: text `n s` then `x y deg` per square (centres).
