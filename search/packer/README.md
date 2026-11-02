# packer: searching for square packings

Packings of n unit squares in the smallest square, by search; results go through the store and, for records, through
`candidate.py` and `SUBMITTING.md`.  How it works: `DESIGN.md`.  What we have learned, with evidence: `CLAIMS.md`.
History (the lab log to 10-24): `archive/PACKER-2026-10.md`.

## Run something

| want | command |
|---|---|
| search a set of n | write `<plugin>/goals/NAME.yaml` (format: `goals/example.yaml`), `python3 -m pk.events add NAME`, then `python3 pk/tasks.py repeat NAME --budget CPU_H --prio low -- python3 -m pk.events work --campaign NAME --max-units 3`; `python3 -m pk.events report NAME` |
| add an idea | a stage (`DESIGN.md` "Adding an idea"), named as a move in a campaign |
| packings at n | `pk.py ls N`, `pk.py get ID -o FILE`, `pk.py frontier` |
| keep the store current | `pk.py ingest --apply` (register, watched repos, fresh public packings) |
| which generators feed records | `pk.py yield` |
| a record candidate | `candidate.py FILE --name NAME`, then `submission.py NAME` (`SUBMITTING.md`) |
| check a packing independently | `satcheck_f64.py ID ...`; exact: `dilate_cert.py`, `../exact/` |

## Files

- `pk/`: the library (store, sink, converge, engine, stages, moves, fill, schedules, events, families, parents, jobs,
  ingestion).
- `explore.py`: the within-n local search (a future source of the event loop).  `compose2.py`, `mondrian.py`, `rectlib.py`: piece
  composition.  `dgen.py`, `bandflip.py`, `harden.py`, `lcut.py`, `ladd.py`, `crop.py`: generators and moves still to be
  ported as stages.  `mcmin.py`, `hop.py`, `layout.py`, `gen.py`: helpers (file I/O, quench wrapper, full-line test,
  exact upright fill).
- `src/`: Rust (`fq` the quench / polish engine, `anneal`, `harden`).  Build: `cargo build --release`.
- `goals/example.yaml`: the campaign format.  Real campaign files (each is its run's pre-registration and, after the
  run, its result) are research results: they live in a plugin's `goals/`, found by name.
- Formats: text `n s` then `x y deg` per square (centres).
