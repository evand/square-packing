#!/usr/bin/env bash
# s60_cert_runs.sh -- the certified runs shipped in certificates/s60/ (s(60) = 8), with manifests.
# The s(60) twin of search/s45_cert_runs.sh.  Run from anywhere (detached is fine):
#   setsid nohup bash search/s60_cert_runs.sh [zmx2|zm_mixed|import RUNDIR|all] > runs/s60cert.out 2>&1 < /dev/null &
# Env: CORES (taskset list, default 6-11), NPROC (default 6).
#
#   zmx2      verify2 zmx2 on the cover: --d4 (6,400 roots) and --full (no symmetry, 51,200 roots); ~3 CPU-min.
#             Writes certificates/s60/zmx2_{d4,full}/{roots.log,run.out,manifest.txt}.
#   zm_mixed  search/zm_mixed.py --d4 --cert-mode on the cover (102,400 roots, ~19.6 CPU-h, ~3.3 h on 6 processes).
#             Resumable (roots.jsonl, header-checked).  Writes certificates/s60/zm_mixed_d4/{roots.jsonl,run.log,
#             manifest.json} and checker/ (the exact checker files that ran).
#   import RUNDIR
#             how the shipped zm_mixed_d4/ was made: RUNDIR is an earlier complete run of the same checker files on
#             the same bytes with the same settings, under another path (runs/zm_mixed_s60/, made by
#             `zm_mixed.py cert runs/s60_mixed_candidate_8.txt ... --resume runs/zm_mixed_s60/roots.jsonl`).  Its
#             records are copied verbatim, its header with only the `input` path changed to the bundled cover;
#             then zm_mixed.py itself is run with --resume on them: it recomputes the header, refuses unless shas,
#             settings and total match, finds all 102,400 roots done and writes the manifest.  run.log = the
#             original run's log followed by the resume's log.
set -eu
cd "$(dirname "$0")/.."            # s12/
WHAT=${1:-all}
CORES=${CORES:-6-11}
NPROC=${NPROC:-6}
B=certificates/s60
C=$B/s60_mixed_cover_8.txt
ZMARGS="--d4 --cert-mode --disj --chain-from 0 --depth 24 --pitch 1/20 --ubins 16"

if [ "$WHAT" = zmx2 ] || [ "$WHAT" = all ]; then
  (cd verify2 && cargo build --release --bin zmx2 2>&1 | tail -1)
  Z=verify2/target/release/zmx2
  for mode in d4 full; do
    O=$B/zmx2_$mode; mkdir -p $O
    rm -f $O/roots.log            # a fresh run, not a resumed one
    t0=$(date +%s)
    taskset -c $CORES $Z cert $C --$mode --threads $NPROC --log $O/roots.log > $O/run.out 2>&1 || true
    t1=$(date +%s)
    {
      echo "# zmx2 run manifest (search/s60_cert_runs.sh; certificates/s60/README.md)"
      echo "date:            $(date -u +%Y-%m-%dT%H:%M:%SZ)   wall $((t1 - t0)) s"
      echo "host:            $(uname -m), taskset -c $CORES, $NPROC threads"
      echo "git HEAD:        $(git rev-parse HEAD) ($(git status --porcelain -- verify2/src/bin/zmx2.rs | wc -l) local change(s) to zmx2.rs)"
      echo "checker source:  verify2/src/bin/zmx2.rs sha256 $(sha256sum verify2/src/bin/zmx2.rs | cut -d' ' -f1)"
      echo "checker binary:  $Z sha256 $(sha256sum $Z | cut -d' ' -f1)"
      echo "rustc:           $(rustc --version)"
      echo "input:           $C sha256 $(sha256sum $C | cut -d' ' -f1)"
      echo "input summary:   $($Z info $C | tr '\n' ' ')"
      echo "D4 check:        $($Z d4 $C)"
      echo "command:         taskset -c $CORES $Z cert $C --$mode --threads $NPROC --log $O/roots.log"
      echo "log header:      $(head -1 $O/roots.log)"
      echo "log sha256:      $(sha256sum $O/roots.log | cut -d' ' -f1)  ($(grep -c '^ROOT ' $O/roots.log) ROOT lines, $(grep -c '^UNCERT ' $O/roots.log || true) UNCERT lines)"
      grep -E '^(roots:|done in|VERIFIED|NOT VERIFIED|INCOMPLETE|REGION)' $O/run.out | sed 's/^/result:          /'
    } > $O/manifest.txt
    cat $O/manifest.txt
  done
fi

if [ "$WHAT" = zm_mixed ] || [ "$WHAT" = all ]; then
  O=$B/zm_mixed_d4; mkdir -p $O/checker
  for f in zm_mixed.py mixed_cover.py zeromargin.py; do cp search/$f $O/checker/$f; done
  taskset -c $CORES python3 search/zm_mixed.py cert $C $ZMARGS --nproc $NPROC --progress 1000 \
      --resume $O/roots.jsonl --manifest $O/manifest.json >> $O/run.log 2>&1 || true
  for f in zm_mixed.py mixed_cover.py zeromargin.py; do
    cmp -s search/$f $O/checker/$f || { echo "search/$f changed during the run"; exit 1; }
  done
  tail -5 $O/run.log
fi

if [ "$WHAT" = import ]; then
  R=${2:?usage: $0 import RUNDIR}
  O=$B/zm_mixed_d4
  [ -e $O/roots.jsonl ] && { echo "$O/roots.jsonl exists; remove it first"; exit 1; }
  mkdir -p $O/checker
  for f in zm_mixed.py mixed_cover.py zeromargin.py; do cp search/$f $O/checker/$f; done
  python3 - "$R/roots.jsonl" $O/roots.jsonl $C <<'EOF'
import json, sys, hashlib
src, dst, cover = sys.argv[1:]
lines = open(src).read().splitlines(keepends=True)
h = json.loads(lines[0])
assert h['kind'] == 'zm_mixed cert header', "no header"
assert h['sha256']['input'] == hashlib.sha256(open(cover, 'rb').read()).hexdigest(), "different cover bytes"
print(f"import: {src}: {len(lines) - 1} records, header input {h['input']!r} -> {cover!r}")
h['input'] = cover
with open(dst, 'w') as fh:
    fh.write(json.dumps(h) + "\n")
    fh.writelines(lines[1:])
EOF
  { cat $R/run.log
    echo "=== search/s60_cert_runs.sh import $R ($(date -u +%Y-%m-%dT%H:%M:%SZ)): records copied into $O/roots.jsonl" \
         "(header input path -> $C), then zm_mixed.py --resume on them ==="
  } > $O/run.log
  taskset -c $CORES python3 search/zm_mixed.py cert $C $ZMARGS --nproc 1 --progress 1000 \
      --resume $O/roots.jsonl --manifest $O/manifest.json >> $O/run.log 2>&1 || true
  for f in zm_mixed.py mixed_cover.py zeromargin.py; do
    cmp -s search/$f $O/checker/$f || { echo "search/$f changed during the run"; exit 1; }
  done
  tail -6 $O/run.log
fi
