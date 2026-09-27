#!/usr/bin/env bash
# zmx2_run.sh -- the official zmx2 runs on the s(21) candidate, resumable, with a manifest.
# Run from s12/ (detached is fine):
#   setsid nohup bash search/zmx2_run.sh > runs/zmx2/run.out 2>&1 < /dev/null &
# Each sweep keeps a per-root log (runs/zmx2/*.log); rerunning skips finished roots.
# Writes search/zmx2_manifest.txt (checker sha, settings, roots, result).
set -eu
cd "$(dirname "$0")/.."
C=${C:-runs/line-cover_m5_candidate_x1003.txt}
TH=${TH:-10}
CORES=${CORES:-0-9}
OUT=runs/zmx2; mkdir -p $OUT
MAN=search/zmx2_manifest.txt
(cd verify2 && cargo build --release --bin zmx2 2>&1 | tail -1)
Z=verify2/target/release/zmx2
sweep() { # name, flags...
  local name=$1; shift
  local t0=$(date +%s)
  taskset -c $CORES $Z cert $C "$@" --threads $TH --log $OUT/$name.log > $OUT/$name.out 2>&1 || true
  local t1=$(date +%s)
  echo "$name: $((t1 - t0)) s wall this invocation"
}
sweep cand_d4 --d4
sweep cand_full --full
{
  echo "# zmx2 run manifest (written by search/zmx2_run.sh; see search/ZMX2.md sec 9)"
  echo "date:            $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host:            $(uname -n) $(uname -m), cores $CORES ($TH threads)"
  echo "git HEAD:        $(git rev-parse HEAD) ($(git status --porcelain -- verify2/src/bin/zmx2.rs | wc -l) local change(s) to zmx2.rs)"
  echo "checker source:  verify2/src/bin/zmx2.rs sha256 $(sha256sum verify2/src/bin/zmx2.rs | cut -d' ' -f1)"
  echo "checker binary:  $Z sha256 $(sha256sum $Z | cut -d' ' -f1)"
  echo "rustc:           $(rustc --version)"
  echo "input:           $C sha256 $(sha256sum $C | cut -d' ' -f1)"
  echo "input summary:   $($Z info $C | tr '\n' ' ')"
  echo "D4 check:        $($Z d4 $C)"
  for name in cand_d4 cand_full; do
    echo
    echo "[$name]"
    echo "command:         taskset -c $CORES $Z cert $C $( [ $name = cand_d4 ] && echo --d4 || echo --full ) --threads $TH --log $OUT/$name.log"
    echo "log header:      $(head -1 $OUT/$name.log)"
    echo "log sha256:      $(sha256sum $OUT/$name.log | cut -d' ' -f1)  ($(grep -c '^ROOT ' $OUT/$name.log) ROOT lines, $(grep -c '^UNCERT ' $OUT/$name.log) UNCERT lines)"
    grep -E '^(roots:|done in|VERIFIED|NOT VERIFIED|INCOMPLETE|REGION)' $OUT/$name.out | sed 's/^/                 /'
  done
} > $MAN
cat $MAN
