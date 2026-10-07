#!/usr/bin/env bash
# s21_cert_runs.sh -- the certified runs shipped in certificates/s21/ (s(21) = 5), with manifests.
# Run from anywhere (detached is fine):
#   setsid nohup bash search/s21_cert_runs.sh [zmx2|zm_mixed|all] > runs/s21cert.out 2>&1 < /dev/null &
# Env: CORES (taskset list, default 0-9), NPROC (default 10).
#
#   zmx2      verify2 zmx2 on the cover: --d4 (2,500 roots) and --full (no symmetry, 20,000 roots); ~2 CPU-min.
#             Writes certificates/s21/zmx2_{d4,full}/{roots.log,run.out,manifest.txt}.
#   zm_mixed  search/zm_mixed.py --d4 --cert-mode on the cover (40,000 roots, ~20 CPU-h, ~2 h on 10 processes).
#             Resumable (roots.jsonl, header-checked).  Writes certificates/s21/zm_mixed_d4/{roots.jsonl,run.log,
#             manifest.json} and checker/ (the exact checker files that ran).
set -eu
cd "$(dirname "$0")/.."            # repo root
WHAT=${1:-all}
CORES=${CORES:-0-9}
NPROC=${NPROC:-10}
B=certificates/s21
C=$B/s21_mixed_cover_5.txt

if [ "$WHAT" = zmx2 ] || [ "$WHAT" = all ]; then
  (cd verify2 && cargo build --release --bin zmx2 -q) || { echo "cargo build of zmx2 failed"; exit 1; }
  Z=verify2/target/release/zmx2
  for mode in d4 full; do
    O=$B/zmx2_$mode; mkdir -p $O
    rm -f $O/roots.log            # a fresh run, not a resumed one
    t0=$(date +%s)
    taskset -c $CORES $Z cert $C --$mode --threads $NPROC --log $O/roots.log > $O/run.out 2>&1 || true
    t1=$(date +%s)
    {
      echo "# zmx2 run manifest (search/s21_cert_runs.sh; certificates/s21/README.md)"
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
  taskset -c $CORES python3 search/zm_mixed.py cert $C --d4 --cert-mode --disj --chain-from 0 --depth 24 \
      --pitch 1/20 --ubins 16 --nproc $NPROC --progress 1000 \
      --resume $O/roots.jsonl --manifest $O/manifest.json >> $O/run.log 2>&1 || true
  for f in zm_mixed.py mixed_cover.py zeromargin.py; do
    cmp -s search/$f $O/checker/$f || { echo "search/$f changed during the run"; exit 1; }
  done
  tail -5 $O/run.log
fi
