#!/usr/bin/env bash
# zmx2_sym_run.sh -- the zmx2 --sym-atoms runs (search/ZMX2.md sec 4.9, sec 12), resumable, with a manifest.
# Run from the repo root (detached is fine):
#   setsid nohup bash search/zmx2_sym_run.sh > runs/zmx2_sym/run.out 2>&1 < /dev/null &
# Sweeps (each with a per-root log runs/zmx2_sym/<name>.log; a rerun skips finished roots):
#   s32_full   s(32) cover, --full --pair-points --sym-atoms   (28,800 roots: the goal, no symmetry used)
#   s32_d4     s(32) cover, --d4   --pair-points --sym-atoms   (3,600 roots)
#   s21_d4/s21_full, s60_d4/s60_full   the s(21), s(60) bundle covers, --sym-atoms   (regression; there the
#              mirrored assignment equals the default one, so the census must equal the shipped roots.log)
#   s21_pp_d4, s60_pp_d4   --d4 --pair-points --sym-atoms, where the two assignments differ (--full with
#              --pair-points costs ~8x more; not run)
# Writes search/zmx2_sym_manifest.txt (checker sha, settings, census, verdict per sweep).
# ONLY="s32_full s21_d4" restricts the sweeps (the manifest then lists those only).
set -eu
cd "$(dirname "$0")/.."
TH=${TH:-6}
CORES=${CORES:-0-5}
OUT=runs/zmx2_sym; mkdir -p $OUT
MAN=${MAN:-search/zmx2_sym_manifest.txt}
(cd verify2 && cargo build --release --bin zmx2 -q) || { echo "cargo build of zmx2 failed"; exit 1; }
Z=verify2/target/release/zmx2
S32=certificates/s32/s32_closed_cover_6.txt
S21=certificates/s21/s21_mixed_cover_5.txt
S60=certificates/s60/s60_mixed_cover_8.txt
ALL="s32_full s32_d4 s21_d4 s21_full s60_d4 s60_full s21_pp_d4 s60_pp_d4"
ONLY=${ONLY:-$ALL}
declare -A CMD
CMD[s32_full]="$S32 --full --pair-points --sym-atoms"
CMD[s32_d4]="$S32 --d4 --pair-points --sym-atoms"
CMD[s21_d4]="$S21 --d4 --sym-atoms"
CMD[s21_full]="$S21 --full --sym-atoms"
CMD[s60_d4]="$S60 --d4 --sym-atoms"
CMD[s60_full]="$S60 --full --sym-atoms"
CMD[s21_pp_d4]="$S21 --d4 --pair-points --sym-atoms"
CMD[s60_pp_d4]="$S60 --d4 --pair-points --sym-atoms"
for name in $ONLY; do
  t0=$(date +%s)
  taskset -c $CORES $Z cert ${CMD[$name]} --threads $TH --log $OUT/$name.log > $OUT/$name.out 2>&1 || true
  t1=$(date +%s)
  echo "$name: $((t1 - t0)) s wall this invocation"
done
{
  echo "# zmx2 --sym-atoms run manifest (written by search/zmx2_sym_run.sh; see search/ZMX2.md sec 12)"
  echo "date:            $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host:            $(uname -m), taskset -c $CORES ($TH threads)"
  echo "git HEAD:        $(git rev-parse HEAD) ($(git status --porcelain -- verify2/src/bin/zmx2.rs | wc -l) local change(s) to zmx2.rs)"
  echo "checker source:  verify2/src/bin/zmx2.rs sha256 $(sha256sum verify2/src/bin/zmx2.rs | cut -d' ' -f1)"
  echo "checker binary:  $Z sha256 $(sha256sum $Z | cut -d' ' -f1)"
  echo "rustc:           $(rustc --version)"
  for name in $ONLY; do
    set -- ${CMD[$name]}
    echo
    echo "[$name]"
    echo "input:           $1 sha256 $(sha256sum $1 | cut -d' ' -f1)"
    echo "command:         taskset -c $CORES $Z cert ${CMD[$name]} --threads $TH --log $OUT/$name.log"
    echo "log header:      $(head -1 $OUT/$name.log)"
    echo "log sha256:      $(sha256sum $OUT/$name.log | cut -d' ' -f1)  ($(grep -c '^ROOT ' $OUT/$name.log) ROOT lines, $(grep -c '^UNCERT ' $OUT/$name.log || true) UNCERT lines)"
    grep -E '^(roots:|done in|VERIFIED|NOT VERIFIED|INCOMPLETE|REGION)' $OUT/$name.out | sed 's/^/                 /'
    case $name in s21_d4|s21_full|s60_d4|s60_full)
      echo "vs shipped:      $(python3 search/zmx2_census_cmp.py $OUT/$name.log certificates/${name%_*}/zmx2_${name#*_}/roots.log)" ;;
    esac
  done
} > $MAN
cat $MAN
