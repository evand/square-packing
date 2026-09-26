#!/usr/bin/env bash
# zmcheck performance benchmark: builds several commits of verify2 and times them on fixed slices
# of the s(32) candidate, checking that every build gives the identical ROOT census.
#
#   verify2/bench/bench.sh [-n REPS] [-q] [-o OUTDIR] [COMMIT ...]
#
#   COMMIT   commits to compare (default: every commit touching verify2 from the first one with
#            ZM_UBINS, 4d49633, to HEAD; the first one listed is the baseline for the census diff)
#   -n REPS  repetitions of every (build, benchmark, threads) run (default 3)
#   -q       quick: skip the long germ root `germ0` (~6 min per run at baseline)
#   -o DIR   output directory (default /tmp/zmbench.<date>)
#
# Benchmarks (all on runs/s32-close_candidate.txt; ZM_ROOTLOG=1 for the census):
#   tile   binding edge tile cell x in [3.4,3.6], y in [5.4,5.6], default settings, 32 roots;
#          at 1 and 4 threads.  Census: 678 boxes, depth 11, ADM 99 / DISJ 184 / EMPTY 72.
#   germB  one root near the (3.5,2.5) germ, germ settings (--branch-cap 640 --node-cap 4000000
#          --depth 22), x3.4 y2.3 u[2/8,3/8]: ranking-bound (weight_of); 1 thread.
#   germ0  the same cell's u[0/8,1/8] root: forced two-chain path, cover bookkeeping; 1 thread.
#
# Per run it records user+sys CPU time, wall time, max RSS, and (if `perf` is on PATH or in
# $PERF) user-mode cycles, instructions, IPC, cache references/misses, L1d load misses.  Needs
# GNU time (`command -v time` is often a shell keyword; set $GNUTIME if it is not found).
#
# Expected runtime on an idle 16-core Zen 5 with the default 6 commits and -n 3: roughly
# 1.5-2 h (dominated by the baseline builds' germ0 runs); with -q, about 45 min.  Each run uses
# at most 4 threads, one run at a time.
#
# Output: DIR/results.csv (one row per run), DIR/summary.txt (per-build medians and speedups
# against the baseline), DIR/census.diff (empty = all builds agree), DIR/<commit>/... raw files.
set -euo pipefail

REPS=3; QUICK=0; OUT=/tmp/zmbench.$(date +%Y%m%d-%H%M%S)
while getopts "n:qo:" o; do
  case $o in n) REPS=$OPTARG ;; q) QUICK=1 ;; o) OUT=$OPTARG ;; *) sed -n 2,30p "$0"; exit 2 ;; esac
done
shift $((OPTIND - 1))

HERE=$(cd "$(dirname "$0")" && pwd)
V2=$(dirname "$HERE")                       # .../s12/verify2
S12=$(dirname "$V2")
REPO=$(git -C "$V2" rev-parse --show-toplevel)
CERT=${CERT:-$S12/runs/s32-close_candidate.txt}    # runs/ is gitignored: set CERT in a worktree
[ -f "$CERT" ] || { echo "missing $CERT (set CERT=path/to/s32-close_candidate.txt)"; exit 2; }
GNUTIME=${GNUTIME:-$(for t in /usr/bin/time /run/current-system/sw/bin/time; do [ -x "$t" ] && echo "$t" && break; done)}
[ -n "$GNUTIME" ] || { echo "GNU time not found; set GNUTIME"; exit 2; }
PERF=${PERF:-$(command -v perf || true)}

if [ $# -eq 0 ]; then
  set -- $(git -C "$REPO" rev-list --reverse 4d49633^..HEAD -- s12/verify2/src s12/verify2/Cargo.toml)
fi
COMMITS=()
for c in "$@"; do COMMITS+=("$(git -C "$REPO" rev-parse --short "$c")"); done
BASE=${COMMITS[0]}
mkdir -p "$OUT"
echo "commits: ${COMMITS[*]} (baseline $BASE); reps $REPS; quick $QUICK; out $OUT"
echo "load at start: $(cut -d' ' -f1-3 /proc/loadavg)" | tee "$OUT/load.txt"

# ---- build: each commit from `git archive` into its own source tree and CARGO_TARGET_DIR
for c in "${COMMITS[@]}"; do
  d=$OUT/$c
  if [ ! -x "$d/zmcheck" ]; then
    rm -rf "$d/src"; mkdir -p "$d/src"
    git -C "$REPO" archive "$c" s12/verify2 | tar -x -C "$d/src"
    ( cd "$d/src/s12/verify2" && CARGO_TARGET_DIR="$d/target" cargo build --release ) > "$d/build.log" 2>&1 \
      || { cat "$d/build.log"; echo "build of $c failed"; exit 1; }
    cp "$d/target/release/zmcheck" "$d/zmcheck"
  fi
done

# ---- runs
declare -A ARGS ENVV THREADS
ARGS[tile]="--xlo 3.4 --xhi 3.5 --ylo 5.4 --yhi 5.5";                                ENVV[tile]="";          THREADS[tile]="1 4"
GERM="--branch-cap 640 --node-cap 4000000 --depth 22 --xlo 3.4 --xhi 3.4 --ylo 2.3 --yhi 2.3"
ARGS[germB]=$GERM;                                                                    ENVV[germB]="ZM_UBINS=2"; THREADS[germB]="1"
ARGS[germ0]=$GERM;                                                                    ENVV[germ0]="ZM_UBINS=0"; THREADS[germ0]="1"
BENCHES="tile germB"; [ "$QUICK" = 1 ] || BENCHES="$BENCHES germ0"
EV=cycles:u,instructions:u,cache-references:u,cache-misses:u,L1-dcache-load-misses:u

echo "commit,bench,threads,rep,user_s,sys_s,wall_s,maxrss_kb,cycles,instructions,ipc,cache_refs,cache_misses,l1d_misses,load1" > "$OUT/results.csv"
for r in $(seq 1 "$REPS"); do
  for c in "${COMMITS[@]}"; do
    for b in $BENCHES; do
      for t in ${THREADS[$b]}; do
        o=$OUT/$c/$b.t$t.r$r
        pre=()
        [ -n "$PERF" ] && pre=("$PERF" stat -x, -e "$EV" -o "$o.perf" --)
        env ZM_ROOTLOG=1 ${ENVV[$b]} "$GNUTIME" -f "%U,%S,%e,%M" -o "$o.time" \
          "${pre[@]}" "$OUT/$c/zmcheck" cert "$CERT" --threads "$t" ${ARGS[$b]} > "$o.out" 2> "$o.err"
        # census: stdout without timings, plus the ROOT lines without their seconds, sorted
        { grep -v -E '^done in|progress:|threads$' "$o.out" | sed -E 's/; [0-9]+ threads//'
          grep '^ROOT' "$o.err" | sed -E 's/ [0-9.]+s$//' | sort; } > "$o.census"
        pv() { [ -f "$o.perf" ] && awk -F, -v e="$1" '$3 ~ "^"e {print $1}' "$o.perf" | head -1 || true; }
        cy=$(pv cycles); in=$(pv instructions)
        ipc=$( [ -n "$cy" ] && [ -n "$in" ] && awk -v a="$in" -v b="$cy" 'BEGIN{printf "%.3f", a/b}' || true)
        echo "$c,$b,$t,$r,$(cat "$o.time"),$cy,$in,$ipc,$(pv cache-references),$(pv cache-misses),$(pv L1-dcache-load-misses),$(cut -d' ' -f1 /proc/loadavg)" \
          | tee -a "$OUT/results.csv"
      done
    done
  done
done

# ---- census: every build and rep against the baseline's first rep
: > "$OUT/census.diff"
for c in "${COMMITS[@]}"; do for b in $BENCHES; do for t in ${THREADS[$b]}; do for r in $(seq 1 "$REPS"); do
  diff "$OUT/$BASE/$b.t$t.r1.census" "$OUT/$c/$b.t$t.r$r.census" > /dev/null \
    || { echo "== $c $b t$t r$r differs from $BASE"; diff "$OUT/$BASE/$b.t$t.r1.census" "$OUT/$c/$b.t$t.r$r.census"; } >> "$OUT/census.diff"
done; done; done; done
[ -s "$OUT/census.diff" ] && echo "CENSUS DIFFERS: see $OUT/census.diff" || echo "census: all builds identical"

# ---- summary: medians of CPU time (user+sys), wall and instructions; speedup vs baseline
python3 - "$OUT/results.csv" "$BASE" > "$OUT/summary.txt" <<'EOF'
import csv, sys, statistics as st
rows = list(csv.DictReader(open(sys.argv[1]))); base = sys.argv[2]
key = lambda r: (r["bench"], r["threads"])
med = {}
for r in rows:
    k = (r["commit"],) + key(r)
    cpu = float(r["user_s"]) + float(r["sys_s"])
    med.setdefault(k, []).append((cpu, float(r["wall_s"]), float(r["instructions"] or "nan"),
                                  float(r["ipc"] or "nan")))
print(f"{'commit':10} {'bench':6} {'thr':>3} {'cpu_s':>8} {'wall_s':>8} {'Ginstr':>8} {'IPC':>5} {'cpu x':>6} {'instr x':>7}")
for k in sorted(med, key=lambda k: (k[1], k[2], [r['commit'] for r in rows].index(k[0]))):
    v = med[k]; c = st.median(x[0] for x in v); w = st.median(x[1] for x in v)
    i = st.median(x[2] for x in v); p = st.median(x[3] for x in v)
    b = med.get((base,) + k[1:]); bc = st.median(x[0] for x in b); bi = st.median(x[2] for x in b)
    print(f"{k[0]:10} {k[1]:6} {k[2]:>3} {c:8.1f} {w:8.1f} {i/1e9:8.1f} {p:5.2f} {bc/c:6.2f} {bi/i:7.2f}")
EOF
cat "$OUT/summary.txt"
