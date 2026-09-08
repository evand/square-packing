#!/bin/sh
# Restart driver for search/t4screen.py.
#
# The screening loop (rows + clique separation + column pricing) does not converge in one sitting
# at t = 4 -- the same degeneracy search/CLIQUE_CEILING.md records -- and no single process here is
# allowed to run longer than ~15 minutes.  So run it in chunks: each chunk writes
# `runs/t4_<TAG>_poses.txt` (its pose set) and `runs/t4_<TAG>_cliques.txt` (its cliques as exact
# integer anchors), and the next chunk warm-starts from both.  Cliques come back as CANDIDATES,
# so a clique the LP has stopped using is only re-added if it is violated again.
#
#   sh search/t4screen_loop.sh T TAG NCHUNKS SECS  [extra t4screen args...]
#
# Deterministic given (T, TAG, NCHUNKS, SECS, args) only up to the wall-clock cut inside a chunk;
# for a bit-reproducible run use one chunk with --stages/--rowloops as the stopping rule and
# --time large.
set -e
T=$1; TAG=$2; N=$3; SECS=$4
shift 4
for i in $(seq 1 "$N"); do
  EXTRA=""
  [ -f "runs/t4_${TAG}_poses.txt" ]   && EXTRA="$EXTRA --load-poses runs/t4_${TAG}_poses.txt"
  [ -f "runs/t4_${TAG}_cliques.txt" ] && EXTRA="$EXTRA --cq-load runs/t4_${TAG}_cliques.txt"
  echo "=== $TAG chunk $i/$N (t=$T, ${SECS}s) ==="
  # shellcheck disable=SC2086
  timeout $((SECS + 240)) python3 search/t4screen.py "$T" "$TAG" --time "$SECS" $EXTRA "$@" || true
done
