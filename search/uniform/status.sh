#!/bin/sh
# one-line status of every polish and ILP run under runs/
cd "$(dirname "$0")/../.." || exit 1
echo "== polish runs (last improvement line)"
for f in runs/plog_*.txt; do
  last=$(grep "D_crit" "$f" | tail -n 1 | cut -c1-100)
  echo "$(basename "$f" .txt): $last"
done
echo "== ILP runs"
grep -h "RESULT" runs/log_*.txt 2>/dev/null | sort
echo "== still running:"
pgrep -fc "ilp_uniform.py"; pgrep -fc "polish.py"
uptime
