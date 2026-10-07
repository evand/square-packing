#!/bin/sh
# wait_polish.sh COUNT PATTERN... : block until COUNT of the given polish logs contain "final:"
cd "$(dirname "$0")/../.." || exit 1
want=$1; shift
while :; do
  n=0
  for f in "$@"; do
    if grep -q "final:" "$f" 2>/dev/null; then n=$((n+1)); fi
  done
  if [ "$n" -ge "$want" ]; then echo "done: $n finished"; exit 0; fi
  sleep 20
done
