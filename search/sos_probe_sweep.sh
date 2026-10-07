#!/usr/bin/env bash
# epsilon sweep: sos_probe_sweep.sh "<orders>" "<eps list>" <sdp args...>
orders=$1; eps=$2; shift 2
for o in $orders; do for e in $eps; do
  T=$(python3 -c "print(3-$e)")
  "$(dirname "$0")/sos_probe_run.sh" --order $o --T $T "$@"
done; done
