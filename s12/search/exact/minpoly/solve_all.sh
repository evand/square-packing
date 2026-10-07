#!/usr/bin/env bash
# exactsolve on every register input (batch/inputs), writing NAME.contacts.json for minpoly.py; run inside ./env.sh
cd "$(dirname "$0")/.."
ls batch/inputs/n-*.txt | sort -t- -k2 -n | xargs -P "${PROCS:-6}" -I{} sh -c 'b=$(basename {} .txt); [ -f minpoly/solve/$b.contacts.json ] || timeout 3600 python3 exactsolve.py {} --out minpoly/solve > minpoly/solve/$b.log 2>&1'
