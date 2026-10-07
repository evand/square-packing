#!/usr/bin/env bash
# Replay verify_exact.py (stdlib only) on every committed exact form: minpoly/data/n-N.minpoly.json.gz.
# Writes minpoly/data/verify.txt (one line per n) and prints the count of VALID.  PROCS parallel (default 4).
cd "$(dirname "$0")/.."
ls minpoly/data/n-*.minpoly.json.gz | sort -t- -k2 -n | xargs -P "${PROCS:-4}" -I{} sh -c \
  'timeout 3600 python3 verify_exact.py {} 2>&1 | tail -1' > minpoly/data/verify.txt
sort -t- -k2 -n -o minpoly/data/verify.txt minpoly/data/verify.txt
echo "VALID: $(grep -c ': VALID\.' minpoly/data/verify.txt) of $(ls minpoly/data/n-*.minpoly.json.gz | wc -l)"
