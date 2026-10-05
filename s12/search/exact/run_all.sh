#!/usr/bin/env bash
# Re-run the test set: exact solve + certificate for every input; logs/JSON/certificates in results/, then both verifiers.
set -u
cd "$(dirname "$0")"
mkdir -p results
for f in site5 site10; do python3 exactsolve.py inputs/$f.txt --algdeg 4 > results/$f.log 2>&1; done
python3 exactsolve.py inputs/site11.txt --dps 300 --algdeg 16 > results/site11.log 2>&1
for f in site17 site71 rec110 cand110; do python3 exactsolve.py inputs/$f.txt > results/$f.log 2>&1; done
for c in results/*.cert; do
  echo "== $c: $(python3 verify_cert.py "$c" | tail -1); verify_cert2: $(python3 verify_cert2.py "$c" > /dev/null && echo VALID || echo INVALID)"
done
