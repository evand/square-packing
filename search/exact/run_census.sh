#!/usr/bin/env bash
# LOCAL (needs ../packer): census minima through pipeline.py (exact solve; if not a local minimum, kick + slp2 + solve) -> descended/.
set -u
cd "$(dirname "$0")"
mkdir -p descended
for f in inputs/cen3_*.txt; do
  python3 pipeline.py "$f" > descended/$(basename "$f" .txt).pipeline.txt 2>&1
done
for c in descended/*.cert; do echo "== $c"; python3 verify_cert.py "$c" | tail -1; done
