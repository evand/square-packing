#!/bin/sh
# finalize.sh CERT NAME : install a uniform certificate as certificates/NAME.txt (+ .json) and run
# every check the repo requires: verify at N=2000 and N=8000, xcheck (every bin, N=2000),
# JSON export / roundtrip, and scale_to_critical.
cd "$(dirname "$0")/../.." || exit 1
CERT=$1; NAME=$2
OUT=certificates/$NAME.txt
cp "$CERT" "$OUT"
echo "== $OUT"
head -n 4 "$OUT" | tr '\n' ' '; echo
for N in 2000 8000; do
  printf "verify N=%s: " "$N"
  ( time verify/target/release/verify "$OUT" 12 "$N" "$(nproc)" 0 ) 2>&1 | grep -E "^(VERIFIED|NOT VERIFIED|ERROR|real)" | tr '\n' ' '; echo
done
printf "xcheck N=2000 --all: "
( time python3 xcheck.py "$OUT" 2000 --all --n 12 ) 2>&1 | grep -E "VERIF|minimum|real" | tr '\n' ' '; echo
python3 search/export_points.py export "$OUT" -o "certificates/$NAME.json" --n 12 && \
python3 search/export_points.py --roundtrip "$OUT" "certificates/$NAME.json" | tail -n 1
D=$(sed -n 2p "$OUT")
LO=$((D - 40)); if [ "$LO" -lt 1 ]; then LO=1; fi
python3 search/scale_to_critical.py "$OUT" --n 12 --lo $LO | grep -E "^(critical|best|==>)"
