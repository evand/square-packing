#!/bin/sh
# s(21) = 5: re-check the bundle (certificates/s21/README.md).  Run from anywhere.
#
#   certificates/s21/verify.sh          about a minute: hashes; the cover's exact total and D4 invariance (own
#                                       parser); the Lean data file is this cover; the shipped zm_mixed.py and
#                                       zmx2 run records re-summarised from scratch; and a FRESH zmx2 run over the
#                                       whole pose space with no symmetry assumed (20,000 roots, ~100 CPU-s).
#   certificates/s21/verify.sh --full   also re-runs zm_mixed.py --d4 --cert-mode over the whole D4 region
#                                       (40,000 roots, ~20 CPU-h; about 2 h on 10 processes) and compares its
#                                       census with the shipped records.
set -e
cd "$(dirname "$0")/../.."          # s12/
B=certificates/s21
C=$B/s21_mixed_cover_5.txt
FULL=0
case "${1:-}" in --full) FULL=1 ;; "") ;; *) echo "usage: $0 [--full]"; exit 2 ;; esac
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
run() {   # run a check, indent its output, stop on failure (a pipe into sed would hide the exit status under sh)
  "$@" > "$TMP/o" 2>&1 && st=0 || st=$?
  sed 's/^/    /' "$TMP/o"
  [ "$st" = 0 ] || { echo "FAILED ($st): $*"; exit 1; }
}
NP=${NPROC:-$(nproc)}

echo "--- hashes"
( cd $B && sha256sum -c SHA256SUMS )

echo "--- the cover: well-formed, exact total < 21, D4-invariant measure (search/s21_records.py, own parser)"
python3 search/s21_records.py cover $C

echo "--- Lean data (lean/Sqpack/S21Data.lean) is s21_mixed_cover_5.txt"
python3 lean/scripts/gen_s21_data.py $C --check

echo "--- shipped zm_mixed.py run (zm_mixed_d4/): header shas = checker/ + cover, settings, all 40,000 roots certified"
run python3 search/s21_records.py zm_mixed $B/zm_mixed_d4 $C
cmp $B/zm_mixed_d4/checker/zeromargin.py $B/../s32/zeromargin_d4/checker/zeromargin.py   # the pinned zeromargin.py

echo "--- shipped zmx2 runs (zmx2_d4/, zmx2_full/): every root of each region certified"
run python3 search/s21_records.py zmx2 $B/zmx2_d4/roots.log $C d4
run python3 search/s21_records.py zmx2 $B/zmx2_full/roots.log $C full

echo "--- fresh zmx2 run: the whole pose space, no symmetry assumed (--full: cover + mirror, 20,000 roots)"
( cd verify2 && cargo build --release --bin zmx2 2>&1 | tail -1 )
Z=verify2/target/release/zmx2
$Z d4 $C
out=$($Z cert $C --full --threads "$NP" --log "$TMP/full.log" 2>&1) || { echo "$out" | tail -20; echo "zmx2 failed"; exit 1; }
echo "$out" | grep -E '^(cover:|done in|VERIFIED|NOT VERIFIED|INCOMPLETE)'
echo "$out" | grep -q '^VERIFIED: ' || { echo "REJECTED by zmx2: $C"; exit 1; }
run python3 search/s21_records.py zmx2 "$TMP/full.log" $C full
python3 - "$TMP/full.log" $B/zmx2_full/roots.log <<'EOF'
import sys
def census(p):   # ROOT id pass P root R boxes b cert c empty e uncert u maxdepth m capped k ms t
    return {(f[3], f[5]): tuple(f[i] for i in (7, 9, 11, 13, 15, 17))
            for f in (l.split() for l in open(p) if l.startswith("ROOT "))}
a, b = census(sys.argv[1]), census(sys.argv[2])
assert a == b, "fresh zmx2 census differs from the shipped zmx2_full/roots.log"
print("    fresh census identical to the shipped zmx2_full/roots.log, root for root")
EOF

if [ "$FULL" = 1 ]; then
  echo "--- full re-run: zm_mixed.py --d4 --cert-mode over the D4 region (40,000 roots; ~20 CPU-h)"
  O="$TMP/zmm"; mkdir -p "$O/checker"; cp search/zm_mixed.py search/mixed_cover.py search/zeromargin.py "$O/checker/"
  for f in zm_mixed.py mixed_cover.py zeromargin.py; do
    cmp -s search/$f $B/zm_mixed_d4/checker/$f || echo "    note: search/$f differs from the shipped checker/$f"
  done
  python3 search/zm_mixed.py cert $C --d4 --cert-mode --disj --chain-from 0 --depth 24 --pitch 1/20 --ubins 16 \
      --nproc "$NP" --progress 5000 --resume "$O/roots.jsonl" --manifest "$O/manifest.json" > "$O/run.log" 2>&1 || true
  tail -4 "$O/run.log"
  run python3 search/s21_records.py zm_mixed "$O" $C
  python3 - "$O/roots.jsonl" $B/zm_mixed_d4/roots.jsonl <<'EOF'
import json, sys
def census(p):
    L = [json.loads(l) for l in open(p) if l.strip()][1:]
    return {tuple(r["root"]): {k: v for k, v in r["st"].items() if k != "cpu"} for r in L}
a, b = census(sys.argv[1]), census(sys.argv[2])
same = sum(1 for k in a if a[k] == b.get(k))
print(f"    census identical to the shipped records for {same} / {len(b)} roots")
EOF
fi
echo "s(21) bundle: OK"
