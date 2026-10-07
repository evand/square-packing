#!/bin/sh
# s(32) = 6: re-check the bundle (certificates/s32/README.md).  Run from anywhere.
#
#   certificates/s32/verify.sh          fast (seconds): hashes; the shipped zeromargin.py run records
#                                       re-summarised from scratch (all 7,200 roots present and
#                                       certified, hashes match, D4 invariance of the weighted cover
#                                       re-checked exactly); the Lean data file is this certificate.
#   certificates/s32/verify.sh --full   also re-runs the whole D4 sweep with zeromargin.py for both
#                                       covers (~2.8 + 2.4 CPU-h; about 12 min on 28 processes).
set -e
cd "$(dirname "$0")/../.."          # repo root
B=certificates/s32
FULL=0
case "${1:-}" in --full) FULL=1 ;; "") ;; *) echo "usage: $0 [--full]"; exit 2 ;; esac
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

echo "--- hashes"
( cd $B && sha256sum -c SHA256SUMS )
cmp $B/zeromargin_d4/checker/zeromargin.py $B/zeromargin_d4_shift_v1/checker/zeromargin.py

echo "--- Lean data (lean/Sqpack/S32Data.lean) is s32_closed_cover_6.txt"
python3 lean/scripts/gen_s32_data.py $B/s32_closed_cover_6.txt --check

summ() {   # $1 = records dir; summarised in a scratch copy (summary rewrites SUMMARY.txt)
  mkdir -p "$TMP/$1"; cp -r "$B/$1/." "$TMP/$1/"
  out=$(python3 search/zm_d4_sweep.py summary --out "$TMP/$1") || { echo "$out"; echo "NOT CLEAN: $1"; exit 1; }
  echo "$out" | grep -E '^(D4 invariance|region:|totals:|D4 RECHECK)'
}
for d in zeromargin_d4 zeromargin_d4_shift_v1; do echo "--- shipped zeromargin.py run: $d"; summ $d; done

echo "--- fast re-execution: two interior tile-germ roots of s32_closed_cover_6.txt, census must equal the shipped records"
for r in '3/2,8/5,7/5,3/2,1/16,1/8' '7/5,3/2,7/5,3/2,0,1/16'; do
  python3 search/zm_d4_sweep.py run --cert $B/s32_closed_cover_6.txt --out "$TMP/sample" --depth 24 --nproc 2 --only "$r" > /dev/null 2>&1
done
python3 - "$TMP/sample/roots.jsonl" $B/zeromargin_d4/roots.jsonl <<'EOF'
import json, sys
new = [json.loads(l) for l in open(sys.argv[1])]
old = {}
for l in open(sys.argv[2]):
    r = json.loads(l); old.setdefault(r["root"], r)
assert len(new) == 2, new
for r in new:
    assert r["stats"]["UNCERT"] == 0 and r["stats"] == old[r["root"]]["stats"], (r, old[r["root"]])
    print(f"    {r['root']}: {r['stats']['boxes']} boxes, certified, census identical")
EOF

echo "--- shipped zmx2 --full --pair-points --sym-atoms run (zmx2_full_sym/): all 28,800 roots, none uncertified, no D4 fold"
xz -dc certificates/s32/zmx2_full_sym/roots.log.xz | python3 -c '
import sys
r = [l.split() for l in sys.stdin if l.startswith("ROOT ")]   # ROOT id pass P root R boxes b cert c empty e uncert u ...
keys = {(f[3], f[5]) for f in r}
assert len(r) == 28800 and len(keys) == 28800, (len(r), len(keys))
assert {f[3] for f in r} == {"0", "1"} and all(f[13] == "0" for f in r), "uncertified or unexpected pass"
print("    28,800 distinct roots (both passes), 0 uncertified")'

if [ "$FULL" = 1 ]; then
  for c in s32_closed_cover_6 s32_shift_v1; do
    echo "--- full re-run: zeromargin.py over the D4 region of $c.txt (depth 24, then 30 on any leftover root)"
    O="$TMP/rerun_$c"
    python3 search/zm_d4_sweep.py run --cert $B/$c.txt --out "$O" --depth 24 --nproc "$(nproc)" > "$O.log" 2>&1
    python3 search/zm_d4_sweep.py run --cert $B/$c.txt --out "$O" --depth 24 --deepen 30 --nproc "$(nproc)" >> "$O.log" 2>&1
    out=$(python3 search/zm_d4_sweep.py summary --out "$O") || { echo "$out"; tail -20 "$O.log"; echo "NOT CLEAN: $c"; exit 1; }
    echo "$out" | grep -E '^(checker|region:|totals:|D4 RECHECK)'
  done
fi
echo "s(32) bundle: OK"
