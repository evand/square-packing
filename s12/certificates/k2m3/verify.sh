#!/bin/sh
# s(k² − 3) = k for all k ≥ 6: re-check the bundle (certificates/k2m3/README.md).  Run from anywhere.
#
#   certificates/k2m3/verify.sh          a few seconds: hashes; the box cover's exact total and D4 invariance (own
#                                        parser); the family rebuilds the box cover exactly (qx2_family_check.py);
#                                        Lemma Z (the θ = 0 face) re-run with the shipped checker and compared with the
#                                        shipped output; the shipped run V3 record re-checked from scratch
#                                        (search/qx2_records.py: header, settings, roots, every leaf, coverage,
#                                        census); the Lean data file BentzData.lean regenerated from the shipped
#                                        files and compared byte for byte.
#   certificates/k2m3/verify.sh --full   also re-runs qx2_zm.py with the run V3 settings over the whole D4 region
#                                        (9,800 roots, ~81,000 CPU-s; about 3 h on 8 processes), re-checks the new
#                                        record the same way and compares its census with the shipped one root for
#                                        root (every leaf count and every leaf box).
# NPROC (default: nproc) sets the number of processes for --full.
set -e
cd "$(dirname "$0")/../.."          # s12/
export PYTHONDONTWRITEBYTECODE=1
B=certificates/k2m3
C=$B/L4_k02_box7.txt
K=$B/qx2_zm/checker
R=$B/qx2_zm/runV3_6294052a_leaves
FULL=0
case "${1:-}" in --full) FULL=1 ;; "") ;; *) echo "usage: $0 [--full]"; exit 2 ;; esac
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
NP=${NPROC:-$(nproc)}
run() {   # run a check, indent its output, stop on failure (a pipe into sed would hide the exit status)
  "$@" > "$TMP/o" 2>&1 && st=0 || st=$?
  sed 's/^/    /' "$TMP/o"
  [ "$st" = 0 ] || { echo "FAILED ($st): $*"; exit 1; }
}

echo "--- hashes"
( cd $B && sha256sum -c SHA256SUMS )

echo "--- the box cover: well-formed, exact total < 46, D4-invariant measure (search/qx2_records.py, own parser)"
run python3 search/qx2_records.py cover $C 46

echo "--- the family file rebuilds the box cover exactly; sigma = 0, D > 3/4 (search/qx2_family_check.py)"
run python3 search/qx2_family_check.py $B/L4_k02_family.txt $C

echo "--- Lemma Z: the theta = 0 face, exact (shipped checker, qx2_zm.py axis)"
python3 $K/qx2_zm.py axis $C > "$TMP/lemmaZ.out"
sed 's/^/    /' "$TMP/lemmaZ.out"
cmp "$TMP/lemmaZ.out" $B/qx2_zm/lemmaZ.out
echo "    identical to the shipped qx2_zm/lemmaZ.out"

echo "--- shipped run V3 (qx2_zm/): header shas = checker/ + cover, settings, 9,800 D4 roots, every leaf, coverage"
run python3 search/qx2_records.py record $R.jsonl.gz $R.out $K $C

echo "--- Lean data (lean/Sqpack/BentzData.lean) regenerated from the shipped cover and family"
mkdir -p "$TMP/lean/qx2_data"; cp $C $B/L4_k02_family.txt "$TMP/lean/qx2_data/"
run python3 lean/scripts/gen_bentz_data.py --search "$TMP/lean" --out "$TMP/lean/BentzData.lean"
cmp "$TMP/lean/BentzData.lean" lean/Sqpack/BentzData.lean
echo "    identical to lean/Sqpack/BentzData.lean"

if [ "$FULL" = 1 ]; then
  echo "--- full re-run: qx2_zm.py (shipped checker/) with the run V3 settings over the D4 region (~81,000 CPU-s)"
  for f in qx2_zm.py zm_mixed.py mixed_cover.py zeromargin.py; do
    cmp -s search/$f $K/$f || echo "    note: search/$f differs from the shipped checker/$f"
  done
  O="$TMP/full"; mkdir -p "$O"
  python3 $K/qx2_zm.py $C --depth 18 --nproc "$NP" --exact-umax 1/2 --exact-from 3 --progress 1000 \
      --resume "$O/V3.jsonl" --dump-leaves > "$O/run.out" 2>&1 || true
  tail -3 "$O/run.out" | sed 's/^/    /'
  gzip -k "$O/V3.jsonl"
  run python3 search/qx2_records.py record "$O/V3.jsonl.gz" "$O/run.out" $K $C
  python3 - "$O/V3.jsonl" $R.jsonl.gz <<'EOF'
import gzip, json, sys
def census(p):
    op = gzip.open if p.endswith('.gz') else open
    L = [json.loads(l) for l in op(p, 'rt') if l.strip()][1:]
    return {tuple(r["root"]): ({k: v for k, v in r["st"].items() if k != "cpu"}, sorted(map(str, r["leaves"])))
            for r in L}
a, b = census(sys.argv[1]), census(sys.argv[2])
same = sum(1 for k in b if a.get(k) == b[k])
print(f"    census and leaves identical to the shipped record for {same} / {len(b)} roots")
assert same == len(b) == len(a), "the fresh run differs from the shipped record"
EOF
fi
echo "s(k^2 - 3) = k bundle: OK"
