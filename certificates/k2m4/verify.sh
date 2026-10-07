#!/bin/sh
# s(k² − 4) = k for all k ≥ 8 (the family half of s(k² − 4) = k, k ≥ 5): re-check the bundle
# (certificates/k2m4/README.md).  Run from anywhere.
#
#   certificates/k2m4/verify.sh          ~20 s: hashes; the box cover's exact total < 77 (i.e. D > 1) and D4
#                                        invariance (own parser); the family rebuilds the box cover exactly
#                                        (qx2_family_check.py --k 9 --R 3 --w 3); Lemma Z (the θ = 0 face) re-run
#                                        with the shipped checker and compared with the shipped output; the shipped
#                                        record re-checked from scratch (search/qx2_records.py: header, settings,
#                                        16,200 roots, every leaf, coverage, census); the Lean data file
#                                        Bentz4Data.lean regenerated from the shipped files and compared byte for
#                                        byte.  NOT done in this mode: recomputing any theta > 0 leaf's mass bound.
#                                        The record check is structure and coverage only, so this is not a fresh
#                                        geometric check (--full is).
#   certificates/k2m4/verify.sh --full   also re-runs qx2_zm.py with the run's settings over the whole D4 region
#                                        (16,200 roots, ~815,000 CPU-s ≈ 226 CPU-h; about 16 h on 15 processes),
#                                        re-checks the new record the same way and compares it with the shipped one
#                                        root for root (every leaf count and every leaf box).
# NPROC (default: nproc) sets the number of processes for --full.
set -e
cd "$(dirname "$0")/../.."          # repo root
export PYTHONDONTWRITEBYTECODE=1
B=certificates/k2m4
C=$B/K4_k008_box9.txt
K=$B/qx2_zm/checker
R=$B/qx2_zm/run_k4x_k008_leaves
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

echo "--- the box cover: well-formed, exact total < 77 (= 81 - 4, i.e. D > 1), D4-invariant (search/qx2_records.py, own parser)"
run python3 search/qx2_records.py cover $C 77

echo "--- the family file rebuilds the box cover exactly; sigma = 0 (search/qx2_family_check.py, R = w = 3, k = 9)"
echo "    (its printed '4D > 3' and '< k^2 - 3' are the k^2 - 3 thresholds; D > 1 is the 'total < 77' above)"
run python3 search/qx2_family_check.py $B/K4_k008_family.txt $C --k 9 --R 3 --w 3

echo "--- Lemma Z: the theta = 0 face, exact (shipped checker, qx2_zm.py axis)"
python3 $K/qx2_zm.py axis $C > "$TMP/lemmaZ.out"
sed 's/^/    /' "$TMP/lemmaZ.out"
cmp "$TMP/lemmaZ.out" $B/qx2_zm/lemmaZ.out
echo "    identical to the shipped qx2_zm/lemmaZ.out"

echo "--- shipped run qx2_k4x_k008 (qx2_zm/): header shas = checker/ + cover, settings, 16,200 D4 roots, every leaf, coverage"
run python3 search/qx2_records.py record $R.jsonl.gz $R.out $K $C

echo "--- Lean data (lean/Sqpack/Bentz4Data.lean) regenerated from the shipped cover and family"
mkdir -p "$TMP/lean/qx2_data"; cp $C $B/K4_k008_family.txt "$TMP/lean/qx2_data/"
run python3 lean/scripts/gen_bentzfam_data.py --search "$TMP/lean" --out "$TMP/lean/Bentz4Data.lean"
cmp "$TMP/lean/Bentz4Data.lean" lean/Sqpack/Bentz4Data.lean
echo "    identical to lean/Sqpack/Bentz4Data.lean"

if [ "$FULL" = 1 ]; then
  echo "--- full re-run: qx2_zm.py (shipped checker/) with the run's settings over the D4 region (~815,000 CPU-s)"
  for f in qx2_zm.py zm_mixed.py mixed_cover.py zeromargin.py; do
    cmp -s search/$f $K/$f || echo "    note: search/$f differs from the shipped checker/$f"
  done
  O="$TMP/full"; mkdir -p "$O"
  python3 $K/qx2_zm.py $C --depth 18 --nproc "$NP" --exact-umax 1/2 --exact-from 3 --progress 1000 \
      --resume "$O/run.jsonl" --dump-leaves > "$O/run.out" 2>&1 || true
  tail -3 "$O/run.out" | sed 's/^/    /'
  gzip -k "$O/run.jsonl"
  run python3 search/qx2_records.py record "$O/run.jsonl.gz" "$O/run.out" $K $C
  python3 - "$O/run.jsonl" $R.jsonl.gz <<'EOF'
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
if [ "$FULL" = 1 ]; then echo "s(k^2 - 4) = k (k >= 8) bundle: OK (full: theta > 0 re-run with qx2_zm.py)"
else echo "s(k^2 - 4) = k (k >= 8) bundle: OK (fast: theta = 0 re-run; theta > 0 record checked for structure and coverage, leaves NOT recomputed; use --full)"; fi
