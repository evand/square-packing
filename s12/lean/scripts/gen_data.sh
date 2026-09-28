#!/usr/bin/env bash
# Regenerate the Lean certificate data that is too big to commit (gitignored).
# Run from anywhere; then build the opt-in module, e.g. `lake build Sqpack.S11Lower`.
# The generator is deterministic and untrusted: the kernel checks whatever it writes.
#   usage: lean/scripts/gen_data.sh [S11 S13 ...]     (default: all)
set -euo pipefail
cd "$(dirname "$0")/../.."          # the s12/ directory
G="python3 lean/scripts/gen_boxtree.py"

gen() {
  case "$1" in
    S11) $G certificates/s11_lower_3.8143.txt --n 11 --look 0 --parts 24 --name S11 --outdir lean/Sqpack/S11 ;;  # ~4 min
    S13) python3 lean/scripts/gen_zmtree.py certificates/rung2/s13_closed_cover_4.txt --n 13 --name S13 --outdir lean/Sqpack/S13 --nproc 4 ;;  # ~2.5 min on 4 cores
    *) echo "unknown data set: $1" >&2; exit 2 ;;
  esac
}

for d in ${@:-S11 S13}; do gen "$d"; done
