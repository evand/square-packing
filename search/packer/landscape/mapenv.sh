#!/usr/bin/env bash
# Run a command with python 3.13 + umap-learn, scikit-learn, matplotlib, scipy, numpy from nixpkgs (landscape map work).
#   ./mapenv.sh python3 map.py ...
set -euo pipefail
exec nix --extra-experimental-features 'nix-command flakes' shell --impure --expr \
  'with import (builtins.getFlake "nixpkgs") {}; [ (python313.withPackages (p: with p; [ umap-learn scikit-learn matplotlib scipy numpy ])) ]' \
  -c env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 "$@"
