#!/usr/bin/env bash
# python 3.13 + cyipopt (IPOPT/MUMPS), scipy, numpy, mpmath from nixpkgs
set -euo pipefail
exec nix --extra-experimental-features 'nix-command flakes' shell --impure --expr \
  'with import (builtins.getFlake "nixpkgs") {}; [ (python313.withPackages (p: with p; [ cyipopt scipy numpy mpmath ])) ]' \
  -c env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$@"
