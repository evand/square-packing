# source search/packer/activate.sh     (from anywhere; bash)
# The project's own Python environment (11-10): ~/math/square-packing/venv, built from the repo's requirements.lock
# (exact versions; requirements.txt lists the direct dependencies).  Builds the venv if it is missing.
# Start the pool, jobs and sessions from a shell that sourced this, so every unit runs on the pinned versions.
_PK_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_PK_REPO="$(cd "$_PK_HERE/../.." && pwd)"
_PK_VENV="${PK_VENV:-$(cd "$_PK_REPO/.." && pwd)/venv}"
# NixOS: pip wheels with C extensions (numpy, scipy, highspy) need nix-ld's libstdc++ / zlib.
_NIX_LD_LIB=/run/current-system/sw/share/nix-ld/lib
if [[ -d "$_NIX_LD_LIB" && ":$LD_LIBRARY_PATH:" != *":$_NIX_LD_LIB:"* ]]; then
  export LD_LIBRARY_PATH="$_NIX_LD_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
if [[ ! -x "$_PK_VENV/bin/python3" ]]; then
  echo "activate.sh: building $_PK_VENV from requirements.lock" >&2
  /run/current-system/sw/bin/python3 -m venv "$_PK_VENV" && "$_PK_VENV/bin/pip" install -q -r "$_PK_REPO/requirements.lock"
fi
[[ -n "$VIRTUAL_ENV" && "$VIRTUAL_ENV" != "$_PK_VENV" ]] && type deactivate >/dev/null 2>&1 && deactivate
. "$_PK_VENV/bin/activate"
unset _PK_HERE _PK_REPO _NIX_LD_LIB
