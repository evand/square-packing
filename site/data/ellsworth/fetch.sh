#!/usr/bin/env bash
# Polite sequential fetch of Ellsworth's SVGs; skips files already present.
cd "$(dirname "$0")"
ok=0; fail=0
while read -r f; do
  [ -z "$f" ] && continue
  out="svg/$f"
  [ -s "$out" ] && continue
  if curl -sS -f -A "square-packing-explorer (+https://github.com/evand/square-packing)" -o "$out" "https://kingbird.myphotos.cc/packing/$f"; then ok=$((ok+1)); else fail=$((fail+1)); echo "FAIL $f"; rm -f "$out"; fi
  sleep 0.4
done < svglist.txt
echo "done ok=$ok fail=$fail total=$(ls svg | wc -l)"
