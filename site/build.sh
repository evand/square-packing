#!/usr/bin/env bash
# Rebuild the site's data from scratch (or incrementally).  Steps:
#   1. fetch    Ellsworth's HTML pages + every SVG they link (polite, sequential, skips files already present)
#   2. scrape   pages -> data/records_raw.json  (labels, s, attribution prose, polynomials)
#   3. parse    SVGs -> data/packings.json      (exact squares, validated: count / containment / overlap)
#   4. export   analysis -> www/data/p/*.json + www/data/index.json  (angles, contacts, rigidity, free regions)
#   5. timeline prose dates -> www/data/timeline.json
#   6. bounds   data/lower_bounds.json -> www/data (hand-curated with citations; see notes/lower-bounds-notes.md)
#   7. lists    the record packings n <= 324 from search/regularize/lists (certified, regularized) + records posted since
#               the catalogue (data/records2026.json) -> www/data/p/{best,alt,step}-*.json, index.json, timeline.json
#   provenance  (not in `all`: needs the packing store and register clone under search/packer/runs/) who posted each
#               record newer than the catalogue, and when -> data/records2026.json (committed)
# Usage: ./build.sh [fetch|scrape|parse|export|timeline|bounds|lists|provenance|all]   (default: all but fetch)
set -euo pipefail
cd "$(dirname "$0")"
step=${1:-site}
run() { echo "== $1"; shift; "$@"; }
case "$step" in
  fetch)    run "fetch pages" bash -c 'cd data/ellsworth && for f in squares_in_squares.html squares_in_squares__compared.html squares_in_squares__rigid.html "squares_in_squares__n^2-n-1.html" "squares_in_squares__Göbel_squares.html" "squares_in_squares__Göbel_strips.html" squares_in_squares__triangular_table.html squares_in_squares__analytic_minimization.html; do curl -sS -A "square-packing-explorer (+https://github.com/evand/square-packing)" -o "$f" "https://kingbird.myphotos.cc/packing/$f"; done; grep -oh "href=\"[^\"]*\.svg\"" *.html | sed "s/.*href=\"//;s/\"$//" | sort -u > svglist.txt; ./fetch.sh' ;;
  scrape)   run scrape python3 tools/scrape_pages.py ;;
  parse)    run parse python3 tools/parse_all.py ;;
  export)   run export python3 tools/export.py "${@:2}" ;;
  timeline) run timeline python3 tools/timeline.py && cp data/timeline.json www/data/ ;;
  bounds)   run bounds cp data/lower_bounds.json www/data/ ;;
  lists)    run lists python3 tools/import_lists.py "${@:2}" ;;
  provenance) run provenance python3 tools/provenance.py ;;
  site|all) [ "$step" = all ] && "$0" fetch; "$0" scrape; "$0" parse; "$0" export; "$0" timeline; "$0" bounds; "$0" lists ;;
  *) echo "unknown step $step"; exit 2 ;;
esac
