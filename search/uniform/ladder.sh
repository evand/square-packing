#!/bin/sh
# Launch a ladder of ILP runs: ladder.sh FINE TL "k list" "s list" [extra args]
# logs go to runs/log_<tag>.txt, tag = k<k>_s<s>_f<fine>
cd "$(dirname "$0")/../.." || exit 1
FINE=$1; TL=$2; KS=$3; SS=$4; shift 4
for k in $KS; do
  for s in $SS; do
    if [ "$k" = 1 ]; then SYM=none; else SYM=d4; fi
    tag=k${k}_s${s}_f${FINE}${TAGSUF}
    nohup python3 search/uniform/ilp_uniform.py "$k" "$s" --fine "$FINE" --sym $SYM --tag "$tag" --tl "$TL" --rounds 40 --threads 2 "$@" > "runs/log_$tag.txt" 2>&1 &
  done
done
echo launched
