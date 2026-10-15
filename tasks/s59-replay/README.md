# jlevy#279: s(59) exact replay (paused 10-07; moved out of TODO.md 10-15)

Evan 10-07: poor RoI on CPU for now; resume when cores are idle.  One `zm_mixed` sweep at depth 34 over the whole D4
region (wand125 needed two runs whose joint coverage jlevy couldn't check).  33,725 / 102,400 roots done, ~127 CPU-h
total; the rest is the expensive part.

Resume:

    python3 search/zm_mixed.py cert runs/s59_wand125/n59_mixed_cover_8.txt --d4 --cert-mode --disj --chain-from 0 \
        --depth 34 --pitch 1/20 --ubins 16 --nproc 15 --progress 5000 \
        --resume runs/s59_wand125/zm_mixed_d4_depth34.jsonl --manifest runs/s59_wand125/zm_mixed_d4_depth34_manifest.json

Then reply with draft (D) in `outreach/draft-jlevy-replies-2026-10-07.md` (posting needs Evan's OK).
