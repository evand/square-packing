#!/usr/bin/env bash
# bisect_all.sh V1 V2 ... -- zbisect.sh on search/golf/covers/V.txt for each V, in [0.99, 1.07], 8 steps.
cd "$(dirname "$0")/../.."
for v in "$@"; do
  bash search/golf/zbisect.sh search/golf/covers/$v.txt ${LO:-0.990} ${HI:-1.070} ${STEPS:-8} 2>&1 | tee search/golf/zbis/$v.bisect | tail -1
done
