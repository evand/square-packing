#!/usr/bin/env bash
# zbisect.sh COVER LO HI [STEPS] -- smallest factor f in [LO,HI] at which zmx2 --d4 VERIFIES COVER x f
# (search/m8_zbisect.sh with golf's scaler and cores 6-11).  Run from anywhere; scratch in search/golf/zbis/.
set -eu
cd "$(dirname "$0")/../.."
C=$1; LO=$2; HI=$3; N=${4:-8}
Z=verify2/target/release/zmx2; B=$(basename $C .txt); mkdir -p search/golf/zbis
for i in $(seq $N); do
  F=$(python3 -c "print(f'{($LO+$HI)/2:.7f}')")
  O=search/golf/zbis/${B}_x$F.txt
  python3 search/golf/golf.py scale $C $O --factor $F > /dev/null
  R=$(taskset -c ${ZCORES:-6-11} $Z cert $O --d4 --threads ${ZTHREADS:-6} --uncert-cap 50 --log search/golf/zbis/${B}_x$F.log 2>&1 | grep -E "^(VERIFIED|NOT VERIFIED|INCOMPLETE)" | head -1)
  T=$($Z info $O | grep '^total' | sed 's/.*= //')
  echo "f=$F total=$T : $R"
  if [[ $R == VERIFIED* ]]; then HI=$F; else LO=$F; fi
  rm -f $O
done
echo "RESULT $C: verifies at f=$HI, refused at f=$LO"
