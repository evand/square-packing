#!/usr/bin/env bash
# m8_zbisect.sh COVER LO HI [STEPS] -- bisect the smallest scale factor f in [LO,HI] at which zmx2 --d4 VERIFIES
# COVER x f (search/S60_COVER.md; copy of m7_zbisect.sh with ZCORES/ZTHREADS).  Prints each probe; the exact certified cost is total(COVER) x f_hi.
# Cores: $ZCORES (default 2-13), threads $ZTHREADS (default 12).  Scaled files and logs go to runs/m8zmx/.
set -eu
cd "$(dirname "$0")/.."
C=$1; LO=$2; HI=$3; N=${4:-10}
Z=verify2/target/release/zmx2; B=$(basename $C .txt); mkdir -p runs/m8zmx
for i in $(seq $N); do
  F=$(python3 -c "print(f'{($LO+$HI)/2:.7f}')")
  O=runs/m8zmx/${B}_x$F.txt
  python3 search/line_cover.py scale $C $O --factor $F > /dev/null
  R=$(taskset -c ${ZCORES:-2-13} $Z cert $O --d4 --threads ${ZTHREADS:-12} --uncert-cap 50 --log runs/m8zmx/${B}_x$F.log 2>&1 | grep -E "^(VERIFIED|NOT VERIFIED|INCOMPLETE)" | head -1)
  T=$($Z info $O | grep '^total' | sed 's/.*= //')
  echo "f=$F total=$T : $R"
  if [[ $R == VERIFIED* ]]; then HI=$F; rm -f $O; else LO=$F; rm -f $O; fi
done
echo "RESULT $C: verifies at f=$HI (total x f), refused at f=$LO"
