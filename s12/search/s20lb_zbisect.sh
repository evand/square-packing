#!/usr/bin/env bash
# s20lb_zbisect.sh COVER LO HI [STEPS] -- bisect the smallest scale factor f in [LO,HI] at which zmx2 --d4 VERIFIES
# COVER x f (search/S20_LB.md; m8_zbisect.sh without core pinning; extra zmx2 flags in $ZFLAGS).
# Prints each probe; the certified cost is total(COVER) x f_hi.  Threads $ZTHREADS (default 8).  Logs: runs/s20lb_zmx/.
set -eu
cd "$(dirname "$0")/.."
C=$1; LO=$2; HI=$3; N=${4:-10}
Z=verify2/target/release/zmx2; B=$(basename $C .txt); mkdir -p runs/s20lb_zmx
for i in $(seq $N); do
  F=$(python3 -c "print(f'{($LO+$HI)/2:.7f}')")
  O=runs/s20lb_zmx/${B}_x$F.txt
  rm -f runs/s20lb_zmx/${B}_x$F.log
  python3 search/line_cover.py scale $C $O --factor $F > /dev/null
  R=$($Z cert $O --d4 ${ZFLAGS:-} --threads ${ZTHREADS:-8} --uncert-cap 50 --log runs/s20lb_zmx/${B}_x$F.log 2>&1 | grep -E "^(VERIFIED|NOT VERIFIED|INCOMPLETE)" | head -1)
  T=$($Z info $O | grep '^total' | sed 's/.*= //')
  echo "f=$F total=$T : $R"
  if [[ $R == VERIFIED* ]]; then HI=$F; rm -f $O; else LO=$F; rm -f $O; fi
done
echo "RESULT $C: verifies at f=$HI (total x f), refused at f=$LO"
