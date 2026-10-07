#!/usr/bin/env bash
# run1.sh N exactfile rattlers-json tol-exp
H=/tmp/claude-1000/-home-evand-math/e6f049c2-2e76-4498-80cb-37969d16b554/scratchpad; cd $H/sym
n=$1; f=$2; r=$3; te=$4
python3 build2.py $f n$n.ms "$r" '{}' $te 2> n$n.build.log || { tail -3 n$n.build.log; exit 1; }
python3 reduce.py n$n.ms n${n}r.ms 2> n$n.reduce.log; head -1 n$n.reduce.log
python3 - <<PY
L=open('n${n}r.ms').read().split('\n',2); v=L[0].split(','); v=['S']+[x for x in v if x!='S']
open('n${n}s.ms','w').write(','.join(v[::-1])+'\n'+L[1]+'\n'+L[2])
PY
timeout ${TMO:-600} nix --extra-experimental-features 'nix-command flakes' shell nixpkgs#msolve -c msolve -t 2 -P 1 -f n${n}s.ms -o n${n}s.out > n$n.msolve.log 2>&1; tail -1 n$n.msolve.log
S=$(head -1 $f | awk '{print substr($2,1,60)}')
$H/venv/bin/python param.py n${n}s.out $S 2>&1 | grep -v "^MINPOLY" ; grep -o "MINPOLY.*" <($H/venv/bin/python param.py n${n}s.out $S 2>/dev/null) | cut -c1-200
