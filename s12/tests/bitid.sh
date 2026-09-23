#!/bin/sh
# Bit-identity of the verifier against a pre-change binary.  Everything without an anchor block
# must be byte-identical; a file WITH one must be identical up to the credited-clique columns the
# witness file gained on 2026-09-07 (search/WITNESS.md).  Usage: sh tests/bitid.sh <old-binary>
cd "$(dirname "$0")/.."
B=$1
V=verify/target/release/verify
T=runs/bitid; mkdir -p $T
fail=0
one() {   # name, args...
  n=$1; shift
  $B "$@" > $T/$n.base 2>&1; rcb=$?
  $V "$@" > $T/$n.new  2>&1; rcn=$?
  if [ "$rcb" != "$rcn" ]; then echo "  EXIT DIFFERS $n ($rcb vs $rcn)"; fail=1; return; fi
  if diff -q $T/$n.base $T/$n.new >/dev/null; then echo "  identical  $n"; else echo "  DIFFERS    $n"; diff $T/$n.base $T/$n.new | head -6; fail=1; fi
}
for c in s12_lower_3.9686 s12_lower_3.931795_sparse s12_56points_3.8 s12_boxclique_demo_3.9318_N2000 s12_lower_3.9676; do
  one "$c" certificates/$c.txt 12 2000 1 0
done
# witness mode (topk = 6): stdout and the witness file
for c in s12_lower_3.9686 s12_boxclique_demo_3.9318_N2000; do
  $B certificates/$c.txt 12 2000 1 6 $T/$c.sepbase > $T/$c.wbase 2>&1
  $V certificates/$c.txt 12 2000 1 6 $T/$c.sepnew  > $T/$c.wnew  2>&1
  if diff -q $T/$c.wbase $T/$c.wnew >/dev/null && diff -q $T/$c.sepbase $T/$c.sepnew >/dev/null
  then echo "  identical  $c (witness mode, $(wc -l < $T/$c.sepnew) witnesses)"
  else echo "  DIFFERS    $c (witness mode)"; fail=1; fi
done
# a branch trailer with per-box lambdas, in witness mode
cat certificates/s12_56points_3.8.txt > $T/br.txt
printf 'region corner 6 5\nlambda -5 -5 0 0\nk 0 0 1 1\n' >> $T/br.txt
$B $T/br.txt 12 2000 1 6 $T/br.sepbase > $T/br.wbase 2>&1
$V $T/br.txt 12 2000 1 6 $T/br.sepnew  > $T/br.wnew  2>&1
if diff -q $T/br.wbase $T/br.wnew >/dev/null && diff -q $T/br.sepbase $T/br.sepnew >/dev/null
then echo "  identical  per-box branch trailer (witness mode)"; else echo "  DIFFERS    per-box branch trailer"; fail=1; fi
# a real branch certificate
one "branch_k1" certificates/branch/s12_t3.98_corner_k1.txt 12 120 1 0
# An ANCHOR-CLIQUE certificate in witness mode is the one place the output is meant to change
# (search/WITNESS.md): stdout must still be identical, and the witness file must be the old one
# with the credited-clique columns appended, so cutting it back to five fields reproduces it.
awk 'NR==5{print $1, $2, 0; next} {print}' certificates/s12_56points_3.8.txt > $T/anc.txt
printf 'anchors 1 1\nanchorP 375 395 400\n0 1\npiece 0 0\n' >> $T/anc.txt
$B $T/anc.txt 12 200 1 6 $T/anc.sepbase > $T/anc.wbase 2>&1
$V $T/anc.txt 12 200 1 6 $T/anc.sepnew  > $T/anc.wnew  2>&1
cut -d' ' -f1-5 $T/anc.sepnew > $T/anc.sepcut
ncred=$(awk '{n+=$6} END{print n+0}' $T/anc.sepnew)
# a reference binary from after 2026-09-07 writes the two columns itself: then the files must be
# identical in full (the cut comparison is for a reference from before that date)
if diff -q $T/anc.wbase $T/anc.wnew >/dev/null && diff -q $T/anc.sepbase $T/anc.sepnew >/dev/null
then echo "  identical  anchor clique (stdout and the full witness file, $ncred credited-clique columns)"
elif diff -q $T/anc.wbase $T/anc.wnew >/dev/null && diff -q $T/anc.sepbase $T/anc.sepcut >/dev/null
then echo "  identical  anchor clique (stdout, and witnesses up to the $ncred credited-clique columns)"
else echo "  DIFFERS    anchor clique (witness mode)"; diff $T/anc.sepbase $T/anc.sepnew | head -4; fail=1; fi
# TIGHT_DUMP mode
TIGHT_DUMP=$T/td.base TIGHT_THRESH=10000000 TIGHT_MAX=20000 $B certificates/s12_56points_3.8.txt 12 200 1 0 > $T/td.wbase 2>&1
TIGHT_DUMP=$T/td.new  TIGHT_THRESH=10000000 TIGHT_MAX=20000 $V certificates/s12_56points_3.8.txt 12 200 1 0 > $T/td.wnew  2>&1
# the stdout names the dump file, so normalise the two paths before comparing
sed "s|$T/td.base|DUMP|" $T/td.wbase > $T/td.wbase2; sed "s|$T/td.new|DUMP|" $T/td.wnew > $T/td.wnew2
if diff -q $T/td.wbase2 $T/td.wnew2 >/dev/null && diff -q $T/td.base $T/td.new >/dev/null
then echo "  identical  TIGHT_DUMP ($(grep -c '^c ' $T/td.new) cells dumped)"; else echo "  DIFFERS    TIGHT_DUMP"; fail=1; fi
[ $fail -eq 0 ] && echo "ALL IDENTICAL" || echo "BIT-IDENTITY FAILED"
exit $fail
