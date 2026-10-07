#!/usr/bin/env bash
# minpoly.py + verify_exact.py on every solved register packing (minpoly/solve/n-*.contacts.json); run inside ./env.sh.
# Per n: minpoly/solve/n-N.minpoly.log (minpoly.py) and n-N.verify.txt (verify_exact.py).  Skips n already done
# (a .minpoly.log exists) unless FORCE=1.
cd "$(dirname "$0")/.."
ls minpoly/solve/n-*.contacts.json | sed 's/.contacts.json//' | sort -t- -k2 -n | xargs -P "${PROCS:-6}" -I{} sh -c '
  [ -z "$FORCE" ] && [ -f {}.minpoly.log ] && exit 0
  rm -f {}.minpoly.json {}.verify.txt; timeout ${TMO:-600} python3 minpoly.py {} > {}.minpoly.log 2>&1; echo "exit $?" >> {}.minpoly.log
  [ -f {}.minpoly.json ] && timeout ${TMO:-600} python3 verify_exact.py {}.minpoly.json > {}.verify.txt 2>&1'
