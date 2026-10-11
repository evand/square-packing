#!/usr/bin/env python3
"""cvc5_run.py FILE.smt2 [--timeout s]: decide an SMT-LIB QF_NRA problem with cvc5's cylindrical algebraic
coverings (nl-cov), as an engine independent of z3's NLSAT.  Prints sat/unsat/unknown and seconds."""
import sys, time, argparse
import cvc5

a = argparse.ArgumentParser()
a.add_argument('file')
a.add_argument('--timeout', type=float, default=600)
o = a.parse_args()
text = open(o.file).read()
text = text.replace('(set-info :status unknown)', '(set-logic QF_NRA)', 1)
tm = cvc5.TermManager()
slv = cvc5.Solver(tm)
for k, v in [('produce-models', 'true'), ('nl-cov', 'true'), ('tlimit-per', str(int(o.timeout * 1000)))]:
    slv.setOption(k, v)
p = cvc5.InputParser(slv)
p.setStringInput(cvc5.InputLanguage.SMT_LIB_2_6, text, 'q')
sm = p.getSymbolManager()
t0 = time.time()
res = None
while True:
    cmd = p.nextCommand()
    if cmd.isNull():
        break
    out = cmd.invoke(slv, sm)
    if 'check-sat' in str(cmd):
        res = out.strip()
if res is None:
    r = slv.checkSat()
    res = str(r)
print(dict(file=o.file, verdict=res, secs=round(time.time() - t0, 2)), flush=True)
if res == 'sat':
    for name in ('X', 'Y', 'c', 's', 'lam'):
        try:
            t = sm.getDeclaredTerms()
        except Exception:
            t = []
    for t in sm.getDeclaredTerms():
        print(' ', t, '=', slv.getValue(t))
