"""Table of minpoly.py / verify_exact.py results over the register batch, with cross-checks.

  python3 minpoly/summary.py            (from search/exact; writes minpoly/results.md and minpoly/results.json)

Cross-checks: the integer-relation polynomials found from 68 digits (tasks/exact-minpoly/results/findpoly68_tilted.txt,
equal up to sign) and Ellsworth's published degrees (tasks/exact-minpoly/results/kingbird_crosscheck.json, 'match').
"""
import json, glob, os, re

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.dirname(here))
task = os.path.join('..', '..', 'tasks', 'exact-minpoly', 'results')
fp = {}
for line in open(os.path.join(task, 'findpoly68_tilted.txt')):
    m = re.match(r'(\d+) .*?(\[[-\d, ]+\]|None)\s*$', line.strip())
    if m and m.group(2) != 'None':
        fp[int(m.group(1))] = json.loads(m.group(2))
kb = {n: d for n, d in json.load(open(os.path.join(task, 'kingbird_crosscheck.json')))['match']}

rows = []
for n in range(1, 325):
    b = f'minpoly/solve/n-{n}'
    r = {'n': n}
    if not os.path.exists(b + '.contacts.json'):
        r['status'] = 'no exact point' if os.path.exists(b + '.log') else 'not run'
        rows.append(r)
        continue
    log = open(b + '.minpoly.log').read() if os.path.exists(b + '.minpoly.log') else ''
    if not os.path.exists(b + '.minpoly.json'):
        m = re.search(r'not handled: (.*)', log)
        r['status'] = 'open: ' + (m.group(1)[:90] if m else (log.strip().splitlines() or ['?'])[-2 if 'exit' in log else -1][:90])
        rows.append(r)
        continue
    d = json.load(open(b + '.minpoly.json'))
    v = open(b + '.verify.txt').read() if os.path.exists(b + '.verify.txt') else ''
    st = d['stats']
    p = d['S']['poly']
    r.update(status='VALID' if 'VALID.' in v and 'NOT' not in v else 'verify failed', classes=st['n_classes_tilted'],
             free_classes=st.get('n_classes_free'), field_degree=st['field_degree'], S_degree=st['S_degree'],
             height_digits=len(str(max(abs(x) for x in p))), stationary=d.get('stationary'), seconds=st['seconds'],
             S=d['S']['value'][:20], poly=p if len(p) <= 9 else None)
    x = []
    if n in fp:
        q = fp[n]
        x.append('findpoly ok' if (q == p[::-1] or [-c for c in q] == p[::-1]) else 'findpoly DIFFERS')
    if n in kb:
        x.append('Ellsworth deg ok' if kb[n] == st['S_degree'] else f'Ellsworth deg {kb[n]} DIFFERS')
    r['crosscheck'] = ', '.join(x)
    rows.append(r)

json.dump(rows, open('minpoly/results.json', 'w'), indent=0)
from collections import Counter
c = Counter(r['status'].split(':')[0] for r in rows)
with open('minpoly/results.md', 'w') as fh:
    fh.write('# minpoly results (register batch)\n\n' + ', '.join(f'{k}: {v}' for k, v in sorted(c.items())) + '\n\n')
    fh.write('| n | status | tilted classes (free) | field deg | S deg | height digits | stationary | s | cross-check |\n')
    fh.write('|---|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        if 'classes' in r:
            fh.write(f"| {r['n']} | {r['status']} | {r['classes']} ({r['free_classes']}) | {r['field_degree']} | "
                     f"{r['S_degree']} | {r['height_digits']} | {'yes' if r['stationary'] else ''} | {r['seconds']} | "
                     f"{r['crosscheck']} |\n")
        else:
            fh.write(f"| {r['n']} | {r['status']} | | | | | | | |\n")
print(', '.join(f'{k}: {v}' for k, v in sorted(c.items())))
bad = [r['n'] for r in rows if 'DIFFERS' in r.get('crosscheck', '') or r['status'] == 'verify failed']
print('cross-check differences / verify failures:', bad)
