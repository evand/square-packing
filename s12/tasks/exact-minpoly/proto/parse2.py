import re, json, sympy as sp
from mpmath import mp, mpf, polyroots
mp.dps = 80
h = open('main.html').read()
boxes = h.split('<div class="box">')[1:]
ours = {r['n']: r for r in json.load(open('/home/evand/math/square-packing/public/s12/search/exact/batch/results.json'))}
mine = {}
for l in open('/home/evand/math/square-packing/public/s12/tasks/exact-minpoly/results/findpoly68_tilted.txt'):
    p = l.split(None, 4)
    mine[int(p[0])] = p[4].strip().split(']')[0] + ']' if '[' in p[4] else None
s = sp.Symbol('s')
out = []
for b in boxes:
    m = re.match(r'\s*<font size="\+3">(\d+)<br>', b)
    if not m: continue
    n = int(m.group(1))
    fr = re.search(r'<span class="frames">\$(.*?)\$</span>', b, re.S)
    f1 = re.search(r'<span class="frame1">\$s = \{\}\^\{(\d+)\}', b)
    num = re.search(r'\\Nn\{([\d.]+)\}', b)
    rec = {'n': n, 'deg': int(f1.group(1)) if f1 else None, 'num': num.group(1) if num else None}
    if fr:
        txt = fr.group(1).split('=')[0]
        txt = re.sub(r's\^\{?(\d+)\}?', r's**\1', txt)
        txt = txt.replace('\\sqrt{2}', '*sqrt(2)').replace('\\sqrt2', '*sqrt(2)')
        txt = re.sub(r'(\d)s', r'\1*s', txt); txt = txt.replace(')s', ')*s').replace('(*', '(').replace('+*', '+').replace('-*', '-')
        try:
            E = sp.expand(sp.sympify(txt, locals={'s': s}))
            if E.has(sp.sqrt(2)):
                E = sp.expand(E * E.subs(sp.sqrt(2), -sp.sqrt(2)))     # norm to Q
                rec['sqrt2'] = True
            P = sp.Poly(E, s)
            fl = sp.factor_list(P)[1]
            rec['factors'] = [[int(c) for c in q.all_coeffs()] for q, m in fl]
            rec['poly'] = [int(c) for c in P.all_coeffs()]
        except Exception as e:
            rec['err'] = str(e)[:80]
    out.append(rec)
json.dump(out, open('kingbird.json', 'w'))
print(len(out), 'entries;', sum('poly' in r for r in out), 'with polynomial;', sum('err' in r for r in out), 'parse errors')
# compare
for r in out:
    if 'poly' not in r: continue
    n = r['n']; o = ours.get(n)
    if not o or not o.get('S_exact'): print(n, 'deg', len(r['poly'])-1, 'no exact S of ours'); continue
    S = mpf(o['S_exact'])
    best = None
    for c in r['factors']:
        v = sum(ci * S**(len(c)-1-i) for i, ci in enumerate(c)); sc = sum(abs(ci) * abs(S)**(len(c)-1-i) for i, ci in enumerate(c))
        if best is None or abs(v)/sc < best[0]: best = (abs(v)/sc, c)
    v, c = best[0]*1, best[1]; sc = 1
    print(n, 'Q-deg of matching factor', len(c)-1, 'sqrt2' if r.get('sqrt2') else '', 'height', len(str(max(map(abs, c)))), 'rel residual at our S', mp.nstr(abs(v)/sc, 3), 'ours:', mine.get(n, '-'))
