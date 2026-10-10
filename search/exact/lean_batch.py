"""Batch Lean proofs of `Packs n S*` (S* exact) for the minpoly-verified register records.

For each n: lean_cert.py --split (certificate module + chunk modules + glue), then `lake build` of the data module,
the chunks (in parallel, each under a systemd memory cap), the glue, and an axiom check.  Generated Lean files go to
lean/Sqpack/Exact/Packs/ (committed since 2026-10-10, so each theorem can be cited at a commit; not in the default
build: `lake build Sqpack.Exact.Packs.N<n>`).  Input: the committed exact forms minpoly/data/n-N.minpoly.json.gz.  Results are appended to
lean_batch/results.tsv (and, for quadratic S*, the closed form to lean_batch/closed.tsv: theorems `packs_closed`,
`minSide_le_closed`); `--summary` writes lean_batch/results.md.

  ./env.sh python3 lean_batch.py --maxdeg 20 --maxn 324 [--jobs 3] [--chunk-jobs 2] [--mem 14G] [--only 5,11,294]
  ./env.sh python3 lean_batch.py --summary
"""
import argparse, json, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
LEAN = os.path.normpath(os.path.join(HERE, '..', '..', 'lean'))
OUTD = os.path.join(HERE, 'lean_batch')
TSV = os.path.join(OUTD, 'results.tsv')
CLOSED = os.path.join(OUTD, 'closed.tsv')   # n, closed form of S* (quadratic S*): packs_closed, minSide_le_closed
STD_AXIOMS = '[propext, Classical.choice, Quot.sound]'


def records(maxdeg, maxn):
    R = json.load(open(os.path.join(HERE, 'minpoly', 'results.json')))
    out = []
    for r in R:
        n = r['n']
        if r.get('status') != 'VALID' or n > maxn:
            continue
        if (r.get('field_degree') or 1) <= maxdeg:
            out.append((n, r.get('field_degree') or 1))
    return sorted(out)


def run(cmd, mem, timeout, cwd, log):
    cmd = ['systemd-run', '--user', '--scope', '-q', '-p', f'MemoryMax={mem}', '-p', 'MemorySwapMax=0',
           'timeout', str(timeout)] + cmd
    with open(log, 'a') as f:
        f.write(f'$ {" ".join(cmd)}\n')
        f.flush()
        return subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT).returncode


def one(n, deg, A):
    t0 = time.time()
    log = os.path.join(OUTD, 'logs', f'n-{n}.log')
    open(log, 'w').close()
    mod = f'N{n}'
    out = os.path.join(LEAN, 'Sqpack', 'Exact', 'Packs', f'{mod}.lean')
    status, chunks = 'ok', 0
    rc = run([sys.executable, os.path.join(HERE, 'lean_cert.py'), os.path.join(HERE, 'minpoly', 'data', f'n-{n}.minpoly.json.gz'),
              out, '--split', str(A.chunk)], '12G', 3000, HERE, log)
    tgen = time.time() - t0
    if rc != 0:
        status = f'gen failed ({rc})'
    else:
        chunks = len([f for f in os.listdir(out[:-5]) if re.fullmatch(r'C\d+\.lean', f)])
        base = f'Sqpack.Exact.Packs.{mod}'
        if run(['lake', 'build', f'{base}.Data'], A.mem, 3600, LEAN, log) != 0:
            status = 'data failed'
        else:
            with ThreadPoolExecutor(A.chunk_jobs) as ex:
                rcs = list(ex.map(lambda c: run(['lake', 'build', f'{base}.C{c}'], A.mem, 3600, LEAN, log), range(chunks)))
            bad = [c for c, r in enumerate(rcs) if r != 0]
            if bad:
                status = f'chunks failed {bad} (rc {[rcs[c] for c in bad]})'
            elif run(['lake', 'build', base], A.mem, 3600, LEAN, log) != 0:
                status = 'glue failed'
            else:
                ax = os.path.join(OUTD, 'logs', f'ax-{n}.lean')
                src = open(out).read()
                m = re.search(r'theorem minSide_le_closed : minSide \d+ ≤ (.*) :=', src)
                thms = ['packs'] + (['packs_closed', 'minSide_le_closed'] if m else [])
                open(ax, 'w').write(f'import {base}\n' + ''.join(f'#print axioms UnitSquarePacking.EC.{mod}.{t}\n'
                                                                  for t in thms))
                r = subprocess.run(['lake', 'env', 'lean', ax], cwd=LEAN, capture_output=True, text=True)
                if r.returncode != 0 or r.stdout.count(STD_AXIOMS) != len(thms):
                    status = 'AXIOMS: ' + r.stdout.strip().replace('\n', ' ')[:200]
                elif m:
                    with open(CLOSED, 'a') as f:
                        f.write(f'{n}\t{m.group(1)}\n')
    ttot = time.time() - t0
    line = f'{n}\t{deg}\t{chunks}\t{tgen:.0f}\t{ttot:.0f}\t{status}\n'
    with open(TSV, 'a') as f:
        f.write(line)
    print(line, end='', flush=True)


def summary():
    rows = {}
    for line in open(TSV):
        n, deg, ch, tg, tt, st = line.rstrip('\n').split('\t')
        rows[int(n)] = (deg, ch, tg, tt, st)              # latest run wins
    ok = [n for n, r in rows.items() if r[4] == 'ok']
    L = ['# Lean `Packs n S*` batch (lean_batch.py)', '',
         f'Proved (kernel-checked, standard axioms): **{len(ok)}** of {len(rows)} attempted.', '',
         '| n | field deg | chunks | gen s | total s | status |', '|---|---|---|---|---|---|']
    for n in sorted(rows):
        L.append('| ' + ' | '.join([str(n)] + list(rows[n])) + ' |')
    if os.path.exists(CLOSED):
        cl = {}
        for line in open(CLOSED):
            n, e = line.rstrip('\n').split('\t')
            cl[int(n)] = e                                # latest run wins
        cl = {n: e for n, e in cl.items() if n in rows and rows[n][4] == 'ok'}
        L += ['', '## Closed forms (quadratic `S*`)', '',
              f'**{len(cl)}** records with `UnitSquarePacking.EC.N<n>.packs_closed : Packs n (S*)` and '
              '`minSide_le_closed : minSide n ≤ S*`, `S*` written as `a + b √d` (`Sqpack/ExactQuad.lean`).', '',
              '| n | S* |', '|---|---|']
        L += [f'| {n} | `{cl[n]}` |' for n in sorted(cl)]
    open(os.path.join(OUTD, 'results.md'), 'w').write('\n'.join(L) + '\n')
    print(f'{len(ok)} ok of {len(rows)}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--maxdeg', type=int, default=8)
    ap.add_argument('--maxn', type=int, default=324)
    ap.add_argument('--only', default='')
    ap.add_argument('--jobs', type=int, default=3, help='records in parallel')
    ap.add_argument('--chunk-jobs', type=int, default=2, help='chunks in parallel per record')
    ap.add_argument('--chunk', type=int, default=25, help='squares per chunk')
    ap.add_argument('--mem', default='14G')
    ap.add_argument('--skip-done', action='store_true')
    ap.add_argument('--summary', action='store_true')
    A = ap.parse_args()
    os.makedirs(os.path.join(OUTD, 'logs'), exist_ok=True)
    if A.summary:
        return summary()
    R = records(A.maxdeg, A.maxn)
    if A.only:
        want = {int(x) for x in A.only.split(',')}
        R = [r for r in R if r[0] in want]
    if A.skip_done and os.path.exists(TSV):
        done = {int(l.split('\t')[0]) for l in open(TSV) if l.rstrip('\n').endswith('\tok')}
        R = [r for r in R if r[0] not in done]
    print(f'{len(R)} records', flush=True)
    with ThreadPoolExecutor(A.jobs) as ex:
        list(ex.map(lambda r: one(r[0], r[1], A), R))
    summary()


if __name__ == '__main__':
    main()
