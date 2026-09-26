#!/usr/bin/env python3
"""D4-reduced exact sweep of the s(32) candidate with `zmcheck cert --d4` (search/S32_EXACT.md sec 7).

The fundamental region is the 30 x 30 centre cells of the quadrant [0,3]^2 times u-bins 0..3
(u in [0,1/2]): 3,600 of the 28,800 root boxes.  It is split into jobs, one per x-cell column
and block of y-cells; each job is one `zmcheck cert --d4 --xlo a --xhi a --ylo b --yhi c` run
(ZM_ROOTLOG=1).  Every job re-checks the D4 invariance of the certificate itself.

    python3 search/s32_sweep.py plan    [--out DIR] [--yblock B]      # list the jobs, cost-ordered
    python3 search/s32_sweep.py run     [--out DIR] [--jobs J] [--threads T] [--only REGEX]
    python3 search/s32_sweep.py summary [--out DIR]                   # -> DIR/SUMMARY.txt

Run from s12/.  Resumable: a finished job leaves DIR/<job>.log (+ .err); a job whose .log is
present is skipped, an interrupted one (only .part files) is rerun from scratch.  The job
blocking (--yblock) is fixed per DIR at the first `run`/`plan` (DIR/manifest.json).
Checker settings are fixed: --branch-cap 640 --node-cap 4000000 --depth 22.
"""
import argparse, hashlib, json, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

CERT = 'runs/s32-close_candidate.txt'
ZM = 'verify2/target/release/zmcheck'
SETTINGS = ['--branch-cap', '640', '--node-cap', '4000000', '--depth', '22']
NCELL, NBIN = 30, 4            # quadrant [0,3]^2 at pitch 1/10; u-bins 0..3
GERM = {14, 15, 24, 25}        # cells touching the interior tile-centre lines x or y = 1.5, 2.5
INTL = {9, 10, 19, 20, 29}     # cells touching the lines x or y = 1, 2, 3


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def jobs(yblock):
    out = []
    for i in range(NCELL):
        for j0 in range(0, NCELL, yblock):
            j1 = min(j0 + yblock, NCELL) - 1
            # cost heuristic, only used for ordering (longest first); cells < 5 are EMPTY
            est = 0.0
            for j in range(j0, j1 + 1):
                if i < 5 or j < 5:
                    continue
                est += 1 + 8 * ((i in GERM) + (j in GERM)) + 3 * ((i in INTL) + (j in INTL)) + 40 * (i in GERM and j in GERM)
            out.append(dict(name=f'c{i:02d}_y{j0:02d}-{j1:02d}', i=i, j0=j0, j1=j1,
                            roots=(j1 - j0 + 1) * NBIN, est=est))
    out.sort(key=lambda r: -r['est'])
    return out


def manifest(a, write):
    os.makedirs(a.out, exist_ok=True)
    mp = os.path.join(a.out, 'manifest.json')
    if os.path.exists(mp):
        m = json.load(open(mp))
        if a.yblock is not None and a.yblock != m['yblock']:
            sys.exit(f'{mp} has yblock {m["yblock"]}; not {a.yblock} (use another --out)')
        return m
    m = dict(yblock=a.yblock or 5, cert=CERT, cert_sha256=sha(CERT), settings=SETTINGS)
    if write:
        json.dump(m, open(mp, 'w'), indent=1)
    return m


def cmd(job, threads, zm):
    return [zm, 'cert', CERT, '--d4', '--threads', str(threads), *SETTINGS,
            '--xlo', f'{job["i"] / 10:.1f}', '--xhi', f'{job["i"] / 10:.1f}',
            '--ylo', f'{job["j0"] / 10:.1f}', '--yhi', f'{job["j1"] / 10:.1f}']


def run_job(job, a):
    base = os.path.join(a.out, job['name'])
    c = cmd(job, a.threads, a.zmcheck)
    t0 = time.time()
    env = dict(os.environ, ZM_ROOTLOG='1')
    with open(base + '.log.part', 'w') as fo, open(base + '.err.part', 'w') as fe:
        fo.write('# ' + ' '.join(c) + '\n')
        fo.flush()
        rc = subprocess.run(c, stdout=fo, stderr=fe, env=env).returncode
    txt = open(base + '.log.part').read()
    ok = rc == 0 and 'done in' in txt and 'NOT VERIFIED' in txt
    if ok:
        os.replace(base + '.err.part', base + '.err')
        os.replace(base + '.log.part', base + '.log')   # the .log is the completion marker
    m = re.search(r'UNCERTIFIED (\d+)', txt)
    print(f'{time.strftime("%H:%M:%S")} {job["name"]}: {"done" if ok else f"FAILED rc={rc}"} '
          f'{time.time() - t0:.0f}s uncert {m[1] if m else "?"}', flush=True)
    return ok


def parse(a, job):
    """status of one job from its files: None if not finished, else a dict."""
    base = os.path.join(a.out, job['name'])
    if not os.path.exists(base + '.log'):
        return None
    txt = open(base + '.log').read()
    err = open(base + '.err').read() if os.path.exists(base + '.err') else ''
    g = lambda pat: (re.search(pat, txt) or [None, None])[1]
    r = dict(roots=int(g(r'(\d+) root boxes') or -1), depth=g(r'depth limit (\d+)'),
             uncert=int(g(r'UNCERTIFIED (\d+)') or -1), boxes=int(g(r'boxes (\d+), max depth') or 0),
             maxdepth=int(g(r'max depth (\d+)') or 0), wall=float(g(r'done in (\d+)s') or 0),
             d4=bool(re.search(r'^D4: weight function invariant', txt, re.M)),
             reduced='D4-REDUCED' in txt, refused=int(g(r'(\d+) boxes REFUSED') or 0),
             cpu=sum(float(x) for x in re.findall(r'^ROOT .* ([\d.]+)s$', err, re.M)),
             uroots=[l[5:] for l in err.splitlines() if l.startswith('ROOT ') and
                     not l.split(' uncert ')[1].startswith('0 ')],
             ulines=[l[7:] for l in err.splitlines() if l.startswith('UNCERT ')])
    r['problems'] = [p for p, bad in [
        (f'root count {r["roots"]} != {job["roots"]}', r['roots'] != job['roots']),
        ('no D4 invariance line', not r['d4']), ('not a --d4 run', not r['reduced']),
        (f'depth limit {r["depth"]} != 22', r['depth'] != '22'),
        ('no census', r['uncert'] < 0)] if bad]
    return r


def cmd_plan(a):
    m = manifest(a, write=False)
    js = jobs(m['yblock'])
    for j in js:
        print(f'{j["name"]}  roots {j["roots"]:3d}  est {j["est"]:5.1f}  ' + ' '.join(cmd(j, a.threads, a.zmcheck)[3:]))
    print(f'{len(js)} jobs, {sum(j["roots"] for j in js)} roots')


def cmd_run(a):
    m = manifest(a, write=True)
    if sha(CERT) != m['cert_sha256']:
        sys.exit(f'{CERT} changed since {a.out}/manifest.json was written')
    js = [j for j in jobs(m['yblock']) if not a.only or re.search(a.only, j['name'])]
    todo = [j for j in js if not os.path.exists(os.path.join(a.out, j['name'] + '.log'))]
    zsha = sha(a.zmcheck)
    with open(os.path.join(a.out, 'binaries.txt'), 'a') as f:
        f.write(f'{time.strftime("%Y-%m-%d %H:%M:%S")} {a.zmcheck} sha256 {zsha} jobs {len(todo)}\n')
    print(f'{len(js)} jobs selected, {len(js) - len(todo)} already done, {len(todo)} to run '
          f'({a.jobs} at a time x {a.threads} threads); zmcheck sha256 {zsha[:16]}', flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(lambda j: run_job(j, a), todo))
    print(f'{sum(res)} finished, {len(res) - sum(res)} failed')
    cmd_summary(a)


def cmd_summary(a):
    m = manifest(a, write=False)
    js = sorted(jobs(m['yblock']), key=lambda j: j['name'])
    lines, unc, pend, bad = [], [], [], []
    tot = dict(roots=0, uncert=0, boxes=0, cpu=0.0, maxdepth=0, refused=0)
    for j in js:
        r = parse(a, j)
        if r is None:
            pend.append(j['name'])
            continue
        if r['problems']:
            bad.append(f'{j["name"]}: ' + '; '.join(r['problems']))
        for k in ('roots', 'uncert', 'boxes', 'cpu', 'refused'):
            tot[k] += r[k]
        tot['maxdepth'] = max(tot['maxdepth'], r['maxdepth'])
        if r['uncert']:
            unc.append(f'## {j["name"]}: {r["uncert"]} uncertified')
            unc += ['   root ' + u for u in r['uroots']]
            unc += ['   ' + u for u in r['ulines']]
    need = sum(j['roots'] for j in js)
    lines.append(f'D4-reduced s(32) sweep in {a.out}: {len(js) - len(pend)}/{len(js)} jobs done, '
                 f'roots {tot["roots"]}/{need}')
    lines.append(f'certificate {m["cert"]} sha256 {m["cert_sha256"]}; settings {" ".join(m["settings"])}')
    lines.append(f'boxes {tot["boxes"]}, max depth {tot["maxdepth"]}, CPU {tot["cpu"] / 3600:.1f} h '
                 f'(sum of ROOT times), uncertified {tot["uncert"]}, refused {tot["refused"]}')
    bp = os.path.join(a.out, 'binaries.txt')
    if os.path.exists(bp):
        shas = {l.split()[4] for l in open(bp) if 'sha256' in l}
        lines.append(f'zmcheck binaries used: {len(shas)} ({", ".join(s[:16] for s in shas)})')
    if bad:
        lines += ['PROBLEMS:'] + ['  ' + b for b in bad]
    if pend:
        lines.append(f'pending ({len(pend)}): ' + ' '.join(pend[:12]) + (' ...' if len(pend) > 12 else ''))
    if not pend and not bad and tot['uncert'] == 0 and tot['roots'] == need == NCELL * NCELL * NBIN:
        lines.append('D4 SWEEP CLEAN: all 3600 roots of the fundamental region [0,3]^2 x u in [0,1/2] '
                     'certified by zmcheck --d4 (D4 invariance checked exactly in every job) => every '
                     'closed unit square in [0,6]^2 captures weight >= 1 (S32_EXACT.md sec 7).')
    else:
        lines.append('NOT (YET) CLEAN')
    if unc:
        lines += ['', 'UNCERTIFIED BOXES (per root: at most 80 listed):'] + unc
    txt = '\n'.join(lines) + '\n'
    open(os.path.join(a.out, 'SUMMARY.txt'), 'w').write(txt)
    print(txt, end='')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['plan', 'run', 'summary'])
    ap.add_argument('--out', default='runs/s32d4')
    ap.add_argument('--yblock', type=int, default=None, help='y-cells per job (default 5; fixed per --out)')
    ap.add_argument('--jobs', type=int, default=1, help='concurrent zmcheck processes')
    ap.add_argument('--threads', type=int, default=4, help='threads per zmcheck process')
    ap.add_argument('--only', default=None, help='regex on job names (e.g. "c2[45]_")')
    ap.add_argument('--zmcheck', default=ZM, help='a zmcheck build with --d4')
    a = ap.parse_args()
    if not os.path.exists(CERT):
        sys.exit(f'run from s12/ ({CERT} not found)')
    {'plan': cmd_plan, 'run': cmd_run, 'summary': cmd_summary}[a.mode](a)


if __name__ == '__main__':
    main()
