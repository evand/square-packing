#!/usr/bin/env python3
"""Resumable per-root runner for search/zeromargin.py over the D4 fundamental region
(S32_EXACT.md sec 11).  An exact re-check, independent of verify2/zmcheck.

    python3 search/zm_d4_sweep.py run     --cert runs/s32-close_candidate.txt --out runs/s32py_d4 \\
                                          [--nproc 8 | --nproc-auto 8:28] [--depth 24]
    python3 search/zm_d4_sweep.py summary --out runs/s32py_d4

run
  * On first use of --out it FREEZES the checker: search/zeromargin.py is copied to
    OUT/checker/zeromargin.py and every worker imports that copy, so later edits to search/ cannot
    leak into a running or resumed sweep.  OUT/manifest.json records the certificate's path and
    sha256, the checker copy's sha256, this runner's sha256 and the settings (depth, theta_bias, ...).
    A later `run` refuses a different certificate, checker copy or setting.
  * Each worker re-reads the certificate, checks its sha256, and checks EXACTLY that the weighted
    point set is D4-invariant (Checker.symmetric_d4); if not, the sweep stops.
  * One task per root box (zeromargin.d4_roots: [0,m/2]^2 cells of pitch 1/10 x 8 bins of
    u in [0,1/2]).  A finished root is appended (flushed + fsync'd) to OUT/roots.jsonl with its
    census, uncertified boxes, CPU time, and the checker and certificate sha256.  A rerun skips
    finished roots and redoes the ones in flight.  There is NO time cap: a root either finishes
    (certified, or with uncertified boxes at the depth limit) or is not in the file.
  * --deepen D2 (D2 > the base depth): rerun, at depth D2, every root whose records so far all
    have uncertified boxes.  Its records carry depth D2.  Each record is a complete, independent
    run of its root, so a root is certified if ANY valid record of it (checker and certificate
    sha256 = the manifest's, depth >= the base depth) has no uncertified box.
  * --nproc-auto A:B: at most A roots in flight while another sweep (`pgrep -f` of --busy-pattern)
    is running, B once it has finished.  Heaviest-first order (the germ roots' u-bin 0 first).

summary
  Writes OUT/SUMMARY.txt: one line per root of the region (verdict, boxes, depth, census, CPU,
  checker sha, cert sha), the uncertified boxes, and the totals; prints `D4 RECHECK CLEAN` only if
  every root of the region is present, every record carries the manifest's checker and certificate
  sha256 and depth, and no root has an uncertified box.  Anything missing is INCOMPLETE; anything
  uncertified is NOT CLEAN.
"""
import sys, os, json, time, argparse, hashlib, shutil, subprocess, importlib.util, math
from fractions import Fraction as F
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait

HERE = os.path.dirname(os.path.abspath(__file__))
SETTINGS = dict(depth=24, theta_bias=1, use_chain=True, chain_from=0, use_adm=True, clip=True,
                pitch='1/10', ubins=8, reduction='D4')

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

def load_checker(out):
    p = os.path.join(out, 'checker', 'zeromargin.py')
    spec = importlib.util.spec_from_file_location('zeromargin_frozen', p)
    Z = importlib.util.module_from_spec(spec); spec.loader.exec_module(Z)
    return Z

def rkey(root): return ",".join(str(v) for v in root)

# ------------------------------------------------------------------ worker
G = {}
def winit(out, cert, cert_sha, chk_sha, depth, theta_bias):
    os.nice(0)
    Z = load_checker(out)
    if sha256(Z.__file__) != chk_sha: raise SystemExit("checker copy changed under the sweep")
    if sha256(cert) != cert_sha: raise SystemExit("certificate changed under the sweep")
    m, pts, ws = Z.read_cert(cert)
    chk = Z.Checker(m, pts, ws, max_depth=depth, use_chain=True, chain_from=0, theta_bias=theta_bias)
    chk.fast = True
    if not chk.symmetric_d4(): raise SystemExit("certificate is not D4-invariant")
    G.update(chk=chk, cert_sha=cert_sha, chk_sha=chk_sha, depth=depth)

def wroot(root_s):
    chk = G['chk']
    root = tuple(F(v) for v in root_s.split(','))
    t0 = time.time(); c0 = time.process_time()
    for k in chk.nstat: chk.nstat[k] = 0
    st, unc, _ = chk.run_box(root)
    return dict(root=root_s, stats=st, unc=[[str(v) for v in b] for b in unc],
                sec=round(time.time() - t0, 2), cpu=round(time.process_time() - c0, 2),
                checker_sha=G['chk_sha'], cert_sha=G['cert_sha'], depth=G['depth'],
                nstat=dict(chk.nstat), host=os.uname()[1], pid=os.getpid(),
                finished=time.strftime('%Y-%m-%d %H:%M:%S'))

# ------------------------------------------------------------------ helpers
def region(Z, m):
    return [rkey(r) for r in Z.d4_roots(m, F(SETTINGS['pitch']), SETTINGS['ubins'])]

def heaviness(rs, m):
    """sort key, heaviest first: u-bin 0/1 near a tile centre (k+1/2, l+1/2), then near a grid
    line, then the rest (a scheduling heuristic only)."""
    x0, x1, y0, y1, u0, u1 = (F(v) for v in rs.split(','))
    def dist(a0, a1):          # distance of the cell to the nearest half-integer line
        best = 9
        for k in range(int(m) + 1):
            h = F(2 * k + 1, 2)
            d = 0 if a0 <= h <= a1 else min(abs(a0 - h), abs(a1 - h))
            best = min(best, d)
        return best
    dx, dy = dist(x0, x1), dist(y0, y1)
    germ = (dx == 0 and dy == 0)
    line = (dx == 0 or dy == 0)
    return (0 if germ and u0 == 0 else 1 if germ else 2 if line and u0 == 0 else 3 if line else 4,
            float(u0), rs)

def read_records(out):
    p = os.path.join(out, 'roots.jsonl')
    recs = []
    if os.path.exists(p):
        with open(p) as f:
            for line in f:
                line = line.strip()
                if not line: continue
                try: recs.append(json.loads(line))
                except json.JSONDecodeError: pass      # a torn last line from a kill: redo that root
    return recs

def busy(pattern):
    try:
        r = subprocess.run(['pgrep', '-f', pattern], capture_output=True, text=True)
        me = {str(os.getpid()), str(os.getppid())}
        return any(p and p not in me for p in r.stdout.split())
    except Exception:
        return True

# ------------------------------------------------------------------ run
def cmd_run(a):
    out = a.out; os.makedirs(os.path.join(out, 'checker'), exist_ok=True)
    man_p = os.path.join(out, 'manifest.json')
    cert_sha = sha256(a.cert)
    settings = dict(SETTINGS, depth=a.depth)
    if not os.path.exists(man_p):
        src = os.path.join(HERE, 'zeromargin.py')
        shutil.copy2(src, os.path.join(out, 'checker', 'zeromargin.py'))
        man = dict(cert=os.path.abspath(a.cert), cert_sha=cert_sha,
                   checker_src=src, checker_sha=sha256(os.path.join(out, 'checker', 'zeromargin.py')),
                   runner_sha=sha256(os.path.abspath(__file__)), settings=settings,
                   created=time.strftime('%Y-%m-%d %H:%M:%S'))
        json.dump(man, open(man_p, 'w'), indent=1)
    man = json.load(open(man_p))
    chk_sha = sha256(os.path.join(out, 'checker', 'zeromargin.py'))
    if man['cert_sha'] != cert_sha: sys.exit("REFUSED: certificate sha256 differs from the manifest")
    if man['checker_sha'] != chk_sha: sys.exit("REFUSED: frozen checker copy differs from the manifest")
    if man['settings'] != settings: sys.exit(f"REFUSED: settings {settings} differ from the manifest {man['settings']}")
    Z = load_checker(out)
    m, pts, ws = Z.read_cert(a.cert)
    chk = Z.Checker(m, pts, ws)
    if not chk.symmetric_d4(): sys.exit("REFUSED: certificate is not D4-invariant")
    allr = region(Z, m)
    valid = [r for r in read_records(out)
             if r.get('checker_sha') == chk_sha and r.get('cert_sha') == cert_sha and r.get('depth', 0) >= a.depth]
    rdepth = a.depth
    if a.deepen:
        if a.deepen <= a.depth: sys.exit("--deepen must exceed the base depth")
        cert_ok = {r['root'] for r in valid if r['stats']['UNCERT'] == 0}
        tried = {r['root'] for r in valid if r['depth'] >= a.deepen}
        present = {r['root'] for r in valid}
        todo = sorted((r for r in allr if r in present and r not in cert_ok and r not in tried),
                      key=lambda rs: heaviness(rs, m))
        rdepth = a.deepen
    else:
        done = {r['root'] for r in valid if r['depth'] == a.depth}
        todo = sorted((r for r in allr if r not in done), key=lambda rs: heaviness(rs, m))
        done = done
    if a.only: todo = [r for r in todo if a.only in r]
    lo, hi = (int(v) for v in a.nproc_auto.split(':')) if a.nproc_auto else (a.nproc, a.nproc)
    done = {r['root'] for r in valid}
    print(f"{time.strftime('%H:%M:%S')} region {len(allr)} roots, with records {len(done)}, todo {len(todo)}; "
          f"checker {chk_sha[:12]} cert {cert_sha[:12]} depth {rdepth}; procs {lo}..{hi}", flush=True)
    resf = open(os.path.join(out, 'roots.jsonl'), 'a')
    ctx = mp.get_context('fork')
    t0 = time.time(); nd = 0; cpu = 0.0; unc_roots = 0
    with ProcessPoolExecutor(max_workers=hi, mp_context=ctx, initializer=winit,
                             initargs=(out, a.cert, cert_sha, chk_sha, rdepth, SETTINGS['theta_bias'])) as ex:
        pend = {}; it = iter(todo); exhausted = False; last_hb = 0; last_chk = 0; cap = lo
        while True:
            now = time.time()
            if now - last_chk > 60:
                cap = lo if (a.nproc_auto and busy(a.busy_pattern)) else hi
                last_chk = now
            while not exhausted and len(pend) < cap:
                try: rs = next(it)
                except StopIteration: exhausted = True; break
                pend[ex.submit(wroot, rs)] = rs
            if not pend: break
            dn, _ = wait(list(pend), timeout=30, return_when=FIRST_COMPLETED)
            for fu in dn:
                rs = pend.pop(fu)
                rec = fu.result()                     # a worker exception stops the sweep
                resf.write(json.dumps(rec) + "\n"); resf.flush(); os.fsync(resf.fileno())
                nd += 1; cpu += rec['cpu']
                s = rec['stats']
                if s['UNCERT']: unc_roots += 1
                if s['UNCERT'] or rec['cpu'] > 600:
                    print(f"{time.strftime('%H:%M:%S')} ROOT {rs} boxes {s['boxes']} depth {s['maxdepth']} "
                          f"uncert {s['UNCERT']} cpu {rec['cpu']:.0f}s", flush=True)
            if time.time() - last_hb > 300:
                last_hb = time.time(); el = time.time() - t0
                print(f"{time.strftime('%H:%M:%S')} HB done {nd}/{len(todo)} this run ({len(done)+nd}/{len(allr)} total), "
                      f"cpu {cpu/3600:.2f} h, wall {el/3600:.2f} h, in flight {len(pend)} (cap {cap}), "
                      f"roots with uncertified boxes {unc_roots}", flush=True)
    print(f"{time.strftime('%H:%M:%S')} run finished: {nd} roots, {cpu/3600:.2f} CPU-h, {(time.time()-t0)/3600:.2f} h wall; "
          f"roots with uncertified boxes {unc_roots}", flush=True)

# ------------------------------------------------------------------ summary
def cmd_summary(a):
    out = a.out
    man = json.load(open(os.path.join(out, 'manifest.json')))
    chk_sha = sha256(os.path.join(out, 'checker', 'zeromargin.py'))
    Z = load_checker(out)
    m, pts, ws = Z.read_cert(man['cert'])
    cert_sha_now = sha256(man['cert'])
    chk = Z.Checker(m, pts, ws)
    d4ok = chk.symmetric_d4()
    allr = region(Z, m)
    recs = read_records(out)
    by = {}; bad = []; nrec = {}; cpu_all = 0.0
    for r in recs:
        if (r.get('checker_sha') != man['checker_sha'] or r.get('cert_sha') != man['cert_sha']
                or r.get('depth', 0) < man['settings']['depth']):
            bad.append(r['root']); continue
        nrec[r['root']] = nrec.get(r['root'], 0) + 1; cpu_all += r['cpu']
        o = by.get(r['root'])
        if o is not None and o['depth'] == r['depth'] and o['stats'] != r['stats']:
            bad.append(r['root'] + ' (conflicting duplicate)')
        # the record shown: a certified one if any (smallest depth), else the deepest attempt
        if (o is None or (r['stats']['UNCERT'] == 0 and (o['stats']['UNCERT'] > 0 or r['depth'] < o['depth']))
                or (o['stats']['UNCERT'] > 0 and r['stats']['UNCERT'] > 0 and r['depth'] > o['depth'])):
            by[r['root']] = r
    lines = []; tot = {}; cpu = cpu_all; missing = []; uncl = []; maxd = 0
    for rs in allr:
        r = by.get(rs)
        if r is None:
            missing.append(rs); lines.append(f"{rs:40s} MISSING"); continue
        s = r['stats']; maxd = max(maxd, s['maxdepth'])
        for k in ('ADM', 'P1', 'MIX', 'CHAIN', 'EMPTY', 'UNCERT', 'boxes'): tot[k] = tot.get(k, 0) + s[k]
        v = 'CERTIFIED' if s['UNCERT'] == 0 else f"UNCERTIFIED({s['UNCERT']})"
        if s['UNCERT']: uncl.append(r)
        lines.append(f"{rs:40s} {v:18s} limit {r['depth']} records {nrec[rs]} boxes {s['boxes']:6d} depth {s['maxdepth']:2d} ADM {s['ADM']} P1 {s['P1']} "
                     f"MIX {s['MIX']} CHAIN {s['CHAIN']} EMPTY {s['EMPTY']} cpu {r['cpu']:.1f}s "
                     f"checker {r['checker_sha'][:16]} cert {r['cert_sha'][:16]}")
    extra = set(by) - set(allr)
    clean = (not missing and not uncl and not bad and not extra and d4ok
             and cert_sha_now == man['cert_sha'] and chk_sha == man['checker_sha'])
    head = [
        f"zeromargin.py D4 re-check: {out}",
        f"certificate {man['cert']}  sha256 {man['cert_sha']}  (now {cert_sha_now})",
        f"checker     {out}/checker/zeromargin.py  sha256 {man['checker_sha']}  (now {chk_sha})",
        f"runner sha256 {man['runner_sha']} (now {sha256(os.path.abspath(__file__))}); settings {json.dumps(man['settings'])}",
        f"roots certified only by a --deepen record: {sum(1 for r in by.values() if r['depth'] > man['settings']['depth'] and r['stats']['UNCERT'] == 0)}",
        f"D4 invariance of the weighted point set (exact, x->m-x and x<->y): {d4ok}",
        f"region: {len(allr)} roots = [0,{m/2}]^2 x u in [0,1/2]; present {len(allr)-len(missing)}, "
        f"missing {len(missing)}, with uncertified boxes {len(uncl)}, bad/foreign records {len(bad)}, extra {len(extra)}",
        f"totals: boxes {tot.get('boxes',0)}, max depth {maxd}, leaves ADM {tot.get('ADM',0)} P1 {tot.get('P1',0)} "
        f"MIX {tot.get('MIX',0)} CHAIN {tot.get('CHAIN',0)} EMPTY {tot.get('EMPTY',0)} UNCERTIFIED {tot.get('UNCERT',0)}; "
        f"CPU {cpu/3600:.2f} h (all valid records); census of the shown records",
    ]
    verdict = ("D4 RECHECK CLEAN: every root of the D4 region certified" if clean else
               "INCOMPLETE" if (missing and not uncl) else "NOT CLEAN")
    body = head + [""] + [f"uncertified: {r['root']}: " + "; ".join(",".join(b) for b in r['unc'][:50]) for r in uncl] \
        + [f"bad record: {b}" for b in bad[:50]] + [""] + lines + ["", verdict]
    with open(os.path.join(out, 'SUMMARY.txt'), 'w') as f: f.write("\n".join(body) + "\n")
    print("\n".join(head));
    for r in uncl: print(f"  UNCERTIFIED root {r['root']}: {r['stats']['UNCERT']} boxes")
    print(verdict)
    return 0 if clean else 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['run', 'summary'])
    ap.add_argument('--cert', default='runs/s32-close_candidate.txt')
    ap.add_argument('--out', default='runs/s32py_d4')
    ap.add_argument('--depth', type=int, default=24)
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--deepen', type=int, default=None, help='rerun the uncertified roots at this depth')
    ap.add_argument('--nproc-auto', default=None, help='A:B -- A while --busy-pattern runs, else B')
    ap.add_argument('--busy-pattern', default=r's32d4_v1_sweep|out runs/s32d4_v1( |$)')
    ap.add_argument('--only', default=None, help='substring filter on root keys (testing only)')
    a = ap.parse_args()
    if a.cmd == 'run':
        if a.depth < 24: sys.exit("depth must be >= 24")
        cmd_run(a)
    else:
        sys.exit(cmd_summary(a))

if __name__ == '__main__':
    main()
