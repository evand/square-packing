#!/usr/bin/env python3
"""Record candidate -> registration packet (10-12; written after s(154), jlevy/squares#552, was assembled by hand;
10-14: descent, prior record and distinctness built in after s(146), jlevy/squares#606, was assembled by hand).

  candidate.py FILE --name NAME [--vs FILE ...] [--chain RUNDIR] [--no-regularize] [--no-descend] [--no-prior] [--stage]
               [--import --ref TEXT]
  candidate.py --live N [S] [--no-sync]   just the live check (10-17: the store, synced first when the probe
                                       finds newer issues / commits; --no-sync: read-only, stale items listed)

Steps (work dir runs/candidates/NAME/, summary.json + a facts block for the issue draft):
  1. exactsolve (../exact/env.sh: the pinned batch environment) -> KKT point, rational certificate.  If it says "not
     jammed" (a first-order descent exists once corner-corner touches may slip; 2 of 9 TT2 successes, one of them the
     s(146) record), ../exact/descend.py (MILP kick + slp2) and exactsolve again; the original input is kept as
     input.orig.txt, summary['descent'] has both sides.  --no-descend: fail instead
  2. both exact checkers (verify_cert.py separating axes, verify_cert2.py polygons) on the certificate: both must say VALID
  3. live check (10-17): the store, after pk.watch.probe checks it is current (register main sha, issues / PRs updated
     since the last watch, watched repos' HEADs); anything newer -> `pk.py sync-register` + `pk.py watch --apply`
     first (--no-sync: don't; stale items are listed).  Sources: every frontier row at n (register, each issue's
     claims), the best stored packing (anyone's, and by others), the monotone envelope.  "beaten" = any of them at or
     below S' + 1e-12
  4. matching distance (pk.crossing) to the store's best packing for n and to each --vs file
  5. local-optimality facts from exactsolve's report (min lambda, flat motions, corner-corner jams, rattlers)
  6. --chain RUNDIR: lineage through a chain.py run (chain_lineage.py), for the credit paragraph
  7. regularized view (search/regularize, coincide style, conservative policy) at the same S*, both checkers on its
     certificate.  The input goes through a private REFRESH dir in the work dir, so runs/regularize/exact (the
     shared refresh inputs) is untouched
  8. prior record (--no-prior skips): the best public packing at n in the store by someone else
     (pk.distinct.prior_packing), through the same steps 1, 2, 5, 7 in prior/; then pk.distinct.compare on the two
     regularized views -> vs_prior.json (sides, optimality, contacts, tilt angles of load-bearing squares, matched
     square-by-square diff, verdict 'distinct minimum' / 'same minimum') and vs_prior.md (the issue's comparison table)
  9. --stage: copy to candidates/NAME{.cert,.exact.txt,.input.txt,.json} (+ NAME_regularized.{cert,json},
     NAME.predescent.txt, NAME_vs_prior.json, NAME_prior.json) for a commit to pin.  --import: add to the store
     (source candidate:NAME, finder ours, status certified, --ref e.g. "jlevy/squares#552")
Nothing is posted: the issue text is drafted in outreach/ and posted per Evan's OK (CLAUDE.md).
"""
import argparse, json, os, re, shutil, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
EX = os.path.join(HERE, '..', 'exact')
RG = os.path.join(HERE, '..', 'regularize')


def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def checkers(cert):
    out = {}
    for v in ('verify_cert', 'verify_cert2'):
        r = sh([sys.executable, os.path.join(EX, v + '.py'), cert])
        last = (r.stdout.strip().splitlines() or [''])[-1].strip()
        out[v] = dict(ok=r.returncode == 0 and 'VALID' in last and 'NOT' not in last and 'INVALID' not in last, line=last)
    return out


# ------------------------------------------------------------------ live check (read-only)

def live(n, S, sync=False):
    """10-17: the record check reads the store, after a cheap probe that the store is current (pk.watch.probe:
    register main sha, issues / PRs updated since the last watch, watched repos' HEADs).  sync=True: when the probe
    finds anything newer, run `pk.py sync-register` and `pk.py watch --apply` first, then probe again.  What is still
    stale is listed (L['stale']) and printed.  Sources: every frontier row at n (register + each issue's claims, with
    or without a file), the best stored packing at n (anyone's) and by someone other than us, and the monotone envelope.
    (Until 10-17 this read the register, case page, issue texts (decimal grep) and Ellsworth directly; the grep gave
    n = 155 the s(154) value from #401 and blocked a real record, and caught nothing the store missed.)"""
    from pk import watch
    from pk.store import Store
    L = dict(n=n, S=S, sources={}, errors=[], stale=[], as_of={})
    try:
        P = watch.probe()
        if P['stale'] and sync:
            for cmd in (['sync-register'], ['watch', '--apply']):
                r = sh([sys.executable, os.path.join(HERE, 'pk.py'), *cmd], cwd=HERE, timeout=3600)
                if r.returncode:
                    L['errors'].append(f'pk.py {cmd[0]}: rc {r.returncode}: {(r.stderr or r.stdout)[-300:]}')
            P = watch.probe()
        L['stale'], L['as_of'] = P['stale'], P['as_of']
    except Exception as e:
        L['errors'].append(f'probe: {e!r}'[:300])
    try:
        st = Store()
        for r in st.db.execute('SELECT * FROM frontier WHERE n=?', (n,)):
            L['sources'][f"store frontier: {r['source']} ({r['finder'] or '?'})"] = r['s']
        rows = st.query(n=n, limit=None)
        if rows:
            L['sources'][f"store best any ({(rows[0]['finders'] or '?')[:40]})"] = rows[0]['s']
            oth = next((r for r in rows if any(f.strip() != 'ours' for f in (r['finders'] or 'ours').split(','))), None)
            if oth is not None:
                L['sources'][f"store best public ({(oth['finders'] or '?')[:40]})"] = oth['s']
        # monotone envelope (10-12): s(n) <= s(m) for m > n by dropping squares
        r = st.db.execute("SELECT n, s FROM packing WHERE n > ? AND n <= ? AND status != 'infeasible' AND pen <= 1e-9 "
                          "ORDER BY s LIMIT 1", (n, n + 80)).fetchone()
        if r and r['s'] <= S + 0.05:
            L['sources'][f'store s({r["n"]}) (envelope)'] = r['s']
    except Exception as e:
        L['errors'].append(f'store: {e!r}'[:300])
    L['register_main'] = L['as_of'].get('register')
    L['latest_issue'] = L['as_of'].get('last_issue')
    L['best_other'] = min(L['sources'].values()) if L['sources'] else None
    L['beaten'] = [k for k, v in L['sources'].items() if v <= S + 1e-12]
    return L


def print_live(L):
    print(f'live check n = {L["n"]} vs {L["S"]:.13f}  (store as of: register {L.get("register_main")}, issues updated '
          f'through {L.get("as_of", {}).get("issues_updated")}, last issue #{L.get("latest_issue")})')
    for k, v in sorted(L['sources'].items(), key=lambda kv: kv[1]):
        print(f'  {v:.13f}  {v - L["S"]:+.2e}  {k}' + ('   <-- at or below' if k in L['beaten'] else ''))
    for e in L['errors']:
        print(f'  ERROR {e}')
    for s_ in L.get('stale', []):
        print(f'  STALE (store may be behind): {s_}')
    print('  BEATEN / TIED: ' + ', '.join(L['beaten']) if L['beaten'] else '  nothing at or below')


# ------------------------------------------------------------------ regularize (private REFRESH dir)

def regularize(work, n):
    d = os.path.join(work, 'reg'); os.makedirs(os.path.join(d, 'exact', f'n-{n}'), exist_ok=True)
    shutil.copy(os.path.join(work, 'cand.exact.txt'), os.path.join(d, 'exact', f'n-{n}', 'input.exact.txt'))
    shutil.copy(os.path.join(work, 'cand.json'), os.path.join(d, 'exact', f'n-{n}', 'input.json'))
    code = ('import sys, os, json; sys.path.insert(0, sys.argv[1]); import reg, atlas; d = sys.argv[2]; n = int(sys.argv[3])\n'
            "reg.REFRESH = os.path.join(d, 'exact'); atlas.RUNS = os.path.join(d, 'out'); atlas.COORDS = os.path.join(d, 'out', 'coords')\n"
            "assert reg.refreshed(n), 'candidate input not picked up'\n"
            "print(json.dumps(atlas.compute_one((n, 'coincide', 'conservative'))))\n")
    r = sh([sys.executable, '-c', code, os.path.abspath(RG), d, str(n)], cwd=d, timeout=7200)
    base = os.path.join(d, 'out', 'poster', 'conservative', f'n-{n:03d}')
    if r.returncode or not os.path.exists(base + '.json'):
        return dict(ok=False, error=(r.stderr or r.stdout)[-500:])
    j = json.load(open(base + '.json'))
    shutil.copy(base + '.cert', os.path.join(work, 'regularized.cert')); shutil.copy(base + '.json', os.path.join(work, 'regularized.json'))
    return dict(ok=bool(j['cert_valid'] and j['verify_80']['ok']), S_exact=j['S_exact'], S_cert=j['S_cert'],
                groups=[j['groups_before'], j['groups_after']], level=j['level'], level_text=j['level_text'],
                checkers=checkers(os.path.join(work, 'regularized.cert')))


# ------------------------------------------------------------------ exactsolve (+ corner-slip descent)

def exactsolve(work):
    """exactsolve on work/cand.txt -> (cand.json as a dict or None, the subprocess result)"""
    r = sh([os.path.join(EX, 'env.sh'), 'python3', os.path.join(EX, 'exactsolve.py'), os.path.join(work, 'cand.txt'),
            '--out', work, '-q'], cwd=EX, timeout=7200)
    f = os.path.join(work, 'cand.json')
    return (json.load(open(f)) if os.path.exists(f) else None), r


def not_jammed(J):
    return bool(J) and not J.get('cert_valid') and (J.get('status') == 'not jammed' or
                                                    (J.get('vv_milp_dS') or 0) < -J.get('jam_threshold', 1e-7))


def solve(work, src, rounds=30, descend=True, log=print):
    """src -> work/cand.txt -> exactsolve; if exactsolve says 'not jammed' (a first-order descent once corner-corner
    touches may slip), exact/descend.py (MILP kick + slp2, as for s(146), jlevy#606) and exactsolve again.  The original
    input is kept as work/input.orig.txt.  Returns (J, info) with info['descent'] when it ran; J None = no certificate."""
    from pk.packing import read
    os.makedirs(work, exist_ok=True)
    for f in ('cand.json', 'cand.cert', 'cand.exact.txt', 'cand.contacts.json', 'input.orig.txt', 'descended.txt',
              'descend.log', 'regularized.cert', 'regularized.json', 'vs_prior.json', 'vs_prior.md'):
        if os.path.exists(os.path.join(work, f)):                 # a rerun must not stage an earlier run's files
            os.remove(os.path.join(work, f))
    if os.path.isdir(os.path.join(work, 'prior')):
        shutil.rmtree(os.path.join(work, 'prior'))
    dst = os.path.join(work, 'cand.txt')
    if not (os.path.exists(dst) and os.path.samefile(src, dst)):
        shutil.copy(src, dst)
    J, r = exactsolve(work)
    info = {}
    if descend and not_jammed(J):
        shutil.copy(src, os.path.join(work, 'input.orig.txt'))
        d0 = dict(s_before=float(J['S_input']), vv_milp_dS=J.get('vv_milp_dS'), why=J.get('why'))
        log(f'    not jammed (corner-corner MILP dS* = {J.get("vv_milp_dS"):.3e}): exact/descend.py, then exactsolve again',
            flush=True)
        out = os.path.join(work, 'descended.txt')
        q = sh([os.path.join(EX, 'env.sh'), 'python3', os.path.join(EX, 'descend.py'), os.path.join(work, 'input.orig.txt'),
                out, '--rounds', str(rounds)], cwd=EX, timeout=7200)
        open(os.path.join(work, 'descend.log'), 'w').write(q.stdout + q.stderr)
        if q.returncode or not os.path.exists(out):
            info['descent'] = dict(d0, ok=False, error=(q.stderr or q.stdout)[-300:])
            return None, info
        for f in ('cand.json', 'cand.cert', 'cand.exact.txt', 'cand.contacts.json'):
            if os.path.exists(os.path.join(work, f)):
                os.remove(os.path.join(work, f))
        shutil.copy(out, os.path.join(work, 'cand.txt'))
        J, r = exactsolve(work)
        d0.update(ok=True, s_after=read(out).s, rounds=sum(l.startswith('round ') for l in q.stdout.splitlines()))
        info['descent'] = d0
        log(f'    descended {d0["s_before"]:.13f} -> {d0["s_after"]:.13f}', flush=True)
    if not J or not J.get('cert_valid'):
        info['exactsolve_tail'] = (r.stdout[-2000:], r.stderr[-2000:])
        info['why'] = (J or {}).get('why') or (J or {}).get('status')
        return None, info
    return J, info


def optimality(J):
    so = J.get('second_order', {})
    return dict(kkt=J.get('newton_converged'), min_lambda=J.get('lambda_A_maxmin'), status=so.get('status'),
                flat_motions=so.get('n_zero'), flat_realised=so.get('zero_modes_realised'),
                flat_residual=(so.get('zero_mode_probe') or {}).get('worst_residual'),
                negative_eigs=so.get('n_neg'), corner_corner=J.get('n_vv_exact', J.get('n_vv')),
                vv_milp_dS=J.get('vv_milp_dS'), rattlers=len(J.get('free_squares') or []),
                min_pair=J.get('cert_min_pair'), min_wall=J.get('cert_min_wall'))


def prior_record(work, n, do_reg=True, log=print):
    """The record this candidate beats (pk.distinct.prior_packing: the best public store packing by someone else), through
    the same solve / checkers / regularize path in work/prior/; then pk.distinct.compare -> work/vs_prior.{json,md}"""
    from pk import distinct as D
    got = D.prior_packing(n)
    if not got:
        log('    no public packing at n in the store'); return dict(ok=False, error='no public packing')
    pid, Q, src = got
    pw = os.path.join(work, 'prior'); os.makedirs(pw, exist_ok=True)
    D.write_float(Q, os.path.join(pw, 'input.txt'))
    log(f'    prior: store {pid}, s = {Q.s:.13f}; ' + '; '.join(src[:4]) + (' ...' if len(src) > 4 else ''), flush=True)
    J, info = solve(pw, os.path.join(pw, 'input.txt'), log=log)
    info.pop('exactsolve_tail', None)
    res = dict(store_id=pid, s=Q.s, sources=src, **info)
    if not J:
        log(f'    prior: no certificate ({info.get("why")})'); res['ok'] = False
        json.dump(res, open(os.path.join(pw, 'summary.json'), 'w'), indent=1, default=str)
        return res
    res.update(ok=True, S_exact=J['S_exact'], S_cert=J['S_cert_decimal'], checkers=checkers(os.path.join(pw, 'cand.cert')),
               optimality=optimality(J))
    log(f'    prior S* = {J["S_exact"][:30]}; ' + ', '.join(f'{k} {"VALID" if v["ok"] else "FAILED"}'
                                                       for k, v in res['checkers'].items()), flush=True)
    if do_reg:
        res['regularized'] = R = regularize(pw, n)
        if not R.get('ok'):
            log(f'    prior: regularize failed ({R.get("error", "")[-160:]})')
        elif os.path.exists(os.path.join(work, 'regularized.json')):
            C = D.compare(work, pw)
            json.dump(C, open(os.path.join(work, 'vs_prior.json'), 'w'), indent=1, default=str)
            open(os.path.join(work, 'vs_prior.md'), 'w').write(D.table_md(C, 'this packing', f'prior (store {pid})'))
            for l in D.summary_lines(C):
                log('    ' + l)
            res['verdict'] = C['verdict']
            res['difference'] = {k: v for k, v in C['difference'].items() if k != 'moved_squares'}
    json.dump(res, open(os.path.join(pw, 'summary.json'), 'w'), indent=1, default=str)
    return res


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('file', nargs='?'); ap.add_argument('--name')
    ap.add_argument('--vs', nargs='*', default=[]); ap.add_argument('--chain')
    ap.add_argument('--no-regularize', action='store_true'); ap.add_argument('--stage', action='store_true')
    ap.add_argument('--no-descend', action='store_true', help='no exact/descend.py when exactsolve says not jammed')
    ap.add_argument('--no-prior', action='store_true', help='skip certifying / comparing the record this candidate beats')
    ap.add_argument('--import', dest='imp', action='store_true'); ap.add_argument('--ref')
    ap.add_argument('--live', nargs='+', metavar=('N', 'S'))
    ap.add_argument('--no-sync', action='store_true', help='live check: do not sync a stale store (read-only)')
    a = ap.parse_args()
    if a.live:
        n = int(a.live[0])
        from pk.store import Store
        S = float(a.live[1]) if len(a.live) > 1 else Store().best(n).s
        print_live(live(n, S, sync=not a.no_sync)); return
    if not (a.file and a.name):
        ap.error('FILE and --name are required')
    from pk.packing import read
    from pk.store import Store
    from pk.crossing import match_distance
    P = read(a.file); n = P.n
    work = os.path.join(HERE, 'runs', 'candidates', a.name); os.makedirs(work, exist_ok=True)
    from pk import manifest
    man = manifest.start(work, 'candidate', args=dict(vars(a), n=n), procs=1)     # measured CPU for pk.effort
    out = dict(name=a.name, n=n, input=os.path.abspath(a.file), s_input=P.s)

    print(f'[1] exactsolve n = {n}, s = {P.s:.12f}', flush=True)
    J, info = solve(work, a.file, descend=not a.no_descend)
    out.update({k: v for k, v in info.items() if k != 'exactsolve_tail'})
    if not J:
        if 'exactsolve_tail' in info:
            print(*info['exactsolve_tail'])
        json.dump(out, open(os.path.join(work, 'summary.json'), 'w'), indent=1, default=str)
        sys.exit(f'exactsolve: no valid certificate ({info.get("why")})')
    out['S_exact'], out['S_cert'] = J['S_exact'], J['S_cert_decimal']
    S = float(J['S_cert_decimal'])
    Pc = read(os.path.join(work, 'cand.txt'))                 # the certified f64 input (after the descent, if any)
    print(f'    S* = {J["S_exact"][:30]}   S\' = {J["S_cert_decimal"]}', flush=True)

    print('[2] exact checkers', flush=True)
    out['checkers'] = checkers(os.path.join(work, 'cand.cert'))
    for k, v in out['checkers'].items():
        print(f'    {k}: {"VALID" if v["ok"] else "FAILED"}  ({v["line"][:90]})')
    if not all(v['ok'] for v in out['checkers'].values()):
        sys.exit('a checker rejected the certificate')

    print('[3] live check', flush=True)
    out['live'] = live(n, S, sync=not a.no_sync); print_live(out['live'])

    print('[4] matching distance', flush=True)
    st = Store(); out['distance'] = {}
    b = st.best(n)
    if b:
        out['distance'][f'store best {b.s:.12f}'] = match_distance(Pc, b)
    for f in a.vs:
        q = read(f); out['distance'][f'{os.path.basename(f)} {q.s:.12f}'] = match_distance(Pc, q)
    for k, v in out['distance'].items():
        print(f'    {v:.3f}  vs {k}')

    out['optimality'] = optimality(J)
    print('[5] local optimality: ' + ', '.join(f'{k} {v}' for k, v in out['optimality'].items()))

    if a.chain:
        import chain_lineage
        out['lineage'] = chain_lineage.trace(a.chain, n, P.s)
        moved = [s for s in out['lineage'] if s.get('moves') or s.get('root', '').startswith('seed_k')]
        print(f'[6] lineage: {len(moved)} productive rounds, {sum(len(s.get("moves", [])) for s in moved)} moves; '
              f'origin {out["lineage"][-1]}')

    if not a.no_regularize:
        print('[7] regularize', flush=True)
        out['regularized'] = R = regularize(work, n)
        print(f'    {json.dumps({k: v for k, v in R.items() if k != "checkers"})}')
        for k, v in R.get('checkers', {}).items():
            print(f'    {k}: {"VALID" if v["ok"] else "FAILED"}')

    if not a.no_prior:
        print('[8] prior record: certified and regularized the same way; is ours a different minimum?', flush=True)
        out['prior'] = prior_record(work, n, do_reg=not a.no_regularize)

    if a.stage:
        dst = os.path.join(HERE, 'candidates', a.name)
        for src, ext in (('cand.cert', '.cert'), ('cand.exact.txt', '.exact.txt'), ('cand.txt', '.input.txt'), ('cand.json', '.json'),
                         ('regularized.cert', '_regularized.cert'), ('regularized.json', '_regularized.json'),
                         ('input.orig.txt', '.predescent.txt'), ('vs_prior.json', '_vs_prior.json'),
                         ('prior/summary.json', '_prior.json')):
            if os.path.exists(os.path.join(work, src)):
                shutil.copy(os.path.join(work, src), dst + ext)
        print(f'[9] staged candidates/{a.name}.*')
    json.dump(out, open(os.path.join(work, 'summary.json'), 'w'), indent=1, default=str)
    man.finish(summary=dict(n=n, S_cert=out.get('S_cert'), beaten=bool(out['live']['beaten'])))
    print(f'summary: {os.path.relpath(os.path.join(work, "summary.json"), HERE)}')
    if a.imp:
        if out['live']['beaten']:
            sys.exit('not importing as a record candidate: beaten / tied (see the live check)')
        r = sh([sys.executable, os.path.join(HERE, 'pk.py'), 'import', os.path.join(work, 'cand.exact.txt'), '--source', f'candidate:{a.name}',
                '--finder', 'ours', '--status', 'certified', '--kind', 'record-candidate'] + (['--ref', a.ref] if a.ref else []), cwd=HERE)
        print('[10] store import: ' + (r.stdout.strip() or r.stderr.strip())[-300:])


if __name__ == '__main__':
    main()
