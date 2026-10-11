"""Allowlist scan of the generated V7 root modules (our code; reads text only).

Usage: python3 scan_v7.py OUTDIR [--survey]
Checks every .lean file:
  * forbidden tokens anywhere (axiom, sorry, native_decide, macro, elab, #eval, set_option other than the two
    benign ones, attributes, unsafe, extern, implemented_by, ofReduceBool, debug options, ...);
  * structural lines only from a fixed set (imports of Sqpack modules, the fixed header, the one namespace);
  * every declaration is `theorem`/`def`/`abbrev` with a plain generated name (no dots, no guillemets), matching
    root_NNNNN or rNNNNN<suffix>, so nothing can declare or shadow a name like Bentz.Valid7;
  * file names match the module names they are imported as.
"""
import os, re, sys, collections

FORBID = re.compile(r'\b(axiom|axioms|sorry|sorryAx|native_decide|ofReduceBool|ofReduceNat|implemented_by|extern|unsafe|'
                    r'opaque|partial|macro|macro_rules|elab|elab_rules|syntax|notation|infix|infixl|infixr|prefix|postfix|'
                    r'initialize|builtin_initialize|run_cmd|run_tac|run_elab|instance|attribute|export|'
                    r'Lean|IO|debug|skipKernelTC|csimp|reducible|irreducible|local|scoped|universe|variable|'
                    r'mutual|termination_by|decreasing_by|example)\b|#[a-z_]+|@\[|«|»|`\(|\$\(')
SETOPT_OK = {'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0'}
HEADER_OK = {'noncomputable section  -- checked by the kernel only: no compiled code',
             'namespace SquarePacking.LemmaELeaf', 'end SquarePacking.LemmaELeaf',
             'open ZMTreeM LemmaEVert LemmaEPoly ZMTree LemmaE'}
IMPORT = re.compile(r'^import (Sqpack(?:\.[A-Za-z][A-Za-z0-9_]*)+)$')
DECL = re.compile(r'^(theorem|def|abbrev|lemma|noncomputable def|private def|private theorem) ([^\s:(]+)')
NAME_OK = re.compile(r'^(root_\d{5}|r\d{5}[A-Za-z0-9_]*)$')
KW_START = re.compile(r'^([A-Za-z#@_][A-Za-z0-9_]*)')

def main():
    out = sys.argv[1]; survey = '--survey' in sys.argv
    starts = collections.Counter(); kinds = collections.Counter(); imports_ext = collections.Counter()
    bad = []; nfiles = 0; nbytes = 0; decls = 0; names = set()
    for dp, _, fs in os.walk(out):
        for f in fs:
            if not f.endswith('.lean'):
                bad.append((os.path.join(dp, f), 'non-lean file')); continue
            p = os.path.join(dp, f); rel = os.path.relpath(p, out)
            if not re.fullmatch(r'(R\d{5}\.lean|R\d{5}/r\d{5}[A-Za-z0-9_]*\.lean)', rel):
                bad.append((rel, 'unexpected file name'))
            txt = open(p, encoding='utf-8').read(); nfiles += 1; nbytes += len(txt)
            for i, line in enumerate(txt.split('\n'), 1):
                s = line.rstrip()
                m = KW_START.match(s)
                if m: starts[m.group(1)] += 1
                code = s.split('--', 1)[0] if not s.startswith('noncomputable section') else ''
                if s.startswith('set_option'):
                    if s not in SETOPT_OK: bad.append((rel, i, 'set_option', s[:120]))
                    continue
                if s.startswith('import '):
                    mi = IMPORT.match(s)
                    if not mi: bad.append((rel, i, 'import', s[:120])); continue
                    mod = mi.group(1)
                    if not mod.startswith('Sqpack.V7.R'): imports_ext[mod] += 1
                    continue
                if s in HEADER_OK or s == '': continue
                if s.startswith(('namespace', 'section', 'end ', 'open ')):
                    bad.append((rel, i, 'structure', s[:120])); continue
                fm = FORBID.search(code)
                if fm: bad.append((rel, i, 'forbidden ' + fm.group(0), s[:120]))
                dm = DECL.match(s)
                if dm:
                    decls += 1; kinds[dm.group(1)] += 1; nm = dm.group(2)
                    if not NAME_OK.match(nm): bad.append((rel, i, 'decl name', nm))
                    if nm in names: bad.append((rel, i, 'duplicate decl', nm))
                    names.add(nm)
                elif m and not line.startswith((' ', '\t')) and m.group(1) not in ('theorem', 'def', 'abbrev'):
                    bad.append((rel, i, 'unexpected top-level', s[:120]))
    roots = sorted(n for n in names if n.startswith('root_'))
    print(f'{nfiles} files, {nbytes/1e9:.2f} GB, {decls} declarations ({dict(kinds)}), {len(roots)} root theorems')
    print('imports outside V7.R*:', dict(imports_ext))
    if survey:
        print('line-start tokens:', starts.most_common(40))
    missing = [i for i in range(9800) if f'root_{i:05d}' not in names]
    print('missing root theorems:', len(missing), missing[:10])
    print('problems:', len(bad))
    for b in bad[:40]: print('  ', b)
    sys.exit(1 if bad or missing else 0)

if __name__ == '__main__':
    main()
