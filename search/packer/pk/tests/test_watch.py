"""Ingestion safety + parsing tests (no network): python3 pk/tests/test_watch.py"""
import io, os, sys, tarfile, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pk.packing import parse, Packing
from pk import github, issues, watch

SQCERT = """# sqcert v1
# delta: 0
n 2
s 2
# x y t
1/2 1/2 0
3/2 3/2 0
"""


def raises(f, *a, **k):
    try:
        f(*a, **k)
    except Exception:
        return True
    return False


def test_sqcert_sqpk():
    p = parse(SQCERT, strict=True)
    assert p.n == 2 and p.s == 2 and p.s_str == '2' and p.max_pen() == 0
    q = parse(SQCERT.replace('# sqcert', '# sqpk').replace('3/2 3/2 0', '1.5 1.5 0.0'), strict=True)
    assert q.canon_hash() == p.canon_hash()
    t = parse(SQCERT.replace('3/2 3/2 0', '3/2 3/2 1/3'), strict=True)          # t = tan(theta/2)
    assert abs(t.sq[1, 2] - 36.86989764584402) < 1e-9
    for bad in (SQCERT.replace('n 2', 'n 3'),                     # row count != n
                SQCERT + '1/2 3/2 0\n',                           # extra row
                SQCERT.replace('3/2 3/2 0', '3/2 3/2'),           # short row
                SQCERT.replace('3/2 3/2 0', '1.5 3/2 0'),         # decimal in an exact file
                SQCERT.replace('3/2 3/2 0', '3/0 3/2 0'),         # zero denominator
                SQCERT.replace('3/2 3/2 0', 'nan 3/2 0'),
                SQCERT.replace('3/2 3/2 0', '1' * 5000 + ' 3/2 0'),
                SQCERT.replace('s 2', 's 1')):                    # side below sqrt(n)
        assert raises(parse, bad, strict=True), bad


def test_strict_generic():
    assert raises(parse, '2 2\n0.5 0.5 0\n1.5 1.5 0\n9 9 9\n', strict=True)     # trailing junk
    assert raises(parse, '2 2\n0.5 0.5 0\n1e400 1.5 0\n', strict=True)          # inf
    assert raises(parse, '2 2\n0.5 0.5 0\n50 1.5 0\n', strict=True)             # centre far outside
    assert parse('2 2\n0.5 0.5 0\n1.5 1.5 0\n', strict=True).n == 2
    assert raises(parse, '{"s": 1e999, "squares": [[0.5, 0.5, 0]]}', strict=True)


def test_refs():
    body = ('see https://github.com/alice/sq-certs/tree/0123abcd/certs/n100 and github.com/bob/x.git, '
            'https://github.com/jlevy/squares/issues/5, https://github.com/user-attachments/files/123/a.cert.txt, '
            'https://evil.example/github.com/x/y, commit `deadbeef12` and [r](https://github.com/alice/sq-certs).')
    R, loose = issues.refs(body)
    assert list(R) == ['alice/sq-certs', 'bob/x'], R
    assert R['alice/sq-certs'] == dict(shas=['0123abcd'], links=[('0123abcd', 'certs/n100')])
    assert loose == ['deadbeef12']
    assert [m.group(0) for m in issues.ATTACH_RE.finditer(body)] == ['https://github.com/user-attachments/files/123/a.cert.txt']
    R, _ = issues.refs('https://github.com/a/b/blob/main/../../etc/passwd')
    assert R['a/b']['links'] == []


def _tgz(members):
    b = io.BytesIO()
    with tarfile.open(fileobj=b, mode='w:gz') as t:
        for name, data, kind in members:
            ti = tarfile.TarInfo(name)
            if kind == 'sym':
                ti.type = tarfile.SYMTYPE; ti.linkname = data; t.addfile(ti)
            else:
                ti.size = len(data); t.addfile(ti, io.BytesIO(data))
    return b.getvalue()


def test_extract_caps():
    with tempfile.TemporaryDirectory() as d:
        f = f'{d}/a.tgz'
        open(f, 'wb').write(_tgz([('top/ok.txt', b'1 1\n0.5 0.5 0\n', 'f'), ('top/big.bin', b'x' * 2048, 'f'),
                                  ('top/link', '/etc/passwd', 'sym')]))
        old = github.MAX_MEMBER; github.MAX_MEMBER = 1024
        try:
            assert raises(github._extract_capped, f, f'{d}/x')             # absolute symlink: data filter refuses
            open(f, 'wb').write(_tgz([('top/ok.txt', b'1 1\n', 'f'), ('top/big.bin', b'x' * 2048, 'f')]))
            n, skipped = github._extract_capped(f, f'{d}/y')
            assert skipped == ['top/big.bin'] and os.path.exists(f'{d}/y/top/ok.txt')
            assert not os.path.exists(f'{d}/y/top/big.bin')
        finally:
            github.MAX_MEMBER = old
        open(f, 'wb').write(_tgz([('../escape.txt', b'x', 'f')]))
        assert raises(github._extract_capped, f, f'{d}/z') and not os.path.exists(f'{d}/escape.txt')


def test_download_cap():
    with tempfile.TemporaryDirectory() as d:
        assert raises(github._download_capped, ['head', '-c', '5000', '/dev/zero'], f'{d}/o', 1000)
        assert github._download_capped(['head', '-c', '500', '/dev/zero'], f'{d}/o', 1000) == 500


def test_bad_names():
    assert raises(github.fetch_repo, '../x/y', 'a' * 40)
    assert raises(github.fetch_repo, 'a/b', 'HEAD; rm -rf /')
    assert raises(github.fetch_attachment, 'https://evil.example/user-attachments/files/1/a', 1)


def test_text_claims():
    body = '| n | side |\n|---|---|\n| 155 | 12.952503202605 | x |\n'
    cm = ['update:\n| 155 | 12.952400000000 |\n| 154 | 12.926562245853 |', 'no table here 155 12.0000001']
    assert issues.text_claims('', body) == {155: 12.952503202605}
    assert issues.text_claims('s(271) ≤ 16.948947733714: registration request', body, cm) == \
        {155: 12.9524, 154: 12.926562245853, 271: 16.948947733714}


def test_reconcile():
    from pk.store import Store
    with tempfile.TemporaryDirectory() as d:
        st = Store(f'{d}/s.sqlite')
        a = parse('2 2\n0.5 0.5 0\n1.5 1.5 0\n', strict=True)                       # n = 2 at s = 2: the claim
        b = parse('3 3\n0.5 0.5 0\n1.5 1.5 0\n2.5 2.5 0\n', strict=True)           # n = 3: context (no claim)
        c = parse('4 2\n0.5 0.5 0\n1.5 0.5 0\n0.5 1.5 0\n1.5 1.5 0\n', strict=True)  # n = 4 at 2, claimed 2.5: below
        for p_ in (a, b, c):
            st.add(p_, 'pending:#9', 'x', 'r', 'pending')
        R = watch.reconcile(st, 'pending:#9', {2: 2.0000000000001, 4: 2.5, 5: 3.0})
        kinds = dict(st.db.execute("SELECT p.n, g.kind FROM sighting g JOIN packing p ON p.id = g.packing").fetchall())
        assert kinds == {2: 'claim', 3: 'context', 4: 'context'}, kinds
        assert set(R['matched']) == {2} and R['context'] == 2 and R['unmatched_claims'] == {4: 2.5, 5: 3.0}
        assert [x[0] for x in R['below_claim']] == [4]


def test_frontier_view():
    from pk.store import Store
    with tempfile.TemporaryDirectory() as d:
        st = Store(f'{d}/s.sqlite')
        a = parse('2 2\n0.5 0.5 0\n1.5 1.5 0\n', strict=True)
        pa, _ = st.add(a, 'pending:#9', 'alice', 'r', 'claim')
        st.set_frontier(2, 2.0000000001, '2.0000000001', 'alice', 'pending:#9', 'r', True, pa)   # backed (rounded up)
        st.set_frontier(3, 2.5, '2.5', 'bob', 'pending:#10', 'r', True, None)                  # nothing backs it
        b = parse('3 2.4\n0.5 0.5 0\n1.5 0.5 0\n0.5 1.5 0\n', strict=True)
        st.add(b, 'pool:ours', 'ours', 'r', 'seed-pool')                                       # ours: never public
        F = st.frontier()
        assert F[2]['s'] == 2 and F[2]['verified'] and F[2]['source'] == 'pending:#9'          # the packing, not the claim
        assert F[3]['s'] == 2.5 and not F[3]['verified']                                        # pool:ours is not public
        assert [r['n'] for r in st.unverified()] == [3]
        st.add(b, 'watch:bob/r@1234567', 'bob', 'r', 'watch')                                   # now published
        F = st.frontier()
        assert abs(F[3]['s'] - 2.4) < 1e-12 and F[3]['verified'] and st.unverified() == []


A = 'a' * 40
B = 'b' * 40
H = 'f' * 40


def _res(repo, ref, issue):                         # fake resolver: two commits exist on alice/c; HEAD is H
    return H if ref == 'HEAD' else next((s for s in (A, B) if s.startswith(ref)), None)


def _issue(body, comments=(), who='alice', title='t'):
    return issues.parse(dict(number=7, user=dict(login=who), state='open', title=title, body=body), comments)


def test_parse_and_units():
    body = ('| n | old | new |\n|---|---|---|\n| [17](https://github.com/alice/c/blob/aaaaaaa/n17.cert) | 4.70 | '
            '4.675900000000 |\nhttps://github.com/alice/c/blob/bbbbbbb/n18.cert and https://github.com/bob/r\n'
            'https://github.com/user-attachments/files/1/a.cert.txt')
    P = _issue(body, [('jlevy', '| 18 | 4.8000001 |'), ('alice', '| 18 | 4.8200001 |')])
    assert P['claims'] == {17: 4.6759, 18: 4.8200001}          # the maintainer's table is not alice's claim
    assert P['attachments'] == ['https://github.com/user-attachments/files/1/a.cert.txt']
    U, att, cl, notes = issues.units(P, _res)
    assert U == [dict(repo='alice/c', commit=A, paths=['n17.cert']), dict(repo='alice/c', commit=B, paths=['n18.cert']),
                 dict(repo='bob/r', commit=H, paths=['.'])], U   # one unit per linked commit; a bare repo whole at HEAD
    assert any('bob/r: pin = HEAD' in x for x in notes)
    # a path linked at a commit that no longer exists goes with the repo's pin (named commit, else HEAD)
    U, _, _, _ = issues.units(_issue('https://github.com/alice/c/blob/0123456/x.cert pinned at aaaaaaa'), _res)
    assert U == [dict(repo='alice/c', commit=A, paths=['x.cert'])], U
    # ours only from our own issues
    assert issues.units(_issue('https://github.com/evand/square-packing/tree/aaaaaaa/c'), _res)[0] == []
    assert issues.units(_issue('https://github.com/evand/square-packing/tree/aaaaaaa/c', who='evand'), _res)[0] != []


def test_overrides():
    P = _issue('| 17 | 4.680000000000 |\nhttps://github.com/alice/c/blob/aaaaaaa/n17.cert')
    U, _, cl, _ = issues.units(P, _res, dict(claims={17: 4.6759, 20: 5.1}, units=[dict(repo='x/y', commit=B)]))
    assert cl == {17: 4.6759, 20: 5.1} and U[-1] == dict(repo='x/y', commit=B, paths=['.'])
    U, _, _, _ = issues.units(P, _res, dict(units=[dict(repo='x/y', commit=B, paths=['p'])], only_override_units=True))
    assert U == [dict(repo='x/y', commit=B, paths=['p'])]
    U, att, cl, notes = issues.units(P, _res, dict(skip=True, note='refinement only'))
    assert U == [] and att == [] and 'skipped' in notes[0]


if __name__ == '__main__':
    for k, f in list(globals().items()):
        if k.startswith('test_'):
            f()
    print('ok')
