"""Ingestion safety + parsing tests (no network): python3 pk/tests/test_watch.py"""
import io, os, sys, tarfile, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pk.packing import parse, Packing
from pk import sources, watch

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
    R, loose = watch.refs(body)
    assert list(R) == ['alice/sq-certs', 'bob/x'], R
    assert R['alice/sq-certs'] == dict(shas=['0123abcd'], paths=['certs/n100'])
    assert loose == ['deadbeef12']
    assert [m.group(0) for m in sources.ATTACH_RE.finditer(body)] == ['https://github.com/user-attachments/files/123/a.cert.txt']
    R, _ = watch.refs('https://github.com/a/b/blob/main/../../etc/passwd')
    assert R['a/b']['paths'] == []


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
        old = sources.MAX_MEMBER; sources.MAX_MEMBER = 1024
        try:
            assert raises(sources._extract_capped, f, f'{d}/x')             # absolute symlink: data filter refuses
            open(f, 'wb').write(_tgz([('top/ok.txt', b'1 1\n', 'f'), ('top/big.bin', b'x' * 2048, 'f')]))
            n, skipped = sources._extract_capped(f, f'{d}/y')
            assert skipped == ['top/big.bin'] and os.path.exists(f'{d}/y/top/ok.txt')
            assert not os.path.exists(f'{d}/y/top/big.bin')
        finally:
            sources.MAX_MEMBER = old
        open(f, 'wb').write(_tgz([('../escape.txt', b'x', 'f')]))
        assert raises(sources._extract_capped, f, f'{d}/z') and not os.path.exists(f'{d}/escape.txt')


def test_download_cap():
    with tempfile.TemporaryDirectory() as d:
        assert raises(sources._download_capped, ['head', '-c', '5000', '/dev/zero'], f'{d}/o', 1000)
        assert sources._download_capped(['head', '-c', '500', '/dev/zero'], f'{d}/o', 1000) == 500


def test_bad_names():
    assert raises(sources.fetch_repo, '../x/y', 'a' * 40)
    assert raises(sources.fetch_repo, 'a/b', 'HEAD; rm -rf /')
    assert raises(sources.fetch_attachment, 'https://evil.example/user-attachments/files/1/a', 1)


if __name__ == '__main__':
    for k, f in list(globals().items()):
        if k.startswith('test_'):
            f()
    print('ok')
