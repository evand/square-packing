"""Tests for pk/identity.py and the store's minimum / copy tables (10-11).  python pk/tests/test_identity.py"""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pk import identity as ID


def test_classify_status():
    c = ID.classify_status
    assert c('local minimum of S, strict modulo 17 exact flat motions') == ID.CERTIFIED
    assert c('PD: strict local minimum (second-order sufficient)') == ID.CERTIFIED
    assert c('rigid (null(J_A) = 0)') == ID.CERTIFIED
    assert c('PSD with 21 zero modes (degenerate; ...)') == ID.WEAK
    assert c('inconclusive: negative curvature blocked (first order) by weakly active contacts') == ID.WEAK
    assert c('NOT a local minimum: second-order descent direction') == ID.NOT_MIN
    assert c('not jammed') == ID.NOT_MIN
    assert c('solve failed') == ID.UNRESOLVED
    assert c('local minimum of S', has_S=False) == ID.UNRESOLVED          # no exact S: no verdict
    assert c('local minimum of S', cert_valid=False) == ID.UNRESOLVED
    assert ID.classify(dict(status='not jammed')) == ID.NOT_MIN
    assert ID.classify(dict(S_exact='10.9', cert_valid=True, second_order=dict(status='PSD with 3 zero modes'))) == ID.WEAK


def test_min_key():
    S = '10.99678339663159163950277843399'
    assert ID.min_key(S) == ID.min_key(S[:23]) == ID.min_key('10.996783396631591640')     # truncated / rounded sources
    a, b, c = ID.same_min_keys(S)
    assert a < b < c and b == ID.min_key(S)
    assert ID.min_key('10.99678339663159163') != ID.min_key('10.99678339663159170')       # 7e-17 apart: distinct


def test_store_minimum():
    with tempfile.TemporaryDirectory() as d:
        from pk.store import Store
        st = Store(os.path.join(d, 's.sqlite'))
        k1, new1 = st.add_minimum(110, '10.99678339663159163950', ID.WEAK, 'PSD with 3 zero modes')
        k2, new2 = st.add_minimum(110, '10.99678339663159163950277843399', ID.CERTIFIED, 'local minimum of S')
        k3, new3 = st.add_minimum(110, '10.996783396631591639', ID.WEAK, 'inconclusive')          # 1e-17 neighbour
        assert new1 and not new2 and not new3 and k1 == k2
        m = st.minima(110)
        assert len(m) == 1 and m[0]['cls'] == ID.CERTIFIED                                         # upgraded, never down
        assert m[0]['S'] == '10.99678339663159163950277843399'                                     # longest S kept
        mk = st.record_verdict(110, '10.996783397', dict(cls='certified', S='10.99678339663159163950', status='PD'), 'test')
        assert mk == k1 and st.basin(110, '10.996783397')['cls'] == 'certified'
        assert st.record_verdict(110, '10.998000000', dict(cls='not-min', status='not jammed'), 'test') is None
        assert len(st.db.execute('SELECT * FROM copy').fetchall()) == 1
        st.close()


if __name__ == '__main__':
    for f in (test_classify_status, test_min_key, test_store_minimum):
        f(); print('ok', f.__name__)
