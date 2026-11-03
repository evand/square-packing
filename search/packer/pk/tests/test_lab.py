"""pk/lab.py and Store.record (11-02) on a temp store: one real application (lab.apply: move-log row, child with parent
link and tags, params with arm / rep / seed), Store.record losses, the paired report, the stage table x move log, and
react skipping rows tagged react:hold.  Needs fq (one relax at n = 50) and a store with n = 50 rows; no pool.

  python3 pk/tests/test_lab.py
"""
import json, os, sqlite3, sys, tempfile
import numpy as np
D = tempfile.mkdtemp(prefix='pklab_')
os.environ['PK_STORE'] = os.path.join(D, 'store.sqlite')                  # before pk.store is imported
os.environ['PK_TASKS'] = os.path.join(D, 'tasks')                         # react's cursor / hold files
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pk import lab, react
from pk.packing import Packing
from pk.store import Store


def main():
    st = Store()
    real = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'runs', 'store.sqlite')
    db = sqlite3.connect(f'file:{os.path.abspath(real)}?mode=ro', uri=True)    # a sub-grid n = 50 start, dilated 0.2 %
    s0, c = db.execute("SELECT s, coords FROM packing WHERE n=50 AND status!='infeasible' ORDER BY s LIMIT 1").fetchone()
    db.close()
    sq = np.frombuffer(c, '<f8').reshape(-1, 3).copy(); sq[:, :2] *= 1.002
    pid, _ = st.add(Packing(s0 * 1.002, sq), 'test'); st.close()
    r = lab.apply('t1', 'ctl', 'relax.round.15c.p5', pid, 0, 123, 120.0, ('try:t1', 'react:hold'))
    assert r['s'] is not None and r['s'] < s0 * 1.002 and r['child'], r
    st = Store()
    e = st.db.execute("SELECT mechanism, params, child, run FROM expansion WHERE packing=?", (pid,)).fetchone()
    P = json.loads(e['params'])
    assert e['mechanism'] == 'relax.round.15c.p5' and e['run'] == 'try:t1' and P['arm'] == 'ctl' and P['seed'] == 123, dict(e)
    c = st.db.execute('SELECT parent, tags FROM packing WHERE id=?', (r['child'],)).fetchone()
    assert c['parent'] == pid and 'react:hold' in c['tags']
    # paired rows by hand: arm a1 better than ctl on rep 0, a loss (no output) on rep 1
    st.record(pid, 'x | relax.round.15c.p5', dict(arm='a1', rep=0, seed=123), None, 2.0, run='try:t1', s=r['s'] - 1e-6)
    st.record(pid, 'x | relax.round.15c.p5', dict(arm='a1', rep=1, seed=5), None, 2.0, run='try:t1')
    st.record(pid, 'relax.round.15c.p5', dict(arm='ctl', rep=1, seed=5), None, 2.0, run='try:t1', s=s0 * 1.001)
    st.close()
    rep = lab.report('t1')
    row = [l for l in rep.splitlines() if l.startswith('| a1 | 50 |')][0]
    assert '| 2 | 1 |' in row and '1 / 1 / 0' in row, rep              # tries 2, failed 1; one pair up, one down
    S = {x['name']: x for x in lab.stage_rows()}
    assert S['relax.round.15c.p5']['apps'] == 4 and S['relax.round.15c.p5']['children'] == 1, S['relax.round.15c.p5']
    assert all(x['status'] in ('lab', 'default', 'retired') for x in S.values())
    assert lab.stage_table(lab.stage_rows(), md=True).startswith('| stage |')
    seen = []                                                             # react: the held child is not handed to rules
    react.RULES.clear(); react.RULES['t'] = lambda st_, ids: seen.extend(ids) or len(ids)
    os.makedirs(os.path.dirname(react.CURSOR), exist_ok=True); open(react.CURSOR, 'w').write(json.dumps({'t': 0}))
    react.run_pass()
    assert pid in seen and r['child'] not in seen, (seen, r['child'])
    assert json.load(open(react.CURSOR))['t'] == r['child']              # the cursor moves past held rows
    # 11-02: arm syntax, cross-n rows scored against the target n's reference, determinism, features
    assert lab.arm_spec('dn=-2: relax.round.15c.p5') == (-2, 'relax.round.15c.p5') and lab.arm_spec(' relax.polish ') == (0, 'relax.polish')
    from pk import stages as S, features as FT, events as E
    assert S.deterministic('relax.round.15c.p5') and not S.deterministic('perturb.kick | relax.round.15c.p5')
    assert E.below(1.0, 1.0 + 1e-9) and not E.below(1.0, 1.0 + 1e-12) and not E.below(None, 1.0)
    st = Store()
    st.record(pid, 'shrink.drop | relax.round.15c.p5', dict(arm='a2', rep=0, seed=123, dn=-1), None, 2.0, run='try:t1', s=1.0)
    st.close()
    row = [l for l in lab.report('t1').splitlines() if l.startswith('| a2 | 49 |')][0]
    assert '| 1 | 0 |  | 1 |' in row and '1 / 0 / 0' in row, row        # a win at 49: scored vs 49's reference
    FT.FEATURES.clear()
    FT.feature('nsq')(lambda s, sq: len(sq)); FT.feature('big')(lambda s, sq: len(sq) > 100)
    assert FT.filter_ids([pid], ['nsq>=50', '!big', 'nsq=50']) == [pid] and FT.filter_ids([pid], ['big']) == []
    print('lab ok')


if __name__ == '__main__':
    main()
