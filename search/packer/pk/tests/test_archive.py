"""pk/archive: dedupe removes only byte-identical files; reads write them back."""
import json, os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pk import archive
from pk.store import Store
from pk.packing import Packing, read
import mcmin


def test_roundtrip(monkeypatch):
    with tempfile.TemporaryDirectory() as d:
        db = os.path.join(d, 'store.sqlite'); monkeypatch.setenv('PK_STORE', db)
        run = os.path.join(d, 'run'); os.makedirs(run)
        st = Store(db)
        sq = [(0.5, 0.5, 0.0), (1.5000000000000004, 0.5, 0.0), (0.5, 1.5, 89.99999999999999), (1.6, 1.6, 12.5)]
        E = []
        for i, s in enumerate((2.3, 2.4, 2.5)):
            f = os.path.join(run, f'b{i:05d}.txt'); mcmin.write_deg(f, s, sq); E.append(dict(i=i, s=s, path=f))
        open(os.path.join(run, 'archive.jsonl'), 'w').write(''.join(json.dumps(e) + '\n' for e in E))
        src = 'run:' + os.path.relpath(run, archive.PKG)
        for e in E[:2]:                                   # entry 2 never imported
            st.add(read(e['path']), src, 'ours', str(e['i']), check=False)
        st.db.execute('UPDATE packing SET s=s+1e-12 WHERE s>2.39 AND s<2.41'); st.db.commit()   # entry 1: store differs
        st.close()
        before = {e['path']: open(e['path'], 'rb').read() for e in E}
        r = archive.dedupe_run(run, apply=True, min_age_h=0, store=db)
        assert (r['removed'], r['not_in_store'], r['differs']) == (1, 1, 1), r
        assert not os.path.exists(E[0]['path'])
        assert read(E[0]['path']).n == 4 and open(E[0]['path'], 'rb').read() == before[E[0]['path']]
        os.remove(E[0]['path']); assert archive.restore(run) == 1
        assert all(open(p, 'rb').read() == b for p, b in before.items())
