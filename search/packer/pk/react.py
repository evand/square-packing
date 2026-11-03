#!/usr/bin/env python3
"""Follow-up work for new store rows (11-01): one pass reads the rows added since the last pass and hands them to each
rule, which queues pool units (pk/tasks.py).  Rules are small and independent; a plugin adds its own with @rule.

  python3 pk/tasks.py repeat _react --every 60 --prio high -- python3 -m pk.react      keep it running (once)
  python3 -m pk.react hold | release          a hold keeps experiments isolated: rows that arrive during it get no
                                              follow-ups, ever (release does not catch up).  For one experiment
                                              instead: rows tagged react:hold are skipped (`pk.py try --hold`)

Cursor per rule: runs/tasks/react.json (a new rule starts at the current end of the store, without a backlog).
"""
from __future__ import annotations
import fcntl, json, os, sys
from .store import Store
from . import tasks

CURSOR = os.path.join(tasks.DIR, 'react.json')
HOLD = os.path.join(tasks.DIR, 'react.hold')
HOLD_TAG = '"react:hold"'  # a row tagged react:hold never gets follow-ups (one experiment held, the rest live; 11-02)
PK_DIR = os.path.dirname(tasks.HERE)
MAX_ROWS = 5000          # rows per rule per pass; the rest wait for the next pass
RULES = {}


def rule(name):
    def deco(f):
        RULES[name] = f
        return f
    return deco


@rule('converge')
def _converge(st, ids, batch=8, cap=50):
    """near-record rows -> units of a few packings each (about 5 CPU-s per packing)"""
    from . import converge
    J = [i for i, _, _ in converge.select(st, ids=ids, cap=cap)]
    cmds = [[sys.executable, 'pk.py', 'converge', '--serial', '--ids', *map(str, J[k:k + batch])]
            for k in range(0, len(J), batch)]
    if cmds:
        tasks.add('converge', cmds, cwd=PK_DIR)
    return len(J)


@rule('nearmiss')
def _nearmiss(st, ids):
    from . import nearmiss as NM
    pub = {n: r['s'] for n, r in st.frontier().items()}
    return NM.enqueue(st, NM.scan(st.db, pub.get, ids=ids), note='react')


def run_pass():
    """-> rows handed to rules"""
    from . import stages  # noqa: F401  (loads plugins: their rules)
    st = Store()
    end = st.db.execute('SELECT MAX(id) FROM packing').fetchone()[0]
    cur = json.load(open(CURSOR)) if os.path.exists(CURSOR) else {}
    held, rows = os.path.exists(HOLD), 0
    for name, f in RULES.items():
        c = cur.get(name, end)
        R = st.db.execute('SELECT id, tags FROM packing WHERE id>? ORDER BY id LIMIT ?', (c, MAX_ROWS)).fetchall()
        ids = [i for i, t in R if HOLD_TAG not in t]          # rows of a held experiment (`pk.py try --hold`) skipped
        if ids and not held:
            print(f'{name}: {len(ids)} rows -> {f(st, ids)}', flush=True)
            rows += len(ids)
        cur[name] = R[-1][0] if R else c
        with open(CURSOR + '.tmp', 'w') as o:          # per rule: a failing rule does not repeat the others
            json.dump(cur, o)
        os.replace(CURSOR + '.tmp', CURSOR)
    return rows


def main():
    a = sys.argv[1:]
    if a == ['hold']:
        open(HOLD, 'w').close()
    elif a == ['release']:
        if os.path.exists(HOLD):
            os.remove(HOLD)
    else:
        with open(CURSOR + '.lock', 'w') as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            sys.exit(0 if run_pass() else tasks.EX_IDLE)


if __name__ == '__main__':                     # the imported module, so plugin rules register into the same RULES
    from pk import react as _self
    _self.main()
