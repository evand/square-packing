#!/usr/bin/env python3
"""pk/tasks.py with real processes in a temporary pool: priority order, sharing, idle exits, failure pause, budgets."""
import os, sqlite3, subprocess, sys, tempfile, time

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tasks.py')
ENV = dict(os.environ, PK_TASKS=tempfile.mkdtemp(prefix='pktasks_'), PK_TASKS_NO_OUTSIDE='1')


def t(*args, check=True):
    return subprocess.run([sys.executable, TOOL, *args], env=ENV, capture_output=True, text=True, check=check, timeout=120).stdout


def stamp(name):
    """a unit that records its start"""
    return ['bash', '-c', f'echo {name} $(date +%s.%N) >> {ENV["PK_TASKS"]}/order; sleep 0.3']


def order():
    return [l.split() for l in open(ENV['PK_TASKS'] + '/order')]


def square(x):
    return x * x


def boom(_):
    raise ValueError('boom')


def die(_):
    os._exit(9)


G = 0


def setg(v):
    global G
    G = v


def plusg(x):
    return x + G


def power(b, e=1):
    return b ** e


def hang(_):
    return sum(range(10 ** 12))                 # one C call: no Python signal handler runs before it returns


def slow(x):
    time.sleep(1.5); return x


def main():
    t('limit', '0')                            # nothing starts until all three wait
    for c, prio in (('lo', 'low'), ('no', 'normal'), ('hi', 'high')):
        t('add', c, '--prio', prio, '--', *stamp(c))
    t('limit', '1')
    for c in ('lo', 'no', 'hi'):
        assert t('wait', c).strip().endswith('done'), c
    got = [o[0] for o in order()]
    assert got == ['hi', 'no', 'lo'], got

    # equal priority: two campaigns alternate instead of one running first
    open(ENV['PK_TASKS'] + '/order', 'w').close()
    t('limit', '0')
    for c in ('a', 'b'):
        for _ in range(4):
            t('add', c, '--', *stamp(c))
    t('limit', '2')
    t('wait', 'a'); t('wait', 'b')
    got = [o[0] for o in order()]
    assert sorted(got[:4]) == ['a', 'a', 'b', 'b'], got

    # failures in a row pause a campaign and keep the output
    t('limit', '1')
    for _ in range(8):
        t('add', 'bad', '--', 'bash', '-c', 'echo oops; exit 3')
    r = subprocess.run([sys.executable, TOOL, 'wait', 'bad'], env=ENV, capture_output=True, text=True, timeout=120)
    assert r.returncode == 1, r.stdout
    assert 'oops' in open(os.path.join(ENV['PK_TASKS'], 'bad', sorted(f for f in os.listdir(os.path.join(ENV['PK_TASKS'], 'bad'))
                                                                     if f.startswith('u'))[0])).read()
    assert 'paused' in t('ls')
    t('add', 'bad', '--', 'true'); assert 'paused' in [l for l in t('ls').splitlines() if l.startswith('bad')][0]    # stays paused
    t('cancel', 'bad')

    # a command that cannot start is a failed unit, not a dead dispatcher
    t('add', 'nosuch', '--', '/no/such/binary')
    t('add', 'after', '--', 'true')
    assert t('wait', 'after').strip().endswith('done')

    # repeat runs until its budget (near zero: the first finished unit ends it)
    t('repeat', 'rep', '--budget', '0.0000001', '--', 'bash', '-c', 'for i in $(seq 200000); do :; done')
    assert t('wait', 'rep').strip().endswith('done')
    assert 'budget spent' in t('ls')

    # an idle exit backs off (one unit in about 15 s, not a spin) and costs no failure
    t('repeat', 'idle', '--budget', '1', '--', 'bash', '-c', 'exit 75')
    time.sleep(3)
    db = sqlite3.connect(ENV['PK_TASKS'] + '/tasks.sqlite')
    fails, units, next_t = db.execute("SELECT fails, units, next_t FROM campaign WHERE name='idle'").fetchone()
    assert (fails, units) == (0, 0) and next_t > time.time(), (fails, units, next_t)
    t('cancel', 'idle')

    # --every runs one unit at a time: a unit that outlives the interval does not get a sibling (11 hung ticks held 11 slots, 11-02)
    t('limit', '4')
    t('repeat', 'tick', '--every', '1', '--prio', 'high', '--', 'bash', '-c', 'sleep 4')
    time.sleep(3.5)
    assert db.execute("SELECT COUNT(*) FROM unit WHERE campaign='tick' AND state='running'").fetchone()[0] == 1
    t('cancel', 'tick')

    # many units ending at once: every slot comes back (a Popen dropped too early reaped children behind wait4)
    t('limit', '16')
    subprocess.run([sys.executable, TOOL, 'add', 'burst', '--file', '-'], env=ENV, input='true\n' * 400, text=True, check=True, capture_output=True)
    assert t('wait', 'burst').strip().endswith('done')
    assert db.execute("SELECT units FROM campaign WHERE name='burst'").fetchone()[0] == 400

    # Executor: a function of the calling script, a raised exception, a unit that dies
    from importlib import import_module
    sys.path.insert(0, os.path.join(os.path.dirname(TOOL), '..'))
    os.environ.update(PK_TASKS=ENV['PK_TASKS'], PK_TASKS_NO_OUTSIDE='1')
    tasks = import_module('pk.tasks')
    with tasks.Executor('calls') as ex:
        assert list(ex.map(square, [2, 3, 4])) == [4, 9, 16]
        try:
            ex.submit(boom, 1).result(60); raise SystemExit('no exception')
        except ValueError as e:
            assert 'boom' in str(e)
        try:
            ex.submit(die, 1).result(60); raise SystemExit('no exception')
        except RuntimeError as e:
            assert 'without a result' in str(e)

    # 11-02 (private pools removed): kwargs, chunked map, initializer per unit, a timeout that kills a hang in C code,
    # two executors sharing a campaign, leaving the block waits, a call pickle from an older Executor
    with tasks.Executor('calls2', initializer=setg, initargs=(7,), timeout=3) as ex:
        assert ex.submit(power, 2, e=3).result(60) == 8
        os.environ['PK_T_ENV'] = 'yes'; assert ex.submit(os.getenv, 'PK_T_ENV').result(60) == 'yes'
        assert list(ex.map(plusg, range(10), chunksize=4)) == [x + 7 for x in range(10)]
        time.sleep(1); assert db.execute("SELECT units FROM campaign WHERE name='calls2'").fetchone()[0] == 5   # 2 + 3 chunks
        t1 = time.time()
        try:
            ex.submit(hang, 1).result(60); raise SystemExit('no exception')
        except TimeoutError as e:
            assert 'killed' in str(e) and time.time() - t1 < 15, (e, time.time() - t1)
        time.sleep(0.5)
        assert db.execute("SELECT fails FROM campaign WHERE name='calls2'").fetchone()[0] == 0    # a timeout is no failure
        late = ex.submit(slow, 5)
        other = tasks.Executor('calls2'); other.dir += 'b'; os.makedirs(other.dir)
        f2 = other.submit(slow, 6); other.shutdown(cancel_futures=True)
    assert late.done() and late.result() == 5
    assert tasks.UNIT_CPU > 2, tasks.UNIT_CPU                # the hang's CPU and the others'
    import pickle
    p = os.path.join(ENV['PK_TASKS'], 'legacy'); os.makedirs(p, exist_ok=True); p += '/1'
    with open(p + '.call', 'wb') as o:
        pickle.dump(os.path.abspath(__file__), o); pickle.dump((square, (9,)), o)
    subprocess.run([sys.executable, TOOL, '_call', p], env=ENV, check=True)
    assert pickle.load(open(p + '.done', 'rb')) == (True, 81)

    # a unit that keeps dying pauses the campaign: waiting calls fail instead of hanging; a new run clears it; cancel clears it too
    ex = tasks.Executor('crash')
    fs = [ex.submit(die, 1) for _ in range(12)]
    msgs = []
    for f in fs:
        try:
            f.result(60); raise SystemExit('no exception')
        except RuntimeError as e:
            msgs.append(str(e))
    assert any('paused' in m for m in msgs), msgs
    ex.stop = True; ex.poller.join()
    assert db.execute("SELECT status FROM campaign WHERE name='crash'").fetchone()[0] == 'paused'
    with tasks.Executor('crash') as ex:
        assert ex.submit(square, 5).result(60) == 25
    t('add', 'crash2', '--', 'false')
    for _ in range(6): t('add', 'crash2', '--', 'false')
    assert t('wait', 'crash2', check=False).strip().endswith('paused')
    t('cancel', 'crash2'); assert db.execute("SELECT status, consec FROM campaign WHERE name='crash2'").fetchone() == ('cancelled', 0)

    # drain: limit 0 starts nothing
    t('limit', '0'); t('add', 'held', '--', 'true'); time.sleep(1.5)
    assert ' 1 ' in [l for l in t('ls').splitlines() if l.startswith('held')][0]
    t('limit', '2'); assert t('wait', 'held').strip().endswith('done')
    print('tasks: ok')


if __name__ == '__main__':
    main()
