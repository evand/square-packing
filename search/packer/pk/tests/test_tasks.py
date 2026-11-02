#!/usr/bin/env python3
"""pk/tasks.py with real processes in a temporary pool: priority order, sharing, idle exits, failure pause, budgets."""
import os, sqlite3, subprocess, sys, tempfile, time

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tasks.py')
ENV = dict(os.environ, PK_TASKS=tempfile.mkdtemp(prefix='pktasks_'))


def t(*args, check=True):
    return subprocess.run([sys.executable, TOOL, *args], env=ENV, capture_output=True, text=True, check=check, timeout=120).stdout


def stamp(name):
    """a unit that records its start"""
    return ['bash', '-c', f'echo {name} $(date +%s.%N) >> {ENV["PK_TASKS"]}/order; sleep 0.3']


def order():
    return [l.split() for l in open(ENV['PK_TASKS'] + '/order')]


def main():
    t('limit', '2')
    # a blocker fills both slots, so that everything below waits and then starts in priority order
    subprocess.run([sys.executable, TOOL, 'add', 'block', '--file', '-'], env=ENV, input='sleep 1.5\nsleep 1.5\n', text=True,
                   check=True, capture_output=True)
    time.sleep(0.5)
    for c, prio in (('lo', 'low'), ('no', 'normal'), ('hi', 'high')):
        t('add', c, '--prio', prio, '--', *stamp(c))
    for c in ('block', 'lo', 'no', 'hi'):
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

    # many units ending at once: every slot comes back (a Popen dropped too early reaped children behind wait4)
    t('limit', '16')
    subprocess.run([sys.executable, TOOL, 'add', 'burst', '--file', '-'], env=ENV, input='true\n' * 400, text=True, check=True, capture_output=True)
    assert t('wait', 'burst').strip().endswith('done')
    assert db.execute("SELECT units FROM campaign WHERE name='burst'").fetchone()[0] == 400

    # drain: limit 0 starts nothing
    t('limit', '0'); t('add', 'held', '--', 'true'); time.sleep(1.5)
    assert ' 1 ' in [l for l in t('ls').splitlines() if l.startswith('held')][0]
    t('limit', '2'); assert t('wait', 'held').strip().endswith('done')
    print('tasks: ok')


if __name__ == '__main__':
    main()
