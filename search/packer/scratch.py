"""Run-scoped scratch space.  Many packer scripts call tempfile.mkdtemp() per quench and never remove it (786k leftover
dirs in /tmp by 10-11).  scope() points TMPDIR (env + tempfile.tempdir, inherited by forked / spawned workers and
subprocesses) at one fresh directory and removes it at exit, including on SIGTERM / SIGHUP (pk/jobs.py stop, kill).
A SIGKILL still leaves it behind, but as one directory with a recognisable prefix (/tmp/pk-<name>-*)."""
import atexit, os, shutil, signal, sys, tempfile


def scope(name):
    root = tempfile.mkdtemp(prefix=f'pk-{name}-')
    os.environ['TMPDIR'] = root
    tempfile.tempdir = root
    pid = os.getpid()

    def clean():
        if os.getpid() == pid:              # forked workers inherit atexit handlers: only the owner removes the dir
            shutil.rmtree(root, ignore_errors=True)
    atexit.register(clean)

    def on_signal(sig, frame):
        sys.exit(128 + sig)                 # runs atexit (clean) in the owner; workers just exit
    for sig in (signal.SIGTERM, signal.SIGHUP):
        if signal.getsignal(sig) in (signal.SIG_DFL, None):
            signal.signal(sig, on_signal)
    return root
