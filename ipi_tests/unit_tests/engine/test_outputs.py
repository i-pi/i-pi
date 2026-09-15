import threading
from types import SimpleNamespace

import numpy as np

from ipi.engine.outputs import TrajectoryOutput
from ipi.utils.softexit import softexit


class _BlockingStream:
    """A stream that pauses its first write until the test releases it."""

    def __init__(self):
        self.write_started = threading.Event()
        self.release_write = threading.Event()
        self.close_called = threading.Event()
        self.closed = False

    def write(self, data):
        """Blocks a write and rejects it if the stream was closed meanwhile."""

        if self.closed:
            raise ValueError("I/O operation on closed file")
        self.write_started.set()
        if not self.release_write.wait(timeout=5):
            raise TimeoutError("timed out waiting to release stream write")
        if self.closed:
            raise ValueError("I/O operation on closed file")
        return len(data)

    def close(self):
        """Marks the stream as closed."""

        self.closed = True
        self.close_called.set()


def test_trajectory_close_waits_for_pending_write(monkeypatch):
    """Closing an output must wait for an in-flight trajectory frame."""

    monkeypatch.setattr(softexit, "triggered", False)
    stream = _BlockingStream()
    output = TrajectoryOutput(what="positions", format="xyz", ibead=0, flush=0)
    output._stream_lock = threading.RLock()
    output.out = [stream]
    output.system = SimpleNamespace(
        simul=SimpleNamespace(step=0),
        trajs={
            "positions": (
                np.zeros((1, 1, 3)),
                "length",
                "atomic_unit",
            )
        },
        beads=SimpleNamespace(
            natoms=1,
            nbeads=1,
            names=np.array(["H"]),
        ),
        cell=SimpleNamespace(h=np.eye(3)),
    )
    write_errors = []

    def write_output():
        """Records errors raised by the writer thread."""

        try:
            output.write()
        except BaseException as error:
            write_errors.append(error)

    writer = threading.Thread(target=write_output)
    writer.start()
    assert stream.write_started.wait(timeout=2)

    closer = threading.Thread(target=output.close_stream)
    closer.start()
    try:
        assert not stream.close_called.wait(timeout=0.2)
    finally:
        stream.release_write.set()

    writer.join(timeout=2)
    closer.join(timeout=2)

    assert not writer.is_alive()
    assert not closer.is_alive()
    assert write_errors == []
    assert stream.closed
