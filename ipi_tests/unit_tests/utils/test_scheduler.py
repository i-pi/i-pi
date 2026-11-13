"""The cooperative scheduler runs tasks in one thread, switching between them
while they wait, and polls the forcefields meanwhile."""

import contextvars
import socket

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import FFDirect
from ipi.utils import scheduler


class Counter:
    """A forcefield whose requests are numbers, completed one per poll."""

    name = "counter"
    poll_interval = None

    def __init__(self, idle_only=False):
        self.requests = []
        self.polls = []
        self.idle_only = idle_only

    def queue(self, value):
        request = {"value": value, "done": False}
        self.requests.append(request)
        scheduler.watch(self)
        return request

    def poll(self, idle=False):
        self.polls.append(idle)
        if self.idle_only and not idle:
            return False
        pending = [r for r in self.requests if not r["done"]]
        if pending:
            pending[0]["done"] = True
            self.requests.remove(pending[0])
        return bool(pending)

    def fds(self):
        return []


@pytest.fixture(autouse=True)
def fresh_scheduler(monkeypatch):
    """Each test gets its own scheduler: one whose tasks are stopped leaves
    their requests queued, and the forcefield watched."""
    monkeypatch.setattr(
        scheduler._local, "scheduler", scheduler.Scheduler(), raising=False
    )


def test_wait_in_the_root_polls_without_tasks():
    ff = Counter()
    requests = [ff.queue(i) for i in range(3)]
    scheduler.wait(lambda: all(r["done"] for r in requests))
    assert all(r["done"] for r in requests)
    # nothing else can run, so the forcefield is polled as idle
    assert ff.polls == [True, True, True]
    # a condition that holds returns at once
    scheduler.wait(lambda: True)
    assert len(ff.polls) == 3


def test_tasks_switch_while_they_wait():
    ff = Counter(idle_only=True)
    events = []

    def task(name):
        events.append(f"{name} queues")
        request = ff.queue(name)
        scheduler.wait(lambda: request["done"])
        events.append(f"{name} done")
        return name

    assert scheduler.run([lambda: task("a"), lambda: task("b")]) == ["a", "b"]
    # both tasks queue their request before any is evaluated
    assert events[:2] == ["a queues", "b queues"]
    assert sorted(events[2:]) == ["a done", "b done"]


def test_tasks_make_different_numbers_of_requests():
    ff = Counter()
    counts = {}

    def task(name, n):
        for i in range(n):
            request = ff.queue((name, i))
            scheduler.wait(lambda: request["done"])
        counts[name] = n

    scheduler.run([lambda n=n: task(n, n) for n in [0, 1, 5, 2]])
    assert counts == {0: 0, 1: 1, 5: 5, 2: 2}
    assert not ff.requests


def test_failure_of_a_task_stops_the_others():
    ff = Counter(idle_only=True)
    stopped = []

    def waits():
        request = ff.queue(1)
        try:
            scheduler.wait(lambda: False)
        finally:
            stopped.append(request)

    def fails():
        raise ValueError("task failed")

    with pytest.raises(ValueError, match="task failed"):
        scheduler.run([waits, fails])
    assert len(stopped) == 1


def test_failure_of_a_subtask_reaches_the_task_that_joins_it():
    ff = Counter(idle_only=True)

    def waits():
        ff.queue(0)
        scheduler.wait(lambda: False)

    def fails():
        raise ValueError("subtask failed")

    def outer():
        try:
            scheduler.join([scheduler.spawn(waits), scheduler.spawn(fails)])
        except ValueError:
            return "caught"

    assert scheduler.run([outer, lambda: "other"]) == ["caught", "other"]


def test_tasks_spawned_by_tasks():
    ff = Counter()

    def inner(i):
        request = ff.queue(i)
        scheduler.wait(lambda: request["done"])
        return i * i

    def outer(n):
        return sum(scheduler.join([scheduler.spawn(inner, i) for i in range(n)]))

    assert scheduler.run([lambda: outer(3), lambda: outer(4)]) == [5, 14]


def test_waiters_of_one_resource_check_it_again():
    """A task resumed because a resource was free finds it taken if another
    task resumed before took it."""

    free = [1]
    ff = Counter()
    taken = []

    def take(name):
        scheduler.wait(lambda: free[0] > 0)
        free[0] -= 1
        taken.append(name)
        request = ff.queue(name)
        scheduler.wait(lambda: request["done"])
        free[0] += 1

    scheduler.run([lambda n=n: take(n) for n in range(4)])
    assert sorted(taken) == [0, 1, 2, 3]
    assert free == [1]


def test_tasks_see_the_context_of_the_caller():
    variable = contextvars.ContextVar("variable", default="unset")
    variable.set("caller")
    assert scheduler.run([variable.get]) == ["caller"]


def test_waiting_for_what_cannot_happen_raises():
    with pytest.raises(RuntimeError, match="cannot be computed"):
        scheduler.wait(lambda: False)


def test_waiting_while_polling_raises():
    class Waits(Counter):
        def poll(self, idle=False):
            scheduler.wait(lambda: False)

    ff = Waits()
    request = ff.queue(0)
    with pytest.raises(RuntimeError, match="while polling"):
        scheduler.wait(lambda: request["done"])
    scheduler.unwatch(ff)


def test_errors_of_a_driver_reach_the_task_that_waits(tmp_path):
    path = tmp_path / "failing.py"
    path.write_text(
        "from ipi.pes.dummy import Dummy_driver\n"
        "__DRIVER_NAME__ = 'failing'\n"
        "__DRIVER_CLASS__ = 'Failing'\n"
        "class Failing(Dummy_driver):\n"
        "    def compute_structure(self, cell, pos):\n"
        "        raise ArithmeticError('driver failed')\n"
    )
    ff = FFDirect(pes="custom", pes_path=str(path))

    def task():
        request = ff.queue(Atoms(1), Cell(np.eye(3)))
        scheduler.wait(lambda: request["status"] == "Done")

    with pytest.raises(ArithmeticError, match="driver failed"):
        scheduler.run([task, task])
    ff.stop()


def test_requests_of_all_tasks_are_evaluated_in_one_batch():
    ff = FFDirect(pes="harmonic", pars={"k1": 1.0}, batch_size=4)
    calls = []
    driver = ff.driver

    def counting_driver(cell, pos):
        calls.append(len(pos))
        return driver(cell, pos)

    ff.driver = counting_driver

    def task(x, nrequests):
        atoms = Atoms(1)
        atoms.q_[:] = [x, 0, 0]
        requests = [ff.queue(atoms, Cell(np.eye(3))) for _ in range(nrequests)]
        scheduler.wait(lambda: all(r["status"] == "Done" for r in requests))
        for r in requests:
            ff.release(r)
        return requests[0]["result"][0]

    energies = scheduler.run([lambda x=x: task(x, 1 + x) for x in range(3)])
    np.testing.assert_allclose(energies, [0.5 * x**2 for x in range(3)])
    # 6 requests: a full batch of 4 as soon as possible, and then the rest
    # once all the tasks wait
    assert calls == [4, 2]


class Polled(Counter):
    """A forcefield that cannot be waited on (e.g. MPI), polled every
    `poll_interval`, whose request is done at its third idle poll. It may
    also have sockets to wait on, which never get a reply."""

    def __init__(self, interval, fds=()):
        super().__init__(idle_only=True)
        self.poll_interval = interval
        self.sockets = list(fds)
        self.idle_polls = 0

    def poll(self, idle=False):
        self.polls.append(idle)
        if not idle:
            return False
        self.idle_polls += 1
        if self.idle_polls % 3 and self.requests:
            return False
        return super().poll(idle)

    def fds(self):
        return self.sockets


class RecordingPoller:
    """Stands for select.poll(), recording the time-out of each call."""

    def __init__(self, calls):
        self.calls = calls

    def register(self, fd, mask):
        pass

    def poll(self, timeout=None):
        self.calls.append(timeout)
        return []


@pytest.mark.parametrize("with_socket", [False, True], ids=["no-socket", "socket"])
@pytest.mark.parametrize("interval", [1e-4, 2e-3])
def test_poll_intervals_are_not_rounded_to_milliseconds(
    monkeypatch, interval, with_socket
):
    """select.poll counts in whole milliseconds and rounds up, so a shorter
    `poll_interval` (the latency of MPI drivers, 1e-4 s by default) is
    slept instead, after checking the sockets without waiting."""
    sleeps, polls = [], []
    monkeypatch.setattr(scheduler.time, "sleep", sleeps.append)
    monkeypatch.setattr(scheduler.select, "poll", lambda: RecordingPoller(polls))
    pair = socket.socketpair() if with_socket else ()
    try:
        ff = Polled(interval, fds=pair[:1])
        request = ff.queue(0)
        scheduler.wait(lambda: request["done"])
    finally:
        for s in pair:
            s.close()
    # two idle polls without progress, each followed by a wait
    if not with_socket:
        assert sleeps == [interval, interval] and polls == []
    elif interval < scheduler.MINPOLL:
        assert polls == [0, 0] and sleeps == [interval, interval]
    else:
        assert polls == pytest.approx([interval * 1e3] * 2) and sleeps == []
