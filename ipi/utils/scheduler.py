"""A cooperative scheduler for force evaluations, built on greenlets.

Work that can overlap, such as the steps of the systems of a replica-exchange
run, runs in tasks. A task is a function run by a greenlet (a `Worker`): the
greenlets of a thread take turns, and switch only when a task waits, usually
for forces. Nothing runs in parallel, so the tasks need no locks.

Control moves between the tasks and the root greenlet, the one in which the
thread started:

- A task that calls `wait(condition)` first polls the forcefields, which
  sends its requests to the drivers, or evaluates them in process. If the
  condition is still false, the task is parked in `waiting`, and switches to
  the root.
- The root runs `_loop`. It resumes the parked tasks whose condition holds,
  and starts the new tasks, one at a time; each runs until it waits again or
  ends. Once every task waits, the root polls the forcefields with
  `idle=True`, which evaluates the incomplete batches, and then blocks in one
  `select.poll` over the sockets of all the drivers, until a reply arrives.
- `wait` called in the root runs `_loop` there, so that a run without tasks
  computes its forces where it requests them, without switching.

The loop has no notion of steps or of systems: tasks can make any number of
requests, and each one only waits for its own. Waits also make the soft
exits that have been requested (`ipi.utils.softexit`), which stop the
forcefields, and so end the waits.

A forcefield is polled while it has `requests`, from the moment it queues
one (`watch`). It takes part through three members:

- `poll(idle)` makes all the progress it can without blocking, and returns
  whether anything changed. `idle` is true when every task waits, which is
  when incomplete batches are evaluated. It must neither wait nor read depend
  objects, which could compute forces.
- `fds()` returns the sockets on which a reply or a connection is expected.
- `poll_interval` is None, or the time between polls of drivers that have no
  socket to wait on (e.g. MPI).
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import contextvars
import select
import threading
import time
import weakref
from collections import deque

from greenlet import getcurrent, greenlet

from ipi.utils.softexit import softexit

__all__ = ["Task", "spawn", "join", "run", "wait", "watch", "unwatch"]

# longest wait for a reply before polling again, so that the time-outs of
# the drivers and the requests for a soft exit are noticed
MAXWAIT = 0.1
# shortest wait that poll() can make: it counts in milliseconds, rounding up
MINPOLL = 1e-3


class Task:
    """A function run by the scheduler, with its result or error.

    Attributes:
        context: The context variables of the code that spawned the task,
            which the task sees.
        worker: The worker that runs the task, once it has started.
        done: Whether the task has finished.
    """

    __slots__ = ("func", "args", "kwargs", "context", "worker", "done")
    __slots__ += ("result", "error")

    def __init__(self, func, args, kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs
        self.context = contextvars.copy_context()
        self.worker = None
        self.done = False
        self.result = None
        self.error = None


class Worker(greenlet):
    """A greenlet that runs tasks, one after the other. When its task ends,
    it joins the idle workers and switches back, until the root gives it a
    new task: starting a greenlet costs ten times more than switching to one.

    Attributes:
        scheduler: The scheduler the worker belongs to.
        task: The task being run, or None.
    """

    def __init__(self, scheduler):
        super().__init__(parent=scheduler.root)
        self.scheduler = scheduler
        self.task = None

    def run(self, task):
        while True:
            # switches out inside the function whenever the task waits
            try:
                task.result = task.func(*task.args, **task.kwargs)
            except BaseException as err:
                # also GreenletExit, thrown in by `_kill`
                task.error = err
            task.done = True
            self.task = None
            self.scheduler.idle.append(self)
            # returns to the parent (the root, or the task that stopped this
            # one); the next task arrives as the value returned by the switch
            task = self.parent.switch()


class Scheduler:
    """The tasks and the forcefields of one thread.

    Attributes:
        root: The root greenlet of the thread, which runs the loop.
        ready: Tasks that can start, and workers that can resume their task,
            in the order in which they run.
        waiting: Workers parked until a condition holds, as (worker,
            condition).
        idle: Workers without a task, to be reused.
        watched: Forcefields with pending requests, as an ordered set of weak
            references, which does not keep forcefields alive.
        polling: True while the forcefields are polled.
    """

    def __init__(self):
        # the greenlet that the thread started in has no parent
        root = getcurrent()
        while root.parent is not None:
            root = root.parent
        self.root = root
        self.ready = deque()
        self.waiting = []
        self.idle = []
        self.watched = weakref.WeakKeyDictionary()
        self.polling = False

    def spawn(self, func, *args, **kwargs):
        """Creates a task that runs `func(*args, **kwargs)`. It starts the
        next time the loop runs, i.e. when the caller waits."""

        task = Task(func, args, kwargs)
        self.ready.append(task)
        return task

    def _resume(self, item):
        """Resumes a worker, or starts a task on a free worker. Returns when
        the worker switches back, because its task waits or has ended."""

        if isinstance(item, Worker):
            item.switch()
            return
        # get a worker from the idle pool, or create a new one.
        # workers are recycled to avoid the creation cost
        worker = self.idle.pop() if self.idle else Worker(self)
        worker.parent = getcurrent()
        # the task sees the context variables of the code that spawned it
        worker.gr_context = item.context
        worker.task = item
        item.worker = worker
        # a new worker gets the task as the argument of `run`, an idle one
        # as the return value of its last switch
        worker.switch(item)

    def join(self, tasks):
        """Waits until the tasks are done, and returns their results. If one
        of them fails, the others are stopped and its error is raised."""

        def finished():
            return all(t.done for t in tasks) or any(t.error is not None for t in tasks)

        try:
            self.wait(finished)
        finally:
            # also stops the tasks if the wait is interrupted (e.g. by an
            # error in a forcefield)
            failed = [t for t in tasks if t.error is not None]
            for task in tasks:
                self._kill(task)
        if failed:
            raise failed[0].error
        return [task.result for task in tasks]

    def run(self, funcs):
        """Runs the functions in tasks, and returns their results."""

        return self.join([self.spawn(f) for f in funcs])

    def watch(self, forcefield):
        """Polls `forcefield` until it has no requests. Called when it
        queues a request."""

        self.watched[forcefield] = None

    def unwatch(self, forcefield):
        """Stops polling `forcefield`."""

        self.watched.pop(forcefield, None)

    def wait(self, condition):
        """Returns when `condition()` is true, running the other tasks and
        polling the forcefields meanwhile.

        A task parks itself and switches to the root, which resumes it once
        the condition holds. The root runs the loop itself."""

        if condition():
            return
        # a wait inside a poll would run the loop that is polling
        if self.polling:
            raise RuntimeError("Forces requested while polling the forcefields")
        # makes a requested soft exit, which stops the forcefields
        softexit.check()
        current = getcurrent()
        if isinstance(current, Worker):
            # another task may make the condition false again before this
            # one resumes (e.g. by taking a resource both wait for)
            while not condition():
                # makes the progress that needs no waiting (e.g. sends the
                # requests of the task to a driver, or evaluates them in
                # process), and runs the other tasks only if it is not enough:
                # a switch costs more than a fast evaluation
                if self._poll(idle=False) and condition():
                    return
                entry = (current, condition)
                self.waiting.append(entry)
                try:
                    # parks the task: execution continues in the root loop,
                    # and comes back here when the loop resumes this worker
                    self.root.switch()
                finally:
                    if entry in self.waiting:
                        self.waiting.remove(entry)
        elif current is self.root:
            self._loop(condition)
        else:
            raise RuntimeError("Waiting for forces from a foreign greenlet")

    def _loop(self, until):
        """Runs the tasks and polls the forcefields until `until()` is true."""

        stuck = False
        while not until():
            if not self.ready:
                self._wake()
            if self.ready:
                # runs a task until it waits again or ends
                self._resume(self.ready.popleft())
                stuck = False
                continue
            # every task waits: makes a requested soft exit, evaluates the
            # incomplete batches, or else sleeps until a driver replies
            softexit.check()
            if self._poll(idle=True) or self._block():
                stuck = False
            elif stuck:
                names = [getattr(f, "name", type(f).__name__) for f in self.watched]
                raise RuntimeError(
                    "Waiting for forces that cannot be computed: no task can run, "
                    f"and the forcefields {names} expect no reply"
                )
            else:
                # checks again, in case something changed since the last poll
                # (e.g. a soft exit that stopped the forcefields)
                stuck = True

    def _wake(self):
        """Moves the parked workers whose condition holds to `ready`."""

        waiting = []
        for worker, condition in self.waiting:
            if condition():
                self.ready.append(worker)
            else:
                waiting.append((worker, condition))
        self.waiting = waiting

    def _poll(self, idle):
        """Polls each watched forcefield once, and returns True if any made
        progress. `idle` tells them that every task waits."""

        progress = False
        self.polling = True
        try:
            # a copy, as the forcefields left without requests are removed
            for forcefield in list(self.watched):
                if forcefield.poll(idle):
                    progress = True
                if not forcefield.requests:
                    del self.watched[forcefield]
        finally:
            self.polling = False
        return progress

    def _block(self):
        """Sleeps until a socket of the forcefields can be read, or for the
        shortest `poll_interval`, and at most MAXWAIT. Returns False, without
        sleeping, if there is nothing to wait for."""

        poller = select.poll()
        nfds = 0
        timeouts = []
        for forcefield in self.watched:
            for fd in forcefield.fds():
                try:
                    poller.register(fd, select.POLLIN)
                    nfds += 1
                except (ValueError, OSError):
                    # closed by a soft exit: the next poll sees it stopped
                    timeouts.append(0.0)
            if forcefield.poll_interval is not None:
                timeouts.append(forcefield.poll_interval)
        if nfds == 0 and not timeouts:
            return False
        timeout = min(timeouts + [MAXWAIT])
        if nfds == 0:
            time.sleep(timeout)
        elif timeout >= MINPOLL:
            poller.poll(timeout * 1e3)
        elif not poller.poll(0):
            # poll() counts in whole milliseconds and rounds up, which would
            # turn a shorter interval (e.g. the latency of MPI drivers) into
            # one millisecond: checks the sockets, and sleeps instead
            time.sleep(timeout)
        return True

    def _kill(self, task):
        """Stops a task, if it is not done, by raising GreenletExit where it
        waits."""

        if task.done:
            return
        worker = task.worker
        if worker is None:
            # a task that never started is discarded without running
            self.ready.remove(task)
            task.done = True
            return
        if worker in self.ready:
            self.ready.remove(worker)
        self.waiting = [(w, c) for w, c in self.waiting if w is not worker]
        # the worker records the exit as the error of its task, and returns
        # here, also when the task is stopped by another task
        worker.parent = getcurrent()
        worker.throw()


_local = threading.local()


def _scheduler():
    """Returns the scheduler of the current thread."""

    try:
        return _local.scheduler
    except AttributeError:
        _local.scheduler = Scheduler()
        return _local.scheduler


def spawn(func, *args, **kwargs):
    """Creates a task that runs `func(*args, **kwargs)`."""

    return _scheduler().spawn(func, *args, **kwargs)


def join(tasks):
    """Waits until the tasks are done, and returns their results. If one of
    them fails, the others are stopped and its error is raised."""

    return _scheduler().join(tasks)


def run(funcs):
    """Runs the functions in tasks, and returns their results."""

    return _scheduler().run(funcs)


def wait(condition):
    """Returns when `condition()` is true, running the other tasks and
    polling the forcefields meanwhile."""

    _scheduler().wait(condition)


def watch(forcefield):
    """Polls `forcefield` until it has no requests."""

    _scheduler().watch(forcefield)


def unwatch(forcefield):
    """Stops polling `forcefield`."""

    _scheduler().unwatch(forcefield)
