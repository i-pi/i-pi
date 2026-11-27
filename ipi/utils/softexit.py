"""Soft exits: ending a simulation half way through, in a state from which it
can be restarted.

A soft exit is requested by a kill signal (SIGINT or SIGTERM), by a file
named EXIT, when the wall-clock time runs out, or by the code (`request`). A
request only records a message: the main thread makes the exit at its next
`check`, which it calls between the steps and whenever it waits for forces.
Making the exit runs the cleanup functions, which write RESTART, close the
outputs and stop the forcefields, so that the waits for forces end. During a
run, RESTART holds the state that was last stored (`rollback`), and the
outputs written since then are truncated (`stored` starts a new record of
their positions), so that a restart writes each of them once. A second
signal interrupts the main thread at once.

A watchdog thread deals with a main thread that does not reach `check`, e.g.
because an in-process force evaluation hangs in native code, where Python
does not even run its signal handlers. It learns of the signals through
`signal.set_wakeup_fd`, and checks the EXIT file and the wall-clock time
itself. When an exit is pending and the main thread has not checked for
GRACE seconds, the watchdog writes the last stored checkpoint and ends the
process. It does so only with `hard_exit`, i.e. when i-PI runs as a program.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.


import atexit
import os
import select
import signal
import sys
import threading
import time

from ipi.utils.messages import verbosity, warning

__all__ = ["Softexit", "softexit"]


# time between checks for an EXIT file and for the wall-clock time
SOFTEXITLATENCY = 1.0
# time after which the watchdog makes a pending exit, if the main thread has
# not checked for requests meanwhile
GRACE = 10.0


class Softexit(object):
    """Requests to stop a simulation half way through, and the functions that
    make it restartable.

    Attributes:
       flist: The cleanup functions, as (function, args, kwargs), run once.
       emergency: The functions that the watchdog runs when the main thread
          does not respond. They must not use what the main thread may be
          using, e.g. they write a checkpoint stored beforehand.
       requested: The message of a pending request, or None.
       triggered: Whether the exit has started.
       exiting: True while the cleanup functions run.
       rollback: Whether the exit goes back to the state that was last
          stored, also truncating the outputs written since then.
       nstored: The number of times the state has been stored.
       hard_exit: If True, the process is ended with os._exit() once the
          cleanup is complete, skipping the finalization of the interpreter.
       timeout: The wall-clock time at which to exit, if positive.
    """

    def __init__(self):
        """Initializes SoftExit."""

        self.flist = []
        self.emergency = []
        self.requested = None
        self.triggered = False
        self.exiting = False
        self.hard_exit = False
        self.timeout = -1.0
        self.rollback = False
        self.nstored = 0
        self._started = False
        self._signalled = False
        self._kill = {}
        # the last check of the main thread, and the next one of the files
        self._beat = time.monotonic()
        self._next_check = 0.0
        # taken by whoever makes the exit: the main thread or the watchdog
        self._lock = threading.Lock()
        # held by the main thread while it writes to the outputs or stores the
        # state, and taken for good by the watchdog before it truncates them,
        # so that nothing reaches the files after the truncation
        self.io_lock = threading.RLock()
        self._watchdog = None

    def register_function(self, func, *args, **kwargs):
        """Adds a cleanup function to flist."""

        self.flist.append((func, args, kwargs))

    def register_emergency(self, func, *args, **kwargs):
        """Adds a function for the watchdog to run if the main thread hangs."""

        self.emergency.append((func, args, kwargs))

    def request(self, message=""):
        """Asks for a soft exit, which the main thread makes at its next
        check. Can be called from any thread, and from signal handlers."""

        if self.requested is None and not self.triggered:
            self.requested = message

    def check(self):
        """Makes the exit that has been requested. Called by the main thread
        between the steps and whenever it waits for forces, which also tells
        the watchdog that it responds."""

        self._beat = time.monotonic()
        self.poll()
        if self.requested is not None and not self.exiting:
            self.cleanup(message=self.requested)

    def poll(self, force=False):
        """Requests an exit if there is an EXIT file or the wall-clock time
        has run out, looking at most once per SOFTEXITLATENCY unless `force`."""

        now = time.monotonic()
        if self._started and (force or now > self._next_check):
            self._next_check = now + SOFTEXITLATENCY
            message = self._expired()
            if message is not None:
                self.request(message)

    def stored(self):
        """Records that the state to go back to has been stored: the outputs
        record their positions again, at their next write."""

        self.nstored += 1

    def trigger(self, *, status="restartable", message=""):
        """Halts the simulation.

        Prints out a warning message, then runs all the exit functions in flist
        before terminating the simulation.

        Args:
           status: which kind of stop it is: simulation restartable as is,
                   successful finish or aborted because of some problem.
           message: The message to output to standard output.
        """

        self.cleanup(status=status, message=message)
        self.kill()

    def cleanup(self, *, status="restartable", message=""):
        """Runs the cleanup functions, once, without exiting the process."""

        # reports the reason of a pending request
        message = message or self.requested or ""
        print(" @softexit.trigger:  SOFTEXIT CALLED", message)
        if not self.triggered:
            self.triggered = True

            if status == "restartable":
                message += " Restartable as is: YES."
            elif status == "success":
                message += " I-PI reports success. Restartable as is: NO."
            elif status == "bad":
                message += " I-PI reports a problem. Restartable as is: NO."
            else:
                raise ValueError("Unknown option for softexit status.")

            warning(
                "Soft exit has been requested with message: '"
                + message
                + "'. Cleaning up.",
                verbosity.low,
            )
        self.requested = None

        if self.exiting:
            return
        # waits, if the watchdog is ending the process
        with self._lock:
            self.exiting = True
            try:
                # calls all the registered emergency softexit procedures
                while self.flist:
                    f, a, ka = self.flist.pop(0)
                    # the other functions run also if one fails (e.g. writing
                    # the checkpoint on a full disk), so that the drivers exit
                    try:
                        f(*a, **ka)
                    except Exception as err:
                        print("Error running emergency softexit, ", err)
            finally:
                self.exiting = False

    def kill(self):
        """Terminates the process, or unwinds the main thread.

        With hard_exit set, and once the cleanup functions have run, the
        process is terminated. Outputs and RESTART have been written at that
        point, so nothing is lost by skipping the finalization of the
        interpreter, which can crash when extension libraries (e.g. libtorch
        worker threads) release thread-local state.
        """

        if self.hard_exit and self.triggered and not (self.flist or self.exiting):
            # MPI is otherwise finalized when the interpreter terminates
            mpi = sys.modules.get("mpi4py.MPI")
            if mpi is not None and mpi.Is_initialized() and not mpi.Is_finalized():
                mpi.Finalize()
            sys.stdout.flush()
            sys.stderr.flush()
            os._exit(0)
        sys.exit()

    def reset(self):
        """Stops the watchdog, restores the intercepted signals, and forgets
        the requests and the functions."""

        self._stop_watchdog()
        for signum, handler in self._kill.items():
            signal.signal(signum, handler)
        self.__init__()

    def start(self, timeout=0.0):
        """Starts watching for exit requests: kill signals, the EXIT file, and
        the wall-clock time.

        Args:
           timeout: Number of seconds after which the exit is requested, if
              positive.
        """

        self.timeout = -1.0
        if timeout > 0.0:
            self.timeout = time.time() + timeout
        self._started = True
        self._beat = time.monotonic()
        # an EXIT file that is there from the start is noticed at the end of
        # a short first step, which the run then completes
        self._next_check = self._beat + SOFTEXITLATENCY
        # signals are handled by the main thread
        if threading.current_thread() is not threading.main_thread():
            return
        for signum in (signal.SIGINT, signal.SIGTERM):
            handler = signal.signal(signum, self._kill_handler)
            # keeps the original handler when started more than once
            self._kill.setdefault(signum, handler)
        if self._watchdog is None:
            self._start_watchdog()

    def _expired(self):
        """Returns the message of an exit asked by an EXIT file or by the
        wall-clock time, or None."""

        if os.path.exists("EXIT"):
            return " @SOFTEXIT: EXIT file detected."
        if 0.0 < self.timeout < time.time():
            return " @SOFTEXIT: Maximum wallclock time elapsed."
        return None

    def _kill_handler(self, signum, frame):
        """Requests a soft exit on the first SIGINT or SIGTERM, and falls
        back to the original handler on the second one.

        Python runs this in the main thread, on top of whatever it was doing,
        so the exit is left to the next check.
        """

        warning(
            " @SOFTEXIT:   Kill signal. Trying to make a clean exit.", verbosity.low
        )

        if not self._signalled:
            self._signalled = True
            self.request(" @SOFTEXIT: Kill signal received")
        elif callable(self._kill.get(signum)):
            self._kill[signum](signum, frame)
        else:
            self.kill()

    def _start_watchdog(self):
        """Starts the watchdog thread, and routes the signals to its pipe."""

        read, write = os.pipe()
        os.set_blocking(write, False)
        # the signal handler at the C level writes the signal number to the
        # pipe, even when the main thread is stuck in native code
        previous = signal.set_wakeup_fd(write, warn_on_full_buffer=False)
        stop = threading.Event()
        thread = threading.Thread(
            target=self._watch, args=(read, stop), name="softexit", daemon=True
        )
        self._watchdog = (thread, read, write, stop, previous)
        thread.start()

    def _stop_watchdog(self):
        """Stops the watchdog thread, and restores the signal wake-up."""

        if self._watchdog is None:
            return
        thread, read, write, stop, previous = self._watchdog
        stop.set()
        os.write(write, b"\0")  # wakes the thread up
        thread.join()
        if threading.current_thread() is threading.main_thread():
            signal.set_wakeup_fd(previous)
        os.close(read)
        os.close(write)
        self._watchdog = None

    def _watch(self, read, stop):
        """Waits for signals, checks the EXIT file and the wall-clock time,
        and makes a pending exit if the main thread does not respond."""

        while True:
            ready, _, _ = select.select([read], [], [], SOFTEXITLATENCY)
            if stop.is_set():
                return
            if ready:
                signums = os.read(read, 512)
                if signal.SIGINT in signums or signal.SIGTERM in signums:
                    self.request(" @SOFTEXIT: Kill signal received")
            if self._started:
                message = self._expired()
                if message is not None:
                    self.request(message)
            pending = self.requested is not None or (self.triggered and self.flist)
            if (
                pending
                and self.hard_exit
                and not self.exiting
                and time.monotonic() - self._beat > GRACE
                and self._lock.acquire(blocking=False)
            ):
                self._emergency()

    def _emergency(self):
        """Runs the emergency functions, and ends the process."""

        sys.stderr.write(
            " @SOFTEXIT: The main thread has not responded for %g s: writing the "
            "stored checkpoint, and exiting.\n" % GRACE
        )

        # waits for a write of the main thread in progress, then keeps the
        # outputs to itself: the main thread blocks at its next write
        if not self.io_lock.acquire(timeout=GRACE):
            sys.stderr.write(" @SOFTEXIT: The outputs are busy: exiting as is.\n")
            sys.stderr.flush()
            os._exit(1)

        def run():
            for f, a, ka in self.emergency:
                try:
                    f(*a, **ka)
                except Exception as err:
                    sys.stderr.write("Error running emergency softexit, %s\n" % err)

        # the checkpoint copies arrays from the device, which can be stuck
        # with the main thread: the process ends also if it does not finish
        writer = threading.Thread(target=run, daemon=True)
        writer.start()
        writer.join(GRACE)
        if writer.is_alive():
            sys.stderr.write(" @SOFTEXIT: The checkpoint could not be written.\n")
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(1)

    def _finish(self):
        """Completes a soft exit before the interpreter terminates."""

        if self.requested is not None or (self.triggered and self.flist):
            self.cleanup(message=self.requested or "")
        if self.triggered and self.hard_exit:
            self.kill()

    def _after_fork(self):
        """Detaches a forked child (e.g. a driver) from the soft exit of the
        parent, whose watchdog would otherwise receive its signals."""

        if self._watchdog is not None:
            signal.set_wakeup_fd(-1)
            for fd in self._watchdog[1:3]:
                os.close(fd)
            self._watchdog = None
        for signum, handler in self._kill.items():
            signal.signal(signum, handler)
        self._kill = {}


softexit = Softexit()
atexit.register(softexit._finish)
os.register_at_fork(after_in_child=softexit._after_fork)
