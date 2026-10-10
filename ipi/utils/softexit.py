"""Classes to deal with calls for a soft exit."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import atexit
import sys
import os
import time
import threading
import signal

from ipi.utils.messages import verbosity, warning

__all__ = ["Softexit", "softexit"]


SOFTEXITLATENCY = 1.0  # seconds to sleep between checking for soft exit


class Softexit(object):
    """Class to deal with stopping a simulation half way through.

    Provides a mechanism to end a simulation from any thread that has
    been properly registered, and to call a series of "emergency" functions
    to try as hard as possible to produce a restartable snapshot of
    the simulation.
    Also, provides a loop to check for soft-exit requests and
    trigger termination when necessary.

    Attributes:
       flist: A list of callback functions used to clean up and exit gracefully.
       tlist: A list of threads registered for monitoring
       lock: Held while running the callback functions. A thread can also hold
          it to delay them until it has reached a consistent state: it must
          then check if a soft exit was triggered, and stop if this is the case.
       hard_exit: If True, the process is ended with os._exit() once the
          cleanup is complete, skipping the finalization of the interpreter.
    """

    def __init__(self):
        """Initializes SoftExit."""

        self.flist = []
        self.tlist = []
        self._kill = {}
        self._thread = None
        self.triggered = False
        self.exiting = False
        self._doloop = [False]
        self._killed = False
        self.lock = threading.RLock()
        self.hard_exit = False

    def register_function(self, func, *args, **kwargs):
        """Adds another function to flist.

        Args:
           func: The function to be added to flist.
        """

        self.flist.append((func, args, kwargs))

    def register_thread(self, thread, loop_control=None):
        """Adds a thread to the monitored list.

        Args:
           thread: The thread to be monitored.
           loop_control: the variable that causes the thread to terminate.
        """

        self.tlist.append((thread, loop_control))

    def trigger(self, status="restartable", message=""):
        """Halts the simulation.

        Prints out a warning message, then runs all the exit functions in flist
        before terminating the simulation.

        Args:
           status: which kind of stop it is: simulation restartable as is,
                   successful finish or aborted because of some problem.
           message: The message to output to standard output.
        """

        self.cleanup(status, message)
        self.kill()

    def cleanup(self, status="restartable", message=""):
        """Runs registered cleanup functions without exiting the process.

        The functions are called only once. If another thread holds the lock,
        they are left pending, and run the next time this is called.
        """

        print(
            " @softexit.trigger:  SOFTEXIT CALLED FROM THREAD",
            threading.currentThread(),
            message,
        )
        if not self.triggered:  # avoid double calls from different threads
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

        if not self.lock.acquire(blocking=False):
            return
        try:
            self.exiting = True
            # calls all the registered emergency softexit procedures
            while self.flist:
                f, a, ka = self.flist.pop(0)
                try:
                    f(*a, **ka)
                except RuntimeError as err:
                    print("Error running emergency softexit, ", err)
        finally:
            self.exiting = False  # emergency is over, signal we can be relaxed
            self.lock.release()

        for t, dl in self.tlist:  # set thread exit flag
            dl[0] = False

        # wait for all (other) threads to finish
        for t, dl in self.tlist:
            if not (
                threading.currentThread() is self._thread
                or threading.currentThread() is t
            ):
                t.join()

    def kill(self):
        """Terminates the current thread/process.

        With hard_exit set, and once the cleanup functions have run, the whole
        process is terminated from whichever thread gets here. Outputs and
        RESTART have been written at that point, so nothing is lost by skipping
        the finalization of the interpreter, which can crash when extension
        libraries (e.g. libtorch worker threads) release thread-local state.
        This also stops the main thread from carrying on with a step that uses
        force fields that have been shut down.
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
        """Resets the soft-exit state and restores intercepted signals."""

        for signum, handler in self._kill.items():
            signal.signal(signum, handler)
        self.__init__()

    def start(self, timeout=0.0):
        """Starts the softexit monitoring loop.

        Args:
           timeout: Number of seconds to wait before softexit is triggered.
        """

        self._main = threading.currentThread()
        self.timeout = -1.0
        if timeout > 0.0:
            self.timeout = time.time() + timeout

        self._thread = threading.Thread(target=self._softexit_monitor, name="softexit")
        self._thread.daemon = True
        self._doloop[0] = True
        self._kill[signal.SIGINT] = signal.signal(signal.SIGINT, self._kill_handler)
        self._kill[signal.SIGTERM] = signal.signal(signal.SIGTERM, self._kill_handler)
        self._thread.start()
        self.register_thread(self._thread, self._doloop)

    def _kill_handler(self, signal, frame):
        """Deals with handling a kill call gracefully.

        Intercepts kill signals to trigger softexit.
        Called when signals SIG_INT and SIG_TERM are received.

        Args:
           signal: An integer giving the signal number of the signal received
              from the socket.
           frame: Current stack frame.
        """

        warning(
            " @SOFTEXIT:   Kill signal. Trying to make a clean exit.", verbosity.low
        )

        if not self._killed:
            # the handler runs on top of whatever the main thread was doing,
            # so the soft exit is left to the monitoring thread
            self._killed = True
        elif callable(self._kill.get(signal)):
            # a second signal falls back to the original handler
            self._kill[signal](signal, frame)
        else:
            self.kill()

    def _softexit_monitor(self):
        """Keeps checking for soft exit conditions."""

        while self._doloop[0]:
            time.sleep(SOFTEXITLATENCY)
            if self._killed:
                self.trigger(
                    status="restartable", message=" @SOFTEXIT: Kill signal received"
                )
                break

            if os.path.exists("EXIT"):
                self.trigger(
                    status="restartable", message=" @SOFTEXIT: EXIT file detected."
                )
                break

            if self.timeout > 0 and self.timeout < time.time():
                self.trigger(
                    status="restartable",
                    message=" @SOFTEXIT: Maximum wallclock time elapsed.",
                )
                break

    def _finish(self):
        """Completes a soft exit before the interpreter terminates."""

        if self.triggered:
            with self.lock:  # waits for a cleanup in progress
                if self.flist:
                    self.cleanup()
            if self.hard_exit:
                self.kill()


softexit = Softexit()
atexit.register(softexit._finish)
