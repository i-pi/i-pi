"""Classes to deal with calls for a soft exit."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


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

        print(
            " @softexit.trigger:  SOFTEXIT CALLED FROM THREAD",
            threading.currentThread(),
            message,
        )
        if not self.triggered:  # avoid double calls from different threads
            self.exiting = True
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

            # calls all the registered emergency softexit procedures
            for f, a, ka in self.flist:
                try:
                    f(*a, **ka)
                except RuntimeError as err:
                    print("Error running emergency softexit, ", err)
                    pass

            self.exiting = False  # emergency is over, signal we can be relaxed

            for t, dl in self.tlist:  # set thread exit flag
                dl[0] = False

        # wait for all (other) threads to finish
        for t, dl in self.tlist:
            if not (
                threading.currentThread() is self._thread
                or threading.currentThread() is t
            ):
                t.join()

        sys.exit()

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

        self.trigger(status="restartable", message=" @SOFTEXIT: Kill signal received")

        try:
            self.__del__()
        except AttributeError:
            pass
        if signal in self._kill:
            self._kill[signal](signal, frame)

    def _softexit_monitor(self):
        """Keeps checking for soft exit conditions."""

        while self._doloop[0]:
            time.sleep(SOFTEXITLATENCY)
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


softexit = Softexit()
