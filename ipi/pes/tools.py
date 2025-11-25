import os
import time
from datetime import datetime
from typing import Protocol, Callable, Any, TypeVar
from typing_extensions import ParamSpec
import functools


# --------------------------------------- #
class Timer:
    def __init__(self, enabled=True, file: str = None):
        """
        Independent timer manager.

        Each Timer instance keeps:
        - its own stack for nested timing contexts
        - its own records for reporting

        Usage:
            t = Timer()
            with t.section("outer"):
                ...
                with t.section("inner"):
                    ...
            t.report()
        """
        self._stack = []
        self._records = []
        self.enabled = enabled
        if file is not None and file == "":
            file = None
        self.file = file

    def section(self, name, log=True):
        """Create a context manager for a timed section."""
        if not self.enabled:
            return _DummySection()
        return _TimerSection(self, name, log)

    def report(self):
        """Print or save all recorded timings for this Timer instance."""
        if not self.enabled:
            return

        lines = []
        for t, level, name, elapsed in sorted(self._records):
            timestamp = datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M:%S")
            indent = "    " * level
            slurmid = int(os.environ.get("SLURM_LOCALID", "0"))
            lines.append(
                f"{timestamp} [id:{slurmid:3d}] {indent}{name}: {elapsed:.4f} s"
            )

        if self.file:
            with open(self.file, "a", encoding="utf-8") as f:
                for line in lines:
                    f.write(line + "\n")
                    f.flush()  # optional: flush to disk immediately
        else:
            for line in lines:
                print(line, flush=True)

        self.reset()

    def reset(self):
        """Clear all records for this Timer instance."""
        self._records.clear()


class _TimerSection:
    """Context manager for a single timing block inside a Timer."""

    def __init__(self, timer, name, log):
        self.timer = timer
        self.name = name
        self.log = log

    def __enter__(self):
        self.start = time.time()
        self.level = len(self.timer._stack)
        self.timer._stack.append(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = time.time()
        self.elapsed = self.end - self.start
        self.timer._stack.pop()
        if self.log:
            self.timer._records.append(
                (self.start, self.level, self.name, self.elapsed)
            )


class _DummySection:
    """No-op context manager when timing is disabled."""

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


# --------------------------------------- #
# Define a protocol that requires a `logger` attribute
class HasLogger(Protocol):
    logger: "Timer"  # replace "Timer" with your actual Timer class name


P = ParamSpec("P")
R = TypeVar("R")


def timeit(
    name: str, report: bool = False
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """
    Decorator to measure the execution time of a class method using `self.logger.section`.

    Parameters
    ----------
    name : str
        The name of the timed section; will be shown in the logger.
    report : bool, default=False
        If True, calls `self.logger.report()` after the method finishes.
    """

    def decorator(func: Callable[P, R]) -> Callable[P, R]:
        @functools.wraps(func)
        def wrapper(self: HasLogger, *args, **kwargs) -> Any:
            with self.logger.section(name):
                result = func(self, *args, **kwargs)
            if report:
                self.logger.report()
            return result

        return wrapper

    return decorator
