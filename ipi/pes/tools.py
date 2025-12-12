import os
import json
import time
import numpy as np
from datetime import datetime
from typing import Protocol, Callable, Any, TypeVar, Dict, Tuple, List, Union
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
            timestamp = datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
            indent = "    " * level
            slurmid = int(os.environ.get("SLURM_LOCALID", "-1"))
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


# --------------------------------------- #
class JSONLogger:

    def __init__(self, file: str = None):
        self.file = file
        self.enabled = file is not None

    def _serialize(self, obj: Any) -> Any:
        """Recursively convert objects into JSON-serializable types."""

        # NumPy arrays → lists
        if isinstance(obj, np.ndarray):
            return obj.tolist()

        # NumPy scalar types → Python scalars
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()

        # Standard JSON scalars
        if isinstance(obj, (int, float, str, bool)) or obj is None:
            return obj

        # Dicts → recursively serialize values
        if isinstance(obj, dict):
            return {k: self._serialize(v) for k, v in obj.items()}

        # Lists / tuples / sets → recursively process each element
        if isinstance(obj, (list, tuple, set)):
            return [self._serialize(x) for x in obj]

        # Fallback: store as string if unrecognized
        return str(obj)

    def save(self, results: Dict[str, Any], file: str = None):
        if not self.enabled:
            return

        safe_data = self._serialize(results)

        if file is None:
            file = self.file
        with open(file, "a") as f:
            f.write(json.dumps(safe_data, indent=4))
            f.write("\n")


# --------------------------------------- #
Parent = Dict[str, Union[float, np.ndarray]]


class StructureResults(dict):
    """
    Stores the results produced by one model for one structure.

    Parameters
    ----------
    natoms : int
        Number of atoms in this structure. Used to expand dynamic
        shape entries equal to the string "natoms".
    shapes : Dict[str, Tuple]
        Desired output shapes for each property. Tuples may contain the
        literal "natoms", which is replaced by the actual integer.
    """

    def __init__(self, natoms: int, shapes: Dict[str, Tuple]):
        self.natoms = natoms
        # Expand "natoms" once at initialization
        self.shapes = {k: self._expand(shape) for k, shape in shapes.items()}

    def _expand(self, shape: Tuple) -> Tuple:
        """Replace 'natoms' in shape by the actual integer."""
        if not isinstance(shape, tuple):
            raise TypeError(f"Expected tuple, got {type(shape)}")
        return tuple(self.natoms if d == "natoms" else d for d in shape)

    def store(self, key: str, value: Union[float, np.ndarray]):
        """Store a property value, reshaping it to the expected shape."""
        if key not in self.shapes:
            raise KeyError(f"Unknown property '{key}'")

        shape = self.shapes[key]
        if shape == ():  # scalar
            try:
                self[key] = float(value)
            except Exception:
                raise ValueError(f"Value for '{key}' must be convertible to float")
        else:
            try:
                self[key] = np.asarray(value).reshape(shape)
            except Exception:
                raise ValueError(
                    f"Value for '{key}' has wrong shape: expected {shape}, got {np.asarray(value).shape}"
                )

    def __add__(self, other: "StructureResults") -> "StructureResults":
        out = StructureResults(self.natoms, self.shapes)
        out.update({k: self[k] + other[k] for k in self})
        return out

    def __truediv__(self, divisor: float) -> "StructureResults":
        out = StructureResults(self.natoms, self.shapes)
        out.update({k: self[k] / divisor for k in self})
        return out

    def as_dict(self, copy_arrays: bool = True) -> Parent:
        """Return stored properties as a plain dictionary, optionally copying arrays."""
        return {
            k: (v.copy() if copy_arrays and isinstance(v, np.ndarray) else v)
            for k, v in self.items()
        }


class ModelResults:
    """
    Handle results returned by one model for multiple structures.
    """

    def __init__(self, shapes: Dict[str, Tuple]):
        self._shapes = shapes
        self._results: List[StructureResults] = []

    def store(self, natoms: List[int], results: Dict[str, Any]):
        """
        Store results of one model over multiple structures.
        """
        ptr = np.cumsum([0] + natoms)
        new_structs = [StructureResults(n, self._shapes) for n in natoms]

        for key, value in results.items():
            if key not in self._shapes:
                raise ValueError(f"Unknown property '{key}'")
            if "natoms" in self._shapes[key]:
                value = np.split(value, ptr[1:], axis=0)[:-1]

            for i, s in enumerate(new_structs):
                s.store(key, value[i])

        self._results.extend(new_structs)

    def __len__(self) -> int:
        return len(self._results)

    def natoms(self) -> List[int]:
        return [s.natoms for s in self._results]

    def shapes(self) -> List[Dict[str, Tuple]]:
        return [s.shapes for s in self._results]

    @staticmethod
    def mean(models: List["ModelResults"]) -> "ModelResults":
        if not models:
            raise ValueError("Cannot compute mean of empty list")

        ref_natoms = models[0].natoms()
        ref_shapes = models[0].shapes()
        shapes = models[0]._shapes
        n_structures = len(ref_natoms)
        n_models = len(models)

        # Validate consistency
        for m in models[1:]:
            if m.natoms() != ref_natoms or m.shapes() != ref_shapes:
                raise ValueError(
                    "All ModelResults must have the same natoms and shapes per structure"
                )

        out = ModelResults(shapes)
        for s in range(n_structures):
            # Initialize with a zeroed StructureResults
            summed = StructureResults(ref_natoms[s], shapes)
            for key in models[0]._results[s].keys():
                summed[key] = np.zeros(summed.shapes[key])

            # Sum over models
            for m in models:
                summed += m._results[s]

            out._results.append(summed / n_models)

        return out

    def __getitem__(self, i: int) -> Parent:
        return self._results[i].as_dict()
