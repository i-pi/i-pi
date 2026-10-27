"""Recording of a function that operates on depend objects as CUDA graphs, and
replay of the graphs instead of the function.

With the torch backend on a CUDA device, a step of a simulation launches many
small kernels, and its cost is set by the time the host takes to run the
Python code and to launch them, rather than by the device. A CUDA graph
records the kernels launched by the function, and replays all of them at once
without running the Python code again.

The function runs eagerly for the first calls, which fill the caches of the
depend objects, and is recorded at the next one. The recording can be
interrupted to run some parts of the function eagerly, e.g. those that
communicate with external codes or synchronize with the device: the function
is then recorded in segments, separated by eager functions that run between
their replays. Code that is part of a recorded function marks these points
with `split`, which lets the recording interrupt it (see `GraphReplay.split`),
and `eager`, which calls a function eagerly also when the graphs are replayed,
copying its results into those read by the graphs (see `ipi.utils.jit`, which
also has the recording of the operations on the depend objects).

As the Python code does not run during a replay, the depend objects are kept
consistent by the replay:

- before replaying a segment, the depend objects that it reads, as they were
  before being modified by the segment, are brought up to date. Scalars that
  have been replaced by new objects are copied into the ones read by the
  graph, and host values (such as the time step), which are fixed in the
  graph, must be unchanged;
- after replaying a segment, the taint status of the depend objects and the
  objects held by scalars are restored as they were at the end of the
  recorded segment.

A `JitError` is raised if a host value read by a graph changes, or if the
function does something that cannot be recorded.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import threading
import warnings

import numpy as np

from ipi.utils import depend
from ipi.utils.jit import WARMUP_CALLS, JitError, _Recorder
from ipi.utils.messages import info, verbosity

__all__ = ["GraphReplay"]


class GraphReplay:
    """Runs a function eagerly for the first calls, records the next call as
    CUDA graphs, and replays them at the following calls. The arguments of the
    calls are only passed to the eager calls and to the recorded one.

    Attributes:
        func: The function.
        generator: The torch generator used by the function, if any.
        segments: The recorded segments, or None before the recording.
    """

    def __init__(self, func, generator=None):
        import torch

        self.torch = torch
        self.func = func
        self.generator = generator
        self.stream = torch.cuda.Stream()
        self.pool = torch.cuda.graph_pool_handle()
        self.segments = None
        self.ncalls = 0
        # the threads that have run eager calls
        self.threads = set()

    def __call__(self, *args):
        if self.segments is not None:
            for segment in self.segments:
                segment.replay()
        elif self.ncalls < WARMUP_CALLS or threading.get_ident() not in self.threads:
            # also initializes the libraries on the stream used by the recording,
            # and in this thread (e.g. cuBLAS handles are created per thread)
            default = self.torch.cuda.current_stream()
            self.stream.wait_stream(default)
            with self.torch.cuda.stream(self.stream):
                self.func(*args)
            default.wait_stream(self.stream)
            self.ncalls += 1
            self.threads.add(threading.get_ident())
        else:
            self.record(*args)

    def record(self, *args):
        """Calls the function, recording its kernels as CUDA graphs."""

        self.default = self.torch.cuda.current_stream()
        self.segments = []
        self.recorder = _GraphRecorder(self)
        depend.set_recorder(self.recorder)
        try:
            self.begin()
            self.func(*args)
            self.end(None)
        except Exception as err:
            try:
                self.graph.capture_end()
            except Exception:
                pass
            raise JitError(
                "The step cannot be replayed as CUDA graphs: run i-PI without --jit"
            ) from err
        finally:
            depend.set_recorder(None)
            self.torch.cuda.set_stream(self.default)
        info(
            f" @cudagraphs: step recorded in {len(self.segments)} CUDA graphs",
            verbosity.low,
        )

    def begin(self):
        """Starts recording a segment."""

        self.graph = self.torch.cuda.CUDAGraph()
        if self.generator is not None:
            self.graph.register_generator_state(self.generator)
        self.stream.wait_stream(self.default)
        self.torch.cuda.set_stream(self.stream)
        self.graph.capture_begin(pool=self.pool, capture_error_mode="thread_local")
        self.recorder.reset()

    def end(self, after):
        """Stops recording a segment and replays it. The segment is followed
        by the eager function `after`, if not None."""

        with warnings.catch_warnings():
            # segments between two interruptions may be empty
            warnings.filterwarnings("ignore", message=".*CUDA Graph is empty.*")
            self.graph.capture_end()
        self.torch.cuda.set_stream(self.default)
        self.segments.append(self.recorder.segment(self.graph, after))
        self.graph.replay()

    def interrupt(self, after):
        """Interrupts the recording to run the function `after` eagerly, here
        and after each replay of the segment recorded so far."""

        self.end(after)
        depend.set_recorder(None)
        after()
        depend.set_recorder(self.recorder)
        self.begin()

    def split(self, obj):
        """Called by `split(obj)` while recording, to interrupt the recording
        if the point described by `obj` must run eagerly: does nothing, and is
        overridden by the classes that replay specific functions."""

    def call(self, func, *args):
        """Interrupts the recording to call `func(*args)` eagerly, and returns
        the outputs."""

        call = _Call(func, args)
        self.interrupt(call)
        self.recorder.add_outputs(call.outputs)
        return call.outputs


class _Call:
    """A function called eagerly after a segment. The graphs read the outputs
    of the call made while recording, which are lists: the tensors returned
    by the following calls are copied into them, and the other items replace
    theirs."""

    def __init__(self, func, args):
        import torch

        self.copy = torch._foreach_copy_
        self.func = func
        self.args = args
        self.outputs = None

    def __call__(self):
        outputs = self.func(*self.args)
        if self.outputs is None:
            self.outputs = outputs
            return
        targets, sources = [], []
        _match(self.outputs, outputs, targets, sources)
        if targets:
            self.copy(targets, sources)


def _match(recorded, new, targets, sources):
    """Collects the tensors of the list `new` (and of the lists it contains)
    to be copied into those of `recorded`, and puts its other items in
    `recorded`."""

    if not isinstance(new, list) or len(new) != len(recorded):
        raise JitError("The results of an eager call have changed structure")
    for i, (a, b) in enumerate(zip(recorded, new)):
        if isinstance(a, list):
            _match(a, b, targets, sources)
        elif _is_cuda_tensor(a):
            if not _is_cuda_tensor(b) or a.shape != b.shape or a.dtype != b.dtype:
                raise JitError("The results of an eager call have changed shape")
            targets.append(a)
            sources.append(b)
        else:
            recorded[i] = b


class _Segment:
    """A recorded segment, the state of the depend objects that it reads and
    modifies, and the eager function that follows it."""

    def __init__(self, graph, inputs, statuses, synchros, values, after):
        self.graph = graph
        self.inputs = inputs
        self.statuses = statuses
        self.synchros = synchros
        self.values = values
        self.after = after

    def replay(self):
        for dep, value, copy in self.inputs:
            _stage(dep, value, copy)
        self.graph.replay()
        for status, tainted in self.statuses:
            status.tainted = tainted
        for synchro, manual in self.synchros:
            synchro.manual = manual
        for dep, value in self.values:
            dep._value = value
        if self.after is not None:
            self.after()


def _stage(dep, recorded, copy):
    """Brings `dep` up to date, and makes its value available to a graph that
    reads the object `recorded` (with a copy of its host value `copy`)."""

    current = dep.get()
    if isinstance(dep, depend.depend_array):
        if dep._wvalue is not recorded:
            raise _changed(dep, "has been replaced")
        if copy is not None and not np.array_equal(current, copy):
            raise _changed(dep, "has changed")
    elif current is not recorded:
        if _is_cuda_tensor(recorded):
            if not _is_cuda_tensor(current):
                recorded.fill_(float(current))
            elif current.shape == recorded.shape:
                recorded.copy_(current)
            else:
                raise _changed(dep, "has changed shape")
        elif not _equal(current, recorded):
            raise _changed(dep, "has changed")
    elif copy is not None and not np.array_equal(current, copy):
        raise _changed(dep, "has changed")


def _changed(dep, what):
    return JitError(
        f"The value of <{dep._name}>, which is fixed in the CUDA graphs of the "
        f"step, {what}: the parameters of the simulation that are fixed in the "
        "graphs (e.g. the time step) cannot change when running with --jit"
    )


def _is_cuda_tensor(value):
    return getattr(value, "is_cuda", False)


def _equal(a, b):
    try:
        return bool(np.all(np.asarray(a) == np.asarray(b)))
    except Exception:
        return False


class _GraphRecorder(_Recorder):
    """Recorder of a function replayed as CUDA graphs, which returns the
    recorded segments."""

    def segment(self, graph, after):
        """Returns the recorded segment, with the state of the depend objects
        at its end."""

        if self.host_writes:
            raise JitError(f"The step modifies host arrays: {sorted(self.host_writes)}")
        for dep in self.values:
            if not _is_cuda_tensor(dep._value) and id(dep._value) not in self.outputs:
                raise JitError(f"The step modifies the host value <{dep._name}>")

        statuses, synchros = self.final_state()
        return _Segment(
            graph,
            list(self.inputs.values()),
            statuses,
            synchros,
            [(dep, dep._value) for dep in self.values],
            after,
        )
