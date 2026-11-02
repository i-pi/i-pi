"""Compilation of a function that operates on depend objects with jax.jit, and
replay of the compiled function instead of the function.

This is the counterpart of the CUDA graphs of `ipi.utils.cudagraphs` for the
JAX backend, and it uses the same recording of the operations on the depend
objects (see `ipi.utils.jit`). The function runs eagerly for the first calls,
which are recorded to find the depend objects that it reads before writing
them, i.e. its inputs, until two calls read the same ones (e.g. positions
written before a call make it recompute the forces, which are then not
inputs). The next call traces the function with jax.jit, with the inputs
brought up to date as when the compiled function is called: the inputs that
hold JAX arrays are the arguments of the traced function, which returns the
values of the depend objects that the function writes, i.e. its outputs. The
compiled function then replaces the function: at each call, the inputs are
brought up to date and passed to it, its results are stored in the outputs,
and the taint status of the depend objects is restored as it was at the end
of the traced call.

The other inputs (e.g. the time step, held by a Python float) are fixed in
the compiled function, and a `JitError` is raised if they change. Unlike CUDA
graphs, the compiled function cannot be interrupted to run parts of the
function eagerly: the calls to `eager` are traced, and the function is not
compiled if it must be interrupted at a `split`.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

from ipi.utils import depend
from ipi.utils.jit import WARMUP_CALLS, JitError, Replay, _input, _Recorder
from ipi.utils.messages import info, verbosity, warning

__all__ = ["JaxReplay"]

# eager calls after which the function is not compiled, if no two consecutive
# ones have read the same depend objects
MAX_EAGER_CALLS = 8


class JaxReplay(Replay):
    """Runs a function eagerly for the first calls, until two calls read the
    same depend objects, compiles it with jax.jit at the next call, and calls
    the compiled function at the following calls. The arguments of the calls
    are only passed to the eager calls and to the traced one.

    Attributes:
        func: The function.
        prng: The random number generator whose JAX key is used by the
            function, if any.
        compiled: The compiled function, or None before the compilation.
    """

    method = "compiled with jax.jit"

    def __init__(self, func, prng=None):
        super().__init__(func)
        self.prng = prng
        self.ncalls = 0
        # the depend objects read before being written by the last observed
        # call, and by the function once two observed calls agree
        self.observed = None
        self.inputs = None
        # why the function cannot be compiled, if it cannot
        self.reason = None

    def ready(self):
        return self.reason is None and self.inputs is not None

    def warm(self, *args):
        if self.reason is not None or self.ncalls < WARMUP_CALLS - 1:
            self.func(*args)
        else:
            inputs = self.observe(*args)
            if inputs == self.observed:
                self.inputs = inputs
            elif self.ncalls >= MAX_EAGER_CALLS:
                self.reason = (
                    "the depend objects that the steps read change at every step, "
                    "e.g. as the state is written between the steps"
                )
            self.observed = inputs
            if self.reason is not None:
                warning(f"The steps are not compiled: {self.reason}", verbosity.low)
        self.ncalls += 1

    def observe(self, *args):
        """Calls the function eagerly, and returns the depend objects that it
        reads before writing them."""

        recorder = _Recorder(self)
        depend.set_recorder(recorder)
        try:
            self.func(*args)
        finally:
            depend.set_recorder(None)
        return [dep for dep, value, copy in recorder.inputs.values()]

    def record(self, *args):
        """Traces and compiles the function, and calls the compiled function."""

        import jax

        # the inputs are brought up to date, as when the compiled function is called
        for dep in self.inputs:
            dep.get()
        values = [dep._load() for dep in self.inputs]
        arrays = [v for v in values if _is_jax(v)]
        # the inputs passed to the compiled function, and the fixed ones
        self.arrays = [dep for dep, v in zip(self.inputs, values) if _is_jax(v)]
        self.fixed = [
            _input(dep, v) for dep, v in zip(self.inputs, values) if not _is_jax(v)
        ]
        key = None if self.prng is None else self.prng.jax_key

        def traced(key, arrays):
            for dep, value in zip(self.arrays, arrays):
                dep._store(value)
            if key is not None:
                self.prng.jax_key = key
            self.recorder = _TraceRecorder(self, {dep._status for dep in self.inputs})
            depend.set_recorder(self.recorder)
            try:
                self.func(*args)
            finally:
                depend.set_recorder(None)
            key = None if key is None else self.prng.jax_key
            return key, self.outputs_of(self.recorder)

        try:
            self.compiled = jax.jit(traced).lower(key, arrays).compile()
        except Exception as err:
            raise self.failed(err)
        finally:
            for dep, value in zip(self.arrays, arrays):
                dep._store(value)
            if key is not None:
                self.prng.jax_key = key
        self.statuses, self.synchros = self.recorder.final_state()
        self.run(key, arrays)
        info(
            f" @jaxjit: step compiled, with {len(arrays)} input and "
            f"{len(self.outputs)} output arrays",
            verbosity.low,
        )

    def outputs_of(self, recorder):
        """Returns the arrays held by the depend objects written during the
        traced call, and stores how to put them back in the objects."""

        import jax

        if recorder.host_writes:
            raise JitError(
                f"The step modifies host arrays: {sorted(recorder.host_writes)}"
            )
        self.outputs = list(recorder.arrays) + list(recorder.values)
        leaves, self.treedef = jax.tree_util.tree_flatten(
            [dep._load() for dep in self.outputs]
        )
        # the leaves that are not arrays are fixed in the compiled function
        self.leaves = [None if _is_jax(leaf) else leaf for leaf in leaves]
        for dep in recorder.values:
            if not any(_is_jax(leaf) for leaf in jax.tree_util.tree_leaves(dep._value)):
                raise JitError(f"The step modifies the host value <{dep._name}>")
        return [leaf for leaf in leaves if _is_jax(leaf)]

    def run(self, key, arrays):
        """Calls the compiled function, and stores its results."""

        import jax

        try:
            key, results = self.compiled(key, arrays)
        except TypeError as err:
            # the compiled function checks the shapes and types of its inputs
            raise JitError(
                "The shape or type of an input of the compiled step has changed: "
                "the parameters of the simulation that are fixed in the compiled "
                "step cannot change when running with --jit"
            ) from err
        if key is not None:
            self.prng.jax_key = key
        results = iter(results)
        leaves = [next(results) if leaf is None else leaf for leaf in self.leaves]
        for dep, value in zip(
            self.outputs, jax.tree_util.tree_unflatten(self.treedef, leaves)
        ):
            dep._store(value)

    def replay(self):
        """Calls the compiled function on the current values of the inputs, and
        restores the state of the depend objects at the end of the traced call."""

        key = None if self.prng is None else self.prng.jax_key
        self.run(key, [dep.get() for dep in self.arrays])
        for status, tainted in self.statuses:
            status.tainted = tainted
        for synchro, manual in self.synchros:
            synchro.manual = manual

    def cannot(self, reason):
        """Called where the function must run eagerly: the function is then
        not compiled, or cannot be traced."""

        if self.inputs is None:
            self.reason = reason
        else:
            raise JitError(f"The step cannot be compiled: {reason}")

    def call(self, func, *args):
        """Called by `eager(func, *args)`: the call is traced as part of the
        function."""

        return func(*args)


class _TraceRecorder(_Recorder):
    """Recorder of the traced call, which checks that the function reads the
    inputs found when it was observed, and collects the arrays it writes."""

    def __init__(self, replay, inputs):
        super().__init__(replay)
        self.known = inputs
        self.arrays = {}

    def read(self, dep):
        status = dep._status
        if status not in self.written and status not in self.known:
            # its value would be fixed in the compiled function
            raise JitError(
                f"<{dep._name}> is read by the traced step, but not by the eager one"
            )

    def tainted(self, dep, taintme):
        super().tainted(dep, taintme)
        if not taintme and isinstance(dep, depend.depend_array):
            self.arrays[dep._root()] = None


def _is_jax(value):
    import jax

    return isinstance(value, jax.Array)
