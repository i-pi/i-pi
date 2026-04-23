"""Dependency tracking, lazy caching, and automatic update of variables.

This module provides a small data and storage class used throughout i-PI to
represent physical quantities. Two concrete classes wrap values:

- `depend_value`: wraps an arbitrary Python value.
- `depend_array`: wraps an array/tensor.

Both share the machinery defined in `depend_base`: a `_status` bundle
carrying the tainted flag, the list of dependants, an optional
`synchronizer` that links equivalent quantities (e.g. cartesian vs
normal-mode bead coordinates), and the write-side threadlock. The
`_status` instance is shared by reference between a `depend_array` and
its slice-views so tainting and dependant notification propagate
transparently across them.

Values are recomputed lazily: when an access finds the tainted flag set,
`_func` is called and the cached value is refreshed. Writes push: they
taint all downstream dependants eagerly.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2017 i-PI developers
# See the "licenses" directory for full license information.


import weakref
import threading

from copy import deepcopy

import operator as _op

import numpy as np

from ipi.utils.array_backend import xp, xp_size

# NOTE: the info()/warning() calls below are disabled because they
# add measurable cost in the hot path. Re-enable
# manually when debugging tainting/synchronizer issues.
# from ipi.utils.messages import verbosity, warning, info


__all__ = [
    "depend_value",
    "depend_array",
    "depend_status",
    "synchronizer",
    "dpipe",
    "dcopy",
    "dstrip",
    "depraise",
    "dproperties",
    "ddot",
    "noddot",
]


class synchronizer(object):
    """Registry of synchronized peer depend objects.

    Peers share state: setting one makes the others stale (tainted) until
    they are read, at which point they rebuild themselves from the manually
    set peer via a per-peer entry in `_func`. `manual` tracks which peer was
    last written.
    """

    __slots__ = ("synced", "manual")

    def __init__(self, deps=None):
        self.synced = dict() if deps is None else deps
        self.manual = None


class depend_status(object):
    """Shared per-dependency-group state.

    Slice-views of a `depend_array` and their parent point at the same
    `depend_status` instance, so tainting the parent taints the view (and vice
    versa), and both share the same dependant list and threadlock.
    """

    __slots__ = ("tainted", "dependants", "synchro", "threadlock")

    def __init__(self, tainted=True, dependants=None, synchro=None):
        self.tainted = bool(tainted)
        self.dependants = [] if dependants is None else dependants
        self.synchro = synchro
        self.threadlock = threading.RLock()

    def __getstate__(self):
        # threadlock is unpicklable; recreate it on load.
        return (self.tainted, self.dependants, self.synchro)

    def __setstate__(self, state):
        self.tainted, self.dependants, self.synchro = state
        self.threadlock = threading.RLock()


class depend_base(object):
    """Shared tainting, dependency, and synchronization machinery.

    Attributes:
        _name: Human-readable name (used for diagnostics and synchronizer
            lookup).
        _func: Callable (or dict of callables for synchronizer peers) used
            to recompute the value when tainted.
        _status: Shared `depend_status` bundle holding the tainted flag, the
            dependants list, the synchronizer reference, and the write-side
            threadlock. Slice-views share this instance with their parent.
    """

    # `__weakref__` is needed so `weakref.ref(self)` works (used when we
    # register ourselves as a dependant of another node). Without it,
    # `__slots__` would suppress the weakref support slot that Python
    # normally adds to every class.
    __slots__ = ("_name", "_func", "_status", "__weakref__")

    def __init__(
        self,
        name,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        status=None,
    ):
        self._name = name
        self._func = func
        # A fresh Status unless a caller passes one in (e.g. slice-view
        # construction, or `dcopy` sharing state between peers).
        if status is None:
            status = depend_status(
                tainted=True,
                dependants=[] if dependants is None else dependants,
                synchro=synchro,
            )
        else:
            # Caller-provided status takes precedence; seed dependants/synchro
            # only if we're asked to *add* them. Currently nothing passes both.
            if dependants is not None:
                status.dependants = dependants
            if synchro is not None and status.synchro is None:
                status.synchro = synchro
        self._status = status

        self.add_synchro(status.synchro)

        # set up dependencies (reverse edges from upstream sources).
        if dependencies is not None:
            for item in dependencies:
                item.add_dependant(self, tainted=True)

        # Coerce any raw dependants into weakrefs (slice-view dependants
        # come through already-wrapped; user-constructed ones may not).
        deps = self._status.dependants
        for i, item in enumerate(deps):
            if not isinstance(item, weakref.ref):
                deps[i] = weakref.ref(item)

        # Primitive objects start untainted; computed objects stay tainted.
        if self._func is None:
            self._status.tainted = False

    def __deepcopy__(self, memo):
        """Slot-aware deepcopy.

        Mirrors the historical behaviour: rebuild a fresh instance, copy
        every slot, but (a) share `_func` by reference (bound methods
        often carry unpicklable state like threading.Lock on the
        ``__self__``), and (b) give the copy a fresh threadlock inside
        the copied `depend_status`.
        """
        newone = type(self).__new__(type(self))
        for cls in type(self).__mro__:
            for name in getattr(cls, "__slots__", ()):
                # `__weakref__` is a read-only descriptor slot that Python
                # manages itself; setattr on the fresh instance raises.
                if name == "__weakref__":
                    continue
                try:
                    val = getattr(self, name)
                except AttributeError:
                    continue
                if name == "_func":
                    setattr(newone, name, val)
                elif name == "_status":
                    copied = depend_status.__new__(depend_status)
                    copied.tainted = deepcopy(val.tainted, memo)
                    copied.dependants = deepcopy(val.dependants, memo)
                    copied.synchro = deepcopy(val.synchro, memo)
                    copied.threadlock = threading.RLock()
                    setattr(newone, name, copied)
                else:
                    setattr(newone, name, deepcopy(val, memo))
        return newone

    def add_synchro(self, synchro=None):
        if synchro is None:
            return
        status = self._status
        assert (
            status.synchro is None or status.synchro is synchro
        ), "This object must not have a previous synchronizer!"
        status.synchro = synchro
        if self._name not in synchro.synced:
            synchro.synced[self._name] = self
            synchro.manual = self._name

    def add_dependant(self, newdep, tainted=True):
        newdep.add_dependency(self, tainted=tainted)

    def add_dependency(self, newdep, tainted=True):
        newdep._status.dependants.append(weakref.ref(self))
        if tainted:
            self.taint(taintme=True)

    def taint(self, taintme=True):
        """Plain directed acyclic graph walk: taint dependants, optionally self."""
        for item in self._status.dependants:
            item = item()
            if item is None:
                continue
            if not item._status.tainted:
                item.taint()
        self._status.tainted = taintme

    def _taint_synchro(self):
        """Mark peer synchronized quantities as stale.

        Called from setter paths when the synchronizer is present. Sets
        `synchro.manual` to self and taints all other peers (and their
        dependants) while leaving self clean.
        """
        synchro = self._status.synchro
        if synchro is None:
            return
        synchro.manual = self._name
        for v in synchro.synced.values():
            if v is self or v._status.tainted:
                continue
            v._status.tainted = True
            for item in v._status.dependants:
                item = item()
                if item is None:
                    continue
                if not item._status.tainted:
                    item.taint()

    def tainted(self):
        return self._status.tainted

    def update_auto(self):
        """Recompute from `_func` when tainted.

        For synchronized peers, `_func` is a dict keyed by peer name; the
        value for the manually-set peer is invoked.
        """
        synchro = self._status.synchro
        if synchro is not None:
            if self._name != synchro.manual:
                self.set(self._func[synchro.manual](), manual=False)
                # info(f" @depend: Set value for {self._name} (syncro)", verbosity.debug)
            else:
                # warning(self._name + " probably shouldn't be tainted (synchro)", verbosity.low)
                pass
        elif self._func is not None:
            self.set(self._func(), manual=False)
            # info(f" @depend: Set value for {self._name} (auto)", verbosity.debug)
        else:
            # warning(self._name + " probably shouldn't be tainted (value)", verbosity.low)
            pass

    def update_man(self):
        """Post-manual-write hook: propagate synchro state or reject writes
        to auto-computed quantities."""
        if self._status.synchro is not None:
            self._taint_synchro()
            # info(f" @depend: Set value for {self._name} (manual)", verbosity.debug)
        elif self._func is not None:
            raise NameError(
                "Cannot set manually the value of the automatically-computed property <"
                + self._name
                + ">"
            )

    def set(self, value, manual=False):
        raise ValueError("Undefined set function for base depend class")

    def get(self):
        raise ValueError("Undefined get function for base depend class")

    def _refresh(self):
        """Cold-path recompute. Callers are expected to gate with
        `if self._status.tainted: self._refresh()` so the fast path stays
        inline. Double-checks under the lock to make the update thread-safe."""
        status = self._status
        with status.threadlock:
            if status.tainted:
                self.update_auto()
                self.taint(taintme=False)


class depend_value(depend_base):
    """Depend wrapper around a scalar/Python value."""

    __slots__ = ("_value",)

    def __init__(
        self,
        name,
        value=None,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        status=None,
    ):
        self._value = value
        super().__init__(name, synchro, func, dependants, dependencies, status)

    def get(self):
        return self.__get__(self, self.__class__)

    def __get__(self, instance=None, owner=None):
        if self._status.tainted:
            self._refresh()
        return self._value

    def set(self, value, manual=True):
        status = self._status
        with status.threadlock:
            self._value = value
            if manual:
                self.update_man()
            self.taint(taintme=False)

    def __set__(self, instance, value):
        self.set(value)


def _is_scalar_index(index, ndim):
    """True if `array[index]` returns a scalar for an array of rank `ndim`."""
    # 0-d arrays/tensors index like plain integers (len() would raise on torch).
    if getattr(index, "ndim", None) == 0:
        return ndim <= 1
    if hasattr(index, "__len__"):
        if len(index) == ndim and isinstance(index, tuple):
            for i in index:
                if isinstance(i, slice):
                    return False
            return True
        return False
    return ndim <= 1


class depend_array(depend_base):
    """Depend wrapper around an ndarray (composition, not inheritance).

    The underlying array lives in `self._value`. Slice-views share memory
    with the parent's `_value` and share the parent's `_status` (and thus
    the tainted flag, dependants list, synchronizer, and threadlock), so
    writing through a slice taints the parent's downstream consumers.

    Passing ``on_host=True`` pins the storage to the host (`device="cpu"`
    in array-API terms). Under torch the value is still a `torch.Tensor`
    so the surrounding code stays backend-agnostic, but it lives on CPU
    and reading scalars out of it (`.item()`, `int()`, `bool()`, `%`,
    `==`) does not trigger `cudaStreamSynchronize`. Intended for small
    control-flow arrays (MTS counts, index masks) that are only touched
    for Python branching. Do **not** set it for arrays that get combined
    with backend-native tensors on the hot path (e.g. `pdt`, `qdt_on_m`).
    """

    __slots__ = ("_value", "_parent", "_on_host")

    def __init__(
        self,
        value,
        name,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        status=None,
        base=None,
        parent=None,
        on_host=False,
    ):
        self._on_host = bool(on_host)
        asarray_kwargs = {"device": "cpu"} if self._on_host else {}
        if base is not None:
            # Slice-view: wrap an already-allocated view of the parent;
            # storage location is inherited automatically.
            self._value = base
        elif value is None:
            # Sentinel used by __deepcopy__; members will be overwritten.
            self._value = xp.empty(0, **asarray_kwargs)
        else:
            # xp only handles numeric dtypes; strings/objects stay as numpy.
            # Try xp first so torch tensors (possibly on non-CPU devices)
            # pass through without a numpy round-trip that would force a
            # device->host copy and fail on CUDA.
            try:
                self._value = xp.asarray(value, **asarray_kwargs)
            except (TypeError, ValueError, RuntimeError):
                self._value = np.asarray(value)

        # Slice-views keep a reference to the parent so that refreshing a
        # stale slice delegates the recompute to the parent (which owns
        # the full-shape `_value`). Without this, `update_auto` on a
        # slice would try to assign a full-shape result into the slice's
        # narrower view.
        self._parent = parent

        super().__init__(name, synchro, func, dependants, dependencies, status)

    # ------------------------------------------------------------------
    # Propagate common array attributes to allow direct access without unwrapping.
    # ------------------------------------------------------------------

    @property
    def shape(self):
        return self._value.shape

    @property
    def dtype(self):
        return self._value.dtype

    @property
    def ndim(self):
        return self._value.ndim

    @property
    def size(self):
        return xp_size(self._value)

    @property
    def T(self):
        return self._value.T

    @property
    def flat(self):
        return self._value.flat

    @property
    def nbytes(self):
        return self._value.nbytes

    @property
    def itemsize(self):
        return self._value.itemsize

    @property
    def real(self):
        return self._value.real

    @property
    def imag(self):
        return self._value.imag

    # ------------------------------------------------------------------
    # Cover any other attribute (.astype, .copy, .tobytes, .ctypes, .view, ...)
    # ------------------------------------------------------------------

    def __getattr__(self, name):
        # Called only when normal lookup fails. Guard against recursion during
        # partial initialization: _value may not exist yet.
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            value = object.__getattribute__(self, "_value")
        except AttributeError:
            raise AttributeError(name)
        return getattr(value, name)

    # ------------------------------------------------------------------
    # Arithmetic / comparison dunders
    #
    # Delegate to the underlying `_value`'s native operators so the result
    # type matches the backend.
    # ------------------------------------------------------------------

    # Arithmetic / comparison dunders are attached after the class body
    # (see `_install_depend_array_operators`). Each dunder inlines the
    # tainted-check + other-unwrap to avoid an extra Python call per
    # operation on the hot path.

    # Required because we override __eq__.
    __hash__ = None

    # ------------------------------------------------------------------
    # Shape helpers (explicit so they return the intended type)
    # ------------------------------------------------------------------

    def reshape(self, newshape):
        """Return a reshaped depend_array sharing memory with self."""
        return depend_array(
            value=None,
            base=self._value.reshape(newshape),
            name=self._name,
            status=self._status,
            parent=self,
        )

    def flatten(self):
        """Return a 1D depend_array sharing memory with self."""
        return self.reshape(-1)

    # ------------------------------------------------------------------
    # Descriptor + indexing
    # ------------------------------------------------------------------

    def __get__(self, instance=None, owner=None):
        if self._status.tainted:
            self._refresh()
        return self

    def __set__(self, instance, value):
        self.set(value, manual=True)

    def __getitem__(self, index):
        if self._status.tainted:
            self._refresh()
        if _is_scalar_index(index, self._value.ndim):
            return self._value[index]
        return depend_array(
            value=None,
            base=self._value[index],
            name=self._name,
            status=self._status,
            parent=self,
        )

    def __setitem__(self, index, value):
        # Safety net: some backends (e.g. torch) reject a depend_array RHS
        # outright. Callers are expected to dstrip() explicitly; this peel
        # costs ~12 ns and only fires when someone forgot.
        if isinstance(value, depend_base):
            value = value._value
        status = self._status
        with status.threadlock:
            self._value[index] = value
            # Skip update_man() on plain leaf-storage arrays (no synchro,
            # no computed func) — there is nothing for it to do.
            if status.synchro is not None or self._func is not None:
                self.update_man()
            self.taint(taintme=False)

    def set(self, value, manual=True):
        """Full-array assignment. `manual=False` is used by `update_auto`."""
        if isinstance(value, depend_base):
            value = value._value
        status = self._status
        with status.threadlock:
            self._value[...] = value
            if manual and (status.synchro is not None or self._func is not None):
                self.update_man()
            self.taint(taintme=False)

    def get(self):
        if self._status.tainted:
            self._refresh()
        return self

    def _refresh(self):
        """Slice-view aware refresh.

        A slice-view shares `_status` with its parent but its `_value`
        is only a narrow view. Recomputing through `update_auto` must
        happen on the parent, which owns the full-shape storage; the
        slice's `_value` (a view) then observes the refreshed data for
        free, and the shared `tainted` flag is cleared by the parent.
        """
        if self._parent is not None:
            self._parent._refresh()
        else:
            super()._refresh()

    # ------------------------------------------------------------------
    # len / iter / repr / bool
    # ------------------------------------------------------------------

    def __len__(self):
        return len(self._value)

    def __iter__(self):
        if self._status.tainted:
            self._refresh()
        return iter(self._value)

    def __repr__(self):
        return f"depend_array({self._name!r}, {self._value!r})"

    def __bool__(self):
        if self._status.tainted:
            self._refresh()
        return bool(self._value)

    # ------------------------------------------------------------------
    # Pickling: `depend_status.__getstate__` already omits the threadlock, and
    # slot classes pickle naturally via `__reduce_ex__`. No custom methods
    # needed.
    # ------------------------------------------------------------------


# ----------------------------------------------------------------------
# Install arithmetic / comparison dunders on depend_array.
#
# Factored here (not in the class body) so a single template generates
# all ~25 dunders. Each generated method inlines the tainted-check +
# operand-unwrap, so the hot path costs one Python call (the dunder
# itself) plus one bool check.
# ----------------------------------------------------------------------


def _make_binop(op):
    def method(self, other):
        if self._status.tainted:
            self._refresh()
        if isinstance(other, depend_array):
            if other._status.tainted:
                other._refresh()
            other = other._value
        return op(self._value, other)

    method.__name__ = f"__{op.__name__.strip('_')}__"
    return method


def _make_rbinop(op):
    def method(self, other):
        if self._status.tainted:
            self._refresh()
        if isinstance(other, depend_array):
            if other._status.tainted:
                other._refresh()
            other = other._value
        return op(other, self._value)

    method.__name__ = f"__r{op.__name__.strip('_')}__"
    return method


def _make_unop(op, name):
    def method(self):
        if self._status.tainted:
            self._refresh()
        return op(self._value)

    method.__name__ = name
    return method


_BINOPS = [
    ("__add__", "__radd__", _op.add),
    ("__sub__", "__rsub__", _op.sub),
    ("__mul__", "__rmul__", _op.mul),
    ("__truediv__", "__rtruediv__", _op.truediv),
    ("__floordiv__", "__rfloordiv__", _op.floordiv),
    ("__mod__", "__rmod__", _op.mod),
    ("__pow__", "__rpow__", _op.pow),
    ("__matmul__", "__rmatmul__", _op.matmul),
    ("__and__", "__rand__", _op.and_),
    ("__or__", "__ror__", _op.or_),
    ("__xor__", "__rxor__", _op.xor),
]

for _fwd, _rev, _fn in _BINOPS:
    setattr(depend_array, _fwd, _make_binop(_fn))
    setattr(depend_array, _rev, _make_rbinop(_fn))

for _name, _fn in (
    ("__eq__", _op.eq),
    ("__ne__", _op.ne),
    ("__lt__", _op.lt),
    ("__le__", _op.le),
    ("__gt__", _op.gt),
    ("__ge__", _op.ge),
):
    setattr(depend_array, _name, _make_binop(_fn))

setattr(depend_array, "__neg__", _make_unop(_op.neg, "__neg__"))
setattr(depend_array, "__pos__", _make_unop(_op.pos, "__pos__"))
setattr(depend_array, "__abs__", _make_unop(abs, "__abs__"))


# ----------------------------------------------------------------------
# Module-level helpers
# ----------------------------------------------------------------------


# Kept for backward compatibility with callers that imported `noddot`.
# The historical `np.dot = noddot` global patch is removed: `__array_ufunc__`
# now handles `np.dot` correctly for depend_array inputs.
noddot = np.dot


def ddot(da, db):
    """`np.dot` on the unwrapped operands."""
    return np.dot(dstrip(da), dstrip(db))


def dstrip(da):
    """Return the underlying array/value without depend machinery.

    Passes non-depend inputs through unchanged. This is needed because
    a computed `depend_value`/`depend_array` accessed via the
    `dproperties` descriptor returns the raw ndarray/scalar (not the
    wrapper), and several call sites `dstrip` that result a second
    time. The historical low-verbosity warning for the non-depend case
    has been removed.
    """
    if isinstance(da, (depend_array, depend_value)):
        return da._value
    return da


def dpipe(dfrom, dto, item=-1):
    """Synchonizes two depend objects.

    Takes two depend objects, and makes one of them depend on the other in such
    a way that both keep the same value. Used for attributes such as
    temperature that are used in many different modules, and so need different
    depend objects in each, but which should all have the same value.

    Args:
        dfrom: The object that is depend on.
        dto: The object that depends on the former one.
    """

    if item < 0:
        dto._func = lambda: dstrip(dfrom.__get__(dfrom, dfrom.__class__))
    else:
        dto._func = lambda i=item: dfrom.__getitem__(i)
    dto.add_dependency(dfrom)


def dcopy(dfrom, dto):
    """Share dependency state between two depend objects.

    After this call `dto` observes the same tainted flag, dependants list,
    synchronizer and threadlock as `dfrom` (they point to the same
    `depend_status` instance). `_func` is shared by reference. Used for
    attributes such as `fx`/`fy`/`fz` that are views into `_f` and must
    invalidate alongside it.
    """
    dto._status = dfrom._status
    dto._func = dfrom._func
    # Register dto under its own name in the shared synchronizer, if any.
    synchro = dfrom._status.synchro
    if synchro is not None and dto._name not in synchro.synced:
        synchro.synced[dto._name] = dto


def depraise(exception):
    raise exception


def _inject_depend_property(cls, attr_name):
    """Adds a property `<attr_name>` to `cls` that forwards to the depend member
    `_<attr_name>`. The depend member must already exist on the class; this is
    intended to be used with `depend_value` or `depend_array` members defined in
    the class body, and is used by `dproperties` to expose them as attributes."""
    private_name = f"_{attr_name}"

    def getter(self):
        dep = self.__dict__[private_name]
        if dep._status.tainted:
            dep._refresh()
        return dep._value if type(dep) is depend_value else dep

    # `__set__` on both depend_value and depend_array just forwards to
    # `.set(value)`; call it directly.
    def setter(self, value):
        self.__dict__[private_name].set(value)

    setattr(cls, attr_name, property(getter, setter))


def dproperties(cls, attr_names):
    """Expose depend members `_<name>` as attributes `<name>` on `cls`.

    Example: `dproperties(A, ["val"])` on a class with `self._val =
    depend_value(...)` makes `obj.val` read/write through the descriptor.
    """
    if not isinstance(attr_names, list):
        attr_names = [attr_names]
    for attr_name in attr_names:
        _inject_depend_property(cls, attr_name)
