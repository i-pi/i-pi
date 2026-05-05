"""Dependency tracking, lazy caching, and automatic update of variables.

This module provides a small data and storage class used throughout i-PI to
represent physical quantities. Two concrete classes wrap values:

- `depend_value`: wraps an arbitrary Python value.
- `depend_array`: wraps an array/tensor

Both share the machinery defined in `depend_base`: a tainted flag, an
optional recompute function `_func`, an optional `synchronizer` that links
equivalent quantities (e.g. cartesian vs normal-mode bead coordinates), and
a list of dependants to notify on change.

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

# NOTE: the info()/warning() calls below are disabled because they
# add measurable cost in the hot path. Re-enable
# manually when debugging tainting/synchronizer issues.
# from ipi.utils.messages import verbosity, warning, info


__all__ = [
    "depend_value",
    "depend_array",
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

    def __init__(self, deps=None):
        self.synced = dict() if deps is None else deps
        self.manual = None


class depend_base(object):
    """Shared tainting, dependency, and synchronization machinery.

    Attributes:
        _tainted: 1-element bool ndarray; True if value must be recomputed.
        _func: Callable (or dict of callables for synchronizer peers) used
            to recompute the value when tainted.
        _name: Human-readable name (used for diagnostics and synchronizer
            lookup).
        _synchro: Optional `synchronizer` linking peer quantities.
        _dependants: Weak references to objects that must be tainted when
            this one changes.
    """

    def __init__(
        self,
        name,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        tainted=None,
    ):
        if tainted is None:
            tainted = np.array([True], bool)
        if dependants is None:
            dependants = []
        if dependencies is None:
            dependencies = []

        self._tainted = tainted
        self._func = func
        self._name = name
        self._threadlock = threading.RLock()
        self._dependants = []
        self._synchro = None

        self.add_synchro(synchro)

        # set up dependencies and dependants given on initialization
        for item in dependencies:
            item.add_dependant(self, tainted)

        for item in dependants:
            if not isinstance(item, weakref.ref):
                dependants.remove(item)
                dependants.append(weakref.ref(item))
        self._dependants = dependants

        # Primitive objects start untainted; computed objects inherit the
        # incoming tainted flag.
        if tainted[0]:
            if self._func is None:
                self.taint(taintme=False)
            else:
                self.taint(taintme=tainted)

    def __deepcopy__(self, memo):
        # Deepcopy every attribute except `_threadlock` (unpicklable RLock)
        # and `_func` (a bound method whose `__self__` may hold non-picklable
        # state such as a `threading.Lock`; the copy shares the reference).
        newone = type(self).__new__(type(self))
        for member, value in self.__dict__.items():
            if member == "_threadlock":
                newone._threadlock = threading.RLock()
            elif member == "_func":
                newone._func = value
            else:
                setattr(newone, member, deepcopy(value, memo))
        return newone

    def add_synchro(self, synchro=None):
        assert (
            self._synchro is None
        ), "This object must not have a previous synchronizer!"
        self._synchro = synchro
        if self._synchro is not None and self._name not in self._synchro.synced:
            self._synchro.synced[self._name] = self
            self._synchro.manual = self._name

    def add_dependant(self, newdep, tainted=True):
        newdep.add_dependency(self, tainted=tainted)

    def add_dependency(self, newdep, tainted=True):
        newdep._dependants.append(weakref.ref(self))
        if tainted:
            self.taint(taintme=True)

    def taint(self, taintme=True):
        """Plain directed acyclic graph walk: taint dependants, optionally self."""
        for item in self._dependants:
            item = item()
            if item is None:
                continue
            if not item._tainted[0]:
                item.taint()
        self._tainted[0] = taintme

    def _taint_synchro(self):
        """Mark peer synchronized quantities as stale.

        Called from setter paths when `_synchro` is present. Sets
        `_synchro.manual` to self and taints all other peers (and their
        dependants) while leaving self clean.
        """
        if self._synchro is None:
            return
        self._synchro.manual = self._name
        for v in self._synchro.synced.values():
            if v is self or v._tainted[0]:
                continue
            v._tainted[0] = True
            for item in v._dependants:
                item = item()
                if item is None:
                    continue
                if not item._tainted[0]:
                    item.taint()

    def tainted(self):
        return self._tainted[0]

    def update_auto(self):
        """Recompute from `_func` when tainted.

        For synchronized peers, `_func` is a dict keyed by peer name; the
        value for the manually-set peer is invoked.
        """
        if self._synchro is not None:
            if self._name != self._synchro.manual:
                self.set(self._func[self._synchro.manual](), manual=False)
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
        if self._synchro is not None:
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
        `if self._tainted[0]: self._refresh()` so the fast path stays inline.
        Double-checks under the lock to make the update thread-safe."""
        with self._threadlock:
            if self._tainted[0]:
                self.update_auto()
                self.taint(taintme=False)


class depend_value(depend_base):
    """Depend wrapper around a scalar/Python value."""

    def __init__(
        self,
        name,
        value=None,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        tainted=None,
    ):
        self._value = value
        super().__init__(name, synchro, func, dependants, dependencies, tainted)

    def get(self):
        return self.__get__(self, self.__class__)

    def __get__(self, instance=None, owner=None):
        if self._tainted[0]:
            self._refresh()
        return self._value

    def set(self, value, manual=True):
        with self._threadlock:
            self._value = value
            if manual:
                self.update_man()
            self.taint(taintme=False)

    def __set__(self, instance, value):
        self.set(value)


def _is_scalar_index(index, ndim):
    """True if `array[index]` returns a scalar for an array of rank `ndim`."""
    if isinstance(index, (slice, type(Ellipsis))):
        return False
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
    with the parent's `_value` and share the parent's `_tainted` flag and
    `_dependants` list, so writing through a slice taints the parent's
    downstream consumers.
    """

    def __init__(
        self,
        value,
        name,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        tainted=None,
        base=None,
        parent=None,
    ):
        if base is not None:
            # Slice-view: wrap an already-allocated view of the parent.
            self._value = base
        elif value is None:
            # Sentinel used by __deepcopy__; members will be overwritten.
            self._value = np.empty(0)
        elif isinstance(value, np.ndarray):
            self._value = value
        else:
            self._value = np.asarray(value)

        # Slice-views keep a reference to the parent so that refreshing a
        # stale slice delegates the recompute to the parent (which owns
        # the full-shape `_value`). Without this, `update_auto` on a
        # slice would try to assign a full-shape result into the slice's
        # narrower view.
        self._parent = parent

        super().__init__(name, synchro, func, dependants, dependencies, tainted)

    # ------------------------------------------------------------------
    # numpy interop
    # ------------------------------------------------------------------

    def __array__(self, dtype=None, copy=None):
        """Lets np.asarray / np.* accept us transparently.

        Zero-copy fast path when the dtype matches.
        """
        if self._tainted[0]:
            self._refresh()
        v = self._value
        if dtype is None or dtype == v.dtype:
            return v
        return v.astype(dtype)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Unwrap depend_array inputs, run the ufunc on raw arrays."""
        unwrapped = []
        for x in inputs:
            if isinstance(x, depend_array):
                if x._tainted[0]:
                    x._refresh()
                unwrapped.append(x._value)
            else:
                unwrapped.append(x)
        if "out" in kwargs:
            kwargs["out"] = tuple(
                o._value if isinstance(o, depend_array) else o for o in kwargs["out"]
            )
        return getattr(ufunc, method)(*unwrapped, **kwargs)

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
        return self._value.size

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

    def reshape(self, *shape):
        """Return a reshaped depend_array sharing memory with self.

        Accepts both `da.reshape((2, 3))` and `da.reshape(2, 3)`, matching
        numpy's signature.
        """
        if len(shape) == 1:
            shape = shape[0]
        return depend_array(
            value=None,
            base=self._value.reshape(shape),
            name=self._name,
            synchro=self._synchro,
            dependants=self._dependants,
            tainted=self._tainted,
            parent=self,
        )

    def flatten(self):
        """Return a 1D depend_array sharing memory with self."""
        return self.reshape(self._value.size)

    # ------------------------------------------------------------------
    # Descriptor + indexing
    # ------------------------------------------------------------------

    def __get__(self, instance=None, owner=None):
        if self._tainted[0]:
            self._refresh()
        return self

    def __set__(self, instance, value):
        self.set(value, manual=True)

    def __getitem__(self, index):
        if self._tainted[0]:
            self._refresh()
        if _is_scalar_index(index, self._value.ndim):
            return self._value[index]
        return depend_array(
            value=None,
            base=self._value[index],
            name=self._name,
            synchro=self._synchro,
            dependants=self._dependants,
            tainted=self._tainted,
            parent=self,
        )

    def __setitem__(self, index, value):
        with self._threadlock:
            self._value[index] = value
            # Skip update_man() on plain leaf-storage arrays (no synchro,
            # no computed func) — there is nothing for it to do.
            if self._synchro is not None or self._func is not None:
                self.update_man()
            self.taint(taintme=False)

    def set(self, value, manual=True):
        """Full-array assignment. `manual=False` is used by `update_auto`.

        Note: this writes through `_value[...] = value` so child slice-views
        (which hold references to the same underlying buffer) stay valid.
        Do NOT replace `self._value` with a new ndarray — that would
        silently desync any existing slice's `_value` from the parent's.
        """
        with self._threadlock:
            self._value[...] = value
            if manual and (self._synchro is not None or self._func is not None):
                self.update_man()
            self.taint(taintme=False)

    def get(self):
        if self._tainted[0]:
            self._refresh()
        return self

    def _refresh(self):
        """Slice-view aware refresh.

        A slice-view shares `_tainted` with its parent but its `_value`
        is only a narrow view. Recomputing through `update_auto` must
        happen on the parent, which owns the full-shape storage; the
        slice's `_value` (a view) then observes the refreshed data for
        free, and the shared `_tainted` flag is cleared by the parent.
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
        if self._tainted[0]:
            self._refresh()
        return iter(self._value)

    def __repr__(self):
        return f"depend_array({self._name!r}, {self._value!r})"

    def __bool__(self):
        if self._tainted[0]:
            self._refresh()
        return bool(self._value)

    # ------------------------------------------------------------------
    # Support pickling by ignoring the unpicklable threadlock and rebuilding it on unpickle.
    # ------------------------------------------------------------------

    def __getstate__(self):
        state = {k: v for k, v in self.__dict__.items() if k != "_threadlock"}
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._threadlock = threading.RLock()


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
        if self._tainted[0]:
            self._refresh()
        if isinstance(other, depend_array):
            if other._tainted[0]:
                other._refresh()
            other = other._value
        return op(self._value, other)

    method.__name__ = f"__{op.__name__.strip('_')}__"
    return method


def _make_rbinop(op):
    def method(self, other):
        if self._tainted[0]:
            self._refresh()
        if isinstance(other, depend_array):
            if other._tainted[0]:
                other._refresh()
            other = other._value
        return op(other, self._value)

    method.__name__ = f"__r{op.__name__.strip('_')}__"
    return method


def _make_unop(op, name):
    def method(self):
        if self._tainted[0]:
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
        dto._func = lambda: dfrom.__get__(dfrom, dfrom.__class__)
    else:
        dto._func = lambda i=item: dfrom.__getitem__(i)
    dto.add_dependency(dfrom)


def dcopy(dfrom, dto):
    """Copies the dependencies of one depend object to another.

    Args:
        see dpipe.
    """
    dto._dependants = dfrom._dependants
    dto._synchro = dfrom._synchro
    dto.add_synchro(dfrom._synchro)
    dto._tainted = dfrom._tainted
    dto._func = dfrom._func


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
        if dep._tainted[0]:
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
