"""Dependency tracking, lazy caching, and automatic update of variables.

This module provides the storage classes used throughout i-PI to represent
physical quantities. Two concrete classes wrap values:

- `depend_value`: wraps an arbitrary Python value.
- `depend_array`: wraps an array.

Both share the machinery defined in `depend_base`: an optional recompute
function `_func`, and a `depend_status` object holding a tainted flag, the
list of dependants to notify on change, and an optional `depend_synchro`
that links two representations of the same data (e.g. cartesian and
normal-mode bead coordinates), created with `dsynchro`. An object is
computed if it is created with `func`, or later made computed with `dcompute`
(or `dpipe`, to follow the value of another object). A computed object can be
given more dependencies with `ddeps`.

Values are recomputed lazily: when an access finds the tainted flag set,
`_func` is called and the cached value is refreshed. Writes push: they
taint all downstream dependants eagerly.

Classes expose their depend members through `dproperties`. For a member
stored as `obj.x_`, reading `obj.x` returns the current value (recomputed if
needed) and `obj.x = value` overwrites it. The depend object `obj.x_` is used
to declare dependencies, to modify an array in place (`obj.x_[i] = value`,
`obj.x_ += value`), and to obtain views of parts of an array that stay
connected to the dependency graph (`obj.x_[i]`, `obj.x_.flatten()`). As for
numpy arrays, indexing with masks or arrays of indices gives a copy, which is
a plain array: `obj.x_[mask] op= value` writes it back.

The value read from `obj.x` (or with `get()`) must not be modified, and is
valid until the object is written or recomputed: depending on the backend,
later writes may or may not change it. Code that needs the stored value
without refreshing it (e.g. a stale value, to seed a new calculation) reads it
with `_load`.

Arrays are stored using the array backend selected in
`ipi.utils.array_backend`. Arrays stored as numpy arrays are read as
non-writeable views, so modifying them without going through the depend
object raises an error. JAX arrays, which cannot be modified in place, are
stored by `depend_immutable` objects, which replace them when they are written.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2017 i-PI developers
# See the "licenses" directory for full license information.


import functools
import types
import weakref
import threading

from copy import deepcopy

import operator as _op

import array_api_compat
import numpy as np

from ipi.utils.array_backend import xp

__all__ = [
    "depend_value",
    "depend_array",
    "dcompute",
    "ddeps",
    "dpipe",
    "dsynchro",
    "dproperties",
]

# object notified of the reads and writes of depend objects while a function
# is recorded as CUDA graphs (see ipi.utils.cudagraphs), None otherwise
_recorder = None


def set_recorder(recorder):
    """Sets the object notified of the operations on depend objects."""

    global _recorder
    _recorder = recorder


class depend_synchro(object):
    """Link between two depend objects holding two representations of the
    same data, created by `dsynchro`.

    Attributes:
        peers: The two synchronized depend objects.
        manual: Status of the peer that holds the current data, i.e. the one
            that was written last.
    """

    __slots__ = ("peers", "manual")

    def __init__(self, first, second):
        self.peers = (first, second)
        self.manual = first._status


class depend_status(object):
    """State shared by a depend object and the views created from it.

    Attributes:
        tainted: True if the value must be recomputed before being used.
        dependants: Weak references to the objects that must be tainted when
            the value changes.
        synchro: Optional `depend_synchro` linking the object to a peer.
        lock: Lock making updates of the value thread-safe.
    """

    __slots__ = ("tainted", "dependants", "synchro", "lock")

    def __init__(self, tainted=False):
        self.tainted = tainted
        self.dependants = []
        self.synchro = None
        self.lock = threading.RLock()

    def __getstate__(self):
        return (self.tainted, self.dependants, self.synchro)

    def __setstate__(self, state):
        self.tainted, self.dependants, self.synchro = state
        self.lock = threading.RLock()


class depend_base(object):
    """Shared tainting, dependency, and synchronization machinery.

    Attributes:
        _status: `depend_status` object, shared with views of the same data.
        _func: Callable used to recompute the value when tainted. For
            synchronized objects, it computes the value from the peer.
        _name: Human-readable name, used in error messages.
        _value: The stored value.
    """

    __slots__ = ("_status", "_func", "_name", "_value", "__weakref__")

    def __init__(self, name, func=None, dependencies=None):
        self._func = func
        self._name = name
        self._status = depend_status()
        if dependencies is not None:
            _add_dependencies(self, dependencies)
        # computed objects start tainted, primitive objects start clean
        self._status.tainted = func is not None

    def __deepcopy__(self, memo):
        # `_func` is shared rather than copied: it is often a bound method
        # whose `__self__` may hold non-copyable state such as a lock.
        cls = type(self)
        newone = cls.__new__(cls)
        memo[id(self)] = newone
        for member in _slot_names(cls):
            if hasattr(self, member):
                value = getattr(self, member)
                if member != "_func":
                    value = deepcopy(value, memo)
                setattr(newone, member, value)
        return newone

    def __getstate__(self):
        return {
            m: getattr(self, m) for m in _slot_names(type(self)) if hasattr(self, m)
        }

    def __setstate__(self, state):
        for member, value in state.items():
            setattr(self, member, value)

    def taint(self, taintme=True):
        """Plain directed acyclic graph walk: taint dependants, optionally self."""

        for item in self._status.dependants:
            item = item()
            if item is not None and not item._status.tainted:
                item.taint()
        self._status.tainted = taintme
        if _recorder is not None:
            _recorder.tainted(self, taintme)

    def _taint_synchro(self):
        """Marks self as holding the current data of a synchronized pair, and
        taints the peer (and its dependants)."""

        synchro = self._status.synchro
        synchro.manual = self._status
        for peer in synchro.peers:
            if peer._status is not self._status and not peer._status.tainted:
                peer.taint()
        if _recorder is not None:
            _recorder.synchronized(synchro)

    def tainted(self):
        return self._status.tainted

    def update_auto(self):
        """Recomputes the value from `_func`. A synchronized object is
        recomputed only if its peer holds the current data."""

        synchro = self._status.synchro
        if self._func is not None and (
            synchro is None or synchro.manual is not self._status
        ):
            self.set(self._func(), manual=False)

    def update_man(self):
        """Post-manual-write hook: propagates synchro state or rejects writes
        to auto-computed quantities."""

        if self._status.synchro is not None:
            self._taint_synchro()
        elif self._func is not None:
            raise NameError(
                "Cannot set manually the value of the automatically-computed property <"
                + self._name
                + ">"
            )

    def set(self, value, manual=True):
        raise ValueError("Undefined set function for base depend class")

    def get(self):
        """Returns the current value, recomputing it if needed."""

        if self._status.tainted:
            self._refresh()
        if _recorder is not None:
            _recorder.read(self)
        return self._value

    def _refresh(self):
        """Cold-path recompute. Callers are expected to gate with
        `if self._status.tainted: self._refresh()` so the fast path stays inline.
        Double-checks under the lock to make the update thread-safe."""

        status = self._status
        with status.lock:
            if status.tainted:
                self.update_auto()
                self.taint(taintme=False)

    def _root(self):
        """Returns the object that holds the stored value."""

        return self

    def _load(self):
        """Returns the stored value, without recomputing it."""

        return self._value


def _slot_names(cls):
    """Names of all the slots defined by `cls` and its bases, except `__weakref__`."""

    return [
        s
        for c in cls.__mro__
        for s in getattr(c, "__slots__", ())
        if s != "__weakref__"
    ]


def _add_dependencies(dep, dependencies):
    """Makes `dep` a dependant of each object in `dependencies`."""

    for item in dependencies:
        item._status.dependants.append(weakref.ref(dep))


class depend_value(depend_base):
    """Depend wrapper around a scalar/Python value."""

    __slots__ = ()

    def __init__(
        self,
        name,
        value=None,
        func=None,
        dependencies=None,
    ):
        self._value = value
        super().__init__(name, func, dependencies)

    def set(self, value, manual=True):
        with self._status.lock:
            self._value = value
            if manual:
                self.update_man()
            self.taint(taintme=False)


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


def _is_basic_index(index):
    """True if `array[index]` is a view of the array, i.e. if the index only holds
    integers, slices, Ellipsis and None."""
    for item in index if type(index) is tuple else (index,):
        if not (
            item is None
            or item is Ellipsis
            or type(item) in (int, slice)
            or isinstance(item, np.integer)
        ):
            return False
    return True


def _readonly(array):
    """Returns a non-writeable view of a numpy array, and other arrays as they are."""

    if not isinstance(array, np.ndarray):
        return array
    view = array.view()
    view.flags.writeable = False
    return view


# index of the depend_arrays that are not views obtained by indexing
_NO_INDEX = object()


class depend_array(depend_base):
    """Depend wrapper around an array.

    The array is stored using the selected array backend, or as a numpy array
    if `host` is True (for arrays of strings, or small arrays of metadata that
    are only used on the host). The buffer `_wvalue` is modified by the write
    methods, and the current value is read from `_value`, which for numpy
    arrays is a non-writeable view of `_wvalue`, and is `_wvalue` itself
    otherwise. Views created by indexing or reshaping share memory and status
    with the parent, so writing through a view taints the parent's
    dependants, and a stale view is refreshed by recomputing the parent.
    """

    __slots__ = ("_wvalue", "_parent", "_index")

    def __init__(
        self,
        value,
        name,
        func=None,
        dependencies=None,
        host=False,
    ):
        if isinstance(value, np.ndarray) and not value.flags.writeable:
            value = value.copy()
        wvalue = np.asarray(value) if host else xp.asarray(value)
        if array_api_compat.is_jax_array(wvalue):
            self.__class__ = depend_immutable
        self._wvalue = wvalue
        self._value = _readonly(wvalue)
        self._parent = None
        self._index = _NO_INDEX
        super().__init__(name, func, dependencies)

    def _view(self, value, wvalue, index=_NO_INDEX):
        """Creates a depend_array sharing status and memory with self."""

        view = depend_array.__new__(depend_array)
        view._status = self._status
        view._func = None
        view._name = self._name
        view._value = value
        view._wvalue = wvalue
        view._parent = self
        view._index = index
        return view

    def __deepcopy__(self, memo):
        newone = super().__deepcopy__(memo)
        newone._value = _readonly(newone._wvalue)
        return newone

    def __setstate__(self, state):
        super().__setstate__(state)
        self._value = _readonly(self._wvalue)

    def __array__(self, dtype=None, copy=None):
        raise TypeError(
            f"depend_array <{self._name}> is not an array: read its value "
            "through the corresponding property, or with get()"
        )

    def __repr__(self):
        return f"depend_array({self._name!r}, {self._value!r})"

    def __getitem__(self, index):
        if self._status.tainted:
            self._refresh()
        if _recorder is not None:
            _recorder.read(self)
        value = self._value[index]
        if _is_scalar_index(index, self._value.ndim) or not _is_basic_index(index):
            # scalars, and the copies selected by masks or arrays of indices, are
            # returned as they are: `obj.x_[i] op= v` then writes them back
            return value
        wvalue = value if self._wvalue is self._value else self._wvalue[index]
        return self._view(value, wvalue, index)

    def __setitem__(self, index, value):
        status = self._status
        if isinstance(value, depend_base):
            if getattr(value, "_parent", None) is self and value._index is index:
                # `obj.x_[i] op= v` stores back the view obtained from
                # `obj.x_[i]`, which has already been updated in place
                if status.synchro is None and self._func is not None:
                    self.update_man()  # raises, as the value is computed
                return
            value = value._value
        with status.lock:
            self._assign(index, value)
            if status.synchro is not None or self._func is not None:
                self.update_man()
            self.taint(taintme=False)

    def _assign(self, index, value):
        self._wvalue[index] = value

    def _apply(self, op, other):
        op(self._wvalue, other)

    def set(self, value, manual=True):
        """Full-array assignment. `manual=False` is used by `update_auto`.

        Writes into the existing buffer, so that views stay valid.
        """

        if isinstance(value, depend_base):
            value = value._value
        status = self._status
        with status.lock:
            self._assign(Ellipsis, value)
            if manual and (status.synchro is not None or self._func is not None):
                self.update_man()
            self.taint(taintme=False)

    def reshape(self, *shape):
        """Returns a reshaped view. Accepts `reshape((2, 3))` and `reshape(2, 3)`."""

        if len(shape) == 1:
            shape = shape[0]
        return self._view(self._value.reshape(shape), self._wvalue.reshape(shape))

    def flatten(self):
        """Returns a 1D view."""

        return self.reshape(-1)

    def _refresh(self):
        """Refreshes a view by recomputing the root array, which owns the full data."""

        depend_base._refresh(self._root())

    def _root(self):
        dep = self
        while dep._parent is not None:
            dep = dep._parent
        return dep

    def _load(self):
        return self._wvalue


def _make_inplace(op):
    def method(self, other):
        if isinstance(other, depend_base):
            other = other._value
        status = self._status
        with status.lock:
            if status.tainted:
                self._refresh()
            self._apply(op, other)
            if status.synchro is not None or self._func is not None:
                self.update_man()
            self.taint(taintme=False)
        return self

    method.__name__ = op.__name__
    return method


for _name, _fn in (
    ("__iadd__", _op.iadd),
    ("__isub__", _op.isub),
    ("__imul__", _op.imul),
    ("__itruediv__", _op.itruediv),
):
    setattr(depend_array, _name, _make_inplace(_fn))


class _reshaped:
    """Index of the views obtained by reshaping an immutable array."""

    __slots__ = ("shape",)

    def __init__(self, shape):
        self.shape = shape


class _static:
    """Hashable wrapper of a basic index, given to jitted functions as a static
    argument (slices are hashable only from Python 3.12)."""

    __slots__ = ("index", "key")

    def __init__(self, index):
        self.index = index
        self.key = repr(index)

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return self.key == other.key


@functools.cache
def _jax_indexing():
    """Jitted functions that read, update and fill JAX arrays, given basic
    indices, shapes or types (as JAX is much slower at indexing eagerly)."""

    import jax
    import jax.numpy as jnp

    def put(a, index, v):
        # like numpy, ignores the leading dimensions of size one of v
        v = jnp.asarray(v)
        v = v.reshape(v.shape[max(v.ndim - a[index].ndim, 0) :])
        return a.at[index].set(v)

    def put_static(a, index, v):
        return put(a, index.index, v)

    def fill(v, shape, dtype):
        return jnp.broadcast_to(jnp.asarray(v, dtype), shape)

    return types.SimpleNamespace(
        take=jax.jit(lambda a, index: a[index.index], static_argnums=1),
        put=jax.jit(put_static, static_argnums=1),
        put_eagerly=put,
        reshape=jax.jit(lambda a, shape: a.reshape(shape), static_argnums=1),
        fill=jax.jit(fill, static_argnums=(1, 2)),
    )


def _take(value, index):
    """Returns the part of `value` selected by `index`."""

    jx = _jax_indexing()
    if type(index) is _reshaped:
        return jx.reshape(value, index.shape)
    if _is_basic_index(index):
        return jx.take(value, _static(index))
    return value[index]


def _put(value, index, part):
    """Returns `value` with the part selected by `index` replaced by `part`."""

    jx = _jax_indexing()
    if type(index) is _reshaped:
        return jx.reshape(part, value.shape)
    if _is_basic_index(index):
        return jx.put(value, _static(index), part)
    return jx.put_eagerly(value, index, part)


def _fill(like, part):
    """Returns an array with the shape and type of `like`, whose values are not
    read, filled with `part`."""

    if (
        array_api_compat.is_jax_array(part)
        and part.shape == like.shape
        and part.dtype == like.dtype
    ):
        return part  # arrays are immutable, so it can be shared
    return _jax_indexing().fill(part, like.shape, like.dtype)


class depend_immutable(depend_array):
    """depend_array storing an array that cannot be modified in place (JAX).

    Writes replace the stored array with an updated copy. Views do not share
    memory with their parent: they read their value from the parent at each
    access, and write by replacing the value of the parent.
    """

    __slots__ = ()

    @property
    def _value(self):
        if self._parent is None:
            # also the reads made by the writes, which use the stored array
            if _recorder is not None:
                _recorder.read(self)
            return self._wvalue
        return _take(self._parent._value, self._index)

    @_value.setter
    def _value(self, value):
        pass  # held by `_wvalue`, or by the parent of a view

    def _store(self, value):
        if self._parent is None:
            self._wvalue = value
        else:
            parent = self._parent
            parent._store(_put(parent._value, self._index, value))

    def _load(self):
        return self._value

    def _assign(self, index, value):
        if index is Ellipsis and self._parent is None:
            # replaces the stored array without reading it
            self._store(_fill(self._wvalue, value))
        else:
            self._store(_put(self._value, index, value))

    def _apply(self, op, other):
        self._store(op(self._value, other))

    def _view(self, index):
        view = depend_immutable.__new__(depend_immutable)
        view._status = self._status
        view._func = None
        view._name = self._name
        view._parent = self
        view._index = index
        return view

    def __getitem__(self, index):
        if self._status.tainted:
            self._refresh()
        value = self._value
        if _is_scalar_index(index, value.ndim) or not _is_basic_index(index):
            # scalars, and the copies selected by masks or arrays of indices
            return _take(value, index)
        return self._view(index)

    def reshape(self, *shape):
        return self._view(_reshaped(shape[0] if len(shape) == 1 else shape))

    def __deepcopy__(self, memo):
        return depend_base.__deepcopy__(self, memo)

    def __setstate__(self, state):
        depend_base.__setstate__(self, state)


def dpipe(dfrom, dto, item=-1):
    """Synchonizes two depend objects.

    Takes two depend objects, and makes one of them depend on the other in such
    a way that both keep the same value. Used for attributes such as
    temperature that are used in many different modules, and so need different
    depend objects in each, but which should all have the same value.

    Args:
        dfrom: The object that is depend on.
        dto: The object that depends on the former one. If it already follows
            another object, it follows `dfrom` instead (e.g. the time step of
            the normal modes, which follows the one of the motion until a
            dynamics is bound).
        item: If non-negative, `dto` takes the value of `dfrom[item]`.
    """

    dto._func = None
    if item < 0:
        dcompute(dto, dfrom.get, [dfrom])
    else:
        dcompute(dto, lambda: dfrom.get()[item], [dfrom])


def dcompute(dep, func, dependencies):
    """Makes a depend object computed, as if it had been created with `func`
    and `dependencies`. Its current value is discarded, and it can no longer
    be written.

    Args:
        dep: The depend object. It must not be computed already, be
            synchronized with another object, or be a view of an array.
        func: Function computing the value.
        dependencies: The depend objects the value is computed from.
    """

    assert (
        dep._func is None and dep._status.synchro is None and dep._root() is dep
    ), f"<{dep._name}> is computed, synchronized, or a view"
    dep._func = func
    _add_dependencies(dep, dependencies)
    dep.taint(taintme=True)


def ddeps(dep, dependencies):
    """Adds dependencies to a computed depend object, e.g. ones that are only
    known when it is bound, so that it is also recomputed when they change.

    Args:
        dep: The depend object. It must be computed.
        dependencies: The depend objects to add to its dependencies.
    """

    assert dep._func is not None, f"<{dep._name}> is not computed"
    _add_dependencies(dep, dependencies)
    dep.taint(taintme=True)


def dsynchro(first, second, first_from_second, second_from_first):
    """Links two depend objects holding two representations of the same data.

    Writing one of the two marks the other as stale, and reading a stale one
    recomputes it from the other. Used e.g. for the cartesian and normal-mode
    coordinates of the beads, that can both be modified by the integrators.

    Args:
        first: The depend object that holds the current data.
        second: The depend object that is computed from `first` when read.
        first_from_second: Function computing the value of `first` from `second`.
        second_from_first: Function computing the value of `second` from `first`.
    """

    assert (
        first._status.synchro is None and second._status.synchro is None
    ), "Objects can only be synchronized with one peer"
    synchro = depend_synchro(first, second)
    first._func = first_from_second
    second._func = second_from_first
    first._status.synchro = synchro
    second._status.synchro = synchro
    first._taint_synchro()


def _inject_depend_property(cls, attr_name):
    """Adds a property `<attr_name>` to `cls` that reads and writes the value of
    the depend member `<attr_name>_`."""

    member_name = f"{attr_name}_"

    def getter(self):
        dep = self.__dict__[member_name]
        if dep._status.tainted:
            dep._refresh()
        if _recorder is not None:
            _recorder.read(dep)
        return dep._value

    def setter(self, value):
        self.__dict__[member_name].set(value)

    setattr(cls, attr_name, property(getter, setter))


def dproperties(cls, attr_names):
    """Expose depend members `<name>_` as attributes `<name>` on `cls`.

    Example: `dproperties(A, ["val"])` on a class with `self.val_ =
    depend_value(...)` makes `obj.val` read/write the value.
    """
    if not isinstance(attr_names, list):
        attr_names = [attr_names]
    for attr_name in attr_names:
        _inject_depend_property(cls, attr_name)
