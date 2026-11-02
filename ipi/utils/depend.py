"""Dependency tracking, lazy caching, and automatic update of variables.

i-PI stores physical quantities in depend objects, which know the objects
they are computed from, and recompute their value only when it is read after
one of those has changed.

Classes:

- `depend_value` holds any Python object (typically a scalar).
- `depend_array` holds an array. Indexing it with integers, slices, Ellipsis
  and None, or reshaping it, gives a `depend_array_view`: a depend object for
  part of the array, which stays connected to the dependency graph.
- `depend_immutable` and `depend_immutable_view` take their place for JAX
  arrays, which cannot be modified in place.

State. A depend object holds its value, and a `depend_status` holding
everything else: the function that computes the value, a tainted flag (the
value is out of date), the dependants to taint when the value changes, the
synchronization with a peer, and a lock. The views of an array share the
status of the array they belong to (their root), so they are always in the
same state as the root.

Kinds of objects:

- Primitive objects hold data, and can be written.
- Computed objects are created with `func`, or made computed later with
  `dsetfunc` (or `dpipe`, to follow the value of another object). Their value
  is a cache of the function, which cannot be written. `ddeps` gives them more
  dependencies.
- Synchronized objects come in pairs, linked by `dsynchro`, that hold two
  representations of the same data (e.g. the cartesian and normal-mode
  coordinates of the beads). Both can be written: the one written last holds
  the data, and the other is computed from it.

Reads: reading an object that is out of date recomputes it (or its root, for
a view), from its function or from its peer, and marks it up to date.

Writes: all the writes (`set`, `obj.x = v`, `obj.x_[i] = v`, in-place
operators, also on views) take the same steps, numbered as here in the code
of `depend_value.set`, `depend_array.set`, `depend_array.__setitem__` and
`_inplace`:

1. raise before changing anything if the object is computed;
2. bring the value up to date first, if the write keeps part of it (a
   partial index, an in-place operator, or any write to a view);
3. modify the stored value;
4. make a synchronized object the one that holds the data, and taint its
   peer;
5. taint the dependants, recursively.

`set(value, manual=False)` stores the value of a computed object (when it is
recomputed, or to restore a value computed earlier) with steps 3 and 5 only.

Copies: deep copies and pickles hold a copy of the value, and are isolated:
they have no dependants and no synchronization (a synchronized object that is
out of date is brought up to date first). Computed copies share the function
of the original. Views cannot be copied.

Usage: classes expose their depend members through `dproperties`. For a
member stored as `obj.x_`, reading `obj.x` returns the current value, and
`obj.x = v` replaces it. The depend object `obj.x_` is used to declare
dependencies, to modify an array in place (`obj.x_[i] = v`, `obj.x_ += v`), and
to obtain views (`obj.x_[i]`, `obj.x_.flatten()`).

Arrays are stored with the array backend selected in
`ipi.utils.array_backend`, and behave as numpy arrays with all the backends:

- `obj.x_[i]` is a view if `i` holds only integers, slices, Ellipsis and None,
  and does not select a single element. Elements, and the copies selected by
  masks or arrays of indices, are plain values: `obj.x_[mask] op= v` writes
  them back. As the values read, they must not be modified (with torch,
  elements share memory with the stored array).
- `reshape` and `flatten` give views, and raise if that needs a copy.
- In-place operators keep the shape and the type of the array (types are cast
  if numpy would cast them in place), or raise.
- The value read from `obj.x` (or with `get()`) must not be modified. With
  numpy it is read-only, so modifying it raises. It is valid until the object
  is written or recomputed: depending on the backend, later writes may or may
  not change it. `_load` returns the stored value without recomputing it
  (e.g. a stale value, to seed a new calculation).
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2017 i-PI developers
# See the "licenses" directory for full license information.


import operator
import threading
import weakref

import array_api_compat
import numpy as np

from ipi.utils import jaxindexing
from ipi.utils.array_backend import xp

__all__ = [
    "depend_value",
    "depend_array",
    "dsetfunc",
    "ddeps",
    "dpipe",
    "dsynchro",
    "dproperties",
]

# object notified of the reads and writes of depend objects while a function
# is recorded to be compiled (see ipi.utils.jit), None otherwise
_recorder = None


def set_recorder(recorder):
    """Sets the object notified of the operations on depend objects."""

    global _recorder
    _recorder = recorder


# --------------------------------------------------------------------------
# Depend objects
#
# The link between two synchronized objects, the status shared by an object
# and its views, and the base class that reads, recomputes and taints all
# depend objects. `depend_value` holds any Python object.
# --------------------------------------------------------------------------


class depend_synchro(object):
    """Link between two depend objects holding two representations of the
    same data, created by `dsynchro`.

    Attributes:
        peers: The two synchronized objects.
        funcs: The functions computing the value of each peer from the other.
        manual: Status of the peer that holds the data, i.e. the one that was
            written last.
    """

    __slots__ = ("peers", "funcs", "manual")

    def __init__(self, first, second, first_from_second, second_from_first):
        self.peers = (first, second)
        self.funcs = (first_from_second, second_from_first)
        self.manual = first._status

    def compute(self, status):
        """Returns the value of the peer with `status`, computed from the other."""

        return self.funcs[0 if self.peers[0]._status is status else 1]()


class depend_status(object):
    """State of a depend object other than its value, shared with its views.

    Attributes:
        func: Function computing the value of a computed object, None for the
            other objects.
        tainted: True if the value is out of date.
        dependants: Weak references to the objects to taint when the value
            changes.
        synchro: The `depend_synchro` linking the object to a peer, or None.
        lock: Lock making the updates of the value thread-safe.
    """

    __slots__ = ("func", "tainted", "dependants", "synchro", "lock")

    def __init__(self, func=None, tainted=False):
        self.func = func
        self.tainted = tainted
        self.dependants = []
        self.synchro = None
        self.lock = threading.RLock()

    # copies have no dependants and no synchronization. Deep copies share
    # `func`, which is often a bound method whose `__self__` holds state that
    # cannot be copied, such as a lock
    def __deepcopy__(self, memo):
        return depend_status(self.func, self.tainted)

    def __getstate__(self):
        return self.func, self.tainted

    def __setstate__(self, state):
        self.__init__(*state)


class depend_base(object):
    """Value, status and name of a depend object, and the reads shared by all
    the depend objects. The writes (`set`) are defined by the subclasses.

    Attributes:
        _status: The `depend_status` of the object, shared with its views.
        _name: Name of the object, used in error messages.
        _value: The value returned by reads.
    """

    __slots__ = ("_status", "_name", "_value", "__weakref__")

    # True for the views of arrays
    _is_view = False

    def __init__(self, name, func=None, dependencies=None):
        self._name = name
        # computed objects start out of date, primitive objects up to date
        self._status = depend_status(func, tainted=func is not None)
        if dependencies is not None:
            _add_dependencies(self, dependencies)

    # deep copies and pickles: Python's default protocol makes the copy with
    # `__getstate__` and `__setstate__`
    def __getstate__(self):
        status = self._status
        if status.tainted and status.func is None:
            # a synchronized object computed from its peer, which the copy
            # does not have
            self._refresh()
        return {
            m: getattr(self, m) for m in _slot_names(type(self)) if hasattr(self, m)
        }

    def __setstate__(self, state):
        for member, value in state.items():
            setattr(self, member, value)

    def tainted(self):
        """Returns True if the value is out of date."""

        return self._status.tainted

    def taint(self, taintme=True):
        """Taints the dependants, recursively, and sets the tainted flag of
        the object to `taintme`."""

        for item in self._status.dependants:
            item = item()
            if item is not None and not item._status.tainted:
                item.taint()
        self._status.tainted = taintme
        if _recorder is not None:
            _recorder.tainted(self, taintme)

    def _taint_synchro(self):
        """Makes the object the one of a synchronized pair that holds the
        data, and taints its peer (and the dependants of the peer)."""

        synchro = self._status.synchro
        synchro.manual = self._status
        for peer in synchro.peers:
            if peer._status is not self._status and not peer._status.tainted:
                peer.taint()
        if _recorder is not None:
            _recorder.synchronized(synchro)

    def get(self):
        """Returns the current value, recomputing it if needed."""

        if self._status.tainted:
            self._refresh()
        if _recorder is not None:
            _recorder.read(self)
        return self._value

    def _refresh(self):
        """Recomputes the value if it is out of date, and marks it up to date.

        Callers check `self._status.tainted` first, so that reading a value
        that is up to date does not call it. The check is repeated under the
        lock, so that only one thread recomputes the value.
        """

        status = self._status
        with status.lock:
            if status.tainted:
                synchro = status.synchro
                if synchro is not None and synchro.manual is not status:
                    self.set(synchro.compute(status), manual=False)
                elif synchro is None and status.func is not None:
                    self.set(status.func(), manual=False)
                else:
                    # nothing to recompute (e.g. the peer that holds the data)
                    self.taint(taintme=False)

    def _root(self):
        """Returns the object that stores the value: the object itself, or the
        array that a view belongs to."""

        return self

    def _load(self):
        """Returns the stored value, without recomputing it."""

        return self._value

    def _store(self, value):
        """Replaces the stored value, without tainting the dependants."""

        self._value = value


def _slot_names(cls):
    """Names of all the slots defined by `cls` and its bases, except
    `__weakref__` and those replaced by a property (the value of
    `depend_immutable`)."""

    return [
        s
        for c in cls.__mro__
        for s in getattr(c, "__slots__", ())
        if s != "__weakref__" and not isinstance(getattr(cls, s), property)
    ]


def _add_dependencies(dep, dependencies):
    """Makes `dep` a dependant of each object in `dependencies`."""

    for item in dependencies:
        item._status.dependants.append(weakref.ref(dep))


def _computed_error(dep):
    """The error raised when the computed object `dep` is written."""

    return NameError(
        "Cannot set manually the value of the automatically-computed property "
        f"<{dep._name}>"
    )


class depend_value(depend_base):
    """Depend object holding any Python object (typically a scalar)."""

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
        """Replaces the value, with the steps of the writes listed in the
        docstring of the module.

        Args:
            value: The new value.
            manual: False to store the value of a computed object, e.g. when
                it is recomputed, or to restore a value computed earlier.
        """

        status = self._status
        with status.lock:
            # 1. computed objects cannot be written
            if manual and status.func is not None:
                raise _computed_error(self)
            # 3. (a value is written whole, so there is no step 2)
            self._value = value
            # 4. the written peer holds the data
            if manual and status.synchro is not None:
                self._taint_synchro()
            # 5. taints the dependants
            self.taint(taintme=False)


# --------------------------------------------------------------------------
# Arrays and their views
#
# A `depend_array` holds a numpy or torch array, which writes modify in place.
# Basic indices and `reshape` give views (`depend_array_view`) that share the
# memory and the status of the array; elements, masks and arrays of indices
# give plain values.
# --------------------------------------------------------------------------


def _is_basic_index(index):
    """True if `array[index]` is a view of the array, i.e. if the index only
    holds integers, slices, Ellipsis and None."""

    for item in index if type(index) is tuple else (index,):
        if not (
            item is None
            or item is Ellipsis
            or type(item) in (int, slice)
            or isinstance(item, np.integer)
        ):
            return False
    return True


def _is_element_index(index, ndim):
    """True if the basic index `index` selects a single element of an array
    of rank `ndim`, i.e. if it is made of `ndim` integers."""

    items = index if type(index) is tuple else (index,)
    if len(items) != ndim:
        return False
    for item in items:
        if not (type(item) is int or isinstance(item, np.integer)):
            return False
    return True


def _readonly(array):
    """Returns a non-writeable view of a numpy array, and other arrays as
    they are."""

    if not isinstance(array, np.ndarray):
        return array
    view = array.view()
    view.flags.writeable = False
    return view


class _reshaped:
    """Index of the views obtained by reshaping an array."""

    __slots__ = ("shape",)

    def __init__(self, shape):
        self.shape = shape


def _inplace(op):
    """Returns the in-place operator of depend arrays that updates the array
    with `op`, with the steps of the writes listed in the docstring of the
    module."""

    def method(self, other):
        if isinstance(other, depend_base):
            other = other.get()
        status = self._status
        with status.lock:
            # 1. computed objects cannot be written
            if status.func is not None:
                raise _computed_error(self)
            # 2. the operator reads the value it updates
            if status.tainted:
                self._refresh()
            if _recorder is not None:
                _recorder.read(self)
            # 3. updates the stored array
            self._apply(op, other)
            # 4. the written peer holds the data
            if status.synchro is not None:
                self._taint_synchro()
            # 5. taints the dependants
            self.taint(taintme=False)
        return self

    method.__name__ = f"__{op.__name__}__"
    return method


class depend_array(depend_base):
    """Depend object holding an array.

    The array is stored in `_data`, with the selected array backend, or with
    numpy if `host` is True (for arrays of strings, or small arrays of
    metadata only used on the host). Writes modify `_data` in place, so that
    the views of the array stay valid.

    Reads return `_value`: with numpy, a non-writeable view of `_data`, so
    that modifying the value read raises. It is made when the array is stored
    rather than at each read, which would cost more than the read itself.
    With the other backends, which have no read-only arrays, `_value` is
    `_data`.

    Attributes:
        _data: The stored array.
    """

    __slots__ = ("_data",)

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
        data = np.asarray(value) if host else xp.asarray(value)
        if array_api_compat.is_jax_array(data):
            self.__class__ = depend_immutable
        super().__init__(name, func, dependencies)
        self._store(data)

    def __getstate__(self):
        state = super().__getstate__()
        # the value read is made again from `_data` by `__setstate__`
        state.pop("_value", None)
        return state

    def __setstate__(self, state):
        super().__setstate__(state)
        self._store(self._data)

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
        if _is_element_index(index, self._data.ndim) or not _is_basic_index(index):
            # elements, and the copies selected by masks or arrays of indices,
            # are plain values, which `obj.x_[i] op= v` writes back
            return self._value[index]
        data = self._data[index]
        # with numpy, indexing the read-only value gives a read-only view
        value = data if self._value is self._data else self._value[index]
        return depend_array_view(self, index, data, value)

    def set(self, value, manual=True):
        """Replaces the value, as `depend_value.set`. Writes into the stored
        array, so that the views stay valid."""

        if isinstance(value, depend_base):
            value = value.get()
        status = self._status
        with status.lock:
            if manual:
                # 1. computed objects cannot be written
                if status.func is not None:
                    raise _computed_error(self)
                # 2. a view is part of its root, whose other elements are kept
                if self._is_view:
                    if status.tainted:
                        self._refresh()
                    if _recorder is not None:
                        _recorder.read(self)
            # 3. writes the stored array
            self._assign(Ellipsis, value)
            # 4. the written peer holds the data
            if manual and status.synchro is not None:
                self._taint_synchro()
            # 5. taints the dependants
            self.taint(taintme=False)

    def __setitem__(self, index, value):
        if isinstance(value, depend_base):
            if (
                isinstance(value, depend_view)
                and value._parent is self
                and _is_basic_index(index)
                and value._index == index
            ):
                # `obj.x_[i] op= v` updates the view `obj.x_[i]` in place, and
                # then assigns it to `obj.x_[i]`: there is nothing left to do
                return
            value = value.get()
        status = self._status
        with status.lock:
            # 1. computed objects cannot be written. With torch, the elements
            # given by `obj.x_[i]` share memory with the stored array, so
            # `obj.x_[i] op= v` has already modified it: it is recomputed
            if status.func is not None:
                self.taint()
                raise _computed_error(self)
            # 2. the elements that are not assigned are kept
            if index is not Ellipsis or self._is_view:
                if status.tainted:
                    self._refresh()
                if _recorder is not None:
                    _recorder.read(self)
            # 3. writes the stored array
            self._assign(index, value)
            # 4. the written peer holds the data
            if status.synchro is not None:
                self._taint_synchro()
            # 5. taints the dependants
            self.taint(taintme=False)

    __iadd__ = _inplace(operator.iadd)
    __isub__ = _inplace(operator.isub)
    __imul__ = _inplace(operator.imul)
    __itruediv__ = _inplace(operator.itruediv)

    def _assign(self, index, value):
        self._data[index] = value

    def _apply(self, op, other):
        # an in-place operator that falls back to an out-of-place one (e.g. for
        # a torch tensor and a numpy array) leaves the stored array unchanged
        if op(self._data, other) is not self._data:
            raise TypeError(
                f"Cannot update <{self._name}> in place with a "
                f"{type(other).__name__}: use an array of the backend"
            )

    def reshape(self, *shape):
        """Returns a view with another shape, given as `reshape((2, 3))` or
        `reshape(2, 3)`. Raises if the array cannot be reshaped without a copy
        (e.g. a view of non-contiguous elements that is flattened), as writing
        the copy would not modify the array."""

        if len(shape) == 1:
            shape = shape[0]
        if isinstance(self._data, np.ndarray):
            data = self._data.reshape(shape, copy=False)
        else:
            # torch tensors: `view` raises if a copy is needed
            data = self._data.view(shape)
        return depend_array_view(self, _reshaped(shape), data, _readonly(data))

    def flatten(self):
        """Returns a 1D view (see `reshape`)."""

        return self.reshape(-1)

    def _load(self):
        return self._data

    def _store(self, data):
        self._data = data
        self._value = _readonly(data)


class depend_view(object):
    """Methods shared by the views of depend arrays.

    A view is part of the array it belongs to (its root), and has no state of
    its own: it shares the status of the root, and reading it when it is out
    of date recomputes the root. Writing a view writes part of the root, so
    the root is brought up to date first, and the dependants of the root are
    tainted. As for the root, a view of a computed array cannot be written.
    Views cannot be copied, made computed or synchronized: copy or modify the
    root instead.

    The subclasses define the slots:
        _parent: The array or view that the view was obtained from.
        _index: The index that selected the view from `_parent`.
    """

    __slots__ = ()

    _is_view = True

    def _attach(self, parent, index):
        """Makes the view the part of `parent` selected by `index`."""

        self._status = parent._status
        self._name = parent._name
        self._parent = parent
        self._index = index

    def __getstate__(self):
        raise TypeError(
            f"<{self._name}> is a view, which cannot be copied: copy the array "
            "it belongs to"
        )

    def _root(self):
        return self._parent._root()

    def _refresh(self):
        self._root()._refresh()


class depend_array_view(depend_view, depend_array):
    """Part of a depend_array (see `depend_view`), whose `_data` shares memory
    with the array. `value` is the value read, as in `depend_array._store`."""

    __slots__ = ("_parent", "_index")

    def __init__(self, parent, index, data, value):
        self._attach(parent, index)
        self._data = data
        self._value = value


# --------------------------------------------------------------------------
# JAX arrays
#
# JAX arrays cannot be modified in place, so they are held by a
# `depend_immutable`, whose writes replace the stored array with an updated
# copy, and whose views read and write through their parent. `_jax_take` and
# `_jax_put` read and update JAX arrays with the jitted functions of
# `ipi.utils.jaxindexing`.
# --------------------------------------------------------------------------


def _jax_take(value, index):
    """Returns the part of `value` selected by `index`."""

    if type(index) is _reshaped:
        return jaxindexing.reshape(value, index.shape)
    if _is_basic_index(index):
        return jaxindexing.take(value, index)
    return value[index]


def _jax_put(value, index, part):
    """Returns `value` with the part selected by `index` replaced by `part`."""

    if type(index) is _reshaped:
        return jaxindexing.reshape(part, value.shape)
    if _is_basic_index(index):
        return jaxindexing.put(value, index, part)
    return jaxindexing.put_eagerly(value, index, part)


class depend_immutable(depend_array):
    """depend_array holding a JAX array, which cannot be modified in place.

    Writes replace `_data` with an updated copy. Reads return `_data` itself,
    as JAX arrays are read-only. The views (`depend_immutable_view`) cannot
    share memory with the array: they read and write through their parent.
    """

    __slots__ = ()

    @property
    def _value(self):
        # notifies all the reads of the stored array, also those made by
        # writes and by the views
        if _recorder is not None:
            _recorder.read(self)
        return self._data

    def __getitem__(self, index):
        if self._status.tainted:
            self._refresh()
        value = self._value
        if _is_element_index(index, value.ndim) or not _is_basic_index(index):
            return _jax_take(value, index)
        return depend_immutable_view(self, index)

    def reshape(self, *shape):
        return depend_immutable_view(
            self, _reshaped(shape[0] if len(shape) == 1 else shape)
        )

    def _assign(self, index, value):
        if index is Ellipsis:
            # replaces the array without reading it
            self._store(jaxindexing.fill(self._data, value))
        else:
            self._store(_jax_put(self._value, index, value))

    def _apply(self, op, other):
        value = self._value
        result = op(value, other)
        # as numpy does in place, keeps the shape, and casts the type
        if result.shape != value.shape or not np.can_cast(
            result.dtype, value.dtype, "same_kind"
        ):
            raise TypeError(
                f"Cannot update <{self._name}> in place with a result of shape "
                f"{result.shape} and type {result.dtype}"
            )
        if result.dtype != value.dtype:
            result = result.astype(value.dtype)
        self._store(result)

    def _load(self):
        return self._value

    def _store(self, data):
        self._data = data


class depend_immutable_view(depend_view, depend_immutable):
    """Part of a depend_immutable (see `depend_view`). It holds no array: its
    value is taken from the value of its parent at each read, and writing it
    replaces the value of its parent."""

    __slots__ = ("_parent", "_index")

    def __init__(self, parent, index):
        self._attach(parent, index)

    @property
    def _value(self):
        return _jax_take(self._parent._value, self._index)

    def _assign(self, index, value):
        self._store(_jax_put(self._value, index, value))

    def _store(self, value):
        parent = self._parent
        parent._store(_jax_put(parent._value, self._index, value))


# --------------------------------------------------------------------------
# Functions acting on depend objects
#
# `dsetfunc`, `dpipe`, `ddeps` and `dsynchro` change how existing depend
# objects are computed from each other, e.g. when a module is bound to
# another. `dproperties` exposes the values of the depend members of a class
# as attributes.
# --------------------------------------------------------------------------


def dpipe(dfrom, dto, item=-1):
    """Makes a depend object follow the value of another.

    Used for quantities such as the temperature, that are used in many
    modules, each with its own depend object, which should all have the same
    value.

    Args:
        dfrom: The object to follow.
        dto: The object that follows `dfrom`. If it already follows another
            object, it follows `dfrom` instead (e.g. the time step of the
            normal modes follows the one of the motion, until a dynamics is
            bound to them).
        item: If non-negative, `dto` takes the value of `dfrom[item]`.
    """

    status = dto._status
    assert (
        status.synchro is None and not dto._is_view
    ), f"<{dto._name}> is synchronized, or a view"
    status.func = None
    if item < 0:
        dsetfunc(dto, dfrom.get, [dfrom])
    else:
        dsetfunc(dto, lambda: dfrom.get()[item], [dfrom])


def dsetfunc(dep, func, dependencies):
    """Sets the function computing a depend object, as if it had been created
    with `func` and `dependencies`. Its current value is discarded, and it can
    no longer be written.

    Args:
        dep: The depend object. It must not be computed already, be
            synchronized with another object, or be a view.
        func: Function computing the value.
        dependencies: The depend objects the value is computed from.
    """

    status = dep._status
    assert (
        status.func is None and status.synchro is None and not dep._is_view
    ), f"<{dep._name}> is computed, synchronized, or a view"
    status.func = func
    _add_dependencies(dep, dependencies)
    dep.taint(taintme=True)


def ddeps(dep, dependencies):
    """Adds dependencies to a computed depend object, e.g. ones that are only
    known when it is bound, so that it is also recomputed when they change.

    Args:
        dep: The depend object. It must be computed, and not be a view.
        dependencies: The depend objects to add to its dependencies.
    """

    assert (
        dep._status.func is not None and not dep._is_view
    ), f"<{dep._name}> is not computed, or is a view"
    _add_dependencies(dep, dependencies)
    dep.taint(taintme=True)


def dsynchro(first, second, first_from_second, second_from_first):
    """Links two depend objects holding two representations of the same data.

    Writing one of the two marks the other as out of date, and reading the
    other then recomputes it. Used e.g. for the cartesian and normal-mode
    coordinates of the beads, which can both be modified by the integrators.

    Args:
        first: The object that holds the data.
        second: The object computed from `first` when it is read.
        first_from_second: Function computing the value of `first` from `second`.
        second_from_first: Function computing the value of `second` from `first`.
    """

    for dep in (first, second):
        status = dep._status
        assert (
            status.func is None and status.synchro is None and not dep._is_view
        ), f"<{dep._name}> is computed, synchronized, or a view"
    synchro = depend_synchro(first, second, first_from_second, second_from_first)
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
    """Exposes the depend members `<name>_` as attributes `<name>` of `cls`.

    Example: `dproperties(A, ["val"])` on a class with `self.val_ =
    depend_value(...)` makes `obj.val` read and write the value.
    """

    if not isinstance(attr_names, list):
        attr_names = [attr_names]
    for attr_name in attr_names:
        _inject_depend_property(cls, attr_name)
