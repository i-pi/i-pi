"""Indexing of JAX arrays with jitted functions.

JAX arrays cannot be modified in place, and JAX is much slower than numpy at
indexing them eagerly. The depend objects that hold JAX arrays
(`depend.depend_immutable`) read, update, reshape and fill them with the
functions of this module, which are jitted for each basic index (integers,
slices, Ellipsis and None), shape and type.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import functools
import types

import array_api_compat

__all__ = ["take", "put", "put_eagerly", "reshape", "fill"]


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
def _jitted():
    """The jitted functions, and the eager update of an array."""

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


def take(value, index):
    """Returns the part of `value` selected by the basic index `index`."""

    return _jitted().take(value, _static(index))


def put(value, index, part):
    """Returns `value` with the part selected by the basic index `index`
    replaced by `part`."""

    return _jitted().put(value, _static(index), part)


def put_eagerly(value, index, part):
    """Returns `value` with the part selected by any index `index` (e.g. a mask
    or an array of indices) replaced by `part`."""

    return _jitted().put_eagerly(value, index, part)


def reshape(value, shape):
    """Returns `value` with the shape `shape`."""

    return _jitted().reshape(value, shape)


def fill(like, part):
    """Returns an array with the shape and type of `like`, whose values are not
    read, filled with `part`."""

    if (
        array_api_compat.is_jax_array(part)
        and part.shape == like.shape
        and part.dtype == like.dtype
    ):
        return part  # arrays are immutable, so it can be shared
    return _jitted().fill(part, like.shape, like.dtype)
