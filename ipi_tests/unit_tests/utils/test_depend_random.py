"""Randomized test of the depend objects.

Random sequences of operations (views, writes, reads, in-place operators,
dsetfunc, dpipe, copies) are run on a small graph of depend objects, and on a
model of the same graph made of numpy arrays, in which the computed values are
recomputed after each write. The values read from the depend objects must be
those of the model, and the writes to computed objects, or that numpy refuses,
must raise and change nothing.

A view is modelled by the flat indices of the elements of its root that it
selects, so that writing the view writes these elements. A failing sequence is
shortened before it is reported, as the list of the operations that reproduce
the failure.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import operator
import pickle
from copy import deepcopy

import numpy as np
import pytest

import ipi.utils.array_backend as ab
import ipi.utils.depend as dp

SHAPE = (4, 3)
INPLACE = {"+=": operator.iadd, "-=": operator.isub, "*=": operator.imul}
# relative frequencies of the operations
WEIGHTS = {
    "view": 4,
    "setitem": 4,
    "set": 3,
    "iop": 2,
    "index_iop": 2,
    "bad_iop": 0.5,
    "rhs_dep": 1,
    "value_set": 1,
    "read": 4,
    "dsetfunc": 0.2,
    "dpipe": 0.2,
    "copy": 0.4,
}
# seeds and lengths of the sequences that are run (JAX compiles a function for
# each new index, so it is slow at indexing eagerly)
SEQUENCES = {"numpy": (range(200), 60), "torch": (range(40), 60), "jax": (range(4), 40)}


class Failure(Exception):
    """A difference between the depend objects and the model."""


def flip(x):
    return ab.xp.flip(x, axis=0)


class Graph:
    """The depend objects: primitive (p, s, d, e, and the integers i),
    computed (c, t, k), and synchronized (a, b, whose conversions reverse the
    order of the rows). `d` and `e` are made computed by `dsetfunc` and `dpipe`
    during the sequences."""

    def __init__(self, init):
        xp = ab.xp
        self.p_ = dp.depend_array(name="p", value=xp.asarray(init["p"]))
        self.s_ = dp.depend_value(name="s", value=init["s"])
        self.c_ = dp.depend_array(
            name="c",
            value=xp.zeros(SHAPE),
            func=lambda: 2 * self.p + self.s,
            dependencies=[self.p_, self.s_],
        )
        self.a_ = dp.depend_array(name="a", value=xp.asarray(init["a"]))
        self.b_ = dp.depend_array(name="b", value=xp.zeros(SHAPE))
        dp.dsynchro(
            self.a_, self.b_, lambda: 0.5 * flip(self.b), lambda: 2 * flip(self.a)
        )
        self.t_ = dp.depend_value(
            name="t",
            func=lambda: self.a.sum() + self.s,
            dependencies=[self.a_, self.s_],
        )
        self.k_ = dp.depend_array(
            name="k",
            value=xp.zeros(4),
            func=lambda: self.b.sum(axis=1),
            dependencies=[self.b_],
        )
        self.d_ = dp.depend_array(name="d", value=xp.asarray(init["d"]))
        self.e_ = dp.depend_array(name="e", value=xp.asarray(init["e"]))
        self.i_ = dp.depend_array(name="i", value=xp.asarray(init["i"]))


PROPERTIES = ["p", "s", "c", "a", "b", "t", "k", "d", "e", "i"]
dp.dproperties(Graph, PROPERTIES)


class Root:
    """An object of the graph that is not a view, and its value in the model.

    Attributes:
        kind: "primitive", "computed" or "synchronized".
        func: Function computing the value in the model (from the peer, for a
            synchronized object).
        peer: Name of the synchronized peer.
        flat: The flat index of each element of the value.
    """

    def __init__(self, name, dep, value, kind="primitive", func=None, peer=None):
        self.name = name
        self.dep = dep
        self.value = np.array(value)
        self.kind = kind
        self.func = func
        self.peer = peer
        self.flat = np.arange(self.value.size).reshape(self.value.shape)
        self.root = self


class View:
    """A view, and the flat indices of the elements of its root it selects."""

    def __init__(self, name, dep, root, flat):
        self.name = name
        self.dep = dep
        self.root = root
        self.flat = flat


def is_int(obj):
    return obj.root.value.dtype.kind == "i"


def basic_index(rng, shape):
    """A random index made of integers, slices, Ellipsis and None."""

    def integer(n):
        return int(rng.integers(-n, n))

    def piece(n):
        start = None if rng.random() < 0.3 else int(rng.integers(0, n + 1))
        stop = None if rng.random() < 0.3 else int(rng.integers(0, n + 1))
        return slice(start, stop, [None, 1, 2][rng.integers(3)])

    choices = [lambda: Ellipsis, lambda: None]
    if len(shape) >= 1:
        choices += [lambda: piece(shape[0])] * 2
        if shape[0]:
            choices += [lambda: integer(shape[0]), lambda: (integer(shape[0]), None)]
        if shape[-1]:
            choices += [lambda: (Ellipsis, integer(shape[-1]))]
    if len(shape) >= 2:
        choices += [lambda: (piece(shape[0]), piece(shape[1]))]
        if shape[0]:
            choices += [lambda: (integer(shape[0]), piece(shape[1]))]
        if shape[1]:
            choices += [lambda: (piece(shape[0]), integer(shape[1]))]
        if shape[0] and shape[1]:
            choices += [lambda: (integer(shape[0]), integer(shape[1]))]
    return choices[rng.integers(len(choices))]()


def random_index(rng, shape):
    """A random basic index or, with probability 0.3, a mask, an array of
    indices or a 0-d array of the first dimension, or a mask of the second."""

    if not (len(shape) and shape[0] and rng.random() < 0.3):
        return basic_index(rng, shape)
    n = shape[0]
    kind = rng.integers(4 if len(shape) >= 2 else 3)
    if kind == 0:
        return rng.random(n) < 0.5
    if kind == 1:
        return np.sort(rng.permutation(n)[: rng.integers(0, n + 1)])
    if kind == 2:
        return np.asarray(rng.integers(0, n))
    return (slice(None), rng.random(shape[1]) < 0.5)


def is_basic(index):
    return not any(
        isinstance(i, np.ndarray) for i in (index if type(index) is tuple else (index,))
    )


def to_backend(index):
    """The index with its arrays converted to the array backend."""

    if type(index) is tuple:
        return tuple(to_backend(i) for i in index)
    return ab.xp.asarray(index) if isinstance(index, np.ndarray) else index


def rhs(rng, shape, integer):
    """A value to write in a part of the given shape, as given to the model
    and to the depend object: a scalar, or an array of the shape (or of its
    trailing dimensions)."""

    if rng.random() < 0.3:
        v = int(rng.integers(-4, 5))
        v = v if integer else float(v)
        return v, v
    if len(shape) and rng.random() < 0.2:
        shape = shape[1:]
    v = np.asarray(rng.integers(-4, 5, size=shape))
    v = v if integer else v.astype(float)
    return v, ab.xp.asarray(v)


def operand(rng, shape, integer, symbol):
    """The right-hand side of an in-place operator (for `*=`, factors that
    keep the values small)."""

    if symbol != "*=":
        return rhs(rng, shape, integer)
    choices = np.array([-1, 2] if integer else [-2.0, -1.0, -0.5, 0.5, 2.0])
    if rng.random() < 0.5:
        v = choices[rng.integers(len(choices))].item()
        return v, v
    v = choices[rng.integers(len(choices), size=shape)]
    return v, ab.xp.asarray(v)


class Model:
    """The graph of depend objects and its model, on which operations run."""

    def __init__(self, seed):
        rng = np.random.default_rng(seed)
        init = {k: rng.integers(-4, 5, size=SHAPE).astype(float) for k in "pad"}
        init["s"] = float(rng.integers(-3, 4))
        init["e"] = rng.integers(-4, 5, size=4).astype(float)
        init["i"] = rng.integers(-4, 5, size=5)
        self.graph = g = Graph(init)
        m = self.value
        self.roots = {
            "p": Root("p", g.p_, init["p"]),
            "s": Root("s", g.s_, init["s"]),
            "c": Root(
                "c", g.c_, np.zeros(SHAPE), "computed", lambda: 2 * m("p") + m("s")
            ),
            "a": Root(
                "a", g.a_, init["a"], "synchronized", lambda: 0.5 * m("b")[::-1], "b"
            ),
            "b": Root(
                "b",
                g.b_,
                np.zeros(SHAPE),
                "synchronized",
                lambda: 2 * m("a")[::-1],
                "a",
            ),
            "t": Root("t", g.t_, 0.0, "computed", lambda: m("a").sum() + m("s")),
            "k": Root("k", g.k_, np.zeros(4), "computed", lambda: m("b").sum(1)),
            "d": Root("d", g.d_, init["d"]),
            "e": Root("e", g.e_, init["e"]),
            "i": Root("i", g.i_, init["i"]),
        }
        self.views = {}
        # the operations run, reported on failures
        self.lines = []
        # dsetfunc and dpipe, which run once
        self.done = set()
        self.written(self.roots["a"])

    def value(self, name):
        return self.roots[name].value

    def written(self, root):
        """Updates the model after `root` has been written: its synchronized
        peer, then the computed objects, each after those it depends on."""

        if root.kind == "synchronized":
            peer = self.roots[root.peer]
            peer.value[...] = peer.func()
        for name in ("c", "t", "k", "d", "e"):
            computed = self.roots[name]
            if computed.kind == "computed":
                computed.value[...] = computed.func()

    def model_value(self, obj):
        return obj.root.value.reshape(-1)[obj.flat]

    def update(self, obj, index, func):
        """Replaces the elements of `obj[index]` in the model with `func` of
        their current values."""

        flat = obj.flat[index]
        values = obj.root.value.reshape(-1)
        values[flat] = np.broadcast_to(func(values[flat]), np.shape(flat))

    # checks

    def check(self, what, value, expected):
        value = np.asarray(ab.to_numpy(value))
        if value.shape != np.shape(expected) or not np.allclose(
            value, expected, rtol=1e-12, atol=1e-12
        ):
            raise Failure(
                f"{what} is {value.tolist()}, expected {np.asarray(expected).tolist()}"
            )

    def runs(self, what, func):
        """Calls `func`, which must not raise."""

        try:
            return func()
        except Exception as err:
            raise Failure(f"{what} raised {type(err).__name__}: {err}") from err

    def raises(self, what, func, errors=Exception):
        """Calls `func`, which must raise one of `errors`."""

        try:
            func()
        except errors:
            return
        except Exception as err:
            raise Failure(f"{what} raised {type(err).__name__}: {err}") from err
        raise Failure(f"{what} did not raise")

    def read_back(self, root):
        code = f"{root.name}.get()"
        self.lines.append(code)
        self.check(code, self.runs(code, root.dep.get), root.value)

    def write(self, rng, obj, code, run, func, index=Ellipsis):
        """Runs the write `run` of `obj[index]`, which must raise if `obj`
        belongs to a computed object, and otherwise replaces those elements in
        the model with `func` of their values. Then, half of the times, reads
        the root back."""

        root = obj.root
        self.lines.append(code)
        if root.kind == "computed":
            self.raises(code, run, NameError)
        else:
            self.runs(code, run)
            self.update(obj, index, func)
            self.written(root)
        if rng.random() < 0.5:
            self.read_back(root)

    # operations: each gets its random generator, the object it targets, the
    # other object it uses (or None) and its key

    def apply(self, op):
        kind, target, other, key = op
        kind_of, name = target
        obj = (self.roots if kind_of == "root" else self.views).get(name)
        # the targets made by operations removed by the shrinking do not exist
        if obj is not None:
            getattr(self, "op_" + kind)(np.random.default_rng(key), obj, other, key)

    def op_view(self, rng, obj, other, key):
        """Makes a view of an array or of a view, by indexing or reshaping it."""

        if isinstance(obj.dep, dp.depend_value):
            return
        name = f"v{len(self.views)}"
        if rng.random() < 0.25:
            size = obj.flat.size
            shapes = [(-1,), (size,), (1, size)]
            shapes += [(2, size // 2)] if size % 2 == 0 else []
            shape = shapes[rng.integers(len(shapes))]
            self.lines.append(f"{name} = {obj.name}.reshape({shape})")
            try:
                dep = obj.dep.reshape(shape)
            except (ValueError, RuntimeError):
                # a view of non-contiguous elements cannot be reshaped
                self.lines[-1] += "  # raised"
                return
            flat = obj.flat.reshape(shape)
        else:
            index = basic_index(rng, obj.flat.shape)
            code = f"{name} = {obj.name}[{index!r}]"
            self.lines.append(code)
            dep = self.runs(code, lambda: obj.dep[index])
            flat = obj.flat[index]
            if not isinstance(flat, np.ndarray):
                # a single element is a plain value
                if isinstance(dep, dp.depend_base):
                    raise Failure(f"{code} gives a depend object")
                self.check(code, dep, obj.root.value.reshape(-1)[flat])
                return
        if not isinstance(dep, dp.depend_array):
            raise Failure(f"{name} is a {type(dep).__name__}")
        self.views[key] = View(name, dep, obj.root, flat)

    def op_setitem(self, rng, obj, other, key):
        """`obj[index] = value`, with a basic or an advanced index."""

        if isinstance(obj.dep, dp.depend_value):
            return
        index = random_index(rng, obj.flat.shape)
        value, bvalue = rhs(rng, np.shape(obj.flat[index]), is_int(obj))
        bindex = to_backend(index)
        self.write(
            rng,
            obj,
            f"{obj.name}[{index!r}] = {value!r}",
            lambda: obj.dep.__setitem__(bindex, bvalue),
            lambda current: value,
            index,
        )

    def op_set(self, rng, obj, other, key):
        """Replaces the whole value of an array or a view, with `set`,
        `[...] =`, or the property of a root."""

        if isinstance(obj.dep, dp.depend_value):
            return
        value, bvalue = rhs(rng, obj.flat.shape, is_int(obj))
        how = rng.integers(3)
        if how == 2 and obj is obj.root and obj.name in PROPERTIES:
            code = f"g.{obj.name} = {value!r}"
            run = lambda: setattr(self.graph, obj.name, bvalue)  # noqa: E731
        elif how == 1:
            code = f"{obj.name}[...] = {value!r}"
            run = lambda: obj.dep.__setitem__(Ellipsis, bvalue)  # noqa: E731
        else:
            code = f"{obj.name}.set({value!r})"
            run = lambda: obj.dep.set(bvalue)  # noqa: E731
        self.write(rng, obj, code, run, lambda current: value)

    def op_iop(self, rng, obj, other, key):
        """`obj op= value`, on an array or a view."""

        if isinstance(obj.dep, dp.depend_value):
            return
        symbol = list(INPLACE)[rng.integers(3)]
        value, bvalue = operand(rng, obj.flat.shape, is_int(obj), symbol)
        self.write(
            rng,
            obj,
            f"{obj.name} {symbol} {value!r}",
            lambda: INPLACE[symbol](obj.dep, bvalue),
            lambda current: INPLACE[symbol](current, value),
        )

    def op_index_iop(self, rng, obj, other, key):
        """`obj[index] op= value`, run as Python runs it."""

        if isinstance(obj.dep, dp.depend_value):
            return
        index = random_index(rng, obj.flat.shape)
        symbol = list(INPLACE)[rng.integers(3)]
        value, bvalue = operand(rng, np.shape(obj.flat[index]), is_int(obj), symbol)
        bindex = to_backend(index)

        def run():
            part = obj.dep[bindex]
            part = INPLACE[symbol](part, bvalue)
            obj.dep[bindex] = part

        self.write(
            rng,
            obj,
            f"{obj.name}[{index!r}] {symbol} {value!r}",
            run,
            lambda current: INPLACE[symbol](current, value),
            index,
        )

    def op_bad_iop(self, rng, obj, other, key):
        """An in-place operator that numpy refuses, as it would change the
        shape or the type of the array: it must raise and change nothing."""

        if isinstance(obj.dep, dp.depend_value):
            return
        if is_int(obj) and rng.random() < 0.5:
            value, bvalue = 0.5, 0.5
        else:
            v = np.ones((2,) + obj.flat.shape, dtype=obj.root.value.dtype)
            value, bvalue = v, ab.xp.asarray(v)
        code = f"{obj.name} += {value!r}"
        self.lines.append(code)
        self.raises(code, lambda: obj.dep.__iadd__(bvalue))
        self.read_back(obj.root)

    def op_rhs_dep(self, rng, obj, other, key):
        """Writes the value of another depend object, which is read up to
        date (`s` if `other` is None)."""

        if isinstance(obj.dep, dp.depend_value):
            return
        source = self.roots["s"] if other is None else None
        if other is not None:
            kind_of, name = other
            source = (self.roots if kind_of == "root" else self.views).get(name)
        if source is None or source.root is obj.root or is_int(obj) or is_int(source):
            return
        if isinstance(source.dep, dp.depend_value):
            how = "+="
        elif source.flat.shape == obj.flat.shape:
            how = ["set", "[...] =", "+="][rng.integers(3)]
        else:
            return
        value = np.array(self.model_value(source))
        if how == "set":
            code, run = f"{obj.name}.set({source.name})", lambda: obj.dep.set(
                source.dep
            )
        elif how == "[...] =":
            code = f"{obj.name}[...] = {source.name}"
            run = lambda: obj.dep.__setitem__(Ellipsis, source.dep)  # noqa: E731
        else:
            code = f"{obj.name} += {source.name}"
            run = lambda: obj.dep.__iadd__(source.dep)  # noqa: E731
        if how == "+=":
            self.write(rng, obj, code, run, lambda current: current + value)
        else:
            self.write(rng, obj, code, run, lambda current: value)

    def op_value_set(self, rng, obj, other, key):
        """Writes a depend_value: `s` (primitive) or `t` (computed)."""

        root = self.roots[["s", "t"][rng.integers(2)]]
        value = float(rng.integers(-3, 4))
        if rng.random() < 0.5:
            code, run = f"{root.name}.set({value})", lambda: root.dep.set(value)
        else:
            code = f"g.{root.name} = {value}"
            run = lambda: setattr(self.graph, root.name, value)  # noqa: E731
        self.write(rng, root, code, run, lambda current: value)

    def op_read(self, rng, obj, other, key):
        """Reads an object, with `get`, a property, or an index (which gives a
        view, or a plain value)."""

        how = rng.integers(3)
        if how == 1 and obj is obj.root and obj.name in PROPERTIES:
            code = f"g.{obj.name}"
            result = self.runs(code, lambda: getattr(self.graph, obj.name))
            expected = obj.value
        elif how == 2 and isinstance(obj.dep, dp.depend_array):
            index = random_index(rng, obj.flat.shape)
            code = f"{obj.name}[{index!r}]"
            result = self.runs(code, lambda: obj.dep[to_backend(index)])
            flat = obj.flat[index]
            if is_basic(index) and isinstance(flat, np.ndarray):
                if not isinstance(result, dp.depend_array):
                    raise Failure(f"{code} gives a {type(result).__name__}")
                result = self.runs(code + ".get()", result.get)
            elif isinstance(result, dp.depend_base):
                raise Failure(f"{code} gives a depend object")
            expected = obj.root.value.reshape(-1)[flat]
        else:
            code = f"{obj.name}.get()"
            result = self.runs(code, obj.dep.get)
            expected = self.model_value(obj)
        self.lines.append(code)
        self.check(code, result, expected)

    def op_dsetfunc(self, rng, obj, other, key):
        """Makes `d` computed from `p` and `c` (once)."""

        if "dsetfunc" in self.done:
            return
        self.done.add("dsetfunc")
        g = self.graph
        code = "dsetfunc(d, lambda: g.p - g.c, [p, c])"
        self.lines.append(code)
        self.runs(code, lambda: dp.dsetfunc(g.d_, lambda: g.p - g.c, [g.p_, g.c_]))
        d = self.roots["d"]
        d.kind, d.func = "computed", lambda: self.value("p") - self.value("c")
        self.written(d)

    def op_dpipe(self, rng, obj, other, key):
        """Makes `e` follow `k` (once)."""

        if "dpipe" in self.done:
            return
        self.done.add("dpipe")
        code = "dpipe(k, e)"
        self.lines.append(code)
        self.runs(code, lambda: dp.dpipe(self.graph.k_, self.graph.e_))
        e = self.roots["e"]
        e.kind, e.func = "computed", lambda: self.value("k")
        self.written(e)

    def op_copy(self, rng, obj, other, key):
        """Copies an object with deepcopy or pickle: views cannot be copied,
        and the copy of an array is a primitive array with its current value."""

        copy = (
            deepcopy
            if rng.random() < 0.5
            else (lambda x: pickle.loads(pickle.dumps(x)))
        )
        if isinstance(obj, View):
            code = f"copy({obj.name})"
            self.lines.append(code)
            self.raises(code, lambda: copy(obj.dep), TypeError)
            return
        if obj.kind == "computed" or isinstance(obj.dep, dp.depend_value):
            return
        name = f"copy{len(self.roots)}"
        code = f"{name} = copy({obj.name})"
        self.lines.append(code)
        self.roots[key] = Root(name, self.runs(code, lambda: copy(obj.dep)), obj.value)
        self.read_back(self.roots[key])


# the targets of the operations (synchronized objects, which have the most
# states, twice as often)
ROOTS = PROPERTIES + ["a", "b"]


def generate(seed, n):
    """A random sequence of `n` operations (kind, target, other, key). The
    targets refer to roots by name, or to the views and copies made by
    earlier operations, by the key of the operation that made them."""

    rng = np.random.default_rng(seed)
    kinds = list(WEIGHTS)
    weights = np.array([WEIGHTS[k] for k in kinds])
    made, ops = [], []
    for _ in range(n):
        kind = kinds[rng.choice(len(kinds), p=weights / weights.sum())]
        key = int(rng.integers(2**62))

        def pick():
            if made and rng.random() < 0.5:
                return made[rng.integers(len(made))]
            return ("root", ROOTS[rng.integers(len(ROOTS))])

        target = pick()
        other = pick() if kind == "rhs_dep" and rng.random() < 0.8 else None
        ops.append((kind, target, other, key))
        if kind == "view":
            made.append(("view", key))
        elif kind == "copy":
            made.append(("root", key))
    return ops


def run(seed, ops):
    """Runs the operations on a new graph. Returns the first failure (with the
    index of its operation) or None, and the operations run."""

    model = Model(seed)
    for i, op in enumerate(ops):
        try:
            model.apply(op)
        except Failure as failure:
            failure.index = i
            return failure, model.lines
    return None, model.lines


def shrink(seed, ops, failure):
    """Removes operations from a failing sequence while it still fails."""

    ops = ops[: failure.index + 1]
    for chunk in (16, 8, 4, 2, 1):
        i = 0
        while i < len(ops):
            trial = ops[:i] + ops[i + chunk :]
            found, _ = run(seed, trial)
            if found is not None:
                ops, failure = trial[: found.index + 1], found
            else:
                i += chunk
    return ops, failure


@pytest.mark.parametrize("backend", ["numpy", "torch", "jax"])
def test_random_operations_match_the_model(select_backend, backend):
    """Depend: random sequences of operations on depend objects of all kinds
    give the values of a model made of numpy arrays recomputed eagerly"""
    select_backend(backend)
    seeds, n = SEQUENCES[backend]
    for seed in seeds:
        ops = generate(seed, n)
        failure, _ = run(seed, ops)
        if failure is not None:
            ops, failure = shrink(seed, ops, failure)
            _, lines = run(seed, ops)
            pytest.fail(f"seed {seed}: {failure}\n" + "\n".join(lines))
