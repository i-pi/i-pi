"""Graph of the depend objects of a simulation, i.e. of which quantities are
computed from which, written in the dot format of Graphviz.

`DependGraph(root)` collects the depend objects that can be reached from
`root` through attributes, lists and dicts, and `DependGraph.dot` writes
them as a graph:

- each node is a depend object, named by the path that leads to it from
  `root` (e.g. `syslist[0].beads.q`), and grouped with the other members of
  the object that holds it. Views share the node of the array they belong
  to. Objects that cannot be reached, but are computed from one that can,
  are named by their `_name` in angle brackets.
- an arrow goes from each object to the objects computed from it.
- a dashed line links the two objects of a synchronized pair.
- primitive objects are boxes, computed objects ellipses, and synchronized
  objects rounded boxes.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2017 i-PI developers
# See the "licenses" directory for full license information.

import collections

from ipi.utils.depend import depend_base

__all__ = ["DependGraph"]

_SHAPES = {
    "primitive": "shape=box",
    "computed": "shape=ellipse",
    "synchronized": 'shape=box, style="rounded"',
}


class DependGraph:
    """Depend objects reachable from an object, and their dependencies.

    Attributes:
        names: The path of each node from the root, by status (an object and
            its views share their status, so they are one node).
        edges: The pairs (dependency, dependant) of statuses.
        pairs: The pairs of statuses of synchronized objects.
    """

    def __init__(self, root):
        self.names = {}
        self.edges = set()
        self.pairs = set()
        self._collect(root)
        self._connect()

    def _collect(self, root):
        """Names the depend objects found walking breadth-first from `root`,
        so that each is named by its shortest path, and an array by its own
        path rather than by that of one of its views."""

        view_names = {}
        seen = {id(root)}
        queue = collections.deque([(root, "")])
        while queue:
            obj, path = queue.popleft()
            for name, value in _members(obj, path):
                if id(value) in seen:
                    continue
                seen.add(id(value))
                if isinstance(value, depend_base):
                    names = view_names if value._is_view else self.names
                    names.setdefault(value._status, name)
                elif _holds_members(value):
                    queue.append((value, name))
        for status, name in view_names.items():
            self.names.setdefault(status, name)

    def _connect(self):
        """Collects the dependencies and the synchronized pairs of the nodes,
        adding the objects computed from them that were not reached."""

        queue = list(self.names)
        while queue:
            status = queue.pop()
            # dependants that have been deleted leave a dead reference
            linked = [ref() for ref in status.dependants]
            linked = [dep for dep in linked if dep is not None]
            for dep in linked:
                self.edges.add((status, dep._status))
            if status.synchro is not None:
                first, second = status.synchro.peers
                self.pairs.add((first._status, second._status))
                linked += [first, second]
            for dep in linked:
                if dep._status not in self.names:
                    self.names[dep._status] = f"<{dep._name}>"
                    queue.append(dep._status)

    def around(self, patterns, depth=None):
        """Returns the nodes that match one of `patterns`, and those they are
        computed from or that are computed from them (through at most `depth`
        links, or any number if None). A node matches a pattern if its path
        is the pattern or ends with `.pattern`, e.g. `beads.q`."""

        selected = set()
        for pattern in patterns:
            matches = {s for s, name in self.names.items() if _matches(name, [pattern])}
            if not matches:
                raise ValueError(f"No depend object matches '{pattern}'")
            selected |= matches
        forward = collections.defaultdict(set)
        backward = collections.defaultdict(set)
        for a, b in self.edges:
            forward[a].add(b)
            backward[b].add(a)
        for a, b in self.pairs:
            # each object of a pair is computed from the other
            for links in (forward, backward):
                links[a].add(b)
                links[b].add(a)
        kept = set()
        for links in (forward, backward):
            reached, front, steps = set(selected), selected, 0
            while front and (depth is None or steps < depth):
                front = {b for a in front for b in links[a]} - reached
                reached |= front
                steps += 1
            kept |= reached
        return kept

    def dot(self, around=(), depth=None):
        """Returns the graph in the dot format of Graphviz, restricted to the
        nodes given by `self.around(around, depth)` if `around` is given."""

        nodes = self.around(around, depth) if around else set(self.names)
        ids = {
            status: f"n{i}"
            for i, status in enumerate(
                sorted(nodes, key=lambda status: self.names[status])
            )
        }
        # nodes grouped by the path of the object that holds them
        clusters = collections.defaultdict(list)
        for status in ids:
            owner, _, label = self.names[status].rpartition(".")
            style = _SHAPES[_kind(status)]
            if around and _matches(self.names[status], around):
                style += ", penwidth=3"
            clusters[owner].append(f"{ids[status]} [label={_quote(label)}, {style}];")

        lines = ["digraph depend {", "  rankdir=LR;", "  node [fontsize=10];"]
        for i, (owner, members) in enumerate(sorted(clusters.items())):
            if owner:
                lines += [f"  subgraph cluster{i} {{", f"    label={_quote(owner)};"]
                lines += [f"    {member}" for member in members] + ["  }"]
            else:
                lines += [f"  {member}" for member in members]
        lines += sorted(
            f"  {ids[a]} -> {ids[b]};" for a, b in self.edges if a in ids and b in ids
        )
        lines += sorted(
            f"  {ids[a]} -> {ids[b]} [style=dashed, dir=none];"
            for a, b in self.pairs
            if a in ids and b in ids
        )
        return "\n".join(lines + ["}"]) + "\n"


def _members(obj, path):
    """The pairs (path, value) of the items of a list, tuple or dict, or of
    the attributes of an object. Depend members are named without their
    trailing underscore, as their values (`beads.q` for `beads.q_`)."""

    if isinstance(obj, dict):
        return [(f"{path}[{key!r}]", value) for key, value in obj.items()]
    if isinstance(obj, (list, tuple)):
        return [(f"{path}[{i}]", value) for i, value in enumerate(obj)]
    prefix = path + "." if path else ""
    members = []
    for key, value in vars(obj).items():
        if isinstance(value, depend_base) and key.endswith("_"):
            key = key[:-1]
        members.append((prefix + key, value))
    return members


def _holds_members(value):
    """True for the containers and the i-PI objects that are walked."""

    if isinstance(value, (list, tuple, dict)):
        return True
    return type(value).__module__.startswith("ipi.") and hasattr(value, "__dict__")


def _kind(status):
    """The kind of the depend objects with the given status."""

    if status.func is not None:
        return "computed"
    if status.synchro is not None:
        return "synchronized"
    return "primitive"


def _matches(name, patterns):
    """True if the path `name` is one of `patterns`, or ends with `.pattern`."""

    return any(name == p or name.endswith("." + p) for p in patterns)


def _quote(text):
    """`text` as a quoted string of the dot format."""

    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'
