"""DependGraph and i-pi-depend-graph: the graph of the depend objects of a
simulation."""

import runpy
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

from ipi.engine.simulation import Simulation
from ipi.utils import depend as dp
from ipi.utils.tools import DependGraph
from ipi_tests.unit_tests.engine.inputs import harmonic_input

SCRIPT = Path(__file__).parents[3] / "tools" / "py" / "depend-graph.py"


class Holder:
    """Object holding the depend objects of a small graph."""


def small_graph():
    """A graph with views, a synchronized pair, an array reached only through
    a view, a computed object that cannot be reached from the root, and a
    deleted dependant. Returns the root, and the unreachable object, which
    must be kept alive."""

    root = Holder()
    root.a_ = dp.depend_array(name="a", value=np.zeros(4))
    root.b_ = dp.depend_value(
        name="b", func=lambda: root.a_.get().sum(), dependencies=[root.a_]
    )
    root.parts = [root.a_[:2], {"tail": root.a_[2:]}]
    root.x_ = dp.depend_array(name="x", value=np.zeros(2))
    root.y_ = dp.depend_array(name="y", value=np.zeros(2))
    dp.dsynchro(root.x_, root.y_, lambda: root.y_.get() / 2, lambda: root.x_.get() * 2)
    root.only_views = (dp.depend_array(name="v", value=np.zeros(3))[1:],)
    hidden = dp.depend_value(
        name="hidden", func=lambda: root.b_.get() + 1, dependencies=[root.b_]
    )
    dp.depend_value(name="gone", func=lambda: 0, dependencies=[root.a_])
    return root, hidden


def test_nodes_edges_and_pairs():
    """DependGraph: views share the node of their array, which is named by its
    own path, and unreachable dependants are named by their _name"""

    root, hidden = small_graph()
    graph = DependGraph(root)
    status = {name: s for s, name in graph.names.items()}

    assert sorted(status) == ["<hidden>", "a", "b", "only_views[0]", "x", "y"]
    assert status["a"] is root.a_._status and status["<hidden>"] is hidden._status
    assert graph.edges == {
        (status["a"], status["b"]),
        (status["b"], status["<hidden>"]),
    }
    assert graph.pairs == {(status["x"], status["y"])}


def test_around():
    """DependGraph: around keeps the objects computed from, and used to compute,
    the selected ones, within the given number of links"""

    root, hidden = small_graph()
    graph = DependGraph(root)

    def around(*patterns, depth=None):
        return sorted(graph.names[s] for s in graph.around(patterns, depth))

    assert around("b") == ["<hidden>", "a", "b"]
    assert around("a") == ["<hidden>", "a", "b"]
    assert around("a", depth=1) == ["a", "b"]
    assert around("b", depth=0) == ["b"]
    # each object of a synchronized pair is computed from the other
    assert around("y") == ["x", "y"]
    assert around("b", "x", depth=0) == ["b", "x"]
    with pytest.raises(ValueError):
        graph.around(["c"])


def test_dot():
    """DependGraph: the dot file has a node per object, with its kind, an arrow
    per dependency and a dashed line per synchronized pair"""

    root, hidden = small_graph()
    graph = DependGraph(root)
    text = graph.dot()
    lines = text.splitlines()

    assert lines[0] == "digraph depend {" and lines[-1] == "}"
    assert sum(" -> " in line for line in lines) == 3
    assert sum("style=dashed" in line for line in lines) == 1
    assert 'label="b", shape=ellipse' in text
    assert 'label="a", shape=box' in text
    assert 'label="x", shape=box, style="rounded"' in text
    assert 'label="<hidden>"' in text
    assert "penwidth" not in text

    selected = graph.dot(around=["b"], depth=0)
    assert 'label="b", shape=ellipse, penwidth=3' in selected
    assert " -> " not in selected


def test_graph_of_a_simulation():
    """DependGraph: the objects of a bound simulation are named by their path,
    and the positions of each replica, views of those of the beads, are not
    separate nodes"""

    simulation = Simulation.load_from_xml(harmonic_input(nbeads=4), read_only=True)
    graph = DependGraph(simulation)
    status = {name: s for s, name in graph.names.items()}

    q, qnm = status["syslist[0].beads.q"], status["syslist[0].nm.qnm"]
    assert {q, qnm} in [set(pair) for pair in graph.pairs]
    assert status["syslist[0].forces.pot"] in graph.around(["beads.q"])
    assert "syslist[0].beads._blist[0].q" not in status
    assert "syslist[0].beads._blist[0].kin" in status


@pytest.mark.parametrize("output", ["graph.dot", "graph.svg"])
def test_script(tmp_path, monkeypatch, output):
    """i-pi-depend-graph: writes the graph of an input file, without leaving
    other files, and renders it with Graphviz"""

    if output.endswith(".svg") and shutil.which("dot") is None:
        pytest.skip("Graphviz is not installed")
    monkeypatch.chdir(tmp_path)
    (tmp_path / "input.xml").write_text(harmonic_input(nbeads=2))
    argv = ["input.xml", "-o", output, "--around", "beads.q", "--depth", "1"]
    monkeypatch.setattr(sys, "argv", ["i-pi-depend-graph", *argv])
    runpy.run_path(str(SCRIPT), run_name="__main__")

    assert sorted(path.name for path in tmp_path.iterdir()) == [output, "input.xml"]
    text = (tmp_path / output).read_text()
    assert "penwidth=3" in text if output.endswith(".dot") else "<svg" in text
