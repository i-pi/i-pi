#!/usr/bin/env python3

"""Writes the graph of the depend objects of a simulation, i.e. of which
quantities are computed from which (see ipi.utils.tools.depend_graph).

The input file is loaded as i-PI does, without connecting to the clients or
writing outputs, so it should be run where i-PI would be (the input may read
other files). The graph is written in the dot format of Graphviz, or rendered
with the `dot` program if the output file ends with another extension.

Syntax:
   i-pi-depend-graph input.xml [-o graph.svg] [--around beads.q] [--depth 2]
"""

import argparse
import os
import subprocess
import sys

from ipi.engine.simulation import Simulation
from ipi.utils.tools import DependGraph


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Writes the graph of the depend objects of a simulation."
    )
    parser.add_argument("input", help="i-PI input file (XML or JSON)")
    parser.add_argument(
        "-o",
        "--output",
        help="output file: .dot, or another format of Graphviz (.svg, .pdf, "
        ".png) rendered with `dot`. Prints the dot file if not given.",
    )
    parser.add_argument(
        "--around",
        action="append",
        default=[],
        metavar="PATH",
        help="keeps only the objects computed from, or used to compute, the "
        "object at PATH (e.g. beads.q or syslist[0].forces.pot). Can be given "
        "more than once.",
    )
    parser.add_argument(
        "--depth",
        type=int,
        help="with --around, keeps only the objects at most DEPTH links away",
    )
    args = parser.parse_args(argv)

    with open(args.input) as xml:
        simulation = Simulation.load_from_xml(
            xml, custom_verbosity="quiet", read_only=True
        )
    text = DependGraph(simulation).dot(args.around, args.depth)

    if args.output is None:
        sys.stdout.write(text)
        return
    fmt = os.path.splitext(args.output)[1][1:]
    if fmt == "dot":
        with open(args.output, "w") as out:
            out.write(text)
    else:
        subprocess.run(
            ["dot", f"-T{fmt}", "-o", args.output], input=text, text=True, check=True
        )


if __name__ == "__main__":
    main()
