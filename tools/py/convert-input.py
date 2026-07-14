#!/usr/bin/env python3
"""
Converts between XML and JSON input formats for i-PI.
"""

import sys
import os
import argparse
import json

# Check that we have the import path for this i-PI set and if not, add it.
dir_root = os.path.realpath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "../../")
)
if not dir_root in sys.path:
    sys.path.insert(0, dir_root)

from ipi.utils.io.inputs.io_xml import xml_parse_file, xml_write
from ipi.utils.io.inputs.io_json import (
    json_parse_file,
    xmlnode_to_json,
)


def main():
    parser = argparse.ArgumentParser(
        description="Convert between i-PI XML and JSON input formats."
    )
    parser.add_argument("input_file", help="Input file path")
    parser.add_argument("output_file", help="Output file path")
    parser.add_argument(
        "--indent", type=int, default=2, help="Indentation for JSON output"
    )

    args = parser.parse_args()

    input_ext = os.path.splitext(args.input_file)[1].lower()
    output_ext = os.path.splitext(args.output_file)[1].lower()

    if input_ext == ".xml" and output_ext == ".json":
        # XML to JSON
        try:
            with open(args.input_file, "r") as f:
                root = xml_parse_file(f)

            # root contains "simulation" as a child.
            # xmlnode_to_json(root) will give {"simulation": ...}
            data = xmlnode_to_json(root)

            with open(args.output_file, "w") as f:
                json.dump(data, f, indent=args.indent)
            print(f"Converted {args.input_file} to {args.output_file}")

        except Exception as e:
            print(f"Error converting XML to JSON: {e}")
            sys.exit(1)

    elif input_ext == ".json" and output_ext == ".xml":
        # JSON to XML
        try:
            with open(args.input_file, "r") as f:
                # json_parse_file returns a root node containing "simulation"
                root = json_parse_file(f)

            # xml_write expects the node to write.
            # If we pass 'root', it writes <root>...</root>
            # But we want <simulation>...</simulation>
            # The root returned by json_parse_file has one child: simulation.

            simulation_node = root.fields[0][1]

            with open(args.output_file, "w") as f:
                f.write(xml_write(simulation_node, name=simulation_node.name))
            print(f"Converted {args.input_file} to {args.output_file}")

        except Exception as e:
            print(f"Error converting JSON to XML: {e}")
            sys.exit(1)

    else:
        print("Unsupported conversion. Supported: .xml -> .json, .json -> .xml")
        sys.exit(1)


if __name__ == "__main__":
    main()
