#!/usr/bin/env python3
import argparse
from xml.etree import ElementTree as ET
import xml.dom.minidom


def prettify(elem):
    """
    https://stackoverflow.com/questions/17402323/use-xml-etree-elementtree-to-print-nicely-formatted-xml-files
    https://stackoverflow.com/questions/14479656/empty-lines-while-using-minidom-toprettyxml
    """
    rough_string = ET.tostring(elem)
    reparsed = xml.dom.minidom.parseString(rough_string)
    return "\n".join(
        [
            line
            for line in reparsed.toprettyxml(indent=" " * 2).split("\n")
            if line.strip()
        ]
    )


def set_simulation_version_if_not_present(input_path):
    tree = ET.parse(input_path)

    simulation_element = tree.getroot()
    if "version" in simulation_element.attrib:
        raise Exception("Input file already has a version, leaving unchanged")
    simulation_element.set("version", "3.1")

    with open(input_path, "wt") as f:
        f.write(prettify(tree.getroot()))


def main():
    parser = argparse.ArgumentParser(
        description="Adding version field to input file, with version 3.1"
    )
    parser.add_argument(
        "-ifile",
        "--input_file",
        required=True,
        type=str,
        default=None,
        help="the relative path to the input file, edited in place",
    )

    args = parser.parse_args()
    input_path = args.input_file

    set_simulation_version_if_not_present(input_path)


if __name__ == "__main__":
    main()
