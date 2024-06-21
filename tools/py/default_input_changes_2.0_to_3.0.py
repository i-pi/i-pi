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


def set_all_ffsocket_pbc_if_default(input_path):
    """
    In version 2.0, the default value of the pbc flag in ffsocket was True.
    In version 3.0, the default value of the pbc flag in ffsocket was False.
    We recommend using the new default in mose cases, because it should not
    change your simulations. A few clients only wored with False and had that
    flag already explicit. However, if you, for some (unknown :-)) reason really
    need to reproduce the behavior of version 2.0, this will write "pbc=True" in
    the right place. Note that this script does not check whether your input was
    for v 2.0 or 3.0 - as we have no way to check at the moment. It is your
    responsibility to know.
    """
    tree = ET.parse(input_path)

    # making the default value of pbc in v2.0 explicit in the input (changed in v3.0)
    ffsocket_elements = tree.getroot().iter("ffsocket")
    for ffsocket_element in ffsocket_elements:
        if "pbc" not in ffsocket_element.attrib:
            ffsocket_element.set("pbc", "True")

    with open(input_path, "wt") as f:
        f.write(prettify(tree.getroot()))


def main():
    parser = argparse.ArgumentParser(
        description="Compensating for input default changes from v2.0 to v3.0."
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

    set_all_ffsocket_pbc_if_default(input_path)


if __name__ == "__main__":
    main()
