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
        description="Upgrade input file from v2.0 to v3.0."
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
