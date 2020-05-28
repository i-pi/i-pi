import xml.etree.ElementTree as ET
import sys


def modify_xml(
    input_name,
    output_name,
    socket_mode="unix",
    port_number=33333,
    address_name="localhost",
    nsteps=5,
):
    """ Modify xml to run dummy tests """
    tree = ET.parse(input_name)
    root = tree.getroot()
    clients = list()

    for ffsocket in root.findall("ffsocket"):
        name = ffsocket.attrib["name"]
        ffsocket.attrib["mode"] = socket_mode

        for element in ffsocket:
            port = port_number + len(clients)
            if element.tag == "port":
                element.text = str(port)
            elif element.tag == "address":
                element.text = address_name + "_" + str(port)
                address = address_name + "_" + str(port)
        clients.append((port, address))

    element = root.find("total_steps")
    if element is not None:
        element.text = str(nsteps)
    else:
        new_element = ET.SubElement(root, "total_steps")
        new_element.text = str(nsteps)

    tree.write(open(output_name, "wb"))


if __name__ == "__main__":
    name = sys.argv[1]
    out = "output.xml"
    modify_xml(name, out)
