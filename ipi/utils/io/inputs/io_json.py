"""Functions used to read input files in JSON format.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import json
from ipi.utils.io.inputs.io_xml import xml_node

__all__ = [
    "json_parse_file",
    "json_parse_string",
]


def json_parse_file(stream):
    """Parses a JSON input file.

    Args:
        stream: A file object.

    Returns:
        A xml_node for the root node of the file.
    """

    try:
        data = json.load(stream)
    except json.JSONDecodeError as e:
        raise ValueError(f"Error parsing JSON: {e}")

    return json_to_xmlnode(data)


def json_parse_string(string):
    """Parses a JSON input string.

    Args:
        string: A string in JSON format.

    Returns:
        A xml_node for the root node of the file.
    """

    try:
        data = json.loads(string)
    except json.JSONDecodeError as e:
        raise ValueError(f"Error parsing JSON string: {e}")

    return json_to_xmlnode(data)


def json_to_xmlnode(data, name="root"):
    """
    Recursively converts a dictionary (from JSON) into an xml_node structure.
    """

    # Attributes and fields containers
    attribs = {}
    fields = []

    # If data is a primitive, it's just text content
    if not isinstance(data, (dict, list)):
        return xml_node(name=name, fields=[("_text", str(data))])

    # If data is a list, this function should probably not be called directly on it
    # unless it's the root? Root should be a dict.
    if isinstance(data, list):
        # Should not happen for a single node context, unless the JSON structure is
        # "tag": [ ... ]
        # In that case, the parent iterates the list and calls json_to_xmlnode for each item.
        # If it happens, we treat it as text content?
        return xml_node(name=name, fields=[("_text", str(data))])

    # Process attributes
    if "attributes" in data:
        for k, v in data["attributes"].items():
            attribs[k] = str(v)

    # Process value (text content)
    if "value" in data:
        val = data["value"]
        # Convert value to string representation expected by i-PI parsers
        if isinstance(val, list):
            # For arrays/lists, we want "[v1, v2, ...]"
            # json.dumps would produce "[v1, v2, ...]" which is compatible
            # but str(val) is also "[v1, v2, ...]" in Python
            fields.append(("_text", str(val)))
        else:
            fields.append(("_text", str(val)))

    # Process children
    for k, v in data.items():
        if k in ("attributes", "value"):
            continue

        if isinstance(v, list):
            for item in v:
                fields.append((k, json_to_xmlnode(item, name=k)))
        else:
            fields.append((k, json_to_xmlnode(v, name=k)))

    return xml_node(name=name, attribs=attribs, fields=fields)
