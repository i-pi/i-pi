""" Utility functions to parse i-PI output files.

These are meant to be used in Python post-processing pipelines, so
trajectory files are read as ASE objects (assuming units to be 
Angstrom and eV), and output files are read in a dictionary of 
numpy arrays.
"""

import re
import numpy as np
from ipi.utils.units import unit_to_user

try:
    import ase
except ImportError:
    ase = None

__all__ = ["read_output", "read_trajectory"]


def read_output(filename):
    """Reads an i-PI output file and returns a dictionary with the properties in a tidy order,
    and information on units and descriptions of the content.

    Usage:
        read_output("filename")

    Returns:
        values, info

        values: a dictionary with the property names as keys, and the values as numpy arrays
        info: a dictionary with the property names as keys and as values tuples of (units, description)
    """

    # Regex pattern to match header lines and capture relevant parts
    header_pattern = re.compile(
        r"#\s*(column|cols\.)\s+(\d+)(?:-(\d+))?\s*-->\s*([^\s\{]+)(?:\{([^\}]+)\})?\s*:\s*(.*)"
    )

    # Reading the file
    with open(filename, "r") as file:
        lines = file.readlines()

    header_lines = [line for line in lines if line.startswith("#")]
    data_lines = [line for line in lines if not line.startswith("#") and line.strip()]

    # Interprets properties
    properties = {}
    for line in header_lines:
        match = header_pattern.match(line)
        if match:
            # Extracting matched groups
            col_type, start_col, end_col, property_name, units, description = (
                match.groups()
            )
            col_info = f"{start_col}-{end_col}" if end_col else start_col
            properties[col_info] = {
                "name": property_name,
                "units": units,
                "description": description,
            }

    # Parse data
    values_dict = {}
    info_dict = {}
    for col_info, prop_info in properties.items():
        # Initialize list to hold values for each property
        values_dict[prop_info["name"]] = []
        # Save units and description
        info_dict[prop_info["name"]] = (prop_info["units"], prop_info["description"])

    for line in data_lines:
        values = line.split()
        for column_info, prop_info in properties.items():
            if "-" in column_info:  # Multi-column property
                start_col, end_col = map(
                    int, column_info.split("-")
                )  # 1-based indexing
                prop_values = values[
                    start_col - 1 : end_col
                ]  # Adjust to 0-based indexing
            else:  # Single column property
                col_index = int(column_info) - 1  # Adjust to 0-based indexing
                prop_values = [values[col_index]]

            values_dict[prop_info["name"]].append([float(val) for val in prop_values])

    for prop_name, prop_values in values_dict.items():
        values_dict[prop_name] = np.array(
            prop_values
        ).squeeze()  # make 1-col into a flat array

    return values_dict, info_dict


def read_trajectory(
    filename,
    format=None,
    dimension="automatic",
    units="automatic",
    cell_units="automatic",
):
    """Reads a file in i-PI format and returns it in ASE format.

    `format` can be `xyz` (i-PI formatted), `pdb`, `binary`, `json`, `ase`, and if not specified it'll
    be inferred from the filename extension.
    units can be inferred from the file content, or specified with `units` and `cell_units`
    """

    if ase is None:
        raise ImportError(
            "read_trajectory requires the `ase` package to return the structure in ASE format"
        )

    from ipi.utils.io import read_file

    if format is None:
        # tries to infer the file format
        format = filename[filename.rfind(".") + 1 :]
        if format not in ["xyz", "pdb", "binary", "json", "ase"]:
            raise ValueError(f"Unrecognized file format: {format}")

    file_handle = open(filename, "r")
    bohr2angstrom = unit_to_user("length", "angstrom", 1.0)
    comment_regex = re.compile(r"(\w+)\{([^}]+)\}")
    step_regex = re.compile(r"Step:\s+(\d+)")

    frames = []
    while True:
        try:
            ret = read_file(
                format,
                file_handle,
                dimension=dimension,
                units=units,
                cell_units=cell_units,
            )

            frame = ase.Atoms(
                ret["atoms"].names,
                positions=ret["atoms"].q.reshape((-1, 3)) * bohr2angstrom,
                cell=ret["cell"].h.T * bohr2angstrom,
                pbc=True,
            )

            # parse comment to get the property
            matches = comment_regex.findall(ret["comment"])

            # get what we have found
            if len(matches) >= 2:
                what = matches[-2][0]
            else:  # defaults to reading positions
                what = "positions"

            # ... and the step
            matches = step_regex.findall(ret["comment"])
            if len(matches) >= 1:
                frame.info["step"] = int(matches[-1][0])

            # if we have forces, set positions to zero (that data is missing!) and set forces instead
            if what == "forces":
                # set forces and convert to eV/angstrom
                frame.positions *= 0
                frame.arrays["forces"] = ret["atoms"].q.reshape((-1, 3)) * unit_to_user(
                    "force", "ev/ang", 1.0
                )

            frames.append(frame)

        except EOFError:
            break
        except:
            raise

    return frames
