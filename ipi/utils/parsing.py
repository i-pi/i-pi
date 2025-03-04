"""Utility functions to parse i-PI output files.

These are meant to be used in Python post-processing pipelines, so
trajectory files are read as ASE objects (assuming units to be
Angstrom and eV), and output files are read in a dictionary of
numpy arrays.
"""

import re
import numpy as np
from ipi.utils.units import unit_to_user
from ipi.utils.messages import warning, verbosity

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
        r"#\s*(column|cols\.)\s+(\d+)(?:-(\d+))?\s*-->\s*([^\s\{\(]+)(?:\{([^\}]+)\})?(?:\(([^\)]+)\))?(?:\{([^\}]+)\})?\s*:\s*(.*)"
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
            (
                col_type,
                start_col,
                end_col,
                property_name,
                units_before,
                args,
                units_after,
                description,
            ) = match.groups()
            col_info = f"{start_col}-{end_col}" if end_col else start_col
            units = units_before
            if units_after is not None:
                units = units_after
            if args is not None:
                property_name += f"({args})"

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
    be inferred from the filename extension. `extras` will read a trajectory containing the additional
    data returned from a calculator; will try to read as a float array, and if it fails it'll resort
    to returning a list of raw strings.
    units can be inferred from the file content, or specified with `dimension`, `units` and `cell_units`
    """

    if ase is None:
        raise ImportError(
            "read_trajectory requires the `ase` package to return the structure in ASE format"
        )

    from ipi.utils.io import read_file

    if format is None:
        # tries to infer the file format
        format = filename[filename.rfind(".") + 1 :]
        if format == "extxyz":
            format = "ase"
        if format not in ["xyz", "pdb", "binary", "json", "ase"]:
            raise ValueError(f"Unrecognized file format: {format}")

    file_handle = open(filename, "r")
    comment_regex = re.compile(r"(\w+)\{([^}]+)\}")
    step_regex = re.compile(r"Step:\s+(\d+)")

    frames = []
    while True:
        try:
            if format == "extras":
                file = open(filename, "r")
                step_regex = re.compile(r"#EXTRAS[^(]*\(([^)]+)\)# *Step:\s+(\d+)")
                step_list = []
                data_list = []
                data_frame = []
                is_array = True
                property_name = ""
                for line in file:
                    matches = step_regex.findall(line)
                    if len(matches) > 0:
                        if len(data_frame) > 0:
                            try:
                                data_processed = np.loadtxt(data_frame)
                            except:
                                is_array = False
                                data_processed = "\n".join(data_frame)
                            data_list.append(data_processed)
                            data_frame = []

                        # Found a new step, update current_step and initialize the list for data
                        step_list.append(int(matches[0][1]))
                        if property_name == "":
                            property_name = matches[0][0]
                        elif property_name != matches[0][0]:
                            raise ValueError(
                                f"Inconsistent property {matches[0][0]} found in extras containing {property_name}"
                            )
                    else:
                        data_frame.append(line)

                if len(data_frame) > 0:
                    try:
                        data_processed = np.loadtxt(data_frame)
                    except:
                        is_array = False
                        data_processed = "\n".join(data_frame)
                    data_list.append(data_processed)
                if is_array:
                    data_list = np.array(data_list).squeeze()
                return {
                    "step": np.asarray(step_list, dtype=int),
                    property_name: data_list,
                }
            else:
                ret = read_file(
                    format,
                    file_handle,
                    dimension=dimension,
                    units=units,
                    cell_units=cell_units,
                )

                frame = ase.Atoms(
                    ret["atoms"].names,
                    # will apply conversion later!
                    positions=ret["atoms"].q.reshape((-1, 3)),
                    cell=ret["cell"].h.T * unit_to_user("length", "ase"),
                    pbc=True,
                )

                # parse comment to get the property
                matches = comment_regex.findall(ret["comment"])
                # get what we have found
                if len(matches) == 2:
                    what = matches[0][0]
                else:  # defaults to reading positions
                    what = "positions"

                # ... and the step
                matches = step_regex.findall(ret["comment"])
                if len(matches) >= 1:
                    frame.info["step"] = int(matches[-1][0])

                # fetch the list of known traj types, cf. `engine/properties.py``
                from ipi.engine.properties import Trajectories  # avoids circular import

                traj_types = Trajectories().traj_dict
                if not what in traj_types:
                    warning(
                        f"{what} is not a known trajectory type. Will apply no units conversion",
                        verbosity.low,
                    )
                elif traj_types[what]["dimension"] == "length":
                    # positions is the right place to store, and we just need to convert
                    frame.positions *= unit_to_user("length", "ase")
                else:
                    # if we have another type of value, set positions to zero
                    # (that data is missing!) and set an array instead
                    frame.positions *= 0.0
                    frame.arrays[what] = ret["atoms"].q.reshape((-1, 3)) * unit_to_user(
                        traj_types[what]["dimension"], "ase"
                    )

                frames.append(frame)

        except EOFError:
            break
        except:
            raise

    return frames
