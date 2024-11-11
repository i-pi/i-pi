""" Utility functions to parse i-PI output files.

These are meant to be used in Python post-processing pipelines, so
trajectory files are read as ASE objects (assuming units to be 
Angstrom and eV), and output files are read in a dictionary of 
numpy arrays.
"""

import re
import numpy as np
from ipi.utils.units import unit_to_user

import glob
import os
from warnings import warn
from typing import List
from ase import Atoms

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
            (
                col_type,
                start_col,
                end_col,
                property_name,
                units,
                description,
            ) = match.groups()
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
    bohr2angstrom = unit_to_user("length", "angstrom", 1.0)
    comment_regex = re.compile(r"(\w+)\{([^}]+)\}")
    step_regex = re.compile(r"Step:\s+(\d+)")

    frames = []
    while True:
        try:
            if format == "extras":
                file = open(filename, "r")
                step_regex = re.compile(r"#EXTRAS\((\w+)\)# *Step:\s+(\d+)")
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
                    frame.arrays["forces"] = ret["atoms"].q.reshape(
                        (-1, 3)
                    ) * unit_to_user("force", "ev/ang", 1.0)

                frames.append(frame)

        except EOFError:
            break
        except:
            raise

    return frames

def merge_beads(prefix:str,folder:str,bead:int):
    
    #------------------#  
    # if bead is not None:
    #     pattern = f"{folder}/{prefix}.*_*.xyz"
    # else:
    pattern = f"{folder}/{prefix}.*_{bead}.xyz"
    print(f"\n\t Looking for files mathing the following pattern: '{pattern}'")
    files = glob.glob(pattern)
    N = len(files)
    print(f"\t {N} files found.")
    if N == 0:
        raise ValueError("No files found.")
    
    #------------------#
    pattern_extract = r'([^/]+)\.(.*?)_(.*?)\.(\w+)$'
    beads = set()
    names = set()
    for n,file in enumerate(files):
        matched = re.search(pattern_extract, os.path.basename(file))
        name = matched.group(2)    # Captures the first `*`
        bead = matched.group(3)   # Captures the second `*`
        names.add(name)
        beads.add(bead)  
    Nn = len(names)        
    Nb = len(beads)  
    print(f"\t {Nb} beads found.")    
    print(f"\t {Nn} arrays found.")    
    assert N == Nb*Nn, f"Found {N} files but {Nb} beads and {Nn} arrays."
    
    if "positions" not in names:
        raise ValueError("No 'positions' array found.")
        
    #------------------#
    print("\n\t Found the following files:")
    # index = integer_to_slice_string(args.index)
    instructions = {}
    files = sorted(files, key=lambda x: 'positions' not in x)
    traj = None
    for n,file in enumerate(files):
        matched = re.search(pattern_extract, os.path.basename(file))
        if n == 0:
            traj:List[Atoms] = read_trajectory(file)
            if "positions" not in file:
                raise ValueError(f"File {file} is not a position file but is the first file.")
        else:
            if traj is None:
                raise ValueError(f"File {file} is not a position file but is not the first file.")
            array:List[Atoms] = read_trajectory(file)
            for i in range(len(traj)):
                traj[i].arrays[name] = array[i].positions
        
        
    #     if matched:
    #         prefix_value = matched.group(1)   # Captures the prefix (variable part)
    #         name = matched.group(2)    # Captures the first `*`
    #         bead = matched.group(3)   # Captures the second `*`
    #         extension = matched.group(4)      # Captures the file extension
    #         assert prefix_value == prefix, f"File {file} does not match the expected pattern"
    #         if extension != "xyz":
    #             warn(f"File {file} was expected to have extension 'xyz' but has extension '{extension}'")
    #         print(f"\t - bead={bead}, name={name}: {file}")
    #         if name not in instructions:
    #             instructions[name] = {}
    #         instructions[name][bead] = Instruction(name=name, bead=bead, filename=file, format="i-pi", index=index, fixed_cell=args.fixed_cell,pbc=args.pbc)
    #     else:
    #         raise ValueError(f"File {file} does not match the expected pattern")
        
    # # #------------------#
    # # print("\n\t Creating ASE trajectories (one for each bead):")
    # # trajs:List[List[Atoms]] = [None]*Nb
    # # for n,bead in enumerate(beads):
    # #     print(f"\t - bead={bead} ... ",end="")
    # #     instr:Instruction = instructions["positions"][bead]
    # #     trajs[n] = instr.get()
    # #     print("done")
    
    # print("\n\t Creating ASE trajectories (one for each bead):")
    # for n,bead in enumerate(beads):
    #     print(f"\t - bead={bead}:")
    #     instr:Instruction = instructions["positions"][bead]
    #     traj = instr.get()
    #     for name in names:
    #         print(f"\t\t - {name} ... ",end="")
    #         if name == "positions":
    #             continue
    #         instr = instructions[name][bead]
    #         arr = instr.get()
    #         for i in range(len(traj)):
    #             traj[i].arrays[name] = arr[i].positions
    #         print("done")
    #     pass