""" Utility functions to parse i-PI output files.

These are meant to be used in Python post-processing pipelines, so
trajectory files are read as ASE objects (assuming units to be 
Angstrom and eV), and output files are read in a dictionary of 
numpy arrays.
"""

import re
import numpy as np
from ipi.utils.units import unit_to_user
import sys
import os

try:
    import ase
    from ase.io import write
except ImportError:
    ase = None

__all__ = ["read_output", "read_trajectory"]


class SuppressOutput:
    """
    A context manager to suppress print output by redirecting sys.stdout to os.devnull.

    Attributes:
        enable (bool): If True, suppresses output. Defaults to True.
    """

    def __init__(self, enable=True):
        """
        Initializes SuppressOutput with the option to enable or disable suppression.

        Args:
            enable (bool): Whether to suppress output. Defaults to True.
        """
        self.enable = enable

    def __enter__(self):
        """
        Redirects sys.stdout to os.devnull if suppression is enabled.
        """
        if self.enable:
            self._original_stdout = sys.stdout
            sys.stdout = open(os.devnull, "w")
        return self

    def __exit__(self, *argv, **kwargs):
        """
        Restores sys.stdout to its original state after exiting the context.
        """
        if self.enable:
            sys.stdout.close()
            sys.stdout = self._original_stdout


class AppendableList:
    """
    A class that can append numpy arrays to a pre-allocated array.
    If the array is full, it doubles its size.
    """

    def __init__(self, size: int = 100000) -> None:
        """
        Initialize a new AppendableArray.

        Args:
            size (int): The initial size of the array. Defaults to 100000.
        """
        self._arr = [None] * size
        self._size = 0
        self._max_size = size
        self._n_update = 0

    def append(self, x) -> None:
        """
        Append a variable to the list.

        Args:
            x (object): The variable to append.

        """
        # If the array is full, double its size
        if self._size + 1 > self._max_size:
            self._expand()
        # Append the float to the array
        self._arr[self._size] = x
        # Increment the size counter
        self._size += 1

    def _expand(self) -> None:
        """
        Double the size of the array.
        """
        new_size = self._max_size * 2
        new_arr = [None] * new_size
        new_arr[: self._size] = self.finalize()
        self._arr = new_arr
        self._max_size = new_size
        self._n_update += 1

    def finalize(self):
        """
        Return the final array.

        Returns:
            T: The final array.
        """
        return self._arr[: self._size]


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
    bohr2angstrom = unit_to_user("length", "angstrom", 1.0)
    comment_regex = re.compile(r"([^)]+)\{([^}]+)\}")
    step_regex = re.compile(r"Step:\s+(\d+)")

    # Do not use `list.append`
    # frames = []

    # Let's optimize this thing
    frames = AppendableList()

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

    return frames.finalize()


def merge_files(prefix, folder, bead):
    """
    Reads a set of i-PI output files and merges all the trajectory information
    from files with the same bead number. The function will raise an error if
    not all the beads have the same number of arrays or if no 'positions' array
    is found.

    Parameters
    ----------
    prefix : str
        The prefix of the i-PI output files.
    folder : str
        The folder where the i-PI output files are stored.
    bead : int
        The bead number to be merged.

    Returns
    -------
    traj : List[Atoms]
        The merged trajectory information.
    """

    # author: Elia Stocco
    # comments:
    # The script cycle over the arrays and then over then over the atoms.
    # This is slower than the other way around but it is more readable
    # and it is less memory expensive, especially for large files.

    if ase is None:
        raise ImportError(
            "read_trajectory requires the `ase` package to return the structure in ASE format"
        )
    # pattern = f"{folder}/{prefix}.*_{bead}.xyz"
    # files = glob.glob(pattern)

    # Construct the regex pattern
    regex_pattern = re.compile(
        rf"^{re.escape(prefix)}\.\w+_{re.escape(str(bead))}\.xyz$"
    )
    # Filter files using regex
    files = [
        os.path.join(folder, filename)
        for filename in os.listdir(folder)
        if regex_pattern.match(filename)
    ]

    N = len(files)
    if N == 0:
        raise ValueError("No files found.")

    # ------------------#
    pattern_extract = r"([^/]+)\.(.*?)_(.*?)\.(\w+)$"
    names = set()
    for n, file in enumerate(files):
        matched = re.search(pattern_extract, os.path.basename(file))
        name = matched.group(2)
        names.add(name)

    if "positions" not in names:
        raise ValueError("No 'positions' array found.")

    # ------------------#
    files = sorted(files, key=lambda x: "positions" not in x)
    traj = None
    for n, file in enumerate(files):  # cycle over arrays
        matched = re.search(pattern_extract, os.path.basename(file))
        name = matched.group(2)
        if n == 0:
            traj = read_trajectory(file)
            if "positions" not in file:
                raise ValueError(
                    f"File {file} is not a position file but is the first file."
                )
        else:
            if traj is None:
                raise ValueError(
                    f"File {file} is not a position file but is not the first file."
                )
            array = read_trajectory(file)
            for i in range(len(traj)):  # cycle over atoms
                traj[i].arrays[name] = array[i].positions

    return traj


def merge_beads(prefix, folder, nbeads, ofile):
    """
    Convert i-PI output files to ASE trajectory files.

    Parameters
    ----------
    prefix : str
        The prefix of the input files.
    folder : str
        The folder where the input files are located.
    nbeads : int
        The number of beads in the i-PI simulation.
    ofile : str
        The output ASE trajectory file.

    Notes
    -----
    The input files are expected to be named as follows:
    `<folder>/<prefix>.*_<bead>.xyz`, where `<prefix>` is the given prefix,
    `<bead>` is the bead number (starting from 0), and `*` is a wildcard
    matching any characters. The output ASE trajectory file will be named
    `<ofile>.<bead>.extxyz`, where `<bead>` is the bead number.
    """
    if ase is None:
        raise ImportError(
            "read_trajectory requires the `ase` package to return the structure in ASE format"
        )
    for n in range(nbeads):  # cycle over arrays
        traj = merge_files(prefix, folder, n)
        base_name, ext = os.path.splitext(ofile)
        assert (
            ext == ".extxyz"
        ), f"Output file must have extension '.extxyz' but has extension '{ext}'"
        tmp_file = f"{base_name}.{n}{ext}"
        write(tmp_file, traj)


def merge_trajectories(files, names, strides, formats):
    """
    Merge multiple i-PI output files into a single ASE trajectory object.

    Parameters
    ----------
    files : list[str]
        A list of file paths to be merged. The first file must contain the positions.
    names : list[str]
        A list of array names corresponding to the data arrays in the input files.
        The first name must be `'positions'`.

    strides : list[int]
        A list of integers specifying the stride (sampling rate) to apply to each
        trajectory file. If `stride[i] = 1`, all frames from file `i` are included.
        A larger stride will select every `stride[i]`-th frame.

    Returns
    -------
    traj : ASE Atoms object
        The merged trajectory, represented as an ASE Atoms object containing the
        atomic positions and any other arrays (e.g., forces) from the input files.

    Notes
    -----
    - The input files are expected to have the format `<prefix>.*_0.xyz`, where:
        - `<prefix>` is the given prefix for the files,
        - `*` is a wildcard matching any characters,
        - `<bead>` is the bead number (starting from 0).
    - The output ASE trajectory will be stored in the file `<prefix>.extxyz`, where
      `<prefix>` is the given prefix.

    - The function assumes that the first file corresponds to atomic positions and
      subsequent files correspond to other data arrays, such as forces or velocities.
      All files are expected to follow this convention.

    - The function is not highly optimized, as it reads the entire trajectory for each
      file and then slices it based on the provided stride. In the future, it would be
      more efficient to directly pass the stride to the trajectory reader.

    """
    # author: Elia Stocco
    # comments:
    #   The script loops over the arrays and then over then over the atoms.
    #   This is slower than the other way around but it is more readable
    #   and it is less memory expensive, especially for large files.

    # ------------------#
    if ase is None:
        raise ImportError(
            "merge_trajectories requires the `ase` package to return the structure in ASE format"
        )

    # ------------------#
    if "positions" not in names:
        raise ValueError("No 'positions' array found.")
    if names[0] != "positions":
        raise ValueError("The first elements must be `positions`.")

    # ------------------#
    pattern_extract = r"([^/]+)\.(.*?)_(.*?)\.(\w+)$"
    # check_names = set()
    for n, file in enumerate(files):
        matched = re.search(pattern_extract, os.path.basename(file))
        name = matched.group(2)
        assert name == names[n], "some error"

    # ------------------#
    files = sorted(files, key=lambda x: "positions" not in x)
    traj = None
    for n, file in enumerate(files):  # cycle over arrays
        matched = re.search(pattern_extract, os.path.basename(file))
        name = matched.group(2)
        print(f"\t - reading file '{file}'")
        if n == 0:

            # Read the trajectory to get the positions
            with SuppressOutput():
                traj = read_trajectory(filename=file, format=formats[n])

            # ATTENTION:
            # This is highly inefficient because the codes reads the whole trajectory
            # and then takes only the snapshots of interest.
            # In the future one should provide `stride` directly to `read_trajectory`.

            # slice the trajectory using the stride
            stride = strides[n]
            if stride != 1:
                traj = traj[::stride]

            if "positions" not in file:
                raise ValueError(
                    f"File {file} is not a position file but is the first file."
                )
        else:
            if traj is None:
                raise ValueError(
                    f"File {file} is not a position file but is not the first file."
                )

            # Read the trajectory to get other arrays
            with SuppressOutput():
                array = read_trajectory(file)

            # ATTENTION:
            # This is highly inefficient because the codes reads the whole trajectory
            # and then takes only the snapshots of interest.
            # In the future one should provide `stride` directly to `read_trajectory`.

            # slice the trajectory using the stride
            stride = strides[n]
            if stride != 1:
                array = array[::stride]

            # add the read trajectory to the output as `ase.Atoms.arrays`
            for i in range(len(traj)):  # cycle over atoms
                traj[i].arrays[name] = array[i].positions

    return traj


def create_classical_trajectory(input_file, trajectories, properties):

    # author: Elia Stocco
    # date: November 12th, 2024
    # comments:
    #   the code could be vastly improved by avoiding reading the whole trajectory from file
    #   Have a look at the `ATTENTION` comments.

    if len(trajectories) > 0:
        # add more keys if you need
        available = [
            "positions",
            "forces",
            "velocities",
            "momenta",
            "becx",
            "becy",
            "becz",
            "Eforces",
        ]
        okay = [a in available for a in trajectories]
        assert all(okay), "Some provided trajectories are not available."

    if "positions" not in trajectories:
        trajectories.append("positions")
    trajectories = sorted(trajectories, key=lambda x: "positions" not in x)

    # import these classes and modules here to avoid circular import errors.
    from ipi.utils.io.inputs import io_xml
    from ipi.inputs.simulation import InputSimulation
    from ipi.engine.outputs import TrajectoryOutput, PropertyOutput

    # suppress the output fo messages to screen
    with SuppressOutput():
        xml = io_xml.xml_parse_file(input_file)
        # Initializes the simulation class.
        isimul = InputSimulation()
        isimul.parse(xml.fields[0][1])
        simul = isimul.fetch()

    ### Beads
    # check the number of beads
    nbeads = simul.syslist[0].beads.nbeads
    assert nbeads == 1, "A classical trajectory can only have `nbeads` == 1."

    ### Prefix
    prefix = simul.outtemplate.prefix

    ### Properties
    # It's gonna be quick to check if everything is alright.
    # Let's prepare the necessary variables and used them and the end of this function.

    ipi_props = [a for a in simul.outtemplate if isinstance(a, PropertyOutput)]
    assert len(ipi_props) <= 1, "Only one (or zero) property output is supported"
    ipi_props = ipi_props[0] if len(ipi_props) == 1 else None

    if ipi_props is not None:

        # read the properties
        prop_file = f"{prefix}.{ipi_props.filename}"
        if len(properties) > 0 and not os.path.exists(prop_file):
            raise ValueError(f"File {prop_file} does not exists.")
        props, _ = read_output(prop_file)

        # keep only the properties of interest
        props = {key: props[key] for key in properties if key in props}

    ### Trajectories files
    # Time and memory consuming.

    # check that positions have been saved to file
    ipi_trajs = [a for a in simul.outtemplate if isinstance(a, TrajectoryOutput)]
    whats = [a.what for a in ipi_trajs]
    assert "positions" in whats, "Positions have not beed saved to file."

    # keep only the trajectories of interest
    ipi_trajs = [a for a in ipi_trajs if a.what in trajectories]

    # move `positions` to the first place
    ipi_trajs = sorted(ipi_trajs, key=lambda x: "positions" not in x.what)
    # overwrite `trajectories` so that it will have the same ordering of `ipi_trajs`
    trajectories = [a.what for a in ipi_trajs]

    ### Stride
    # check the stride (of both properties and trajectories)
    ipi_trajs_props = ipi_trajs if ipi_props is None else ipi_trajs + [ipi_props]
    strides = [a.stride for a in ipi_trajs_props]
    max_stride = max(strides)
    okay = [max_stride % a == 0 for a in strides]
    assert all(okay), "Strides are not compatibles."
    strides = [int(max_stride / s) for s in strides]
    assert all(
        [a * b.stride == max_stride for a, b in zip(strides, ipi_trajs)]
    ), "Some errors occurred."
    traj_strides = strides[:-1]

    # check that all the trajectory files exist
    traj_files = [f"{prefix}.{traj.filename}_0.xyz" for traj in ipi_trajs]
    for file in traj_files:
        assert os.path.exists(file), f"File '{file}' does not exist."

    ### Trajectory building
    # build the output trajectory

    formats = [a.format for a in ipi_trajs]

    # build the trajectory
    traj = merge_trajectories(traj_files, trajectories, traj_strides, formats)

    # check that all the trajectories of interest have been correctly read
    keys = traj[0].arrays.keys()
    for name in trajectories:
        assert name in keys, f"'{name}' is not in the output trajectory."

    ### Properties adding

    # ATTENTION:
    # This is highly inefficient because the codes reads all the properties
    # and then takes only the snapshots of interest.
    # In the future one should provide `stride` directly to `read_output`.

    strided_properties = {}
    stride = strides[-1]
    for p in props.keys():
        strided_properties[p] = props[p][::stride]

    # Let's loop over the atoms and then over the properties
    # Maybe the other way around might be faster.
    for n, snapshot in enumerate(traj):
        for p in props.keys():
            snapshot.info[p] = strided_properties[p][n]

    return traj


def create_centroid_trajectory(input_file, trajectories):

    # author: Elia Stocco
    # date: November 12th, 2024
    # comments:
    #   the code could be vastly improved by avoiding reading the whole trajectory from file
    #   Have a look at the `ATTENTION` comments.

    if len(trajectories) > 0:
        # add more keys if you need
        available = ["x_centroid", "v_centroid", "f_centroid", "p_centroid"]
        okay = [a in available for a in trajectories]
        assert all(okay), "Some provided trajectories are not available."

    if "x_centroid" not in trajectories:
        trajectories.append("x_centroid")
    trajectories = sorted(trajectories, key=lambda x: "x_centroid" not in x)

    # import these classes and modules here to avoid circular import errors.
    from ipi.utils.io.inputs import io_xml
    from ipi.inputs.simulation import InputSimulation
    from ipi.engine.outputs import TrajectoryOutput

    # suppress the output fo messages to screen
    with SuppressOutput():
        xml = io_xml.xml_parse_file(input_file)
        # Initializes the simulation class.
        isimul = InputSimulation()
        isimul.parse(xml.fields[0][1])
        simul = isimul.fetch()

    ### Beads
    # check the number of beads
    nbeads = simul.syslist[0].beads.nbeads
    assert nbeads == 1, "A classical trajectory can only have `nbeads` == 1."

    ### Prefix
    prefix = simul.outtemplate.prefix

    ### Trajectories files
    # Time and memory consuming.

    # check that positions have been saved to file
    ipi_trajs = [a for a in simul.outtemplate if isinstance(a, TrajectoryOutput)]
    whats = [a.what for a in ipi_trajs]
    assert "x_centroid" in whats, "`x_centroid` have not beed saved to file."

    # keep only the trajectories of interest
    ipi_trajs = [a for a in ipi_trajs if a.what in trajectories]

    # move `x_centroid` to the first place
    ipi_trajs = sorted(ipi_trajs, key=lambda x: "x_centroid" not in x.what)
    # overwrite `trajectories` so that it will have the same ordering of `ipi_trajs`
    trajectories = [a.what for a in ipi_trajs]

    ### Stride
    # check the stride (of both properties and trajectories)
    strides = [a.stride for a in ipi_trajs]
    max_stride = max(strides)
    okay = [max_stride % a == 0 for a in strides]
    assert all(okay), "Strides are not compatibles."
    strides = [int(max_stride / s) for s in strides]
    assert all(
        [a * b.stride == max_stride for a, b in zip(strides, ipi_trajs)]
    ), "Some errors occurred."

    # check that all the trajectory files exist
    traj_files = [f"{prefix}.{traj.filename}.xyz" for traj in ipi_trajs]
    for file in traj_files:
        assert os.path.exists(file), f"File '{file}' does not exist."

    ### Trajectory building
    # build the output trajectory

    formats = [a.format for a in ipi_trajs]

    # build the trajectory
    traj = merge_trajectories(traj_files, trajectories, strides, formats)

    # check that all the trajectories of interest have been correctly read
    keys = traj[0].arrays.keys()
    for name in trajectories:
        assert name in keys, f"'{name}' is not in the output trajectory."

    return traj
