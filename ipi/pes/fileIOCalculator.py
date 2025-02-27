import os
import time
import json
import shutil
import re
import numpy as np
from datetime import datetime
from filelock import FileLock
from ase import Atoms
from ase.io import write
from ase.calculators.calculator import Calculator, all_changes
from typing import Dict, Union, Any
from .ase import ASEDriver


# ---------------------------------------#
def convert_lists_to_arrays(data: Any) -> Any:
    """
    Recursively convert lists in a data structure to NumPy arrays where possible.

    Args:
        data (Any): The input data, typically a dictionary or list.

    Returns:
        Any: The data with lists converted to NumPy arrays where applicable.

    Example:
        >>> data = {"scores": [95, 85, 75], "info": {"ages": [21, 22, 23]}}
        >>> convert_lists_to_arrays(data)
        {'scores': array([95, 85, 75]), 'info': {'ages': array([21, 22, 23])}}
    """
    if isinstance(data, list):
        try:
            return np.array(data)  # Attempt to convert the list to a NumPy array
        except ValueError:
            return data  # If conversion fails, keep it as a list
    elif isinstance(data, dict):
        return {key: convert_lists_to_arrays(value) for key, value in data.items()}
    elif isinstance(data, (tuple, set)):
        return type(data)(convert_lists_to_arrays(value) for value in data)
    return data  # Return other data types unchanged


# ---------------------------------------#
def read_json(file: str) -> dict:
    """
    Load a JSON file into a dictionary, converting lists to NumPy arrays where possible.

    Args:
        file (str): The path to the JSON file.

    Returns:
        dict: The loaded data with lists converted to NumPy arrays.

    Example:
        >>> data = {"scores": [95, 85, 75]}
        >>> with open("data.json", "w") as f:
        ...     json.dump(data, f)
        >>> json2dict("data.json")
        {'scores': array([95, 85, 75])}
    """
    with open(file, "r") as f:
        data = json.load(f)
    return convert_lists_to_arrays(data)


# ---------------------------------------#
class Logger:
    def __init__(
        self, log_file: str = None, debug: bool = False, warning: bool = True
    ) -> None:
        """
        Set up a custom logger for the class.

        Args:
            log_file (str): Path to the log file. If None, logs to standard output (stdout).
        """
        self.log_file = log_file
        self._debug = debug
        self._warning = warning

        flags = extract_flags(log_file)
        if flags["debug"] is not None:
            self._debug = flags["debug"]
        if flags["warning"] is not None:
            self._warning = flags["warning"]

    def _write_log(self, premessage: str, message: str) -> None:
        """
        Writes the log message to the desired output (file or stdout).

        Args:
            message (str): The log message to write.
        """
        # Get the current timestamp
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Format the log message
        log_message = f"{timestamp} - {premessage:<7}: {message:<50}"

        # Write to the log file or stdout
        if self.log_file:
            with open(self.log_file, "a") as f:
                f.write(log_message + "\n")
        else:
            print(log_message)

    def info(self, message: str) -> None:
        """Log an info message."""
        self._write_log("INFO", message)

    def debug(self, message: str) -> None:
        """Log a debug message."""
        if self._debug:
            self._write_log("DEBUG", message)

    def error(self, message: str) -> None:
        """Log an error message."""
        self._write_log("ERROR", message)

    def warning(self, message: str) -> None:
        if self._warning:
            self._write_log("WARNING", message)


def check_exit(logger: Logger = None) -> None:
    """
    Check if an 'EXIT' file exists in the current directory.
    If it exists, terminate the program.
    """
    if os.path.exists("EXIT"):
        if logger is not None:
            logger.info("EXIT file detected. Kill this process.")
        time.sleep(0.001)
        exit(0)


# ---------------------------------------#
def extract_flags(input_string: str) -> dict:
    flags = {"debug": None, "warning": None}
    if input_string is None or len(input_string) == 0:
        return flags
    # Define a regex pattern to find both 'debug' and 'warning' and their values (True or False)
    pattern = r"(debug|warning)=(True|False)"

    # Search for the pattern in the string
    matches = re.findall(pattern, input_string)

    # Populate the dictionary with the values from the matches
    for match in matches:
        flag, value = match
        flags[flag] = value == "True"

    return flags


# ---------------------------------------#
# Mandatory properties to be returned by the calculator
MANDATORY = {
    "energy": (float, 1),
    "free_energy": (float, 1),
    "forces": (float, ("natoms", 3)),
    "stress": (float, (3, 3)),
}


# ---------------------------------------#
def get_impl_prop(file: str) -> Dict[str, Any]:
    if file is None:
        return MANDATORY
    else:
        with open(file, "r") as f:
            data: Dict[str, Any] = json.load(f)
        # Convert lists back to tuples
        impl_prop = {}
        for k, v in data.items():
            if not isinstance(v, list):
                v = [v]
            try:
                v[0] = eval(v[0])
            except:
                pass
            v = tuple(v)
            impl_prop[k] = v
        return impl_prop


# ---------------------------------------#
def prepare_folder(folder: str) -> None:
    """
    Prepare the folder by clearing its contents or creating it.

    Args:
        folder (str): Directory to prepare.
    """
    if os.path.exists(folder):
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)  # Remove file or symbolic link
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)  # Remove directory
    else:
        os.makedirs(folder)


__DRIVER_NAME__ = "FileIODriver"
__DRIVER_CLASS__ = "FileIODriver"


# ---------------------------------------#
class FileIODriver(ASEDriver):

    def __init__(self, template, verbose, *args, **kwargs):
        super().__init__(template, *args, **kwargs)
        self.ase_calculator = FileIOCalculator(verbose=verbose, *args, **kwargs)


class FileIOCalculator(Calculator):
    """
    A custom calculator for handling file-based I/O operations with logging.

    Attributes:
        folder (str): Directory for input/output files.
        logger (Logger): Logger instance for logging operations.
    """

    folder: str
    logger: Logger
    implemented_properties: Dict[str, Any]
    verbose: bool

    def __init__(
        self, folder: str, log_file: str, impl_prop: str, verbose: bool
    ) -> None:
        """
        Initialize the calculator and prepare the folder.

        Args:
            folder (str): Path to the working directory.
            log_file (str): Path to the log file.
        """
        super().__init__()
        self.folder = folder
        self.logger = Logger(log_file)
        prepare_folder(self.folder)
        self.implemented_properties = MANDATORY
        impl_prop = get_impl_prop(impl_prop)
        self.implemented_properties.update(impl_prop)

    def calculate(
        self, atoms: Atoms = None, properties=None, system_changes=all_changes
    ) -> None:
        """
        Perform calculations by handling input/output files and processing results.

        Args:
            atoms (Atoms): Atomic configuration for the calculation.
            properties: Requested properties (not used explicitly).
            system_changes: System changes to account for (default: all changes).
        """
        start_time = time.time()
        self.logger.debug("Starting calculation.")
        super().calculate(atoms, properties, system_changes)

        # Write input file
        ifile = os.path.join(self.folder, "input.extxyz")
        if os.path.exists(ifile):
            self.logger.error(f"Input file already exists: {ifile}")
            raise ValueError(f"Error: {ifile} should not exist.")
        with FileLock(f"{ifile}.lock"):
            write(ifile, atoms, format="extxyz")
        self.logger.debug(f"Input file written: {ifile}")

        # Wait for output file
        ofile = os.path.join(self.folder, "output.json")
        timeout_seconds = None  # Set timeout for output file creation
        start_wait = time.time()
        self.logger.debug(f"Waiting for output file: {ofile}")
        while not os.path.exists(ofile):
            check_exit(self.logger)
            if timeout_seconds and (time.time() - start_wait > timeout_seconds):
                self.logger.error(f"Timeout waiting for output file: {ofile}")
                raise TimeoutError(
                    f"Output file not created within {timeout_seconds} seconds."
                )
            time.sleep(0.001)

        # Validate input file removal
        if os.path.exists(ifile):
            self.logger.error(f"Input file still exists: {ifile}")
            raise ValueError(f"Error: {ifile} should not exist anymore.")

        # Read and validate output file
        with FileLock(f"{ofile}.lock"):
            data: Dict[str, Any] = read_json(ofile)
            self.logger.debug(f"Output file read: {ofile}")
            os.remove(ofile)

        # Process and validate results
        self.results: Dict[str, Union[float, np.ndarray, str]] = {}
        for key in MANDATORY.keys():
            if key not in data:
                self.logger.warning(f"Missing key {key} in output data.")
                if key == "energy":
                    data["energy"] = 0.0
                    self.results["energy"] = 0.0
                elif key == "forces":
                    self.results["forces"] = np.zeros(atoms.positions.shape)
                elif key == "free_energy":
                    self.results["free_energy"] = (
                        data["energy"] if "energy" in data else 0.0
                    )
                elif key == "stress":
                    self.results["stress"] = np.zeros((3, 3))
                else:
                    # raise KeyError(f"Missing key {key} in data. Found keys: {list(data.keys())}.")
                    raise KeyError(f"Coding error for key '{key}'.")
            else:
                self.results[key] = data[key]
                self.logger.debug(f"Processed key {key}")

        # Validate forces and stress shapes
        if self.results["forces"].shape != atoms.positions.shape:
            self.logger.error("Invalid forces shape.")
            raise ValueError("Forces shape mismatch.")
        if self.results["stress"].shape != (3, 3):
            self.logger.error("Invalid stress shape.")
            raise ValueError("Stress shape mismatch.")

        # Log total calculation time
        total_time = time.time() - start_time
        self.logger.debug("Calculation completed successfully.")
        self.logger.info(f"Calculation took {total_time:.2f} seconds.")
