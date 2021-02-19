import subprocess as sp
import os
from pathlib import Path
from distutils.dir_util import copy_tree
import xml.etree.ElementTree as ET
import tempfile
import time

fortran_driver_models = [
    "dummy",
    "lj",
    "sg",
    "harm",
    "harm3d",
    "morse",
    "zundel",
    "qtip4pf",
    "pswater",
    "eckart",
    "ch4hcbe",
    "MB",
]

# YL should do this  automatically but fFor now I do it explicitly
python_driver_models = ["dummy", "harmonic"]


def get_test_settings(
    example_folder,
    settings_file="test_settings.dat",
    driver_code="fortran",
    driver_model="dummy",
    socket_mode="unix",
    port_number=33333,
    address_name="localhost",
    flags=[],
    nsteps="2",
):
    """This function looks for the existence of test_settings.txt file.
    This file can contain instructions like number of steps or driver name.
    If the file doesn't  exist, the driver dummy is assigned."""
    try:
        with open(Path(example_folder) / settings_file) as f:
            flags = list()
            while True:
                line = f.readline()
                if not line:
                    break
                elif "driver_model" in line:
                    driver_model = line.split()[1]
                elif "socket_mode" in line:
                    socket_mode = line.split()[1]
                elif "port" in line:
                    port_number = line.split()[1]
                elif "address" in line:
                    address_name = line.split()[1]
                elif "flags" in line:
                    flags.append({line.split()[1]: line.split()[2:]})
                elif "nsteps" in line:
                    nsteps = line.split()[1]
                elif "driver_code" in line:
                    driver_code = line.split()[1]
    except:
        pass

    if driver_code == "fortran":
        if driver_model not in fortran_driver_models:
            driver = "dummy"
    elif driver_code == "python":
        if driver_model not in python_driver_models:
            driver = "dummy"
    else:
        raise ValueError(
            "Drive code not available. Valid options are 'fortran' and 'python'"
        )

    driver_info = {
        "driver_model": driver_model,
        "socket_mode": socket_mode,
        "port_number": port_number,
        "address_name": address_name,
        "driver_code": driver_code,
        "flag": flags,
    }

    test_settings = {"nsteps": nsteps}

    return driver_info, test_settings
