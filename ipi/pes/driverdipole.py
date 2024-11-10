"""Interface with e3nn to run machine learned electric dipole"""

import sys
import json
import numpy as np
import warnings
import importlib
from .dummy import Dummy_driver

__DRIVER_NAME__ = "driverdipole"
__DRIVER_CLASS__ = "driverdipole_driver"

fmt = "%26.20f"  # output format for dipole and BEC

# Some comments:
# The 'driverdipole_driver' allows sending to i-PI the dipole and their derivatives w.r.t. nuclear positions to i-PI
# You can "plug in" (almost) any kind of neural network by simply providing:
# - a JSON file with the following keys (have a look at the 'get_model' function):
#   - "module": the module where the class is provided;
#   - "class": the class name;
#   - "kwargs": the arguments to properly initialize the network;
# - a *.pth file with the parameters of the network.
# - an optional JSON to modify the parameters contained in 'driverdipole_driver.opts_default'
#
# The only requirements that your network has to satisfy are:
# - it should have the 'get' and 'get_value_and_jac' methods defined,
#       which provide the dipole, and the dipole+its jacobian w.r.t. nuclear positions respectively;
# - since a '_symbols' attribute will be added to the network by 'get_model', no conflict should occur due to that;
# - the network should be aware of the atomic species "inside" the methods 'get' and 'get_value_and_jac' through the '_symbols' attribute.


def recursive_copy(source_dict: dict, target_dict: dict) -> dict:
    """Recursively copy keys and values from a source dictionary to a target dictionary, if they are not present in the target."""
    for key, value in source_dict.items():
        if (
            isinstance(value, dict)
            and key in target_dict
            and isinstance(target_dict[key], dict)
        ):
            recursive_copy(value, target_dict[key])
        else:
            if key not in target_dict:
                target_dict[key] = value
    return target_dict


def add_default(dictionary: dict = None, default: dict = None) -> dict:
    """Add default key-value pairs to a dictionary if they are not present."""
    if dictionary is None:
        dictionary = {}
    if not isinstance(dictionary, dict):
        raise ValueError("'dictionary' has to be of 'dict' type")
    return recursive_copy(source_dict=default, target_dict=dictionary)


def get_class(module_name, class_name):
    try:
        # Import the module dynamically
        module = importlib.import_module(module_name)
        # Get the class from the module
        class_obj = getattr(module, class_name)
        # Create an instance of the class
        # instance = class_obj()
        return class_obj
    except ImportError:
        raise ValueError(f"Module '{module_name}' not found.")
    except AttributeError:
        raise ValueError(f"Class '{class_name}' not found in module '{module_name}'.")
    except Exception as e:
        raise ValueError(f"An error occurred: {e}")


def get_model(instructions, parameters: str):
    import torch

    if isinstance(instructions, str):
        with open(instructions, "r") as json_file:
            _instructions = json.load(json_file)
        instructions = _instructions

    # extract values for the instructions
    kwargs = instructions["kwargs"]
    cls = instructions["class"]
    mod = instructions["module"]

    # get the class to be instantiated
    # Call the function and suppress the warning
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")  # , category=UserWarning)
        class_obj = get_class(mod, cls)

    # instantiate class
    model = class_obj(**kwargs)
    if not model:
        raise ValueError(
            "Error instantiating class '{:s}' from module '{:s}'".format(cls, mod)
        )

    try:
        N = model.n_parameters()
        print("\tLoaded model has {:d} parameters".format(N))
    except:
        pass  # print("\tCannot count parameters")

    # Load the parameters from the saved file
    checkpoint = torch.load(parameters)

    # Update the model's state dictionary with the loaded parameters
    model.load_state_dict(checkpoint)
    model.eval()

    # Store the chemical species that will be used during the simulation.
    if "chemical-symbols" not in instructions:
        raise ValueError("'instructions' should contain the key 'chemical-symbols'")
    model._symbols = instructions["chemical-symbols"]

    return model


class driverdipole_driver(Dummy_driver):
    opts_default = {
        "dipole": {
            "send": True,  # whether to send the dipole to i-PI
            "file": None,  # the file where the dipole will be saved (using np.savetxt)
        },
        "BEC": {
            "compute": True,  # whether to compute the BEC
            "send": True,  # whether to send the BEC to i-PI
            "file": None,  # the file where the BEC will be saved (using np.savetxt)
        },
        "restart": False,  # whether remove the files (if already existing) where the dipole and BEC will be saved.
    }

    def __init__(self, args=None):
        self.error_msg = """The parameters of 'driverdipole_driver' are not correctly formatted. \
            They should be two or three strings, separated by a comma."""
        self.opts = dict()
        self.count = 0
        super().__init__(args)

    def check_arguments(self):
        """Check the arguments required to run the driver."""
        try:
            arglist = self.args.split(",")
        except ValueError:
            sys.exit(self.error_msg)

        if len(arglist) >= 2:
            info_file = arglist[0]  # json file to properly allocate a 'model' object
            parameters_file = arglist[1]  # *.pth file with the model parameters
            try:
                opts_file = arglist[2]  # json file with some parameters for this class
            except:
                print("\tNo options file provided: using the default values")
                opts_file = None
        else:
            sys.exit(self.error_msg)  # to be modified

        print("\n\tThe driver is 'driverdipole_driver'")
        print("\tLoading model ...")
        self.model = get_model(info_file, parameters_file)

        if opts_file is not None:
            try:
                # Open and read the JSON file
                with open(opts_file, "r") as json_file:
                    # Parse the JSON data and save it to a variable
                    self.opts = json.load(json_file)
                self.opts = add_default(self.opts, self.opts_default)
            except:
                print("\tError reading the options file '{:s}'".format(opts_file))
                print("\tThe default values will be used")
                self.opts = add_default(None, self.opts_default)
        else:
            self.opts = add_default(None, self.opts_default)

        print(
            "\tComputing BECs: {:s}".format(
                "yes" if self.opts["BEC"]["compute"] else "no"
            )
        )

        if self.opts["BEC"]["file"] is not None:
            file = self.opts["BEC"]["file"]
            print("\tBECs will be saved to file '{:s}'".format(file))
            if self.opts["restart"]:
                print("\t'restart' is true: removing old file '{:s}'".format(file))

        if self.opts["dipole"]["file"] is not None:
            file = self.opts["BEC"]["file"]
            print("\tdipole will be saved to file '{:s}'".format(file))
            if self.opts["restart"]:
                print("\t'restart' is true: removing old file '{:s}'".format(file))

        self.count = 0
        print("\tInitialization completed.")
        pass

    def __call__(self, cell, pos):
        """Get energies, forces, stresses and extra quantities"""

        print("\n@calling 'driverdipole_driver': step {:d}".format(self.count + 1))

        # Check that if 'cell' has some np.inf values (isolated system)
        # the all the other elements are zero
        has_inf_values = np.any(np.isinf(cell))
        if has_inf_values:
            # print("The array contains inf values.")
            non_inf_mask = np.logical_not(np.isinf(cell))
            # Check if all non-inf values are zero
            all_non_inf_are_zero = np.all(cell[non_inf_mask] == 0)
            if not all_non_inf_are_zero:
                raise ValueError(
                    "Error with 'cell': the the are both inf and non-zero values.\nIs this the cell of an isolated system?"
                )

        # For isolated systems the diagonal elements of the cell are np.inf
        # we need to replace them to be sure that 'vir' will be zero
        cell[np.isinf(cell)] = 0.0

        # Get vanishing pot, forces and vir
        pot, force, vir, extras = super().__call__(cell, pos)
        extras = {}

        # computing BEC tensors and dipole
        if self.opts["BEC"]["compute"]:
            dipole, bec, _ = self.model.get_value_and_jac(cell=cell, pos=pos)

            # saving dipole to txt file
            print("dipole:")
            print(str(dipole.numpy().flatten()))

            # add BEC to extras
            if self.opts["BEC"]["send"]:
                extras["BEC"] = bec.tolist()

            # saving BEC to txt file
            if self.opts["BEC"]["file"] is not None:
                with open(self.opts["BEC"]["file"], "a") as f:
                    np.savetxt(
                        f, bec.numpy(), header="step {:d}".format(self.count), fmt=fmt
                    )

            # saving BEC to screen
            print("BEC:")
            print(str(bec.numpy()))

        else:
            # computing only dipole
            dipole, _ = self.model.get(cell=cell, pos=pos)

            # saving dipole to screen
            print("dipole:")
            print(str(dipole.numpy().flatten()))

        # saving dipole to txt file
        if self.opts["dipole"]["file"] is not None:
            with open(self.opts["dipole"]["file"], "a") as f:
                np.savetxt(f, dipole.numpy().reshape((1, 3)), fmt=fmt)

        # add dipole to extras
        if self.opts["dipole"]["send"]:
            extras["dipole"] = dipole.tolist()

        # increment the counter
        self.count += 1

        # return dipole (and BEC)
        return pot, force, vir, json.dumps(extras)
