"""Interface with e3nn to run machine learned electric dipole"""

import sys
import json
import numpy as np
import warnings
import importlib
from .dummy import Dummy_driver


def recursive_copy(source_dict: dict, target_dict: dict) -> dict:
    """
    Recursively copy keys and values from a source dictionary to a target dictionary, if they are not present in the target.

    This function takes two dictionaries, 'source_dict' and 'target_dict', and copies keys and values from 'source_dict' to 'target_dict'. If a key exists in both dictionaries and both values are dictionaries, the function recursively calls itself to copy nested keys and values. If a key does not exist in 'target_dict', it is added along with its corresponding value from 'source_dict'.

    Args:
        source_dict (dict): The source dictionary containing keys and values to be copied.
        target_dict (dict): The target dictionary to which keys and values are copied if missing.

    Returns:
        dict: The modified 'target_dict' with keys and values copied from 'source_dict'.

    Example:
        >>> dict_A = {"a": 1, "b": {"b1": 2, "b2": {"b2_1": 3}}, "c": 4}
        >>> dict_B = {"a": 10, "b": {"b1": 20, "b2": {"b2_2": 30}}, "d": 40}
        >>> result = recursive_copy(dict_A, dict_B)
        >>> print(result)
        {'a': 10, 'b': {'b1': 20, 'b2': {'b2_1': 3, 'b2_2': 30}}, 'd': 40}
    """
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
    """
    Add default key-value pairs to a dictionary if they are not present.

    This function takes two dictionaries: 'dictionary' and 'default'. It checks each key in the 'default' dictionary, and if the key is not already present in the 'dictionary', it is added along with its corresponding value from the 'default' dictionary. If 'dictionary' is not provided, an empty dictionary is used as the base.

    Args:
        dictionary (dict, optional): The input dictionary to which default values are added. If None, an empty dictionary is used. Default is None.
        default (dict): A dictionary containing the default key-value pairs to be added to 'dictionary'.

    Returns:
        dict: The modified 'dictionary' with default values added.

    Raises:
        ValueError: If 'dictionary' is not of type 'dict'.

    Example:
        >>> existing_dict = {'a': 1, 'b': 2}
        >>> default_values = {'b': 0, 'c': 3}
        >>> result = add_default(existing_dict, default_values)
        >>> print(result)
        {'a': 1, 'b': 2, 'c': 3}
    """
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

    if type(instructions) == str:
        with open(instructions, "r") as json_file:
            _instructions = json.load(json_file)
        instructions = _instructions

    # instructions['kwargs']["normalization"] = None

    # wxtract values for the instructions
    kwargs = instructions["kwargs"]
    cls = instructions["class"]
    mod = instructions["module"]

    # get the class to be instantiated
    # Call the function and suppress the warning
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")  # , category=UserWarning)
        class_obj = get_class(mod, cls)

    # instantiate class
    # try :
    model = class_obj(**kwargs)
    if not model:
        raise ValueError(
            "Error instantiating class '{:s}' from module '{:s}'".format(cls, mod)
        )

    try:
        N = model.n_parameters()
        print("\tLoaded model has {:d} parameters".format(N))
    except:
        print("\tCannot count parameters")

    # Load the parameters from the saved file
    checkpoint = torch.load(parameters)

    # Update the model's state dictionary with the loaded parameters
    model.load_state_dict(checkpoint)
    model.eval()

    # Store the chemical species that will be used during the simulation.
    model._symbols = instructions["chemical-symbols"]

    return model


class e3nn_pol(Dummy_driver):
    opts_default = {
        "compute-BEC": True,
        "print": True,
    }

    def __init__(self, args=None):
        self.error_msg = """The parameters of 'e3nn_pol' are not correctly formatted. \
            They should be two or three strings, separated by a comma."""
        super().__init__(args)

    def check_arguments(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in librascal
        """
        try:
            arglist = self.args.split(",")
        except ValueError:
            sys.exit(self.error_msg)

        if len(arglist) >= 2:
            # self.model_file = arglist[0] # file with the torch.nn.Module
            info_file = arglist[0]  # json file to properly allocate a 'model' object
            parameters_file = arglist[1]  # *.pth file with the model parameters
            try:
                opts_file = arglist[2]  # json file with some parameters for this class
            except:
                print("\tNo options file provided: using the default values")
                opts_file = None
        else:
            sys.exit(self.error_msg)  # to be modified

        print("\tThe driver is 'e3nn_pol'")
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
            "\tComputing BECs: {:s}".format("yes" if self.opts["compute-BEC"] else "no")
        )
        print("\tInitialization completed")

        self.count = 0

        pass

    def __call__(self, cell, pos):
        """Get energies, forces, stresses and extra quantities"""

        self.count += 1

        if self.opts["print"]:
            print(" @calling 'e3nn_pol' for the {:d}th time".format(self.count))

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

        if self.opts["compute-BEC"]:
            dipole, bec, X = self.model.get_value_and_jac(cell=cell, pos=pos)
            extras["BEC"] = bec.tolist()
            # print("sum rule:",bec.sum(dim=0))
        else:
            dipole, X = self.model.get(cell=cell, pos=pos)

        extras["dipole"] = dipole.tolist()

        return pot, force, vir, json.dumps(extras)
