# this is how the driver will be referred to in the input files
__DRIVER_NAME__ = "dummy"
__DRIVER_CLASS__ = "Zero_driver"

import json
from ipi.utils.messages import warning, verbosity


class Dummy_driver(object):
    """Base class providing the structure of a PES for the python driver.

    Command line:
        i-pi-py_driver -m dummy [...]
    Init arguments:
        :param verbose: bool to determine whether the PES should output verbose info.
    """

    # whether `compute` also accepts the arrays of the i-PI array backend (e.g.
    # torch tensors on a GPU), returning results of the same kind, rather than
    # only numpy arrays. Such inputs must not be modified.
    backend_arrays = False

    def __init__(self, verbose=False, *args, **kwargs):
        """Initialized dummy drivers"""
        self.verbose = verbose
        self.args = args
        self.kwargs = kwargs

        self.check_parameters()

    def check_parameters(self):
        """Dummy function that checks the arguments required to run the driver"""
        pass

    def compute_structure(self, cell, pos):
        """Evaluate a single structure"""
        # just evaluates zeros
        pot = 0.0
        force = pos * 0.0  # makes a zero force with same shape as pos
        vir = cell * 0.0  # makes a zero virial with same shape as cell
        extras = json.dumps(
            {"dipole": [0.0, 0.0, 0.0]}
        )  # have json formatting to potentially work with some test examples. meaningless value
        return pot, force, vir, extras

    def compute(self, cell, pos):
        """Does nothing, but returns properties that can be used by the driver loop."""

        if isinstance(cell, list):
            if not isinstance(pos, list) or len(cell) != len(pos):
                raise ValueError(
                    "Both position and cell should be given as lists to run in batched mode"
                )
            warning("Batched execution will execute in serial.", verbosity.high)
            return [self.compute_structure(cell, pos) for cell, pos in zip(cell, pos)]
        else:
            return self.compute_structure(cell, pos)

    def __call__(self, cell, pos):
        """Function interface"""

        return self.compute(cell, pos)


class Zero_driver(Dummy_driver):
    """The "dummy" PES: zero energy, forces and virial, computed with the same
    kind of arrays as the positions, i.e. also for the arrays of the i-PI array
    backend.

    Command line:
        i-pi-py_driver -m dummy [...]
    """

    backend_arrays = True

    def compute_structure(self, cell, pos):
        """Evaluate a single structure, with a zero energy of the same kind as pos"""
        pot, force, vir, extras = super().compute_structure(cell, pos)
        return force.sum(), force, vir, extras
