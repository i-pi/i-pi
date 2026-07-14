"""Ideal-gas (non-interacting) potential: a vectorized no-op driver.

Returns zero potential, force and virial for every structure. Useful as a
minimal-cost benchmark target — e.g. to measure socket/dispatch overhead,
including batched (batch_size>1) evaluation — where the driver does no real
work and the wall time is dominated by communication.
"""

from .dummy import Dummy_driver
import numpy as np

__DRIVER_NAME__ = "gas"
__DRIVER_CLASS__ = "Gas_driver"


class Gas_driver(Dummy_driver):
    """Non-interacting particles: zero potential, force and virial.

    Command-line:
        i-pi-py_driver -m gas -u [...]
    """

    def compute_structure(self, cell, pos):
        return 0.0, np.zeros_like(pos), np.zeros((3, 3)), ""

    def compute(self, cell, pos):
        """Vectorized no-op for a batch of structures (single allocation), with
        a fallback to the single-structure path."""

        if not isinstance(cell, list):
            return self.compute_structure(cell, pos)
        if not isinstance(pos, list) or len(cell) != len(pos):
            raise ValueError(
                "Both position and cell should be given as lists to run in batched mode"
            )

        # nothing to compute: hand back zero-valued views of one allocation
        nbatch = len(pos)
        forces = np.zeros((nbatch,) + np.asarray(pos[0]).shape)
        vir = np.zeros((3, 3))
        return [(0.0, forces[i], vir, "") for i in range(nbatch)]
