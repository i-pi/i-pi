"""Python driver for the Fortran q-TIP4P/f dipole-polarizability surface.

This driver deliberately calls the original ``h2o_dipole`` Fortran routine
instead of maintaining a second implementation of its Ewald and polarization
model.  Consequently, it returns the same dipole, dipole derivatives, and
polarizability as ``i-pi-driver -m water_dip_pol`` when linked against the
same source and compiler.
"""

import ctypes
import json
import os
from pathlib import Path

import numpy as np

from ipi.pes.dummy import Dummy_driver

__DRIVER_NAME__ = "h2o_dip_pol"
__DRIVER_CLASS__ = "H2ODipPolDriver"


class H2ODipPolDriver(Dummy_driver):
    """Expose ``drivers/f90/pes/h2o_dip_pol.f90`` through the Python driver.

    The Fortran source must first be compiled as a shared library, e.g.::

        gfortran -O3 -ffree-line-length-none -fPIC -shared \
            drivers/f90/pes/h2o_dip_pol.f90 -o libh2o_dip_pol.so

    Set ``IPI_H2O_DIP_POL_LIBRARY`` to that library path, or pass it as the
    ``library`` argument.  The driver returns zero energy, forces, and virial,
    exactly like the Fortran ``water_dip_pol`` driver; its physical outputs are
    provided in the JSON extras ``dipole``, ``dipole_derivative``, and
    ``polarizability``.

    Positions must contain water molecules ordered ``O H H O H H ...`` and the
    simulation cell must be orthorhombic.
    """

    def __init__(self, library=None, compute_der=True, *args, **kwargs):
        self.library = library
        self.compute_der = bool(compute_der)
        super().__init__(*args, **kwargs)

    def check_parameters(self):
        library = self.library or os.environ.get("IPI_H2O_DIP_POL_LIBRARY")
        if library is None:
            source_tree_library = (
                Path(__file__).resolve().parents[2]
                / "drivers"
                / "f90"
                / "libh2o_dip_pol.so"
            )
            if source_tree_library.is_file():
                library = source_tree_library
            else:
                raise ValueError(
                    "The h2o_dip_pol Fortran library was not found. Compile "
                    "drivers/f90/pes/h2o_dip_pol.f90 as a shared library and "
                    "set IPI_H2O_DIP_POL_LIBRARY to its path."
                )

        library = Path(library).expanduser().resolve()
        if not library.is_file():
            raise ValueError(f"h2o_dip_pol library does not exist: {library}")

        self._library = ctypes.CDLL(str(library))
        try:
            self._h2o_dipole = self._library.h2o_dipole_
        except AttributeError as exc:
            raise ValueError(
                f"{library} does not export the gfortran symbol 'h2o_dipole_'."
            ) from exc

        double_ptr = ctypes.POINTER(ctypes.c_double)
        self._h2o_dipole.argtypes = [
            double_ptr,
            ctypes.POINTER(ctypes.c_int32),
            double_ptr,
            ctypes.POINTER(ctypes.c_int32),
            double_ptr,
            double_ptr,
            double_ptr,
        ]
        self._h2o_dipole.restype = None

    def compute_structure(self, cell, pos):
        cell = np.asarray(cell, dtype=np.float64)
        pos = np.asarray(pos, dtype=np.float64)
        if pos.ndim != 2 or pos.shape[1] != 3 or len(pos) % 3:
            raise ValueError("h2o_dip_pol expects an O-H-H ordered (3N, 3) array")
        if cell.shape != (3, 3) or not np.allclose(cell, np.diag(np.diag(cell))):
            raise ValueError("h2o_dip_pol supports orthorhombic cells only")

        nat = ctypes.c_int32(len(pos))
        compute_der = ctypes.c_int32(int(self.compute_der))
        box = np.ascontiguousarray(np.diag(cell), dtype=np.float64)
        atoms = np.asfortranarray(pos, dtype=np.float64)
        dipole = np.zeros(3, dtype=np.float64)
        dipole_derivative = np.zeros((len(pos), 9), dtype=np.float64, order="F")
        polarizability = np.zeros((3, 3), dtype=np.float64, order="F")

        self._h2o_dipole(
            box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(nat),
            atoms.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(compute_der),
            dipole.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            dipole_derivative.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            polarizability.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )

        BEC = dipole_derivative.ravel(order="C").reshape((len(pos), 3, 3))
        extras = {
            "dipole": dipole.tolist(),
            # Keep the flat Fortran-driver layout expected by FFCavPhSocket:
            # nine dipole derivatives for each atom, in atom order.
            "BEC": BEC.tolist(),
            "polarizability": polarizability.ravel(order="F").tolist(),
        }
        return 0.0, np.zeros_like(pos), np.zeros((3, 3)), json.dumps(extras)
