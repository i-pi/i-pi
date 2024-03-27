"""Creates objects that handle normal mode transformations."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from copy import copy

import numpy as np

from ipi.engine.normalmodes import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *


__all__ = ["InputNormalModes", "InputNMFrequencies", "InputBosons"]


class InputNMFrequencies(InputArray):
    """Storage class for NormalModes engine."""

    attribs = copy(InputArray.attribs)
    attribs["style"] = (
        InputAttribute,
        {
            "dtype": str,
            "default": "rpmd",
            "help": """
Specifies the technique to be used to calculate the dynamical masses.
'rpmd' simply assigns the bead masses the physical mass.
'manual' sets all the normal mode frequencies except the centroid normal mode manually.
'pa-cmd' takes an argument giving the frequency to set all the non-centroid normal modes to.
'wmax-cmd' is similar to 'pa-cmd', except instead of taking one argument it takes two
([wmax,wtarget]). The lowest-lying normal mode will be set to wtarget for a
free particle, and all the normal modes will coincide at frequency wmax. 
""",
            "options": ["pa-cmd", "wmax-cmd", "manual", "rpmd"],
        },
    )

    default_label = "NMFREQUENCIES"
    default_help = "Provides a compact way of specifying the ring polymer frequencies"

    def __init__(self, help=None, dimension=None, default=None, dtype=None):
        """Initializes InputNormalModes.

        Just calls the parent initialization function with appropriate arguments.
        """

        super(InputNMFrequencies, self).__init__(
            help=help, default=default, dtype=float, dimension="frequency"
        )

    def store(self, mf):
        """Takes a modes and frequencies ans store them
        of it.

        Args:
            mf: A tuple containing a string and an array, ("MODE", [FREQS]).
        """
        mode, freqs = mf
        super(InputNMFrequencies, self).store(freqs)
        self.style.store(mode)

    def fetch(self):
        """Creates a normal modes object.

        Returns:
            A normal modes object.
        """

        super(InputNMFrequencies, self).check()
        return (self.style.fetch(), super(InputNMFrequencies, self).fetch())


class InputBosons(InputArray):
    """Storage class for the input of which atoms have bosonic idistinguishability.

    Attributes:
       shape: The shape of the array.
    """

    attribs = copy(InputArray.attribs)
    attribs["id"] = (
        InputAttribute,
        {
            "dtype": str,
            "default": "index",
            "options": ["index", "label"],
            "help": "If 'id' is 'index', then bosonic atoms are specified a list of indices (zero-based). If 'id' is 'label' then specify a list of labels.",
        },
    )

    default_label = "BOSONS"
    default_help = "Deals with the specification of which atoms are have bosonic indistinguishability. The specified atoms participate in exchange interaction, which forms ring polymers that combine several particles together. The algorithm scales quadratically with the number of atoms and linearly with the number of beads. The implementation is based on Hirshberg et al.'s doi:10.1073/pnas.1913365116 and Feldman and Hirshberg's doi:10.1063/5.0173749."

    def __init__(self, help=None, default=None):
        """Initialises InputBosons.

        Args:
           help: A help string.
           dimension: The dimensionality of the value.
           default: A default value.
           dtype: An optional data type. Defaults to None.
        """

        super(InputBosons, self).__init__(help=help, default=default, dtype=str)

    def store(self, value):
        """Converts the data to the appropriate data type, shape and units and
        stores it.

        Args:
           value: The raw data to be stored.
           id: The specification mode of the atoms.
        """

        self.id.store("index")
        super(InputBosons, self).store(value)

    def fetch(self):
        """Returns the stored data in the user defined units."""

        value = super(InputBosons, self).fetch()

        return (value, self.id.fetch())


class InputNormalModes(Input):
    """Storage class for NormalModes engine.

    Describes how normal-modes transformation and integration should be
    performed.

    Attributes:
        frequencies: Specifies how the frequencies given should be interpreted
            when creating the mass matrix.
        transform: Specifies whether the normal mode calculation will be
            done using a FFT transform or a matrix multiplication.
    """

    attribs = {
        "transform": (
            InputAttribute,
            {
                "dtype": str,
                "default": "fft",
                "help": "Specifies whether to calculate the normal mode transform using a fast Fourier transform or a matrix multiplication. For small numbers of beads the matrix multiplication may be faster.",
                "options": ["fft", "matrix"],
            },
        ),
        "propagator": (
            InputAttribute,
            {
                "dtype": str,
                "default": "exact",
                "help": "How to propagate the free ring polymer dynamics. Cayley transform is not exact but is strongly stable and avoid potential resonance issues. A bab scheme performs numerical verlet type propagation. All three options work for distinguishable particles. Only the bab propagator can be used with bosonic particles.",
                "options": ["exact", "cayley", "bab"],
            },
        ),
        "fft_threads": (
            InputValue,
            {
                "dtype": int,
                "default": 1,
                "help": "The number of threads to be used for the FFT.",
                "dimension": None,
            },
        ),
        "fft_float32": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "Whether to use single precision FFT.",
                "dimension": None,
            },
        ),
    }

    fields = {
        "frequencies": (
            InputNMFrequencies,
            {
                "default": ("rpmd", np.zeros(0)),
                "help": "Specifies normal mode frequencies for a (closed path) calculation",
            },
        ),
        "open_paths": (
            InputArray,
            {
                "dtype": int,
                "default": np.zeros(0, int),
                "help": "Indices of the atoms whose path should be opened (zero-based).",
            },
        ),
        "bosons": (
            InputBosons,
            {
                "default": np.zeros(0, str),
                "help": "Specify which atoms are bosons.",
            },
        ),
        "nmts": (
            InputValue,
            {
                "dtype": int,
                "default": 10,
                "help": "The number of iterations to perform one bab step.",
                "dimension": None,
            },
        ),
    }

    default_label = "NORMALMODES"
    default_help = "Deals with the normal mode transformations, including the adjustment of bead masses to give the desired ring polymer normal mode frequencies if appropriate. Takes as arguments frequencies, of which different numbers must be specified and which are used to scale the normal mode frequencies in different ways depending on which 'mode' is specified."

    def store(self, nm):
        self.transform.store(nm.transform_method)
        self.fft_threads.store(nm.fft_threads)
        self.fft_float32.store(nm.fft_float32)
        self.frequencies.store((nm.mode, nm.nm_freqs))
        self.propagator.store(nm.propagator)
        self.open_paths.store(nm.open_paths)
        self.bosons.store(nm.bosons)
        self.nmts.store(nm.nmts)

    def fetch(self):
        mode, freqs = self.frequencies.fetch()
        return NormalModes(
            mode,
            self.transform.fetch(),
            self.propagator.fetch(),
            freqs,
            open_paths=self.open_paths.fetch(),
            bosons=self.bosons.fetch(),
            nmts=self.nmts.fetch(),
            fft_threads=self.fft_threads.fetch(),
            fft_float32=self.fft_float32.fetch(),
        )
