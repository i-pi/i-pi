"""Creates objects that handle normal mode transformations."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from copy import copy

import numpy as np

from ipi.engine.normalmodes import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *


__all__ = ['InputNormalModes', 'InputNMFrequencies']


class InputNMFrequencies(InputArray):

    """ Storage class for NormalModes engine. """

    attribs = copy(InputArray.attribs)
    attribs["style"] = (InputAttribute, {"dtype": str,
                                         "default": "rpmd",
                                         "help": """Specifies the technique to be used to calculate the dynamical masses.
                                                'rpmd' simply assigns the bead masses the physical mass.
                                                'manual' sets all the normal mode frequencies except the centroid normal mode manually.
                                                'pa-cmd' takes an argument giving the frequency to set all the non-centroid normal modes to.
                                                'wmax-cmd' is similar to 'pa-cmd', except instead of taking one argument it takes two
                                                      ([wmax,wtarget]). The lowest-lying normal mode will be set to wtarget for a
                                                      free particle, and all the normal modes will coincide at frequency wmax. """,
                                         "options": ['pa-cmd', 'wmax-cmd', 'manual', 'rpmd']})

    default_label = "NMFREQUENCIES"
    default_help = "Provides a compact way of specifying the ring polymer frequencies"

    def __init__(self, help=None, dimension=None, default=None, dtype=None):
        """ Initializes InputNormalModes.

        Just calls the parent initialization function with appropriate arguments.
        """

        super(InputNMFrequencies, self).__init__(help=help, default=default, dtype=float, dimension="frequency")

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


class InputNormalModes(Input):

    """ Storage class for NormalModes engine.

    Describes how normal-modes transformation and integration should be
    performed.

    Attributes:
        frequencies: Specifies how the frequencies given should be interpreted
            when creating the mass matrix.
        transform: Specifies whether the normal mode calculation will be
            done using a FFT transform or a matrix multiplication.
    """

    attribs = {
        "transform": (InputValue, {"dtype": str,
                                   "default": "fft",
                                   "help": "Specifies whether to calculate the normal mode transform using a fast Fourier transform or a matrix multiplication. For small numbers of beads the matrix multiplication may be faster.",
                                   "options": ['fft', 'matrix']})
    }

    fields = {
        "frequencies": (InputNMFrequencies, {
            "default": ("rpmd", np.zeros(0)),
                     "help": "Specifies normal mode frequencies for a (closed path) calculation"
        }),
        "open_paths": (InputArray, {"dtype": int,
                                    "default": np.zeros(0, int),
                                    "help": "Indices of the atmoms whose path should be opened (zero-based)."})
    }

    default_label = "NORMALMODES"
    default_help = "Deals with the normal mode transformations, including the adjustment of bead masses to give the desired ring polymer normal mode frequencies if appropriate. Takes as arguments frequencies, of which different numbers must be specified and which are used to scale the normal mode frequencies in different ways depending on which 'mode' is specified."

    def store(self, nm):
        self.transform.store(nm.transform_method)
        self.frequencies.store((nm.mode, nm.nm_freqs))
        self.open_paths.store(nm.open_paths)

    def fetch(self):
        mode, freqs = self.frequencies.fetch()
        return NormalModes(mode, self.transform.fetch(), freqs, open_paths=self.open_paths.fetch())
