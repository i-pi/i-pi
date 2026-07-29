"""Simple vector-field functions for :class:`FFDielectric` inputs.

All functions receive time in the units selected by the XML ``time_units``
attribute and return a Cartesian vector in the XML ``units``. Frequencies are
angular frequencies (radians per selected time unit), phases are in radians,
and Gaussian ``sigma`` and ``peak`` values use the selected time unit.
"""

import numpy as np


def _gaussian_envelope(time, sigma, peak):
    if sigma <= 0:
        raise ValueError("The Gaussian sigma must be positive.")
    return np.exp(-0.5 * ((time - peak) / sigma) ** 2)


def static(time, amplitude):
    """Returns a time-independent vector with the requested amplitude."""

    del time
    return np.asarray(amplitude)


def plane_wave(time, amplitude, frequency, phase=0.0):
    """Returns ``amplitude * cos(frequency * time + phase)``."""

    return np.asarray(amplitude) * np.cos(frequency * time + phase)


def gaussian(time, amplitude, sigma, peak=0.0):
    """Returns a vector with a Gaussian time profile."""

    return np.asarray(amplitude) * _gaussian_envelope(time, sigma, peak)


def plane_wave_gaussian(time, amplitude, frequency, sigma, peak=0.0, phase=0.0):
    """Returns a plane wave multiplied by a Gaussian envelope."""

    return plane_wave(time, amplitude, frequency, phase) * _gaussian_envelope(
        time, sigma, peak
    )
