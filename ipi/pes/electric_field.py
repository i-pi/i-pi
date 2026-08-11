"""Simple vector-field functions for :class:`FFDielectric` inputs.

The callable receives simulation time in atomic units and returns a Cartesian
vector in the units selected by the XML ``units`` attribute. Frequencies are
specified as cyclic frequencies and converted to angular atomic frequencies
inside this module. Phases are in radians, while Gaussian ``sigma`` and
``peak`` values are converted from their respective ``*_units`` arguments.
"""

import numpy as np

from ipi.utils.units import unit_to_internal


def _gaussian_envelope(
    time,
    sigma,
    peak,
    sigma_units="atomic_unit",
    peak_units="atomic_unit",
):
    sigma = unit_to_internal("time", sigma_units, sigma)
    peak = unit_to_internal("time", peak_units, peak)
    if sigma <= 0:
        raise ValueError("The Gaussian sigma must be positive.")
    return np.exp(-0.5 * ((time - peak) / sigma) ** 2)


def _angular_frequency(frequency, frequency_units):
    """Convert a cyclic frequency to angular frequency in atomic units."""

    cycles = unit_to_internal("frequency-cyclic", frequency_units, frequency)
    return 2.0 * np.pi * cycles


def static(time, amplitude):
    """Return a time-independent vector with the requested amplitude."""

    del time
    return np.asarray(amplitude)


def plane_wave(
    time,
    amplitude,
    frequency,
    phase=0.0,
    frequency_units="atomic_unit",
):
    """Return ``amplitude * cos(omega * time + phase)``.

    ``frequency`` is a cyclic frequency, such as ``100`` with
    ``frequency_units="GHz"``. The returned vector remains in the units
    specified by the enclosing ``<electric_field units="...">`` element.
    """

    omega = _angular_frequency(frequency, frequency_units)
    return np.asarray(amplitude) * np.cos(omega * time + phase)


def gaussian(
    time,
    amplitude,
    sigma,
    peak=0.0,
    sigma_units="atomic_unit",
    peak_units="atomic_unit",
):
    """Return a vector with a Gaussian time profile."""

    return np.asarray(amplitude) * _gaussian_envelope(
        time,
        sigma,
        peak,
        sigma_units=sigma_units,
        peak_units=peak_units,
    )


def plane_wave_gaussian(
    time,
    amplitude,
    frequency,
    sigma,
    peak=0.0,
    phase=0.0,
    frequency_units="atomic_unit",
    sigma_units="atomic_unit",
    peak_units="atomic_unit",
):
    """Return a plane wave multiplied by a Gaussian envelope."""

    return plane_wave(
        time,
        amplitude,
        frequency,
        phase=phase,
        frequency_units=frequency_units,
    ) * _gaussian_envelope(
        time,
        sigma,
        peak,
        sigma_units=sigma_units,
        peak_units=peak_units,
    )
