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


def static(time, amplitude):
    """Return a time-independent vector with the requested amplitude."""

    del time
    return np.asarray(amplitude)


def plane_wave(
    time,
    amplitude,
    omega,
    phase=0.0,
    omega_units="atomic_unit",
):
    """Return ``amplitude * sin(omega * time + phase)``.

    ``omega`` is a angular omega, such as ``100`` with
    ``omega_units="GHz"``. The returned vector remains in the units
    specified by the enclosing ``<electric_field units="...">`` element.
    """
    omega = unit_to_internal("frequency", omega_units, omega)
    return np.asarray(amplitude) * np.sin(omega * time + phase)


def ramp(
    time,
    amplitude,
    period,
    phase=0.0,
    period_units="atomic_unit",
):
    """Return a periodic triangular ramp with a prescribed period.

    Over one cycle the scalar profile follows ``0 -> 1 -> 0 -> -1 -> 0``:
    it ramps up during the first quarter-cycle, down during the next two, and
    up during the final quarter-cycle. ``period`` is the duration of one
    cycle, interpreted in ``period_units``. ``phase`` is accepted but does
    not affect the ramp.
    """
    del phase
    period = unit_to_internal("time", period_units, period)
    if period <= 0:
        raise ValueError("The ramp period must be positive.")
    angle = 2.0 * np.pi * time / period
    profile = 2.0 / np.pi * np.arcsin(np.sin(angle))
    return np.asarray(amplitude) * profile


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
    omega,
    sigma,
    peak=0.0,
    phase=0.0,
    omega_units="atomic_unit",
    sigma_units="atomic_unit",
    peak_units="atomic_unit",
):
    """Return a plane wave multiplied by a Gaussian envelope."""

    return plane_wave(
        time,
        amplitude,
        omega,
        phase=phase,
        omega_units=omega_units,
    ) * _gaussian_envelope(
        time,
        sigma,
        peak,
        sigma_units=sigma_units,
        peak_units=peak_units,
    )
