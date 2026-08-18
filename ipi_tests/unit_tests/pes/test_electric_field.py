"""Tests for the built-in FFDielectric vector fields."""

import numpy as np

from ipi.pes.electric_field import gaussian, plane_wave
from ipi.utils.units import unit_to_internal


def test_plane_wave_frequency_in_ghz():
    cycles_per_au = unit_to_internal("frequency-cyclic", "GHz", 1.0)
    quarter_period = 1.0 / (4.0 * cycles_per_au)
    value = plane_wave(
        quarter_period,
        [0.0, 0.0, 1.0],
        1.0,
        frequency_units="GHz",
    )
    np.testing.assert_allclose(value, [0.0, 0.0, 1.0], atol=1.0e-12)


def test_gaussian_converts_time_parameters_to_atomic_units():
    one_fs = unit_to_internal("time", "femtosecond", 1.0)
    value = gaussian(
        one_fs,
        [1.0, 2.0, 3.0],
        1.0,
        peak=1.0,
        sigma_units="femtosecond",
        peak_units="femtosecond",
    )
    np.testing.assert_allclose(value, [1.0, 2.0, 3.0])
