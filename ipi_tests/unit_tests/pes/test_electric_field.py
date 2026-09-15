"""Tests for the built-in FFDielectric vector-field functions."""

import numpy as np

from ipi.engine.motion.driven_dynamics import PythonVectorField
from ipi.pes.electric_field import plane_wave, ramp


def test_plane_wave_starts_at_zero_without_a_phase():
    """The default plane wave has zero field at the start of the simulation."""
    assert np.allclose(plane_wave(0.0, [1.0, -2.0, 0.5], omega=1.0), 0.0)


def test_ramp_has_four_linear_segments_per_cycle():
    """The periodic ramp follows up, down, down, then up."""
    amplitude = np.array([1.0, -2.0, 0.5])
    period = 4.0
    times = period / 4 * np.arange(5)

    values = np.array([ramp(time, amplitude, period=period) for time in times])

    expected_profile = np.array([0.0, 1.0, 0.0, -1.0, 0.0])
    assert np.allclose(values, expected_profile[:, None] * amplitude)


def test_ramp_accepts_period_parameters_via_ffdielectric_field_loader():
    """The built-in ramp can be selected with a period and period units."""
    field = PythonVectorField(
        file="",
        name="ramp",
        family="electric-field",
        units="atomic_unit",
        parameters={
            "amplitude": [0.0, 0.0, 2.0],
            "period": 4.0,
            "phase": np.pi / 2,
        },
    )

    assert np.allclose(field.get(0.0), [0.0, 0.0, 0.0])
    assert np.allclose(field.get(1.0), [0.0, 0.0, 2.0])
