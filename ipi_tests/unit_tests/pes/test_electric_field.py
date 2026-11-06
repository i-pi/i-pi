"""Tests for the built-in FFDielectric vector-field functions."""

import numpy as np
import pytest

from ipi.engine.motion.driven_dynamics import PythonVectorField
from ipi.inputs.motion.driven_dynamics import (
    InputFieldParameters,
    InputPythonVectorField,
)
from ipi.pes.electric_field import gaussian, plane_wave, plane_wave_gaussian, ramp


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


# Codecov: covers Gaussian and Gaussian-modulated plane-wave profiles.
def test_gaussian_and_plane_wave_gaussian_profiles():
    amplitude = np.array([1.0, -2.0, 0.5])
    assert np.allclose(gaussian(1.0, amplitude, sigma=2.0, peak=1.0), amplitude)
    expected = plane_wave(0.5, amplitude, omega=2.0) * np.exp(-0.5 * (0.5 / 2.0) ** 2)
    np.testing.assert_allclose(
        plane_wave_gaussian(0.5, amplitude, omega=2.0, sigma=2.0), expected
    )


@pytest.mark.parametrize(
    "function, kwargs, message",
    [
        (gaussian, {"amplitude": [1, 0, 0], "sigma": 0.0}, "sigma"),
        (ramp, {"amplitude": [1, 0, 0], "period": 0.0}, "period"),
    ],
)
# Codecov: covers validation of nonpositive Gaussian/ramp time scales.
def test_profiles_reject_nonpositive_time_scales(function, kwargs, message):
    with pytest.raises(ValueError, match=message):
        function(0.0, **kwargs)


# Codecov: covers custom vector-field import and callable/shape errors.
def test_python_vector_field_reports_invalid_user_functions(tmp_path):
    missing = tmp_path / "missing.py"
    with pytest.raises(ValueError, match="does not exist"):
        PythonVectorField(missing, "field", "electric-field")

    source = tmp_path / "field.py"
    source.write_text("field = 1\n")
    with pytest.raises(ValueError, match="not callable"):
        PythonVectorField(source, "field", "electric-field")

    source_valid = tmp_path / "valid.py"
    source_valid.write_text("def field(time):\n    return [time, 0]\n")
    field = PythonVectorField(source_valid, "field", "electric-field")
    with pytest.raises(ValueError, match="shape"):
        field.get(0.0)


# Codecov: covers InputPythonVectorField store/fetch round trips.
def test_python_vector_field_input_round_trip():
    field = PythonVectorField(
        file="",
        name="static",
        family="electric-field",
        units="atomic_unit",
        parameters={"amplitude": [1.0, 2.0, 3.0]},
    )
    input_field = InputPythonVectorField.specialize("electric-field")()
    input_field.store(field)
    restored = input_field.fetch()
    assert restored.name == field.name
    assert restored.family == field.family
    assert restored.parameters == field.parameters


# Codecov: covers vector-field family and JSON-parameter validation.
def test_python_vector_field_input_validates_family_and_parameters():
    with pytest.raises(ValueError, match="Unknown unit family"):
        InputPythonVectorField.specialize("not-a-family")

    parameters = InputFieldParameters()
    parameters.store({"amplitude": [1, 2, 3]})
    assert parameters.fetch() == {"amplitude": [1, 2, 3]}
    parameters.store("[1, 2, 3]")
    with pytest.raises(ValueError, match="JSON dictionary"):
        parameters.fetch()
