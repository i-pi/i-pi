"""Tests for server-side fixed-electric-field FFDielectric support."""

import threading
from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.forcefields import FFDielectric, ForceField, ForceRequest
from ipi.engine.properties import Properties
from ipi.inputs.forcefields import InputFFDielectric
from ipi.utils.io.inputs.io_xml import xml_parse_string
from ipi.utils.inputvalue import ArrayFromDict, InputValueFromDict
from ipi.utils.units import unit_to_internal


def _server_dielectric(**kwargs):
    """Create the smallest server-side dielectric wrapper for unit tests."""
    electric_fields = kwargs.pop("electric_fields", [])
    return FFDielectric(
        name="dielectric",
        dipole={"family": "electric-dipole", "units": "atomic_unit", "key": "dipole"},
        bec={"family": "charge", "units": "e", "key": "BEC"},
        piezo={
            "family": "electric-polarization",
            "units": "atomic_unit",
            "key": "piezoelectric",
        },
        electric_fields=electric_fields,
        forcefield=ForceField(),
        **kwargs,
    )


def test_dielectric_is_server_only_for_now():
    """Client-side coupling belongs to the follow-up PR."""
    with pytest.raises(ValueError, match="where='server'"):
        _server_dielectric(where="client")


def test_dielectric_rejects_wrapped_forcefield():
    """The wrapper must receive unwrapped coordinates."""
    with pytest.raises(ValueError, match="does not support pbc='True'"):
        FFDielectric(
            name="dielectric",
            dipole={
                "family": "electric-dipole",
                "units": "atomic_unit",
                "key": "dipole",
            },
            bec={"family": "charge", "units": "e", "key": "BEC"},
            piezo={
                "family": "electric-polarization",
                "units": "atomic_unit",
                "key": "piezoelectric",
            },
            electric_fields=[],
            forcefield=ForceField(dopbc=True),
        )


def test_input_ffdielectric_parses_server_e_configuration():
    xml = """
    <ffdielectric name="mu+E" where="server">
      <electric_field name="static">
        <parameters>{"amplitude": [1.0, 0.0, 0.0]}</parameters>
      </electric_field>
      <ffdirect name="mace"><pes>dummy</pes></ffdirect>
    </ffdielectric>
    """
    input_ff = InputFFDielectric()
    input_ff.parse(xml_parse_string(xml).fields[0][1])
    dielectric = input_ff.fetch()
    assert dielectric.name == "mu+E"
    assert dielectric.where == "server"
    assert len(dielectric.electric_fields) == 1
    np.testing.assert_allclose(dielectric.electric_fields[0].get(0.0), [1.0, 0.0, 0.0])
    assert dielectric.forcefield.pes == "dummy"


def test_response_extractors_round_trip_units_and_defaults():
    extractor_type = InputValueFromDict.specialize("electric-dipole", "eang", "dipole")
    extractor = extractor_type()
    extractor.store({"family": "electric-dipole", "units": "eang", "key": "dipole"})
    extractor._explicit = True
    assert extractor.fetch() == {
        "family": "electric-dipole",
        "units": "eang",
        "key": "dipole",
    }

    array = ArrayFromDict("electric-field", "V/ang", "field")
    np.testing.assert_allclose(
        array.get({"field": [1.0, 0.0, 0.0]}),
        [unit_to_internal("electric-field", "V/ang", 1.0), 0.0, 0.0],
    )
    assert repr(array).startswith("ArrayFromDict(")


def test_dielectric_applies_fixed_e_energy_force_and_virial():
    """Fixed E uses the driver dipole, BEC, and proper piezoelectric tensor."""
    dielectric = _server_dielectric()
    field = np.array([2.0, 0.0, 0.0])
    piezoelectric = np.zeros((3, 3, 3))
    piezoelectric[0, 0, 0] = 4.0
    request = ForceRequest(
        {
            "electric_field": field,
            "cell": (np.diag([2.0, 1.0, 1.0]), np.eye(3)),
            "result": (
                10.0,
                np.zeros(3),
                np.zeros((3, 3)),
                {
                    "dipole": np.array([0.5, 0.0, 0.0]),
                    "BEC": np.array(
                        [[[3.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]
                    ),
                    "piezoelectric": piezoelectric,
                },
            ),
        }
    )

    energy, forces, virial, _ = dielectric.fixed_E(request)

    assert energy == pytest.approx(9.0)
    np.testing.assert_allclose(forces, [6.0, 0.0, 0.0])
    np.testing.assert_allclose(
        virial, [[16.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    )

    dielectric.electric_fields = [object()]
    dielectric._field_cache = (field, None)
    properties = Properties()
    properties.simul = SimpleNamespace(fflist={"dielectric": dielectric})
    properties._threadlock = threading.Lock()
    np.testing.assert_allclose(properties["electric_field(dielectric)"][0], field)
    np.testing.assert_allclose(
        properties["electric_displacement(dielectric)"][0],
        [2.0 + np.pi, 0.0, 0.0],
    )


def test_dielectric_property_rejects_invalid_forcefield_names():
    properties = Properties()
    properties.simul = SimpleNamespace(fflist={})
    properties._threadlock = threading.Lock()
    with pytest.raises(ValueError, match="Specify the FFDielectric"):
        properties.get_dielectric_electric_field()
    with pytest.raises(ValueError, match="No forcefield named"):
        properties.get_dielectric_electric_field("missing")

    properties.simul.fflist["plain"] = ForceField()
    with pytest.raises(ValueError, match="not an FFDielectric"):
        properties.get_dielectric_electric_field("plain")


def test_dielectric_field_evaluation_is_cached_and_sums_fields():
    class Field:
        def __init__(self, value):
            self.value = np.asarray(value)
            self.calls = 0

        def get(self, time):
            self.calls += 1
            return self.value + time

    first = Field([1.0, 0.0, 0.0])
    second = Field([0.0, 2.0, 0.0])
    dielectric = _server_dielectric(electric_fields=[first, second])
    np.testing.assert_allclose(dielectric._evaluate_fields(0.5)[0], [2.0, 3.0, 1.0])
    dielectric._evaluate_fields(0.5)
    assert first.calls == second.calls == 1
    np.testing.assert_allclose(dielectric._evaluate_fields(1.0)[0], [3.0, 4.0, 2.0])
    assert first.calls == second.calls == 2
    dielectric._field_cache = (np.array([1.0, 2.0, 3.0]), None)
    np.testing.assert_allclose(dielectric.get_electric_field(), [1.0, 2.0, 3.0])


def test_dielectric_field_configuration_errors():
    dielectric = _server_dielectric()
    with pytest.raises(ValueError, match="has no electric field"):
        dielectric.get_electric_field()
    with pytest.raises(ValueError, match="requires an electric_field"):
        dielectric.apply_ensemble({"result": None})
    with pytest.raises(NotImplementedError):
        dielectric.update()
    with pytest.raises(NotImplementedError):
        dielectric.poll()


def test_dielectric_expands_voigt_piezoelectric_tensor():
    """Voigt piezoelectric data expands to symmetric Cartesian strain axes."""
    voigt = np.arange(18).reshape(3, 6)
    tensor = FFDielectric._coerce_piezoelectric_tensor(voigt)
    np.testing.assert_array_equal(tensor[:, 1, 2], voigt[:, 3])
    np.testing.assert_array_equal(tensor[:, 0, 2], voigt[:, 4])
    np.testing.assert_array_equal(tensor[:, 0, 1], voigt[:, 5])
    np.testing.assert_array_equal(tensor, tensor.swapaxes(1, 2))


def test_dielectric_rejects_nonsymmetric_piezoelectric_tensor():
    """The proper piezoelectric tensor must be symmetric in strain axes."""
    tensor = np.zeros((3, 3, 3))
    tensor[0, 1, 2] = 1.0
    with pytest.raises(ValueError, match="symmetric in its two strain indices"):
        FFDielectric._coerce_piezoelectric_tensor(tensor)
