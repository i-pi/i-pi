import numpy as np
import gc
import weakref
from types import SimpleNamespace

import pytest

from ipi.engine.forcefields import FFPlumed, ForceRequest


class _Plumed:
    def cmd(self, _command, *_args):
        pass


def test_mtd_update_fallback_uses_force_request():
    ff = FFPlumed.__new__(FFPlumed)
    ff.lastq = np.zeros(3)
    ff.plumed_step = 0
    ff.compute_work = False
    ff.plumed = _Plumed()
    captured = []

    def evaluate(request):
        captured.append(request)
        ff.lastq[:] = request["pos"]

    ff.evaluate = evaluate

    assert ff.mtd_update(np.ones(3), np.eye(3)) == 0.0
    assert len(captured) == 1
    assert isinstance(captured[0], ForceRequest)


def test_evaluate_retains_plumed_force_buffer_until_update():
    class PointerProbe(_Plumed):
        def cmd(self, command, *args):
            if command == "setForces":
                self.buffer = weakref.ref(args[0])

    ff = FFPlumed.__new__(FFPlumed)
    ff.natoms = 1
    ff.lastq = np.zeros(3)
    ff.plumed_step = 0
    ff.charges = np.zeros(1)
    ff.masses = np.ones(1)
    ff.system_force = None
    ff.plumed_data = {}
    ff.plumed = PointerProbe()
    request = ForceRequest(
        {"pos": np.zeros(3), "cell": (np.eye(3), None), "result": None}
    )
    ff.evaluate(request)
    gc.collect()
    retained = ff.plumed.buffer()
    assert retained is not None
    retained[:] = 1.0
    # Subsequent PLUMED work must not mutate the already returned MD force.
    assert np.array_equal(request["result"][1], np.zeros(3))


@pytest.mark.parametrize("has_system_force", [False, True])
@pytest.mark.parametrize("bias_scale", [0.0, 0.5, -2.0])
def test_evaluate_returns_only_bias_correction(has_system_force, bias_scale):
    """Remove physical contributions after converting PLUMED's virial sign."""
    physical_force = np.array([[1.0, -2.0, 3.0], [-4.0, 5.0, -6.0]])
    physical_virial = np.array(
        [[2.0, -1.0, 0.5], [-1.0, 3.0, -0.25], [0.5, -0.25, 4.0]]
    )
    bias_force = bias_scale * np.array([[2.0, 1.0, -3.0], [-2.0, -1.0, 3.0]])
    bias_virial = bias_scale * np.array(
        [[1.0, 0.25, -0.5], [0.25, -2.0, 0.75], [-0.5, 0.75, 3.0]]
    )
    bias_energy = bias_scale * 1.25

    class IncrementProbe(_Plumed):
        """Emulate PLUMED's additive force and opposite-sign virial contract."""

        def cmd(self, command, *args):
            if command == "setForces":
                self.force = args[0]
            elif command == "setVirial":
                self.virial = args[0]
            elif command == "setEnergy":
                assert args[0] == 7.0
            elif command == "performCalcNoUpdate":
                np.testing.assert_array_equal(
                    self.force, physical_force if has_system_force else 0.0
                )
                np.testing.assert_array_equal(
                    self.virial, -physical_virial if has_system_force else 0.0
                )
                self.force[:] += bias_force
                self.virial[:] -= bias_virial
            elif command == "getBias":
                args[0][:] = bias_energy

    ff = FFPlumed.__new__(FFPlumed)
    ff.natoms = 2
    ff.lastq = np.zeros(6)
    ff.plumed_step = 0
    ff.charges = np.zeros(2)
    ff.masses = np.ones(2)
    ff.system_force = (
        SimpleNamespace(f=physical_force.copy(), vir=physical_virial.copy(), pot=7.0)
        if has_system_force
        else None
    )
    ff.plumed_data = {}
    ff.plumed = IncrementProbe()
    request = ForceRequest(
        {"pos": np.zeros(6), "cell": (np.eye(3), None), "result": None}
    )
    ff.evaluate(request)
    energy, force, virial, _ = request["result"]
    assert energy == bias_energy
    np.testing.assert_array_equal(force, bias_force.flatten())
    np.testing.assert_array_equal(virial, bias_virial)
    assert request["status"] == "Done"
    assert request._event_done.is_set()
    if has_system_force:
        np.testing.assert_array_equal(ff.system_force.f, physical_force)
        np.testing.assert_array_equal(ff.system_force.vir, physical_virial)


@pytest.mark.parametrize("stiffness", [0.5, 2.0])
def test_energy_bias_correction_matches_coordinate_and_strain_derivatives(stiffness):
    """Check an energy-dependent bias against independent energy differences."""
    positions = np.array([[0.5, -1.0, 1.5], [-0.5, 1.0, -1.5]])
    physical_force = -positions
    physical_virial = positions.T @ physical_force
    physical_energy = 0.5 * np.sum(positions**2)

    def bias_energy(q):
        """Harmonic restraint on the physical harmonic potential energy."""
        energy = 0.5 * np.sum(q**2)
        return 0.5 * stiffness * (energy - 1.0) ** 2

    class EnergyProbe(_Plumed):
        """Apply the force rescaling used by PLUMED's ENERGY action."""

        def cmd(self, command, *args):
            if command == "setForces":
                self.force = args[0]
            elif command == "setVirial":
                self.virial = args[0]
            elif command == "setEnergy":
                self.energy = args[0]
            elif command == "performCalcNoUpdate":
                factor = 1.0 + stiffness * (self.energy - 1.0)
                self.force[:] *= factor
                self.virial[:] *= factor
            elif command == "getBias":
                args[0][:] = bias_energy(positions)

    ff = FFPlumed.__new__(FFPlumed)
    ff.natoms = 2
    ff.lastq = np.zeros(6)
    ff.plumed_step = 0
    ff.charges = np.zeros(2)
    ff.masses = np.ones(2)
    ff.system_force = SimpleNamespace(
        f=physical_force.copy(), vir=physical_virial.copy(), pot=physical_energy
    )
    ff.plumed_data = {}
    ff.plumed = EnergyProbe()
    request = ForceRequest(
        {"pos": positions.flatten(), "cell": (np.eye(3), None), "result": None}
    )
    ff.evaluate(request)
    energy, force, virial, _ = request["result"]
    assert energy == bias_energy(positions)
    step = 1.0e-5
    force_fd = np.zeros_like(positions)
    virial_fd = np.zeros((3, 3))
    for index in np.ndindex(positions.shape):
        delta = np.zeros_like(positions)
        delta[index] = step
        force_fd[index] = -(
            bias_energy(positions + delta) - bias_energy(positions - delta)
        ) / (2.0 * step)
    for index in np.ndindex(virial_fd.shape):
        strain = np.zeros((3, 3))
        strain[index] = step
        virial_fd[index] = -(
            bias_energy(positions @ (np.eye(3) + strain))
            - bias_energy(positions @ (np.eye(3) - strain))
        ) / (2.0 * step)
    np.testing.assert_allclose(force.reshape((-1, 3)), force_fd, rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(virial, virial_fd, rtol=1e-8, atol=1e-9)
