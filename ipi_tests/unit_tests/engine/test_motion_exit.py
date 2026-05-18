import pytest

from ipi.engine.motion import Motion
from ipi.engine.simulation import Simulation


class _System:
    def __init__(self, motion):
        self.motion = motion


class _FinishingMotion(Motion):
    def step(self, step=None):
        self.finish(status="success", message="finished cleanly")


def _make_simulation(motion):
    simulation = Simulation.__new__(Simulation)
    simulation.syslist = [_System(motion)]
    simulation.threading = False
    simulation.smotion = None
    simulation.finished = False
    simulation.exit_status = "success"
    simulation.exit_message = " @ SIMULATION: Exiting cleanly."
    return simulation


def test_motion_finish_marks_simulation_finished():
    simulation = _make_simulation(_FinishingMotion())

    simulation.run_step(0)

    assert simulation.finished
    assert simulation.exit_status == "success"
    assert simulation.exit_message == "finished cleanly"


def test_motion_finish_rejects_unknown_status():
    motion = Motion()

    with pytest.raises(ValueError):
        motion.finish(status="unknown")
