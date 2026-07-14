from pathlib import Path

import pytest

try:
    from ipi_tests.regression_tests.runstools import Runner_regression
except ImportError:
    from runstools import Runner_regression


RESTART_TESTS_FOLDER = Path(__file__).resolve().parent / "restart_tests"
RESTART_CASES = [
    pytest.param(
        0,
        "dynamics",
        RESTART_TESTS_FOLDER / "NVE" / "NVE_1" / "harmonic_f90",
        id="dynamics",
    ),
    pytest.param(
        1,
        "minimize",
        RESTART_TESTS_FOLDER / "GEOP" / "CG" / "harmonic",
        id="minimize",
    ),
    pytest.param(
        2,
        "vibrations",
        RESTART_TESTS_FOLDER / "PHONONS" / "fd_phonons" / "ch4hcbe",
        id="vibrations",
    ),
    pytest.param(
        3,
        "neb",
        RESTART_TESTS_FOLDER / "NEB" / "FIRE" / "zundel",
        id="neb",
    ),
    pytest.param(
        4,
        "instanton",
        RESTART_TESTS_FOLDER / "INSTANTON" / "070K_D",
        id="instanton",
    ),
    pytest.param(
        5,
        "driven_dynamics",
        RESTART_TESTS_FOLDER / "EDANVE",
        id="driven_dynamics",
    ),
]


class RestartRunner(Runner_regression):
    def __init__(self, *args, request_exit=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.request_exit = request_exit

    def create_client_list(self, driver_info, nid, test_settings):
        clients = super().create_client_list(driver_info, nid, test_settings)
        if self.request_exit:
            (self.tmp_dir / "EXIT").touch()
        return clients


def _run(runner, test_input, nid):
    error_msg = runner.run(test_input, nid)
    if error_msg is not None:
        raise RuntimeError(error_msg)


@pytest.mark.parametrize(("case_id", "motion_mode", "test_input"), RESTART_CASES)
def test_restart_matches_regression_reference(case_id, motion_mode, test_input):
    nid = 1000 + case_id
    interrupted = RestartRunner(
        Path("."),
        request_exit=True,
        check_numpy_output=False,
        check_xyz_output=False,
    )
    _run(interrupted, test_input, nid)
    assert (interrupted.tmp_dir / "RESTART").exists()
    exit_file = interrupted.tmp_dir / "EXIT"
    if exit_file.exists():
        exit_file.unlink()

    restarted = RestartRunner(Path("."), call_ipi="i-pi RESTART")
    _run(restarted, interrupted.tmp_dir, nid)
