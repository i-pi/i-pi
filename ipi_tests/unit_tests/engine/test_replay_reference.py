"""Replay must preserve each bead's coordinates, units, and frame ordering."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.engine.atoms import Atoms as IAtoms
from ipi.engine.cell import Cell
from ipi.utils.io import print_file
from ipi_tests.unit_tests.engine.inputs import harmonic_input


@pytest.mark.parametrize(
    "mode,separate",
    [("xyz", False), ("pdb", False), ("ase", False), ("xyz", True), ("pdb", True)],
)
def test_replay_reads_known_frames_in_bead_order(simulation, tmp_path, mode, separate):

    q = np.arange(2 * 3 * 6, dtype=float).reshape(2, 3, 6) / 10
    paths = (
        [tmp_path / f"trajectory_{b}.{mode}" for b in range(3)]
        if separate
        else [tmp_path / f"trajectory.{mode}"]
    )
    for i, path in enumerate(paths):
        with path.open("w") as f:
            for frame in range(2):
                for bead in ([i] if separate else range(3)):
                    atoms = IAtoms(2)
                    atoms._names[:] = ["H", "O"]
                    atoms._q[:] = q[frame, bead]
                    print_file(
                        mode,
                        atoms,
                        Cell(np.eye(3) * 10),
                        filedesc=f,
                        title=f"Step: {frame} Bead: {bead} ",
                        key="positions",
                    )
    root = ET.fromstring(harmonic_input(nbeads=3))
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "replay")
    ET.SubElement(motion, "file", mode=mode, units="atomic_unit").text = (
        f"trajectory_*.{mode}" if separate else str(paths[0])
    )
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        system = sim.syslist[0]
        system.motion.exit_on_convergence = False
        for frame in range(2):
            system.motion.step(frame)
            np.testing.assert_allclose(system.beads.q, q[frame], atol=2e-3)
        system.motion.step(2)
        assert system.motion.finished
        for stream in system.motion.rfile if separate else [system.motion.rfile]:
            stream.close()
