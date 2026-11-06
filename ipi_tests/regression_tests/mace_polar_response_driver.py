"""Small MACE-POLAR socket driver used by the EDA/FFDielectric comparison test.

The response tensors are deliberately supplied from a file.  This keeps the
comparison focused on the dynamics: both implementations receive the same
MACE energy/forces and the same Born charges.
"""

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

from ipi.pes._mace import MACE_driver


def _load_socket_driver():
    path = Path(__file__).resolve().parents[2] / "bin" / "i-pi-py_driver"
    spec = importlib.util.spec_from_file_location("ipi_py_driver", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run_driver


class MACEPolarResponseDriver(MACE_driver):
    """MACE-POLAR with fixed response tensors added to socket extras."""

    def __init__(self, bec_path, *args, **kwargs):
        self._bec = np.loadtxt(bec_path).reshape((-1, 3, 3))
        super().__init__(*args, **kwargs)

    def post_process(self, properties, structure):
        potential, forces, virial, extras = super().post_process(properties, structure)
        data = json.loads(extras) if extras else {}
        data.setdefault("dipole", [0.0, 0.0, 0.0])
        data["BEC"] = self._bec.tolist()
        data["piezoelectric"] = np.zeros((3, 3, 3)).tolist()
        return potential, forces, virial, json.dumps(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--address", required=True)
    parser.add_argument("--template", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--mace-kwargs", required=True)
    parser.add_argument("--bec", required=True)
    args = parser.parse_args()

    driver = MACEPolarResponseDriver(
        bec_path=args.bec,
        template=args.template,
        model=args.model,
        device="cpu",
        mace_kwargs=args.mace_kwargs,
    )
    _load_socket_driver()(unix=True, address=args.address, driver=driver)


if __name__ == "__main__":
    main()
