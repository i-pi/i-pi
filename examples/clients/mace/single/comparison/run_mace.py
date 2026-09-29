#!/usr/bin/env python3
"""Serve ASE's MACECalculator through the i-PI UNIX socket."""

import json

from ase.calculators.socketio import SocketClient
from ase.io import read
from mace.calculators import MACECalculator


atoms = read("init.xyz")
with open("mace_kwargs.json", encoding="utf-8") as handle:
    mace_kwargs = json.load(handle)

atoms.calc = MACECalculator(
    model_paths="../mace.model",
    device="cpu",
    **mace_kwargs,
)
SocketClient(unixsocket="mace-comparison-calculator").run(
    atoms,
    use_stress=True,
)
