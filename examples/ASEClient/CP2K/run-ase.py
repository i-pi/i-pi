import os
import sys

from ase.calculators.socketio import SocketClient
from ase.io import read
from ase.calculators.cp2k import CP2K

# Define atoms object

atoms = read("init.xyz", 0)

# Set CP2K calculator #################

workdir = "CP2Ktest"
aux_settings = {"label": workdir}

if "ASE_CP2K_COMMAND" not in os.environ:
    print("$ASE_CP2K_COMMAND not defined")
    print('You can try with:\n "export ASE_CP2K_COMMAND=/usr/bin/cp2k_shell" \n\n')
    sys.exit()

calc = CP2K(**aux_settings)

atoms.set_calculator(calc)

# Create Client
# inet
port = 10200
host = "localhost"
client = SocketClient(host=host, port=port)

client.run(atoms)
