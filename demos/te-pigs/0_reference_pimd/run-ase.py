from ase.io import read
from mace.calculators import mace_off
from ase.calculators.socketio import SocketClient

# Define atoms object
atoms = read("init.pdb", 0)

# Set ASE calculator #################
calcs = []
calc = mace_off(model="small", device="cuda", default_dtype="float32")
atoms.set_calculator(calc)

# Create Client
host = "driver"
client = SocketClient(unixsocket=host)
client.run(atoms, use_stress=True)
