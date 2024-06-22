from ase.io import read
from mace.calculators import MACECalculator
from ase.calculators.socketio import SocketClient

# Define atoms object
atoms = read("init.pdb", 0)

# Set ASE calculator #################
calcs = []
calc = MACECalculator(
    "../2_training/TePIGS_model.model", device="cuda", default_dtype="float32"
)
atoms.set_calculator(calc)

# Create Client
host = "driver-pigs"
client = SocketClient(unixsocket=host)
client.run(atoms, use_stress=True)
