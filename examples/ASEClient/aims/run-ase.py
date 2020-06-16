from ase.calculators.aims import Aims

# from ipi.interfaces.clients import ClientASE
from ase.calculators.socketio import SocketClient
from ase.io import read

# Define geometry to initialize AIMS   ################

atoms = read("molecule.in", 0, "aims")

# Define settings for AIMS calculator  ################
usr_settings = {
    "aims_command": "<AIMS-COMMAND>",
    "species_dir": "<SPECIES-DIR>",
    "outfilename": "aims.out",
}

dft_settings = {
    # system settings
    "xc": "pbe",
    "spin": "none",
    "compute_forces": True,
}

# Set aims calculator ######################
workdir = "aims_rundir"
aux_settings = {"label": workdir}
calc = Aims(**usr_settings, **dft_settings, **aux_settings)

atoms.set_calculator(calc)

# Create Client ############################
# inet
port = 10200
host = "localhost"
client = SocketClient(host=host, port=port)

# 'unix'
# client = SocketClient(unixsocket=host)

client.run(atoms)
