from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketClient, SocketIOCalculator
from ase.io import read

# ###### Define geometry to initialize AIMS   ################

atoms = read("molecule.in", 0, "aims")

# ###### Define settings for AIMS calculator  ################
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

# ################# Set aims calculator ######################
workdir = "aims_rundir"
port_aims = 12345
aux_settings = {"label": workdir, "use_pimd_wrapper": ("localhost", port_aims)}

calc = Aims(**usr_settings, **dft_settings, **aux_settings)

# atoms.set_calculator(calc)

# ################# Create Client ############################
# inet
port_ipi = 10200
host_ipi = "localhost"
client = SocketClient(host=host_ipi, port=port_ipi)

# ################# Create ASE SERVER ############################

with SocketIOCalculator(calc, log="socketio.log", port=port_aims) as io_calc:
    atoms.set_calculator(io_calc)
    client.run(atoms)
