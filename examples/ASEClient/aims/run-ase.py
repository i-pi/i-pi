import os
import sys

import ase.io
from ase.calculators.aims import Aims
#from ipi.interfaces.clients import ClientASE
from ase.calculators.socketio import SocketClient
from ase.io import read 

####### Define geometry to initialize AIMS   ################

atoms = read("molecule.in", 0, "aims")

####### Define settings for AIMS calculator  ################
usr_settings = {
    "aims_command":"/home/litman/Yair/Installation/OpenMpi/openmpi-4.0.1/orte/tools/orterun/orterun -np 4  /home/litman/Codes/bin/aims.x",
    "species_dir":"/home/litman/Codes/FHIaims/species_defaults/light",
    "outfilename": "aims.out",
}

dft_settings = {
# system settings
    "xc": "pbe",
    "spin": "none",
    "compute_forces": True,
}

################## Set aims calculator ######################
workdir = "aims_rundir"
aux_settings = {"label": workdir}
calc = Aims(**usr_settings, **dft_settings, **aux_settings)

atoms.set_calculator(calc)

################## Create Client ############################
#inet
port=10200
host='localhost'
client = SocketClient(host=host, port=port)

#'unix'
#client = SocketClient(unixsocket=host)

client.run(atoms)


