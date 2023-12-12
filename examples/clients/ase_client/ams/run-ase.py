#!/usr/bin/env amspython
from ase.calculators.socketio import SocketClient
from ase.io import read
from scm.plams import Settings, finish, init
from scm.plams.interfaces.adfsuite.ase_calculator import AMSCalculator

"""
Example illustrating how to use AMS as a client together with i-PI.

The i-PI configuration in input.xml sets up a short ring-polymer molecular
dynamics (RPMD) simulation. 

The connection between server and client is of type 'unixsocket' (a file in
/tmp). This example will only run on Linux.

The settings below set up AMS to use the UFF force field to calculate the
energy, forces, and stress tensor of whatever system the i-PI server requests
and return them to i-PI.

AMS is run in "AMSWorker" mode (interactive mode), which means that AMS does
**not** shut down and start up again between calculations.

This example can also be found in the AMS documentation. See
https://www.scm.com/doc/plams/examples/i-PI-AMS.html

Make sure to launch this script using the Python interpreter included with AMS:
$AMSBIN/amspython run-ase.py

"""


def main():
    init()

    use_stress = False
    atoms = read("firstframe.xyz")

    sett = Settings()
    sett.input.ams.Task = "SinglePoint"
    sett.input.ams.Properties.Gradients = "True"
    sett.input.ams.Properties.StressTensor = str(use_stress)
    sett.input.forcefield.type = "UFF"
    sett.runscript.nproc = 1

    with AMSCalculator(settings=sett, amsworker=True) as calc:
        atoms.set_calculator(calc)
        client = SocketClient(
            unixsocket="driver-irpmd-16"
        )  # socket should match the one given in input.xml
        client.run(atoms, use_stress=use_stress)

    finish()


if __name__ == "__main__":
    main()
