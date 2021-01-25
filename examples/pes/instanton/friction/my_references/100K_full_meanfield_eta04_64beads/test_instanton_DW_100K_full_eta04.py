from sysbath.engine.simulation import Simulation
from sysbath.utils.units import Constants, Convert
from sysbath.i_o import init_from_file, Position_report, NM_report
import pytest
import numpy as np
import sys
from sysbath.utils import Report
from pathlib import Path

test = False

folder = Path(__file__).parent
pos_file_name = folder / "init_pos.dat"
folder_pos = Path(__file__).parent
nfreq = 1

temperature = 100 * Convert.K2au
sys_pes = "Double_well"
optimizer = "Hessian_optimizator"
mapper = "Instanton_mapper"
asr = "tete"
inst_mode = "full"
update_hessian = "recompute"
tol_gradient = 1e-8
nsteps = 100
nbeads = 64

bath_type = "Ohmic"
omega_c = 500 * Convert.invcm2au
eta = 4.185819  # eta/mw_b ~1
eta *= 0.4
epsilon = +0.00
delta = 0.00
deltaQ = 0.5


def test_instanton_DW_100K_full_SD_meanfield_04():
    RP = init_from_file(pos_file_name=pos_file_name)
    RP = RP.RPC(nbeads)
    RP.set_temperature(temperature)
    RP.set_SD_mean_field(bath_type, eta, omega_c, epsilon, delta, deltaQ)
    simul = Simulation(RP)

    simul.init_calculator(sys_pes)
    simul.RP.set_hessian(dx=0.001)

    simul.init_mapper(
        mapper, update_hessian=update_hessian, asr=asr, inst_mode=inst_mode
    )
    simul.init_optimizator(optimizer, tol_gradient=tol_gradient)

    pos_rep = Position_report(simul=simul, folder=folder_pos, nfreq=nfreq)
    simul.add_report(pos_rep)

    simul.run(nsteps)

    np.savetxt(folder / "pot.dat", RP.potentials * Convert.au2eV)
    np.savetxt(folder / "hessian.dat", simul.RP.hessian)
    if test:

        RP1 = init_from_file(pos_file_name=folder / "positions.dat")
        RP2 = init_from_file(pos_file_name=folder / "reference.dat")
        np.testing.assert_allclose(RP1.get_positions(), RP2.get_positions())

        pot_test = np.loadtxt(folder / "pot.dat")
        pot_ref = np.loadtxt(folder / "reference_potential.dat")
        np.testing.assert_allclose(pot_test, pot_ref)

        pot_test = np.loadtxt(folder / "hessian.dat")
        pot_ref = np.loadtxt(folder / "reference_hessian.dat")
        np.testing.assert_allclose(pot_test, pot_ref)


if __name__ == "__main__":

    test_instanton_DW_100K_full_SD_meanfield_04()
