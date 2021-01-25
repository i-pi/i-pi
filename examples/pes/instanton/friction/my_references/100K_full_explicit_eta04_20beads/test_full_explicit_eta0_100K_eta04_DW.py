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

temperature = 100.00 * Convert.K2au
sys_pes = "Double_well"
optimizer = "Hessian_optimizator"
mapper = "Instanton_mapper"
asr = "tete"
inst_mode = "full"
update_hessian = "recompute"
tol_gradient = 3e-8
nsteps = 30
nbeads = 20
eta = 4.185819  # eta/mw_b ~1
eta *= 0.4
epsilon = 0.00
delta = 0.00
deltaQ = 1.00


def test_full_explicit_eta04_100K_DW():
    RP = init_from_file(pos_file_name=pos_file_name)
    RP = RP.RPC(nbeads)
    RP.set_temperature(temperature)
    simul = Simulation(RP)

    simul.init_calculator(sys_pes, eta=eta, epsilon=epsilon, delta=delta, deltaQ=deltaQ)
    simul.RP.set_hessian(dx=0.01)

    simul.init_mapper(
        mapper, update_hessian=update_hessian, asr=asr, inst_mode=inst_mode
    )
    simul.init_optimizator(
        optimizer, tol_gradient=tol_gradient, big_step=0.10 * np.sqrt(1837)
    )

    pos_rep = Position_report(simul=simul, folder=folder_pos, nfreq=nfreq)
    simul.add_report(pos_rep)

    simul.run(nsteps)

    # np.savetxt(folder/ 'SD.dat',np.asarray(SD))
    np.savetxt(folder / "pot.dat", RP.potentials * Convert.au2eV)
    np.savetxt(folder / "hessian.dat", simul.RP.hessian)

    with open(folder / "pot_sys.dat", "w") as f:
        sys_pos = RP.get_positions()
        # sys_pot = simul.RP.calculator.pes[0].v_sys(sys_pos )
        for p in sys_pos:
            f.write("{}\n".format(simul.RP.calculator.pes[0].v_sys(p) * Convert.au2eV))

    if test:

        RP1 = init_from_file(pos_file_name=folder / "positions.dat")
        RP2 = init_from_file(pos_file_name=folder / "reference.dat")
        np.testing.assert_allclose(RP1.get_positions(), RP2.get_positions())

        pot_test = np.loadtxt(folder / "pot.dat")
        pot_ref = np.loadtxt(folder / "reference_potential.dat")
        np.testing.assert_allclose(pot_test, pot_ref)

        pot_sys_test = np.loadtxt(folder / "pot_sys.dat")
        pot_sys_ref = np.loadtxt(folder / "reference_sys_potential.dat")
        np.testing.assert_allclose(pot_test, pot_ref)

        pot_test = np.loadtxt(folder / "hessian.dat")
        pot_ref = np.loadtxt(folder / "reference_hessian.dat")
        np.testing.assert_allclose(pot_test, pot_ref)


if __name__ == "__main__":

    test_full_explicit_eta04_100K_DW()
