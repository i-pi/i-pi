import numpy as np
import sys

import argparse
from ipi.utils.messages import verbosity, info

""" Reads all the information needed from a i-pi RESTART file and compute the partition functions of the reactant, transition state (TS) or
instanton according to J. Phys. Chem. Lett. 7, 437(2016) (Instanton Rate calculations) or J. Chem. Phys. 134, 054109 (2011) (Tunneling Splitting)


Syntax:    python  Instanton_postproc.py  <checkpoint_file> -c <case> -t  <temperature (K)>  (-n <nbeads(full polymer)>) (-freq <freq_reactant.dat>)

Examples for rate calculation:
           python  Instanton_postproc.py   RESTART  -c  instanton    -t   300
           python  Instanton_postproc.py   RESTART  -c  reactant     -t   300            -n 50
           python  Instanton_postproc.py   RESTART  -c    TS         -t   300

Examples for splitting  calculation (2 steps):
         i)   python  Instanton_postproc.py   RESTART  -c  reactant   -t   10  -n 32 --->this generate the 'eigenvalues_reactant.dat' file
         ii)  python  Instanton_postproc.py   RESTART  -c  instanton  -t   10  -freq eigenvalues_reactant.dat



Type python Instanton_postproc.py -h for more information


Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.
"""

# Y. Litman, 2017.

# You can insert the i-pi path with the following lines.
# Uncomment them and adjust the ipi_path variable

# ipi_path='/home/litman/Yair/Instanton/I-PI-mc/i-pi-mc'

# if not (os.path.exists(ipi_path)):
#    print 'We can not find ipi in %s' %ipi_path
#    print 'Please correct the path'
#    sys.exit()
# sys.path.insert(0, ipi_path)

from ipi.engine.simulation import Simulation
from ipi.utils.units import unit_to_internal, Constants
from ipi.utils.instools import red2comp
from ipi.utils.hesstools import clean_hessian
from ipi.engine.motion.instanton import SpringMapper

# UNITS
K2au = unit_to_internal("temperature", "kelvin", 1.0)
kb = Constants.kb
hbar = Constants.hbar
eV2au = unit_to_internal("energy", "electronvolt", 1.0)
cal2au = unit_to_internal("energy", "cal/mol", 1.0)
cm2au = unit_to_internal("frequency", "hertz", 1.0) * 3e10

# INPUT
parser = argparse.ArgumentParser(
    description="""Post-processing routine in order to obtain different quantities from an instanton (or instanton related) calculation. These quantities can be used for the calculation of rates or tunneling splittings in the instanton approximation."""
)
parser.add_argument("input", help="Restart file")
parser.add_argument(
    "-c",
    "--case",
    default=False,
    help="Type of the calculation to analyse. Options: 'instanton', 'reactant' or 'TS'.",
)
parser.add_argument(
    "-t", "--temperature", type=float, default=0.0, help="Temperature in K."
)
parser.add_argument(
    "-asr",
    "--asr",
    default="poly",
    help="Removes the zero frequency vibrational modes depending on the symmerty of the system",
)
parser.add_argument(
    "-e", "--energy_shift", type=float, default=0.0, help="Zero of energy in eV"
)
parser.add_argument(
    "-f",
    "--filter",
    default=[],
    help="List of atoms indexes to filter (i.e. eliminate its componentes in the position,mass and hessian arrays. It is 0 based.",
    type=int,
    action="append",
)
parser.add_argument(
    "-n",
    "--nbeadsR",
    default=0,
    help="Number of beads (full polymer) to compute the approximate partition function (only reactant case)",
    type=int,
)
parser.add_argument(
    "-freq",
    "--freq_reac",
    default=None,
    help="List of frequencies of the minimum. Required for splitting calculation.",
)
parser.add_argument(
    "-q",
    "--quiet",
    default=False,
    action="store_true",
    help="Avoid the Qvib and Qrot calculation in the instanton case.",
)

args = parser.parse_args()
inputt = args.input
case = args.case
temp = args.temperature * K2au
asr = args.asr
V00 = args.energy_shift
filt = args.filter
nbeadsR = args.nbeadsR
input_freq = args.freq_reac
quiet = args.quiet
Verbosity = verbosity
Verbosity.level = "quiet"

if case not in list(["reactant", "TS", "instanton"]):
    raise ValueError(
        "We can not indentify the case. The valid cases are: 'reactant', 'TS' and 'instanton'"
    )
if asr not in list(["poly", "linear", "crystal", "none"]):
    raise ValueError(
        "We can not indentify asr case. The valid cases are: 'poly', 'crystal' , 'linear' and 'none'"
    )

if asr == "poly":
    nzeros = 6
elif asr == "crystal":
    nzeros = 3
else:
    nzeros = 0
if asr == "linear":
    raise NotImplementedError("Sum rules for linear molecules is not implemented")

if args.temperature == 0.0:
    raise ValueError("The temperature must be specified.'")


# ----- Define some functions-----------------


def get_double(q0, nbeads0, natoms, h0):
    """Takes nbeads, positions and hessian (only the 'physcal part') of the half polymer and
    returns the equivalent for the full ringpolymer."""
    q = np.concatenate((q0, np.flipud(q0)), axis=0)
    nbeads = 2 * nbeads0
    ii = 3 * natoms
    iii = 3 * natoms * nbeads0

    h = np.zeros((iii * 2, iii * 2))
    h[0:iii, 0:iii] = h0

    # diagonal block
    for i in range(nbeads0):
        x = i * ii + iii
        y = ((nbeads0 - 1) - i) * ii
        h[x : x + ii, x : x + ii] = h0[y : y + ii, y : y + ii]

    return q, nbeads, h


def spring_pot(nbeads, q, omega2, m3):
    e = 0.0
    for i in range(nbeads - 1):
        dq = q[i + 1, :] - q[i, :]
        e += omega2 * 0.5 * np.dot(m3[0] * dq, dq)
    return e


def Filter(pos, h, natoms, m, m3, filt):
    filt3 = []
    for i in filt:
        filt3.append(3 * i)
        filt3.append(3 * i + 1)
        filt3.append(3 * i + 2)
    pos = np.delete(pos, filt3, axis=1)
    aux = np.delete(h, filt3, axis=1)
    h = np.delete(aux, filt3, axis=0)
    m = np.delete(m, filt, axis=0)
    m3 = np.delete(m3, filt3, axis=1)
    natoms = natoms - len(filt)
    return pos, h, natoms, m, m3


def get_rp_freq(w0, nbeads, temp, mode="rate"):
    """Compute the ring polymer frequencies for multidimensional harmonic potential
    defined by the frequencies w0."""
    hbar = 1.0
    kb = 1
    betaP = 1 / (kb * nbeads * temp)
    factor = betaP * hbar
    w = 0.0
    ww = []

    if np.amin(w0) < 0.0:
        print("@get_rp_freq: We have a negative frequency, something is going wrong.")
        sys.exit()

    if mode == "rate":
        for n in range(w0.size):
            for k in range(nbeads):
                if w0[n] == 0 and k == 0:
                    continue
                w += np.log(
                    factor
                    * np.sqrt(
                        4.0
                        / (betaP * hbar) ** 2
                        * np.sin(np.absolute(k) * np.pi / nbeads) ** 2
                        + w0[n]
                    )
                )
                # note the w0 is the eigenvalue ( the square of the frequency )
        return w

    elif mode == "splitting":
        for n in range(w0.size):
            for k in range(nbeads):
                # note the w0 is the eigenvalue ( the square of the frequency )
                ww = np.append(
                    ww,
                    np.sqrt(
                        4.0
                        / (betaP * hbar) ** 2
                        * np.sin((k + 1) * np.pi / (2 * nbeads + 2)) ** 2
                        + w0[n]
                    ),
                )
        return np.array(ww)
    else:
        print("We can't indentify the mode")
        sys.exit()


def save_frequencies(d, nzeros, filename="freq.dat"):
    """Small function to save the frequencies in a file
    d: array with the eigenvalues of the (extended) hessian"""

    outfile = open(filename, "w")
    aux = np.zeros(nzeros)
    freq = np.sign(d) * np.absolute(d) ** 0.5 / cm2au
    dd = np.concatenate((aux, freq))
    np.savetxt(outfile, dd.reshape(dd.size, 1), header="Frequencies (cm^-1)")
    outfile.close()
    print("We saved the frequencies in {}".format(filename))


# -----END of some functions-----------------

# -----READ---------------------------------
print("\nWe are ready to start. Reading {} ... (This can take a while)".format(inputt))

simulation = Simulation.load_from_xml(
    open(inputt), custom_verbosity="quiet", request_banner=False, read_only=True
)


beads = simulation.syslist[0].motion.beads.clone()
m = simulation.syslist[0].motion.beads.m.copy()
nbeads = simulation.syslist[0].motion.beads.nbeads
natoms = simulation.syslist[0].motion.beads.natoms

if case == "reactant":
    if nbeadsR == 0:
        print(
            "We have to specify number of beads for computing the partition function in the reactant case"
        )
        sys.exit()

if case != "instanton" and nbeads > 1:
    print(("Incompatibility between case and nbeads in {}.".format(inputt)))
    print(("case {} , beads {}".format(case, nbeads)))
    sys.exit()


# Depending the case we read from the restart file different things:
if case == "reactant":
    dynmat = simulation.syslist[0].motion.dynmatrix.copy()

    ism = beads.m3 ** (0.5)
    ismm = np.outer(ism, ism)
    h = np.multiply(dynmat, ismm)
    pos = beads.q
    m3 = beads.m3
    if len(filt) > 0:
        pos, h, natoms, m, m3 = Filter(pos, h, natoms, m, beads.m3, filt)

elif case == "TS":
    pos = beads.q
    h = simulation.syslist[0].motion.optarrays["hessian"].copy()
    m3 = beads.m3
    pots = simulation.syslist[0].motion.optarrays["old_u"]
    V0 = simulation.syslist[0].motion.optarrays["energy_shift"]

    if V00 != 0.0:
        print("Overwriting energy shift with the provided values")
        V0 = V00 * eV2au
elif case == "instanton":
    hessian = simulation.syslist[0].motion.optarrays["hessian"].copy()
    mode = simulation.syslist[0].motion.options["mode"]
    temp2 = simulation.syslist[0].ensemble.temp
    pots = simulation.syslist[0].motion.optarrays["old_u"]
    grads = -simulation.syslist[0].motion.optarrays["old_f"]
    V0 = simulation.syslist[0].motion.optarrays["energy_shift"]

    if V00 != 0.0:
        print("Overwriting energy shift with the provided values")
        V0 = V00 * eV2au

    if np.absolute(temp - temp2) / K2au > 2:
        print(
            "\n Mismatch between provided temperature and temperature in the calculation"
        )
        sys.exit()

    if mode == "rate":
        h0 = red2comp(hessian, nbeads, natoms)
        pos, nbeads, hessian2 = get_double(beads.q, nbeads, natoms, h0)
        hessian = hessian2
        m3 = np.concatenate((beads.m3, beads.m3), axis=0)
        omega2 = (temp * nbeads * kb / hbar) ** 2
        if not quiet:
            spring = SpringMapper.spring_hessian(
                natoms, nbeads, beads.m3[0], omega2, mode="full"
            )
            h = np.add(hessian, spring)
    elif mode == "splitting":
        if input_freq is None:
            print(
                'Please provide a name of the file containing the list of the frequencies for the minimum using "-freq" flag'
            )
            print(" You can generate that file using this script in the case reactant.")
            sys.exit()

        print(("Our linear polymer has  {}".format(nbeads)))
        pos = beads.q
        m3 = beads.m3
        omega2 = (temp * nbeads * kb / hbar) ** 2

        if not quiet:
            h0 = red2comp(hessian, nbeads, natoms)
            spring = SpringMapper.spring_hessian(
                natoms, nbeads, beads.m3[0], omega2, mode="splitting"
            )
            h = np.add(h0, spring)
            if asr != "none":
                print(
                    "We are changing asr to none since we consider a fixed ended linear polimer for the post-processing"
                )
                asr = "none"
    else:
        print("We can not recognize the mode. STOP HERE")
        sys.exit()

# ----------------------------------------------------------START----------------------------------------------
beta = 1.0 / (kb * temp)
betaP = 1.0 / (kb * (nbeads) * temp)

print(("\nTemperature: {} K".format(temp / K2au)))
print(("NBEADS: {}".format(nbeads)))
print(("atoms:  {}".format(natoms)))
print(("ASR:    {}".format(asr)))
print(("1/(betaP*hbar) = {:8.5f}".format((1 / (betaP * hbar)))))

if not quiet or case == "reactant" or case == "TS":
    print("Diagonalization ... \n\n")
    d, w, detI = clean_hessian(h, pos, natoms, nbeads, m, m3, asr, mofi=True)
    print("Lowest 10 frequencies (cm^-1)")
    d10 = np.array2string(
        np.sign(d[0:10]) * np.absolute(d[0:10]) ** 0.5 / cm2au,
        precision=2,
        max_line_width=100,
        formatter={"float_kind": lambda x: "%.2f" % x},
    )
    print(("{}".format(d10)))

    save_frequencies(d, nzeros)

if case == "reactant":
    Qtras = ((np.sum(m)) / (2 * np.pi * beta * hbar**2)) ** 1.5

    if asr == "poly":
        Qrot = (8 * np.pi * detI / ((hbar) ** 6 * (beta) ** 3)) ** 0.5
    else:
        Qrot = 1.0

    # logQvib    = -np.sum( np.log( 2*np.sinh( (beta*hbar*np.sqrt(d)/2.0) )  ))   #Limit n->inf
    logQvib_rp = -get_rp_freq(d, nbeadsR, temp)

    # This file is created to use afterwards in the splitting calculations
    outfile = open("eigenvalues_reactant.dat", "w")
    aux = np.zeros(nzeros)
    dd = np.concatenate((aux, d))
    np.savetxt(outfile, dd.reshape(1, dd.size))
    outfile.close()

    print(("\nWe are done. Reactants. Nbeads {}".format(nbeadsR)))
    print(("{:14s} | {:8s} | {:8s}".format("Qtras(bohr^-3)", "Qrot", "logQvib_rp")))
    print(("{:14.3f} | {:8.3f} |{:8.3f}\n".format(Qtras, Qrot, logQvib_rp)))
    print("A file with the eigenvalues in atomic units was generated\n")

elif case == "TS":
    Qtras = ((np.sum(m)) / (2 * np.pi * beta * hbar**2)) ** 1.5

    if asr == "poly":
        Qrot = (8 * np.pi * detI / ((hbar) ** 6 * (beta) ** 3)) ** 0.5
    else:
        Qrot = 1.0

    logQvib = -np.sum(
        np.log(2 * np.sinh((beta * hbar * np.sqrt(np.delete(d, 0)) / 2.0)))
    )

    U = pots.sum() - V0

    print("\nWe are done. TS")
    print(("Partition functions at {} K".format(temp / K2au)))
    print(("\nQtras: {}".format(Qtras)))
    print(("Qrot: {}".format(Qrot)))
    print(("logQvib: {}".format(logQvib)))
    print(
        (
            "Potential energy at TS:  {} eV, V/kBT {}\n".format(
                U / eV2au, U / (kb * temp)
            )
        )
    )

elif case == "instanton":
    if mode == "rate":
        Qtras = ((np.sum(m)) / (2 * np.pi * beta * hbar**2)) ** 1.5

        if asr == "poly" and not quiet:
            Qrot = (8 * np.pi * detI / ((hbar) ** 6 * (betaP) ** 3)) ** 0.5
            Qrot /= nbeads**3
        else:
            Qrot = 1.0

        if not quiet:
            del_freq = np.sign(d[1]) * np.absolute(d[1]) ** 0.5 / cm2au
            print("Deleted frequency: {:8.3f} cm^-1".format(del_freq))

            if asr != "poly":
                print("WARNING asr != poly")
                print("First 10 eigenvalues")
                ten_eigv = np.sign(d[0:10]) * np.absolute(d[0:10]) ** 0.5 / cm2au
                print("{}".format(ten_eigv))
                print(
                    "Please check that this you don't have any unwanted zero frequency"
                )
            logQvib = (
                -np.sum(np.log(betaP * hbar * np.sqrt(np.absolute(np.delete(d, 1)))))
                + nzeros * np.log(nbeads)
                + np.log(nbeads)
            )
        else:
            logQvib = 0.0

        BN = 2 * np.sum(beads.m3[1:, :] * (beads.q[1:, :] - beads.q[:-1, :]) ** 2)
        factor = 1.0000  # default
        action1 = (2 * pots.sum() * factor - nbeads * V0) * 1.0 / (temp * nbeads * kb)
        action2 = spring_pot(nbeads, pos, omega2, m3) / (temp * nbeads * kb)

        print(
            "\nWe are done. Instanton rate. Nbeads {} (diff only {})".format(
                nbeads, nbeads / 2
            )
        )
        print(
            "   {:8s} {:8s}  | {:11s} | {:11s} | {:11s} | {:8s} ( {:8s},{:8s} ) |".format(
                "BN",
                "(BN*N)",
                "Qt(bohr^-3)",
                "Qrot",
                "log(Qvib*N)",
                "S/hbar",
                "S1/hbar",
                "S2/hbar",
            )
        )
        print(
            "{:8.3f} ( {:8.3f} ) | {:11.3f} | {:11.3f} | {:11.3f} | {:8.3f} ( {:8.3f} {:8.3f} ) |".format(
                BN,
                BN * nbeads,
                Qtras,
                Qrot,
                logQvib,
                (action1 + action2),
                action1,
                action2,
            )
        )
        print("\n\n")

    elif mode == "splitting":
        out = open(input_freq, "r")
        d_min = np.zeros(natoms * 3)
        aux = out.readline().split()
        if len(aux) != (natoms * 3):
            print(("We are expecting {} frequencies.".format((natoms * 3 - 6))))
            print(("instead we have read  {}".format(len(aux))))
        for i in range((natoms * 3)):
            d_min[i] = float(aux[i])
        d_min = d_min.reshape((natoms * 3))
        out.close()
        ww = get_rp_freq(np.sign(d_min) * d_min**2, nbeads, temp, mode="splitting")
        react = np.sum(np.log(ww))

        action1 = (pots.sum() - nbeads * V0) * 1 / (temp * nbeads * kb)
        action2 = spring_pot(nbeads, pos, omega2, m3) / (temp * nbeads * kb)
        action = action1 + action2
        if action / hbar > 5.0:
            print(
                "WARNING, S/h seems to big. Probably a proper energy shift is missing."
            )

        BN = np.sum(beads.m3[1:, :] * (beads.q[1:, :] - beads.q[:-1, :]) ** 2)

        if not quiet:
            inst = np.sum(np.log(np.sqrt(np.absolute(np.delete(d, [1])))))
            phi = np.exp(inst - react)
        else:
            phi = 1

        tetaphi = (
            betaP * hbar * np.sqrt(action / (2 * hbar * np.pi)) * np.exp(-action / hbar)
        )
        teta = tetaphi / phi
        h = -teta / betaP
        # cm2au= (2 * np.pi * 3e10 * 2.4188843e-17)  # Handy for debugging

        print("\n\nWe are done")
        print("Nbeads {}, betaP {} a.u.,hbar {} a.u".format(nbeads, betaP, hbar))
        print("")
        print("V0  {} eV ( {} Kcal/mol) ".format(V0 / eV2au, V0 / cal2au / 1000))
        print(
            "S1/hbar {} ,S2/hbar {} ,S/hbar {}".format(
                action1 / hbar, action2 / hbar, action / hbar
            )
        )
        print("BN {} a.u.".format(BN))
        print(
            "BN/(hbar^2 * betaN)  {}  (should be same as S/hbar) ".format(
                (BN / ((hbar**2) * betaP))
            )
        )
        print("")
        if quiet:
            print("phi is not computed because you specified the quiet option")
            print(
                ("We can provied only Tetaphi which value is {} a.u. ".format(tetaphi))
            )
        else:
            print(("phi {} a.u.   Teta {} a.u. ".format(phi, tetaphi / phi)))
            print(
                "Tunnelling splitting matrix element (h)  {} a.u ({} cm^-1)".format(
                    h, h / cm2au
                )
            )
    else:
        print("We can not recongnize the mode.")
        sys.exit()

info("\n\n", Verbosity.medium)
info(
    "Remember that the output obtained from this script simply gives you components\
      that you can use in order to calculate a rate or a tunneling splitting in the \
      instanton approximation.",
    verbosity.medium,
)
info(
    "Use, for example, the references below in order to obtain final desired results.",
    verbosity.medium,
)
info("Instanton Rate: J. Phys. Chem. Lett.  7, 4374(2016)", verbosity.medium)
info("Tunneling Splitting: J. Chem. Phys. 134, 054109 (2011)", verbosity.medium)
info("\n\n", verbosity.medium)
sys.exit(0)
