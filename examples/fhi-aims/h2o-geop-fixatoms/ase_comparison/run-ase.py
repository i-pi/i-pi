#! /usr/bin/env python2
""" This script performs geometry optimization for 1 water molecule 
    using LAMMPS and fixing some atoms to compare with the new 
    'fixatoms' feature in the i-PI optimization
"""
import io, os, sys
from socket import gethostname
from warnings import warn
import numpy as np
from ase import Atoms, Atom
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator
from ase.io import read, write
from ase.optimize import BFGS
from ase.io.trajectory import Trajectory
from ase.constraints import *

# Environment-dependent parameters -- please configure according to machine
# Note that FHI-aims support for the i-PI socket needs smth like
# `make -f Makefile.ipi ipi.mpi`
species_dir = '/home/fidanyan/soft/fhi-aims.190214/species_defaults/light'
command = 'mpirun -np 4 ipi.aims.190214.mpi.x'
hostname = gethostname()
port = 31415

def Optimize_And_Print(atoms=None, energies=None, n=None):
    """ Performs geometry optimization and outputs the geometries.
        Also appends energy of this displacement
    """
    if atoms == None or energies == None or n == None:
        raise RuntimeError("All arguments of 'Optimize_And_Print()' are mandatory.")

    print(n)
    opt = BFGS(atoms, trajectory='opt.traj', logfile='opt.log')
    opt.run(fmax=0.01, steps=60)
    energies.append(atoms.get_potential_energy(force_consistent=False, apply_constraint=True))

    #reading in the trajectory file created during optimization
    traj = Trajectory("opt.traj", 'r')
    nsteps = opt.get_number_of_steps()

    outFileName_geop = "relaxation_%i.xyz" % n
    if os.path.exists(outFileName_geop):
        os.remove(outFileName_geop)
    #write each structure from the .traj file in .xyz format
    if nsteps == 0:
        write(outFileName_geop, atoms, format='xyz', append=True)
    else:
        for i in range(-nsteps, 0):
            atoms2 = traj[i]
            write(outFileName_geop, atoms2, format='xyz', append=True)

    atoms2 = traj[-1]
    write(outFileName, atoms2, format='xyz', append=True)


# === HERE THE ALGORITHM STARTS ===
energies = []
constraints = []  # the list for all constraints applied to a system

atoms = read('geometry.initial.in', format='aims')
initial_pos = atoms.positions.copy()

global outFileName
outFileName = 'relaxed.xyz'
#if os.path.exists(outFileName): os.remove(outFileName)

aims = Aims(command=command,
            species_dir=species_dir,
            use_pimd_wrapper=(hostname, port),
            compute_forces=True,
            xc='PBE',
            spin='none',
#            relativistic='atomic_zora scalar',
            k_grid=[1,1,1],
            # trying to reach convergence
            sc_iter_limit='200',
            # re-initialize Pulay mixer sometimes:
            sc_init_iter='200',
            # small mixing needed for metal
            charge_mix_param='0.02',
            # big blur also needed for metal
            occupation_type='gaussian 0.1',
            sc_accuracy_forces='1E-3',
#            use_dipole_correction=True,
#            vdw_correction_hirshfeld=True,
#            vdw_pair_ignore='Rh Rh',
            output_level='MD_light'
           )

calc = SocketIOCalculator(aims, log='log.ase', port=port)
atoms.calc = calc

# 0 fixed atoms
n = 0
constraints=[]
atoms.set_constraint(constraints)
Optimize_And_Print(atoms, energies, n)

# 1 fixed atom
n = 1
atoms.set_positions(initial_pos, apply_constraint=False)
constraints=FixAtoms(indices=[1])
atoms.set_constraint(constraints)
Optimize_And_Print(atoms, energies, n)

# 2 fixed atoms
n = 2
atoms.set_positions(initial_pos, apply_constraint=False)
constraints=FixAtoms(indices=[1,2])
atoms.set_constraint(constraints)
Optimize_And_Print(atoms, energies, n)

# 3 fixed atoms
n = 3
atoms.set_positions(initial_pos, apply_constraint=False)
constraints=FixAtoms(indices=[0,1,2])
atoms.set_constraint(constraints)
Optimize_And_Print(atoms, energies, n)

