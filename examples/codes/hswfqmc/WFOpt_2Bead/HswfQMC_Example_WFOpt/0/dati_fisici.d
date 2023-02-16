&dati_fisici
r_s=1.33                            !Density in r_s=a/a_0 units
crystal_cell='datex'                 !Crystal cell that specifies the initial protonic positions: 'bcc', 'fcc', 'hcp', 'hcp_w', 'mhcpo', 'sc_', 'mol', 'dat', 'datex', 'grp__', 'quadr'
file_reticolo='reticolo/rp_now.d'    !Path to the file containing the protonic initial position, for the cases 'dat__' and 'datex'
flag_molecular=T                     !Molecule in the crystal cell (implemented only for 'hcp__' and 'mhcpo')
strecthing_cov_bond=1.               !Distance between H atoms in the H2 molecule, in 0.74 Angstrom units (the equilibrium distance)
N_cell_side=8                        !Number of replica of the crystal cell on each direction (x,y,z). In the cases 'dat__' and 'datex', then N_cell_side specifies the number of atoms
/

!Distances are in a_0 units, where a_0 is the Bohr radius
!In order to express them in Angstrom, multiply for 0.53
!Energies are in Rydberg units per atom (multiply for 2 in order to have them in Hartree units)
!'.datex' files contain non-normalized positions, and the first line contains the L(1:3) (i.e. the x,y,z lengths of the simulation box)
!'.dat' files contain normalized positions
