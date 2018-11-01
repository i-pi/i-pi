!name of the folder where the DFT orbitals will be stored inside orbitals/ . If nothing is provided it will be used the default name "qespresso"
ORBITALS_FOLDER="simple"
!want to overwrite a previous orbitals' folder? If so set equal to T, otherwise it will be assumed false
OVERWRITE=T
!specify the path to the folder where the pseudo potential that is to be used is stored
PSEUDO_DIR="/storage/Codes/espresso-5.1.2/pseudo"
!pseudo potential to be used
PSEUDO="H.coulomb-ae.UPF"
!Energy cutoff for the plane waves basis set [Rydberg]
ECUTWFC=3.0
!type of K-points to use (gamma if no TABC are used, automatic otherwise)
K_POINTS="automatic"
!if K_points=="gamma" then the K-points' grid must be specified (example: 3 3 3) - offset is automatically set to zero
K_GRID=3 3 3
!Input wave function to be set for using the generated orbitals
WF="dati_funzione_onda.d"
