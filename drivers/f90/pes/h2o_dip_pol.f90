!
!=============================================================================
!
! File: h2o_dip_pol.f90
! Original author: Tomislav Begusic
! Date: 03 April 2022
! Description: Simple test program to call stand-alone
!              Implements water dipole and polarizability 
!              functions for qTIP4P/f model described in
!              T. Begusic, G. A. Blake, Nat. Commun. 14, 1950 (2023)
!
! NOTES:
! (1) ATOMIC UNITS ARE USED THROUGHOUT
! (2) WE ASSUME THAT THE ATOMS ARE IN A LIST ARRANGED WITH MOLECULES
!     AS [ {O,H,H},{O,H,H},..... ].
!
! Disclaimer: please thoroughly check the output of this before 
! using it for anything interesting!
!
!=============================================================================
!


!
!====================================================================================
!
! The routines below implement dipole, its derivative, and polarizability of water. 
!
! Expects atoms in the order O H H O H H O H H ... 
!====================================================================================
!

! Dipole moment.
SUBROUTINE h2o_dipole(box, nat, atoms, compute_der, dip, dip_der, pol)

  IMPLICIT NONE
  
  !====================================================
  ! Input arguments and output values.
  !====================================================
  DOUBLE PRECISION, INTENT(IN) :: box(3)              !Lattice constants for PBC.
  INTEGER, INTENT(IN) :: nat                          !Number of atoms.
  DOUBLE PRECISION, INTENT(IN) :: atoms(nat, 3)       !Atom coordinates.
  LOGICAL, INTENT(IN) :: compute_der                  !Switch for computing dipole gradients.
  DOUBLE PRECISION, INTENT(INOUT) :: dip(3)           !Output dipole.
  DOUBLE PRECISION, INTENT(INOUT) :: dip_der(nat, 3)  !Output dipole gradient.
  DOUBLE PRECISION, INTENT(INOUT) :: pol(3, 3)        !Output polarizability.
  !====================================================

  !====================================================
  ! Useful constants.
  !====================================================
  DOUBLE PRECISION, PARAMETER :: pi = DACOS(-1.0d0), twopi = 2 * pi, sqrtpi = SQRT(pi), sqrt2 = SQRT(2.0d0)
  DOUBLE COMPLEX, PARAMETER :: IU = CMPLX(0.0d0, 1.0d0, KIND=8)
  DOUBLE PRECISION, PARAMETER :: angtoau = 1.88973d0
  !====================================================

  !========================================================================================
  ! Parameters related to the TIP4P force fields, needed for computing the M site.
  !========================================================================================
  DOUBLE PRECISION, PARAMETER :: qm = -1.1128d0 ! Charge of M site.
  DOUBLE PRECISION, PARAMETER :: qh = - qm / 2  ! Charge of H atom.
  DOUBLE PRECISION, PARAMETER :: charges(3) = (/ qm, qh, qh /) !(n_charged_atom_per_mol)
  DOUBLE PRECISION, PARAMETER :: gam = 0.73612d0, gam2 = 0.5d0 * (1-gam) ! Parameter controlling M site position.
  !========================================================================================

  !========================================================================================
  ! Ewald parameter specifying the degree of convergence of real and reciprocal sums.
  ! Larger number implies longer and more accurate calculation.
  ! In other codes, this parameter is multiplied by pi, here pi is separately
  ! implemented in the equations (i.e., ewald_param = 1 corresponds to pi in other codes).
  !========================================================================================
  DOUBLE PRECISION, PARAMETER :: ewald_param = 1.0d0 ! Ewald sum parameter.
  !========================================================================================

  !================================================
  ! Polarizability-related parameters.
  !================================================
  !Coefficients for the flexible polarizability from G. Avila, JCP 122, 144310 (2005).
  !All units in angstrom, last constant is 1/ (4 pi eps_0) * 10^-10.
  !Axes defined as in the original paper, with a rotation x->y, y->z, z->x.
  DOUBLE PRECISION, PARAMETER :: r_eq = 0.95843d0, phi_eq = 104.44d0 / 180 * pi ! Distance in angstrom, angle in degrees.
  !Parameters from tables 6 and 7 of the JCP paper cited above. a_diag corresponds 
  !to the xx, yy, zz components, and a_off_diag is for the xy component.
  !a_index_diag and a_index_off_diag are the exponents for the s1, s2, and s3 variables (i, j, k in tables 6 and 7).
  INTEGER, PARAMETER :: ncoeff_diag = 22, ncoeff_off_diag = 13 ! Number of coefficients for the diagonal and off-diagonal terms.
  DOUBLE PRECISION, PARAMETER, DIMENSION(3, ncoeff_diag) :: a_diag = RESHAPE( (/ &
    1.64720d0, 1.58787d0, 1.53818d0, & !000
    0.24730d0,-0.07614d0, 0.08361d0, & !010
    2.39394d0, 1.60710d0, 0.82863d0, & !100
    0.89549d0, 0.71119d0,-0.04232d0, & !002
    0.03954d0, 0.13331d0, 0.05666d0, & !020
    1.99483d0, 0.61842d0, 0.08172d0, & !200
    0.99680d0,-0.70521d0, 0.08557d0, & !110
   -0.08900d0, 0.05328d0,-0.00271d0, & !030
    0.53619d0,-0.20959d0,-0.11892d0, & !300
    0.41179d0,-1.08201d0,-0.06926d0, & !012
    0.41179d0, 0.35624d0,-0.10806d0, & !102
   -0.08290d0, 0.00892d0, 0.08190d0, & !120
    1.62467d0,-0.35495d0, 0.06999d0, & !210
    0.04290d0,-0.22179d0,-0.05890d0, & !004
   -0.02203d0,-0.03182d0,-0.00800d0, & !040
   -0.62916d0,-0.51709d0,-0.11111d0, & !400
   -0.20222d0, 0.34644d0, 0.12995d0, & !022
    3.29320d0, 0.01374d0,-0.09119d0, & !112
   -2.27559d0,-1.25979d0,-0.48361d0, & !202
    0.15971d0,-0.03825d0,-0.03088d0, & !220
    1.20526d0, 0.18326d0,-0.09411d0, & !310
   -0.18405d0, 0.17426d0, 0.05992d0  & !130
    /), (/3, ncoeff_diag/) ) * 0.89875517923d0 
  INTEGER, PARAMETER, DIMENSION(3, ncoeff_diag) :: a_index_diag = RESHAPE( (/ &
    0, 0, 0, &
    0, 1, 0, &
    1, 0, 0, &
    0, 0, 2, &
    0, 2, 0, &
    2, 0, 0, &
    1, 1, 0, &
    0, 3, 0, &
    3, 0, 0, &
    0, 1, 2, &
    1, 0, 2, &
    1, 2, 0, &
    2, 1, 0, &
    0, 0, 4, &
    0, 4, 0, &
    4, 0, 0, &
    0, 2, 2, &
    1, 1, 2, &
    2, 0, 2, &
    2, 2, 0, &
    3, 1, 0, &
    1, 3, 0  &
    /), (/3, ncoeff_diag/) )
  DOUBLE PRECISION, PARAMETER, DIMENSION(ncoeff_off_diag) :: a_off_diag = (/ &
     1.17904d0, & !001
    -0.37700d0, & !011
     1.85559d0, & !101
     0.00698d0, & !003
    -0.36099d0, & !021
    -0.00202d0, & !111
     0.62070d0, & !201
     0.53070d0, & !013
    -0.64252d0, & !103
    -0.63496d0, & !121
    -0.12823d0, & !031
     1.40356d0, & !211
    -1.24044d0  & !301
    /) * 0.89875517923d0 
  INTEGER, PARAMETER, DIMENSION(3, ncoeff_off_diag) :: a_index_off_diag = RESHAPE( (/ &
    0, 0, 1, &
    0, 1, 1, &
    1, 0, 1, &
    0, 0, 3, &
    0, 2, 1, &
    1, 1, 1, &
    2, 0, 1, &
    0, 1, 3, &
    1, 0, 3, &
    1, 2, 1, &
    0, 3, 1, &
    2, 1, 1, &
    3, 0, 1  &
    /), (/3, ncoeff_off_diag/) )
  !================================================

  !================================================
  ! Useful arrays needed globally.
  !================================================
  INTEGER :: nmol

  DOUBLE PRECISION :: alpha(nat/3, 3, 3) !Array of permanent molecular polarizabilities. (nmol, xyz, xyz)
  DOUBLE PRECISION :: atoms_charged(nat, 3) !List of atoms, with O replaced by M site ('charged atoms'). (natoms, xyz)
  DOUBLE PRECISION :: dip_der_full(nat, 3, 3) !Gradient of dipole moment. (natoms, xyz, xyz)
  DOUBLE PRECISION :: dalpha_dr(nat, 3, 3, 3) !Gradient of permanent molecular polarizabilities. (natoms, xyz, xyz, xyz)

  DOUBLE PRECISION :: dip_ind(3) !Induced dipole moment. (nmol, xyz)
  DOUBLE PRECISION :: dip_ind_der(nat, 3, 3) !Gradient of the induced dipole moment. (nmol, nat, xyz, xyz)
  DOUBLE PRECISION :: pol_ind(3, 3) !Induced polarizability. (nmol, nmol, xyz, xyz)
  
  !Number of molecules.
  nmol = nat / 3

  !---------------------------------------------------
  ! Molecular polarizability.
  !---------------------------------------------------
  ! Computes molecular polarizability in the molecular frame (or sets it to constant for rigid models)
  ! and rotates it to the lab frame. Performed on each molecule.
  ! Also evaluates the analytical gradient of the polarizability, which is needed for the induced dipole gradients.
  !---------------------------------------------------
  CALL calc_alpha
  !---------------------------------------------------
  
  !-----------------------------------------------------------------------------------------------------
  ! Charged atoms with site M instead of oxygen, as defined by TIP4P force field.
  !-----------------------------------------------------------------------------------------------------
  CALL tip4p_particles
  !------------------------------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------------------------------
  ! Coordinate-dependent induced dipoles, derivatives of dipoles, and polarizabilities 
  ! (dip_ind, dip_ind_der, and pol_ind).
  ! This is where we perform Ewald summation.
  !------------------------------------------------------------------------------------------------------
  CALL calc_induced_part(dip_ind, dip_ind_der, pol_ind)
  !------------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------------
  ! Calculate total dipole moment. Using scaled charges for permanent dipole according to Hamm.
  !------------------------------------------------------------------------------------------------------
  dip(:) = calc_dipole(atoms_charged, charges / 1.3d0, dip_ind)
  !------------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------------
  !Calculate total dipole moment derivative. Using scaled charges for permanent dipole according to Hamm.
  !------------------------------------------------------------------------------------------------------
  IF (compute_der) THEN
     dip_der_full(:, : ,:) = calc_dipole_derivative(charges / 1.3d0, dip_ind_der)
  ELSE
     dip_der_full(:, :, :) = 0.0d0
  ENDIF
  ! Only z-component of the dipole gradient is used (arbitrary choice).
  dip_der = dip_der_full(:, 3, :)
  !------------------------------------------------------------------------------------------------------

  !--------------------------------------
  !Calculate total polarizability.
  !--------------------------------------
  pol = calc_polarizability(pol_ind)
  !--------------------------------------

 CONTAINS

  !************************************************************************************************************
  ! Charged atom positions as defined for TIP4P force field: two H atoms and the M site that 
  ! is computed below from O and H coordinates.
  !************************************************************************************************************
  SUBROUTINE tip4p_particles

    INTEGER :: i, j
    
    atoms_charged = atoms
    DO i=1, nat, 3
       !Apply PBC to hydrogens so that the molecule is not split over the box.
       !Probably not needed if we do not expect PBC applied to atoms before passed from i-pi.
       DO j = 1, 2
          atoms_charged(i + j, :) = atoms(i + j, :) - NINT((atoms(i + j, :) - atoms(i, :))/ box(:)) * box(:) 
       ENDDO
       !Compute M site position.
       atoms_charged(i, :) = gam * atoms(i, :) + gam2 * (atoms_charged(i+1, :) + atoms_charged(i+2, :))
    ENDDO

  END SUBROUTINE tip4p_particles
  !************************************************************************************************************
  !************************************************************************************************************

  !************************************************************************************************************
  ! Computation of permanent molecular polarizabilities, first in the molecular frame and then rotated 
  ! into the lab frame. 
  ! Other procedures that are needed are right below calc_alpha and are enclosed with comment lines to 
  ! separate from other routines.
  !************************************************************************************************************
  SUBROUTINE calc_alpha

    INTEGER :: i, iatom

    DO i = 1, nmol
       iatom = 3 * (i - 1) + 1
       CALL calc_alpha_mol(alpha(i, :, :), dalpha_dr(iatom : iatom + 2, :, :, :), atoms(iatom : iatom + 2, :))
       CALL rotate_alpha(alpha(i, :, :), dalpha_dr(iatom : iatom + 2, :, :, :), atoms(iatom : iatom + 2, :))
    ENDDO

  END SUBROUTINE calc_alpha

  !************************************************************************************************************
  ! Compute molecular polarizability in the molecular frame. Several options
  ! available: isotropic rigid, anisotropic rigid, and anisotropic flexible. For
  ! the flexible polarizability, see routine below named calc_alpha_Avila.
  !************************************************************************************************************
  SUBROUTINE calc_alpha_mol(alpha_mol, dalpha_dr_mol, atoms_mol)

    DOUBLE PRECISION, INTENT(INOUT) :: alpha_mol(3, 3), dalpha_dr_mol(3, 3, 3, 3)
    DOUBLE PRECISION, INTENT(IN) :: atoms_mol(3, 3)

    DOUBLE PRECISION :: a_Avila(4), da_dr(3, 3, 4)
    INTEGER :: k

    alpha_mol = 0.0d0
    dalpha_dr_mol = 0.0d0
    CALL calc_alpha_Avila(a_Avila, da_dr, atoms_mol) !Anisotropic molecular polarizability.
    DO k = 1, 3
       alpha_mol(k, k) = a_Avila(k)               !Anisotropic flexible molecular polarizability (diagonal terms).
       dalpha_dr_mol(:, :, k, k) = da_dr(:, :, k) !Gradient of the anisotropic flexible molecular polarizability (diagonal terms).
    ENDDO
    alpha_mol(1,2) = a_Avila(4)                !Anisotropic flexible molecular polarizability (xy component).
    alpha_mol(2,1) = a_Avila(4)                !Anisotropic flexible molecular polarizability (xy component).
    dalpha_dr_mol(:, :, 1, 2) = da_dr(:, :, 4) !Gradient of the anisotropic flexible molecular polarizability (xy component).
    dalpha_dr_mol(:, :, 2, 1) = da_dr(:, :, 4) !Gradient of the anisotropic flexible molecular polarizability (xy component).

  END SUBROUTINE calc_alpha_mol

  !************************************************************************************************************
  ! Molecular polarizability in the molecular frame from G. Avila, JCP 122, 144310 (2005).
  ! Analytical gradient available, needed for computing the gradient of the induced dipole moment.
  !************************************************************************************************************
  SUBROUTINE calc_alpha_Avila(a_mol, da_dr_mol, atoms_mol)

    DOUBLE PRECISION, INTENT(INOUT) :: a_mol(4), da_dr_mol(3, 3, 4)
    DOUBLE PRECISION, INTENT(IN) :: atoms_mol(3, 3)

    DOUBLE PRECISION :: r1(3), r2(3), r1_norm, r2_norm, tnsr1(3, 3), tnsr2(3, 3)
    DOUBLE PRECISION :: delta_r1, delta_r2, cos_phi, phi
    DOUBLE PRECISION :: dalpha_ds(3, 4), dacos_dx
    DOUBLE PRECISION :: s1, s2, s3
    INTEGER :: i, j, k, l

    CALL normalized_vec_norm(atoms_mol(1,:) - atoms_mol(2,:), r1, r1_norm)
    CALL normalized_vec_norm(atoms_mol(1,:) - atoms_mol(3,:), r2, r2_norm)
    delta_r1 = r1_norm / angtoau - r_eq !In angstrom
    delta_r2 = r2_norm / angtoau - r_eq !In angstrom
    cos_phi = DOT_PRODUCT(r1, r2)
    phi = ACOS(cos_phi) !In radian.

    !All derived coordinates in angstrom, to match the units of the parameters defined in the paper.
    s1 = (delta_r1 + delta_r2) / sqrt2
    s2 = phi - phi_eq
    s3 = (delta_r2 - delta_r1) / sqrt2

    ! Inefficient evaluation of the polynomial, but clean code. Could be improved further for efficiency.
    a_mol = 0.0d0
    DO l = 1, ncoeff_diag
       a_mol(1:3) = a_mol(1:3) + a_diag(:, l) * s1**a_index_diag(1, l) * s2**a_index_diag(2, l) * s3**a_index_diag(3, l)
    ENDDO
    DO l = 1, ncoeff_off_diag
       a_mol(4) = a_mol(4) + a_off_diag(l) * s1**a_index_off_diag(1, l) * s2**a_index_off_diag(2, l) * s3**a_index_off_diag(3, l)
    ENDDO
    a_mol = a_mol * angtoau**3 ! Convert back to atomic units.

    IF (compute_der) THEN
       !Helper variables, all in atomic units.
       tnsr1 = rot_grad_tnsr(r1, r1_norm)
       tnsr2 = rot_grad_tnsr(r2, r2_norm)
       dacos_dx = - angtoau / SQRT(1.0d0 - cos_phi**2)

       ! Gradient with respect to internal coordinates. Inefficient code, but readable.
       dalpha_ds = 0.0d0
       DO l = 1, ncoeff_diag
          i = a_index_diag(1, l); j = a_index_diag(2, l); k = a_index_diag(3, l);
          IF (i > 0) dalpha_ds(1, 1:3) = dalpha_ds(1, 1:3) + a_diag(:, l) * i * s1**(i-1) * s2**j * s3**k
          IF (j > 0) dalpha_ds(2, 1:3) = dalpha_ds(2, 1:3) + a_diag(:, l) * j * s1**i * s2**(j-1) * s3**k
          IF (k > 0) dalpha_ds(3, 1:3) = dalpha_ds(3, 1:3) + a_diag(:, l) * k * s1**i * s2**j * s3**(k-1)
       ENDDO
       DO l = 1, ncoeff_off_diag
          i = a_index_off_diag(1, l); j = a_index_off_diag(2, l); k = a_index_off_diag(3, l);
          IF (i > 0) dalpha_ds(1, 4) = dalpha_ds(1, 4) + a_off_diag(l) * i * s1**(i-1) * s2**j * s3**k
          IF (j > 0) dalpha_ds(2, 4) = dalpha_ds(2, 4) + a_off_diag(l) * j * s1**i * s2**(j-1) * s3**k
          IF (k > 0) dalpha_ds(3, 4) = dalpha_ds(3, 4) + a_off_diag(l) * k * s1**i * s2**j * s3**(k-1)
       ENDDO

       ! Gradient with respect to Cartesian coordinates:
       ! Elements are da_dr_mol(i, j, k) = (ds_l / dr_ij) * (dalpha_k / ds_l), 
       ! where r_ij corresponds to j-th coordinate of i-th atom, alpha_k is k-th component of the polarizability tensor's diagonal,
       ! and s_l are s1, s2, and s3 coordinates (over which we sum).
       !Oxygen coordinates.
       da_dr_mol(1, :, :) =  outer((r1 + r2) / sqrt2, dalpha_ds(1, :)) &
                            + outer((MATMUL(tnsr1, r2) + MATMUL(tnsr2, r1)) * dacos_dx, dalpha_ds(2, :)) &
                            + outer((r2 - r1) / sqrt2, dalpha_ds(3, :))
       !Hydrogen 1 coordinates.
       da_dr_mol(2, :, :) = - outer(r1 / sqrt2, dalpha_ds(1, :)) &
                            - outer(MATMUL(tnsr1, r2) * dacos_dx, dalpha_ds(2, :)) &
                            + outer(r1 / sqrt2, dalpha_ds(3, :))
       !Hydrogen 2 coordinates.
       da_dr_mol(3, :, :) = - outer(r2 / sqrt2, dalpha_ds(1, :)) &
                            - outer(MATMUL(tnsr2, r1) * dacos_dx, dalpha_ds(2, :)) &
                            - outer(r2 / sqrt2, dalpha_ds(3, :))
       da_dr_mol = da_dr_mol * angtoau**2 ! Convert to atomic units.
    ENDIF

  END SUBROUTINE calc_alpha_Avila

  !************************************************************************************************************
  ! Rotation of all permanent molecular polarizabilities into the lab frame.
  ! Definition of axes from G. Avila, JCP 122, 144310 (2005).
  ! y is the axis that bisects the H-O-H angle and molecules lies in the xy plane.
  ! Following this definition, x and y are computed from r1 = normalized(rO - rH1) and r2 = normalized(rO-rH2) as
  ! x = normalized(r1 -r2), y = normalized(r1+r2).
  ! z is perpendicular to the plane of the molecule and can be computed as the cross product of x and y.
  !************************************************************************************************************
  SUBROUTINE rotate_alpha(alpha_mol, dalpha_dr_mol, atoms_mol)

    DOUBLE PRECISION, INTENT(INOUT) :: alpha_mol(3, 3), dalpha_dr_mol(3, 3, 3, 3)
    DOUBLE PRECISION, INTENT(IN) :: atoms_mol(3, 3)

    DOUBLE PRECISION :: O(3, 3), dO_dr(3, 3, 3, 3)
    DOUBLE PRECISION :: r1(3), r2(3), r_sum(3), r_diff(3), r1_norm, r2_norm, r_sum_norm, r_diff_norm
    DOUBLE PRECISION :: tnsr1(3, 3), tnsr2(3, 3), tnsr_sum(3, 3), tnsr_diff(3, 3), tnsr_tmp(3, 3)
    INTEGER :: j, k

    CALL normalized_vec_norm(atoms_mol(1,:) - atoms_mol(2,:), r1, r1_norm)
    CALL normalized_vec_norm(atoms_mol(1,:) - atoms_mol(3,:), r2, r2_norm)
    CALL normalized_vec_norm(r1 + r2, r_sum, r_sum_norm)
    CALL normalized_vec_norm(r1 - r2, r_diff, r_diff_norm)

    !x, y, and z molecular axes in the laboratory frame.
    O(:, 1) = r_diff
    O(:, 2) = r_sum
    O(:, 3) = cross_product(O(:, 1), O(:, 2))

    IF (compute_der) THEN
       tnsr1 = rot_grad_tnsr(r1, r1_norm)
       tnsr2 = rot_grad_tnsr(r2, r2_norm)
       tnsr_sum = rot_grad_tnsr(r_sum, r_sum_norm)
       tnsr_diff = rot_grad_tnsr(r_diff, r_diff_norm)
       !Derivatives of the x axis w.r.t. O, H1, and H2 atomic coordinates.
       dO_dr(1, :, :, 1) =  MATMUL(tnsr1-tnsr2, tnsr_diff)     
       dO_dr(2, :, :, 1) = -MATMUL(tnsr1, tnsr_diff)       
       dO_dr(3, :, :, 1) =  MATMUL(tnsr2, tnsr_diff)       
       !Derivatives of the y axis w.r.t. O, H1, and H2 atomic coordinates.
       dO_dr(1, :, :, 2) =  MATMUL(tnsr1 + tnsr2, tnsr_sum )
       dO_dr(2, :, :, 2) = -MATMUL(tnsr1, tnsr_sum )        
       dO_dr(3, :, :, 2) = -MATMUL(tnsr2, tnsr_sum )        
       !Derivatives of the z axis w.r.t. atomic coordinates.
       DO j = 1, 3
          DO k = 1, 3
             dO_dr(j, k, :, 3) = cross_product(dO_dr(j, k, :, 1), O(:, 2)) + cross_product(O(:, 1), dO_dr(j, k, :, 2)) 
          ENDDO
       ENDDO

      !Rotated polarizability gradient. Uses molecular-frame alpha_mol, so must be evaluated before rotating polarizability.
       DO j = 1, 3
          DO k = 1, 3
             tnsr_tmp = multiply_matrices(dO_dr(j, k, :, :), alpha_mol(:, :), TRANSPOSE(O))
             dalpha_dr_mol(j, k, :, :) = tnsr_tmp + TRANSPOSE(tnsr_tmp) + multiply_matrices(O, dalpha_dr_mol(j, k, :, :), TRANSPOSE(O))
          ENDDO
       ENDDO
    ENDIF

    !Rotated polarizability.
    alpha_mol(:, :) = multiply_matrices(O, alpha_mol(:, :), TRANSPOSE(O))

  END SUBROUTINE rotate_alpha

  ! Computes matrix (Id - outer(vec, vec)) / norm, which appears in the computation of the gradient of the rotation matrix.
  FUNCTION rot_grad_tnsr(vec, norm) RESULT(tnsr)

    DOUBLE PRECISION, INTENT(IN) :: vec(3), norm

    DOUBLE PRECISION :: tnsr(3, 3)

    INTEGER :: k

    tnsr = - outer(vec, vec)
    DO k = 1, 3
       tnsr(k, k) = tnsr(k, k) + 1.0d0
    ENDDO
    tnsr = tnsr / norm

  END FUNCTION rot_grad_tnsr
  !************************************************************************************************************
  !************************************************************************************************************

  !************************************************************************************************************
  ! Induced dipoles, dipole gradients, and polarizabilities. We use the first-order (not self-consistently solved) 
  ! dipole-induced dipole model, accounting for the intermolecular (two-body) interactions.
  ! This computation involves multiple routines that are all gathered below and separated from the rest 
  ! of the code by comment lines.
  !************************************************************************************************************
  SUBROUTINE calc_induced_part(dip_ind, dip_ind_der, pol_ind)

    DOUBLE PRECISION, INTENT(INOUT) :: dip_ind(3), dip_ind_der(nat, 3, 3), pol_ind(3, 3)

    DOUBLE PRECISION :: a, rcut
    DOUBLE PRECISION :: ro(nat/3, 3) !List of O atoms. (nmol, xyz)

    !-------------------------------------------------------
    ! It is practical to have a list of O atom coordinates.
    !-------------------------------------------------------
    ro(:, :) = atoms(1:nat:3, :)

    !----------------------------------------------
    ! Parameters for Ewald (calculated only once).
    !----------------------------------------------
    rcut = ewald_param * MINVAL(box) * MIN(0.5d0,1.2d0*nat**(-1.d0/6.d0))
    a = ewald_param * pi/rcut

    !----------------------------------------------
    ! Short-range part - sum over pairs.
    !----------------------------------------------
    CALL short_range_ew(atoms_charged, ro, a, rcut, dip_ind, dip_ind_der, pol_ind)

    !----------------------------------------------
    ! Long-range part - performs sum in k space.
    !----------------------------------------------
    CALL long_range_ew(atoms_charged, ro, a, dip_ind, dip_ind_der, pol_ind)

  END SUBROUTINE calc_induced_part

  !************************************************************************************************************
  ! Short-range part, solved by direct summation in the real space.
  !************************************************************************************************************
  SUBROUTINE short_range_ew(r, ro, a, rcut, dip_ind, dip_ind_der, pol_ind)

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nmol, 3), a, rcut 
    DOUBLE PRECISION, INTENT(INOUT) :: dip_ind(3), dip_ind_der(nat, 3, 3), pol_ind(3, 3)

    INTEGER :: i, j, k, l, m, iatom, jatom, latom, nx, ny, nz, nmax
    DOUBLE PRECISION :: r_ij_0(3), r_ij(3), dr, dr2, dr3, rcut2, screen, a_dr, gauss_part
    DOUBLE PRECISION :: T_tnsr(3, 3), T_tnsr_sum(3, 3), E_ij(3), a3_self_term, prefac
    LOGICAL :: self_term

    nmax=NINT(rcut/MINVAL(box))
    dip_ind = 0.0d0
    dip_ind_der = 0.0d0
    rcut2 = rcut**2
    prefac = 2.0d0 / sqrtpi
    DO i = 1, nmol
       iatom = 3 * (i - 1) + 1
       T_tnsr_sum = 0.0d0
       E_ij = 0.0d0
       DO j = 1, nmol
          jatom = 3 * (j - 1)
          DO k = 1, 3
             T_tnsr = 0.0d0
             !Nearest neighbor.
             r_ij_0(:) = ro(i,:) - r(jatom + k, :)
             r_ij_0(:) = r_ij_0(:) - box(:) * NINT(r_ij_0(:)/box(:))
             DO nx = -nmax, nmax
                DO ny = -nmax, nmax
                   DO nz = -nmax, nmax
                      r_ij(:) = r_ij_0(:) - box(:) * (/nx, ny, nz/) 
                      dr2 = SUM(r_ij**2)
                      IF (dr2 .LT. rcut2) THEN
                         ! Helper variables that are used more than once or convenient to compute separately.
                         self_term = (i .EQ. j) .AND. (nx .EQ. 0) .AND. (ny .EQ. 0).AND. (nz .EQ. 0)
                         CALL compute_helper_vars(a, dr2, prefac, self_term, dr, dr3, a_dr, gauss_part, screen)
                         ! Contribution to the electric field on molecule i induced by atom k in molecule j of cell n.
                         E_ij(:) = E_ij(:) + charges(k) * r_ij(:) / dr3 * screen
                         ! Contribution to the T tensor for given i, j, k, n (used for dipole derivative).
                         IF (compute_der) T_tnsr = T_tnsr + short_range_T_tnsr(r_ij, dr2, dr3, a_dr, gauss_part, screen)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             !--------------------------------------------------------------------------------------------------------
             ! Compute the derivative of the electric field of molecule i w.r.t. coordinates of atom k of molecule j.
             !--------------------------------------------------------------------------------------------------------
             IF (compute_der) THEN
                T_tnsr =  charges(k) * MATMUL(alpha(i, :, :), T_tnsr)
                !Derivative of induced dipole at molecule i w.r.t. coordinates of molecule j (includes self-term i==j).
                CALL dip_ind_der_ij(dip_ind_der(jatom + 1: jatom + 3, :, :), T_tnsr(:, :), k)
                !Sum for derivative w.r.t. i-th molecule atoms.
                T_tnsr_sum = T_tnsr_sum + T_tnsr
             ENDIF
             !--------------------------------------------------------------------------------------------------------
             !--------------------------------------------------------------------------------------------------------
          ENDDO
       ENDDO
       !----------------------------------------------------------------------------------------------------------
       ! Contribution to the dipole moment of molecule i induced by atom k in molecule j summed over all cells.
       !----------------------------------------------------------------------------------------------------------
       dip_ind(:) = dip_ind(:) + MATMUL(alpha(i, :, :), E_ij(:))
       !----------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------

       !--------------------------------------------------------------------------------------------------------
       ! Compute the derivative of the electric field of molecule i w.r.t. coordinates of atoms in molecule i.
       !--------------------------------------------------------------------------------------------------------
       IF (compute_der) THEN
          !Oxygen of molecule i.
          dip_ind_der(iatom, :, :) = dip_ind_der(iatom, :, :) + T_tnsr_sum(:, :)
          !Derivative of alpha w.r.t. atoms of i-th molecule (loop over oxygen and hydrogens, and over xyz).
          DO l = 1, 3
             latom = iatom + l - 1
             DO m = 1, 3
                dip_ind_der(latom, m, :) = dip_ind_der(latom, m, :) + MATMUL(dalpha_dr(latom, :, m, :), E_ij)
             ENDDO
          ENDDO
       ENDIF
       !--------------------------------------------------------------------------------------------------------
       !--------------------------------------------------------------------------------------------------------
    ENDDO

    !--------------------------------------------------------------------------------------------------------
    ! Compute induced part of the polarizability.
    !--------------------------------------------------------------------------------------------------------
    a3_self_term = 4.0d0 * a**3 / (3.0d0 * sqrtpi)
    pol_ind = 0.0d0
    DO i = 1, nmol
       DO j = i, nmol
          r_ij_0(:) = ro(i,:) - ro(j, :)
          r_ij_0(:) = r_ij_0(:) - box(:) * NINT(r_ij_0(:)/box(:))
          T_tnsr = 0.0d0
          !Sum over cells and consider only non-self terms.
          DO nx = -nmax, nmax
             DO ny = -nmax, nmax
                DO nz = -nmax, nmax
                   r_ij(:) = r_ij_0(:) - box(:) * (/nx, ny, nz/) 
                   dr2 = SUM(r_ij**2)
                   self_term = (i .EQ. j) .AND. (nx .EQ. 0) .AND. (ny .EQ. 0).AND. (nz .EQ. 0)
                   IF (dr2 .LT. rcut2 .AND. (.NOT. self_term)) THEN
                      CALL compute_helper_vars(a, dr2, prefac, .FALSE., dr, dr3, a_dr, gauss_part, screen)
                      T_tnsr = T_tnsr + short_range_T_tnsr(r_ij, dr2, dr3, a_dr, gauss_part, screen)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          !Multiply computed tensor by alpha and add to system's polarizability.
          IF (i .EQ. j) THEN
             !Self term, include one for each value of i. 
             pol_ind = pol_ind - a3_self_term * MATMUL(alpha(i, :, :), alpha(j, :, :))
          ELSE
             !Not a self term, so we have to add also the transpose, to correct for including only j>i.
             ! alpha_i * T_ij * alpha_j = (alpha_j * T_ji * alpha_i)^T
             T_tnsr = MATMUL(MATMUL(alpha(i, :, :), T_tnsr), alpha(j, :, :))
             pol_ind = pol_ind + T_tnsr + TRANSPOSE(T_tnsr)
          ENDIF
       ENDDO
    ENDDO
    !--------------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------

  END SUBROUTINE short_range_ew

  !************************************************************************************************************
  ! Used for simplifying the short-range Ewald code above.
  !************************************************************************************************************
  SUBROUTINE compute_helper_vars(a, dr2, prefac, self_term, dr, dr3, a_dr, gauss_part, screen)

    DOUBLE PRECISION, INTENT(IN) :: a, dr2, prefac
    LOGICAL, INTENT(IN) :: self_term
    DOUBLE PRECISION, INTENT(OUT) :: dr, dr3, a_dr, gauss_part, screen

    ! Helper variables for real-space summation:
    dr = SQRT(dr2)                                                  !Absolute value of distance.
    dr3 = dr*dr2                                                    !dr cubed.
    a_dr = a * dr                                                   !a * dr.
    gauss_part = prefac * a_dr * EXP(-a_dr**2)                      !Gaussian multiplied by a_dr and the prefactor.
    screen = short_range_ew_screen(a_dr, gauss_part, self_term)     !screen term equal to erfc(x) + gauss_part (-1 if self term).

  END SUBROUTINE compute_helper_vars

  !************************************************************************************************************
  ! Screening term from the short-range Ewald summation. 
  !************************************************************************************************************
  FUNCTION short_range_ew_screen(a_r, gauss_part, self_term) RESULT(sc)

    DOUBLE PRECISION :: sc

    DOUBLE PRECISION, INTENT(IN) :: a_r, gauss_part
    LOGICAL, INTENT(IN) :: self_term

    sc = erfc(a_r) + gauss_part
    IF (self_term) sc = sc - 1.0d0 !Self-interaction term.

  END FUNCTION short_range_ew_screen

  !************************************************************************************************************
  ! Standard dipole-induced dipole interaction tensor with short-range Ewald screening.
  !************************************************************************************************************
  FUNCTION short_range_T_tnsr(r_ij, r2, r3, a_r, gauss_part, screen) RESULT(T_ij)

    DOUBLE PRECISION :: T_ij(3, 3)

    DOUBLE PRECISION, INTENT(IN) :: r_ij(3), r2, r3, a_r, gauss_part, screen

    INTEGER :: l

    T_ij = - outer(r_ij, r_ij) / (r3 * r2) * (3.0d0 * screen + 2.0d0 * a_r**2 * gauss_part)
    DO l = 1, 3
       T_ij(l, l) = T_ij(l, l) + screen / r3
    ENDDO

  END FUNCTION short_range_T_tnsr

  !************************************************************************************************************
  ! Long-range part solved by performing operations in the k-space. 
  !************************************************************************************************************
  SUBROUTINE long_range_ew(r, ro, a, dip_ind, dip_ind_der, pol_ind)

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nmol, 3), a    
    DOUBLE PRECISION, INTENT(INOUT) :: dip_ind(3), dip_ind_der(nat, 3, 3), pol_ind(3, 3)

    INTEGER :: i, l, m, iatom, latom, kx, ky, kz, kmax, k
    DOUBLE PRECISION :: b, f, f_mod, rk(3), rk_i, rk_out(3, 3), rk2, rkmax2, lat(3), q(nat), tnsr_tmp(3, 3), prefac, sum_alpha(3, 3) 
    DOUBLE COMPLEX :: sk_i(nat), exp_ikr(nmol), sk, sk_o(3, 3), sk_o_times_rk(3), exp_ikr_i_times_sk, sk_o_times_rk_out(3, 3)
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:, :, :) :: exp_ikr_charged, exp_ikr_o

    kmax = INT(ewald_param * a * MAXVAL(box))
    lat(:) = twopi/box(:) 
    rkmax2 = (ewald_param * twopi * a)**2
    b = 0.25d0/(a**2)
    f = 4.0d0 * pi / PRODUCT(box) ! 4pi / V !
    f_mod = f
    q(:) = PACK(SPREAD(charges(:), 2, nmol), .TRUE.)

    ALLOCATE(exp_ikr_charged(1:nat, -kmax:kmax, 1:3), exp_ikr_o(1:nmol, -kmax: kmax, 1:3))

    ! Evaluate exp(i k_i r_i) terms for i = x, y, z. 
    ! Note: It is much faster to precompute and use in the computationally 
    ! intensive part below than to re-evaluate the exponentials kmax^3 * nat times.
    DO k = -kmax, kmax
       DO i = 1, 3
          rk_i = lat(i) * k
          exp_ikr_charged(:, k, i) = EXP(-IU * rk_i * r(:, i))
          exp_ikr_o(:, k, i) = EXP(IU * rk_i * ro(:, i))
       ENDDO
    ENDDO

    DO kx = 0,kmax
       IF (kx .EQ. 1) f_mod = 2.d0*f
       rk(1) = lat(1)*kx
       DO ky = -kmax,kmax
          rk(2) = lat(2)*ky
          DO kz = -kmax,kmax
             rk(3) = lat(3)*kz
             rk2 = SUM(rk(:)**2)
             IF (rk2 .LT. rkmax2 .AND. rk2 .GT. EPSILON(0.d0)) THEN
                !Helper variables:
                prefac = f_mod * EXP(-b*rk2) / rk2         !Prefactor common to all properties.
                !charge(i) * exp(i k r_i), where r_i are charged sites.
                sk_i = q * exp_ikr_charged(:, kx, 1) * exp_ikr_charged(:, ky, 2) * exp_ikr_charged(:, kz, 3)
                sk = SUM(sk_i)                              !Structure factor for the charged sites.
                !exp(i k r_o), where r_o are oxygen atoms.
                exp_ikr = exp_ikr_o(:, kx, 1) * exp_ikr_o(:, ky, 2) * exp_ikr_o(:, kz, 3)
                sk_o = 0.0d0                               !Structure factor for the oxygen atoms.
                DO i = 1, nmol
                   sk_o = sk_o + alpha(i, :, :) * exp_ikr(i) 
                ENDDO

                !--------------------------------------------------------------------------------------------
                ! Contribution to the induced dipole moment of the system.
                !--------------------------------------------------------------------------------------------
                dip_ind = dip_ind +  prefac * AIMAG(sk * MATMUL(sk_o, rk))
                !--------------------------------------------------------------------------------------------
                !--------------------------------------------------------------------------------------------

                !--------------------------------------------------------------------------------------------
                ! Contribution to the derivative of the induced dipole moment of the system.
                !--------------------------------------------------------------------------------------------
                IF (compute_der) THEN
                   rk_out = prefac * outer(rk, rk)
                   sk_o_times_rk_out = MATMUL(sk_o, rk_out)
                   DO i = 1, nmol
                      exp_ikr_i_times_sk = exp_ikr(i) * sk
                      iatom = 3 * (i-1) + 1
                      !Sum for i-th oxygen atom.
                      dip_ind_der(iatom, :, :) = dip_ind_der(iatom, :, :) + MATMUL(alpha(i, :, :), rk_out) * REAL(exp_ikr_i_times_sk, KIND=KIND(1.0d0))
                      DO l = 1, 3
                         latom = iatom + l - 1
                         !Derivatives with respect to atom l of molecule j.
                         tnsr_tmp = REAL(sk_i(latom) * sk_o_times_rk_out, KIND=KIND(1.0d0))
                         CALL dip_ind_der_ij(dip_ind_der(iatom : iatom + 2, :, :), tnsr_tmp, l)
                         !Contribution from the derivative of alpha.
                         DO m = 1, 3
                            dip_ind_der(latom, m, :) = dip_ind_der(latom, m, :) + prefac * MATMUL(dalpha_dr(latom, :, m, :), rk) * AIMAG(exp_ikr_i_times_sk)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                !--------------------------------------------------------------------------------------------
                !--------------------------------------------------------------------------------------------

                !--------------------------------------------------------------------------------------------
                ! Contribution to the induced polarizability of the system.
                !--------------------------------------------------------------------------------------------
                sk_o_times_rk = MATMUL(sk_o, rk)
                pol_ind = pol_ind + prefac * REAL(outer_cmplx(sk_o_times_rk, CONJG(sk_o_times_rk)), KIND=KIND(1.0d0))
                !--------------------------------------------------------------------------------------------
                !--------------------------------------------------------------------------------------------
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !--------------------------------------------------------------------------------------------------
    ! Non-zero constant contribution to the induced polarizability of the system, from kx=ky=kz=0
    !--------------------------------------------------------------------------------------------------
    sum_alpha = SUM(alpha, DIM=1)
    pol_ind = pol_ind + f * MATMUL(sum_alpha, sum_alpha) / 3.0d0
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------

  END SUBROUTINE long_range_ew

  !************************************************************************************************************
  ! Derivative of induced dipole at molecule i w.r.t. coordinates of molecule j (includes self-term i==j).
  ! This code is common to short-range and long-range Ewald components, used in both subroutines.
  !************************************************************************************************************
  SUBROUTINE dip_ind_der_ij(dmui_drj, tnsr, k)

    DOUBLE PRECISION, INTENT(INOUT) :: dmui_drj(3, 3, 3)
    DOUBLE PRECISION, INTENT(IN) :: tnsr(3, 3)
    INTEGER, INTENT(IN) :: k

    INTEGER :: h

    IF (k .EQ. 1) THEN
       !j-th molecule oxygen atom. Includes i==j term.
       dmui_drj(k, :, :) = dmui_drj(k, :, :) - gam * tnsr(:, :)
       !j-th molecule hydrogen atoms (part that comes from M site).
       DO h = 1, 2
          dmui_drj(k + h, :, :) = dmui_drj(k + h, :, :) - gam2 * tnsr(:, :)
       ENDDO
    ELSE
       !j-th molecule hydrogen atoms (part that comes from hydrogen atom positions).
       dmui_drj(k, :, :) = dmui_drj(k, :, :) - tnsr(:, :)
    ENDIF   

  END SUBROUTINE dip_ind_der_ij
  !************************************************************************************************************
  !************************************************************************************************************

  !************************************************************************************************************
  ! Full dipole moment. One-body term is here computed only from the point
  ! charges, but could be replaced by a more sophisticated dipole surface.
  !************************************************************************************************************
  FUNCTION calc_dipole(atoms, charges, dip_ind)

    DOUBLE PRECISION, INTENT(IN) :: atoms(nat,3), charges(3), dip_ind(3)

    DOUBLE PRECISION :: calc_dipole(3)

    INTEGER :: i, iatom, k

    calc_dipole = 0.0d0
    !Permanent dipoles.
    DO i = 1, nmol
       iatom = 3 * (i - 1)
       DO k = 1, 3
          calc_dipole(:) = calc_dipole(:) + charges(k) * atoms(iatom + k, :) 
       ENDDO
    ENDDO
    !Induced dipoles from electrostatic interaction with charges on other water molecules.
    calc_dipole = calc_dipole + dip_ind

  END FUNCTION calc_dipole
  !************************************************************************************************************
  !************************************************************************************************************

  !************************************************************************************************************
  ! Full derivative of the dipole moment. One-body term is here computed only from the point
  ! charges, but could be replaced by a more sophisticated dipole surface.
  !************************************************************************************************************
  FUNCTION calc_dipole_derivative(charges, dip_ind_der) RESULT(dip_der)

    DOUBLE PRECISION, INTENT(IN) :: charges(3)
    DOUBLE PRECISION, INTENT(IN) :: dip_ind_der(nat, 3, 3)

    DOUBLE PRECISION :: dip_der(nat, 3, 3)

    INTEGER :: iatom, i, k

    dip_der = 0.0d0
    DO i=1, nmol
       iatom = 3 * (i - 1) + 1
       DO k = 1, 3
          !Gradient of the permanent molecular dipole moment.
          dip_der(iatom, k, k) = dip_der(iatom, k, k) + gam * charges(1) !i-th oxygen.
          dip_der(iatom+1:iatom+2, k, k) = dip_der(iatom+1:iatom+2, k, k) + charges(2) + gam2 * charges(1) !i-th hydrogen.
       ENDDO
    ENDDO
    dip_der = dip_der + dip_ind_der

  END FUNCTION calc_dipole_derivative
  !************************************************************************************************************
  !************************************************************************************************************

  !************************************************************************************************************
  ! Full polarizability as a sum of permanent and induced components. The
  ! negative sign in front of pol_ind is due to the definition of the dipole tensor used in the code above.
  !************************************************************************************************************
  FUNCTION calc_polarizability(pol_ind) RESULT(pol)

    DOUBLE PRECISION, INTENT(IN) :: pol_ind(3,3)

    DOUBLE PRECISION :: pol(3, 3)

    pol = 0.0d0
    pol = SUM(alpha, DIM=1)
    pol = pol - pol_ind 

  END FUNCTION calc_polarizability
  !************************************************************************************************************
  !************************************************************************************************************

  !************************************************************************************************************
  ! Below are some helper procedures that are used in different places.
  !************************************************************************************************************

  !--------------------------------------
  ! Evaluates array k.r_i (i = 1, n),
  ! where k and r_i are 3-dim vectors.
  !--------------------------------------
  FUNCTION k_dot_r(n, k, r)

    DOUBLE PRECISION :: k_dot_r(n)

    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: k(3), r(n, 3)

    INTEGER :: i

    k_dot_r = 0.0d0
    DO i = 1, 3
       k_dot_r(:) = k_dot_r(:) + k(i) * r(:, i) 
    ENDDO

  END FUNCTION k_dot_r

  !--------------------------------------
  ! Outer product of two vectors.
  ! c_ij = a_i * b_j
  !--------------------------------------
  FUNCTION outer(a, b)

    DOUBLE PRECISION, INTENT(IN) :: a(:), b(:)

    DOUBLE PRECISION :: outer(SIZE(a,1),SIZE(b,1))

    INTEGER :: i

    DO i=1, SIZE(a, 1)
       outer(i, :)  = a(i) * b(:)
    ENDDO

  END FUNCTION outer

  !--------------------------------------
  ! Outer product of two complex vectors.
  ! c_ij = a_i * b_j
  !--------------------------------------
  FUNCTION outer_cmplx(a, b)

    DOUBLE COMPLEX, INTENT(IN) :: a(3), b(3)

    DOUBLE COMPLEX :: outer_cmplx(3,3)

    INTEGER :: i

    DO i=1, 3
       outer_cmplx(i, :)  = a(i) * b(:)
    ENDDO

  END FUNCTION outer_cmplx

  !--------------------------------------
  ! Cross product of two 3-dim vectors.
  !--------------------------------------
  FUNCTION cross_product(a, b) RESULT(c)

    DOUBLE PRECISION, INTENT(IN) :: a(3), b(3)

    DOUBLE PRECISION :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  END FUNCTION cross_product

  !--------------------------------------
  ! Product of three matrices.
  !--------------------------------------
  FUNCTION multiply_matrices(A, B, C) RESULT(ABC)

    DOUBLE PRECISION, INTENT(IN) :: A(3, 3), B(3, 3), C(3, 3)

    DOUBLE PRECISION :: ABC(3, 3)

    ABC = MATMUL(A, MATMUL(B, C))

  END FUNCTION multiply_matrices

  !-----------------------------------------
  ! Takes a general vector x and computes
  ! its norm (x_norm) and normalized vector 
  ! x_vec = x / x_norm.
  !-----------------------------------------
  SUBROUTINE normalized_vec_norm(x, x_vec, x_norm)

    DOUBLE PRECISION, INTENT(IN) :: x(:)
    DOUBLE PRECISION, INTENT(OUT) :: x_vec(:), x_norm

    x_norm = NORM2(x)
    x_vec = x / x_norm

  END SUBROUTINE normalized_vec_norm
  !************************************************************************************************************
  !************************************************************************************************************

END SUBROUTINE h2o_dipole

