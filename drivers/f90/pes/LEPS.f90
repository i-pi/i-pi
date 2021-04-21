! These models are taken from 
! H. Jonsson, G. Mills and K.W. Jacobsen, in Classical and Quantum Dynamics in 
! Condensed Phase Simulations, edited by B. J. Berne, G. Ciccotti and D. F. Coker
! (World Scientific, Singapore 1998). p. 385
! 
! This probably the first occurrence to the "modern" NEB in the literature.


SUBROUTINE LEPS_M1(na,atoms,V,F)
! See the reference above for more details.
! You must specify three atoms lying on the x axis.
! The atom 1st is fixed.
! The second atom must be within the other two:
! 1st ------ 2nd ------ 3rd.
  IMPLICIT NONE
  REAL(8),PARAMETER :: acb(3) = (/ .05D0, .05D0, .3D0 /), &
       & d(3) = (/ 4.746D0, 4.746D0, 3.445D0 /), r_0=.742D0, alpha=1.942D0, &
       & step=1D-4
  REAL(8), INTENT(IN) :: atoms(3,3) ! A=> atoms(1), B=> atoms(2) C=> atoms(3)
  REAL(8), INTENT(OUT) :: V,F(3,3)
  REAL(8) :: V_plus, V_mins, moved_atoms(3,3)
  INTEGER, INTENT(IN) :: na
  INTEGER :: i

  CALL V_LEPS(atoms, acb, d, r_0, alpha, V)

  ! Compute numerical forces
  F = 0.0D0
  do i=2,na
     moved_atoms = atoms
     moved_atoms(i,1) = atoms(i,1) + step
     CALL V_LEPS(moved_atoms, acb, d, r_0, alpha, V_plus)
     moved_atoms(i,1) = atoms(i,1) - step
     CALL V_LEPS(moved_atoms, acb, d, r_0, alpha, V_mins)
     F(i,1) = 1 * (V_plus - V_mins)/ (4*step)
  end do


END SUBROUTINE LEPS_M1

SUBROUTINE LEPS_M2(na,atoms,V,F)
! See the reference above for more details.
! You must specify three atoms lying on the x axis 
! plus another atoms lying in the xy plane. See the
! picture below and set atomic coordinates accordingly.
! Atoms 1st and 3rd are fixed.
! The second atom must be within the other two:
! 1st ------ 2nd ------ 3rd
!             |
!             |
!            4th
  IMPLICIT NONE
  REAL(8),PARAMETER :: acb(3) = (/ .05D0, .05D0, .8D0 /), &
       & d(3) = (/4.746D0, 4.746D0, 3.445D0 /), r_0=.742D0, alpha=1.942D0, &
       & kc=0.2025D0, c=1.154, step=1D-4
  REAL(8), INTENT(IN) :: atoms(4,3) ! A=> atoms(1), B=> atoms(2) C=> atoms(3)
  REAL(8), INTENT(OUT) :: V,F(4,3)
  INTEGER, INTENT(IN) :: na
  REAL(8) :: rab, rac, x, dist, moved_atoms(4,3), V_plus, V_mins
  INTEGER :: i

  CALL V_LEPS(atoms(1:3,:), acb, d, r_0, alpha, V)
  rab = dist(atoms(1,:), atoms(2,:))
  rac = dist(atoms(1,:), atoms(3,:))
  x = atoms(2,1) - atoms(4,1)
  V = V + 2 * kc * (rab - (rac/2 - x/c))**2

  ! Compute numerical forces
  F = 0.0D0
  do i=2,na,2
     moved_atoms = atoms
     moved_atoms(i,1) = atoms(i,1) + step
     CALL V_LEPS(moved_atoms, acb, d, r_0, alpha, V_plus)
     moved_atoms(i,1) = atoms(i,1) - step
     CALL V_LEPS(moved_atoms, acb, d, r_0, alpha, V_mins)
     F(i,1) = -1 * (V_plus - V_mins)/ (4 * step)
  end do

  
END SUBROUTINE LEPS_M2

SUBROUTINE V_LEPS(atoms, acb, d, r_0, alpha, V)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: atoms(3,3), acb(3), d(3), r_0, alpha
  REAL(8), INTENT(OUT) :: V
  REAL(8) :: r, Jiii(3)=0D0, J_temp=0D0, dist
  INTEGER :: i, ii
  
  V = 0.0D0
  J_temp = 0.0D0
  do i = 1,3
     ii = mod(i,3)+1
     r = dist(atoms(i,:), atoms(ii,:))
     V = V + Q(r, d(i)) / (1 + acb(i))
     Jiii(i) = J(r, d(i)) / (1 + acb(i))
  end do
  do i = 1,3
     ii = mod(i,3)+1
     J_temp = J_temp + Jiii(i)**2 - Jiii(i)*Jiii(ii)
  end do
  
  V = V - J_temp**0.5
  
CONTAINS
  REAL(8) FUNCTION Q(r, d)
    REAL(8), INTENT(IN) :: d, r
    Q = d/2.0D0 * (3.0D0/2.0D0 * exp(-2.0D0 * alpha * (r-r_0)) &
         - exp(-1.0D0 * alpha * (r-r_0)))
  END FUNCTION Q

  REAL(8) FUNCTION J(r, d)
    REAL(8), INTENT(IN) :: r, d
    J = d/4.0D0 * (exp(-2.0D0 * alpha * (r-r_0)) &
         - 6.0D0 * exp(-1.0D0 * alpha * (r-r_0)))
  END FUNCTION J

END SUBROUTINE V_LEPS


REAL(8) FUNCTION DIST(p1,p2)
  IMPLICIT NONE
  REAL(8),INTENT(IN),DIMENSION(3) :: p1, p2

  dist = ((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2D0)**0.5D0

  
END FUNCTION DIST




! PROGRAM test
!   IMPLICIT NONE
!   real(8), ALLOCATABLE, DIMENSION(:,:) :: atoms, step1, step2, atoms_orig,F
!   real(8) :: V, dist
!   integer :: na, i, j
!   logical :: M1 = .true.
!   logical :: only_force = .false.

!   ! print*, 'asdasd'
!   ! stop
  
!   if (M1) then
!      na = 3
!      ALLOCATE(atoms(na,3), step1(na,3), step2(na,3), atoms_orig(na,3), F(na,3))
!      atoms = reshape((/ 1.D0, 1.D0, 1.D0, .0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /),(/ 3,3 /)) * .0D0

!      if (only_force) then
!         atoms = reshape((/ 0.D0, .5934D0, 40.1D0, .0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /),(/ 3,3 /))
!         CALL LEPS_M1(na, atoms, V, F)
!         write(*,*) V
!         do i = 1,3
!            write(*,'(3F12.6)')(F(i,j), j=1,3)
!         end do
!         stop
!      else
!         atoms_orig = atoms
!         step1 = reshape((/ 0.D0, 1D0, 1D0, .0D0, 0.D0, 0.D0, .0D0, 0.D0, 0.D0 /),(/ 3,3 /)) * 4D-2
!         step2 = reshape((/ 0.D0, 0D0, 1D0, .0D0, 0.D0, 0.D0, .0D0, 0.D0, 0.D0 /),(/ 3,3 /)) * 4D-2
        
!         ! atoms = reshape((/ 1.D0, 2.D0, 3.D0, .0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /),(/ 3,3 /))
!         ! CALL LEPS_M1_V(na,atoms,V,F)
        
!         do i=1,101
!            do j = 1,101
!               CALL LEPS_M1(na,atoms,V,F)
!               ! V = 0D0
!               ! write(*,*) atoms(1,1), atoms(2,1), atoms(3,1)
!               ! write(*,*) dist(atoms(3,:), atoms(2,:)), dist(atoms(1,:), atoms(2,:)), dist(atoms(1,:),atoms(3,:))
!               write(*,*) dist(atoms(3,:), atoms(2,:)), dist(atoms(1,:), atoms(2,:)), F(2,1), F(3,1), V
!               atoms = atoms + step1
!            end do
!            atoms = atoms_orig + DBLE(i)*step2
!            write(*,*)
!         end do
!      end if
!   else
!      na = 4
!      ALLOCATE(atoms(na,3), step1(na,3), step2(na,3), atoms_orig(na,3), F(na,3))
!      ! xA,xB,xC,Xd,yA,yB,yC,yD,...
!      atoms = reshape((/ 0.D0, 0.D0, 3.742D0, -2.D0, .0D0, .0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /),(/ 4,3 /)) !* .0D0
     
!      if (only_force) then

!         atoms = reshape((/ 0.D0, 3.0D0, 3.742D0, -2.0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0D0, 0D0 /),(/ 4,3 /))
!         CALL LEPS_M2(na, atoms, V, F)
!         write(*,*) V
!         do i = 1,4
!            write(*,'(3F12.6)')(F(i,j), j=1,3)
!         end do
!         stop

!      else
!         atoms_orig = atoms
!         step1 = reshape((/ 0.D0, 1D0, 0D0, 1D0, 0.D0, 0.D0, .0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /),(/ 4,3 /)) * 4D-2
!         step2 = reshape((/ 0.D0, 0D0, 0D0, 1D0, 0.D0, 0.D0, .0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /),(/ 4,3 /)) * 4D-2
        
!         ! atoms = reshape((/ 1.D0, 2.D0, 3.D0, 4.D0, 5D0, 6D0, 7.D0, 8.D0, 9.D0, 10.D0, 11.D0, 12.D0 /),(/ 4,3 /)) !* .0D0
!         ! write(*,*) atoms(2,:)
!         ! stop
        
!         do i=1,101
!            do j = 1,101
!               CALL LEPS_M2(na,atoms,V,F)
!               ! V = 0D0
!               ! write(*,*) atoms(1,1), atoms(2,1), atoms(3,1)
!               ! write(*,*) dist(atoms(3,:), atoms(2,:)), dist(atoms(1,:), atoms(2,:)), dist(atoms(1,:),atoms(3,:))
!               write(*,*) dist(atoms(1,:), atoms(2,:)), atoms(2,1)-atoms(4,1), F(2,1) !V
!               atoms = atoms + step1
!            end do
!            atoms = atoms_orig + DBLE(i)*step2
!            write(*,*)
!         end do
!      end if
     
!   end if
  
  
! 1000 format('3F12.6')
  
! END PROGRAM test

