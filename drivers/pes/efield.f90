!
!=============================================================================
!
! File: efield.f90
!
! NOTES:
! (1) ATOMIC UNITS ARE USED THROUGHOUT
! (2) WE ASSUME THAT THE ATOMS ARE IN A LIST ARRANGED WITH MOLECULES
!     AS [ {O,H,H},{O,H,H},..... ].
! (3) PBCs are NOT included - clusters only!
! (4) q-TIP4P/F parameters are hard-coded in the "potential" subroutine.
!
! Disclaimer: please thoroughly check the output of this before 
! using it for anything interesting!
!
!=============================================================================
!

subroutine efield_v(rt,na,dvdrt,v, vir,efield)
  implicit none  
  integer :: i,j,na,nm,ic,k
  real(8) :: r(3,na), dvdr(3,na), rt(na,3), dvdrt(na,3),efield(3)
  real(8) :: oo_eps, oo_sig, rcut,v
  real(8) :: apot, bpot, alp, alpha, alpha2
  real(8) :: qo,qh,theta,reoh,vir(3,3)
  real(8), allocatable :: ro(:,:), z(:)

  ! Set and check number of molecules.
  ! 
  nm = na/3
  if (3*nm .ne. na) stop 'ERROR 1 in POTENTIAL'

  
  !
  ! Set potential parameters - ATOMIC UNITS!
  ! NOTE: HERE WE HAVE THE PARAMETERS FOR q-TIP4P/F.
  !
  rcut = 9.d0 / 0.5291772108d0               ! 9 Angstrom cutoff for LJ
  alpha = 0.73612d0
  qo = -1.1128d0
  qh = +0.5d0 * 1.1128d0
  oo_sig = 5.96946d0
  oo_eps = 2.95147d-4
  theta = 107.4d0 * (dacos(-1.d0) / 180.d0)
  reoh = 1.78d0
  apot = 0.185d0
  bpot = 0.07d0
  alp = 1.21d0
  alpha2 = 0.5d0 * (1.d0 - alpha)


  ! Zero-out potential and derivatives.
  !
  v = 0.d0
  vir = 0.d0
  dvdr(:,:) = 0.d0

  do i=1, na 
    r(:,i)=rt(i,:)
  enddo


  !
  ! *** COULOMB CALCULATION ***
  !
  ! Allocate space for oxygen position storage.
  !
  allocate (ro(3,nm), z(na) )


  ! Determine the positions of the m-sites.
  !
  ic = 0
  do i = 1, na, 3
     ic = ic + 1
     do j = 1, 3
        ro(j,ic) = r(j,i)
        r(j,i) = alpha * r(j,i) + alpha2*(r(j,i+1)+r(j,i+2))
     enddo
  enddo


  ! Allocate atomic charges and molecular identities.
  !
  ic = 0
  do i = 1, na, 3
     z(i) = qo
     z(i+1) = qh
     z(i+2) = qh
   enddo

  ! Calculate forces, potential and virial

  do i=1,na
   ! Potential
   v = v - z(i)*dot_product(r(:,i),efield(:))
   ! Force
   dvdr(:,i) = z(i)*efield(:)
   ! Virial
   do j=1,3
    do k=1,3
     vir(j,k) = vir(j,k) - r(j,i)*dvdr(k,i)
    enddo
   enddo
  enddo
 
  !
  ! Use the chain rule to calculate the correct forces on the atoms. 
  ! Also, replace the m-site with the original oxygen atoms.
  !
  ic = 0
  do j = 1, na, 3
     ic = ic + 1
     do i = 1, 3
        dvdr(i,j+1) = dvdr(i,j+1) + alpha2 * dvdr(i,j) 
        dvdr(i,j+2) = dvdr(i,j+2) + alpha2 * dvdr(i,j) 
        dvdr(i,j) = alpha * dvdr(i,j) 
        r(i,j) = ro(i,ic)
     enddo
  enddo  
  deallocate (ro, z)
  vir = -1.0d0 *vir
  do i=1, na 
    dvdrt(i,:) = -dvdr(:,i)
  enddo
  return
end Subroutine efield_v
