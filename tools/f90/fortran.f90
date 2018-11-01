subroutine F2divM(f2, fatxyz, masses, nat, nob)
!
IMPLICIT NONE
!
! Input
!
REAL(8), INTENT(in), DIMENSION(:,:) :: fatxyz    ! forces
REAL(8), INTENT(in), DIMENSION(:)   :: masses    ! masses
INTEGER, INTENT(in)                 :: nat, nob  ! number of atoms and number of beads
!
! Output
!
REAL(8), INTENT(out)              :: f2
!
! Local variables
!
INTEGER                           :: ia, ib
!
! Square forces divided by mass calculations
!
f2 = 0
DO ib=1,nob
  DO ia=1,nat
    f2 = f2 + (fatxyz(ib, 3*ia-2)*fatxyz(ib, 3*ia-2) + fatxyz(ib, 3*ia-1)*fatxyz(ib, 3*ia-1) + fatxyz(ib, 3*ia)*fatxyz(ib, 3*ia))/masses(ia)
  END DO ! ia
END DO !ib
!
END subroutine F2divM
!
!
!
!
subroutine findCentroidVirialKineticEnergy(Kcv, atxyz, fatxyz, nat, nob)
!
IMPLICIT NONE
!
! Input
!
INTEGER, INTENT(in)                  :: nat, nob                    ! number of atoms and number of beads
REAL(8), INTENT(in), DIMENSION(:,:)  :: atxyz, fatxyz               ! coordinates and forces
!
! Output
!
REAL(8), INTENT(out)                 :: Kcv
!
! Local variables
!
INTEGER                              :: ia, ib, ic
REAL(8)                              :: temp, xc(3)
!
!
Kcv = 0.0
DO ia=1,nat
  DO ic=1,3
    xc(ic) = 0.0
    DO ib=1,nob
      xc(ic) = xc(ic) + atxyz(ib, ic + 3*ia - 3)
    END DO !ib
    xc(ic) = xc(ic)/nob
  END DO !ic
  temp = 0.0
  DO ic=1,3
    DO ib=1,nob
      temp = temp + (xc(ic) - atxyz(ib, ic + 3*ia - 3))*fatxyz(ib, ic + 3*ia - 3)
    END DO !ib
  END DO !ic
  Kcv = Kcv + temp
END DO !ia
Kcv = Kcv/(2.0*nob)
!
END SUBROUTINE findCentroidVirialKineticEnergy
!
!
!
!
subroutine findCoupling(coupling, atxyz, mass, temperature, nat, nob)
!
IMPLICIT NONE
!
! Input
!
INTEGER, INTENT(in)                 :: nat, nob      ! number of atoms and number of beads
REAL(8), INTENT(in)                 :: temperature   ! temperature
REAL(8), INTENT(in), DIMENSION(:)   :: mass          ! masses
REAL(8), INTENT(in), DIMENSION(:,:) :: atxyz         ! coordinates
!
! Output
!
REAL(8), INTENT(out)                :: coupling
!
! Local variables
!
INTEGER :: ia, ib
REAL(8) :: temp
!
! Computing coupling term
!
coupling = 0.0
DO ia=1,nat
  temp = 0.0
  DO ib=1,nob-1
    temp = temp + (atxyz(ib, 3*ia - 2) - atxyz(ib+1, 3*ia - 2))**2 + (atxyz(ib, 3*ia - 1) - atxyz(ib+1, 3*ia - 1))**2 + (atxyz(ib, 3*ia) - atxyz(ib+1, 3*ia))**2
  END DO !ib
  temp = temp + (atxyz(nob, 3*ia - 2) - atxyz(1, 3*ia - 2))**2 + (atxyz(nob, 3*ia - 1) - atxyz(1, 3*ia - 1))**2 + (atxyz(nob, 3*ia) - atxyz(1, 3*ia))**2
  coupling = coupling + temp*mass(ia)
END DO !ia
coupling = -coupling*0.5*nob*temperature**2
!
END SUBROUTINE findCoupling
!
!
!
!
subroutine UpdateQRDF(gOOr, f2gOOr, fgOOr, atxyz1, atxyz2, fatxyz1, fatxyz2, nat1, nat2, nbins, r_min, r_max, h, ainv, nbeads, f2, mass1, mass2)
!
IMPLICIT NONE
!
! Input
!
INTEGER, INTENT(in) :: nbeads, nbins, nat1, nat2
REAL(8), INTENT(in) :: mass1, mass2
REAL(8), INTENT(in) :: r_min, r_max, f2
REAL(8), INTENT(in), DIMENSION(:,:) :: atxyz1, atxyz2, fatxyz1, fatxyz2, h, ainv
!
! Output
!
REAL(8), INTENT(inout), DIMENSION(:)    :: f2gOOr, fgOOr
REAL(8), INTENT(inout), DIMENSION(:,:)  :: gOOr
!
! Local variables
!
INTEGER :: ia, ib, ig, ih
REAL(8) :: deltar, dAB, vdAB(3), temp, norm
REAL(8), PARAMETER :: tol = 0.00001
!
! Histogram step initialization
!
deltar = gOOr(2,1) - gOOr(1,1)
!
! Start computing g(r) from MD configurations
!
! Normalization constant
!
IF (mass1.EQ.mass2) THEN
  norm = 1.0/(nat1*(nat2-1))
ELSE
  norm = 1.0/(nat1*nat2)
END IF
!
! Populate histogram bins for gOO(r)...
!
Do ih=1,nbeads
  DO ia=1,nat1
    DO ib=1,nat2
      ! Compute the distance of the closest image of atom B to atom A using minimum image convention...
      CALL CalcMinDist(h,ainv,atxyz1(ih,3*ia-2),atxyz1(ih,3*ia-1),atxyz1(ih,3*ia),atxyz2(ih,3*ib-2),atxyz2(ih,3*ib-1),atxyz2(ih,3*ib),dAB,vdAB)
      ! Screen distances that are outside desired range
      IF (dAB.LT.r_max.AND.dAB.GT.r_min) THEN
        ig=INT((dAB-r_min)/deltar)+1  !bin/histogram position
        gOOr(ig,2)=gOOr(ig,2)+1*norm
        f2gOOr(ig)=f2gOOr(ig)+f2*norm
        !
        ! calc PPI corrections
        IF (dAB > tol) THEN
          temp = ((fatxyz1(ih,3*ia-2)/mass1 - fatxyz2(ih,3*ib-2)/mass2)*vdAB(1) + (fatxyz1(ih,3*ia-1)/mass1 - fatxyz2(ih,3*ib-1)/mass2)*vdAB(2) + (fatxyz1(ih,3*ia)/mass1 - fatxyz2(ih,3*ib)/mass2)*vdAB(3))/(dAB*deltar)
          IF(ig.GT.1) fgOOr(ig-1) = fgOOr(ig-1) - temp*norm
          IF(ig.LT.nbins) fgOOr(ig+1) = fgOOr(ig+1) + temp*norm
        END IF
      END IF
    END DO !ib
  END DO !ia
END DO !ih
!
END SUBROUTINE UpdateQRDF
!
!
!
!
subroutine CalcMinDist(h,ainv,xA,yA,zA,xB,yB,zB,dAB,rAB)
!
IMPLICIT NONE
!
REAL(8) :: xA,yA,zA,xB,yB,zB,dAB
REAL(8) :: rAB(3),rAB2(3),ainv(3,3),h(3,3)
!Initialization of distance
dAB=0.0
!
! Compute distance between atom A and atom B (according to the minimum
! image convention)...
!
rAB(1)=xA-xB   ! r_AB = r_A - r_B
rAB(2)=yA-yB   ! r_AB = r_A - r_B
rAB(3)=zA-zB   ! r_AB = r_A - r_B
!
rAB2(1)=ainv(1,1)*rAB(1)+ainv(1,2)*rAB(2)+ainv(1,3)*rAB(3)   ! s_AB =h^-1 r_AB
rAB2(2)=ainv(2,1)*rAB(1)+ainv(2,2)*rAB(2)+ainv(2,3)*rAB(3)   ! s_AB =h^-1 r_AB
rAB2(3)=ainv(3,1)*rAB(1)+ainv(3,2)*rAB(2)+ainv(3,3)*rAB(3)   ! s_AB =h^-1 r_AB
!
rAB2(1)=rAB2(1)-IDNINT(rAB2(1))   ! impose MIC on s_AB in range:[-0.5,+0.5]
rAB2(2)=rAB2(2)-IDNINT(rAB2(2))   ! impose MIC on s_AB in range:[-0.5,+0.5]
rAB2(3)=rAB2(3)-IDNINT(rAB2(3))   ! impose MIC on s_AB in range:[-0.5,+0.5]
!
rAB(1)=h(1,1)*rAB2(1)+h(1,2)*rAB2(2)+h(1,3)*rAB2(3)   ! r_AB = h s_AB(MIC)
rAB(2)=h(2,1)*rAB2(1)+h(2,2)*rAB2(2)+h(2,3)*rAB2(3)   ! r_AB = h s_AB(MIC)
rAB(3)=h(3,1)*rAB2(1)+h(3,2)*rAB2(2)+h(3,3)*rAB2(3)   ! r_AB = h s_AB(MIC)
!
dAB=DSQRT(rAB(1)*rAB(1)+rAB(2)*rAB(2)+rAB(3)*rAB(3))   ! |r_A -r_B| (MIC)
!
END SUBROUTINE CalcMinDist