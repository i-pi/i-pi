+subroutine UpdateRDF(gOOr, atxyz1, atxyz2, nat1, nat2, nbins, r_min, r_max, h, ainv, nbeads, mass1, mass2)!, forma)
+!
+IMPLICIT NONE
+!
+! Input
+!
+INTEGER, INTENT(in) :: nbeads, nbins, nat1, nat2
+REAL(8), INTENT(in) :: mass1, mass2
+REAL(8), INTENT(in) :: r_min, r_max
+REAL(8), INTENT(in), DIMENSION(:,:) :: atxyz1, atxyz2, h, ainv
+!
+! Output
+!
+!INTEGER, INTENT(inout), DIMENSION(:,:)  :: gOOr
+REAL(8), INTENT(inout), DIMENSION(:,:)  :: gOOr
+!
+! Local variables
+!
+INTEGER :: ia, ib, ig, ih
+REAL(8) :: deltar, dAB, vdAB(3), temp, norm 
+REAL(8), PARAMETER :: tol = 0.00001
+!REAL(8), DIMENSION(nbins,nbins),SAVE  :: forma
+!
+! Histogram step initialization
+!
+deltar = gOOr(2,1) - gOOr(1,1)
+!
+! Start computing g(r) from MD configurations
+!
+! Normalization constant
+!
+!IF (mass1.EQ.mass2) THEN
+!  norm = 2.d0/(nat1*(nat2-1))
+!ELSE
+!  norm = 1.d0/(nat1*nat2)
+!END IF
+!
+! Populate histogram bins for gOO(r)...
+!
+IF (mass1.EQ.mass2) THEN
+  Do ih=1,nbeads
+    DO ia=1,nat1
+      DO ib=ia,nat2

