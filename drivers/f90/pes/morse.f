
      SUBROUTINE getmorse(r0,D,a,q,pot,force)
         IMPLICIT NONE
         DOUBLE PRECISION r0, D, a, q(1,3), pot, force(1,3)
         DOUBLE PRECISION nq, dr
         nq=dsqrt(sum(q*q))
         dr=nq-r0
         pot=D*(dexp(-2.*a*dr)-2.*dexp(-a*dr))
         force=-2*a*D*(dexp(-a*dr)-dexp(-2.0*a*dr))*q/nq
      END SUBROUTINE

      SUBROUTINE getmorsedia(r0,D,a,q,pot,force)
         IMPLICIT NONE
         DOUBLE PRECISION r0, D, a, q(2,3), pot, force(2,3)
         DOUBLE PRECISION nq, dr, diff(3)
         diff(:)=(q(1,:)-q(2,:))
         nq=dsqrt(dot_product(diff,diff))
         dr=nq-r0
         pot=D*(dexp(-2.*a*dr)-2.*dexp(-a*dr))
         force(1,:)=-2*a*D*(dexp(-a*dr)-dexp(-2.0*a*dr))*diff/nq
         force(2,:)=2*a*D*(dexp(-a*dr)-dexp(-2.0*a*dr))*diff/nq
      END SUBROUTINE
