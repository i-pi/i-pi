!Compute the eckart potential for the first dimension of the first atoms. 
!For the other 2 dimensions and the rest of the atoms compute an harmonic potential
!with frequency 3800 cm^-1 (for H)

!Eckart potential
!           AA                     BB
!V(x)=  -------------     +   --------------
!        1+exp(-2x/A)          cosh^2(x/A)
!Default values
!
!AA  = 0.0d0
!A   = 0.66047
!BB  = (6*12)/( 1836 * (A**2) *( (4.D0 * ATAN(1.0d0) )**2 ) )        ---->9.10864e-3
!k   = 1836*(3800.0d0/219323d0)**2                                   ---->3800 cm^-1 with proton mass


      SUBROUTINE geteckart(nat,AA,A,BB,k,q,pot,force)
        IMPLICIT NONE
        integer      :: nat, i
        real(kind=8) :: AA,A,BB,k
        real(kind=8) :: q(nat,3),pot,force(nat,3)
        real(kind=8) :: x,y,z,fx,fy,fz
        real(kind=8) :: V1,V2,EC

        DO i = 1,nat 
           x = q(i,1)
           y = q(i,2)
           z = q(i,3)

           !---------------------------------------
           !Eckart---------------------------------
           IF (i == 1) THEN
               V1         = AA / ( 1 + EXP(-2.0d0*x/A) )
               V2         = BB / ( cosh( x/A )**2 )
               pot        = V1+V2

               EC         = EXP(-2.0d0*X/A)
               V1         = AA*EC / ( (1+EC)**2 )
               V2         = -BB*sinh(x/A) / ( cosh( x/A )**3 )
               fx         = -(V1+V2)*(2.0d0/A)
           ELSE
               pot        = pot + 0.5*k*x**2
               fx         = -k*x
           ENDIF
           !Harmonic
           pot            = pot + 0.5*k*y**2
           fy             = -k*y
           pot            = pot + 0.5*k*z**2
           fz             = -k*z
           !---------------------------------------

           force(i,1) = fx
           force(i,2) = fy
           force(i,3) = fz

        ENDDO       
        return
        end





