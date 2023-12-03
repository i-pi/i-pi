!Harmonic bath

!  V(q,x1,..,xn) = sum_1^(nat-1) [ 0.5 m omega_i^2(q - (c_i s(q))/(m omega_i)^2)^2 ]
!      s(q) = q *sd(q)
!      sd(q) = [1+eps exp( (q-0)^2 / (2deltaQ^2) ) ] + delta tanh(q/deltaQ)
!      If eps=delta=0 then sd(q) =1 and s(q) = q 


      SUBROUTINE get_harmonic_bath(nat,bath_type,friction,omega_c,eps,delta,deltaQ,q,pot,force)
        IMPLICIT NONE
        integer      :: nat, i 
        real(kind=8),PARAMETER ::  pi= 3.141592653589793D0 
        real(kind=8) :: bath_type,friction,omega_c,mass 
        real(kind=8) :: eps,delta,deltaQ,S,dSD_dq
        real(kind=8) :: A,B
        real(kind=8) :: q(nat*3),pot,force(nat*3),x(nat*3-1)
        real(kind=8) :: c(nat*3-1),omega(nat*3-1), omega2(nat*3-1),aux
        mass = 1837.36223469
        pot            =  0.0
        A = -0.00476705894242374
        B = 0.0005980249683218661

        DO i = 1, 3*nat -1 
          x(i) = q(i+1)
        ENDDO

        !Get omega and c_i
        IF (IDINT(bath_type) .EQ. 1 ) THEN !Ohmic
           DO i = 1, 3*nat -1 
              omega(i) = - omega_c * LOG(  (i - 0.5 ) /  (3*nat-1) )  
              omega2(i) = omega(i)**2
              c(i) = omega(i) * ( ( 2 * friction * mass *omega_c) / ( (3*nat-1) * pi ) )**0.5
           ENDDO 
        END IF
         
        !SYSTEM
         pot =  A * (q(1) ** 2) + B * (q(1) ** 4)
         force(1) = - ( 2 * A * (q(1)) + 4 * B * (q(1) ** 3) )

        !BATH
        DO i = 1,3*nat-1 
           aux =  x(i) - ( c(i)*S(q(1),eps,delta,deltaQ) / ( mass*omega2(i) ) )

           pot = pot + 0.5 * mass * omega2(i) * aux **2
           force(1)  = force(1) + aux * c(i)*dSD_dq(q(1),eps,delta,deltaQ)
           force(i+1) = -mass * omega2(i) * aux 


        ENDDO  

        return
        end

        REAL*8 FUNCTION   SD(q,eps,delta,deltaQ)
        IMPLICIT NONE
        real(kind=8), intent(in) :: q,eps,delta,deltaQ  ! q system position 
        real(kind=8)             :: dx       

           dx =  q / deltaQ
           SD = 1.0 + eps*DEXP(-0.5 * (dx**2) ) + delta*DTANH(dx)
           !write(6,*) 'ALBERTO', 1,eps, dx,  DEXP(-0.5 * (dx**2) ),  eps*DEXP(-0.5 * (dx**2) ) 
        END FUNCTION SD

        REAL*8 FUNCTION   S(q,eps,delta,deltaQ)
        IMPLICIT NONE
        real(kind=8), intent(in) :: q,eps,delta,deltaQ  
        real(kind=8)             :: SD 
           S = q* SD(q,eps,delta,deltaQ) 
        END FUNCTION S

        REAL*8 FUNCTION   dSD_dq(q,eps,delta,deltaQ)
        IMPLICIT NONE
        real(kind=8), intent(in) :: q,eps,delta,deltaQ  
        real(kind=8)             :: dx,dsddq1,dsddq2,SD 
           dx  =  q / deltaQ
           dsddq1 = eps*DEXP(-0.5 * (dx**2))*(-dx/deltaQ)
           dsddq2 = delta* (1- DTANH(DX)**2) / deltaQ
           dSD_dq  =  q *(dsddq1 + dsddq2) +  SD(q,eps,delta,deltaQ) 

        END FUNCTION dSD_dq

      SUBROUTINE get_meanfield_harmonic_bath(nat,eta,eps,delta,deltaQ,q,pot,force,friction)
        IMPLICIT NONE
        integer      :: nat, i 
        real(kind=8) :: A,B,k,eta
        real(kind=8) :: eps,delta,deltaQ,dSD_dq
        real(kind=8) :: x,y,z,fx,fy,fz
        real(kind=8) :: q(nat,3),pot,force(nat,3)
        real(kind=8) :: friction(3*nat,3*nat)

        k  = 1836*(3800.0d0/219323d0)**2

        A = -0.00476705894242374
        B = 0.0005980249683218661

        DO i = 1,nat
           x = q(i,1)
           y = q(i,2)
           z = q(i,3)

           pot            =  0.0

           pot            = pot+ 0.5*k*y**2
           fy             = -k*y
           pot            = pot+ 0.5*k*z**2
           fz             = -k*z

           pot = pot + A * (x ** 2) + B * (x ** 4)


           fx = - ( 2 * A * (x) + 4 * B * (x ** 3) )

           force(i,1) = fx
           force(i,2) = fy
           force(i,3) = fz

           friction(nat*3*(i-1)+1,nat*3*(i-1)+1) = eta * (dSD_dq(x,eps,delta,deltaQ) )**2  
           friction(nat*3*(i-1)+2,nat*3*(i-1)+2) = 0 !eta * (dSD_dq(y,eps,delta,deltaQ) )**2 
           friction(nat*3*(i-1)+3,nat*3*(i-1)+3) = 0 !eta * (dSD_dq(z,eps,delta,deltaQ) )**2 

        ENDDO

        return
        end

