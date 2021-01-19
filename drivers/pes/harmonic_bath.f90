!Harmonic bath
!JCP 122, 084106(2005)

!  V(q,x1,..,xn) = sum_1^(nat-1) [ 0.5 m omega_i^2(q - (c_i q)/m omega_i)^2 ]

      SUBROUTINE get_harmonic_bath(nat,bath_type,friction,omega_c,q,pot,force)
        IMPLICIT NONE
        integer      :: nat, i 
        real(kind=8),PARAMETER ::  pi= 3.141592653589793D0 
        real(kind=8) :: bath_type,friction,omega_c,mass 
        real(kind=8) :: A,B
        real(kind=8) :: q(nat*3),pot,force(nat*3),x(nat*3-1)
        real(kind=8) :: c(nat*3-1),omega(nat*3-1), omega2(nat*3-1),aux
        mass = 1837
        pot            =  0.0
        A = -0.00476705894242374
        B = 0.0005980249683218661

        DO i = 1, 3*nat -1 
          x(i) = q(i+1)
        ENDDO
        IF (IDINT(bath_type) .EQ. 1 ) THEN !Ohmic
           DO i = 1, 3*nat -1 
              omega(i) = - omega_c * LOG(  (i - 0.5 ) /  (3*nat-1) )  
              omega2(i) = omega(i)**2
              c(i) = omega(i) * ( ( 2 * friction * mass *omega_c) / ( (3*nat-1) * pi ) )**0.5
              !write(6,*)'Ohmic',i,omega(i),c(i)           
           ENDDO 
        END IF

        !SYSTEM
         pot =  A * (q(1) ** 2) + B * (q(1) ** 4)
         force(1) = - ( 2 * A * (q(1)) + 4 * B * (q(1) ** 3) )

        !BATH
        DO i = 1,3*nat-1 
           aux =  x(i) - c(i)*q(1) / ( mass*omega2(i) ) 

           pot = pot + 0.5 * mass * omega2(i) * aux **2
           
           force(1)  = force(1) + aux * c(i)
           force(i+1) = -mass * omega2(i) * aux 
           !WRITE(6,*)i,mass,omega2(i),aux
           !WRITE(6,*)i ,x(i), c(i),q(1) , mass, omega2(i) 
        ENDDO  

        DO i = 1,3*nat
        WRITE(6,*) i, q(i), force(i)
        !WRITE(6,*) i, force(i),x(i) - (c(i)*q(1)) / ( mass*omega2(i) ),x(i+1) - (c(i+1)*q(1)) / ( mass*omega2(i+1) )
        ENDDO
        return
        end


      SUBROUTINE get_meanfield_harmonic_bath(nat,eta,q,pot,force,friction)
        IMPLICIT NONE
        integer      :: nat, i 
        real(kind=8),PARAMETER ::  pi= 3.141592653589793D0 
        real(kind=8) :: A,B,k,eta
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
           friction = eta 

        ENDDO

        return
        end

