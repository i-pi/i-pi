!doublewell model
!Y. Litman 2019

      SUBROUTINE getdoublewell_1D(nat,q,pot,force)
        IMPLICIT NONE
        integer      :: nat, i
        real(kind=8) :: A,B
        real(kind=8) :: q(nat,3),pot,force(nat,3)
        real(kind=8) :: x,y,z,fx,fy,fz
        real(kind=8) :: k


        k  = 1836*(3800.0d0/219323d0)**2
        pot            =  0.0
        DO i = 1,nat 
           x = q(i,1)
           y = q(i,2)
           z = q(i,3)

           pot            =  0.0
           
           pot            = pot+ 0.5*k*y**2
           fy             = -k*y
           pot            = pot+ 0.5*k*z**2
           fz             = -k*z
           
           A = -0.00476705894242374 
           B = 0.0005980249683218661
           pot = pot + A * (x ** 2) + B * (x ** 4)

 
           fx = - ( 2 * A * (x) + 4 * B * (x ** 3) )

           force(i,1) = fx
           force(i,2) = fy
           force(i,3) = fz

        ENDDO  
                
        return
        end

