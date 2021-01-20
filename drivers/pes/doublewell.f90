!doublewell model
!Y. Litman 2019

      SUBROUTINE getdoublewell(nat,q,pot,force)
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

           A = -0.00476705894242374 
           B = 0.0005980249683218661
           pot = pot + A * (x ** 2) + B * (x ** 4)
           pot = pot + A * (y ** 2) + B * (y ** 4)
           pot = pot + A * (z ** 2) + B * (z ** 4)
 
           fx = - ( 2 * A * (x) + 4 * B * (x ** 3) )
           fy = - ( 2 * A * (y) + 4 * B * (y ** 3) )
           fz = - ( 2 * A * (z) + 4 * B * (z ** 3) )

           force(i,1) = fx
           force(i,2) = fy
           force(i,3) = fz

        ENDDO  
                
        return
        end


    SUBROUTINE dw_friction(nat,q,friction)
        IMPLICIT NONE
        integer      :: nat, i
        real(kind=8) :: q(nat,3),friction(nat,6)
        real(kind=8) :: x,y,z


        DO i = 1,nat
            x = q(i,1)
            y = q(i,2)
            z = q(i,3)

            friction(i,1) = x**2
            friction(i,2) = y**2
            friction(i,3) = z**2
            friction(i,4) = x*y
            friction(i,5) = x*z
            friction(i,6) = y*z
        ENDDO

        return
    end
