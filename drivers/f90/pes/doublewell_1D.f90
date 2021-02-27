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

    SUBROUTINE dw1d_dipole(nat,q,dip)
        IMPLICIT NONE
        integer                   :: nat, i
        real(kind=8),intent(inout):: q(nat,3)
        real(kind=8),intent(out)  :: dip(3)
        real(kind=8)              :: x,y,z

        dip = 0.0d0
        DO i = 1,nat
            x = q(i,1)
            y = q(i,2)
            z = q(i,3)

            dip(1) = dip(1)+x**3
            dip(2) = dip(2)-y**3
            dip(3) = dip(3)+z**3
        ENDDO

        return
    END

    SUBROUTINE dw1d_friction(nat,q,friction)
        IMPLICIT NONE
        integer                  :: nat, i, j, qi, qj, coori, coorj
        real(kind=8),intent(in)  :: q(nat,3)
        real(kind=8),intent(out) :: friction(3*nat,3*nat)

        DO i = 1,nat
            DO j = 1,nat
                DO qi = 1,3
                    coori = (i-1)*nat+qi
                    DO qj = 1,3
                        coorj = (j-1)*nat+qj
                        friction(coori,coorj) = q(i,qi)*q(j,qj)
                    END DO
                END DO
            END DO
        END DO

        return
    END
