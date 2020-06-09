!Muller Brown surface
!Theoret. Chim. Acta (Berl.) 53, 75-93 (1979)
!J. Chem. Theory Comput. 2018, 14, 5489−5498 Article
!Y. Litman 2019

      SUBROUTINE get_MB(nat,scaling,q,pot,force)
      IMPLICIT NONE
      integer      :: nat,i  
      real(kind=8) :: scaling
      real(kind=8) :: q(nat,3),pot,force(nat,3)
      real(kind=8) :: x,y,z,fx,fy,fz
      real(kind=8) :: k

      real(kind=8),dimension(4) :: A0,a,b,c,x0,y0,V,gx,gy

      scaling = 0.004737803248674678
      k  = 1836*(3800.0d0/219323d0)**2

      A0  = (/ -200.0 , -100.0   , -170.0 , 15.0  /)
      a  = (/  -1.0   ,   -1.0   ,   -6.5 ,  0.7  /)
      b  = (/   0.0   ,    0.0   ,   11.0 ,  0.6  /)
      c  = (/ -10.0   ,  -10.0   ,   -6.5 ,  0.7  /)
      x0 = (/   1.0   ,    0.0   ,   -0.5 , -1.0  /)
      y0 = (/   0.0   ,    0.5   ,    1.5 ,  1.0  /)
     
      A0 = A0*scaling 
      DO i = 1,nat
         x = q(i,1)
         y = q(i,2)
         z = q(i,3)

         !z coordinate
         pot            =  0.5*k*z**2
         fz             = -k*z

         !xy coordinate
         V(1) = A0(1)*DEXP( a(1)*(x - x0(1))**2 + b(1)*( x - x0(1) )*( y - y0(1) ) +  c(1)* (y - y0(1))**2  )
         V(2) = A0(2)*DEXP( a(2)*(x - x0(2))**2 + b(2)*( x - x0(2) )*( y - y0(2) ) +  c(2)* (y - y0(2))**2  )
         V(3) = A0(3)*DEXP( a(3)*(x - x0(3))**2 + b(3)*( x - x0(3) )*( y - y0(3) ) +  c(3)* (y - y0(3))**2  )
         V(4) = A0(4)*DEXP( a(4)*(x - x0(4))**2 + b(4)*( x - x0(4) )*( y - y0(4) ) +  c(4)* (y - y0(4))**2  )

         gx(1)  = V(1)* ( 2*a(1)*(x - x0(1)) +  b(1)*( y - y0(1) ) )
         gx(2)  = V(2)* ( 2*a(2)*(x - x0(2)) +  b(2)*( y - y0(2) ) )
         gx(3)  = V(3)* ( 2*a(3)*(x - x0(3)) +  b(3)*( y - y0(3) ) )
         gx(4)  = V(4)* ( 2*a(4)*(x - x0(4)) +  b(4)*( y - y0(4) ) )


         gy(1)  = V(1)* ( b(1)*( x - x0(1) ) + 2 * c(1) * (y - y0(1)) )
         gy(2)  = V(2)* ( b(2)*( x - x0(2) ) + 2 * c(2) * (y - y0(2)) )
         gy(3)  = V(3)* ( b(3)*( x - x0(3) ) + 2 * c(3) * (y - y0(3)) )
         gy(4)  = V(4)* ( b(4)*( x - x0(4) ) + 2 * c(4) * (y - y0(4)) )


         pot = pot + V(1) + V(2) + V(3) + V(4)
         fx  = -gx(1) - gx(2) - gx(3) - gx(4) 
         fy  = -gy(1) - gy(2) - gy(3) - gy(4)
     
         force(i,1) = fx
         force(i,2) = fy
         force(i,3) = fz
        
      ENDDO
      ! write(6,*) 'ALBERTO'
      return
      end
 
