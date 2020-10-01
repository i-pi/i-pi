!********************************************************************************
!------------Note by Xinchuan Huang on 2004-02-16----------------
! Potential and dipole moment surface for H5O2+
! version 4B, capable of H2O...H3O+ dissociation
!
! Citation for this PES&DMS:
!  Xinchuan Huang, Bastiaan J. Braams, and Joel M. Bowman, J. Chem. Phys. 122, 044308 (2005)
!
!  Fitting codes were written by Bastiaan J. Bramms on 2004-01-03, 2004-01-04 and 2004-01-09
!  Codes were modified and applied to fit by Xinchuan Huang    ---- Finished on 2004-01-14
!
!  pgf90 & ifc/ifort compatible
!  Two data file are required :  
!    h5o2.dms4B.coeff.com.dat  h5o2.pes4B.coeff.dat
!
!  The C2-symmetry minimum geometry on PES-4B in the order of O O H H H H H is: (in bohr)
!          X                  Y                 Z 
! O  2.254868053895938 0.000000000000000 0.000000000000000
! O -2.254868053895938 0.000000000000000 0.000000000000000
! H  0.000000000000000 0.000000000000000 0.1235103451781876
! H -3.013194720920241 1.489183788836780 -0.7465980020137635
! H -3.201587564645578 -0.4544810799466731 1.498136372134866
! H  3.013194720920241 -1.489183788836780 -0.7465980020137635
! H  3.201587564645578 0.4544810799466731 1.498136372134866
!  
!  prepot() should be called before first calling calcpot(V,xx)
!    xx(7,3) is cartesian array containing O O H H H H H coor, in bohr
!    returned V is potential in hartree, taking C2-sym minimum as zero potential reference
!
!  predip() should be calld before first calling calcdip(dip,xx)
!    xx(7,3) is cartesian array containing O O H H H H H coor, in bohr
!    returned dip(1),dip(2),dip(3) are the x-, y-, z- dipole components
!              relative to the center-of-mass of H5O2+, in a.u.
!
!  For D5O2+ calculations, you may need to change calcdip() to user D5O2+ center-of-mass frame
!
!
!  Permanent contact:  Joel M. Bowman    E-mail: jmbowma@emory.edu 
!  (Technical questions may send to : xinchuan@hotmail.com)
!------------End of Note----------------------------------------------
!********************************************************************************

        subroutine prezundelpot()

        implicit double precision (a-h,o-z)
        implicit integer (i-n)

        double precision dc0(0:6,0:6),dw0(0:6,0:6), coef(0:7961)

        integer i1

        common/NCOE/ms,mr
        common/h5o2coef/dc0,dw0,coef

        ms=7785 ; mr=59

        open(20,file='h5o2.pes4B.coeff.dat',status='old')
        read(20,*)
        read(20,*)dc0
        read(20,*)
        read(20,*)dw0
        read(20,*)
        read(20,*)
        read(20,*)(coef(i1),i1=0,ms+3*mr-1)
        close(20)

        return
        end

!***************************************************************

        subroutine zundelpot(V,cart_in)

        implicit none

        integer ms,mr
        double precision dc0(0:6,0:6),dw0(0:6,0:6),coef(0:7961)

        common/NCOE/ms,mr
        common/h5o2coef/dc0,dw0,coef
        double precision V, cart_in(7,3)

        integer j


        double precision vec(0:ms+3*mr-1)
        double precision xnuc(0:2,0:6)

         do j=1,3
          xnuc(j-1,5:6)=cart_in(1:2,j)
          xnuc(j-1,0:4)=cart_in(3:7,j)
        end do

        call getvec (ms, mr, xnuc(0:2,0:6), dc0, dw0, vec)

        V = dot_product(coef,vec)

! set C1 min as potential zero-point   ---- PES-4B
        V = V + 153.012245695813d0

! set C2 linear SP as potential zero-point  --- PES-4B  (OH not equal, 7 variable)
!        V = V + 153.0119260317857d0

        return
        end subroutine
!********************************************************
        subroutine getvec (ms, mr, xn, dc0, dw0, vec)                                
        implicit none                                                                
        integer, parameter :: wp=selected_real_kind(12,300)                          
! version for X5Y2                                                           
        integer nk, ms, mr                                                           
        parameter (nk=7)                                                             
        double precision  xn(0:2,0:nk-1), dc0(0:nk-1,0:nk-1),
     $  dw0(0:nk-1,0:nk-1), vec(0:ms+3*mr-1)                                                           
        integer k, l                                                                 
        double precision  rvec(0:3),d0(0:nk-1,0:nk-1),r0(0:nk-1,0:nk-1)               
!-----------------------------------------------------------------------     
        vec=0                                                                      
        call getr0 (nk, xn, r0)                                                      
        call getd0 (nk, r0, dc0, dw0, d0)                                            
        call getv_x5y2 (ms, d0, vec(0:ms-1))                                         
        call getrvec (4, r0, rvec)                                                   
        do l=0, mr-1                                                               
         do k=0, 2                                                                 
          vec(ms+3*l+k)=rvec(k+1)*vec(l)                                           
         enddo                                                                       
        enddo                                                                        
        return                                                                       
        end                                                                          
        subroutine getv_x5y2 (m, d, vec)                                             
        implicit none                                                                
        integer, parameter :: wp=selected_real_kind(12,300)                          
        integer, parameter :: nk=7, mxd=12                                           
        integer m                                                                    
        double precision  d(0:nk-1,0:nk-1), vec(0:m-1)                                  
! version for molecule X5Y2.                                                 
! MolienSeries(0:*): 1 3 12 43 159 554 1879 6066 18755 ...                   
! #Primaries(1:*):   3 5 3 4 3 2 0 0 0 1                                     
! #Secondaries(1:*): 0 1 12 39 113 338 932 2402 ...                          
        integer, parameter :: l0=1, l1=l0+3, l2=l1+12, l3=l2+43,
     $ l4=l3+159,l5=l4+554, l6=l5+1879, l7=l6+6066-932                                      
!! We haven't incorporated the secondaries at degree 7                       
        integer, parameter :: np1=3, np2=5, np3=3, np4=4, np5=3, np6=2, 
     $   np7=0                                                                      
        integer, parameter :: ns1=0, ns2=1, ns3=12, ns4=39, ns5=113,                
     $     ns6=338, ns7=932                                                           
        double precision  x(0:np1-1), y(0:np2-1),z(0:np3-1),u(0:np4-1),              
     $     v(0:np5-1), w(0:np6-1)                                        
        double precision  x2(0:np1-1), x3(0:np1-1), x4(0:np1-1), 
     $ x5(0:np1-1),x6(0:np1-1), x7(0:np1-1), y2(0:np2-1), y3(0:np2-1), 
     $ z2(0:np2-1), ys(0:ns2-1), zs(0:ns3-1), us(0:ns4-1),         
     $     vs(0:ns5-1), ws(0:ns6-1), w7s(0:ns7-1)                                     
        integer mdeg, i, j, i0, i1, i2, i3, i4, j0, j1                         
        double precision  d2(0:nk-1,0:nk-1), d3(0:nk-1,0:nk-1), 
     $  d4(0:nk-1,0:nk-1), d5(0:nk-1,0:nk-1), d6(0:nk-1,0:nk-1),
     $ d7(0:nk-1,0:nk-1)                    
        double precision  t0, di0j1                                                            
        double precision  her2, her3, her4, her5, her6, her7                            
!        her2(t0)=(4*t0**2-2)/dsqrt(dble(8*2))                               
!        her3(t0)=(8*t0**2-12)*t0/dsqrt(dble(16*6))                          
!        her4(t0)=((16*t0**2-48)*t0**2+12)/dsqrt(dble(32*24))                
!        her5(t0)=((32*t0**2-160)*t0**2+120)*t0/dsqrt(dble(64*120))          
!        her6(t0)=(((64*t0**2-480)*t0**2+720)*t0**2-120)/                          
!     $     dsqrt(dble(128*720))                                                
!        her7(t0)=(((128*t0**2-1344)*t0**2+3360)*t0**2-1680)*t0/                   
!     $     dsqrt(dble(256*5040))                                               
!-----------------------------------------------------------------------     
! Test for compatibility, set mdeg                                           
        select case (m)                                                              
        case (l0)                                                                    
         mdeg=0                                                                    
        case (l1)                                                                    
         mdeg=1                                                                    
        case (l2)                                                                    
         mdeg=2                                                                    
        case (l3)                                                                    
         mdeg=3                                                                    
        case (l4)                                                                    
         mdeg=4                                                                    
        case (l5)                                                                    
         mdeg=5                                                                    
        case (l6)                                                                    
         mdeg=6                                                                    
        case (l7)                                                                    
         mdeg=7                                                                    
        case default                                                                 
         stop 'getv - wrong dimension'                                               
        endselect                                                                    
! auxiliary distances                                                        
        do i=0, nk-1                                                               
         do j=i+1, nk-1                                                            
          d2(i,j)=her2(d(i,j))                                                     
          d2(j,i)=d2(i,j)                                                          
          d3(i,j)=her3(d(i,j))                                                     
          d3(j,i)=d3(i,j)                                                          
          d4(i,j)=her4(d(i,j))                                                     
          d4(j,i)=d4(i,j)                                                          
          d5(i,j)=her5(d(i,j))                                                     
          d5(j,i)=d5(i,j)                                                          
          d6(i,j)=her6(d(i,j))                                                     
          d6(j,i)=d6(i,j)                                                          
          d7(i,j)=her7(d(i,j))                                                     
          d7(j,i)=d7(i,j)                                                          
         enddo                                                                       
        enddo                                                                        
! Primary Invariants                                                         
        x=0 ; y=0 ; z=0 ; u=0 ; v=0 ; w=0                                
        do i0=0, 4                                                                 
         t0=0                                                                      
         do i1=0, 4                                                                
         if (i1.ne.i0) then                                                          
          t0=t0+d(i0,i1)/4                                                         
          x(0)=x(0)+d(i0,i1)/20                                                    
          y(0)=y(0)+d2(i0,i1)/20                                                   
          z(0)=z(0)+d3(i0,i1)/20                                                   
          u(0)=u(0)+d4(i0,i1)/20                                                   
          v(0)=v(0)+d5(i0,i1)/20                                                   
          w(0)=w(0)+d6(i0,i1)/20                                                   
         endif                                                                       
         enddo                                                                       
         y(1)=y(1)+her2(t0)/5                                                      
         z(1)=z(1)+her3(t0)/5                                                      
         u(1)=u(1)+her4(t0)/5                                                      
         v(1)=v(1)+her5(t0)/5                                                      
        enddo                                                                        
        x(1)=sum(d(0:4,5:6))/10                                                    
        do i0=0, 4                                                                 
         y(2)=y(2)+her2(sum(d(i0,5:6))/2)/5                                        
         z(2)=z(2)+her3(sum(d(i0,5:6))/2)/5                                        
         u(2)=u(2)+her4(sum(d(i0,5:6))/2)/5                                        
         v(2)=v(2)+her5(sum(d(i0,5:6))/2)/5                                        
        enddo                                                                        
        do j0=5, 6                                                                 
         y(3)=y(3)+her2(sum(d(0:4,j0))/5)/2                                        
        enddo                                                                        
        do i0=0, 4                                                                 
         y(4)=y(4)+sum(d2(i0,5:6))/10                                              
         u(3)=u(3)+sum(d4(i0,5:6))/10                                              
         w(1)=w(1)+sum(d6(i0,5:6))/10                                              
        enddo                                                                        
        x(2)=d(5,6)                                                                
! Required powers                                                            
        do i=0, np1-1                                                              
         x2(i)=her2(x(i))                                                          
         x3(i)=her3(x(i))                                                          
         x4(i)=her4(x(i))                                                          
         x5(i)=her5(x(i))                                                          
         x6(i)=her6(x(i))                                                          
         x7(i)=her7(x(i))                                                          
        enddo                                                                        
        do i=0, np2-1                                                              
         y2(i)=her2(y(i))                                                          
         y3(i)=her3(y(i))                                                          
        enddo                                                                        
        do i=0, np3-1                                                              
         z2(i)=her2(z(i))                                                          
        enddo                                                                        
! Secondary Invariants                                                       
!! ys(0:0), zs(0:11), us(0:38), vs(0:112), ws(0:337), w7s(0:931)             
!! reducible: us(0:0), vs(0:11), ws(0:113), w7s(0:541)                       
        ys=0 ; zs=0 ; us=0 ; vs=0 ; ws=0 ; w7s=0                         
! Irreducible secondaries                                                    
        do i0=0, 4                                                                 
         do i1=0, 4                                                                
         if (i1.ne.i0) then                                                          
          if (3.le.mdeg) then                                                        
          do i2=0, 4                                                               
          if (i2.ne.i1.and.i2.ne.i0) then                                            
           zs(0)=zs(0)+d2(i0,i1)*d(i0,i2)/60                                       
           zs(1)=zs(1)+d(i0,i1)*d(i0,i2)*d(i1,i2)/60                               
           us(1)=us(1)+d3(i0,i1)*d(i0,i2)/60                                       
           us(2)=us(2)+d2(i0,i1)*d2(i0,i2)/60                                      
           us(3)=us(3)+d2(i0,i1)*d(i0,i2)*d(i1,i2)/60                              
        if (5.le.mdeg) then                                                          
           vs(12)=vs(12)+d4(i0,i1)*d(i0,i2)/60                                     
           vs(13)=vs(13)+d3(i0,i1)*d2(i0,i2)/60                                    
           vs(14)=vs(14)+d3(i0,i1)*d(i0,i2)*d(i1,i2)/60                            
           vs(15)=vs(15)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)/60                           
        endif                                                                        
        if (6.le.mdeg) then                                                          
           ws(114)=ws(114)+d5(i0,i1)*d(i0,i2)/60                                   
           ws(115)=ws(115)+d4(i0,i1)*d2(i0,i2)/60                                  
           ws(116)=ws(116)+d3(i0,i1)*d3(i0,i2)/60                                  
           ws(117)=ws(117)+d4(i0,i1)*d(i0,i2)*d(i1,i2)/60                          
           ws(118)=ws(118)+d3(i0,i1)*d2(i0,i2)*d(i1,i2)/60                         
           ws(119)=ws(119)+d2(i0,i1)*d2(i0,i2)*d2(i1,i2)/60                        
        endif                                                                        
           do i3=0, 4                                                              
           if (i3.ne.i0.and.i3.ne.i1.and.i3.ne.i2) then                              
            us(4)=us(4)+d2(i0,i1)*d(i0,i2)*d(i0,i3)/120                            
            us(5)=us(5)+d2(i0,i1)*d(i1,i2)*d(i0,i3)/120                            
        if (5.le.mdeg) then                                                          
            vs(16)=vs(16)+d3(i0,i1)*d(i0,i2)*d(i0,i3)/120                          
            vs(17)=vs(17)+d2(i0,i1)*d2(i0,i2)*d(i0,i3)/120                         
            vs(18)=vs(18)+d3(i0,i1)*d(i1,i2)*d(i0,i3)/120                          
            vs(19)=vs(19)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)/120                 
        endif                                                                        
        if (6.le.mdeg) then                                                          
            ws(120)=ws(120)+d4(i0,i1)*d(i0,i2)*d(i0,i3)/120                        
            ws(121)=ws(121)+d3(i0,i1)*d2(i0,i2)*d(i0,i3)/120                       
            ws(122)=ws(122)+d4(i0,i1)*d(i1,i2)*d(i0,i3)/120                        
            ws(123)=ws(123)+d3(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)/120               
            ws(124)=ws(124)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i0,i3)/120              
            ws(125)=ws(125)+d3(i0,i1)*d2(i1,i2)*d(i0,i3)/120                       
        endif                                                                        
           endif                                                                     
           enddo                                                                     
          endif                                                                      
          enddo                                                                      
          endif                                                                      
         endif                                                                       
         enddo                                                                       
        enddo                                                                        
        do j0=5, 6                                                                 
         do i0=0, 4                                                                
          zs(2)=zs(2)+d3(i0,j0)/10                                                 
          vs(20)=vs(20)+d5(i0,j0)/10                                               
          do i1=0, 4                                                               
          if (i1.ne.i0) then                                                         
           ys(0)=ys(0)+d(i0,i1)*d(i0,j0)/40                                        
           if (3.le.mdeg) then                                                       
           zs(3)=zs(3)+d2(i0,i1)*d(i0,j0)/40                                       
           zs(4)=zs(4)+d(i0,i1)*d2(i0,j0)/40                                       
           zs(5)=zs(5)+d(i0,i1)*d(i0,j0)*d(i1,j0)/40                               
           zs(6)=zs(6)+d2(i0,j0)*d(i1,j0)/40                                       
           us(6)=us(6)+d3(i0,i1)*d(i0,j0)/40                                       
           us(7)=us(7)+d2(i0,i1)*d2(i0,j0)/40                                      
           us(8)=us(8)+d(i0,i1)*d3(i0,j0)/40                                       
           us(9)=us(9)+d3(i0,j0)*d(i1,j0)/40                                       
           us(10)=us(10)+d2(i0,j0)*d2(i1,j0)/40                                    
           us(11)=us(11)+d2(i0,i1)*d(i0,j0)*d(i1,j0)/40                            
           us(12)=us(12)+d(i0,i1)*d2(i0,j0)*d(i1,j0)/40                            
        if (5.le.mdeg) then                                                          
           vs(21)=vs(21)+d4(i0,i1)*d(i0,j0)/40                                     
           vs(22)=vs(22)+d3(i0,i1)*d2(i0,j0)/40                                    
           vs(23)=vs(23)+d2(i0,i1)*d3(i0,j0)/40                                    
           vs(24)=vs(24)+d(i0,i1)*d4(i0,j0)/40                                     
           vs(25)=vs(25)+d3(i0,i1)*d(i0,j0)*d(i1,j0)/40                            
           vs(26)=vs(26)+d2(i0,i1)*d2(i0,j0)*d(i1,j0)/40                           
           vs(27)=vs(27)+d(i0,i1)*d3(i0,j0)*d(i1,j0)/40                            
           vs(28)=vs(28)+d4(i0,j0)*d(i1,j0)/40                                     
           vs(29)=vs(29)+d(i0,i1)*d2(i0,j0)*d2(i1,j0)/40                           
           vs(30)=vs(30)+d3(i0,j0)*d2(i1,j0)/40                                    
        endif                                                                        
        if (6.le.mdeg) then                                                          
           ws(126)=ws(126)+d5(i0,i1)*d(i0,j0)/40                                   
           ws(127)=ws(127)+d4(i0,i1)*d2(i0,j0)/40                                  
           ws(128)=ws(128)+d3(i0,i1)*d3(i0,j0)/40                                  
           ws(129)=ws(129)+d2(i0,i1)*d4(i0,j0)/40                                  
           ws(130)=ws(130)+d(i0,i1)*d5(i0,j0)/40                                   
           ws(131)=ws(131)+d4(i0,i1)*d(i0,j0)*d(i1,j0)/40                          
           ws(132)=ws(132)+d3(i0,i1)*d2(i0,j0)*d(i1,j0)/40                         
           ws(133)=ws(133)+d2(i0,i1)*d3(i0,j0)*d(i1,j0)/40                         
           ws(134)=ws(134)+d(i0,i1)*d4(i0,j0)*d(i1,j0)/40                          
           ws(135)=ws(135)+d5(i0,j0)*d(i1,j0)/40                                   
           ws(136)=ws(136)+d2(i0,i1)*d2(i0,j0)*d2(i1,j0)/40                        
           ws(137)=ws(137)+d4(i0,j0)*d2(i1,j0)/40                                  
           ws(138)=ws(138)+d(i0,i1)*d3(i0,j0)*d2(i1,j0)/40                         
        endif                                                                        
           do i2=0, 4                                                              
           if (i2.ne.i1.and.i2.ne.i0) then                                           
            zs(7)=zs(7)+d(i0,i1)*d(i0,i2)*d(i0,j0)/120                             
            zs(8)=zs(8)+d(i0,i1)*d(i1,i2)*d(i0,j0)/120                             
            zs(9)=zs(9)+d(i0,i2)*d(i0,j0)*d(i1,j0)/120                             
            us(13)=us(13)+d2(i0,i1)*d(i0,i2)*d(i0,j0)/120                          
            us(14)=us(14)+d2(i0,i1)*d(i1,i2)*d(i0,j0)/120                          
            us(15)=us(15)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)/120                  
            us(16)=us(16)+d(i0,i1)*d2(i1,i2)*d(i0,j0)/120                          
            us(17)=us(17)+d(i0,i1)*d(i0,i2)*d2(i0,j0)/120                          
            us(18)=us(18)+d(i0,i1)*d(i1,i2)*d2(i0,j0)/120                          
            us(19)=us(19)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)/120                  
            us(20)=us(20)+d2(i0,i2)*d(i0,j0)*d(i1,j0)/120                          
            us(21)=us(21)+d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,j0)/120                  
            us(22)=us(22)+d(i0,i2)*d2(i0,j0)*d(i1,j0)/120                          
            us(23)=us(23)+d(i1,i2)*d2(i0,j0)*d(i1,j0)/120                          
            us(24)=us(24)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i2,j0)/120                  
        if (5.le.mdeg) then                                                          
            vs(31)=vs(31)+d3(i0,i1)*d(i0,i2)*d(i0,j0)/120                          
            vs(32)=vs(32)+d2(i0,i1)*d2(i0,i2)*d(i0,j0)/120                         
            vs(33)=vs(33)+d3(i0,i1)*d(i1,i2)*d(i0,j0)/120                          
            vs(34)=vs(34)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)/120                 
            vs(35)=vs(35)+d2(i0,i1)*d2(i1,i2)*d(i0,j0)/120                         
            vs(36)=vs(36)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,j0)/120                 
            vs(37)=vs(37)+d(i0,i1)*d3(i1,i2)*d(i0,j0)/120                          
            vs(38)=vs(38)+d2(i0,i1)*d(i0,i2)*d2(i0,j0)/120                         
            vs(39)=vs(39)+d2(i0,i1)*d(i1,i2)*d2(i0,j0)/120                         
            vs(40)=vs(40)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,j0)/120                 
            vs(41)=vs(41)+d(i0,i1)*d2(i1,i2)*d2(i0,j0)/120                         
            vs(42)=vs(42)+d(i0,i1)*d(i0,i2)*d3(i0,j0)/120                          
            vs(43)=vs(43)+d(i0,i1)*d(i1,i2)*d3(i0,j0)/120                          
            vs(44)=vs(44)+d2(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)/120                 
            vs(45)=vs(45)+d(i0,i1)*d2(i0,i2)*d(i0,j0)*d(i1,j0)/120                 
            vs(46)=vs(46)+d3(i0,i2)*d(i0,j0)*d(i1,j0)/120                          
            vs(47)=vs(47)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,
     $j0)/120         
            vs(48)=vs(48)+d2(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,j0)/120                 
            vs(49)=vs(49)+d(i0,i1)*d(i0,i2)*d2(i0,j0)*d(i1,j0)/120                 
            vs(50)=vs(50)+d2(i0,i2)*d2(i0,j0)*d(i1,j0)/120                         
            vs(51)=vs(51)+d(i0,i1)*d(i1,i2)*d2(i0,j0)*d(i1,j0)/120                 
            vs(52)=vs(52)+d(i0,i2)*d(i1,i2)*d2(i0,j0)*d(i1,j0)/120                 
            vs(53)=vs(53)+d2(i1,i2)*d2(i0,j0)*d(i1,j0)/120                         
            vs(54)=vs(54)+d(i0,i2)*d3(i0,j0)*d(i1,j0)/120                          
            vs(55)=vs(55)+d(i1,i2)*d3(i0,j0)*d(i1,j0)/120                          
            vs(56)=vs(56)+d(i0,i2)*d2(i0,j0)*d2(i1,j0)/120                         
            vs(57)=vs(57)+d2(i0,i1)*d(i0,j0)*d(i1,j0)*d(i2,j0)/120                 
            vs(58)=vs(58)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)*d(i2,
     $j0)/120         
            vs(59)=vs(59)+d(i0,i1)*d2(i0,j0)*d(i1,j0)*d(i2,j0)/120                 
        endif                                                                        
        if (6.le.mdeg) then                                                          
            ws(139)=ws(139)+d4(i0,i1)*d(i0,i2)*d(i0,j0)/120                        
            ws(140)=ws(140)+d3(i0,i1)*d2(i0,i2)*d(i0,j0)/120                       
            ws(141)=ws(141)+d4(i0,i1)*d(i1,i2)*d(i0,j0)/120                        
            ws(142)=ws(142)+d3(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)/120               
            ws(143)=ws(143)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i0,j0)/120              
            ws(144)=ws(144)+d3(i0,i1)*d2(i1,i2)*d(i0,j0)/120                       
            ws(145)=ws(145)+d2(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,j0)/120              
            ws(146)=ws(146)+d2(i0,i1)*d3(i1,i2)*d(i0,j0)/120                       
            ws(147)=ws(147)+d(i0,i1)*d(i0,i2)*d3(i1,i2)*d(i0,j0)/120               
            ws(148)=ws(148)+d(i0,i1)*d4(i1,i2)*d(i0,j0)/120                        
            ws(149)=ws(149)+d3(i0,i1)*d(i0,i2)*d2(i0,j0)/120                       
            ws(150)=ws(150)+d2(i0,i1)*d2(i0,i2)*d2(i0,j0)/120                      
            ws(151)=ws(151)+d3(i0,i1)*d(i1,i2)*d2(i0,j0)/120                       
            ws(152)=ws(152)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,j0)/120              
            ws(153)=ws(153)+d2(i0,i1)*d2(i1,i2)*d2(i0,j0)/120                      
            ws(154)=ws(154)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d2(i0,j0)/120              
            ws(155)=ws(155)+d(i0,i1)*d3(i1,i2)*d2(i0,j0)/120                       
            ws(156)=ws(156)+d2(i0,i1)*d(i0,i2)*d3(i0,j0)/120                       
            ws(157)=ws(157)+d2(i0,i1)*d(i1,i2)*d3(i0,j0)/120                       
            ws(156)=ws(156)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d3(i0,j0)/120               
            ws(159)=ws(159)+d(i0,i1)*d2(i1,i2)*d3(i0,j0)/120                       
            ws(160)=ws(160)+d(i0,i1)*d(i0,i2)*d4(i0,j0)/120                        
            ws(161)=ws(161)+d(i0,i1)*d(i1,i2)*d4(i0,j0)/120                        
            ws(162)=ws(162)+d3(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)/120               
            ws(163)=ws(163)+d2(i0,i1)*d2(i0,i2)*d(i0,j0)*d(i1,j0)/120              
            ws(164)=ws(164)+d(i0,i1)*d3(i0,i2)*d(i0,j0)*d(i1,j0)/120               
            ws(165)=ws(165)+d4(i0,i2)*d(i0,j0)*d(i1,j0)/120                        
            ws(166)=ws(166)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,
     $j0)/120      
            ws(167)=ws(167)+d(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,
     $j0)/120      
            ws(168)=ws(168)+d3(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,j0)/120               
            ws(169)=ws(169)+d2(i0,i2)*d2(i1,i2)*d(i0,j0)*d(i1,j0)/120              
            ws(170)=ws(170)+d2(i0,i1)*d(i0,i2)*d2(i0,j0)*d(i1,j0)/120              
            ws(171)=ws(171)+d(i0,i1)*d2(i0,i2)*d2(i0,j0)*d(i1,j0)/120              
            ws(172)=ws(172)+d3(i0,i2)*d2(i0,j0)*d(i1,j0)/120                       
            ws(173)=ws(173)+d2(i0,i1)*d(i1,i2)*d2(i0,j0)*d(i1,j0)/120              
            ws(174)=ws(174)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,j0)*d(i1,
     $j0)/120      
            ws(175)=ws(175)+d2(i0,i2)*d(i1,i2)*d2(i0,j0)*d(i1,j0)/120              
            ws(176)=ws(176)+d(i0,i1)*d2(i1,i2)*d2(i0,j0)*d(i1,j0)/120              
            ws(177)=ws(177)+d(i0,i2)*d2(i1,i2)*d2(i0,j0)*d(i1,j0)/120              
            ws(178)=ws(178)+d3(i1,i2)*d2(i0,j0)*d(i1,j0)/120                       
            ws(179)=ws(179)+d(i0,i1)*d(i0,i2)*d3(i0,j0)*d(i1,j0)/120               
            ws(180)=ws(180)+d2(i0,i2)*d3(i0,j0)*d(i1,j0)/120                       
            ws(181)=ws(181)+d(i0,i1)*d(i1,i2)*d3(i0,j0)*d(i1,j0)/120               
            ws(182)=ws(182)+d(i0,i2)*d(i1,i2)*d3(i0,j0)*d(i1,j0)/120               
            ws(183)=ws(183)+d2(i1,i2)*d3(i0,j0)*d(i1,j0)/120                       
            ws(184)=ws(184)+d(i0,i2)*d4(i0,j0)*d(i1,j0)/120                        
            ws(185)=ws(185)+d(i1,i2)*d4(i0,j0)*d(i1,j0)/120                        
            ws(186)=ws(186)+d(i0,i1)*d(i0,i2)*d2(i0,j0)*d2(i1,j0)/120              
            ws(187)=ws(187)+d2(i0,i2)*d2(i0,j0)*d2(i1,j0)/120                      
            ws(188)=ws(188)+d(i0,i2)*d(i1,i2)*d2(i0,j0)*d2(i1,j0)/120              
            ws(189)=ws(189)+d(i0,i2)*d3(i0,j0)*d2(i1,j0)/120                       
            ws(190)=ws(190)+d3(i0,i1)*d(i0,j0)*d(i1,j0)*d(i2,j0)/120               
            ws(191)=ws(191)+d2(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)*d(i2,
     $j0)/120      
            ws(192)=ws(192)+                                                      
     $       d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,j0)*d(i2,j0)/120              
            ws(193)=ws(193)+d2(i0,i1)*d2(i0,j0)*d(i1,j0)*d(i2,j0)/120              
            ws(194)=ws(194)+d(i0,i1)*d(i0,i2)*d2(i0,j0)*d(i1,j0)*d(i2,
     $j0)/120      
            ws(195)=ws(195)+d(i0,i1)*d(i1,i2)*d2(i0,j0)*d(i1,j0)*d(i2,
     $j0)/120      
            ws(196)=ws(196)+d(i0,i1)*d3(i0,j0)*d(i1,j0)*d(i2,j0)/120               
            ws(197)=ws(197)+d(i1,i2)*d3(i0,j0)*d(i1,j0)*d(i2,j0)/120               
        endif                                                                        
            do i3=0, 4                                                             
            if (i3.ne.i0.and.i3.ne.i1.and.i3.ne.i2) then                             
             us(25)=us(25)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)/120                 
             us(26)=us(26)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)/120                 
             us(27)=us(27)+d(i0,i2)*d(i0,i3)*d(i0,j0)*d(i1,j0)/120                 
             us(28)=us(28)+d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,j0)/120                 
        if (5.le.mdeg) then                                                          
             vs(60)=vs(60)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)/240                
             vs(61)=vs(61)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)/240                
             vs(62)=vs(62)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,
     $j0)/240        
             vs(63)=vs(63)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,j0)/240                
             vs(64)=vs(64)+d(i0,i1)*d(i1,i2)*d2(i0,i3)*d(i0,j0)/240                
             vs(65)=vs(65)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d2(i0,j0)/240                
             vs(66)=vs(66)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d2(i0,j0)/240                
             vs(67)=vs(67)+d2(i0,i2)*d(i0,i3)*d(i0,j0)*d(i1,j0)/240                
             vs(68)=vs(68)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240        
             vs(69)=vs(69)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240        
             vs(70)=vs(70)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240        
             vs(71)=vs(71)+d2(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,j0)/240                
             vs(72)=vs(72)+d2(i0,i2)*d(i2,i3)*d(i0,j0)*d(i1,j0)/240                
             vs(73)=vs(73)+d(i0,i2)*d(i1,i2)*d(i2,i3)*d(i0,j0)*d(i1,
     $j0)/240        
             vs(74)=vs(74)+d(i0,i2)*d(i0,i3)*d(i2,i3)*d(i0,j0)*d(i1,
     $j0)/240        
             vs(75)=vs(75)+d(i0,i2)*d(i0,i3)*d2(i0,j0)*d(i1,j0)/240                
             vs(76)=vs(76)+d(i1,i2)*d(i0,i3)*d2(i0,j0)*d(i1,j0)/240                
             vs(77)=vs(77)+d(i1,i2)*d(i1,i3)*d2(i0,j0)*d(i1,j0)/240                
             vs(78)=vs(78)+d(i0,i2)*d(i2,i3)*d2(i0,j0)*d(i1,j0)/240                
        endif                                                                        
        if (6.le.mdeg) then                                                          
             ws(198)=ws(198)+d3(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)/240              
             ws(199)=ws(199)+d2(i0,i1)*d2(i0,i2)*d(i0,i3)*d(i0,j0)/240             
             ws(200)=ws(200)+d3(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)/240              
             ws(201)=ws(201)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,
     $j0)/240     
             ws(202)=ws(202)+d2(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,j0)/240             
             ws(203)=ws(203)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i0,
     $j0)/240     
             ws(204)=ws(204)+d2(i0,i1)*d(i1,i2)*d2(i0,i3)*d(i0,j0)/240             
             ws(205)=ws(205)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d2(i0,j0)/240             
             ws(206)=ws(206)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d2(i0,j0)/240             
             ws(207)=ws(207)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d2(i0,
     $j0)/240     
             ws(208)=ws(208)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d2(i0,j0)/240             
             ws(209)=ws(209)+d(i0,i1)*d(i1,i2)*d2(i0,i3)*d2(i0,j0)/240             
             ws(210)=ws(210)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d3(i0,j0)/240              
             ws(211)=ws(211)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d3(i0,j0)/240              
             ws(212)=ws(212)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(213)=ws(213)+d(i0,i1)*d2(i0,i2)*d(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(214)=ws(214)+d3(i0,i2)*d(i0,i3)*d(i0,j0)*d(i1,j0)/240              
             ws(215)=ws(215)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(216)=ws(216)+                                                     
     $       d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,j0)/240             
             ws(217)=ws(217)+d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(218)=ws(218)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(219)=ws(219)+d(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(220)=ws(220)+d3(i1,i2)*d(i0,i3)*d(i0,j0)*d(i1,j0)/240              
             ws(221)=ws(221)+d2(i0,i2)*d2(i0,i3)*d(i0,j0)*d(i1,j0)/240             
             ws(222)=ws(222)+d(i0,i2)*d(i1,i2)*d2(i0,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(223)=ws(223)+d2(i1,i2)*d2(i0,i3)*d(i0,j0)*d(i1,j0)/240             
             ws(224)=ws(224)+                                                     
     $        d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,j0)*d(i1,j0)/240             
             ws(225)=ws(225)+d(i0,i1)*d2(i0,i2)*d(i2,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(226)=ws(226)+d3(i0,i2)*d(i2,i3)*d(i0,j0)*d(i1,j0)/240              
             ws(227)=ws(227)+                                                     
     $        d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i2,i3)*d(i0,j0)*d(i1,j0)/240             
             ws(228)=ws(228)+d2(i0,i2)*d(i1,i2)*d(i2,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(229)=ws(229)+                                                     
     $        d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i2,i3)*d(i0,j0)*d(i1,j0)/240             
             ws(230)=ws(230)+d2(i0,i2)*d(i0,i3)*d(i2,i3)*d(i0,j0)*d(i1,
     $j0)/240     
             ws(231)=ws(231)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d2(i0,j0)*d(i1,
     $j0)/240     
             ws(232)=ws(232)+d2(i0,i2)*d(i0,i3)*d2(i0,j0)*d(i1,j0)/240             
             ws(233)=ws(233)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d2(i0,j0)*d(i1,
     $j0)/240     
             ws(234)=ws(234)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d2(i0,j0)*d(i1,
     $j0)/240     
             ws(235)=ws(235)+d2(i1,i2)*d(i0,i3)*d2(i0,j0)*d(i1,j0)/240             
             ws(236)=ws(236)+d(i1,i2)*d2(i0,i3)*d2(i0,j0)*d(i1,j0)/240             
             ws(237)=ws(237)+d(i0,i1)*d(i1,i2)*d(i1,i3)*d2(i0,j0)*d(i1,
     $j0)/240     
             ws(238)=ws(238)+d(i0,i2)*d(i1,i2)*d(i1,i3)*d2(i0,j0)*d(i1,
     $j0)/240     
             ws(239)=ws(239)+d2(i1,i2)*d(i1,i3)*d2(i0,j0)*d(i1,j0)/240             
             ws(240)=ws(240)+d(i0,i1)*d(i0,i2)*d(i2,i3)*d2(i0,j0)*d(i1,
     $j0)/240     
             ws(241)=ws(241)+d2(i0,i2)*d(i2,i3)*d2(i0,j0)*d(i1,j0)/240             
             ws(242)=ws(242)+d(i0,i2)*d(i1,i2)*d(i2,i3)*d2(i0,j0)*d(i1,
     $j0)/240     
             ws(243)=ws(243)+d2(i1,i2)*d(i2,i3)*d2(i0,j0)*d(i1,j0)/240             
             ws(244)=ws(244)+d(i0,i2)*d(i0,i3)*d(i2,i3)*d2(i0,j0)*d(i1,
     $j0)/240     
             ws(245)=ws(245)+d(i1,i2)*d(i0,i3)*d(i2,i3)*d2(i0,j0)*d(i1,
     $j0)/240     
             ws(246)=ws(246)+d(i0,i2)*d2(i2,i3)*d2(i0,j0)*d(i1,j0)/240             
             ws(248)=ws(248)+d(i0,i2)*d(i0,i3)*d3(i0,j0)*d(i1,j0)/240              
             ws(249)=ws(249)+d(i1,i2)*d(i0,i3)*d3(i0,j0)*d(i1,j0)/240              
             ws(250)=ws(250)+d(i1,i2)*d(i1,i3)*d3(i0,j0)*d(i1,j0)/240              
             ws(251)=ws(251)+d(i0,i2)*d(i2,i3)*d3(i0,j0)*d(i1,j0)/240              
             ws(252)=ws(252)+d(i0,i2)*d(i0,i3)*d2(i0,j0)*d2(i1,j0)/240             
             ws(253)=ws(253)+d2(i0,i1)*d(i0,i3)*d(i0,j0)*d(i1,j0)*d(i2,
     $j0)/240     
             ws(254)=ws(254)+d(i1,i2)*d(i0,i3)*d2(i0,j0)*d(i1,j0)*d(i2,
     $j0)/240     
             do i4=0, 4                                                            
             if (i4.ne.i0.and.i4.ne.i1.and.i4.ne.i2.and.i4.ne.i3) then               
              ws(247)=ws(247)+d(i0,i2)*d(i0,i3)*d(i0,i4)*d2(i0,j0)*d(i1,
     $j0)/240    
             endif                                                                   
             enddo                                                                   
        endif                                                                        
            endif                                                                    
            enddo                                                                    
           endif                                                                     
           enddo                                                                     
           endif                                                                     
          endif                                                                      
          enddo                                                                      
         enddo                                                                       
        enddo                                                                        
        if (3.le.mdeg) then                                                          
        do j0=5, 6                                                                 
         do j1=5, 6                                                                
         if (j1.ne.j0) then                                                          
          do i0=0, 4                                                               
        if (6.le.mdeg) then                                                          
           ws(255)=ws(255)+d5(i0,j0)*d(i0,j1)/10                                   
           ws(256)=ws(256)+d4(i0,j0)*d2(i0,j1)/10                                  
        endif                                                                        
           do i1=0, 4                                                              
           if (i1.ne.i0) then                                                        
            zs(10)=zs(10)+d(i0,i1)*d(i0,j0)*d(i0,j1)/40                            
            zs(11)=zs(11)+d(i0,i1)*d(i1,j0)*d(i0,j1)/40                            
            us(29)=us(29)+d2(i0,i1)*d(i0,j0)*d(i0,j1)/40                           
            us(30)=us(30)+d(i0,i1)*d2(i0,j0)*d(i0,j1)/40                           
            us(31)=us(31)+d3(i0,j0)*d(i0,j1)/40                                    
            us(32)=us(32)+d2(i0,i1)*d(i1,j0)*d(i0,j1)/40                           
            us(33)=us(33)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,j1)/40                   
            us(34)=us(34)+d2(i0,j0)*d(i1,j0)*d(i0,j1)/40                           
            us(35)=us(35)+d(i0,i1)*d2(i1,j0)*d(i0,j1)/40                           
        if (5.le.mdeg) then                                                          
            vs(79)=vs(79)+d3(i0,i1)*d(i0,j0)*d(i0,j1)/40                           
            vs(80)=vs(80)+d2(i0,i1)*d2(i0,j0)*d(i0,j1)/40                          
            vs(81)=vs(81)+d(i0,i1)*d3(i0,j0)*d(i0,j1)/40                           
            vs(82)=vs(82)+d4(i0,j0)*d(i0,j1)/40                                    
            vs(83)=vs(83)+d3(i0,i1)*d(i1,j0)*d(i0,j1)/40                           
            vs(84)=vs(84)+d2(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,j1)/40                  
            vs(85)=vs(85)+d(i0,i1)*d2(i0,j0)*d(i1,j0)*d(i0,j1)/40                  
            vs(86)=vs(86)+d3(i0,j0)*d(i1,j0)*d(i0,j1)/40                           
            vs(87)=vs(87)+d2(i0,i1)*d2(i1,j0)*d(i0,j1)/40                          
            vs(88)=vs(88)+d(i0,i1)*d(i0,j0)*d2(i1,j0)*d(i0,j1)/40                  
            vs(89)=vs(89)+d2(i0,j0)*d2(i1,j0)*d(i0,j1)/40                          
            vs(90)=vs(90)+d(i0,i1)*d3(i1,j0)*d(i0,j1)/40                           
            vs(91)=vs(91)+d(i0,i1)*d2(i0,j0)*d2(i0,j1)/40                          
            vs(92)=vs(92)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d2(i0,j1)/40                  
            vs(93)=vs(93)+d(i0,i1)*d2(i1,j0)*d2(i0,j1)/40                          
        endif                                                                        
        if (6.le.mdeg) then                                                          
            ws(257)=ws(257)+d4(i0,i1)*d(i0,j0)*d(i0,j1)/40                         
            ws(258)=ws(258)+d3(i0,i1)*d2(i0,j0)*d(i0,j1)/40                        
            ws(259)=ws(259)+d2(i0,i1)*d3(i0,j0)*d(i0,j1)/40                        
            ws(260)=ws(260)+d(i0,i1)*d4(i0,j0)*d(i0,j1)/40                         
            ws(261)=ws(261)+d4(i0,i1)*d(i1,j0)*d(i0,j1)/40                         
            ws(262)=ws(262)+d3(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,j1)/40                
            ws(263)=ws(263)+d2(i0,i1)*d2(i0,j0)*d(i1,j0)*d(i0,j1)/40               
            ws(264)=ws(264)+d(i0,i1)*d3(i0,j0)*d(i1,j0)*d(i0,j1)/40                
            ws(265)=ws(265)+d4(i0,j0)*d(i1,j0)*d(i0,j1)/40                         
            ws(266)=ws(266)+d3(i0,i1)*d2(i1,j0)*d(i0,j1)/40                        
            ws(267)=ws(267)+d2(i0,i1)*d(i0,j0)*d2(i1,j0)*d(i0,j1)/40               
            ws(268)=ws(268)+d(i0,i1)*d2(i0,j0)*d2(i1,j0)*d(i0,j1)/40               
            ws(269)=ws(269)+d3(i0,j0)*d2(i1,j0)*d(i0,j1)/40                        
            ws(270)=ws(270)+d2(i0,i1)*d3(i1,j0)*d(i0,j1)/40                        
            ws(271)=ws(271)+d(i0,i1)*d(i0,j0)*d3(i1,j0)*d(i0,j1)/40                
            ws(272)=ws(272)+d(i0,i1)*d4(i1,j0)*d(i0,j1)/40                         
            ws(273)=ws(273)+d2(i0,i1)*d2(i0,j0)*d2(i0,j1)/40                       
            ws(274)=ws(274)+d2(i0,i1)*d(i0,j0)*d(i1,j0)*d2(i0,j1)/40               
            ws(275)=ws(275)+d(i0,i1)*d2(i0,j0)*d(i1,j0)*d2(i0,j1)/40               
            ws(276)=ws(276)+d3(i0,j0)*d(i1,j0)*d2(i0,j1)/40                        
            ws(277)=ws(277)+d2(i0,i1)*d2(i1,j0)*d2(i0,j1)/40                       
            ws(278)=ws(278)+d(i0,i1)*d(i0,j0)*d2(i1,j0)*d2(i0,j1)/40               
            ws(279)=ws(279)+d(i0,i1)*d3(i1,j0)*d2(i0,j1)/40                        
        endif                                                                        
            do i2=0, 4                                                             
            if (i2.ne.i1.and.i2.ne.i0) then                                          
             us(36)=us(36)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,j1)/120                 
             us(37)=us(37)+d(i0,i1)*d(i1,i2)*d(i0,j0)*d(i0,j1)/120                 
             us(38)=us(38)+d(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,j1)/120                 
        if (5.le.mdeg) then                                                          
             vs(94)=vs(94)+d2(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,j1)/120                
             vs(95)=vs(95)+d2(i0,i1)*d(i1,i2)*d(i0,j0)*d(i0,j1)/120                
             vs(96)=vs(96)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i0,j1
     $)/120        
             vs(97)=vs(97)+d(i0,i1)*d2(i1,i2)*d(i0,j0)*d(i0,j1)/120                
             vs(98)=vs(98)+d(i0,i1)*d(i0,i2)*d2(i0,j0)*d(i0,j1)/120                
             vs(99)=vs(99)+d(i0,i1)*d(i1,i2)*d2(i0,j0)*d(i0,j1)/120                
             vs(100)=vs(100)+d2(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,j1)/120              
             vs(101)=vs(101)+d(i0,i1)*d2(i0,i2)*d(i1,j0)*d(i0,j1)/120              
             vs(102)=vs(102)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,j0)*d(i0,j1
     $)/120      
             vs(103)=vs(103)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)*d(i0,j1
     $)/120      
             vs(104)=vs(104)+d(i0,i1)*d(i1,i2)*d(i0,j0)*d(i1,j0)*d(i0,j1
     $)/120      
             vs(105)=vs(105)+d(i0,i2)*d2(i0,j0)*d(i1,j0)*d(i0,j1)/120              
             vs(106)=vs(106)+d(i1,i2)*d2(i0,j0)*d(i1,j0)*d(i0,j1)/120              
             vs(107)=vs(107)+d(i0,i1)*d(i0,i2)*d2(i1,j0)*d(i0,j1)/120              
             vs(108)=vs(108)+d(i0,i1)*d(i0,i2)*d(i1,j0)*d(i2,j0)*d(i0,j1
     $)/120      
             vs(109)=vs(109)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i2,j0)*d(i0,j1
     $)/120      
        endif                                                                        
        if (6.le.mdeg) then  

             di0j1=d(i0,j1)/120

             ws(280)=ws(280)+d3(i0,i1)*d(i0,i2)*d(i0,j0)*di0j1              
             ws(281)=ws(281)+d2(i0,i1)*d2(i0,i2)*d(i0,j0)*di0j1             
             ws(282)=ws(282)+d3(i0,i1)*d(i1,i2)*d(i0,j0)*di0j1              
             ws(283)=ws(283)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)*di0j1     
             ws(284)=ws(284)+d2(i0,i1)*d2(i1,i2)*d(i0,j0)*di0j1             
             ws(285)=ws(285)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,j0)*di0j1     
             ws(286)=ws(286)+d(i0,i1)*d3(i1,i2)*d(i0,j0)*di0j1              
             ws(287)=ws(287)+d2(i0,i1)*d(i0,i2)*d2(i0,j0)*di0j1             
             ws(288)=ws(288)+d2(i0,i1)*d(i1,i2)*d2(i0,j0)*di0j1             
             ws(289)=ws(289)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,j0)*di0j1     
             ws(290)=ws(290)+d(i0,i1)*d2(i1,i2)*d2(i0,j0)*di0j1             
             ws(291)=ws(291)+d(i0,i1)*d(i0,i2)*d3(i0,j0)*di0j1              
             ws(292)=ws(292)+d(i0,i1)*d(i1,i2)*d3(i0,j0)*di0j1              
             ws(293)=ws(293)+d3(i0,i1)*d(i0,i2)*d(i1,j0)*di0j1              
             ws(294)=ws(294)+d(i0,i1)*d3(i0,i2)*d(i1,j0)*di0j1              
             ws(295)=ws(295)+d2(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)*di0j1     
             ws(296)=ws(296)+d(i0,i1)*d2(i0,i2)*d(i0,j0)*d(i1,j0)*di0j1     
             ws(297)=ws(297)+d2(i0,i1)*d(i1,i2)*d(i0,j0)*d(i1,j0)*di0j1     
             ws(298)=ws(298)+                                                     
     $          d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,j0)*di0j1             
             ws(299)=ws(299)+d2(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,j0)*di0j1     
             ws(300)=ws(300)+d(i0,i1)*d2(i1,i2)*d(i0,j0)*d(i1,j0)*di0j1     
             ws(301)=ws(301)+d(i0,i1)*d(i0,i2)*d2(i0,j0)*d(i1,j0)*di0j1     
             ws(302)=ws(302)+d2(i0,i2)*d2(i0,j0)*d(i1,j0)*di0j1             
             ws(303)=ws(303)+d(i0,i1)*d(i1,i2)*d2(i0,j0)*d(i1,j0)*di0j1     
             ws(304)=ws(304)+d(i0,i2)*d(i1,i2)*d2(i0,j0)*d(i1,j0)*di0j1     
             ws(305)=ws(305)+d2(i1,i2)*d2(i0,j0)*d(i1,j0)*di0j1             
             ws(306)=ws(306)+d(i0,i2)*d3(i0,j0)*d(i1,j0)*di0j1              
             ws(307)=ws(307)+d(i1,i2)*d3(i0,j0)*d(i1,j0)*di0j1              
             ws(308)=ws(308)+d2(i0,i1)*d(i0,i2)*d2(i1,j0)*di0j1             
             ws(309)=ws(309)+d(i0,i1)*d2(i0,i2)*d2(i1,j0)*di0j1             
             ws(310)=ws(310)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d2(i1,j0)*di0j1     
             ws(311)=ws(311)+d(i0,i1)*d(i1,i2)*d(i0,j0)*d2(i1,j0)*di0j1     
             ws(312)=ws(312)+d(i0,i2)*d2(i0,j0)*d2(i1,j0)*di0j1             
             ws(313)=ws(313)+d(i0,i1)*d(i0,i2)*d3(i1,j0)*di0j1              
             ws(314)=ws(314)+d(i0,i1)*d(i1,i2)*d3(i1,j0)*di0j1              
             ws(315)=ws(315)+d2(i0,i1)*d(i0,i2)*d(i1,j0)*d(i2,j0)*di0j1     
             ws(316)=ws(316)+d2(i0,i1)*d(i1,i2)*d(i1,j0)*d(i2,j0)*di0j1     
             ws(317)=ws(317)+d2(i0,i1)*d(i0,j0)*d(i1,j0)*d(i2,j0)*di0j1     
             ws(318)=ws(318)+                                                     
     $         d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)*d(i2,j0)*di0j1             
             ws(319)=ws(319)+                                                     
     $         d(i0,i1)*d(i1,i2)*d(i0,j0)*d(i1,j0)*d(i2,j0)*di0j1             
             ws(320)=ws(320)+d(i0,i1)*d2(i0,j0)*d(i1,j0)*d(i2,j0)*di0j1     
             ws(321)=ws(321)+d(i1,i2)*d2(i0,j0)*d(i1,j0)*d(i2,j0)*di0j1     
             ws(322)=ws(322)+d(i0,i1)*d(i0,i2)*d2(i1,j0)*d(i2,j0)*di0j1     
             ws(323)=ws(323)+d(i0,i1)*d(i0,j0)*d2(i1,j0)*d(i2,j0)*di0j1     
             ws(324)=ws(324)+d(i0,i1)*d(i1,i2)*d2(i0,j0)*d2(i0,j1)/120             
             ws(325)=ws(325)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)*d2(i0,
     $j1)/120     
        endif                                                                        
             do i3=0, 4                                                            
             if (i3.ne.i0.and.i3.ne.i1.and.i3.ne.i2) then                            
        if (5.le.mdeg) then                                                          
              vs(110)=vs(110)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)*d(i0,
     $j1)/240     
              vs(111)=vs(111)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i0,
     $j1)/240     
              vs(112)=vs(112)+d(i0,i1)*d(i0,i3)*d(i1,j0)*d(i2,j0)*d(i0,
     $j1)/240     
        endif                                                                        
        if (6.le.mdeg) then                                                          
              ws(326)=ws(326)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)*d(i0,
     $j1)/240    
              ws(327)=ws(327)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i0,
     $j1)/240    
              ws(328)=ws(328)+                                                    
     $        d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i0,j1)/240            
              ws(329)=ws(329)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,j0)*d(i0,
     $j1)/240    
              ws(330)=ws(330)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d2(i0,j0)*d(i0,
     $j1)/240    
              ws(331)=ws(331)+                                                    
     $        d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)*d(i1,j0)*d(i0,j1)/240            
              ws(332)=ws(332)+d(i0,i2)*d(i0,i3)*d2(i0,j0)*d(i1,j0)*d(i0,
     $j1)/240    
              ws(333)=ws(333)+d(i1,i2)*d(i0,i3)*d2(i0,j0)*d(i1,j0)*d(i0,
     $j1)/240    
              ws(334)=ws(334)+d(i0,i2)*d(i2,i3)*d2(i0,j0)*d(i1,j0)*d(i0,
     $j1)/240    
              ws(335)=ws(335)+d2(i0,i1)*d(i0,i3)*d(i1,j0)*d(i2,j0)*d(i0,
     $j1)/240    
              ws(336)=ws(336)+                                                    
     $        d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i1,j0)*d(i2,j0)*d(i0,j1)/240            
              ws(337)=ws(337)+                                                    
     $        d(i0,i1)*d(i0,i3)*d(i0,j0)*d(i1,j0)*d(i2,j0)*d(i0,j1)/240            
        endif                                                                        
             endif                                                                   
             enddo                                                                   
            endif                                                                    
            enddo                                                                    
           endif                                                                     
           enddo                                                                     
          enddo                                                                      
         endif                                                                       
         enddo                                                                       
        enddo                                                                        
        endif                                                                        
! Reducible secondaries                                                      
        us(0)=her2(ys(0))                                                          
        vs(0:11)=ys(0)*zs(0:11)                                                    
        ws(0)=her3(ys(0))                                                          
        ws(1)=her2(zs(0))                                                          
        ws(2)=zs(0)*zs(1)                                                          
        ws(3)=her2(zs(1))                                                          
        ws(4:5)=zs(0:1)*zs(2)                                                      
        ws(6)=her2(zs(2))                                                          
        ws(7:9)=zs(0:2)*zs(3)                                                      
        ws(10)=her2(zs(3))                                                         
        ws(11:14)=zs(0:3)*zs(4)                                                    
        ws(15)=her2(zs(4))                                                         
        ws(16:20)=zs(0:4)*zs(5)                                                    
        ws(21)=her2(zs(5))                                                         
        ws(22:27)=zs(0:5)*zs(6)                                                    
! her2(zs(6)) is not independent                                             
        ws(28:34)=zs(0:6)*zs(7)                                                    
        ws(35)=her2(zs(7))                                                         
        ws(36:43)=zs(0:7)*zs(8)                                                    
        ws(44)=her2(zs(8))                                                         
        ws(45:50)=zs(0:5)*zs(9)                                                    
! zs(6)*zs(9) is not independent                                             
        ws(51:52)=zs(7:8)*zs(9)                                                    
! her2(zs(9)) is not independent                                             
        ws(53:62)=zs(0:9)*zs(10)                                                   
        ws(63)=her2(zs(10))                                                        
        ws(64:74)=zs(0:10)*zs(11)                                                  
        ws(75)=her2(zs(11))                                                        
        ws(76:113)=ys(0)*us(1:38)                                                  
! Compute vec(0:*).  The code was created using these parameters:            
! MolienSeries(0:*): 1 3 12 43 159 554 1879 6066 18755                       
! #Primaries(1:*):   3 5 3 4 3 2 0 0                                         
! #Secondaries(1:*): 0 1 12 39 113 338 932 2402                              
! constant term                                                              
        vec(0)=1                                                                   
! first degree terms                                                         
        if (1.le.mdeg) then                                                          
         vec(1)=x(0)                                                               
         vec(2)=x(1)                                                               
         vec(3)=x(2)                                                               
        endif                                                                        
! second degree terms                                                        
        if (2.le.mdeg) then                                                          
         vec(4)=x2(0)                                                              
         vec(5:6)=x(0)*vec(2:3)                                                    
         vec(7)=x2(1)                                                              
         vec(8:8)=x(1)*vec(3:3)                                                    
         vec(9)=x2(2)                                                              
         vec(10)=y(0)                                                              
         vec(11)=y(1)                                                              
         vec(12)=y(2)                                                              
         vec(13)=y(3)                                                              
         vec(14)=y(4)                                                              
         vec(15:15)=ys(0:0)                                                        
        endif                                                                        
! third degree terms                                                         
        if (3.le.mdeg) then                                                          
         vec(16)=x3(0)                                                             
         vec(17:18)=x2(0)*vec(2:3)                                                 
         vec(19:27)=x(0)*vec(7:15)                                                 
         vec(28)=x3(1)                                                             
         vec(29:29)=x2(1)*vec(3:3)                                                 
         vec(30:36)=x(1)*vec(9:15)                                                 
         vec(37)=x3(2)                                                             
         vec(38:43)=x(2)*vec(10:15)                                                
         vec(44)=z(0)                                                              
         vec(45)=z(1)                                                              
         vec(46)=z(2)                                                              
         vec(47:58)=zs(0:11)                                                       
        endif                                                                        
! fourth degree terms                                                        
        if (4.le.mdeg) then                                                          
         vec(59)=x4(0)                                                             
         vec(60:61)=x3(0)*vec(2:3)                                                 
         vec(62:70)=x2(0)*vec(7:15)                                                
         vec(71:101)=x(0)*vec(28:58)                                               
         vec(102)=x4(1)                                                            
         vec(103:103)=x3(1)*vec(3:3)                                               
         vec(104:110)=x2(1)*vec(9:15)                                              
         vec(111:132)=x(1)*vec(37:58)                                              
         vec(133)=x4(2)                                                            
         vec(134:139)=x2(2)*vec(10:15)                                             
         vec(140:154)=x(2)*vec(44:58)                                              
         vec(155)=y2(0)                                                            
         vec(156:160)=y(0)*vec(11:15)                                              
         vec(161)=y2(1)                                                            
         vec(162:165)=y(1)*vec(12:15)                                              
         vec(166)=y2(2)                                                            
         vec(167:169)=y(2)*vec(13:15)                                              
         vec(170)=y2(3)                                                            
         vec(171:172)=y(3)*vec(14:15)                                              
         vec(173)=y2(4)                                                            
         vec(174:174)=y(4)*vec(15:15)                                              
         vec(175)=u(0)                                                             
         vec(176)=u(1)                                                             
         vec(177)=u(2)                                                             
         vec(178)=u(3)                                                             
         vec(179:217)=us(0:38)                                                     
        endif                                                                        
! fifth degree terms                                                         
        if (5.le.mdeg) then                                                          
         vec(218)=x5(0)                                                            
         vec(219:220)=x4(0)*vec(2:3)                                               
         vec(221:229)=x3(0)*vec(7:15)                                              
         vec(230:260)=x2(0)*vec(28:58)                                             
         vec(261:376)=x(0)*vec(102:217)                                            
         vec(377)=x5(1)                                                            
         vec(378:378)=x4(1)*vec(3:3)                                               
         vec(379:385)=x3(1)*vec(9:15)                                              
         vec(386:407)=x2(1)*vec(37:58)                                             
         vec(408:492)=x(1)*vec(133:217)                                            
         vec(493)=x5(2)                                                            
         vec(494:499)=x3(2)*vec(10:15)                                             
         vec(500:514)=x2(2)*vec(44:58)                                             
         vec(515:577)=x(2)*vec(155:217)                                            
         vec(578:592)=y(0)*vec(44:58)                                              
         vec(593:607)=y(1)*vec(44:58)                                              
         vec(608:622)=y(2)*vec(44:58)                                              
         vec(623:637)=y(3)*vec(44:58)                                              
         vec(638:652)=y(4)*vec(44:58)                                              
         vec(653:653)=z(0)*ys(0:0)                                                 
         vec(654:654)=z(1)*ys(0:0)                                                 
         vec(655:655)=z(2)*ys(0:0)                                                 
         vec(656)=v(0)                                                             
         vec(657)=v(1)                                                             
         vec(658)=v(2)                                                             
         vec(659:771)=vs(0:112)                                                    
        endif                                                                        
! sixth degree terms                                                         
        if (6.le.mdeg) then                                                          
         vec(772)=x6(0)                                                            
         vec(773:774)=x5(0)*vec(2:3)                                               
         vec(775:783)=x4(0)*vec(7:15)                                              
         vec(784:814)=x3(0)*vec(28:58)                                             
         vec(815:930)=x2(0)*vec(102:217)                                           
         vec(931:1325)=x(0)*vec(377:771)                                           
         vec(1326)=x6(1)                                                           
         vec(1327:1327)=x5(1)*vec(3:3)                                             
         vec(1328:1334)=x4(1)*vec(9:15)                                            
         vec(1335:1356)=x3(1)*vec(37:58)                                           
         vec(1357:1441)=x2(1)*vec(133:217)                                         
         vec(1442:1720)=x(1)*vec(493:771)                                          
         vec(1721)=x6(2)                                                           
         vec(1722:1727)=x4(2)*vec(10:15)                                           
         vec(1728:1742)=x3(2)*vec(44:58)                                           
         vec(1743:1805)=x2(2)*vec(155:217)                                         
         vec(1806:1999)=x(2)*vec(578:771)                                          
         vec(2000)=y3(0)                                                           
         vec(2001:2005)=y2(0)*vec(11:15)                                           
         vec(2006:2062)=y(0)*vec(161:217)                                          
         vec(2063)=y3(1)                                                           
         vec(2064:2067)=y2(1)*vec(12:15)                                           
         vec(2068:2119)=y(1)*vec(166:217)                                          
         vec(2120)=y3(2)                                                           
         vec(2121:2123)=y2(2)*vec(13:15)                                           
         vec(2124:2171)=y(2)*vec(170:217)                                          
         vec(2172)=y3(3)                                                           
         vec(2173:2174)=y2(3)*vec(14:15)                                           
         vec(2175:2219)=y(3)*vec(173:217)                                          
         vec(2220)=y3(4)                                                           
         vec(2221:2221)=y2(4)*vec(15:15)                                           
         vec(2222:2264)=y(4)*vec(175:217)                                          
         vec(2265)=z2(0)                                                           
         vec(2266:2279)=z(0)*vec(45:58)                                            
         vec(2280)=z2(1)                                                           
         vec(2281:2293)=z(1)*vec(46:58)                                            
         vec(2294)=z2(2)                                                           
         vec(2295:2306)=z(2)*vec(47:58)                                            
         vec(2307:2307)=u(0)*ys(0:0)                                               
         vec(2308:2308)=u(1)*ys(0:0)                                               
         vec(2309:2309)=u(2)*ys(0:0)                                               
         vec(2310:2310)=u(3)*ys(0:0)                                               
         vec(2311)=w(0)                                                            
         vec(2312)=w(1)                                                            
         vec(2313:2650)=ws(0:337)                                                  
        endif                                                                        
! seventh degree terms                                                       
        if (7.le.mdeg) then                                                          
         vec(2651)=x7(0)                                                           
         vec(2652:2653)=x6(0)*vec(2:3)                                             
         vec(2654:2662)=x5(0)*vec(7:15)                                            
         vec(2663:2693)=x4(0)*vec(28:58)                                           
         vec(2694:2809)=x3(0)*vec(102:217)                                         
         vec(2810:3204)=x2(0)*vec(377:771)                                         
         vec(3205:4529)=x(0)*vec(1326:2650)                                        
         vec(4530)=x7(1)                                                           
         vec(4531:4531)=x6(1)*vec(3:3)                                             
         vec(4532:4538)=x5(1)*vec(9:15)                                            
         vec(4539:4560)=x4(1)*vec(37:58)                                           
         vec(4561:4645)=x3(1)*vec(133:217)                                         
         vec(4646:4924)=x2(1)*vec(493:771)                                         
         vec(4925:5854)=x(1)*vec(1721:2650)                                        
         vec(5855)=x7(2)                                                           
         vec(5856:5855)=x6(2)*vec(4:3)                                             
         vec(5856:5861)=x5(2)*vec(10:15)                                           
         vec(5862:5876)=x4(2)*vec(44:58)                                           
         vec(5877:5939)=x3(2)*vec(155:217)                                         
         vec(5940:6133)=x2(2)*vec(578:771)                                         
         vec(6134:6784)=x(2)*vec(2000:2650)                                        
         vec(6785:6799)=y2(0)*vec(44:58)                                           
         vec(6800:6978)=y(0)*vec(593:771)                                          
         vec(6979:6993)=y2(1)*vec(44:58)                                           
         vec(6994:7157)=y(1)*vec(608:771)                                          
         vec(7158:7172)=y2(2)*vec(44:58)                                           
         vec(7173:7321)=y(2)*vec(623:771)                                          
         vec(7322:7336)=y2(3)*vec(44:58)                                           
         vec(7337:7470)=y(3)*vec(638:771)                                          
         vec(7471:7485)=y2(4)*vec(44:58)                                           
         vec(7486:7604)=y(4)*vec(653:771)                                          
         vec(7605:7647)=z(0)*vec(175:217)                                          
         vec(7648:7690)=z(1)*vec(175:217)                                          
         vec(7691:7733)=z(2)*vec(175:217)                                          
         vec(7734:7745)=u(0)*zs(0:11)                                              
         vec(7746:7757)=u(1)*zs(0:11)                                              
         vec(7758:7769)=u(2)*zs(0:11)                                              
         vec(7770:7781)=u(3)*zs(0:11)                                              
         vec(7782:7782)=v(0)*ys(0:0)                                               
         vec(7783:7783)=v(1)*ys(0:0)                                               
         vec(7784:7784)=v(2)*ys(0:0)                                               
!! vec(7785:8716)=w7s(0:931)                                               
!! w7s(*) to follow                                                          
        endif                                                                        
        return                                                                       
        end                                                                          
        subroutine getrvec (m, r, vec)                                               
        implicit none                                                                
        integer, parameter :: wp=selected_real_kind(12,300)                          
! version for X5Y2                                                           
        integer nk, m                                                                
        parameter (nk=7)                                                             
        double precision  r(0:nk-1,0:nk-1), vec(0:m-1)                                  
        integer i, j                                                                 
        double precision  x(0:2), r1(0:nk-1,0:nk-1), t0
!-----------------------------------------------------------------------     
! Test for compatibility                                                     
        if (.not.(m.eq.1.or.m.eq.4)) then                                            
         stop 'getrvec - wrong dimension'                                            
        endif                                                                        
! Computation                                                                
        x=0                                                                        
        do i=0, nk-1                                                               
         do j=0, nk-1                                                              
          if (i.eq.j) then                                                           
           r1(i,j)=0                                                               
          else                                                                       
           r1(i,j)=exp(-r(i,j))/r(i,j)                                             
          endif                                                                      
         enddo                                                                       
        enddo                                                                        
! XY distance                                                                
        x(0)=sum(r1(0:4,5:6))/10                                                   
! XX distance                                                                
        t0=0                                                                       
        do i=0, 4                                                                  
         do j=i+1, 4                                                               
          t0=t0+r1(i,j)/10                                                         
         enddo                                                                       
        enddo                                                                        
        x(1)=t0                                                                    
! YY distance                                                                
        x(2)=r1(5,6)                                                               
! set vec                                                                    
        vec(0)=1                                                                   
        if (4.le.m) then                                                             
         vec(1:3)=x                                                                
        endif                                                                        
        return                                                                       
        end                 
                                                         
        subroutine prezundeldip()

!        implicit none
        implicit double precision (a-h,o-z)
        implicit integer (i-n)

        double precision dip_coef(0:3843)
        double precision dip_dc0(0:6,0:6), dip_dw0(0:6,0:6)

        integer i1

        common/NDIPCOE/ma,mb
        common/h5o2dipcoef/dip_dc0,dip_dw0,dip_coef

        ma=2512 ; mb=1332 

        open(20,file='h5o2.dms4B.coeff.com.dat',status='old')

        read(20,*)
        read(20,*)dip_dc0
        read(20,*)
        read(20,*)dip_dw0
        read(20,*)
        read(20,*)
        read(20,*)(dip_coef(i1),i1=0,ma+mb-1)
!        write(*,*)(dip_coef(i1),i1=0,ma+mb-1)
        close(20)

        return
        end

!*******************************************************************


        subroutine zundeldip(dip,cart_in)

        implicit none

        double precision dip(3), cart_in(7,3), dip_coef(0:3843)
        double precision dip_dc0(0:6,0:6), dip_dw0(0:6,0:6)

        integer ma,mb

        common/NDIPCOE/ma,mb
        common/h5o2dipcoef/dip_dc0,dip_dw0,dip_coef

        integer i,j

        double precision xnuc(0:2,0:6),vec(0:3,0:ma+mb-1)

        double precision dx(0:2),xn(0:2,0:6)

        do j=1,3
          xn(j-1,5:6)=cart_in(1:2,j)
          xn(j-1,0:4)=cart_in(3:7,j)
        end do

!! shift to center-of-charge frame (version for O2H5)
!        dx = (8*xn(0:2,0)+8*xn(0:2,1)+xn(0:2,2)+
!     $   xn(0:2,3)+xn(0:2,4)+xn(0:2,5)+xn(0:2,6))/21

!! shift to h5o2 C.O.M. frame (O2H5)
         dx=(15.9949146d0*xn(0:2,5)+15.9949146d0*xn(0:2,6)+
     $  1.007825d0*(xn(0:2,0)+xn(0:2,1)+xn(0:2,2)+xn(0:2,3)+
     $  xn(0:2,4)))/(15.9949146d0*2+1.007825d0*5)

        do i = 0, 6
         xnuc(0:2,i) = xn(0:2,i)-dx
        enddo

        call getdvec (ma, mb, xnuc(0:2,0:6),dip_dc0,dip_dw0,vec)
        dip = matmul(vec(0:2,0:ma+mb-1),dip_coef)

        dip = dip + dx

!        write(2002,*)dip

        return
        end subroutine
!********************************************************
!*********************************************************************                     
        subroutine getscale (n, xnuc, dc0, dw0)                                      
        implicit none                                                                
        integer, parameter :: wp=selected_real_kind(12,300)                          
! Version for X5Y2                                                           
        integer n, nk                                                                
        parameter (nk=7)                                                             
        double precision xnuc(0:2,0:nk-1,0:n-1), dc0(0:nk-1,0:nk-1),
     $   dw0(0:nk-1,0:nk-1)      
        integer ip, i, j                                                             
        double precision d0(0:nk-1,0:nk-1), dt0(0:nk-1,0:nk-1),
     $ r0(0:nk-1,0:nk-1), dc1(0:nk-1,0:nk-1), dw1(0:nk-1,0:nk-1),t0                            
! evaluate dc0                                                               
        dt0 = 0 ; dc1 = 0 ; dw1 = 1                                                  
        do ip = 0, n-1                                                               
         call getr0 (nk, xnuc(0:2,0:nk-1,ip), r0)                                    
         call getd0 (nk, r0, dc1, dw1, d0)                                           
         dt0 = dt0+d0/n                                                              
        enddo                                                                        
        do i = 0, nk-1                                                               
         dc0(i,i) = 0                                                                
        enddo                                                                        
        t0 = 0                                                                       
        do i = 0, 4                                                                  
         do j = i+1, 4                                                               
          t0 = t0+dt0(i,j)/10                                                        
         enddo                                                                       
        enddo                                                                        
        do i = 0, 4                                                                  
         do j = i+1, 4                                                               
          dc0(i,j) = t0                                                              
          dc0(j,i) = t0                                                              
         enddo                                                                       
        enddo                                                                        
        dc0(5,6) = dt0(5,6)                                                          
        dc0(6,5) = dt0(5,6)                                                          
        t0 = sum(dt0(0:4,5:6))/10                                                    
        dc0(0:4,5:6) = t0                                                            
        dc0(5:6,0:4) = t0                                                            
! evaluate dw0                                                               
        dt0 = 0                                                                      
        do ip = 0, n-1                                                               
         call getr0 (nk, xnuc(0:2,0:nk-1,ip), r0)                                    
         call getd0 (nk, r0, dc0, dw1, d0)                                           
         dt0 = dt0+d0**2/n                                                           
        enddo                                                                        
        t0 = 0                                                                       
        do i = 0, 4                                                                  
         t0 = t0+dt0(i,i)/5                                                          
        enddo                                                                        
        do i = 0, 4                                                                  
         dw0(i,i) = dsqrt(t0)                                                         
        enddo                                                                        
        do i = 5, 6                                                                  
         t0 = t0+dt0(i,i)/2                                                          
        enddo                                                                        
        do i = 5, 6                                                                  
         dw0(i,i) = dsqrt(t0)                                                         
        enddo                                                                        
        t0 = 0                                                                       
        do i = 0, 4                                                                  
         do j = i+1, 4                                                               
          t0 = t0+dt0(i,j)/10                                                        
         enddo                                                                       
        enddo                                                                        
        do i = 0, 4                                                                  
         do j = i+1, 4                                                               
          dw0(i,j) = dsqrt(t0)                                                        
          dw0(j,i) = dsqrt(t0)                                                        
         enddo                                                                       
        enddo                                                                        
        dw0(5,6) = dsqrt(dt0(5,6))                                                    
        dw0(6,5) = dsqrt(dt0(5,6))                                                    
        t0 = sum(dt0(0:4,5:6))/10                                                    
        dw0(0:4,5:6) = dsqrt(t0)                                                      
        dw0(5:6,0:4) = dsqrt(t0)                                                      
        return                                                                       
        end                                                                          
        subroutine getr0 (nk, xn, r0)                                                
        implicit none                                                                
        integer, parameter :: wp=selected_real_kind(12,300)                          
        integer nk                                                                   
        double precision xn(0:2,0:nk-1), r0(0:nk-1,0:nk-1)                             
        integer i, j                                                                 
        do i = 0, nk-1                                                               
         r0(i,i) = 0                                                                 
         do j = i+1, nk-1                                                            
          r0(i,j) = dsqrt((xn(0,j)-xn(0,i))**2+(xn(1,j)-xn(1,i))**2+                 
     $       (xn(2,j)-xn(2,i))**2)                                                    
          r0(j,i) = r0(i,j)                                                          
         enddo                                                                       
        enddo                                                                        
        return                                                                       
        end                                                                          
        subroutine getd0 (nk, r0, dc0, dw0, d0)                                      
        implicit none                                                                
        integer, parameter :: wp=selected_real_kind(12,300)                          
        integer nk                                                                   
        double precision r0(0:nk-1,0:nk-1), dc0(0:nk-1,0:nk-1),                       
     $     dw0(0:nk-1,0:nk-1), d0(0:nk-1,0:nk-1)                                      
        integer i, j                                                                 
        do i = 0, nk-1                                                               
         d0(i,i) = 0                                                                 
         do j = i+1, nk-1                                                            
!          d0(i,j) = (dlog(r0(i,j))-dc0(i,j))/dw0(i,j)                                 

           d0(i,j) = (dexp(-r0(i,j)/3.d0)-dc0(i,j))/dw0(i,j)

          d0(j,i) = d0(i,j)                                                          
         enddo                                                                       
        enddo                                                                        
        return                                                                       
        end                                                                          
        subroutine getdvec (ma, mb, xn, dc0, dw0, vec)                               
        implicit none                                                                
        integer, parameter :: wp=selected_real_kind(12,300)                          
! version for H5O2(+)                                                        
        integer nk, ma, mb                                                           
        parameter (nk=7)                                                             
        double precision xn(0:2,0:nk-1), dc0(0:nk-1,0:nk-1),       
     $  dw0(0:nk-1,0:nk-1),   vec(0:3,0:ma+mb-1)                                                         
        integer i, j, k, ind(0:nk-1)                                                 
        double precision r0(0:nk-1,0:nk-1), d0(0:nk-1,0:nk-1),     
     $  s0(0:nk-1,0:nk-1),   v0(0:ma+mb-1)                                                          
!-----------------------------------------------------------------------     
        vec = 0                                                                      
        call getr0 (nk, xn, r0)                                                      
        call getd0 (nk, r0, dc0, dw0, d0)                                            
! vector factor xn(0:2,0) (species X)                                        
        ind = (/ 1, 2, 3, 4, 5, 6, 0 /)                                              
        do i = 0, nk-1                                                               
         do j = 0, nk-1                                                              
          s0(i,j) = d0(ind(i),ind(j))                                                
         enddo                                                                       
        enddo                                                                        
        call getv_x4y2z (ma, s0, v0)                                                 
        do k = 0, ma-1                                                               
         vec(0:2,k) = vec(0:2,k)+xn(0:2,0)*v0(k)                                     
         vec(3,k) = vec(3,k)+v0(k)                                                   
        enddo                                                                        
! vector factor xn(0:2,1) (species X)                                        
        ind = (/ 0, 2, 3, 4, 5, 6, 1 /)                                              
        do i = 0, nk-1                                                               
         do j = 0, nk-1                                                              
          s0(i,j) = d0(ind(i),ind(j))                                                
         enddo                                                                       
        enddo                                                                        
        call getv_x4y2z (ma, s0, v0)                                                 
        do k = 0, ma-1                                                               
         vec(0:2,k) = vec(0:2,k)+xn(0:2,1)*v0(k)                                     
         vec(3,k) = vec(3,k)+v0(k)                                                   
        enddo                                                                        
! vector factor xn(0:2,2) (species X)                                        
        ind = (/ 0, 1, 3, 4, 5, 6, 2 /)                                              
        do i = 0, nk-1                                                               
         do j = 0, nk-1                                                              
          s0(i,j) = d0(ind(i),ind(j))                                                
         enddo                                                                       
        enddo                                                                        
        call getv_x4y2z (ma, s0, v0)                                                 
        do k = 0, ma-1                                                               
         vec(0:2,k) = vec(0:2,k)+xn(0:2,2)*v0(k)                                     
         vec(3,k) = vec(3,k)+v0(k)                                                   
        enddo                                                                        
! vector factor xn(0:2,3) (species X)                                        
        ind = (/ 0, 1, 2, 4, 5, 6, 3 /)                                              
        do i = 0, nk-1                                                               
         do j = 0, nk-1                                                              
          s0(i,j) = d0(ind(i),ind(j))                                                
         enddo                                                                       
        enddo                                                                        
        call getv_x4y2z (ma, s0, v0)                                                 
        do k = 0, ma-1                                                               
         vec(0:2,k) = vec(0:2,k)+xn(0:2,3)*v0(k)                                     
         vec(3,k) = vec(3,k)+v0(k)                                                   
        enddo                                                                        
! vector factor xn(0:2,4) (species X)                                        
        ind = (/ 0, 1, 2, 3, 5, 6, 4 /)                                              
        do i = 0, nk-1                                                               
         do j = 0, nk-1                                                              
          s0(i,j) = d0(ind(i),ind(j))                                                
         enddo                                                                       
        enddo                                                                        
        call getv_x4y2z (ma, s0, v0)                                                 
        do k = 0, ma-1                                                               
         vec(0:2,k) = vec(0:2,k)+xn(0:2,4)*v0(k)                                     
         vec(3,k) = vec(3,k)+v0(k)                                                   
        enddo                                                                        
! vector factor xn(0:2,5) (species Y)                                        
        ind = (/ 0, 1, 2, 3, 4, 6, 5 /)                                              
        do i = 0, nk-1                                                               
         do j = 0, nk-1                                                              
          s0(i,j) = d0(ind(i),ind(j))                                                
         enddo                                                                       
        enddo                                                                        
        call getv_x5yz (mb, s0, v0)                                                  
        do k = 0, mb-1                                                               
         vec(0:2,ma+k) = vec(0:2,ma+k)+xn(0:2,5)*v0(k)                               
         vec(3,ma+k) = vec(3,ma+k)+v0(k)                                             
        enddo                                                                        
! vector factor xn(0:2,6) (species Y)                                        
        s0 = d0                                                                      
        call getv_x5yz (mb, s0, v0)                                                  
        do k = 0, mb-1                                                               
         vec(0:2,ma+k) = vec(0:2,ma+k)+xn(0:2,6)*v0(k)                               
         vec(3,ma+k) = vec(3,ma+k)+v0(k)                                             
        enddo                                                                        
        return                                                                       
        end                                                                          
        subroutine getv_x4y2z (m, d, vec)                                            
        implicit none                                                                
        integer, parameter :: wp=selected_real_kind(12,300)                          
        integer, parameter :: nk=7                                                   
        integer m                                                                    
        double precision d(0:nk-1,0:nk-1), vec(0:m-1)                                  
! version for molecule X4Y2Z                                                 
! MolienSeries(0:*):                                                         
! #Primaries(1:*):                                                           
! #Secondaries(1:*):                                                         
        integer, parameter :: l0=1, l1=l0+5, l2=l1+26, l3=l2+117, 
     $ l4=l3+501, l5=l4+1975-113, l6=-1, l7=-2                                               
!! Note, we have included only 86 of the 199 secondaries at degree 5.        
        integer, parameter :: np1=5, np2=7, np3=4, np4=4, np5=0, np6=1,            
     $     np7=0                                                                      
        integer, parameter :: ns1=0, ns2=4, ns3=23, ns4=71, ns5=199,                
     $     ns6=548, ns7=1323                                                          
        double precision x(0:np1-1), y(0:np2-1), z(0:np3-1), u(0:np4-1),             
     $     w(0:np6-1)                                        
       double precision x2(0:np1-1),x3(0:np1-1),x4(0:np1-1),x5(0:np1-1),          
     $  x6(0:np1-1), x7(0:np1-1),  y2(0:np2-1), y3(0:np2-1),z2(0:np2-1)                                      
        double precision ys(0:ns2-1), zs(0:ns3-1), 
     $  us(0:ns4-1), vs(0:ns5-1), ws(0:ns6-1), w7s(0:ns7-1)                                     
        integer mdeg, i, j, i0, i1, i2, j0, j1, k0                             
        double precision d2(0:nk-1,0:nk-1), d3(0:nk-1,0:nk-1), 
     $  d4(0:nk-1,0:nk-1), d5(0:nk-1,0:nk-1), d6(0:nk-1,0:nk-1), 
     $  d7(0:nk-1,0:nk-1)                    
        double precision t0                                            
        double precision her2, her3, her4, her5, her6, her7 
!        her2(t0) = (4*t0**2-2)/sqrt(dble(8*2))
!        her3(t0) = (8*t0**2-12)*t0/sqrt(dble(16*6))
!        her4(t0) = ((16*t0**2-48)*t0**2+12)/sqrt(dble(32*24))
!        her5(t0) = ((32*t0**2-160)*t0**2+120)*t0/sqrt(dble(64*120))
!        her6(t0) = (((64*t0**2-480)*t0**2+720)*t0**2-120)/
!     $  sqrt(dble(128*720))
!        her7(t0) = (((128*t0**2-1344)*t0**2+3360)*t0**2-1680)*t0/
!     $  sqrt(dble(256*5040))
!-----------------------------------------------------------------------     
!! We don't have the secondaries at degrees 5 and up yet.                    
!-----------------------------------------------------------------------     
! Test for compatibility, set mdeg                                           
        select case (m)                                                              
        case (l0)                                                                    
         mdeg = 0                                                                    
        case (l1)                                                                    
         mdeg = 1                                                                    
        case (l2)                                                                    
         mdeg = 2                                                                    
        case (l3)                                                                    
         mdeg = 3                                                                    
        case (l4)                                                                    
         mdeg = 4                                                                    
        case (l5)                                                                    
         mdeg = 5                                                                    
        case (l6)                                                                    
         mdeg = 6                                                                    
        case (l7)                                                                    
         mdeg = 7                                                                    
        case default                                                                 
         stop 'getv - wrong dimension'                                               
        endselect                                                                    
! auxiliary distances                                                        
        do i = 0, nk-1                                                               
         do j = i+1, nk-1                                                            
          d2(i,j) = her2(d(i,j))                                                     
          d2(j,i) = d2(i,j)                                                          
          d3(i,j) = her3(d(i,j))                                                     
          d3(j,i) = d3(i,j)                                                          
          d4(i,j) = her4(d(i,j))                                                     
          d4(j,i) = d4(i,j)                                                          
          d5(i,j) = her5(d(i,j))                                                     
          d5(j,i) = d5(i,j)                                                          
          d6(i,j) = her6(d(i,j))                                                     
          d6(j,i) = d6(i,j)                                                          
          d7(i,j) = her7(d(i,j))                                                     
          d7(j,i) = d7(i,j)                                                          
         enddo                                                                       
        enddo                                                                        
! Primary Invariants                                                         
        x = 0 ; y = 0 ; z = 0 ; u = 0 ;  w = 0                                
        do i0 = 0, 3                                                                 
         t0 = 0                                                                      
         do i1 = 0, 3                                                                
         if (i1.ne.i0) then                                                          
          t0 = t0+d(i0,i1)/3                                                         
          x(0) = x(0)+d(i0,i1)/6                                                     
          y(0) = y(0)+d2(i0,i1)/6                                                    
          z(0) = z(0)+d3(i0,i1)/6                                                    
         endif                                                                       
         enddo                                                                       
         y(1) = y(1)+her2(t0)/4                                                      
         z(1) = z(1)+her3(t0)/4                                                      
         u(0) = u(0)+her4(t0)/4                                                      
        enddo                                                                        
        x(1) = sum(d(0:3,4:5))/8                                                     
        do i0 = 0, 3                                                                 
         y(2) = y(2)+her2(sum(d(i0,4:5))/2)/4                                        
         z(2) = z(2)+her3(sum(d(i0,4:5))/2)/4                                        
         u(1) = u(1)+her4(sum(d(i0,4:5))/2)/4                                        
        enddo                                                                        
        do j0 = 4, 5                                                                 
         y(3) = y(3)+her2(sum(d(0:3,j0))/4)/2                                        
        enddo                                                                        
        y(4) = sum(d2(0:3,4:5))/8                                                    
        u(2) = sum(d4(0:3,4:5))/8                                                    
        w(0) = sum(d6(0:3,4:5))/8                                                    
        x(2) = d(4,5)                                                                
        x(3) = sum(d(0:3,6))/4                                                       
        y(5) = sum(d2(0:3,6))/4                                                      
        z(3) = sum(d3(0:3,6))/4                                                      
        u(3) = sum(d4(0:3,6))/4                                                      
        x(4) = sum(d(4:5,6))/2                                                       
        y(6) = sum(d2(4:5,6))/2                                                      
! Required powers                                                            
        do i = 0, np1-1                                                              
         x2(i) = her2(x(i))                                                          
         x3(i) = her3(x(i))                                                          
         x4(i) = her4(x(i))                                                          
         x5(i) = her5(x(i))                                                          
         x6(i) = her6(x(i))                                                          
         x7(i) = her7(x(i))                                                          
        enddo                                                                        
        do i = 0, np2-1                                                              
         y2(i) = her2(y(i))                                                          
         y3(i) = her3(y(i))                                                          
        enddo                                                                        
        do i = 0, np3-1                                                              
         z2(i) = her2(z(i))                                                          
        enddo                                                                        
! Secondary Invariants                                                       
!! reducible: us(0:8), vs(0:85), ws(0:447)                                   
!! Note: We haven't included them all yet                                    
        ys = 0 ; zs = 0 ; us = 0 ; vs = 0 ; ws = 0 ; w7s = 0                         
! Irreducible secondaries                                                    
        do i0 = 0, 3                                                                 
         do i1 = 0, 3                                                                
         if (i1.ne.i0) then                                                          
          us(9) = us(9)+d4(i0,i1)/12                                                 
          do i2 = 0, 3                                                               
          if (i2.ne.i1.and.i2.ne.i0) then                                            
           zs(0) = zs(0)+d2(i0,i1)*d(i0,i2)/24                                       
          endif                                                                      
          enddo                                                                      
         endif                                                                       
         enddo                                                                       
        enddo                                                                        
        do j0 = 4, 5                                                                 
         do i0 = 0, 3                                                                
          zs(1) = zs(1)+d3(i0,j0)/8                                                  
          do i1 = 0, 3                                                               
          if (i1.ne.i0) then                                                         
           ys(0) = ys(0)+d(i0,i1)*d(i0,j0)/24                                        
           zs(2) = zs(2)+d2(i0,i1)*d(i0,j0)/24                                       
           zs(3) = zs(3)+d(i0,i1)*d2(i0,j0)/24                                       
           zs(4) = zs(4)+d(i0,i1)*d(i0,j0)*d(i1,j0)/24                               
           zs(5) = zs(5)+d2(i0,j0)*d(i1,j0)/24                                       
           us(10) = us(10)+d3(i0,i1)*d(i0,j0)/24                                     
           us(11) = us(11)+d2(i0,i1)*d2(i0,j0)/24                                    
           us(12) = us(12)+d(i0,i1)*d3(i0,j0)/24                                     
           us(13) = us(13)+d2(i0,i1)*d(i0,j0)*d(i1,j0)/24                            
           us(14) = us(14)+d(i0,i1)*d2(i0,j0)*d(i1,j0)/24                            
           us(15) = us(15)+d3(i0,j0)*d(i1,j0)/24                                     
           us(16) = us(16)+d2(i0,j0)*d2(i1,j0)/24                                    
           do i2 = 0, 3                                                              
           if (i2.ne.i1.and.i2.ne.i0) then                                           
            zs(6) = zs(6)+d(i0,i1)*d(i0,i2)*d(i0,j0)/48                              
            zs(7) = zs(7)+d(i0,i2)*d(i0,j0)*d(i1,j0)/48                              
            us(17) = us(17)+d2(i0,i1)*d(i0,i2)*d(i0,j0)/48                           
            us(18) = us(18)+d2(i0,i1)*d(i1,i2)*d(i0,j0)/48                           
            us(19) = us(19)+d(i0,i1)*d(i0,i2)*d2(i0,j0)/48                           
            us(20) = us(20)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)/48                   
            us(21) = us(21)+d2(i0,i2)*d(i0,j0)*d(i1,j0)/48                           
            us(22) = us(22)+d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,j0)/48                   
            us(23) = us(23)+d(i0,i2)*d2(i0,j0)*d(i1,j0)/48                           
            us(24) = us(24)+d(i1,i2)*d2(i0,j0)*d(i1,j0)/48                           
           endif                                                                     
           enddo                                                                     
          endif                                                                      
          enddo                                                                      
         enddo                                                                       
        enddo                                                                        
        do j0 = 4, 5                                                                 
         do j1 = 4, 5                                                                
         if (j1.ne.j0) then                                                          
          do i0 = 0, 3                                                               
           us(25) = us(25)+d3(i0,j0)*d(i0,j1)/8                                      
           do i1 = 0, 3                                                              
           if (i1.ne.i0) then                                                        
            zs(8) = zs(8)+d(i0,i1)*d(i0,j0)*d(i0,j1)/24                              
            zs(9) = zs(9)+d(i0,i1)*d(i1,j0)*d(i0,j1)/24                              
            us(26) = us(26)+d2(i0,i1)*d(i0,j0)*d(i0,j1)/24                           
            us(27) = us(27)+d(i0,i1)*d2(i0,j0)*d(i0,j1)/24                           
            us(28) = us(28)+d2(i0,i1)*d(i1,j0)*d(i0,j1)/24                           
            us(29) = us(29)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,j1)/24                   
            us(30) = us(30)+d2(i0,j0)*d(i1,j0)*d(i0,j1)/24                           
            do i2 = 0, 3                                                             
            if (i2.ne.i1.and.i2.ne.i0) then                                          
             us(31) = us(31)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,j1)/48                  
            endif                                                                    
            enddo                                                                    
           endif                                                                     
           enddo                                                                     
          enddo                                                                      
         endif                                                                       
         enddo                                                                       
        enddo                                                                        
        k0 = 6                                                                       
        do i0 = 0, 3                                                                 
         do i1 = 0, 3                                                                
         if (i1.ne.i0) then                                                          
          ys(1) = ys(1)+d(i0,i1)*d(i0,k0)/12                                         
          zs(10) = zs(10)+d2(i0,i1)*d(i0,k0)/12                                      
          zs(11) = zs(11)+d(i0,i1)*d2(i0,k0)/12                                      
          zs(12) = zs(12)+d(i0,i1)*d(i0,k0)*d(i1,k0)/12                              
          us(32) = us(32)+d3(i0,i1)*d(i0,k0)/12                                      
          us(33) = us(33)+d2(i0,i1)*d2(i0,k0)/12                                     
          us(34) = us(34)+d(i0,i1)*d3(i0,k0)/12                                      
          us(35) = us(35)+d2(i0,i1)*d(i0,k0)*d(i1,k0)/12                             
          do i2 = 0, 3                                                               
          if (i2.ne.i1.and.i2.ne.i0) then                                            
           zs(13) = zs(13)+d(i0,i1)*d(i0,i2)*d(i0,k0)/24                             
           us(36) = us(36)+d2(i0,i1)*d(i0,i2)*d(i0,k0)/24                            
           us(37) = us(37)+d2(i0,i1)*d(i1,i2)*d(i0,k0)/24                            
           us(38) = us(38)+d(i0,i1)*d(i0,i2)*d2(i0,k0)/24                            
          endif                                                                      
          enddo                                                                      
         endif                                                                       
         enddo                                                                       
        enddo                                                                        
        do j0 = 4, 5                                                                 
         do i0 = 0, 3                                                                
          ys(2) = ys(2)+d(i0,j0)*d(i0,k0)/8                                          
          ys(3) = ys(3)+d(i0,j0)*d(j0,k0)/8                                          
          zs(14) = zs(14)+d2(i0,j0)*d(i0,k0)/8                                       
          zs(15) = zs(15)+d(i0,j0)*d2(i0,k0)/8                                       
          zs(16) = zs(16)+d(i0,j0)*d(i0,k0)*d(j0,k0)/8                               
          zs(17) = zs(17)+d2(i0,j0)*d(j0,k0)/8                                       
          us(39) = us(39)+d2(i0,j0)*d2(i0,k0)/8                                      
          us(40) = us(40)+d(i0,j0)*d3(i0,k0)/8                                       
          us(41) = us(41)+d3(i0,j0)*d(i0,k0)/8                                       
          us(42) = us(42)+d3(i0,j0)*d(j0,k0)/8                                       
          us(43) = us(43)+d2(i0,j0)*d(i0,k0)*d(j0,k0)/8                              
          us(44) = us(44)+d(i0,j0)*d2(i0,k0)*d(j0,k0)/8                              
          do i1 = 0, 3                                                               
          if (i1.ne.i0) then                                                         
           zs(18) = zs(18)+d(i0,i1)*d(i0,j0)*d(i0,k0)/24                             
           zs(19) = zs(19)+d(i0,i1)*d(i1,j0)*d(i0,k0)/24                             
           zs(20) = zs(20)+d(i0,j0)*d(i1,j0)*d(i0,k0)/24                             
           zs(21) = zs(21)+d(i0,i1)*d(i0,j0)*d(j0,k0)/24                             
           us(45) = us(45)+d2(i0,i1)*d(i0,j0)*d(i0,k0)/24                            
           us(46) = us(46)+d(i0,i1)*d2(i0,j0)*d(i0,k0)/24                            
           us(47) = us(47)+d2(i0,i1)*d(i1,j0)*d(i0,k0)/24                            
           us(48) = us(48)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,k0)/24                    
           us(49) = us(49)+d2(i0,j0)*d(i1,j0)*d(i0,k0)/24                            
           us(50) = us(50)+d(i0,i1)*d2(i1,j0)*d(i0,k0)/24                            
           us(51) = us(51)+d(i0,j0)*d2(i1,j0)*d(i0,k0)/24                            
           us(52) = us(52)+d(i0,i1)*d(i0,j0)*d2(i0,k0)/24                            
           us(53) = us(53)+d(i0,i1)*d(i1,j0)*d2(i0,k0)/24                            
           us(54) = us(54)+d(i0,j0)*d(i1,j0)*d2(i0,k0)/24                            
           us(55) = us(55)+d(i0,j0)*d(i1,j0)*d(i0,k0)*d(i1,k0)/24                    
           us(56) = us(56)+d2(i0,i1)*d(i0,j0)*d(j0,k0)/24                            
           us(57) = us(57)+d(i0,i1)*d2(i0,j0)*d(j0,k0)/24                            
           us(58) = us(58)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(j0,k0)/24                    
           us(59) = us(59)+d(i0,i1)*d(i0,j0)*d(i0,k0)*d(j0,k0)/24                    
           us(60) = us(60)+d(i0,i1)*d(i1,j0)*d(i0,k0)*d(j0,k0)/24                    
           do i2 = 0, 3                                                              
           if (i2.ne.i1.and.i2.ne.i0) then                                           
            us(61) = us(61)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,k0)/48                   
            us(62) = us(62)+d(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,k0)/48                   
            us(63) = us(63)+d(i0,i2)*d(i0,j0)*d(i1,j0)*d(i0,k0)/48                   
            us(64) = us(64)+d(i1,i2)*d(i0,j0)*d(i1,j0)*d(i0,k0)/48                   
            us(65) = us(65)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(j0,k0)/48                   
           endif                                                                     
           enddo                                                                     
          endif                                                                      
          enddo                                                                      
         enddo                                                                       
        enddo                                                                        
        do j0 = 4, 5                                                                 
         do j1 = 4, 5                                                                
         if (j1.ne.j0) then                                                          
          do i0 = 0, 3                                                               
           zs(22) = zs(22)+d(i0,j0)*d(i0,j1)*d(i0,k0)/8                              
           us(66) = us(66)+d2(i0,j0)*d(i0,j1)*d(i0,k0)/8                             
           us(67) = us(67)+d(i0,j0)*d(i0,j1)*d2(i0,k0)/8                             
           us(68) = us(68)+d2(i0,j0)*d(i0,j1)*d(j0,k0)/8                             
           do i1 = 0, 3                                                              
           if (i1.ne.i0) then                                                        
            us(69) = us(69)+d(i0,i1)*d(i0,j0)*d(i0,j1)*d(i0,k0)/24                   
            us(70) = us(70)+d(i0,i1)*d(i1,j0)*d(i0,j1)*d(i0,k0)/24                   
            do i2 = 0, 3                                                             
            if (i2.ne.i1.and.i2.ne.i0) then                                          
            endif                                                                    
            enddo                                                                    
           endif                                                                     
           enddo                                                                     
          enddo                                                                      
         endif                                                                       
         enddo                                                                       
        enddo                                                                        
! Reducible secondaries                                                      
!! Should be us(0:8), vs(0:85), ws(0:447)                                    
!! Note. ws(0:*) is incomplete.                                              
        us(0) = her2(ys(0))                                                          
        us(1) = ys(0)*ys(1)                                                          
        us(2) = her2(ys(1))                                                          
        us(3:4) = ys(0:1)*ys(2)                                                      
        us(5) = her2(ys(2))                                                          
        us(6:8) = ys(0:2)*ys(3)                                                      
! her2(ys(3)) is not independent                                             
        vs(0:3) = ys(0:3)*zs(0)                                                      
        vs(4:7) = ys(0:3)*zs(1)                                                      
        vs(8:11) = ys(0:3)*zs(2)                                                     
        vs(12:15) = ys(0:3)*zs(3)                                                    
        vs(16:19) = ys(0:3)*zs(4)                                                    
        vs(20:22) = ys(0:2)*zs(5)                                                    
        vs(23:26) = ys(0:3)*zs(6)                                                    
        vs(27:29) = ys(0:2)*zs(7)                                                    
        vs(30:33) = ys(0:3)*zs(8)                                                    
        vs(34:37) = ys(0:3)*zs(9)                                                    
        vs(38:41) = ys(0:3)*zs(10)                                                   
        vs(42:45) = ys(0:3)*zs(11)                                                   
        vs(46:49) = ys(0:3)*zs(12)                                                   
        vs(50:53) = ys(0:3)*zs(13)                                                   
        vs(54:57) = ys(0:3)*zs(14)                                                   
        vs(58:61) = ys(0:3)*zs(15)                                                   
        vs(62:64) = ys(0:2)*zs(16)                                                   
        vs(65:67) = ys(0:2)*zs(17)                                                   
        vs(68:71) = ys(0:3)*zs(18)                                                   
        vs(72:75) = ys(0:3)*zs(19)                                                   
        vs(76:78) = ys(0:2)*zs(20)                                                   
        vs(79:81) = ys(0:2)*zs(21)                                                   
        vs(82:85) = ys(0:3)*zs(22)                                                   
! Compute vec(0:*).  The code was created using these parameters:            
! MolienSeries(0:*): 1 5 26 117 501 1975-113 ...                             
! #Primaries(1:*):   5 7 4 4 0 1 0 0                                         
! #Secondaries(1:*): 0 4 23 71 199 548 1323 1                                
! The MolienSeries partial sums (allowed size of vec) are:                   
! Molien Sums(0:*):  1 6 32 149 650 2625-113 ...                             
! constant term                                                              
        vec(0) = 1                                                                   
! first degree terms                                                         
        if (1.le.mdeg) then                                                          
         vec(1) = x(0)                                                               
         vec(2) = x(1)                                                               
         vec(3) = x(2)                                                               
         vec(4) = x(3)                                                               
         vec(5) = x(4)                                                               
        endif                                                                        
! second degree terms                                                        
        if (2.le.mdeg) then                                                          
         vec(6) = x2(0)                                                              
         vec(7:10) = x(0)*vec(2:5)                                                   
         vec(11) = x2(1)                                                             
         vec(12:14) = x(1)*vec(3:5)                                                  
         vec(15) = x2(2)                                                             
         vec(16:17) = x(2)*vec(4:5)                                                  
         vec(18) = x2(3)                                                             
         vec(19:19) = x(3)*vec(5:5)                                                  
         vec(20) = x2(4)                                                             
         vec(21:20) = x(4)*vec(6:5)                                                  
         vec(21) = y(0)                                                              
         vec(22) = y(1)                                                              
         vec(23) = y(2)                                                              
         vec(24) = y(3)                                                              
         vec(25) = y(4)                                                              
         vec(26) = y(5)                                                              
         vec(27) = y(6)                                                              
         vec(28:31) = ys(0:3)                                                        
        endif                                                                        
! third degree terms                                                         
        if (3.le.mdeg) then                                                          
         vec(32) = x3(0)                                                             
         vec(33:36) = x2(0)*vec(2:5)                                                 
         vec(37:57) = x(0)*vec(11:31)                                                
         vec(58) = x3(1)                                                             
         vec(59:61) = x2(1)*vec(3:5)                                                 
         vec(62:78) = x(1)*vec(15:31)                                                
         vec(79) = x3(2)                                                             
         vec(80:81) = x2(2)*vec(4:5)                                                 
         vec(82:95) = x(2)*vec(18:31)                                                
         vec(96) = x3(3)                                                             
         vec(97:97) = x2(3)*vec(5:5)                                                 
         vec(98:109) = x(3)*vec(20:31)                                               
         vec(110) = x3(4)                                                            
         vec(111:110) = x2(4)*vec(6:5)                                               
         vec(111:121) = x(4)*vec(21:31)                                              
         vec(122) = z(0)                                                             
         vec(123) = z(1)                                                             
         vec(124) = z(2)                                                             
         vec(125) = z(3)                                                             
         vec(126:148) = zs(0:22)                                                     
        endif                                                                        
! fourth degree terms                                                        
        if (4.le.mdeg) then                                                          
         vec(149) = x4(0)                                                            
         vec(150:153) = x3(0)*vec(2:5)                                               
         vec(154:174) = x2(0)*vec(11:31)                                             
         vec(175:265) = x(0)*vec(58:148)                                             
         vec(266) = x4(1)                                                            
         vec(267:269) = x3(1)*vec(3:5)                                               
         vec(270:286) = x2(1)*vec(15:31)                                             
         vec(287:356) = x(1)*vec(79:148)                                             
         vec(357) = x4(2)                                                            
         vec(358:359) = x3(2)*vec(4:5)                                               
         vec(360:373) = x2(2)*vec(18:31)                                             
         vec(374:426) = x(2)*vec(96:148)                                             
         vec(427) = x4(3)                                                            
         vec(428:428) = x3(3)*vec(5:5)                                               
         vec(429:440) = x2(3)*vec(20:31)                                             
         vec(441:479) = x(3)*vec(110:148)                                            
         vec(480) = x4(4)                                                            
         vec(481:480) = x3(4)*vec(6:5)                                               
         vec(481:491) = x2(4)*vec(21:31)                                             
         vec(492:518) = x(4)*vec(122:148)                                            
         vec(519) = y2(0)                                                            
         vec(520:529) = y(0)*vec(22:31)                                              
         vec(530) = y2(1)                                                            
         vec(531:539) = y(1)*vec(23:31)                                              
         vec(540) = y2(2)                                                            
         vec(541:548) = y(2)*vec(24:31)                                              
         vec(549) = y2(3)                                                            
         vec(550:556) = y(3)*vec(25:31)                                              
         vec(557) = y2(4)                                                            
         vec(558:563) = y(4)*vec(26:31)                                              
         vec(564) = y2(5)                                                            
         vec(565:569) = y(5)*vec(27:31)                                              
         vec(570) = y2(6)                                                            
         vec(571:574) = y(6)*vec(28:31)                                              
         vec(575) = u(0)                                                             
         vec(576) = u(1)                                                             
         vec(577) = u(2)                                                             
         vec(578) = u(3)                                                             
         vec(579:649) = us(0:70)                                                     
        endif                                                                        
! fifth degree terms                                                         
        if (5.le.mdeg) then                                                          
         vec(650) = x5(0)                                                            
         vec(651:654) = x4(0)*vec(2:5)                                               
         vec(655:675) = x3(0)*vec(11:31)                                             
         vec(676:766) = x2(0)*vec(58:148)                                            
         vec(767:1150) = x(0)*vec(266:649)                                           
         vec(1151) = x5(1)                                                           
         vec(1152:1154) = x4(1)*vec(3:5)                                             
         vec(1155:1171) = x3(1)*vec(15:31)                                           
         vec(1172:1241) = x2(1)*vec(79:148)                                          
         vec(1242:1534) = x(1)*vec(357:649)                                          
         vec(1535) = x5(2)                                                           
         vec(1536:1537) = x4(2)*vec(4:5)                                             
         vec(1538:1551) = x3(2)*vec(18:31)                                           
         vec(1552:1604) = x2(2)*vec(96:148)                                          
         vec(1605:1827) = x(2)*vec(427:649)                                          
         vec(1828) = x5(3)                                                           
         vec(1829:1829) = x4(3)*vec(5:5)                                             
         vec(1830:1841) = x3(3)*vec(20:31)                                           
         vec(1842:1880) = x2(3)*vec(110:148)                                         
         vec(1881:2050) = x(3)*vec(480:649)                                          
         vec(2051) = x5(4)                                                           
         vec(2052:2051) = x4(4)*vec(6:5)                                             
         vec(2052:2062) = x3(4)*vec(21:31)                                           
         vec(2063:2089) = x2(4)*vec(122:148)                                         
         vec(2090:2220) = x(4)*vec(519:649)                                          
         vec(2221:2247) = y(0)*vec(122:148)                                          
         vec(2248:2274) = y(1)*vec(122:148)                                          
         vec(2275:2301) = y(2)*vec(122:148)                                          
         vec(2302:2328) = y(3)*vec(122:148)                                          
         vec(2329:2355) = y(4)*vec(122:148)                                          
         vec(2356:2382) = y(5)*vec(122:148)                                          
         vec(2383:2409) = y(6)*vec(122:148)                                          
         vec(2410:2413) = z(0)*ys(0:3)                                               
         vec(2414:2417) = z(1)*ys(0:3)                                               
         vec(2418:2421) = z(2)*ys(0:3)                                               
         vec(2422:2425) = z(3)*ys(0:3)                                               
         vec(2426:2511) = vs(0:85)                                                   
!! should be vec(2426:2624) = vs(0:198)                                      
        endif                                                                        
        return                                                                       
        end                                                                          
        subroutine getv_x5yz (m, d, vec)                                             
        implicit none                                                                
        integer, parameter :: wp=selected_real_kind(12,300)                          
        integer, parameter :: nk=7                                                   
        integer m                                                                    
        double precision d(0:nk-1,0:nk-1), vec(0:m-1)                                  
! version for molecule X5YZ.                                                 
! MolienSeries(0:*): 1 4 17 68 265 977 3424 11361                            
! #Primaries(1:*):   4 4 4 4 4 1 0 0 ...                                     
! #Secondaries(1:*): 0 3 16 54 165 482 1303                                  
        integer, parameter :: l0=1, l1=l0+4, l2=l1+17, l3=l2+68, 
     $ l4=l3+265, l5=l4+977, l6=l5+3424-482, l7=-1                                           
!! Incomplete at degree 6, and we haven't done degree 7.                     
        integer, parameter :: np1=4, np2=4, np3=4, np4=4, np5=4, np6=1,            
     $     np7=0                                                                      
        integer, parameter :: ns1=0, ns2=3, ns3=16, ns4=54, ns5=165,                
     $     ns6=482, ns7=1303                                                          
        double precision x(0:np1-1), y(0:np2-1), z(0:np3-1), u(0:np4-1),              
     $     v(0:np5-1), w(0:np6-1)                                        
        double precision x2(0:np1-1), x3(0:np1-1), x4(0:np1-1), 
     $ x5(0:np1-1), x6(0:np1-1), x7(0:np1-1),y2(0:np2-1), y3(0:np2-1), 
     $ z2(0:np2-1)                                      
        double precision ys(0:ns2-1), zs(0:ns3-1), 
     $ us(0:ns4-1), vs(0:ns5-1), ws(0:ns6-1), w7s(0:ns7-1)                                     
        integer mdeg, i, j, k, i0, i1, i2, i3, j0, k0                             
        double precision d2(0:nk-1,0:nk-1),d3(0:nk-1,0:nk-1),d4(0:nk-1,
     $ 0:nk-1), d5(0:nk-1,0:nk-1), d6(0:nk-1,0:nk-1), d7(0:nk-1,0:nk-1)                    
        double precision t0                                                            
        double precision her2, her3, her4, her5, her6, her7                            
!        her2(t0) = (4*t0**2-2)/dsqrt(dble(8*2))                               
!        her3(t0) = (8*t0**2-12)*t0/dsqrt(dble(16*6))                          
!        her4(t0) = ((16*t0**2-48)*t0**2+12)/dsqrt(dble(32*24))                
!        her5(t0) = ((32*t0**2-160)*t0**2+120)*t0/dsqrt(dble(64*120))          
!        her6(t0) = (((64*t0**2-480)*t0**2+720)*t0**2-120)/                          
!     $     dsqrt(dble(128*720))                                                
!        her7(t0) = (((128*t0**2-1344)*t0**2+3360)*t0**2-1680)*t0/                   
!     $     dsqrt(dble(256*5040))                                               
!-----------------------------------------------------------------------     
!! We don't have the secondaries at degree 6 yet.                            
!-----------------------------------------------------------------------     
! Test for compatibility, set mdeg                                           
        select case (m)                                                              
        case (l0)                                                                    
         mdeg = 0                                                                    
        case (l1)                                                                    
         mdeg = 1                                                                    
        case (l2)                                                                    
         mdeg = 2                                                                    
        case (l3)                                                                    
         mdeg = 3                                                                    
        case (l4)                                                                    
         mdeg = 4                                                                    
        case (l5)                                                                    
         mdeg = 5                                                                    
        case (l6)                                                                    
         mdeg = 6                                                                    
        case (l7)                                                                    
         mdeg = 7                                                                    
        case default                                                                 
         stop 'getv - wrong dimension'                                               
        endselect                                                                    
! auxiliary distances                                                        
        do i = 0, nk-1                                                               
         do j = i+1, nk-1                                                            
          d2(i,j) = her2(d(i,j))                                                     
          d2(j,i) = d2(i,j)                                                          
          d3(i,j) = her3(d(i,j))                                                     
          d3(j,i) = d3(i,j)                                                          
          d4(i,j) = her4(d(i,j))                                                     
          d4(j,i) = d4(i,j)                                                          
          d5(i,j) = her5(d(i,j))                                                     
          d5(j,i) = d5(i,j)                                                          
          d6(i,j) = her6(d(i,j))                                                     
          d6(j,i) = d6(i,j)                                                          
          d7(i,j) = her7(d(i,j))                                                     
          d7(j,i) = d7(i,j)                                                          
         enddo                                                                       
        enddo                                                                        
! Primary Invariants                                                         
        x = 0 ; y = 0 ; z = 0 ; u = 0 ; v = 0 ; w = 0                                
        do i0 = 0, 4                                                                 
         t0 = 0                                                                      
         do i1 = 0, 4                                                                
         if (i1.ne.i0) then                                                          
          t0 = t0+d(i0,i1)/4                                                         
          x(0) = x(0)+d(i0,i1)/20                                                    
          y(0) = y(0)+d2(i0,i1)/20                                                   
          z(0) = z(0)+d3(i0,i1)/20                                                   
          u(0) = u(0)+d4(i0,i1)/20                                                   
          v(0) = v(0)+d5(i0,i1)/20                                                   
          w(0) = w(0)+d6(i0,i1)/20                                                   
         endif                                                                       
         enddo                                                                       
         y(1) = y(1)+her2(t0)/5                                                      
         z(1) = z(1)+her3(t0)/5                                                      
         u(1) = u(1)+her4(t0)/5                                                      
         v(1) = v(1)+her5(t0)/5                                                      
        enddo                                                                        
        x(1) = sum(d(0:4,5))/5                                                       
        y(2) = sum(d2(0:4,5))/5                                                      
        z(2) = sum(d3(0:4,5))/5                                                      
        u(2) = sum(d4(0:4,5))/5                                                      
        v(2) = sum(d5(0:4,5))/5                                                      
        x(2) = sum(d(0:4,6))/5                                                       
        y(3) = sum(d2(0:4,6))/5                                                      
        z(3) = sum(d3(0:4,6))/5                                                      
        u(3) = sum(d4(0:4,6))/5                                                      
        v(3) = sum(d5(0:4,6))/5                                                      
        x(3) = d(5,6)                                                                
! Required powers                                                            
        do i = 0, np1-1                                                              
         x2(i) = her2(x(i))                                                          
         x3(i) = her3(x(i))                                                          
         x4(i) = her4(x(i))                                                          
         x5(i) = her5(x(i))                                                          
         x6(i) = her6(x(i))                                                          
         x7(i) = her7(x(i))                                                          
        enddo                                                                        
        do i = 0, np2-1                                                              
         y2(i) = her2(y(i))                                                          
         y3(i) = her3(y(i))                                                          
        enddo                                                                        
        do i = 0, np3-1                                                              
         z2(i) = her2(z(i))                                                          
        enddo                                                                        
! Secondary Invariants                                                       
!! ys(0:2), zs(0:15), us(0:53), vs(0:164), ws(0:481), w7s(0:1302)            
!! reducible: us(0:5), vs(0:47), ws(0:..)                                    
!! Note: at this time the secondaries at degree 6 are incomplete.            
        ys = 0 ; zs = 0 ; us = 0 ; vs = 0 ; ws = 0 ; w7s = 0                         
! Irreducible secondaries                                                    
        do i0 = 0, 4                                                                 
         do i1 = 0, 4                                                                
         if (i1.ne.i0) then                                                          
          do i2 = 0, 4                                                               
          if (i2.ne.i1.and.i2.ne.i0) then                                            
           zs(0) = zs(0)+d2(i0,i1)*d(i0,i2)/60                                       
           zs(1) = zs(1)+d(i0,i1)*d(i0,i2)*d(i1,i2)/60                               
           us(6) = us(6)+d3(i0,i1)*d(i0,i2)/60                                       
           us(7) = us(7)+d2(i0,i1)*d2(i0,i2)/60                                      
           us(8) = us(8)+d2(i0,i1)*d(i0,i2)*d(i1,i2)/60                              
           vs(48) = vs(48)+d4(i0,i1)*d(i0,i2)/60                                     
           vs(49) = vs(49)+d3(i0,i1)*d2(i0,i2)/60                                    
           vs(50) = vs(50)+d3(i0,i1)*d(i0,i2)*d(i1,i2)/60                            
           vs(51) = vs(51)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)/60                           
           do i3 = 0, 4                                                              
           if (i3.ne.i0.and.i3.ne.i1.and.i3.ne.i2) then                              
            us(9) = us(9)+d2(i0,i1)*d(i0,i2)*d(i0,i3)/120                            
            us(10) = us(10)+d2(i0,i1)*d(i1,i2)*d(i0,i3)/120                          
            vs(52) = vs(52)+d3(i0,i1)*d(i0,i2)*d(i0,i3)/120                          
            vs(53) = vs(53)+d2(i0,i1)*d2(i0,i2)*d(i0,i3)/120                         
            vs(54) = vs(54)+d3(i0,i1)*d(i1,i2)*d(i0,i3)/120                          
            vs(55) = vs(55)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)/120                 
           endif                                                                     
           enddo                                                                     
          endif                                                                      
          enddo                                                                      
         endif                                                                       
         enddo                                                                       
        enddo                                                                        
        j0 = 5 ; k0 = 6                                                              
        do i0 = 0, 4                                                                 
         do i1 = 0, 4                                                                
         if (i1.ne.i0) then                                                          
          ys(0) = ys(0)+d(i0,i1)*d(i0,j0)/20                                         
          ys(1) = ys(1)+d(i0,i1)*d(i0,k0)/20                                         
          zs(2) = zs(2)+d2(i0,i1)*d(i0,j0)/20                                        
          zs(3) = zs(3)+d(i0,i1)*d2(i0,j0)/20                                        
          zs(4) = zs(4)+d(i0,i1)*d(i0,j0)*d(i1,j0)/20                                
          zs(5) = zs(5)+d2(i0,i1)*d(i0,k0)/20                                        
          zs(6) = zs(6)+d(i0,i1)*d2(i0,k0)/20                                        
          zs(7) = zs(7)+d(i0,i1)*d(i0,k0)*d(i1,k0)/20                                
          us(11) = us(11)+d3(i0,i1)*d(i0,j0)/20                                      
          us(12) = us(12)+d2(i0,i1)*d2(i0,j0)/20                                     
          us(13) = us(13)+d(i0,i1)*d3(i0,j0)/20                                      
          us(14) = us(14)+d2(i0,i1)*d(i0,j0)*d(i1,j0)/20                             
          us(15) = us(15)+d(i0,i1)*d2(i0,j0)*d(i1,j0)/20                             
          us(16) = us(16)+d3(i0,i1)*d(i0,k0)/20                                      
          us(17) = us(17)+d2(i0,i1)*d2(i0,k0)/20                                     
          us(18) = us(18)+d(i0,i1)*d3(i0,k0)/20                                      
          us(19) = us(19)+d2(i0,i1)*d(i0,k0)*d(i1,k0)/20                             
          us(20) = us(20)+d(i0,i1)*d2(i0,k0)*d(i1,k0)/20                             
          vs(56) = vs(56)+d4(i0,i1)*d(i0,j0)/20                                      
          vs(57) = vs(57)+d3(i0,i1)*d2(i0,j0)/20                                     
          vs(58) = vs(58)+d2(i0,i1)*d3(i0,j0)/20                                     
          vs(59) = vs(59)+d(i0,i1)*d4(i0,j0)/20                                      
          vs(60) = vs(60)+d3(i0,i1)*d(i0,j0)*d(i1,j0)/20                             
          vs(61) = vs(61)+d2(i0,i1)*d2(i0,j0)*d(i1,j0)/20                            
          vs(62) = vs(62)+d(i0,i1)*d3(i0,j0)*d(i1,j0)/20                             
          vs(63) = vs(63)+d4(i0,i1)*d(i0,k0)/20                                      
          vs(64) = vs(64)+d3(i0,i1)*d2(i0,k0)/20                                     
          vs(65) = vs(65)+d2(i0,i1)*d3(i0,k0)/20                                     
          vs(66) = vs(66)+d(i0,i1)*d4(i0,k0)/20                                      
          vs(67) = vs(67)+d3(i0,i1)*d(i0,k0)*d(i1,k0)/20                             
          vs(68) = vs(68)+d2(i0,i1)*d2(i0,k0)*d(i1,k0)/20                            
          vs(69) = vs(69)+d(i0,i1)*d3(i0,k0)*d(i1,k0)/20                             
          do i2 = 0, 4                                                               
          if (i2.ne.i1.and.i2.ne.i0) then                                            
           zs(8) = zs(8)+d(i0,i1)*d(i0,i2)*d(i0,j0)/60                               
           zs(9) = zs(9)+d(i0,i1)*d(i1,i2)*d(i0,j0)/60                               
           zs(10) = zs(10)+d(i0,i1)*d(i0,i2)*d(i0,k0)/60                             
           zs(11) = zs(11)+d(i0,i1)*d(i1,i2)*d(i0,k0)/60                             
           us(21) = us(21)+d2(i0,i1)*d(i0,i2)*d(i0,j0)/60                            
           us(22) = us(22)+d2(i0,i1)*d(i1,i2)*d(i0,j0)/60                            
           us(23) = us(23)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)/60                    
           us(24) = us(24)+d(i0,i1)*d2(i1,i2)*d(i0,j0)/60                            
           us(25) = us(25)+d(i0,i1)*d(i0,i2)*d2(i0,j0)/60                            
           us(26) = us(26)+d(i0,i1)*d(i1,i2)*d2(i0,j0)/60                            
           us(27) = us(27)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)/60                    
           us(28) = us(28)+d2(i0,i1)*d(i0,i2)*d(i0,k0)/60                            
           us(29) = us(29)+d2(i0,i1)*d(i1,i2)*d(i0,k0)/60                            
           us(30) = us(30)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,k0)/60                    
           us(31) = us(31)+d(i0,i1)*d2(i1,i2)*d(i0,k0)/60                            
           us(32) = us(32)+d(i0,i1)*d(i0,i2)*d2(i0,k0)/60                            
           us(33) = us(33)+d(i0,i1)*d(i1,i2)*d2(i0,k0)/60                            
           us(34) = us(34)+d(i0,i1)*d(i0,i2)*d(i0,k0)*d(i1,k0)/60                    
           vs(70) = vs(70)+d3(i0,i1)*d(i0,i2)*d(i0,j0)/60                            
           vs(71) = vs(71)+d2(i0,i1)*d2(i0,i2)*d(i0,j0)/60                           
           vs(72) = vs(72)+d3(i0,i1)*d(i1,i2)*d(i0,j0)/60                            
           vs(73) = vs(73)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)/60                   
           vs(74) = vs(74)+d2(i0,i1)*d2(i1,i2)*d(i0,j0)/60                           
           vs(75) = vs(75)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,j0)/60                   
           vs(76) = vs(76)+d(i0,i1)*d3(i1,i2)*d(i0,j0)/60                            
           vs(77) = vs(77)+d2(i0,i1)*d(i0,i2)*d2(i0,j0)/60                           
           vs(78) = vs(78)+d2(i0,i1)*d(i1,i2)*d2(i0,j0)/60                           
           vs(79) = vs(79)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,j0)/60                   
           vs(80) = vs(80)+d(i0,i1)*d2(i1,i2)*d2(i0,j0)/60                           
           vs(81) = vs(81)+d(i0,i1)*d(i0,i2)*d3(i0,j0)/60                            
           vs(82) = vs(82)+d(i0,i1)*d(i1,i2)*d3(i0,j0)/60                            
           vs(83) = vs(83)+d2(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)/60                   
           vs(84) = vs(84)+d(i0,i1)*d2(i0,i2)*d(i0,j0)*d(i1,j0)/60                   
         vs(85) = vs(85)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,j0)/60           
           vs(86) = vs(86)+d(i0,i1)*d(i0,i2)*d2(i0,j0)*d(i1,j0)/60                   
           vs(87) = vs(87)+d3(i0,i1)*d(i0,i2)*d(i0,k0)/60                            
           vs(88) = vs(88)+d2(i0,i1)*d2(i0,i2)*d(i0,k0)/60                           
           vs(89) = vs(89)+d3(i0,i1)*d(i1,i2)*d(i0,k0)/60                            
           vs(90) = vs(90)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,k0)/60                   
           vs(91) = vs(91)+d2(i0,i1)*d2(i1,i2)*d(i0,k0)/60                           
           vs(92) = vs(92)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,k0)/60                   
           vs(93) = vs(93)+d(i0,i1)*d3(i1,i2)*d(i0,k0)/60                            
           vs(94) = vs(94)+d2(i0,i1)*d(i0,i2)*d2(i0,k0)/60                           
           vs(95) = vs(95)+d2(i0,i1)*d(i1,i2)*d2(i0,k0)/60                           
           vs(96) = vs(96)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,k0)/60                   
           vs(97) = vs(97)+d(i0,i1)*d2(i1,i2)*d2(i0,k0)/60                           
           vs(98) = vs(98)+d(i0,i1)*d(i0,i2)*d3(i0,k0)/60                            
           vs(99) = vs(99)+d(i0,i1)*d(i1,i2)*d3(i0,k0)/60                            
           vs(100) = vs(100)+d2(i0,i1)*d(i0,i2)*d(i0,k0)*d(i1,k0)/60                 
           vs(101) = vs(101)+d(i0,i1)*d2(i0,i2)*d(i0,k0)*d(i1,k0)/60                 
           vs(102) = vs(102)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,k0)*d(i1,
     $k0)/60         
           vs(103) = vs(103)+d(i0,i1)*d(i0,i2)*d2(i0,k0)*d(i1,k0)/60                 
           do i3 = 0, 4                                                              
           if (i3.ne.i0.and.i3.ne.i1.and.i3.ne.i2) then                              
            us(35) = us(35)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)/120                  
            us(36) = us(36)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)/120                  
            us(37) = us(37)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,k0)/120                  
            us(38) = us(38)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,k0)/120                  
            vs(104) = vs(104)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)/120               
            vs(105) = vs(105)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)/120               
            vs(106) = vs(106)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,
     $j0)/120       
            vs(107) = vs(107)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,j0)/120               
            vs(108) = vs(108)+d(i0,i1)*d(i1,i2)*d2(i0,i3)*d(i0,j0)/120               
            vs(109) = vs(109)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d2(i0,j0)/120               
            vs(110) = vs(110)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d2(i0,j0)/120               
            vs(111) = vs(111)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,k0)/120               
            vs(112) = vs(112)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,k0)/120               
            vs(113) = vs(113)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,
     $k0)/120       
            vs(114) = vs(114)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,k0)/120               
            vs(115) = vs(115)+d(i0,i1)*d(i1,i2)*d2(i0,i3)*d(i0,k0)/120               
            vs(116) = vs(116)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d2(i0,k0)/120               
            vs(117) = vs(117)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d2(i0,k0)/120               
           endif                                                                     
           enddo                                                                     
          endif                                                                      
          enddo                                                                      
         endif                                                                       
         enddo                                                                       
        enddo                                                                        
        do i0 = 0, 4                                                                 
         ys(2) = ys(2)+d(i0,j0)*d(i0,k0)/5                                           
         zs(12) = zs(12)+d2(i0,j0)*d(i0,k0)/5                                        
         zs(13) = zs(13)+d(i0,j0)*d2(i0,k0)/5                                        
         do i1 = 0, 4                                                                
         if (i1.ne.i0) then                                                          
          zs(14) = zs(14)+d(i0,i1)*d(i0,j0)*d(i0,k0)/20                              
          zs(15) = zs(15)+d(i0,i1)*d(i1,j0)*d(i0,k0)/20                              
          us(39) = us(39)+d2(i0,i1)*d(i0,j0)*d(i0,k0)/20                             
          us(40) = us(40)+d(i0,i1)*d2(i0,j0)*d(i0,k0)/20                             
          us(41) = us(41)+d3(i0,j0)*d(i0,k0)/20                                      
          us(42) = us(42)+d2(i0,i1)*d(i1,j0)*d(i0,k0)/20                             
          us(43) = us(43)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,k0)/20                     
          us(44) = us(44)+d(i0,i1)*d2(i1,j0)*d(i0,k0)/20                             
          us(45) = us(45)+d(i0,i1)*d(i0,j0)*d2(i0,k0)/20                             
          us(46) = us(46)+d2(i0,j0)*d2(i0,k0)/20                                     
          us(47) = us(47)+d(i0,i1)*d(i1,j0)*d2(i0,k0)/20                             
          us(48) = us(48)+d(i0,j0)*d3(i0,k0)/20                                      
          us(49) = us(49)+d(i0,i1)*d(i0,j0)*d(i0,k0)*d(i1,k0)/20                     
          vs(118) = vs(118)+d3(i0,i1)*d(i0,j0)*d(i0,k0)/20                           
          vs(119) = vs(119)+d2(i0,i1)*d2(i0,j0)*d(i0,k0)/20                          
          vs(120) = vs(120)+d(i0,i1)*d3(i0,j0)*d(i0,k0)/20                           
          vs(121) = vs(121)+d4(i0,j0)*d(i0,k0)/20                                    
          vs(122) = vs(122)+d3(i0,i1)*d(i1,j0)*d(i0,k0)/20                           
          vs(123) = vs(123)+d2(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,k0)/20                  
          vs(124) = vs(124)+d(i0,i1)*d2(i0,j0)*d(i1,j0)*d(i0,k0)/20                  
          vs(125) = vs(125)+d2(i0,i1)*d2(i1,j0)*d(i0,k0)/20                          
          vs(126) = vs(126)+d(i0,i1)*d(i0,j0)*d2(i1,j0)*d(i0,k0)/20                  
          vs(127) = vs(127)+d2(i0,i1)*d(i0,j0)*d2(i0,k0)/20                          
          vs(128) = vs(128)+d(i0,i1)*d2(i0,j0)*d2(i0,k0)/20                          
          vs(129) = vs(129)+d3(i0,j0)*d2(i0,k0)/20                                   
          vs(130) = vs(130)+d2(i0,i1)*d(i1,j0)*d2(i0,k0)/20                          
          vs(131) = vs(131)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d2(i0,k0)/20                  
          vs(132) = vs(132)+d(i0,i1)*d2(i1,j0)*d2(i0,k0)/20                          
          vs(133) = vs(133)+d(i0,i1)*d(i0,j0)*d3(i0,k0)/20                           
          vs(134) = vs(134)+d2(i0,j0)*d3(i0,k0)/20                                   
          vs(135) = vs(135)+d(i0,i1)*d(i1,j0)*d3(i0,k0)/20                           
          vs(136) = vs(136)+d(i0,j0)*d4(i0,k0)/20                                    
          vs(137) = vs(137)+d2(i0,i1)*d(i0,j0)*d(i0,k0)*d(i1,k0)/20                  
          vs(138) = vs(138)+d(i0,i1)*d2(i0,j0)*d(i0,k0)*d(i1,k0)/20                  
          vs(139) = vs(139)+d(i0,i1)*d(i0,j0)*d2(i0,k0)*d(i1,k0)/20                  
          do i2 = 0, 4                                                               
          if (i2.ne.i1.and.i2.ne.i0) then                                            
           us(50) = us(50)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,k0)/60                    
           us(51) = us(51)+d(i0,i1)*d(i1,i2)*d(i0,j0)*d(i0,k0)/60                    
           us(52) = us(52)+d(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,k0)/60                    
           us(53) = us(53)+d(i0,i1)*d(i1,i2)*d(i1,j0)*d(i0,k0)/60                    
           vs(140) = vs(140)+d2(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,k0)/60                 
           vs(141) = vs(141)+d2(i0,i1)*d(i1,i2)*d(i0,j0)*d(i0,k0)/60                 
           vs(142) = vs(142)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i0,
     $k0)/60         
           vs(143) = vs(143)+d(i0,i1)*d2(i1,i2)*d(i0,j0)*d(i0,k0)/60                 
           vs(144) = vs(144)+d(i0,i1)*d(i0,i2)*d2(i0,j0)*d(i0,k0)/60                 
           vs(145) = vs(145)+d(i0,i1)*d(i1,i2)*d2(i0,j0)*d(i0,k0)/60                 
           vs(146) = vs(146)+d2(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,k0)/60                 
           vs(147) = vs(147)+d(i0,i1)*d2(i0,i2)*d(i1,j0)*d(i0,k0)/60                 
           vs(148) = vs(148)+d2(i0,i1)*d(i1,i2)*d(i1,j0)*d(i0,k0)/60                 
           vs(149) = vs(149)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,j0)*d(i0,
     $k0)/60         
           vs(150) = vs(150)+d2(i0,i2)*d(i1,i2)*d(i1,j0)*d(i0,k0)/60                 
           vs(151) = vs(151)+d(i0,i1)*d2(i1,i2)*d(i1,j0)*d(i0,k0)/60                 
           vs(152) = vs(152)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)*d(i0,
     $k0)/60         
           vs(153) = vs(153)+d(i0,i1)*d(i1,i2)*d(i0,j0)*d(i1,j0)*d(i0,
     $k0)/60         
           vs(154) = vs(154)+d(i0,i1)*d(i0,i2)*d2(i1,j0)*d(i0,k0)/60                 
           vs(155) = vs(155)+d(i0,i1)*d(i0,i2)*d(i1,j0)*d(i2,j0)*d(i0,
     $k0)/60         
           vs(156) = vs(156)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d2(i0,k0)/60                 
           vs(157) = vs(157)+d(i0,i1)*d(i1,i2)*d(i0,j0)*d2(i0,k0)/60                 
           vs(158) = vs(158)+d(i0,i1)*d(i0,i2)*d(i1,j0)*d2(i0,k0)/60                 
           vs(159) = vs(159)+d(i0,i1)*d(i1,i2)*d(i1,j0)*d2(i0,k0)/60                 
           vs(160) = vs(160)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,k0)*d(i1,
     $k0)/60         
           vs(161) = vs(161)+d(i0,i1)*d(i0,i2)*d(i2,j0)*d(i0,k0)*d(i1,
     $k0)/60         
           do i3 = 0, 4                                                              
           if (i3.ne.i0.and.i3.ne.i1.and.i3.ne.i2) then                              
            vs(162) = vs(162)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,j0)*d(i0,
     $k0)/120       
            vs(163) = vs(163)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,j0)*d(i0,
     $k0)/120       
            vs(164) = vs(164)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i1,j0)*d(i0,
     $k0)/120       
           endif                                                                     
           enddo                                                                     
          endif                                                                      
          enddo                                                                      
         endif                                                                       
         enddo                                                                       
        enddo                                                                        
! Reducible secondaries                                                      
! us(0:5)                                                                    
        k = 0                                                                        
        do i = 0, 2                                                                  
         us(k) = her2(ys(i)) ; k = k+1                                               
         do j = i+1, 2                                                               
          us(k) = ys(i)*ys(j) ; k = k+1                                              
         enddo                                                                       
        enddo                                                                        
! vs(0:47)                                                                   
        do j = 0, 15                                                                 
         vs(3*j:3*j+2) = ys(0:2)*zs(j)                                               
        enddo                                                                        
! ws(..).  We don't have this yet.                                           
! Compute vec(0:*).  The code was created using these parameters:            
! MolienSeries(0:*): 1 4 17 68 265 977 2942 4 4                              
! #Primaries(1:*):   4 4 4 4 4 1 0 3                                         
! #Secondaries(1:*): 0 3 16 54 165 0 1 5                                     
! The MolienSeries partial sums (allowed size of vec) are:                   
! Molien Sums(0:*):  1 5 22 90 355 1332 4274 4278 4282                       
! constant term                                                              
        vec(0) = 1                                                                   
! first degree terms                                                         
        if (1.le.mdeg) then                                                          
         vec(1) = x(0)                                                               
         vec(2) = x(1)                                                               
         vec(3) = x(2)                                                               
         vec(4) = x(3)                                                               
        endif                                                                        
! second degree terms                                                        
        if (2.le.mdeg) then                                                          
         vec(5) = x2(0)                                                              
         vec(6:8) = x(0)*vec(2:4)                                                    
         vec(9) = x2(1)                                                              
         vec(10:11) = x(1)*vec(3:4)                                                  
         vec(12) = x2(2)                                                             
         vec(13:13) = x(2)*vec(4:4)                                                  
         vec(14) = x2(3)                                                             
         vec(15) = y(0)                                                              
         vec(16) = y(1)                                                              
         vec(17) = y(2)                                                              
         vec(18) = y(3)                                                              
         vec(19:21) = ys(0:2)                                                        
        endif                                                                        
! third degree terms                                                         
        if (3.le.mdeg) then                                                          
         vec(22) = x3(0)                                                             
         vec(23:25) = x2(0)*vec(2:4)                                                 
         vec(26:38) = x(0)*vec(9:21)                                                 
         vec(39) = x3(1)                                                             
         vec(40:41) = x2(1)*vec(3:4)                                                 
         vec(42:51) = x(1)*vec(12:21)                                                
         vec(52) = x3(2)                                                             
         vec(53:53) = x2(2)*vec(4:4)                                                 
         vec(54:61) = x(2)*vec(14:21)                                                
         vec(62) = x3(3)                                                             
         vec(63:69) = x(3)*vec(15:21)                                                
         vec(70) = z(0)                                                              
         vec(71) = z(1)                                                              
         vec(72) = z(2)                                                              
         vec(73) = z(3)                                                              
         vec(74:89) = zs(0:15)                                                       
        endif                                                                        
! fourth degree terms                                                        
        if (4.le.mdeg) then                                                          
         vec(90) = x4(0)                                                             
         vec(91:93) = x3(0)*vec(2:4)                                                 
         vec(94:106) = x2(0)*vec(9:21)                                               
         vec(107:157) = x(0)*vec(39:89)                                              
         vec(158) = x4(1)                                                            
         vec(159:160) = x3(1)*vec(3:4)                                               
         vec(161:170) = x2(1)*vec(12:21)                                             
         vec(171:208) = x(1)*vec(52:89)                                              
         vec(209) = x4(2)                                                            
         vec(210:210) = x3(2)*vec(4:4)                                               
         vec(211:218) = x2(2)*vec(14:21)                                             
         vec(219:246) = x(2)*vec(62:89)                                              
         vec(247) = x4(3)                                                            
         vec(248:254) = x2(3)*vec(15:21)                                             
         vec(255:274) = x(3)*vec(70:89)                                              
         vec(275) = y2(0)                                                            
         vec(276:281) = y(0)*vec(16:21)                                              
         vec(282) = y2(1)                                                            
         vec(283:287) = y(1)*vec(17:21)                                              
         vec(288) = y2(2)                                                            
         vec(289:292) = y(2)*vec(18:21)                                              
         vec(293) = y2(3)                                                            
         vec(294:296) = y(3)*vec(19:21)                                              
         vec(297) = u(0)                                                             
         vec(298) = u(1)                                                             
         vec(299) = u(2)                                                             
         vec(300) = u(3)                                                             
         vec(301:354) = us(0:53)                                                     
        endif                                                                        
! fifth degree terms                                                         
        if (5.le.mdeg) then                                                          
         vec(355) = x5(0)                                                            
         vec(356:358) = x4(0)*vec(2:4)                                               
         vec(359:371) = x3(0)*vec(9:21)                                              
         vec(372:422) = x2(0)*vec(39:89)                                             
         vec(423:619) = x(0)*vec(158:354)                                            
         vec(620) = x5(1)                                                            
         vec(621:622) = x4(1)*vec(3:4)                                               
         vec(623:632) = x3(1)*vec(12:21)                                             
         vec(633:670) = x2(1)*vec(52:89)                                             
         vec(671:816) = x(1)*vec(209:354)                                            
         vec(817) = x5(2)                                                            
         vec(818:818) = x4(2)*vec(4:4)                                               
         vec(819:826) = x3(2)*vec(14:21)                                             
         vec(827:854) = x2(2)*vec(62:89)                                             
         vec(855:962) = x(2)*vec(247:354)                                            
         vec(963) = x5(3)                                                            
         vec(964:970) = x3(3)*vec(15:21)                                             
         vec(971:990) = x2(3)*vec(70:89)                                             
         vec(991:1070) = x(3)*vec(275:354)                                           
         vec(1071:1090) = y(0)*vec(70:89)                                            
         vec(1091:1110) = y(1)*vec(70:89)                                            
         vec(1111:1130) = y(2)*vec(70:89)                                            
         vec(1131:1150) = y(3)*vec(70:89)                                            
         vec(1151:1153) = z(0)*ys(0:2)                                               
         vec(1154:1156) = z(1)*ys(0:2)                                               
         vec(1157:1159) = z(2)*ys(0:2)                                               
         vec(1160:1162) = z(3)*ys(0:2)                                               
         vec(1163) = v(0)                                                            
         vec(1164) = v(1)                                                            
         vec(1165) = v(2)                                                            
         vec(1166) = v(3)                                                            
         vec(1167:1331) = vs(0:164)                                                  
        endif                                                                        
! sixth degree terms                                                         
        if (6.le.mdeg) then                                                          
         vec(1332) = x6(0)                                                           
         vec(1333:1335) = x5(0)*vec(2:4)                                             
         vec(1336:1348) = x4(0)*vec(9:21)                                            
         vec(1349:1399) = x3(0)*vec(39:89)                                           
         vec(1400:1596) = x2(0)*vec(158:354)                                         
         vec(1597:2308) = x(0)*vec(620:1331)                                         
         vec(2309) = x6(1)                                                           
         vec(2310:2311) = x5(1)*vec(3:4)                                             
         vec(2312:2321) = x4(1)*vec(12:21)                                           
         vec(2322:2359) = x3(1)*vec(52:89)                                           
         vec(2360:2505) = x2(1)*vec(209:354)                                         
         vec(2506:3020) = x(1)*vec(817:1331)                                         
         vec(3021) = x6(2)                                                           
         vec(3022:3022) = x5(2)*vec(4:4)                                             
         vec(3023:3030) = x4(2)*vec(14:21)                                           
         vec(3031:3058) = x3(2)*vec(62:89)                                           
         vec(3059:3166) = x2(2)*vec(247:354)                                         
         vec(3167:3535) = x(2)*vec(963:1331)                                         
         vec(3536) = x6(3)                                                           
         vec(3537:3543) = x4(3)*vec(15:21)                                           
         vec(3544:3563) = x3(3)*vec(70:89)                                           
         vec(3564:3643) = x2(3)*vec(275:354)                                         
         vec(3644:3904) = x(3)*vec(1071:1331)                                        
         vec(3905) = y3(0)                                                           
         vec(3906:3911) = y2(0)*vec(16:21)                                           
         vec(3912:3984) = y(0)*vec(282:354)                                          
         vec(3985) = y3(1)                                                           
         vec(3986:3990) = y2(1)*vec(17:21)                                           
         vec(3991:4057) = y(1)*vec(288:354)                                          
         vec(4058) = y3(2)                                                           
         vec(4059:4062) = y2(2)*vec(18:21)                                           
         vec(4063:4124) = y(2)*vec(293:354)                                          
         vec(4125) = y3(3)                                                           
         vec(4126:4128) = y2(3)*vec(19:21)                                           
         vec(4129:4186) = y(3)*vec(297:354)                                          
         vec(4187) = z2(0)                                                           
         vec(4188:4206) = z(0)*vec(71:89)                                            
         vec(4207) = z2(1)                                                           
         vec(4208:4225) = z(1)*vec(72:89)                                            
         vec(4226) = z2(2)                                                           
         vec(4227:4243) = z(2)*vec(73:89)                                            
         vec(4244) = z2(3)                                                           
         vec(4245:4260) = z(3)*vec(74:89)                                            
         vec(4261:4263) = u(0)*ys(0:2)                                               
         vec(4264:4266) = u(1)*ys(0:2)                                               
         vec(4267:4269) = u(2)*ys(0:2)                                               
         vec(4270:4272) = u(3)*ys(0:2)                                               
         vec(4273) = w(0)                                                            
!! ws(0:*) to follow                                                         
        endif                                                                        
        return                                                                       
        end                                                                         

!********************************************************


        function her2(t0)
        implicit none
        
        double precision   t0,her2
        her2 = (4*t0**2-2)/sqrt(dble(8*2 ))
        return
        end
        
        function her3(t0)
        implicit none
        
        double precision   t0,her3
        her3 = (8*t0**2-12)*t0/sqrt(dble(16*6 ))
        return
        end
        
        function her4(t0)
        implicit none
        
        double precision   t0,her4
        her4 = ((16*t0**2-48)*t0**2+12)/sqrt(dble(32*24 ))
        return
        end
        
        function her5(t0)
        implicit none
        
        double precision   t0,her5
        her5 = ((32*t0**2-160)*t0**2+120)*t0/sqrt(dble(64*120 ))
        return
        end
        
        function her6(t0)
        implicit none
        
        double precision   t0,her6
        her6=(((64*t0**2-480)*t0**2+720)*t0**2-120)/sqrt(dble(128*720))
        return
        end
        
        function her7(t0)
        implicit none
        
        double precision   t0,her7
        her7 = (((128*t0**2-1344)*t0**2+3360)*t0**2-1680)*t0/sqrt(dble(
     $ 256*5040 ))
        return
        end

!*************************************************
 
