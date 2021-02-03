      SUBROUTINE CH4HPOT_INTER(X1, en1)
C, initialized)
      implicit double precision (A-H,O-Z)
C
      double precision X1(6,3)
      double precision en1
C      logical initialized

      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER

C      IF(initialized.NEQV..TRUE.) THEN
C        call PREPOT
C        initialized = .TRUE.
C      ENDIF

C MR: line to change dimensions for i-PI of X1
      do i = 1,6
         do j = 1,3
            CART(i,j)=X1(i,j)
         enddo
      enddo

      call POT

      en1=PENGYGS

      END SUBROUTINE


C**************************************************************************
C   System:           CH5
C   Number of electronic surfaces: 1
C   Number of derivatives: 0
C   Number of bodies: 6
C   Interface:        potlib2001
C   Common name:      CH5
C   Notes:
C
C   Reference: J.C. Corchado, J.L. Bravo and J. Espinosa-Garcia 
C                  J.Chem.Phys.130,184314 (2009).
C   Fitting to ab initio calculations at the CCSD(T)/cc-pVTZ level.
C
      SUBROUTINE POT
C
C        This potential is written such that:
C
C                       X(1)  - X(3)  : X, Y, Z for H1
C                       X(4)  - X(6)  : X, Y, Z for C
C                       X(7)  - X(9)  : X, Y, Z for H3
C                       X(10) - X(12) : X, Y, Z for H4
C                       X(13) - X(15) : X, Y, Z for H2
C                       X(16) - X(18) : X, Y, Z for H5
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
      CALL CARTOU
      CALL CARTTOR
C
C Changing units from angstrom to bohr and initialize
C
      DO 10 I = 1, 18
         q(I) = R(I)*0.52918d0
         pdot(I)=0.D0
 10   CONTINUE
c
c  calculate relative coordinates and bond lengths
c
       en=0.0d0

       call coorden
c
c  calculate switching functions
c
       call switchf
c
c  calculate reference angles and their derivatives
c
       call refangles
c
c  calculate stretching potential
c
       call stretch(vstr)
c
c  calculate out of plane bending potential
c
       call opbend(vop)
c
c  calculate in plane bending potential
c
       call ipbend(vip)
c
c  total potential energy is vstr+vop+vip
c
       en=vstr+vop+vip
c
c 0.03812 conversion factor from 10(5) j/mol to hartrees
c
       en = en*0.03812D0
       ENGYGS = en
c
      CALL EUNITZERO
      IF(NDER.NE.0) THEN
c
c 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
c
         do i=1,18
            DEGSDR(i)=pdot(i)*0.0201723d0
         enddo
         CALL RTOCART
         CALL DEDCOU
      ENDIF
C
       return
       end
c
c******************************************************
c
       subroutine coorden
c
c  calculates relative coordinates and bond lengths
c
       implicit double precision (a-h,o-z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
c  calculate relative coordinates
c
       do ind=1,3
         tcb(ind)=q(nc(ind))-q(nhb(ind))
         do i=1,4
           tch(i,ind)=q(nc(ind))-q(nh(i,ind))
           tbh(i,ind)=q(nhb(ind))-q(nh(i,ind))
         enddo
       enddo
c
c  calculate bond lengths
c
       rcb=sqrt(tcb(1)*tcb(1)+tcb(2)*tcb(2)+tcb(3)*tcb(3))
       do i=1,4
         rch(i)=sqrt(tch(i,1)*tch(i,1)+tch(i,2)*tch(i,2)+
     *                tch(i,3)*tch(i,3))
         rbh(i)=sqrt(tbh(i,1)*tbh(i,1)+tbh(i,2)*tbh(i,2)+
     *                tbh(i,3)*tbh(i,3))
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine refangles
c
c  subroutine calculates reference angles for the "in-plane" potential
c
       implicit double precision (a-h,o-z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       tau=acos(-1.0d0/3.0d0)
C       pi=4.0d0*atan(1.0d0)
       halfpi=0.5d0*pi
       twopi=2.0d0*pi
c
c  set diagonal elements to zero
c
       do i=1,4
         theta0(i,i)=0.0d0
         do k=1,4
           dtheta0(i,i,k)=0.0d0
         enddo
       enddo
c
c  calculate reference angles
c
       theta0(1,2)=tau+(tau-halfpi)*(sphi(1)*sphi(2)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
       theta0(1,3)=tau+(tau-halfpi)*(sphi(1)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
       theta0(1,4)=tau+(tau-halfpi)*(sphi(1)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
       theta0(2,3)=tau+(tau-halfpi)*(sphi(2)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
       theta0(2,4)=tau+(tau-halfpi)*(sphi(2)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
       theta0(3,4)=tau+(tau-halfpi)*(sphi(3)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
c
c  calculate the derivatives of theta0(i,j) in terms of rch(k)
c  quantity calulated is dtheta0(i,j,k)
c
c  derivatives wrt rch(1)
c
       dtheta0(1,2,1)=(tau-halfpi)*dsphi(1)*sphi(2)
       dtheta0(1,3,1)=(tau-halfpi)*dsphi(1)*sphi(3)
       dtheta0(1,4,1)=(tau-halfpi)*dsphi(1)*sphi(4)
       dtheta0(2,3,1)=(tau-twopi/3.0d0)*dstheta(1)*stheta(4)
       dtheta0(2,4,1)=(tau-twopi/3.0d0)*dstheta(1)*stheta(3)
       dtheta0(3,4,1)=(tau-twopi/3.0d0)*dstheta(1)*stheta(2)
c
c  derivatives wrt rch(2)
c
       dtheta0(1,2,2)=(tau-halfpi)*sphi(1)*dsphi(2)
       dtheta0(1,3,2)=(tau-twopi/3.0d0)*dstheta(2)*stheta(4)
       dtheta0(1,4,2)=(tau-twopi/3.0d0)*dstheta(2)*stheta(3)
       dtheta0(2,3,2)=(tau-halfpi)*dsphi(2)*sphi(3)
       dtheta0(2,4,2)=(tau-halfpi)*dsphi(2)*sphi(4)
       dtheta0(3,4,2)=(tau-twopi/3.0d0)*stheta(1)*dstheta(2)
c
c  derivatives wrt rch(3)
c
       dtheta0(1,2,3)=(tau-twopi/3.0d0)*dstheta(3)*stheta(4)
       dtheta0(1,3,3)=(tau-halfpi)*sphi(1)*dsphi(3)
       dtheta0(1,4,3)=(tau-twopi/3.0d0)*stheta(2)*dstheta(3)
       dtheta0(2,3,3)=(tau-halfpi)*sphi(2)*dsphi(3)
       dtheta0(2,4,3)=(tau-twopi/3.0d0)*stheta(1)*dstheta(3)
       dtheta0(3,4,3)=(tau-halfpi)*dsphi(3)*sphi(4)
c
c  derivatives wrt rch(4)
c
       dtheta0(1,2,4)=(tau-twopi/3.0d0)*stheta(3)*dstheta(4)
       dtheta0(1,3,4)=(tau-twopi/3.0d0)*stheta(2)*dstheta(4)
       dtheta0(1,4,4)=(tau-halfpi)*sphi(1)*dsphi(4)
       dtheta0(2,3,4)=(tau-twopi/3.0d0)*stheta(1)*dstheta(4)
       dtheta0(2,4,4)=(tau-halfpi)*sphi(2)*dsphi(4)
       dtheta0(3,4,4)=(tau-halfpi)*sphi(3)*dsphi(4)
c
c  fill in the other half of the matrix
c
        do i=1,3
          do j=i+1,4
            theta0(j,i)=theta0(i,j)
            do k=1,4
              dtheta0(j,i,k)=dtheta0(i,j,k)
            enddo
          enddo
        enddo
       return
       end
c
c******************************************************
c
c
       subroutine stretch(vstr)
c
c  subroutine to calculate leps-type stretching potential and its
c  derivatives
c
       implicit double precision (a-h,o-z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension vqch(4),vjch(4),vqbh(4),vjbh(4),vq(4),vj(4),
     *           achdc(3),achdh(4,3)
c
c  calculate avergage bond length for the methane moiety
c
       rav=(rch(1)+rch(2)+rch(3)+rch(4))/4.0d0
c
c  initialise:
c
       vstr=0.0d0
c
c  ach:
c
c  in double precision tanh(19.0d0)=1.0d0 and we put the if statement
c  in to avoid overflow/underflow errors
c
       arga=c1ch*(rav-r0ch)
       if(arga.lt.19.0d0)then
         ach=a1ch+b1ch*(tanh(arga)+1.0d0)*0.5d0
         dumach=b1ch*c1ch/(2.0d0*cosh(arga)**2)
       else
         ach=a1ch+b1ch
         dumach=0.0d0
       endif
c
c  calculate singlet: e1, triplet: e3 energies and vq and vj
c  terms for each bond
c
       e1=d1cb*(exp(-2.0d0*acb*(rcb-r0cb))-2.0d0*exp(-acb*(rcb-r0cb)))
       e3=d3cb*(exp(-2.0d0*acb*(rcb-r0cb))+2.0d0*exp(-acb*(rcb-r0cb)))
       vqcb=(e1+e3)*0.5d0
       vjcb=(e1-e3)*0.5d0
       do i=1,4
         e1=d1ch*(exp(-2.0d0*ach*(rch(i)-r0ch))
     *              -2.0d0*exp(-ach*(rch(i)-r0ch)))
         e3=d3ch*(exp(-2.0d0*ach*(rch(i)-r0ch))
     *              +2.0d0*exp(-ach*(rch(i)-r0ch)))
         vqch(i)=(e1+e3)*0.5d0
         vjch(i)=(e1-e3)*0.5d0
         e1=d1hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *              -2.0d0*exp(-ahh*(rbh(i)-r0hh)))
         e3=d3hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *              +2.0d0*exp(-ahh*(rbh(i)-r0hh)))
         vqbh(i)=(e1+e3)*0.5d0
         vjbh(i)=(e1-e3)*0.5d0
c
c  calculate 3 body potential
c
         vq(i)=vqch(i)+vqcb+vqbh(i)
         vj(i)=-sqrt(((vjch(i)-vjcb)**2+(vjcb-vjbh(i))**2
     *                 +(vjbh(i)-vjch(i))**2)*0.5d0)
         vstr=vstr+vq(i)+vj(i)
       enddo
c
c  partial derivatives
c  first we need the derivative of ach:
c
       do ind=1,3
         achdc(ind)=dumach*(tch(1,ind)/rch(1)+tch(2,ind)/rch(2)
     *            +tch(3,ind)/rch(3)+tch(4,ind)/rch(4))/4.0d0
         do i=1,4
           achdh(i,ind)=-dumach*tch(i,ind)/rch(i)/4.0d0
         enddo
       enddo
       dumqcb=-acb*((d1cb+d3cb)*exp(-2.0d0*acb*(rcb-r0cb))-
     *         (d1cb-d3cb)*exp(-acb*(rcb-r0cb)))/rcb
c
c  calculate cartesian derivatives:
c  looping over ch(i) and bh(i)
c
       do i=1,4
         dumqbh=-ahh*((d1hh+d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh))-
     *           (d1hh-d3hh)*exp(-ahh*(rbh(i)-r0hh)))/rbh(i)
         factj=0.5d0/vj(i)
         dumjcb=-acb*((d1cb-d3cb)*exp(-2.0d0*acb*(rcb-r0cb))
     *            -(d1cb+d3cb)*exp(-acb*(rcb-r0cb)))*factj/rcb
         dumjbh=-ahh*((d1hh-d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *            -(d1hh+d3hh)*exp(-ahh*(rbh(i)-r0hh)))*factj/rbh(i)
         do ind=1,3
c
c  deriv wrt hb:
c
                  pdot(nhb(ind))=pdot(nhb(ind))
     *             -tcb(ind)*dumqcb+tbh(i,ind)*dumqbh
     *            +(vjch(i)-vjcb)*(dumjcb*tcb(ind))
     *            +(vjcb-vjbh(i))*(-dumjcb*tcb(ind)-dumjbh*tbh(i,ind))
     *            +(vjbh(i)-vjch(i))*dumjbh*tbh(i,ind)
c
c  dvqch(i)/dc
c
           dumqch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch))
     *              *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                 -(d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
               pdot(nc(ind))=pdot(nc(ind))+dumqch+tcb(ind)*dumqcb
c
c  dvqch(i)/dh(i)
c
           dumqhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch))
     *              *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                 -(d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
              pdot(nh(i,ind))=pdot(nh(i,ind))+dumqhi-tbh(i,ind)*dumqbh
c
c  dvjch(i)/dc
c
           dumjch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch))
     *              *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *               -(d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
c
c  dvj(i)/dnc(ind)
c
           pdot(nc(ind))=pdot(nc(ind))
     *            +(vjch(i)-vjcb)*(dumjch-dumjcb*tcb(ind))
     *            +(vjcb-vjbh(i))*dumjcb*tcb(ind)
     *            -(vjbh(i)-vjch(i))*dumjch
c
c  dvjch(i)/dh(i)
c
           dumjhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch))
     *              *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *               -(d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
c
c  dvj(i)/dnh(i,ind)
c
            pdot(nh(i,ind))=pdot(nh(i,ind))
     *            +(vjch(i)-vjcb)*dumjhi
     *            +(vjcb-vjbh(i))*dumjbh*tbh(i,ind)
     *            +(vjbh(i)-vjch(i))*(-dumjbh*tbh(i,ind)-dumjhi)
c
c  dv(i)/dh(j)
c
           do k=1,3
             j=i+k
             if(j.gt.4)j=j-4
             dumqhj=-achdh(j,ind)*(rch(i)-r0ch)
     *                 *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                    -(d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
             dumjhj=-achdh(j,ind)*(rch(i)-r0ch)
     *                 *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                  -(d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
             pdot(nh(j,ind))=pdot(nh(j,ind))+dumqhj
     *            +(vjch(i)-vjcb)*dumjhj
     *            -(vjbh(i)-vjch(i))*dumjhj
           enddo
         enddo
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine opbend(vop)
c
c  subroutine calculates symmetrized vop potential and derivatives
c
       implicit double precision (a-h,o-z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       double precision norma
       dimension sumd2(4),sumd4(4)
       dimension in(3),a(3),b(3),axb(3),c(4,3),argd(4)
c
c
       vop=0.0d0
c
c  calculate force constants and their derivatives
c
       call opforce
c
c  calculate out-of-plane angle and derivatives
c
       do i=1,4
         j=i+1
         if(j.gt.4)j=j-4
         k=j+1
         if(k.gt.4)k=k-4
         l=k+1
         if(l.gt.4)l=l-4
c
c  modification to ensure that the set of methane CH bond vector
c  (rj,rk,rl) is a right-handed set
c
       in(1)=j
       in(2)=k
       in(3)=l
c
c  vector a is rk-rj, vector b is rl-rj
c
       do ind=1,3
         a(ind)=q(nh(k,ind))-q(nh(j,ind))
         b(ind)=q(nh(l,ind))-q(nh(j,ind))
       enddo
c
c  axb is vector a cross b
c
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
       do ii=1,3
         do ind=1,3
           c(in(ii),ind)=-tch(in(ii),ind)/rch(in(ii))
         enddo
       enddo
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma
c
c  if argd > 0 we need to switch vectors k and l around
c
         if (argd(in(ii)).gt.0.d0) then
             itemp=k
             k=l
             l=itemp
         endif
       enddo
c
c  subroutine performs sum over j, k, l
c  sum2 = sum delta**2
c  sum4 = sum delta**4
c
         call calcdelta(i,j,k,l,sum2,sum4)
         sumd2(i)=sum2
         sumd4(i)=sum4
         vop=vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
       enddo
       do i=1,4
         do j=1,4
c
c  overall derivatives of force constants i wrt the bond-length rch(j)
c
           ddr=dfdelta(i,j)*sumd2(i)+dhdelta(i,j)*sumd4(i)
c
c  calculate derivatives in terms of cartesian coordinates:
c
           do ind=1,3
             pdot(nh(j,ind))=pdot(nh(j,ind))-tch(j,ind)*ddr/rch(j)
             pdot(nc(ind))=pdot(nc(ind))+tch(j,ind)*ddr/rch(j)
           enddo
         enddo
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine ipbend(vip)
c
c  subroutine calculates symmetrised in plane bend term
c  and its derivatives
c
       implicit double precision (a-h,o-z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension costh(4,4),theta(4,4),dth(4,4)
c
c  initialise
c
       vip=0.0d0
c
c  calculate force constants: fk0(i,j), f1(i)
c  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
c
       call ipforce
c
c  calculate theta(i,j) and in plane bend potential
c
       do i=1,3
         do j=i+1,4
           costh(i,j)=tch(i,1)*tch(j,1)+tch(i,2)*tch(j,2)
     *                       +tch(i,3)*tch(j,3)
           costh(i,j)=costh(i,j)/rch(i)/rch(j)
           theta(i,j)=acos(costh(i,j))
           dth(i,j)=theta(i,j)-theta0(i,j)
           vip=vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
c
c  calculate partial derivatives wrt cartesian coordinates
c
c  calculate pdots wrt theta:
c
           termth=-1.0d0/sqrt(1.0d0-costh(i,j)*costh(i,j))
           do ind=1,3
             dthi=-tch(j,ind)/rch(i)/rch(j)
     *                  +costh(i,j)*tch(i,ind)/rch(i)/rch(i)
             dthi=dthi*termth
             dthj=-tch(i,ind)/rch(i)/rch(j)
     *                  +costh(i,j)*tch(j,ind)/rch(j)/rch(j)
             dthj=dthj*termth
             dthc=-(dthi+dthj)
             pdot(nh(i,ind))=pdot(nh(i,ind))
     *                       +fk0(i,j)*f1(i)*f1(j)*dthi*dth(i,j)
             pdot(nh(j,ind))=pdot(nh(j,ind))
     *                       +fk0(i,j)*f1(i)*f1(j)*dthj*dth(i,j)
             pdot(nc(ind))=pdot(nc(ind))
     *                       +fk0(i,j)*f1(i)*f1(j)*dthc*dth(i,j)
             do k=1,4
c
c  calculate pdots wrt force constants and wrt theta0
c
               dth0k=-dtheta0(i,j,k)*tch(k,ind)/rch(k)
               dth0c=-dth0k
               pdot(nh(k,ind))=pdot(nh(k,ind))
     *                -0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k)
     *                -0.5d0*tbh(k,ind)*dfdh(i,j,k)*dth(i,j)**2/rbh(k)
     *                      -fk0(i,j)*f1(i)*f1(j)*dth0k*dth(i,j)
               pdot(nc(ind))=pdot(nc(ind))
     *                +0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k)
     *                      -fk0(i,j)*f1(i)*f1(j)*dth0c*dth(i,j)
               pdot(nhb(ind))=pdot(nhb(ind))
     *                +0.5d0*tbh(k,ind)*dfdh(i,j,k)*dth(i,j)**2/rbh(k)
             enddo
           enddo
         enddo
       enddo
       return
       end
c
c*************************************************************************
c
       subroutine calcdelta(i,j,k,l,sum2,sum4)
c
c  subroutine calculates out of plane angle delta, loops
c  through delta(i,j), delta(i,k), delta(i,l)
c
c   also calculates the derivatives wrt delta
c
       implicit double precision (a-h,o-z)
       double precision norma
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension  delta(4),in(3),a(3),b(3),axb(3),c(4,3),argd(4),
     *            daxb(4,3,3),cdot(4,3,3),atemp2(3)
c
c  initialise
c
       sum2=0.0d0
       sum4=0.0d0
c
c  set j,k,l indices
c
       in(1)=j
       in(2)=k
       in(3)=l
c
c  vector a is rk-rj, vector b is rl-rj
c
       do ind=1,3
         a(ind)=q(nh(k,ind))-q(nh(j,ind))
         b(ind)=q(nh(l,ind))-q(nh(j,ind))
       enddo
c
c  axb is vector a cross b
c
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
       do ii=1,3
         do ind=1,3
           c(in(ii),ind)=-tch(in(ii),ind)/rch(in(ii))
         enddo
       enddo
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma
         delta(in(ii))=acos(argd(in(ii)))-theta0(i,in(ii))
c        write(*,*) 'theta,delta ',theta0(i,in(ii)),delta(in(ii))
         sum2=sum2+delta(in(ii))**2
         sum4=sum4+delta(in(ii))**4
       enddo
c
c  derivatives of axb wrt hj:
c
       daxb(j,1,1)=0.0d0
       daxb(j,1,2)=b(3)-a(3)
       daxb(j,1,3)=-b(2)+a(2)
       daxb(j,2,1)=-b(3)+a(3)
       daxb(j,2,2)=0.0d0
       daxb(j,2,3)=b(1)-a(1)
       daxb(j,3,1)=b(2)-a(2)
       daxb(j,3,2)=-b(1)+a(1)
       daxb(j,3,3)=0.0d0
c
c  derivatives of axb wrt hk:
c
       daxb(k,1,1)=0.0d0
       daxb(k,1,2)=-b(3)
       daxb(k,1,3)=b(2)
       daxb(k,2,1)=b(3)
       daxb(k,2,2)=0.0d0
       daxb(k,2,3)=-b(1)
       daxb(k,3,1)=-b(2)
       daxb(k,3,2)=b(1)
       daxb(k,3,3)=0.0d0
c
c  derivatives of axb wrt hl:
c
       daxb(l,1,1)=0.0d0
       daxb(l,1,2)=a(3)
       daxb(l,1,3)=-a(2)
       daxb(l,2,1)=-a(3)
       daxb(l,2,2)=0.0d0
       daxb(l,2,3)=a(1)
       daxb(l,3,1)=a(2)
       daxb(l,3,2)=-a(1)
       daxb(l,3,3)=0.0d0
c
c   loop over cdot(in(ii),ind,jind) where we consider deriv of c(in(ii))
c   wrt h(in(ii),jind) with components jind
c
       do ii=1,3
c
c  deriv of cdot(in(ii),x) wrt x, y, z
c
         cdot(in(ii),1,1)=1.0d0/rch(in(ii))
     *                   +tch(in(ii),1)*c(in(ii),1)/rch(in(ii))**2
         cdot(in(ii),1,2)=tch(in(ii),2)*c(in(ii),1)/rch(in(ii))**2
         cdot(in(ii),1,3)=tch(in(ii),3)*c(in(ii),1)/rch(in(ii))**2
c
c  deriv of cdot(in(ii),y) wrt x, y, z
c
         cdot(in(ii),2,1)=tch(in(ii),1)*c(in(ii),2)/rch(in(ii))**2
         cdot(in(ii),2,2)=1.0d0/rch(in(ii))
     *                   +tch(in(ii),2)*c(in(ii),2)/rch(in(ii))**2
         cdot(in(ii),2,3)=tch(in(ii),3)*c(in(ii),2)/rch(in(ii))**2
c
c  deriv of cdot(in(ii),z) wrt x, y, z
c
         cdot(in(ii),3,1)=tch(in(ii),1)*c(in(ii),3)/rch(in(ii))**2
         cdot(in(ii),3,2)=tch(in(ii),2)*c(in(ii),3)/rch(in(ii))**2
         cdot(in(ii),3,3)=1.0d0/rch(in(ii))
     *                   +tch(in(ii),3)*c(in(ii),3)/rch(in(ii))**2
       enddo
c
       do ii=1,3
         do ind=1,3
            deldot=-dtheta0(i,in(ii),i)
c
c  derivative wrt h(i,ind)
c  for  rch(i) only terms are from the derivatives of theta0
c
            deldot=-deldot*tch(i,ind)/rch(i)
            pdot(nh(i,ind))=pdot(nh(i,ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
c
c  derivative wrt c(ind)
c
            deldot=-deldot
            pdot(nc(ind))=pdot(nc(ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
           do jj=1,3
c
c  partial derivatives wrt h(in(jj),ind), loop over delta(i,in(ii))
c
c   atemp1 is axb dot daxb wrt h(in(jj))
c
            atemp1=axb(1)*daxb(in(jj),ind,1)
     *            +axb(2)*daxb(in(jj),ind,2)
     *            +axb(3)*daxb(in(jj),ind,3)
            atemp1=atemp1/(norma**3)
c
c  atemp2 is deriv of normalised axb
c
            atemp2(1)=daxb(in(jj),ind,1)/norma-atemp1*axb(1)
            atemp2(2)=daxb(in(jj),ind,2)/norma-atemp1*axb(2)
            atemp2(3)=daxb(in(jj),ind,3)/norma-atemp1*axb(3)
c
c  atemp3 is daxb dot c(in(ii))

            atemp3=atemp2(1)*c(in(ii),1)+atemp2(2)*c(in(ii),2)
     *                             +atemp2(3)*c(in(ii),3)
c
c  atemp4 is axb dot cdot
c
            atemp4=0.0d0
            if(ii.eq.jj)then
c
c  ie deriv of c(in(ii)) wrt h(in(jj)) is non zero only for ii = jj
c
              atemp4=axb(1)*cdot(in(ii),1,ind)
     *                     +axb(2)*cdot(in(ii),2,ind)
     *                     +axb(3)*cdot(in(ii),3,ind)
              atemp4=atemp4/norma
            endif
c
c  atemp5 is deriv of theta0(i,in(ii)) wrt to nh(in(jj),ind)
c
            atemp5=-dtheta0(i,in(ii),in(jj))
c
c  deriv wrt h(in(jj)),ind):
c
            atemp5=-atemp5*tch(in(jj),ind)/rch(in(jj))
            deldot=atemp3+atemp4
            deldot=-1.0d0/sqrt(1.0d0-argd(in(ii))**2)*deldot
            deldot=deldot+atemp5
            pdot(nh(in(jj),ind))=pdot(nh(in(jj),ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
c
c  for carbon the only contributions are from axb dot cdot term and
c  from theta0 and derivative cdot wrt carbon=-cdot wrt hydrogen
c
            deldot=1.0d0/sqrt(1.0d0-argd(in(ii))**2)*atemp4
            deldot=deldot-atemp5
            pdot(nc(ind))=pdot(nc(ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
          enddo
        enddo
       enddo
       return
       end

c******************************************************
c
       subroutine opforce
c
c  calculates the out-of-plane bending force constants
c  and their derivatives
c
       implicit double precision (a-h,o-z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension switch(4),dswitch(4,4)
c
c  calculate switching functions:
c
       switch(1)=(1.0d0-s3(1))*s3(2)*s3(3)*s3(4)
       switch(2)=(1.0d0-s3(2))*s3(3)*s3(4)*s3(1)
       switch(3)=(1.0d0-s3(3))*s3(4)*s3(1)*s3(2)
       switch(4)=(1.0d0-s3(4))*s3(1)*s3(2)*s3(3)
c
c  calculate derivatives:
c  derivative of switch(1) wrt the 4 rch bond lengths
c
       dswitch(1,1)=-ds3(1)*s3(2)*s3(3)*s3(4)
       dswitch(1,2)=(1.0d0-s3(1))*ds3(2)*s3(3)*s3(4)
       dswitch(1,3)=(1.0d0-s3(1))*s3(2)*ds3(3)*s3(4)
       dswitch(1,4)=(1.0d0-s3(1))*s3(2)*s3(3)*ds3(4)
c
c  derivative of switch(2) wrt the 4 rch bond lengths
c
       dswitch(2,1)=(1.0d0-s3(2))*s3(3)*s3(4)*ds3(1)
       dswitch(2,2)=-ds3(2)*s3(3)*s3(4)*s3(1)
       dswitch(2,3)=(1.0d0-s3(2))*ds3(3)*s3(4)*s3(1)
       dswitch(2,4)=(1.0d0-s3(2))*s3(3)*ds3(4)*s3(1)
c
c  derivative of switch(3) wrt the 4 rch bond lengths
c
       dswitch(3,1)=(1.0d0-s3(3))*s3(4)*ds3(1)*s3(2)
       dswitch(3,2)=(1.0d0-s3(3))*s3(4)*s3(1)*ds3(2)
       dswitch(3,3)=-ds3(3)*s3(4)*s3(1)*s3(2)
       dswitch(3,4)=(1.0d0-s3(3))*ds3(4)*s3(1)*s3(2)
c
c  derivative of switch(3) wrt the 4 rch bond lengths
c
       dswitch(4,1)=(1.0d0-s3(4))*ds3(1)*s3(2)*s3(3)
       dswitch(4,2)=(1.0d0-s3(4))*s3(1)*ds3(2)*s3(3)
       dswitch(4,3)=(1.0d0-s3(4))*s3(1)*s3(2)*ds3(3)
       dswitch(4,4)=-ds3(4)*s3(1)*s3(2)*s3(3)
c
c  calculate the force constants and their derivatives
c
       do i=1,4
         fdelta(i)=switch(i)*fch3
         hdelta(i)=switch(i)*hch3
         do j=1,4
           dfdelta(i,j)=dswitch(i,j)*fch3
           dhdelta(i,j)=dswitch(i,j)*hch3
         enddo
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine ipforce
c
c  calculates the symmetrised in plane bend force constants and
c  all partial derivatives involving them
c
       implicit double precision (a-h,o-z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension dfk0(4,4,4),df1dc(4),df1dh(4)
c
c  set force constant at asymptotes
c
       f0=fkinf+ak
       f2=fkinf
c
       fk0(1,2)=f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)*s2(4)-1.0d0)
       fk0(1,3)=f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)*s2(4)-1.0d0)
       fk0(1,4)=f0+f0*(s1(1)*s1(4)-1.0d0)+(f0-f2)*(s2(2)*s2(3)-1.0d0)
       fk0(2,3)=f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)*s2(4)-1.0d0)
       fk0(2,4)=f0+f0*(s1(2)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(3)-1.0d0)
       fk0(3,4)=f0+f0*(s1(3)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(2)-1.0d0)
c
c  derivative of fk0
c
       dfk0(1,2,1)=f0*ds1(1)*s1(2)
       dfk0(1,2,2)=f0*s1(1)*ds1(2)
       dfk0(1,2,3)=(f0-f2)*ds2(3)*s2(4)
       dfk0(1,2,4)=(f0-f2)*s2(3)*ds2(4)
c
       dfk0(1,3,1)=f0*ds1(1)*s1(3)
       dfk0(1,3,2)=(f0-f2)*ds2(2)*s2(4)
       dfk0(1,3,3)=f0*s1(1)*ds1(3)
       dfk0(1,3,4)=(f0-f2)*s2(2)*ds2(4)
c
       dfk0(1,4,1)=f0*ds1(1)*s1(4)
       dfk0(1,4,2)=(f0-f2)*ds2(2)*s2(3)
       dfk0(1,4,3)=(f0-f2)*s2(2)*ds2(3)
       dfk0(1,4,4)=f0*s1(1)*ds1(4)
c
       dfk0(2,3,1)=(f0-f2)*ds2(1)*s2(4)
       dfk0(2,3,2)=f0*ds1(2)*s1(3)
       dfk0(2,3,3)=f0*s1(2)*ds1(3)
       dfk0(2,3,4)=(f0-f2)*s2(1)*ds2(4)
c
       dfk0(2,4,1)=(f0-f2)*ds2(1)*s2(3)
       dfk0(2,4,2)=f0*ds1(2)*s1(4)
       dfk0(2,4,3)=(f0-f2)*s2(1)*ds2(3)
       dfk0(2,4,4)=f0*s1(2)*ds1(4)
c
       dfk0(3,4,1)=(f0-f2)*ds2(1)*s2(2)
       dfk0(3,4,2)=(f0-f2)*s2(1)*ds2(2)
       dfk0(3,4,3)=f0*ds1(3)*s1(4)
       dfk0(3,4,4)=f0*s1(3)*ds1(4)
c
       do i=1,4
c
c  calculate the terms f1(i)
c
         arga1=aa1*rbh(i)*rbh(i)
         arga2=aa4*(rbh(i)-r0hh)*(rbh(i)-r0hh)
         a1=1.0d0-exp(-arga1)
         a2=aa2+aa3*exp(-arga2)
         f1(i)=a1*exp(-a2*(rch(i)-r0ch)**2)
c
c  and calculate the derivatives wrt rch(i) and rbh(i)
c
         duma1=2.0d0*aa1*rbh(i)*exp(-arga1)
         duma2=-2.0d0*aa3*aa4*(rbh(i)-r0hh)*exp(-arga2)
         df1dc(i)=-2.0d0*(rch(i)-r0ch)*a1*a2*exp(-a2*(rch(i)-r0ch)**2)
         df1dh(i)=duma1*exp(-a2*(rch(i)-r0ch)**2)
     *             -duma2*(rch(i)-r0ch)**2*a1*exp(-a2*(rch(i)-r0ch)**2)
       enddo
c
c  derivative of total force constant f(i,j) wrt bond length rch(k)
c  is given by dfdc(i,j,k)
c
      dfdc(1,2,1)=dfk0(1,2,1)*f1(1)*f1(2)+fk0(1,2)*df1dc(1)*f1(2)
      dfdc(1,2,2)=dfk0(1,2,2)*f1(1)*f1(2)+fk0(1,2)*f1(1)*df1dc(2)
      dfdc(1,2,3)=dfk0(1,2,3)*f1(1)*f1(2)
      dfdc(1,2,4)=dfk0(1,2,4)*f1(1)*f1(2)
c
      dfdc(1,3,1)=dfk0(1,3,1)*f1(1)*f1(3)+fk0(1,3)*df1dc(1)*f1(3)
      dfdc(1,3,2)=dfk0(1,3,2)*f1(1)*f1(3)
      dfdc(1,3,3)=dfk0(1,3,3)*f1(1)*f1(3)+fk0(1,3)*f1(1)*df1dc(3)
      dfdc(1,3,4)=dfk0(1,3,4)*f1(1)*f1(3)
c
      dfdc(1,4,1)=dfk0(1,4,1)*f1(1)*f1(4)+fk0(1,4)*df1dc(1)*f1(4)
      dfdc(1,4,2)=dfk0(1,4,2)*f1(1)*f1(4)
      dfdc(1,4,3)=dfk0(1,4,3)*f1(1)*f1(4)
      dfdc(1,4,4)=dfk0(1,4,4)*f1(1)*f1(4)+fk0(1,4)*f1(1)*df1dc(4)
c
      dfdc(2,3,1)=dfk0(2,3,1)*f1(2)*f1(3)
      dfdc(2,3,2)=dfk0(2,3,2)*f1(2)*f1(3)+fk0(2,3)*df1dc(2)*f1(3)
      dfdc(2,3,3)=dfk0(2,3,3)*f1(2)*f1(3)+fk0(2,3)*f1(2)*df1dc(3)
      dfdc(2,3,4)=dfk0(2,3,4)*f1(2)*f1(3)
c
      dfdc(2,4,1)=dfk0(2,4,1)*f1(2)*f1(4)
      dfdc(2,4,2)=dfk0(2,4,2)*f1(2)*f1(4)+fk0(2,4)*df1dc(2)*f1(4)
      dfdc(2,4,3)=dfk0(2,4,3)*f1(2)*f1(4)
      dfdc(2,4,4)=dfk0(2,4,4)*f1(2)*f1(4)+fk0(2,4)*f1(2)*df1dc(4)
c
      dfdc(3,4,1)=dfk0(3,4,1)*f1(3)*f1(4)
      dfdc(3,4,2)=dfk0(3,4,2)*f1(3)*f1(4)
      dfdc(3,4,3)=dfk0(3,4,3)*f1(3)*f1(4)+fk0(3,4)*df1dc(3)*f1(4)
      dfdc(3,4,4)=dfk0(3,4,4)*f1(3)*f1(4)+fk0(3,4)*f1(3)*df1dc(4)
c
c  derivative of total force constant f(i,j) wrt bond length rbh(k)
c  is given by dfdh(i,j,k)
c
c  only non-zero derivatives are those from rbh(i) and rbh(j)
c
       dfdh(1,2,1)=fk0(1,2)*df1dh(1)*f1(2)
       dfdh(1,2,2)=fk0(1,2)*f1(1)*df1dh(2)
       dfdh(1,2,3)=0.0d0
       dfdh(1,2,4)=0.0d0
c
       dfdh(1,3,1)=fk0(1,3)*df1dh(1)*f1(3)
       dfdh(1,3,2)=0.0d0
       dfdh(1,3,3)=fk0(1,3)*f1(1)*df1dh(3)
       dfdh(1,3,4)=0.0d0
c
       dfdh(1,4,1)=fk0(1,4)*df1dh(1)*f1(4)
       dfdh(1,4,2)=0.0d0
       dfdh(1,4,3)=0.0d0
       dfdh(1,4,4)=fk0(1,4)*f1(1)*df1dh(4)
c
       dfdh(2,3,1)=0.0d0
       dfdh(2,3,2)=fk0(2,3)*df1dh(2)*f1(3)
       dfdh(2,3,3)=fk0(2,3)*f1(2)*df1dh(3)
       dfdh(2,3,4)=0.0d0
c
       dfdh(2,4,1)=0.0d0
       dfdh(2,4,2)=fk0(2,4)*df1dh(2)*f1(4)
       dfdh(2,4,3)=0.0d0
       dfdh(2,4,4)=fk0(2,4)*f1(2)*df1dh(4)
c
       dfdh(3,4,1)=0.0d0
       dfdh(3,4,2)=0.0d0
       dfdh(3,4,3)=fk0(3,4)*df1dh(3)*f1(4)
       dfdh(3,4,4)=fk0(3,4)*f1(3)*df1dh(4)
c
       return
       end
c
c******************************************************
c
c
       subroutine switchf
c
c  calculates switching functions: s3,sphi,stheta
c  and their derivatives ds3,dsphi,dstheta
c
       implicit double precision (a-h,o-z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
c jcc valores originales
C      a1s=1.5313681d-7
C      b1s=-4.6696246d0
C      a2s=1.0147402d-7
C      b2s=-12.363798d0
c jcc valores ajustados a CCSDT 
       a1s=1.5132681d-7
       b1s=-4.3792246d0
       a2s=1.9202402d-7
       b2s=-12.323018d0
c
c  use double precision criterion:
c
c  tanh(19.0d0)=1.0d0
c
       argmax=19.0d0
       do i=1,4
         args1=a1s*(rch(i)-r0ch)*(rch(i)-b1s)**8
         if(args1.lt.argmax)then
           s1(i)=1.0d0-tanh(args1)
           ds1(i)=a1s*((rch(i)-b1s)**8
     *                 +8.0d0*(rch(i)-r0ch)*(rch(i)-b1s)**7)
           ds1(i)=-ds1(i)/cosh(args1)**2
         else
           s1(i)=0.0d0
           ds1(i)=0.0d0
         endif
c
         args2=a2s*(rch(i)-r0ch)*(rch(i)-b2s)**6
         if(args2.lt.argmax)then
           s2(i)=1.0d0-tanh(args2)
           ds2(i)=a2s*((rch(i)-b2s)**6
     *                 +6.0d0*(rch(i)-r0ch)*(rch(i)-b2s)**5)
           ds2(i)=-ds2(i)/cosh(args2)**2
         else
           s2(i)=0.0d0
           ds2(i)=0.0d0
         endif
c
c  calculate s3 and ds3
c
         args3=a3s*(rch(i)-r0ch)*(rch(i)-b3s)**2
         if (args3.lt.argmax)then
           s3(i)=1.0d0-tanh(args3)
           ds3(i)=a3s*(3.0d0*rch(i)**2-2.0d0*rch(i)*(r0ch+2.0d0*b3s)
     *          +b3s*(b3s+2.0d0*r0ch))
           ds3(i)=-ds3(i)/cosh(args3)**2
         else
           s3(i)=0.0d0
           ds3(i)=0.0d0
         endif
c
c  calculate sphi and dsphi
c
c  condition here is on the bondlength rch(i)
c  st argsphi is lt approx 19.0d0
c
         if(rch(i).lt.3.8d0)then
           argsphi=aphi*(rch(i)-r0ch)*exp(bphi*(rch(i)-cphi)**3)
           sphi(i)=1.0d0-tanh(argsphi)
           dsphi(i)=aphi*(1.0d0+3.0d0*bphi*(rch(i)-r0ch)
     *                      *(rch(i)-cphi)**2)
           dsphi(i)=dsphi(i)*exp(bphi*(rch(i)-cphi)**3)
           dsphi(i)=-dsphi(i)/cosh(argsphi)**2
         else
           sphi(i)=0.0d0
           dsphi(i)=0.0d0
         endif
c
c  calculate stheta and dstheta
c
         if(rch(i).lt.3.8d0)then
           argstheta=atheta*(rch(i)-r0ch)*exp(btheta*(rch(i)-ctheta)**3)
           stheta(i)=1.0d0-tanh(argstheta)
           dstheta(i)=atheta*(1.0d0+3.0d0*btheta*(rch(i)-r0ch)
     *           *(rch(i)-ctheta)**2)
           dstheta(i)=dstheta(i)*exp(btheta*(rch(i)-ctheta)**3)
           dstheta(i)=-dstheta(i)/cosh(argstheta)**2
         else
           stheta(i)=0.0d0
           dstheta(i)=0.0d0
         endif
       enddo
       return
       end
C
      SUBROUTINE PREPOT
C
C   N3TMMN = 3 * NATOMS
C   NATOMS = the number of atoms represented by this potential function
C
C   The variable N3TMMN is the minimum value of N3TM allowed to be
C   passed by the calling routine for the number of cartesian
C   coordinates needed to represent the full system represented by this
C   potential energy surface routine.
C   N3TM must be greater than or equal to N3TMMN.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
      PARAMETER (N3TMMN = 18)
C
      COMMON/PT3CM/ EZERO(ISURF+1)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
       common /ndx/     nc(3),nhb(3),nh(4,3)
C
C
C  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
C
      IF(NATOMS.GT.25) THEN
         WRITE(NFLAG(18),1111)
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')
         STOP
      END IF
C
       DO I=1,5
          REF(I) = ' '
       END DO
C
       REF(1)='J. Espinosa-Garcia and J. C. Corchado'
       REF(2)='Unpublished results, 1999.'
C
      INDEXES(1) = 1
      INDEXES(2) = 6
      INDEXES(3) = 1
      INDEXES(4) = 1
      INDEXES(5) = 1
      INDEXES(6) = 1
C
C
C
      IRCTNT=6
C
      CALL POTINFO
C
      CALL ANCVRT
c
c  calculate indexes for coordinates
c
       do ind=1,3
         icount=ind-3
         nc(ind)=3*nnc+icount
         nhb(ind)=3*nnb+icount
         do i=1,4
           nh(i,ind)=3*nnh(i)+icount
         enddo
       enddo
c
c  convert to appropriate units:
c
c  energy   in 1.0d+05 j/mol
c  time     in 1.0d-14 s
c  distance in 1.0d-10 m
c  angles   in radians
c  mass     in amu
c
       fact1=0.041840d0
       fact2=6.022045d0
c
       d1ch=d1ch*fact1
       d3ch=d3ch*fact1
       d1cb=d1cb*fact1
       d3cb=d3cb*fact1
       d1hh=d1hh*fact1
       d3hh=d3hh*fact1
       fch3=fch3*fact2
       hch3=hch3*fact2
       fkinf=fkinf*fact2
       ak=ak*fact2
C
      RETURN
      END
C
      BLOCK DATA PTPACM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER(75) REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT3CM/ EZERO(ISURF+1)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
      DATA NASURF /1,35*0/
      DATA NDER /1/
       DATA NFLAG /1,1,15*0,6,0,0/
C
      DATA ANUZERO /-0.70983241091685123/!/0.0D0/ !YL 
      DATA ICARTR,MSURF,MDER/1,0,1/
      DATA NULBL /25*0/
      DATA NATOMS /6/
       DATA nnc     /2/
       DATA nnb     /6/
       DATA nnh     /3,4,5,1/
       DATA r0ch    /   1.08898d0/
       DATA d1ch    / 111.266d0/
       DATA d3ch    /  48.96226d0/
       DATA a1ch    /   1.78374d0/
       DATA b1ch    /   0.14201d0/
       DATA c1ch    /   2.21773d0/
       DATA r0hh    /   0.74239d0/
       DATA d1hh    / 108.382d0/
       DATA d3hh    /  38.42657d0/
       DATA ahh     /   1.9589d0/
       DATA r0cb    /   1.08898d0/
       DATA d1cb    /  56.505d0/
       DATA d3cb    /  19.612d0/
       DATA acb     /   1.4173200d0/
       DATA a3s     /   0.1475300d0/
       DATA b3s     /  -2.9926300d0/
       DATA aphi    /   0.5307000d0/
       DATA bphi    /   0.4012200d0/
       DATA cphi    /   1.9235100d0/
       DATA atheta  /   0.9119800d0/
       DATA btheta  /   0.3537500d0/
       DATA ctheta  /   1.8970500d0/
       DATA fch3    /   0.0693700d0/
       DATA hch3    /   0.1387400d0/
       DATA fkinf   /   0.4291400d0/
       DATA ak      /   0.1353000d0/
       DATA bk      /  10.7132d0/
       DATA aa1     /   1.265960d0/
       DATA aa2     /   0.000710d0/
       DATA aa3     /   0.985920d0/
       DATA aa4     /   2.785060d0/
C
       END

