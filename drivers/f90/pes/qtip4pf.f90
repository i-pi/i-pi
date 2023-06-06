!
!=============================================================================
!
! File: qtip4pf.f90
! Original authors: S. Habershon, D. Manolopoulos, T. Markland
! Date: 23 June 2014
! Description: Simple test program to call stand-alone
!              q-tip4p/f water model to evaluate potential
!              and derivatives for water clusters. Coordinates are read 
!              from 'hex.xyz' for now.
!              We also check the derivatives numerically - seem OK!
!
! NOTES:
! (1) ATOMIC UNITS ARE USED THROUGHOUT
! (2) WE ASSUME THAT THE ATOMS ARE IN A LIST ARRANGED WITH MOLECULES
!     AS [ {O,H,H},{O,H,H},..... ].
! (3) PBCs are NOT included - clusters only!
! (4) q-TIP4P/F parameters are hard-coded in the "potential" subroutine.
!
! Disclaimer: please thoroughly check the output of this before 
! using it for anything interesting!
!
!=============================================================================
!


!
!====================================================================================
!
! The routines below implement q-TIP4P/F with PBCs. 
!
! Unless you want to change the potential, you shouldn't have to change anything 
! below here...unless I've made a mistake...
!
! Expects atoms in the order O H H O H H O H H ... 
!====================================================================================
!
subroutine qtipfpf_setparameters(rcut, alpha, qo, qh, oo_sig, oo_eps, theta, reoh, apot, bpot, alp, alpha2)
  implicit none
  real(8) rcut, alpha, qo, qh, oo_sig, oo_eps, theta, reoh, apot, bpot, alp, alpha2
  !
  ! Set potential parameters - ATOMIC UNITS!
  !
  rcut = 9.d0 / 0.5291772108d0               ! 9 Angstrom cutoff for LJ
  alpha = 0.73612d0
  qo = -1.1128d0
  qh = +0.5d0 * 1.1128d0
  oo_sig = 5.96946d0
  oo_eps = 2.95147d-4
  theta = 107.4d0 * (dacos(-1.d0) / 180.d0)
  reoh = 1.78d0
  apot = 0.185d0
  bpot = 0.07d0
  alp = 1.21d0
  alpha2 = 0.5d0 * (1.d0 - alpha)
  
end subroutine qtipfpf_setparameters

subroutine qtip4pf(box,rt,na,dvdrt,v, vir)
  implicit none  
  integer :: i,j,na,nm,ic,njump
  real(8) :: r(3,na), dvdr(3,na), box(3), rt(na,3), dvdrt(na,3)
  real(8) :: oo_eps, oo_sig, rcut,v,vlj,vint
  real(8) :: apot, bpot, alp, alpha, alpha2
  real(8) :: qo,qh,theta,reoh,vir(3,3), virlj(3,3), virint(3,3)
  real(8), allocatable :: ro(:,:), dvdrlj(:,:), z(:)

  ! Set and check number of molecules.
  ! 
  nm = na/3
  if (3*nm .ne. na) stop 'ERROR 1 in POTENTIAL'

  call qtipfpf_setparameters(rcut, alpha, qo, qh, oo_sig, oo_eps, theta, reoh, apot, bpot, alp, alpha2)

  ! Zero-out potential and derivatives.
  !
  v = 0.d0
  vir = 0.d0
  dvdr(:,:) = 0.d0

  do i=1, na 
    r(:,i)=rt(i,:)
  enddo


  !
  ! *** COULOMB CALCULATION ***
  !
  ! Allocate space for oxygen position storage.
  !
  allocate (ro(3,nm), z(na) )


  ! Determine the positions of the m-sites.
  !
  ic = 0
  do i = 1, na, 3
     ic = ic + 1
     do j = 1, 3
        ro(j,ic) = r(j,i)
        r(j,i) = alpha * r(j,i) + alpha2*(r(j,i+1)+r(j,i+2))
     enddo
  enddo


  ! Allocate atomic charges and molecular identities.
  !
  ic = 0
  do i = 1, na, 3
     z(i) = qo
     z(i+1) = qh
     z(i+2) = qh
   enddo

  call ewald_nc(box,r,dvdr,v,vir,na,z)
 
  !
  ! Use the chain rule to calculate the correct forces on the atoms. 
  ! Also, replace the m-site with the original oxygen atoms.
  !
  ic = 0
  do j = 1, na, 3
     ic = ic + 1
     do i = 1, 3
        dvdr(i,j+1) = dvdr(i,j+1) + alpha2 * dvdr(i,j) 
        dvdr(i,j+2) = dvdr(i,j+2) + alpha2 * dvdr(i,j) 
        dvdr(i,j) = alpha * dvdr(i,j) 
        r(i,j) = ro(i,ic)
     enddo
  enddo  
  deallocate (ro, z)
  
  !
  ! *** LJ CALCULATION ***
  !

  ! Calculate LJ contribution - note that we set njump to 3 because
  ! we only calculate interactions between every 3rd atom (i.e. oxygen).
  !
  allocate (dvdrlj(3,na))
  njump = 3
  Call lj_basic(box,r,dvdrlj,vlj,virlj,na,njump,oo_eps,oo_sig,rcut)
  dvdr(:,:) = dvdr(:,:) + dvdrlj(:,:)
  v = v + vlj
  vir = vir + virlj
  deallocate ( dvdrlj )

  !
  ! *** INTRAMOLECULAR CALCULATION ***
  !
  Call intra_morse_harm(r,dvdr,na,vint,virint,theta,reoh,apot,bpot,alp)  
  v = v + vint
  vir = vir + virint

  ! ....and we're done...
  !
  vir = -1.0d0 *vir
  do i=1, na 
    dvdrt(i,:) = -dvdr(:,i)
  enddo
  return
end Subroutine qtip4pf

subroutine qtip4pf_sr(rt,na, dvdrt,v, vir)
  implicit none  
  integer :: i,na,nm
  real(8) :: r(3,na), dvdr(3,na), rt(na,3), dvdrt(na,3)
  real(8) :: oo_eps, oo_sig, rcut,v
  real(8) :: apot, bpot, alp, alpha, alpha2
  real(8) :: qo,qh,theta,reoh,vir(3,3)

  ! Set and check number of molecules.
  ! 
  nm = na/3
  if (3*nm .ne. na) stop 'ERROR 1 in POTENTIAL'

  call qtipfpf_setparameters(rcut, alpha, qo, qh, oo_sig, oo_eps, theta, reoh, apot, bpot, alp, alpha2)
    
  ! Zero-out potential and derivatives.
  !
  v = 0.d0
  vir = 0.d0
  dvdr(:,:) = 0.d0

  do i=1, na 
    r(:,i)=rt(i,:)
  enddo

  call intra_morse_harm(r,dvdr,na,v,vir,theta,reoh,apot,bpot,alp)  

  vir = -1.0d0 *vir
  do i=1, na 
    dvdrt(i,:) = -dvdr(:,i)
  enddo
  return
end subroutine

!
!==============================================================================
!
! Lennard Jones calculation - taken from original q-tip4p/f code.
!
!==============================================================================
!
Subroutine lj_basic(box,r,dvdr,v,vir,na,njump,oo_eps,oo_sig,rcut)
  implicit none
  integer :: na,i,j,njump
  real(8) :: r(3,na),dvdr(3,na),v,box(3),vir(3,3)
  real(8) :: oo_eps,oo_sig,rcut,ptail
  real(8) :: sigsq,rcutsq,vij
  real(8) :: drsq,onr2,fij,dfx,dfy,dfz,sr2,sr6,wij
  real(8) :: dx,dy,dz,vscale,dscale

  sigsq = oo_sig*oo_sig
  rcutsq = rcut*rcut
  
  v = 0.d0
  dvdr(:,:) = 0.d0
  vir = 0.d0

  do j = 1+njump,na,njump
     do i = 1,j-njump,njump
        dx = r(1,i)-r(1,j)
        dy = r(2,i)-r(2,j)
        dz = r(3,i)-r(3,j)
        dx = dx - box(1)*nint(dx/box(1))
        dy = dy - box(2)*nint(dy/box(2))
        dz = dz - box(3)*nint(dz/box(3))
        drsq = dx*dx + dy*dy + dz*dz
        if (drsq .lt. rcutsq) then
           onr2 = 1.d0/drsq
           sr2 = sigsq * onr2
           sr6 = sr2 * sr2 * sr2
           vij = sr6 * (sr6-1.d0)
           v = v + vij
           wij = sr6 * (sr6-0.5d0)
           fij = wij * onr2
           dfx = fij * dx
           dfy = fij * dy
           dfz = fij * dz
           dvdr(1,i) = dvdr(1,i) - dfx
           dvdr(2,i) = dvdr(2,i) - dfy
           dvdr(3,i) = dvdr(3,i) - dfz
           dvdr(1,j) = dvdr(1,j) + dfx
           dvdr(2,j) = dvdr(2,j) + dfy
           dvdr(3,j) = dvdr(3,j) + dfz
           vir(1,1) = vir(1,1) - dx * dfx
           vir(2,2) = vir(2,2) - dy * dfy
           vir(3,3) = vir(3,3) - dz * dfz
           vir(1,2) = vir(1,2) - dx * dfy
           vir(1,3) = vir(1,3) - dx * dfz
           vir(2,1) = vir(2,1) - dy * dfx
           vir(2,3) = vir(2,3) - dy * dfz
           vir(3,1) = vir(3,1) - dz * dfx
           vir(3,2) = vir(3,2) - dz * dfy
        endif
     enddo
  enddo

  vscale = 4.d0*oo_eps
  dscale = 48.d0*oo_eps

  v = vscale * v
  do i = 1,na,njump
     dvdr(1,i) = dscale*dvdr(1,i)
     dvdr(2,i) = dscale*dvdr(2,i)
     dvdr(3,i) = dscale*dvdr(3,i)
  enddo
  call pres_lj_tail(oo_eps,oo_sig,rcut,ptail,box,na)
  vir(:,:) = dscale * vir(:,:)
  vir(1,1) = vir(1,1)-ptail
  vir(2,3) = vir(2,2)-ptail
  vir(3,3) = vir(3,3)-ptail
 
  return
end subroutine lj_basic

subroutine pres_lj_tail(oo_eps, oo_sig,rcut, ptail,boxlxyz,na)
  implicit none
  ! ------------------------------------------------------------------
  ! LJ tail correction to pressure
  ! ------------------------------------------------------------------
  integer na,nm
  real(8) ptail,boxlxyz(3),vol,pi,rho,prefac
  real(8) oo_eps,oo_sig,rcut
  
  nm = na/3
  pi = dacos(-1.d0)
  vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
  rho = dble(nm) / vol
  prefac = vol*(16.d0*pi*(rho**2)*oo_eps*(oo_sig**3)) / (3.d0)
  ptail  = prefac*( (2.d0/3.d0)*(oo_sig/rcut)**9 &
       - (oo_sig/rcut)**3)
  
  return
end subroutine pres_lj_tail

!
!======================================================================
!
! Intramolecular contribution to q-TIP4P/F energy and derivatives.
!
! Quartic expansion of morse potential for stretch
! with harmonic bend.
!
!======================================================================
!
Subroutine intra_morse_harm(r,dvdr,na,v,vir,theta,reoh,apot,bpot,alp)
  implicit none
  integer :: na,j
  real(8) :: r(3,na), dvdr(3,na), vir(3,3)
  real(8) :: dr1,dr2,dr3,theta,reoh
  real(8) :: dx1,dy1,dz1,v,dx2,dy2,dz2,dr1sq,dr2sq,dr3dx3
  real(8) :: dx3,dy3,dz3,dr3sq,dr1dx1,dr1dy1,dr1dz1
  real(8) :: dr2dx2,dr2dy2,dr2dz2,dr3dy3,dr3dz3,u,vv,v2,arg,ang
  real(8) :: dang,uprime,vprime,grad
  real(8) :: dvdx1,dvdx2,dvdx3,dvdy1,dvdy2,dvdy3
  real(8) :: dvdz1,dvdz2,dvdz3,dvdr1,dvdr2,dvdr3
  real(8) :: dthetadr1,dthetadr2,dthetadr3
  real(8) :: darg,dr1a,dr2a,dr3a,apot,bpot
  real(8) :: de,alp,alp2,alp3,alp4,drasq,drbsq
  real(8) :: f1,f2,deb,a1,a2,a3,xx,yy,zz,xy,xz,yz

  de = apot
  deb = bpot
  alp2  = alp*alp
  alp3  = alp*alp2
  alp4  = alp3*alp
  f1 = 7.d0 / 12.d0
  f2 = 7.d0 / 3.d0

  v = 0.d0
  vir = 0.0d0
  ! Loop over molecules

  do j = 1,na,3
     dx1 = r(1,j+1)-r(1,j)
     dy1 = r(2,j+1)-r(2,j)
     dz1 = r(3,j+1)-r(3,j)
     dr1sq = dx1*dx1 + dy1*dy1 + dz1*dz1
     dr1 = dsqrt(dr1sq)
     dr1a = dr1
     dr1dx1 = dx1/dr1
     dr1dy1 = dy1/dr1
     dr1dz1 = dz1/dr1
     dr1 = dr1-reoh
     drasq=dr1*dr1

     dx2 = r(1,j+2)-r(1,j)
     dy2 = r(2,j+2)-r(2,j)
     dz2 = r(3,j+2)-r(3,j)
     dr2sq = dx2*dx2 + dy2*dy2 + dz2*dz2
     dr2 = dsqrt(dr2sq)
     dr2a = dr2
     dr2dx2 = dx2/dr2
     dr2dy2 = dy2/dr2
     dr2dz2 = dz2/dr2
     dr2 = dr2-reoh
     drbsq=dr2*dr2
     
     dx3 = r(1,j+2)-r(1,j+1)
     dy3 = r(2,j+2)-r(2,j+1)
     dz3 = r(3,j+2)-r(3,j+1)
     dr3sq = dx3*dx3 + dy3*dy3 + dz3*dz3
     dr3 = dsqrt(dr3sq)
     dr3a = dr3
     dr3dx3 = dx3/dr3
     dr3dy3 = dy3/dr3
     dr3dz3 = dz3/dr3
     
     u = (dr1sq + dr2sq - dr3sq)
     vv = (2.d0 * dr1a * dr2a)         
     v2 = 1.d0/(vv * vv)
     arg = u / vv
     darg = -1.d0/dsqrt(1.d0 - arg*arg)
     ang = dacos( arg )
     dang = (ang - theta)
     uprime = 2.d0 * dr1a
     vprime = 2.d0 * dr2a
     grad = (uprime*vv - vprime*u) * v2
     dthetadr1= darg * grad
     uprime = 2.d0 * dr2a
     vprime = 2.d0 * dr1a
     grad = (uprime*vv - vprime*u) * v2
     dthetadr2= darg * grad
     uprime = -2.d0*dr3a
     grad = (uprime*vv) * v2
     dthetadr3= darg * grad

     v = v+de*(alp2*drasq-alp3*dr1*drasq+f1*alp4*drasq*drasq)
     v = v+de*(alp2*drbsq-alp3*dr2*drbsq+f1*alp4*drbsq*drbsq)
     v = v+deb*dang**2

     a1 = de*(2.d0*alp2*dr1-3.d0*alp3*drasq+f2*alp4*dr1*drasq)
     a2 = de*(2.d0*alp2*dr2-3.d0*alp3*drbsq+f2*alp4*dr2*drbsq)
     a3 = 2.d0*deb*dang

     dvdr1 = a1+a3*dthetadr1
     dvdr2 = a2+a3*dthetadr2
     dvdr3 = a3*dthetadr3
     dvdx1 = dvdr1*dr1dx1
     dvdy1 = dvdr1*dr1dy1
     dvdz1 = dvdr1*dr1dz1 
     dvdx2 = dvdr2*dr2dx2
     dvdy2 = dvdr2*dr2dy2
     dvdz2 = dvdr2*dr2dz2
     dvdx3 = dvdr3*dr3dx3
     dvdy3 = dvdr3*dr3dy3
     dvdz3 = dvdr3*dr3dz3
     dvdr(1,j) = dvdr(1,j) - dvdx1 - dvdx2
     dvdr(2,j) = dvdr(2,j) - dvdy1 - dvdy2
     dvdr(3,j) = dvdr(3,j) - dvdz1 - dvdz2
     dvdr(1,j+1) = dvdr(1,j+1) + dvdx1 - dvdx3
     dvdr(2,j+1) = dvdr(2,j+1) + dvdy1 - dvdy3
     dvdr(3,j+1) = dvdr(3,j+1) + dvdz1 - dvdz3
     dvdr(1,j+2) = dvdr(1,j+2) + dvdx2 + dvdx3
     dvdr(2,j+2) = dvdr(2,j+2) + dvdy2 + dvdy3
     dvdr(3,j+2) = dvdr(3,j+2) + dvdz2 + dvdz3

     xx = dx1*dvdx1 + dx2*dvdx2 + dx3*dvdx3
     xy = dx1*dvdy1 + dx2*dvdy2 + dx3*dvdy3
     xz = dx1*dvdz1 + dx2*dvdz2 + dx3*dvdz3
     yy = dy1*dvdy1 + dy2*dvdy2 + dy3*dvdy3
     yz = dy1*dvdz1 + dy2*dvdz2 + dy3*dvdz3
     zz = dz1*dvdz1 + dz2*dvdz2 + dz3*dvdz3
     vir(1,1) = vir(1,1) + xx
     vir(1,2) = vir(1,2) + xy
     vir(1,3) = vir(1,3) + xz
     vir(2,1) = vir(2,1) + xy
     vir(2,2) = vir(2,2) + yy
     vir(2,3) = vir(2,3) + yz
     vir(3,1) = vir(3,1) + xz
     vir(3,2) = vir(3,2) + yz
     vir(3,3) = vir(3,3) + zz
  enddo
  
  return 
end subroutine intra_morse_harm


subroutine ewald_nc(box,r,dvdr,v,vir,n,q)
  implicit none
  ! ------------------------------------------------------------------
  ! Non-Cubic Ewald sum
  ! ------------------------------------------------------------------
  integer n,kmax
  real(8) box(3),dvdr(3,n),r(3,n),q(n),vir(3,3),v
  real(8) cvar,rshort,rlong,pi,wrcut,walpha,rkmax,rkmax2

  pi = dacos(-1.d0)

  
  ! Non-cubic setup - parameters set-up for non-unit box
  cvar = 1.2d0     
  ! Shortest side length

  rshort = min(box(1),box(2),box(3))
  wrcut = rshort*min(0.5d0,cvar*n**(-1.d0/6.d0))
  walpha = pi/wrcut

  ! Longest side length

  rlong = max(box(1),box(2),box(3))
  rkmax = 2.d0*pi*walpha
  kmax = int(walpha*rlong)
  
  v = 0.d0
  vir(:,:) = 0.d0


  ! Point Charge Ewald sum

  call rwald_basic(box,r,dvdr,v,vir,n,q,wrcut,walpha)
  
  ! Reciprocal Space Ewald Sum

  rkmax2 = rkmax*rkmax
  call kwald_nc(box,r,dvdr,v,vir,n,q,walpha,rkmax2,kmax)

  return
end subroutine ewald_nc

subroutine rwald_basic(box,r,dvdr,v,vir,n,q,rcut,alpha) !(n,mol,q,r,v,vir,dvdr,rcut,alpha,boxlxyz)
  implicit none
  ! ------------------------------------------------------------------
  ! Real space part of the Ewald sum for non-cubic systems
  ! Low precision version with an approximation to erfc(x)
  ! ------------------------------------------------------------------
  integer n,i,j
  real(8) q(n),r(3,n),dvdr(3,n),box(3), vir(3,3)
  real(8) pi,rfac,sfac,rcut,rcutsq,alpha,dx,dy,dz,drsq,dr
  real(8) x,e,t,erfc,du,dur,qij,dv,dvr,dvx,dvy,dvz,v
  real(8) p,a1,a2,a3,a4,a5

  ! parameters in the approximation to erfc(x) [A&S, 7.1.26]

  parameter (p = 3.0525860d0)
  parameter (a1 = 0.254829592d0)
  parameter (a2 = -0.284496736d0)
  parameter (a3 = 1.421413741d0)
  parameter (a4 = -1.453152027d0)
  parameter (a5 = 1.061405429d0)

  ! and evaluate the real space sum

  pi = dacos(-1.d0)
  rfac = alpha/dsqrt(pi)
  sfac = -2.d0*rfac
  rcutsq = rcut*rcut

  ! including contributions from within the same cell

  do i = 1,n
     do j = 1,i-1
        dx = r(1,i)-r(1,j)
        dy = r(2,i)-r(2,j)
        dz = r(3,i)-r(3,j)
        dx = dx - box(1)*nint(dx/box(1))
        dy = dy - box(2)*nint(dy/box(2))
        dz = dz - box(3)*nint(dz/box(3))
        drsq = dx*dx+dy*dy+dz*dz
        if (drsq .lt. rcutsq) then
           dr = dsqrt(drsq)
           x = alpha*dr
           e = dexp(-x*x)
           t = p/(p+x)
           erfc = e*t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
           if ((i-1)/3==(j-1)/3) erfc = erfc-1.d0
           du = erfc/dr
           dur = (sfac*e-du)/drsq
           qij = q(i)*q(j)
           dv = qij*du
           dvr = qij*dur
           dvx = dvr*dx
           dvy = dvr*dy
           dvz = dvr*dz
           v = v+dv
           dvdr(1,i) = dvdr(1,i)+dvx
           dvdr(2,i) = dvdr(2,i)+dvy
           dvdr(3,i) = dvdr(3,i)+dvz
           dvdr(1,j) = dvdr(1,j)-dvx
           dvdr(2,j) = dvdr(2,j)-dvy
           dvdr(3,j) = dvdr(3,j)-dvz
           vir(1,1) = vir(1,1) + dx * dvx
           vir(1,2) = vir(1,2) + dx * dvy
           vir(1,3) = vir(1,3) + dx * dvz
           vir(2,1) = vir(2,1) + dy * dvx
           vir(2,2) = vir(2,2) + dy * dvy
           vir(2,3) = vir(2,3) + dy * dvz
           vir(3,1) = vir(3,1) + dz * dvx
           vir(3,2) = vir(3,2) + dz * dvy
           vir(3,3) = vir(3,3) + dz * dvz
        endif
     enddo
  enddo  
  return
end subroutine rwald_basic

subroutine kwald_nc(box,r,dvdr,v,vir,n,q,alpha,rkmax2,kmax) !(r,z,n,v,vir,dvdr,alpha,rkmax2,kmax,box)
  implicit none
  ! ------------------------------------------------------------------
  ! Reciprocal space part of the Ewald sum for non-cubic systems
  ! ------------------------------------------------------------------
  integer n,kmax,i,kx,ky,kz
  real(8) r(3,n),q(n),dvdr(3,n),box(3), vir(3,3)
  real(8) ckx(n,0:kmax),skx(n,0:kmax)
  real(8) cky(n,-kmax:kmax),sky(n,-kmax:kmax)
  real(8) ckz(n,-kmax:kmax),skz(n,-kmax:kmax)
  real(8) cxy(n),sxy(n),dsdr(6,n)
  real(8) v,alpha,rkmax2,pi,rfac,twopi,xbox,ybox,zbox
  real(8) xlat,ylat,zlat,xi,yi,zi,b,f,rkx,rky,rkz,vs
  real(8) c,s,rk2,sr,si,tr,ti,w,et,xx,yy,zz,xy,xz,yz,term

  pi = dacos(-1.d0)
  rfac = alpha/dsqrt(pi)

  twopi = 2.d0*pi
  xbox = box(1)
  ybox = box(2)
  zbox = box(3)

  xlat = twopi/xbox
  ylat = twopi/ybox
  zlat = twopi/zbox

  ! setup the trigonometric arrays

  do i = 1,n
     ckx(i,0) = q(i)
     skx(i,0) = 0.d0
     xi = xlat*r(1,i)
     c = dcos(xi)
     s = dsin(xi)
     do kx = 1,kmax
        ckx(i,kx) = c*ckx(i,kx-1) - s*skx(i,kx-1)
        skx(i,kx) = s*ckx(i,kx-1) + c*skx(i,kx-1)
     enddo
  enddo
  do i = 1,n
     cky(i,0) = 1.d0
     sky(i,0) = 0.d0
     yi = ylat*r(2,i)
     c = dcos(yi)
     s = dsin(yi)
     do ky = 1,kmax
        cky(i,ky) = c*cky(i,ky-1) - s*sky(i,ky-1)
        sky(i,ky) = s*cky(i,ky-1) + c*sky(i,ky-1)
        cky(i,-ky) = cky(i,ky)
        sky(i,-ky) = -sky(i,ky)
     enddo
  enddo
  do i = 1,n
     ckz(i,0) = 1.d0
     skz(i,0) = 0.d0
     zi = zlat*r(3,i)
     c = dcos(zi)
     s = dsin(zi)
     do kz = 1,kmax
        ckz(i,kz) = c*ckz(i,kz-1) - s*skz(i,kz-1)
        skz(i,kz) = s*ckz(i,kz-1) + c*skz(i,kz-1)
        ckz(i,-kz) = ckz(i,kz)
        skz(i,-kz) = -skz(i,kz)
     enddo
  enddo

  ! and evaluate the reciprocal space sum

  b = 0.25d0/(alpha*alpha)
  f = twopi/(xbox*ybox*zbox)      ! 2pi / V !
  do kx = 0,kmax
     if (kx .eq. 1) f = 2.d0*f
     rkx = xlat*kx
     do ky = -kmax,kmax
        rky = ylat*ky
        do i = 1,n
           cxy(i) = ckx(i,kx)*cky(i,ky)-skx(i,kx)*sky(i,ky)
           sxy(i) = ckx(i,kx)*sky(i,ky)+skx(i,kx)*cky(i,ky)
        enddo
        do kz = -kmax,kmax
           rkz = zlat*kz
           rk2 = rkx*rkx + rky*rky + rkz*rkz
           if (rk2.lt.rkmax2 .and. rk2.ne.0.d0) then
              sr = 0.d0
              si = 0.d0
              do i = 1,n
                 tr = cxy(i)*ckz(i,kz) - sxy(i)*skz(i,kz)
                 ti = cxy(i)*skz(i,kz) + sxy(i)*ckz(i,kz)
                 sr = sr + tr
                 si = si + ti
                 dsdr(1,i) = -rkx*ti
                 dsdr(2,i) =  rkx*tr
                 dsdr(3,i) = -rky*ti
                 dsdr(4,i) =  rky*tr
                 dsdr(5,i) = -rkz*ti
                 dsdr(6,i) =  rkz*tr
              enddo
              w = (f/rk2)*dexp(-b*rk2)
              et = w*(sr*sr+si*si)
              v = v + et
              w = 2.d0*w
              sr = w*sr
              si = w*si
              do i = 1,n
                 dvdr(1,i) = dvdr(1,i) + sr*dsdr(1,i) + si*dsdr(2,i)
                 dvdr(2,i) = dvdr(2,i) + sr*dsdr(3,i) + si*dsdr(4,i)
                 dvdr(3,i) = dvdr(3,i) + sr*dsdr(5,i) + si*dsdr(6,i)
              enddo    

              term = 2.d0*(1.d0/rk2 + b)
              xx = et * (term*rkx*rkx-1.d0)
              xy = et * (term*rkx*rky)
              xz = et * (term*rkx*rkz)
              yy = et * (term*rky*rky-1.d0)
              yz = et * (term*rky*rkz)
              zz = et * (term*rkz*rkz-1.d0)
              vir(1,1) = vir(1,1) + xx
              vir(1,2) = vir(1,2) + xy
              vir(1,3) = vir(1,3) + xz
              vir(2,1) = vir(2,1) + xy
              vir(2,2) = vir(2,2) + yy
              vir(2,3) = vir(2,3) + yz
              vir(3,1) = vir(3,1) + xz
              vir(3,2) = vir(3,2) + yz
              vir(3,3) = vir(3,3) + zz          
           endif
        enddo
     enddo
  enddo

  ! ...minus the self term

  vs = 0.d0
  do i = 1,n
     vs = vs+q(i)*q(i)
  enddo
  vs = rfac*vs
  v = v-vs

  return
end subroutine kwald_nc
