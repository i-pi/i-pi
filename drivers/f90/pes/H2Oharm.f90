



        MODULE H2O_Harm

        IMPLICIT NONE  

        double precision :: bohr= 0.529177210903
        
        CONTAINS

           SUBROUTINE Harm_pot(natoms,atoms,pot,forces,vpars,r0inp)
           
           integer, intent (in) ::  natoms       
           double precision, dimension (natoms,3), intent(in) :: atoms        
           double precision, intent(out) :: pot
           double precision, dimension (natoms,3), intent(out) :: forces
           double precision, dimension (3), intent(in),optional :: r0inp
           integer :: vpars
           double precision :: kspring
           double precision, dimension (3),save :: r0=0.0
           logical, save :: first=.true.

           kspring = 0.1 ! Hartree/Bohr 
           if (first) then
            if (present(r0inp)) then
             r0 = r0inp/bohr
            else
             r0 = atoms(vpars,:)
             first = .false.
            endif
           endif 
          ! r0(1) = 1.54087e+00 / bohr  
          ! r0(2) = 1.62009e+00 / bohr 
          ! r0(3) = 1.49775e+00 / bohr
           forces = 0
           pot = 0

           pot = 1.0/2.0 * kspring * ((atoms(vpars,1) - r0(1))**2.0 & 
           + (atoms(vpars,2) - r0(2))**2.0 + (atoms(vpars,3) - r0(3))**2.0)
           forces(vpars,1) = - kspring * (atoms(vpars,1) - r0(1))
           forces(vpars,2) =  - kspring * (atoms(vpars,2) - r0(2))
           forces(vpars,3) =  - kspring * (atoms(vpars,3) - r0(3))
           print*, forces(vpars,:),atoms(vpars,:),vpars,natoms,first,pot
           END SUBROUTINE Harm_pot

           SUBROUTINE Harm2atom_pot(natoms,atoms,pot,forces)

           integer, intent (in) ::  natoms
           double precision, dimension (natoms,3), intent(in) :: atoms
           double precision, intent(out) :: pot
           double precision, dimension (natoms,3), intent(out) :: forces
           integer :: i,j
           double precision :: kspring
           double precision, dimension (:,:),allocatable, save :: r0
           logical, save :: first = .true.

!           kspring = 0.1 ! Hartree/Bohr
           kspring = 0.05 ! Hartree/Bohr
           if (first) then
            allocate(r0(natoms,3))       
            r0 = atoms
            first = .false.
           endif
           forces = 0
           pot = 0

           kspring = 0.1 ! Hartree/Bohr 
           forces = 0
           pot = 0

           do i=1,3
            do j = 1,natoms
             pot = pot + (atoms(j,i) - r0(j,i))**2.0
            enddo
           enddo
           pot = 1.0/2.0 * kspring * pot

           do i=1,3
            do j=1,natoms
             forces(j,i) = - kspring * (atoms(j,i) - r0(j,i))
            enddo
           enddo
           
           print*, natoms,pot
           do i=1,natoms
             print*, forces(i,:),atoms(i,:)
           enddo  

           END SUBROUTINE Harm2atom_pot


           SUBROUTINE Harm_anysotropic_pot(natoms,atoms,pot,forces)

           implicit none

           integer, intent (in) ::  natoms
           double precision :: theta, phi
           integer :: i,j,n
           double precision, dimension (natoms,3), intent(in) :: atoms
           double precision, intent(out) :: pot
           double precision, dimension (natoms,3), intent(out) :: forces

           double precision, dimension(3) :: kspring
           double precision, dimension (:,:),allocatable, save :: r0,Rmat
           double precision, dimension (:), allocatable :: rotvec

           logical, save :: first = .true.

           kspring = 0.1 ! Hartree/Bohr
           if (first) then
            if (natoms .lt. 2) then
              stop
            endif

            allocate(r0(natoms,3))       
            allocate(Rmat(3,3))
            allocate(rotvec(3))

            r0 = atoms
            first = .false.
            rotvec = r0(2,:) - r0(1,:)
            theta = atan(sqrt(rotvec(1)**2+rotvec(2)**2.0)/rotvec(3))
            phi = atan(rotvec(2)/rotvec(1))

            Rmat(1,1) = sin(theta)*cos(phi)
            Rmat(2,1) = -sin(phi)
            Rmat(3,1) = -sin(theta)*cos(phi)
            Rmat(1,2) = sin(theta)*sin(phi)
            Rmat(2,2) = cos(phi)
            Rmat(3,2) = -cos(theta)*sin(phi)
            Rmat(1,3) = cos(theta)
            Rmat(2,3) = 0.0d0
            Rmat(3,3) = sin(theta)
            
           endif
           forces = 0
           pot = 0

           kspring(1) = 1.0
           kspring(2:3) = 0.1 ! Hartree/Bohr 
           forces = 0
           pot = 0

           do j = 1,natoms
             do i=1,3
               pot = pot + kspring(i)*sum(Rmat(i,:)*(atoms(j,:) - r0(j,:)))**2.0
             enddo
           enddo
           pot = 1.0/2.0 * pot

           do j=1,natoms
            do i=1,3
             forces(j,n) = forces(j,n) - kspring(i) * Rmat(i,n)*sum(Rmat(i,:)*(atoms(j,:) - r0(j,:)))
            enddo
           enddo
           
           print*, natoms,pot
           do j=1,natoms
             print*, theta,phi,first,forces(j,:),atoms(j,:),r0,sum(Rmat(i,:)*(atoms(j,:) - r0(j,:))), &
              sum(Rmat(2,:)*(atoms(j,:) - r0(j,:))),sum(Rmat(3,:)*(atoms(j,:) - r0(j,:)))
           enddo  

           END SUBROUTINE Harm_anysotropic_pot

           SUBROUTINE Harm2atom_pot_legacy(natoms,atoms,pot,forces)

           integer, intent (in) ::  natoms
           double precision, dimension (natoms,3), intent(in) :: atoms
           double precision, intent(out) :: pot
           double precision, dimension (natoms,3), intent(out) :: forces
           integer, parameter :: nfix = 2
           integer :: i,j
           integer, dimension(nfix) :: nfixlist=(/ 1,79 /)
           double precision :: kspring
           double precision, dimension (3,nfix) :: r0

           kspring = 0.01 ! Hartree/Bohr 
           r0(:,1) = (/  9.93635e-01 ,  1.13988e+00 , 1.41210e+00 /) /bohr 
           r0(:,2) = (/  5.68719e+00 ,  9.32340e+00-15.4799892 ,  5.88333e+00 /) / bohr
           forces = 0
           pot = 0

           do i=1,3
            do j = 1,nfix
             pot = pot + (atoms(nfixlist(j),i) - r0(i,j))**2.0 
            enddo
           enddo  
           pot = 1.0/2.0 * kspring * pot

           do i=1,3
            do j=1,nfix
             forces(nfixlist(j),i) = - kspring * (atoms(nfixlist(j),i) - r0(i,j))
            enddo
           enddo

           print*, forces(1,:),forces(79,:),atoms(1,:),atoms(79,:),natoms,nfixlist,pot
           END SUBROUTINE Harm2atom_pot_legacy

          SUBROUTINE Harm_z_pot(natoms,atoms,pot,forces,vpars)

           integer, intent (in) ::  natoms
           double precision, dimension (natoms,3), intent(in) :: atoms
           double precision, intent(out) :: pot
           double precision, dimension (natoms,3), intent(out) :: forces
           integer :: i,j
           double precision :: kspring
           double precision, dimension (:,:),allocatable, save :: r0
           double precision, dimension (2) :: vpars
           logical, save :: first = .true.

           if (first) then
            allocate(r0(natoms,3))       
            r0(:,1:2) = atoms(:,1:2)
            r0(1:natoms/2,3) = vpars(1)/bohr 
            r0(natoms/2+1:natoms,3) = vpars(2)/bohr
            first = .false.
           endif
           forces = 0
           pot = 0

           kspring = 0.01 ! Hartree/Bohr 
           forces = 0
           pot = 0

           do j = 1,natoms
             pot = pot + (atoms(j,3) - r0(j,3))**2.0
           enddo
           pot = 1.0/2.0 * kspring * pot

           do j=1,natoms
             forces(j,3) = - kspring * (atoms(j,3) - r0(j,3))
           enddo
           
!           print*, natoms,pot
!           do i=1,natoms
!             print*, forces(i,:),atoms(i,:),r0(i,:),vpars
!           enddo  

          END SUBROUTINE Harm_z_pot

       END MODULE    
