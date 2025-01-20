! The main program which runs our driver test case potentials
!
! Copyright (C) 2013, Joshua More and Michele Ceriotti
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!
! Currently the potentials implemented are the Lennard-Jones
! potential, the Silvera-Goldman para-hydrogen potential and
! the ideal gas (i.e. no interaction at all)

      PROGRAM DRIVER
         USE LJ
         USE LJPolymer
         USE SG
         USE PSWATER
         USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer, f_sleep
         USE DISTANCE, only: CELL_VOLUME
      IMPLICIT NONE

      ! SOCKET VARIABLES
      INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
      INTEGER socket, inet, port        ! socket ID & address of the server
      CHARACTER(LEN=1024) :: host, sockets_prefix="/tmp/ipi_"

      ! COMMAND LINE PARSING
      CHARACTER(LEN=1024) :: cmdbuffer
      INTEGER ccmd, vstyle, vseed
      INTEGER, ALLOCATABLE :: seed(:)
      INTEGER verbose
      INTEGER commas(4), par_count      ! stores the index of commas in the parameter string
      DOUBLE PRECISION vpars(6)         ! array to store the parameters of the potential
      
      ! SOCKET COMMUNICATION BUFFERS
      CHARACTER(LEN=12) :: header
      LOGICAL :: isinit=.false., hasdata=.false.
      INTEGER cbuf, rid, length
      CHARACTER(LEN=65536) :: initbuffer      ! it's unlikely a string this large will ever be passed...      
      CHARACTER(LEN=65536) :: string,string1,string2,string3,trimmed  ! it's unlikely a string this large will ever be passed...
      CHARACTER(LEN=30000) :: longbuffer, longstring ! used in water_dip_pol model to pass dipole-z derivative and polarizability
      DOUBLE PRECISION, ALLOCATABLE :: msgbuffer(:)

      ! PARAMETERS OF THE SYSTEM (CELL, ATOM POSITIONS, ...)
      DOUBLE PRECISION sigma, eps, rc, rn, ks ! potential parameters
      DOUBLE PRECISION stiffness ! lennard-jones polymer
      DOUBLE PRECISION sleep_seconds
      INTEGER n_monomer ! lennard-jones polymer
      INTEGER nat
      DOUBLE PRECISION pot, dpot, dist
      DOUBLE PRECISION, ALLOCATABLE :: atoms(:,:), forces(:,:), datoms(:,:)
      DOUBLE PRECISION cell_h(3,3), cell_ih(3,3), virial(3,3), mtxbuf(9), dip(3), charges(3), dummy(3,3,3), vecdiff(3)
      DOUBLE PRECISION, ALLOCATABLE :: friction(:,:)
      DOUBLE PRECISION volume
      DOUBLE PRECISION, PARAMETER :: fddx = 1.0d-5

      DOUBLE PRECISION, ALLOCATABLE :: dipz_der(:, :) ! Dipole (z-component) derivative (water_dip_pol model)
      DOUBLE PRECISION :: pol(3, 3) !Polarizability (water_dip_pol model)

      ! NEIGHBOUR LIST ARRAYS
      INTEGER, DIMENSION(:), ALLOCATABLE :: n_list, index_list
      DOUBLE PRECISION init_volume, init_rc ! needed to correctly adjust the cut-off radius for variable cell dynamics
      DOUBLE PRECISION, ALLOCATABLE :: last_atoms(:,:) ! Holds the positions when the neighbour list is created
      DOUBLE PRECISION displacement ! Tracks how far each atom has moved since the last call of nearest_neighbours

      ! DMW
      DOUBLE PRECISION efield(3)
      INTEGER i, j
      
      ! parse the command line parameters
      ! intialize defaults
      ccmd = 0
      inet = 1
      host = "localhost"//achar(0)
      port = 31415
      verbose = 0
      par_count = 0
      vstyle = -1
      rc = 0.0d0
      init_rc = 0.0d0
      volume = 0.0d0
      init_volume = 0.0d0

      DO i = 1, COMMAND_ARGUMENT_COUNT()
         CALL GET_COMMAND_ARGUMENT(i, cmdbuffer)
         IF (cmdbuffer == "-u") THEN ! flag for unix socket
            inet = 0
            ccmd = 0
         ELSEIF (cmdbuffer == "-h") THEN ! read the hostname (deprecated)
            ccmd = 1
         ELSEIF (cmdbuffer == "-a") THEN ! read the hostname (address)
            ccmd = 1
         ELSEIF (cmdbuffer == "-p") THEN ! reads the port number
            ccmd = 2
         ELSEIF (cmdbuffer == "-m") THEN ! reads the style of the potential function
            ccmd = 3
         ELSEIF (cmdbuffer == "-o") THEN ! reads the parameters
            ccmd = 4
         ELSEIF (cmdbuffer == "-S") THEN ! reads the socket prefix
            ccmd = 5
         ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            verbose = 1
         ELSEIF (cmdbuffer == "-vv") THEN ! flag for verbose standard output
            verbose = 2
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) " Unrecognized command line argument", ccmd
               CALL helpmessage
               STOP "ENDED"
            ENDIF
            IF (ccmd == 1) THEN
               host = trim(cmdbuffer)//achar(0)
            ELSEIF (ccmd == 2) THEN
               READ(cmdbuffer,*) port
            ELSEIF (ccmd == 3) THEN
               IF (verbose>0) THEN  
                  WRITE(*,*) "Running potential type ", trim(cmdbuffer)
               ENDIF 
               IF (trim(cmdbuffer) == "lj") THEN
                  vstyle = 1
               ELSEIF (trim(cmdbuffer) == "sg") THEN
                  vstyle = 2
               ELSEIF (trim(cmdbuffer) == "harm") THEN
                  vstyle = 3
               ELSEIF (trim(cmdbuffer) == "harm3d") THEN
                  vstyle = 40
               ELSEIF (trim(cmdbuffer) == "morse") THEN
                  vstyle = 4
               ELSEIF (trim(cmdbuffer) == "zundel") THEN
                  vstyle = 5
               ELSEIF (trim(cmdbuffer) == "qtip4pf") THEN
                  vstyle = 6
               ELSEIF (trim(cmdbuffer) == "linear") THEN
                  vstyle = 7
               ELSEIF (trim(cmdbuffer) == "pswater") THEN
                  vstyle = 8
               ELSEIF (trim(cmdbuffer) == "lepsm1") THEN
                  vstyle = 9
               ELSEIF (trim(cmdbuffer) == "lepsm2") THEN
                  vstyle = 10
               ELSEIF (trim(cmdbuffer) == "qtip4pf-efield") THEN
                  vstyle = 11
               ELSEIF (trim(cmdbuffer) == "eckart") THEN
                  vstyle = 20
               ELSEIF (trim(cmdbuffer) == "ch4hcbe") THEN
                  vstyle = 21
               ELSEIF (trim(cmdbuffer) == "ljpolymer") THEN
                  vstyle = 22
               ELSEIF (trim(cmdbuffer) == "MB") THEN
                  vstyle = 23
               ELSEIF (trim(cmdbuffer) == "doublewell") THEN
                  vstyle = 25
               ELSEIF (trim(cmdbuffer) == "doublewell_1D") THEN
                  vstyle = 24
               ELSEIF (trim(cmdbuffer) == "morsedia") THEN
                  vstyle = 26
               ELSEIF (trim(cmdbuffer) == "qtip4pf-sr") THEN
                  vstyle = 27
               ELSEIF (trim(cmdbuffer) == "harmonic_bath") THEN
                  vstyle = 28
               ELSEIF (trim(cmdbuffer) == "meanfield_bath") THEN
                  vstyle = 29
               ELSEIF (trim(cmdbuffer) == "water_dip_pol") THEN
                  vstyle = 31
               ELSEIF (trim(cmdbuffer) == "ljmix") THEN
                  vstyle = 32
               ELSEIF (trim(cmdbuffer) == "qtip4pf-c-1") THEN
                  vstyle = 60
               ELSEIF (trim(cmdbuffer) == "qtip4pf-c-2") THEN
                  vstyle = 61
               ELSEIF (trim(cmdbuffer) == "qtip4pf-c-json") THEN
                  vstyle = 62
               ELSEIF (trim(cmdbuffer) == "qtip4pf-c-1-delta") THEN
                  vstyle = 63
               ELSEIF (trim(cmdbuffer) == "qtip4pf-c-2-delta") THEN
                  vstyle = 64
               ELSEIF (trim(cmdbuffer) == "qtip4pf-c-json-delta") THEN
                  vstyle = 65
               ELSEIF (trim(cmdbuffer) == "noo3-h2o") THEN
                  vstyle = 70
               ELSEIF (trim(cmdbuffer) == "gas") THEN
                  vstyle = 0  ! ideal gas
               ELSEIF (trim(cmdbuffer) == "dummy") THEN
                  vstyle = 99 ! returns non-zero but otherwise meaningless values
               ELSE
                  WRITE(*,*) " Unrecognized potential type ", trim(cmdbuffer)
                  WRITE(*,*) " Use -m [dummy|gas|lj|sg|harm|harm3d|morse|morsedia|zundel|qtip4pf|pswater|lepsm1|lepsm2|qtip4pf-efield|eckart|ch4hcbe|ljpolymer|MB|doublewell|doublewell_1D|water_dip_pol|harmonic_bath|meanfield_bath|ljmix|qtip4pf-sr|qtip4pf-c-1|qtip4pf-c-2|qtip4pf-c-json|qtip4pf-c-1-delta|qtip4pf-c-2-delta|qtip4pf-c-json-delta|noo3-h2o] "
                  STOP "ENDED"
               ENDIF
            ELSEIF (ccmd == 4) THEN
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vpars(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vpars(par_count)
            ELSEIF (ccmd == 5) THEN
               sockets_prefix = trim(cmdbuffer)//achar(0)
            ENDIF
            ccmd = 0
         ENDIF
      ENDDO

      IF (vstyle == -1) THEN
         WRITE(*,*) " Error, type of potential not specified."
         CALL helpmessage
         STOP "ENDED"
      ELSEIF (0 == vstyle) THEN
         IF (par_count == 0) THEN
            sleep_seconds = 0.0
         ELSEIF (par_count == 1) THEN
            sleep_seconds = vpars(1)
         ELSE
            WRITE(*,*) "Error: only an optional delay parameters needed for ideal gas."
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (99 == vstyle) THEN
         IF (par_count == 0) THEN
            sleep_seconds = 0.0
         ELSEIF (par_count == 1) THEN
            sleep_seconds = vpars(1)
         ELSE
            WRITE(*,*) "Error: only an optional delay parameters needed for dummy output."
            STOP "ENDED"
         ENDIF
         CALL RANDOM_SEED(size=vseed)
         ALLOCATE(seed(vseed))
         seed = 12345
         CALL RANDOM_SEED(put=seed)
         isinit = .true.
      ELSEIF (6 == vstyle .OR. 27 == vstyle) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error:  no initialization string needed for qtip4pf or qtip4p-sr."
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (11== vstyle) THEN
         IF (par_count .ne. 3) THEN
            WRITE(*,*) "Error:  incorrect initialization string included for qtip4pf-efield. &
     &             Provide the three components of the electric field in V/nm"
            STOP "ENDED"
         ELSE
            ! We take in an electric field in volts / nm.This must be converted to Eh / (e a0).
            do i=1,3
             efield(i) = vpars(i) / 5.14220652d2
            enddo
         ENDIF
         isinit = .true.
      ELSEIF (5 == vstyle) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for zundel."
            STOP "ENDED"
         ENDIF
         CALL prezundelpot()
         CALL prezundeldip()
         isinit = .true.
      ELSEIF (21 == vstyle) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for CH4+H CBE potential."
            STOP "ENDED"
         ENDIF
         CALL prepot()
         isinit = .true.
      ELSEIF (4 == vstyle) THEN
         IF (par_count == 0) THEN ! defaults (OH stretch)
            vpars(1) = 1.8323926 ! r0
            vpars(2) = 0.18748511263179304 ! D
            vpars(3) = 1.1562696428501682 ! a
         ELSEIF ( 2/= par_count) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For morse potential use -o r0,D,a (in a.u.) "
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (20 == vstyle) THEN !eckart
         IF (par_count == 0) THEN ! defaults values
            vpars(1) = 0.0d0
            vpars(2) = 0.66047
            vpars(3) = (6*12)/( 1836 * (vpars(2)**2) *( (4.D0 * ATAN(1.0d0) )**2 ) )
            vpars(4) = 1836*(3800.0d0/219323d0)**2
         ELSEIF ( 4/= par_count) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For eckart potential use  AA,A,B,k"
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (23 == vstyle) THEN !MB
         IF (par_count == 0) THEN ! defaults values
            vpars(1) = 0.004737803248674678
         ELSEIF ( 1/= par_count) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For MB potential up to 1 param can be specified"
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (28 == vstyle) THEN !harmonic_bath
         WRITE(*,*) "This driver implementation is deprecated. Please use the python driver version "
         STOP "ENDED"
         IF (par_count == 3) THEN ! defaults values 
            vpars(4) = 0
            vpars(5) = 0
            vpars(6) = 1
         ELSEIF (par_count /= 6) THEN 
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For harmonic bath use <bath_type> <friction (atomic units)> <omega_c (invcm)> eps(a.u.) delta (a.u.) deltaQ(a.u.)"
            WRITE(*,*) "Available bath_type are: "
            WRITE(*,*) "1 = Ohmic "
            STOP "ENDED"
         ENDIF
         IF (vpars(1) /= 1) THEN
             WRITE(*,*) "Only Ohmic bath implemented"
             STOP "ENDED"
         END IF
         vpars(3) = vpars(3) * 4.5563353e-06 !Change omega_c from invcm to a.u.
         isinit = .true.
      ELSEIF (29 == vstyle) THEN !meanfield bath
         WRITE(*,*) "This driver implementation is deprecated. Please use the python driver version "
         STOP "ENDED"
         IF (par_count == 3) THEN ! defaults values 
            vpars(2) = 0
            vpars(3) = 0
            vpars(4) = 1
         ELSEIF (par_count /= 4) THEN 
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For harmonic meanfield bath use  <friction (atomic units)> eps(a.u.) delta (a.u.) deltaQ(a.u.)"
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (22 == vstyle .or. 32 == vstyle) THEN !ljpolymer or ljmix
         IF (4/= par_count) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For ljpolymer and ljmix potential use n_monomer,sigma,epsilon,cutoff"
            STOP "ENDED"
         ELSE
            n_monomer = nint(vpars(1))
            sigma = vpars(2)
            eps = vpars(3)
            rc = vpars(4)
            rn = rc * 1.2d0
            stiffness = 36.d0 * (2.d0 ** (2.d0/3.d0))*eps
            isinit = .true.
         ENDIF
      ELSEIF (vstyle == 8) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for Partridge-Schwenke H2O potential."
            STOP "ENDED"
         END IF
      ELSEIF (vstyle == 9) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for LEPSM1."
            STOP "ENDED"
         END IF
      ELSEIF (vstyle == 10) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for LEPSM2."
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (vstyle == 11) THEN
         IF (par_count .ne. 3) THEN
            WRITE(*,*) "Error:  incorrect initialization string included for qtip4pf-efield. &
     &    Provide the three components of the electric field in V/nm"
            STOP "ENDED"
         ELSE
            ! We take in an electric field in volts / nm.This must be converted
            ! to Eh / (e a0).
            do i=1,3
             efield(i) = vpars(i) / 5.14220652d2
            enddo
         ENDIF
         isinit = .true.
      ELSEIF (1 == vstyle) THEN
         IF (par_count /= 3) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For LJ potential use -o sigma,epsilon,cutoff "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         sigma = vpars(1)
         eps = vpars(2)
         rc = vpars(3)
         rn = rc*1.2
         isinit = .true.
      ELSEIF (vstyle == 2) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For SG potential use -o cutoff "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         rc = vpars(1)
         rn = rc*1.2
         isinit = .true.
      ELSEIF (vstyle == 3) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For 1D harmonic potential use -o k "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         ks = vpars(1)
         isinit = .true.
      ELSEIF (vstyle == 40) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For 3D harmonic potential use -o k "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         ks = vpars(1)  !è la k dell'ho, unica chiaramente perché in 1D
         isinit = .true.
      ELSEIF (vstyle == 7) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For a linear potential use -o k "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         ks = vpars(1)
         isinit = .true.

      ELSEIF (25 == vstyle) THEN !doublewell
         IF ( par_count /= 0 ) THEN
                 WRITE(*,*) "Error: no initialization string needed for doublewell."
            STOP "ENDED"
         ENDIF
         isinit = .true.

      ELSEIF (24 == vstyle) THEN !doublewell_1D
         IF ( par_count /= 0 ) THEN
                 WRITE(*,*) "Error: no initialization string needed for 1-dimensional doublewell."
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (26 == vstyle) THEN
         IF (par_count == 0) THEN ! defaults (OH stretch)
             vpars(1) = 1.8323926 ! r0
             vpars(2) = 0.18748511263179304 ! D
             vpars(3) = 1.1562696428501682 ! a
         ELSEIF ( 2/= par_count) THEN
             WRITE(*,*) "Error: parameters not initialized correctly."
             WRITE(*,*) "For morse potential use -o r0,D,a (in a.u.) "
             STOP "ENDED"
         ENDIF
      ELSEIF (vstyle == 31) THEN !water dipole and polarizability
         IF (par_count == 0) THEN
            vpars(1) = 1
         ELSEIF (par_count /= 1 .OR. (vpars(1) /= 0 .AND. vpars(1) /= 1)) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For water_dip_pol use only -o 0 or -o 1."
            vpars(1) = 1
         ENDIF
         isinit = .true.
      ENDIF


      IF (verbose > 0) THEN
         WRITE(*,*) " DRIVER - Connecting to host ", trim(host)
         IF (inet > 0) THEN
            WRITE(*,*) " on port ", port, " using an internet socket."
         ELSE
            WRITE(*,*) " using an UNIX socket."
         ENDIF
      ENDIF

      ! Calls the interface to the POSIX sockets library to open a communication channel
      CALL open_socket(socket, inet, port, host, sockets_prefix)
      nat = -1
      DO WHILE (.true.) ! Loops forever (or until the wrapper ends!)

         ! Reads from the socket one message header
         CALL readbuffer(socket, header, MSGLEN)
         IF (verbose > 0) WRITE(*,*) " Message from server: ", trim(header)

         IF (trim(header) == "STATUS") THEN
            ! The wrapper is inquiring on what we are doing
            IF (.not. isinit) THEN
               CALL writebuffer(socket,"NEEDINIT    ",MSGLEN)  ! Signals that we need initialization data
               IF (verbose > 1) WRITE(*,*) "    !write!=> ", "NEEDINIT    "
            ELSEIF (hasdata) THEN
               CALL writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces
               IF (verbose > 1) WRITE(*,*) "    !write!=> ", "HAVEDATA    "
            ELSE
               CALL writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
               IF (verbose > 1) WRITE(*,*) "    !write!=> ", "READY       "
            ENDIF
         ELSEIF (trim(header) == "INIT") THEN     ! The driver is kindly providing a string for initialization
            CALL readbuffer(socket, rid)
            IF (verbose > 1) WRITE(*,*) "    !read!=> RID: ", rid
            CALL readbuffer(socket, cbuf)
            IF (verbose > 1) WRITE(*,*) "    !read!=> init_length: ", cbuf
            CALL readbuffer(socket, initbuffer, cbuf)
            IF (verbose > 1) WRITE(*,*) "    !read!=> init_string: ", cbuf
            IF (verbose > 0) WRITE(*,*) " Initializing system from wrapper, using ", trim(initbuffer)
            isinit=.true. ! We actually do nothing with this string, thanks anyway. Could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver
         ELSEIF (trim(header) == "POSDATA") THEN  ! The driver is sending the positions of the atoms. Here is where we do the calculation!
            
            ! Parses the flow of data from the socket
            CALL readbuffer(socket, mtxbuf, 9)  ! Cell matrix
            IF (verbose > 1) WRITE(*,*) "    !read!=> cell: ", mtxbuf
            cell_h = RESHAPE(mtxbuf, (/3,3/))
            CALL readbuffer(socket, mtxbuf, 9)  ! Inverse of the cell matrix (so we don't have to invert it every time here)
            IF (verbose > 1) WRITE(*,*) "    !read!=> cell-1: ", mtxbuf
            cell_ih = RESHAPE(mtxbuf, (/3,3/))

            ! The wrapper uses atomic units for everything, and row major storage.
            ! At this stage one should take care that everything is converted in the
            ! units and storage mode used in the driver.
            cell_h = transpose(cell_h)
            cell_ih = transpose(cell_ih)
            ! We compute for a generic cell, just in case (even though usually i-PI passes an upper triangular cell-vector matrix)
            volume = CELL_VOLUME(cell_h) !cell_h(1,1)*cell_h(2,2)*cell_h(3,3)

            CALL readbuffer(socket, cbuf)       ! The number of atoms in the cell
            IF (verbose > 1) WRITE(*,*) "    !read!=> cbuf: ", cbuf
            IF (nat < 0) THEN  ! Assumes that the number of atoms does not change throughout a simulation, so only does this once
               nat = cbuf
               IF (verbose > 0) WRITE(*,*) " Allocating buffer and data arrays, with ", nat, " atoms"
               ALLOCATE(msgbuffer(3*nat))
               ALLOCATE(atoms(nat,3), datoms(nat,3))
               ALLOCATE(forces(nat,3))

               IF (vstyle==24 .or. vstyle==25) THEN
                  ALLOCATE(friction(3*nat,3*nat))
                  friction = 0.0d0
               ENDIF
               atoms = 0.0d0
               datoms = 0.0d0
               forces = 0.0d0
               msgbuffer = 0.0d0
               IF (verbose > 1) WRITE(*,*) " Allocation successful "
            ENDIF            

            CALL readbuffer(socket, msgbuffer, nat*3)
            IF (verbose > 1) WRITE(*,*) "    !read!=> positions: ", msgbuffer(0:2), " ..."
            DO i = 1, nat
               atoms(i,:) = msgbuffer(3*(i-1)+1:3*i)
            ENDDO

            IF (vstyle == 0) THEN   ! ideal gas, so no calculation done
               IF (sleep_seconds > 0) THEN
                  ! artificial delay
                  CALL f_sleep(sleep_seconds)
               ENDIF
               pot = 0
               forces = 0.0d0
               virial = 1.0d-200
               ! returns a tiny but non-zero stress, so it can
               ! bypass the check for zero virial that is used
               ! to avoid running constant-pressure simulations
               ! with a code that cannot compute the virial
            ELSEIF (vstyle == 99) THEN ! dummy output, useful to test that i-PI "just runs"
               IF (sleep_seconds > 0) THEN
                  ! artificial delay
                  CALL f_sleep(sleep_seconds)
               ENDIF
               call random_number(pot)
               pot = pot - 0.5
               call random_number(forces)
               forces = forces - 0.5
               call random_number(virial)
               virial = virial - 0.5
               call random_number(dip)
               dip = dip - 0.5 
            ELSEIF (vstyle == 3) THEN ! 1D harmonic potential, so only uses the first position variable
               pot = 0.5*ks*atoms(1,1)**2
               forces = 0.0d0
               forces(1,1) = -ks*atoms(1,1)
               virial = 0.0d0
               virial(1,1) = forces(1,1)*atoms(1,1)
            ELSEIF (vstyle == 40) THEN ! 3D harmonic potential
                   pot = 0.0d0
                   forces = 0.0d0
                   virial = 0.0d0
               DO i=1,nat
                   pot = pot + 0.5*ks*atoms(i,1)**2 + 0.5*ks*atoms(i,2)**2 + 0.5*ks*atoms(i,3)**2
                   forces(i,1) = -ks*atoms(i,1)
                   forces(i,2) = -ks*atoms(i,2)
                   forces(i,3) = -ks*atoms(i,3)
                   virial(1,1) = virial(1,1) + forces(i,1)*atoms(i,1)
                   virial(1,2) = virial(1,2) + forces(i,1)*atoms(i,2)
                   virial(1,3) = virial(1,3) + forces(i,1)*atoms(i,3)
                   virial(2,1) = virial(2,1) + forces(i,2)*atoms(i,1)
                   virial(2,2) = virial(2,2) + forces(i,2)*atoms(i,2)
                   virial(2,3) = virial(2,3) + forces(i,2)*atoms(i,3)
                   virial(3,1) = virial(3,1) + forces(i,3)*atoms(i,1)
                   virial(3,2) = virial(3,2) + forces(i,3)*atoms(i,2)
                   virial(3,3) = virial(3,3) + forces(i,3)*atoms(i,3)
                enddo
            ELSEIF (vstyle == 7) THEN ! linear potential in x position of the 1st atom
               pot = ks*atoms(1,1)
               forces = 0.0d0
               virial = 0.0d0
               forces(1,1) = -ks
               virial(1,1) = forces(1,1)*atoms(1,1)
            ELSEIF (vstyle == 4) THEN ! Morse potential.
               IF (nat/=1) THEN
                  WRITE(*,*) "Expecting 1 atom for 3D Morse (use the effective mass for the atom mass to get proper frequency!) "
                  STOP "ENDED"
               ENDIF
               CALL getmorse(vpars(1), vpars(2), vpars(3), atoms, pot, forces)
            ELSEIF (vstyle == 26) THEN ! Morse potential with 2 atoms.
               IF (nat/=2) THEN
                   WRITE(*,*) "Expecting 2 atom for 3D Morse "
                   STOP "ENDED"
               ENDIF
               CALL getmorsedia(vpars(1), vpars(2), vpars(3), atoms, pot, forces)
            ELSEIF (vstyle == 5) THEN ! Zundel potential.
               IF (nat/=7) THEN
                  WRITE(*,*) "Expecting 7 atoms for Zundel potential, O O H H H H H "
                  STOP "ENDED"
               ENDIF

               CALL zundelpot(pot,atoms)
               CALL zundeldip(dip,atoms)

               datoms=atoms
               DO i=1,7  ! forces by finite differences
                  DO j=1,3
                     datoms(i,j)=atoms(i,j)+fddx
                     CALL zundelpot(dpot, datoms)
                     datoms(i,j)=atoms(i,j)-fddx
                     CALL zundelpot(forces(i,j), datoms)
                     datoms(i,j)=atoms(i,j)
                     forces(i,j)=(forces(i,j)-dpot)/(2*fddx)
                  ENDDO
               ENDDO
               ! do not compute the virial term

           ELSEIF (vstyle == 21) THEN ! CBE CH4+H potential.
               IF (nat/=6) THEN
                  WRITE(*,*) "Expecting 6 atoms for CH4+H potential, H, C, H, H, H, H "
                  WRITE(*,*) "The expected order is such that atoms 1 to 5 are reactant_1 (CH4)"
                  WRITE(*,*) "and atom 6 is reactant_2 ( H 'free') "
                  STOP "ENDED"
               ENDIF

               CALL ch4hpot_inter(atoms, pot)
               datoms=atoms
               DO i=1,6  ! forces by finite differences
                  DO j=1,3
                     datoms(i,j)=atoms(i,j)+fddx
                     CALL ch4hpot_inter(datoms, dpot)
                     datoms(i,j)=atoms(i,j)-fddx
                     CALL ch4hpot_inter(datoms, forces(i,j))
                     datoms(i,j)=atoms(i,j)
                     forces(i,j)=(forces(i,j)-dpot)/(2*fddx)
                  ENDDO
               ENDDO
               ! do not compute the virial term

            ELSEIF (vstyle == 6) THEN ! qtip4pf potential.
               IF (verbose > 1) WRITE(*,*) "TIP4Pf potential"               
               IF (mod(nat,3)/=0) THEN
                  WRITE(*,*) " Expecting water molecules O H H O H H O H H but got ", nat, "atoms"
                  STOP "ENDED"
               ENDIF
               vpars(1) = cell_h(1,1)
               vpars(2) = cell_h(2,2)
               vpars(3) = cell_h(3,3)
               IF (cell_h(1,2).gt.1d-10 .or. cell_h(1,3).gt.1d-12  .or. cell_h(2,3).gt.1d-12) THEN
                  WRITE(*,*) " qtip4pf PES only works with orthorhombic cells", cell_h(1,2), cell_h(1,3), cell_h(2,3)
                  STOP "ENDED"
               ENDIF
               CALL qtip4pf(vpars(1:3),atoms,nat,forces,pot,virial)
               IF (verbose > 1) WRITE(*,*) "TIP4Pf potential computed"
               dip(:) = 0.0
               DO i=1, nat, 3
                  dip = dip -1.1128d0 * atoms(i,:) + 0.5564d0 * (atoms(i+1,:) + atoms(i+2,:))
               ENDDO
            ELSEIF (vstyle == 27) THEN ! short-range qtip4pf potential.
               IF (mod(nat,3)/=0) THEN
                  WRITE(*,*) " Expecting water molecules O H H O H H O H H but got ", nat, "atoms"
                  STOP "ENDED"
               ENDIF
               CALL qtip4pf_sr(atoms,nat,forces,pot,virial)
            ELSEIF (vstyle .ge. 60 .and. vstyle .le. 65 ) THEN 
               ! qtip4pf committee potential. adds two different types of (small)
               ! LJ potentials just to have variations on a theme

               IF (mod(nat,3)/=0) THEN
                  WRITE(*,*) " Expecting water molecules O H H O H H O H H but got ", nat, "atoms"
                  STOP "ENDED"
               ENDIF
               vpars(1) = cell_h(1,1)
               vpars(2) = cell_h(2,2)
               vpars(3) = cell_h(3,3)
               IF (cell_h(1,2).gt.1d-10 .or. cell_h(1,3).gt.1d-12  .or. cell_h(2,3).gt.1d-12) THEN
                  WRITE(*,*) " qtip4pf PES only works with orthorhombic cells", cell_h(1,2), cell_h(1,3), cell_h(2,3)
                  STOP "ENDED"
               ENDIF
               IF (vstyle == 63 .or. vstyle == 64 .or. vstyle == 65) THEN
                  pot = 0.0
                  forces = 0.0
                  virial = 0.0
               ELSE
                  CALL qtip4pf(vpars(1:3),atoms,nat,forces,pot,virial)
               ENDIF 

               ! adds a small LJ potential to give different committee values
               ! we have to replicate the code to make neighbor lists
               rc = 12.0d0  ! hardcoded cutoff
               IF ((allocated(n_list) .neqv. .true.)) THEN
                  IF (verbose > 0) WRITE(*,*) " Allocating neighbour lists. Cutoff ", rc
                  ALLOCATE(n_list(nat*(nat-1)/2))
                  ALLOCATE(index_list(nat))
                  ALLOCATE(last_atoms(nat,3))
                  last_atoms = 0.0d0
                  rn = rc*1.2
                  CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                  last_atoms = atoms
                  init_volume = volume
                  init_rc = rc
               ENDIF

               ! Checking to see if we need to re-calculate the neighbour list
               rc = init_rc*(volume/init_volume)**(1.0/3.0)
               DO i = 1, nat
                  CALL separation(cell_h, cell_ih, atoms(i,:), last_atoms(i,:), displacement)
                  ! Note that displacement is the square of the distance moved by atom i since the last time the neighbour list was created.
                  IF (4*displacement > (rn-rc)*(rn-rc)) THEN
                     IF (verbose > 0) WRITE(*,*) " Recalculating neighbour lists"
                     CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                     last_atoms = atoms
                     rn = 1.2*rc
                     EXIT
                  ENDIF
               ENDDO

               IF (vstyle == 60 .or. vstyle == 63) THEN ! type 1
                  CALL LJ_getall(rc, 2.5d0, 2d-6, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ELSEIF (vstyle == 61 .or. vstyle == 64) THEN ! type 2
                  CALL LJ_getall(rc, 2.1d0, 24d-6, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ELSEIF (vstyle == 62 .or. vstyle == 65) THEN ! returns both as json
                  ! return both the committee members as a JSON extra string
                  CALL LJ_getall(rc, 2.5d0, 2d-6, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
                  
                  string1=""
                  write(string,'(f15.8)') pot
                  string1 = TRIM(string) // ", "

                  string2="[ "
                  DO i=1,nat-1
                     WRITE(string,'(*(f15.8,","))') forces(i,:)
                     string2 = TRIM(string2) // TRIM(string)
                  END DO
                  write(string,'((f15.8,","), (f15.8,","), f15.8)') forces(nat,:)
                  string2 = TRIM(string2) // TRIM(string) // " ], "                  
               
                  WRITE(string3, '("[", 8(f15.8,",") f15.8, "]")') reshape(virial, (/9/))

                  ! this is a ugly but easy way to compute both terms
                  CALL LJ_getall(rc, 2.5d0, -2d-6, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
                  CALL LJ_getall(rc, 2.1d0, 24d-6, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)

                  write(string,'(f15.8)') pot
                  string1 = TRIM(string1) // TRIM(string)

                  string2 = TRIM(string2) // "[ "
                  DO i=1,nat-1
                     WRITE(string,'(*(f15.8,","))') forces(i,:)
                     string2 = TRIM(string2) // TRIM(string)
                  END DO
                  write(string,'((f15.8,","), (f15.8,","), f15.8)') forces(nat,:)
                  string2 = TRIM(string2) // TRIM(string) // " ] "

                  WRITE(string, '("[", 8(f15.8,",") f15.8, "]")') reshape(virial, (/9/))
                  string3 = TRIM(string3) // ", " // TRIM(string)

                  initbuffer = '{ "committee_pot" : [' 
                  initbuffer = TRIM(initbuffer) // trim(string1)
                  initbuffer = TRIM(initbuffer) // '],  "committee_force" : [ '
                  initbuffer = TRIM(initbuffer) // trim(string2) // ' ], '
                  initbuffer = TRIM(initbuffer) // '"committee_virial" : ['//  trim(string3) // ' ] '

                  initbuffer = TRIM(initbuffer) // '}'

                  ! and now make sure we are returning the ensemble mean
                  CALL LJ_getall(rc, 2.1d0, -12d-6, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
                  CALL LJ_getall(rc, 2.5d0, 1d-6, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ENDIF
            ELSEIF (70 == vstyle) THEN 
               ! potential that can be applied on top of water molecules to add a orientation-dependent term to 
               ! test non-equivariant terms in the potential
               vpars(1) = cell_h(1,1)
               vpars(2) = cell_h(2,2)
               vpars(3) = cell_h(3,3)
               
               pot = 0
               forces = 0
               DO i=1,nat,3
                  dip = 0.5 * (atoms(i+1,:) + atoms(i+2,:)) - atoms(i,:)
                  vpars(4) = sqrt(dip(1)**2+dip(2)**2+dip(3)**2)                  
                  vpars(1) = exp(dip(1)/vpars(4)*1.0)
                  vpars(2) = exp(dip(2)/vpars(4)*0.5)
                  vpars(3) = exp(dip(3)/vpars(4)*2.0)
                  pot = pot + 2e-2*(vpars(1)+vpars(2)+vpars(3))
                  vpars(1) = 2e-2*vpars(1)*1.0/vpars(4)**3
                  vpars(2) = 2e-2*vpars(2)*0.5/vpars(4)**3
                  vpars(3) = 2e-2*vpars(3)*2.0/vpars(4)**3
                  ! gradients on O
                  vpars(4) = -vpars(1)*(dip(2)**2+dip(3)**2) + vpars(2)*dip(1)*dip(2) + vpars(3)*dip(3)*dip(1) 
                  vpars(5) = -vpars(2)*(dip(1)**2+dip(3)**2) + vpars(1)*dip(1)*dip(2) + vpars(3)*dip(2)*dip(3) 
                  vpars(6) = -vpars(3)*(dip(1)**2+dip(2)**2) + vpars(1)*dip(1)*dip(3) + vpars(2)*dip(3)*dip(2) 
                  forces(i,1) = forces(i,1) - vpars(4)
                  forces(i,2) = forces(i,2) - vpars(5)
                  forces(i,3) = forces(i,3) - vpars(6)
                  ! gradients from H are just from Newton's law
                  forces(i+1,1) = forces(i+1,1) + vpars(4)*0.5
                  forces(i+1,2) = forces(i+1,2) + vpars(5)*0.5
                  forces(i+1,3) = forces(i+1,3) + vpars(6)*0.5
                  forces(i+2,1) = forces(i+2,1) + vpars(4)*0.5
                  forces(i+2,2) = forces(i+2,2) + vpars(5)*0.5
                  forces(i+2,3) = forces(i+2,3) + vpars(6)*0.5
               ENDDO
               
            ELSEIF (vstyle == 11) THEN ! efield potential.
               IF (mod(nat,3)/=0) THEN
                  WRITE(*,*) " Expecting water molecules O H H O H H O H H but got ", nat, "atoms"
                  STOP "ENDED"
               ENDIF
               CALL efield_v(atoms,nat,forces,pot,virial,efield)
            ELSEIF (vstyle == 8) THEN ! PS water potential.
               IF (nat/=3) THEN
                  WRITE(*,*) "Expecting 3 atoms for P-S water potential, O H H "
                  STOP "ENDED"
               ENDIF

               dip=0.0
               vecdiff=0.0
               ! lets fold the atom positions back to center in case the water travelled far away.
               ! this avoids problems if the water is splic across (virtual) periodic boundaries
               ! OH_1
               call vector_separation(cell_h, cell_ih, atoms(2,:), atoms(1,:), vecdiff, dist)
               atoms(2,:)=vecdiff(:)
               ! OH_2
               call vector_separation(cell_h, cell_ih, atoms(3,:), atoms(1,:), vecdiff, dist)
               atoms(3,:)=vecdiff(:)
               ! O in center
               atoms(1,:)=0.d0



               atoms = atoms*0.52917721d0    ! pot_nasa wants angstrom
               call pot_nasa(atoms, forces, pot)
               call dms_nasa(atoms, charges, dummy) ! MR: trying to print out the right charges
               dip(:)=atoms(1,:)*charges(1)+atoms(2,:)*charges(2)+atoms(3,:)*charges(3)
               ! MR: the above line looks like it provides correct results in eAngstrom for dipole!
               pot = pot*0.0015946679     ! pot_nasa gives kcal/mol
               forces = forces * (-0.00084329756) ! pot_nasa gives V in kcal/mol/angstrom
               ! do not compute the virial term
            ELSEIF (vstyle == 9) THEN
               IF (nat /= 3) THEN
                  WRITE(*,*) "Expecting 3 atoms for LEPS Model 1  potential, A B C "
                  STOP "ENDED"
               END IF
               CALL LEPS_M1(3, atoms, pot, forces)
            ELSEIF (vstyle == 10) THEN
               IF (nat /= 3) THEN
                  WRITE(*,*) "Expecting 4 atoms for LEPS Model 2  potential, A B C D n"
                  STOP "ENDED"
               END IF
               CALL LEPS_M2(4, atoms, pot, forces)

            ELSEIF (vstyle == 20) THEN ! eckart potential.
               CALL geteckart(nat,vpars(1), vpars(2), vpars(3),vpars(4), atoms, pot, forces)
            ELSEIF (vstyle == 28) THEN ! harmonic_bath.
               CALL get_harmonic_bath(nat,vpars(1),vpars(2),vpars(3),vpars(4),vpars(5),vpars(6),atoms, pot, forces)
            ELSEIF (vstyle == 29) THEN ! meanfield_bath.
               CALL get_meanfield_harmonic_bath(nat,vpars(1),vpars(2),vpars(3),vpars(4),atoms, pot, forces,friction)
            ELSEIF (vstyle == 23) THEN ! MB.
               IF (nat/=1) THEN
                  WRITE(*,*) "Expecting 1 atom for MB"
                  STOP "ENDED"
               ENDIF
               !atoms = atoms*0.52917721d0  !Change to angstrom
               CALL get_MB(nat,vpars(1), atoms, pot, forces)
            ELSEIF (vstyle == 25) THEN ! qQ
               CALL getdoublewell(nat, atoms, pot, forces)
               CALL dw_friction(nat, atoms, friction)

            ELSEIF (vstyle == 24) THEN ! qQ
               CALL getdoublewell_1D(nat, atoms, pot, forces)
               CALL dw1d_friction(nat, atoms, friction)
               CALL dw1d_dipole(nat, atoms, dip)

            ELSEIF (vstyle == 31) THEN   ! Sets force and potential to zero,
                                         ! computes only dipole moment, its gradient, and polarizability.
               pot = 0
               forces = 0.0d0
               DO i = 1, 3
                  vpars(i+1) = cell_h(i, i)
               ENDDO
               IF (.NOT. ALLOCATED(dipz_der)) ALLOCATE (dipz_der(nat, 3))
               CALL h2o_dipole(vpars(2:4), nat, atoms, vpars(1) /= 0, dip, dipz_der, pol)
            ELSE
               IF ((allocated(n_list) .neqv. .true.)) THEN
                  IF (verbose > 0) WRITE(*,*) " Allocating neighbour lists."
                  ALLOCATE(n_list(nat*(nat-1)/2))
                  ALLOCATE(index_list(nat))
                  ALLOCATE(last_atoms(nat,3))
                  last_atoms = 0.0d0
                  CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                  last_atoms = atoms
                  init_volume = volume
                  init_rc = rc
               ENDIF

               ! Checking to see if we need to re-calculate the neighbour list
               rc = init_rc*(volume/init_volume)**(1.0/3.0)
               DO i = 1, nat
                  CALL separation(cell_h, cell_ih, atoms(i,:), last_atoms(i,:), displacement)
                  ! Note that displacement is the square of the distance moved by atom i since the last time the neighbour list was created.
                  IF (4*displacement > (rn-rc)*(rn-rc)) THEN
                     IF (verbose > 0) WRITE(*,*) " Recalculating neighbour lists"
                     CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                     last_atoms = atoms
                     rn = 1.2*rc
                     EXIT
                  ENDIF
               ENDDO

               ! zeroes out forces and al the rest
               forces = 0.0d0
               pot = 0.0d0
               virial = 0.0d0               
               IF (vstyle == 1) THEN
                  CALL LJ_getall(rc, sigma, eps, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ELSEIF (vstyle == 2) THEN
                  CALL SG_getall(rc, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ELSEIF (vstyle == 22) THEN ! ljpolymer potential.
                  CALL ljpolymer_getall(n_monomer,rc,sigma,eps,stiffness,nat,atoms,cell_h,cell_ih,index_list,n_list,pot,forces,virial)
               ELSEIF (vstyle == 32) THEN ! lj mixture.
                  CALL ljmix_getall(n_monomer,rc,sigma,eps,nat,atoms,cell_h,cell_ih,index_list,n_list,pot,forces,virial)
               ENDIF
               IF (verbose > 0) WRITE(*,*) " Calculated energy is ", pot
            ENDIF
            hasdata = .true. ! Signal that we have data ready to be passed back to the wrapper
         ELSEIF (trim(header) == "GETFORCE") THEN  ! The driver calculation is finished, it's time to send the results back to the wrapper

            ! Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
            DO i = 1, nat
               msgbuffer(3*(i-1)+1:3*i) = forces(i,:)
            ENDDO
            virial = transpose(virial)

            CALL writebuffer(socket,"FORCEREADY  ",MSGLEN)
            IF (verbose > 1) WRITE(*,*) "    !write!=> ", "FORCEREADY  "
            CALL writebuffer(socket,pot)  ! Writing the potential
            IF (verbose > 1) WRITE(*,*) "    !write!=> pot: ", pot
            CALL writebuffer(socket,nat)  ! Writing the number of atoms
            IF (verbose > 1) WRITE(*,*) "    !write!=> nat:", nat
            CALL writebuffer(socket,msgbuffer,3*nat) ! Writing the forces
            IF (verbose > 1) WRITE(*,*) "    !write!=> forces:", msgbuffer(0:2), " ..."
            CALL writebuffer(socket,reshape(virial,(/9/)),9)  ! Writing the virial tensor, NOT divided by the volume
            IF (verbose > 1) WRITE(*,*) "    !write!=> strss: ", reshape(virial,(/9/))

 125  format(es21.14,a,es21.14,a,es21.14,a,es21.14,a,es21.14,a,es21.14,a)
 126  format(es21.14,a,es21.14,a,es21.14,a,es21.14,a,es21.14,a)

            IF (vstyle == 29) THEN ! returns meanfield friction
                WRITE(initbuffer,'(a)') "{"
                WRITE(32,'(a)') '{'
                WRITE(string,'(a)') '"friction": ['
                WRITE(32,'(a)') '"friction": ['

                string2 = TRIM(initbuffer) // TRIM(string)
                initbuffer = TRIM(string2)
                DO i=1,3*nat
                    IF(i/=3*nat) THEN
                        WRITE(string,125) ( friction(i,j), "," , j=1,3*nat)
                        WRITE(32,125) ( friction(i,j), "," , j=1,3*nat)
                    ELSE
                        WRITE(string,126) ( friction(i,j), "," , j=1,3*nat-1)
                        WRITE(string2,'(es21.14)') friction(i,3*nat)
                        string3 = TRIM(string) // TRIM(string2)
                        string = string3
                        WRITE(32,126) ( friction(i,j), "," , j=1,3*nat-1)
                        WRITE(32,'(es21.14)') friction(i,3*nat)
                    ENDIF
                    string2 = TRIM(initbuffer) // TRIM(string)
                    initbuffer = TRIM(string2)
                END DO
                string =  TRIM(initbuffer) // ']}'
                initbuffer = TRIM(string)
                WRITE(32,'(a)') "]"
                WRITE(32,'(a)') "}"

                cbuf = LEN_TRIM(initbuffer)
                CALL writebuffer(socket,cbuf)

                IF (verbose > 1) WRITE(*,*) "!write!=> extra_length:", cbuf
                CALL writebuffer(socket,initbuffer,cbuf)
                IF (verbose > 1) WRITE(*,*) "    !write!=> extra: ",  initbuffer

            ELSEIF (vstyle==24 .or. vstyle==25) THEN ! returns fantasy friction
                WRITE(initbuffer,'(a)') "{"
                WRITE(string, '(a,3x,f15.8,a,f15.8,a,f15.8,&
     &          3x,a)') '"dipole": [',dip(1),",",dip(2),",",dip(3),"],"
                string2 = TRIM(initbuffer) // TRIM(string)
                initbuffer = TRIM(string2)

                WRITE(string,'(a)') '"friction": ['
                string2 = TRIM(initbuffer) // TRIM(string)
                initbuffer = TRIM(string2)
                DO i=1,3*nat
                    WRITE(string,'("[ ",*(f15.8,","))') friction(i,:)
                    length = LEN_TRIM(string)
                    trimmed = TRIM(string)
                    IF(i==3*nat) THEN
                        string = TRIM(trimmed(:length-1)) // "]"
                    ELSE
                        string = TRIM(trimmed(:length-1)) // "],"
                    ENDIF
                    string2 = TRIM(initbuffer) // TRIM(string)
                    initbuffer = TRIM(string2)
                END DO
                string =  TRIM(initbuffer) // ']}'
                initbuffer = TRIM(string)
                cbuf = LEN_TRIM(initbuffer)
                CALL writebuffer(socket,cbuf) ! Writes back the fantasy friction
                IF (verbose > 1) WRITE(*,*) "!write!=> extra_length:", &
     &          cbuf
                CALL writebuffer(socket,initbuffer,cbuf)
                IF (verbose > 1) WRITE(*,*) "    !write!=> extra: ",  &
     &          initbuffer(1:cbuf)
            ELSEIF (vstyle==5 .or. vstyle==6 .or. vstyle==8 .or. vstyle==99) THEN ! returns the  dipole through initbuffer
               WRITE(initbuffer, '(a,3x,f15.8,a,f15.8,a,f15.8, &
     &         3x,a)') '{"dipole": [',dip(1),",",dip(2),",",dip(3),"]}"
               cbuf = LEN_TRIM(initbuffer)
               CALL writebuffer(socket,cbuf) ! Writes back the molecular dipole
               IF (verbose > 1) WRITE(*,*)  &
     &         "    !write!=> extra_length: ", cbuf
               CALL writebuffer(socket,initbuffer,cbuf)
               IF (verbose > 1) WRITE(*,*) "    !write!=> extra: ", &
     &         initbuffer(1:cbuf)               
               
            ELSEIF (vstyle==31) THEN ! returns the dipole, dipole derivative, and polarizability through initbuffer
               WRITE(string, '(a,3x,f15.8,a,f15.8,a,f15.8, 3x,a)') '{"dipole": [',dip(1),",",dip(2),",",dip(3),"],"
               longbuffer = TRIM(string)
               WRITE(string2, *) "(a,3x,", 3*nat - 1, '(f15.8, ","),f15.8,3x,a)'
               WRITE(longstring, string2) '"dipole_derivative": [',TRANSPOSE(dipz_der),"],"
               longbuffer = TRIM(longbuffer)//TRIM(longstring)
               WRITE(string, '(a,3x, 8(f15.8, ","),f15.8,3x,a)') '"polarizability": [',pol,"]}"
               longbuffer = TRIM(longbuffer)//TRIM(string)
               cbuf = LEN_TRIM(longbuffer)
               CALL writebuffer(socket,cbuf)
               CALL writebuffer(socket,TRIM(longbuffer),cbuf)
               IF (verbose > 1) WRITE(*,*) "    !write!=> extra: ", &               
     &         initbuffer
            ELSEIF (vstyle==62 .or. vstyle==65) THEN ! returns committee data
               cbuf = LEN_TRIM(initbuffer)
               CALL writebuffer(socket,cbuf)
               CALL writebuffer(socket,initbuffer,cbuf)
            ELSE
               cbuf = 1 ! Size of the "extras" string
               CALL writebuffer(socket,cbuf) ! This would write out the "extras" string, but in this case we only use a dummy string.
               IF (verbose > 1) WRITE(*,*)  &
     &         "    !write!=> extra_length: ", cbuf
               CALL writebuffer(socket,' ',1)
               IF (verbose > 1) WRITE(*,*)  &
     &         "    !write!=> extra: empty"
            ENDIF
            hasdata = .false.
         ELSEIF (trim(header) == "EXIT") THEN
            EXIT
         ELSE
            WRITE(*,*) " Unexpected header ", header
            STOP "ENDED"
         ENDIF
      ENDDO
      IF (nat > 0) DEALLOCATE(atoms, forces, msgbuffer)
      IF (nat>0 .and. (vstyle==24 .or. vstyle==25)) THEN
         DEALLOCATE(friction)
      ENDIF
      STOP

    CONTAINS
      SUBROUTINE helpmessage
         ! Help banner

         WRITE(*,*) " SYNTAX: driver.x [-u] -a address [-p port] -m [dummy|gas|lj|sg|harm|harm3d|morse|morsedia|zundel|qtip4pf|pswater|lepsm1|lepsm2|qtip4p-efield|eckart|ch4hcbe|ljpolymer|MB|doublewell|doublewell_1D|water_dip_pol|harmonic_bath|meanfield_bath|ljmix|qtip4pf-sr|qtip4pf-c-1|qtip4pf-c-2|qtip4pf-c-json|qtip4pf-c-1-delta|qtip4pf-c-2-delta|qtip4pf-c-json-delta|noo3-h2o]"
         WRITE(*,*) "         -o 'comma_separated_parameters' [-S sockets_prefix] [-v] "
         WRITE(*,*) ""
         WRITE(*,*) " For LJ potential use -o sigma,epsilon,cutoff "
         WRITE(*,*) " For SG potential use -o cutoff "
         WRITE(*,*) " For 1D/3D harmonic oscillator use -o k "
         WRITE(*,*) " For 1D morse oscillators use -o r0,D,a"
         WRITE(*,*) " For qtip4pf-efield use -o Ex,Ey,Ez with Ei in V/nm"
         WRITE(*,*) " For ljpolymer or lkmix use -o n_monomer,sigma,epsilon,cutoff "
         WRITE(*,*) " For gas, dummy, use the optional -o sleep_seconds to add a delay"
         WRITE(*,*) " For the ideal, qtip4pf*, zundel, ch4hcbe, nasa, doublewell or doublewell_1D no options are needed! "
       END SUBROUTINE helpmessage

   END PROGRAM
