!F90 ISO_C_BINDING wrapper for socket communication.

!Copyright (C) 2013, Michele Ceriotti

!Permission is hereby granted, free of charge, to any person obtaining
!a copy of this software and associated documentation files (the
!"Software"), to deal in the Software without restriction, including
!without limitation the rights to use, copy, modify, merge, publish,
!distribute, sublicense, and/or sell copies of the Software, and to
!permit persons to whom the Software is furnished to do so, subject to
!the following conditions:

!The above copyright notice and this permission notice shall be included
!in all copies or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!Contains both the functions that transmit data to the socket and read the data
!back out again once finished, and the function which opens the socket initially.

!Functions:
!   open_socket: Opens a socket with the required host server, socket type and
!      port number.
!   write_buffer: Writes a string to the socket.
!   read_buffer: Reads data from the socket.

   MODULE F90SOCKETS
   USE ISO_C_BINDING

   IMPLICIT NONE

   INTEGER,PARAMETER  :: AF_INET = 0, AF_UNIX = 1
   INTEGER,PARAMETER  :: AI_PASSIVE = 1
   INTEGER,PARAMETER  :: SOCK_STREAM = 1
   INTEGER,PARAMETER  :: UNIX_PATH_MAX = 108

   TYPE, BIND(C) :: sockaddr_un
      INTEGER(KIND=C_SHORT) sun_family
      CHARACTER(LEN=1, KIND=C_CHAR) sun_path(UNIX_PATH_MAX)
   END TYPE

   TYPE, BIND(C) :: addrinfo
      INTEGER(KIND=C_INT) :: ai_flags, ai_family, ai_socktype, &
          ai_protocol
      INTEGER(KIND=C_SIZE_T) :: ai_addrlen
      TYPE(C_PTR) :: ai_addr, ai_canonname, ai_next
   END TYPE

   INTERFACE
      FUNCTION getaddrinfo(node, service, hints, res) BIND(C)
         USE ISO_C_BINDING
         TYPE(C_PTR), INTENT(IN), VALUE :: node, service, hints
         TYPE(C_PTR) :: RES
         INTEGER(KIND=C_INT) :: getaddrinfo
      END FUNCTION
   END INTERFACE

   INTERFACE
      SUBROUTINE freeaddrinfo( res) BIND(C)
         USE ISO_C_BINDING
         TYPE(C_PTR), INTENT(IN), VALUE :: RES
      END SUBROUTINE
   END INTERFACE

   INTERFACE
     FUNCTION socket_make(af, type, protocol) BIND(C, name="socket")
      USE ISO_C_BINDING
         INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: af, type, protocol
         INTEGER(KIND=C_INT)  :: socket_make
     END FUNCTION
   END INTERFACE

   INTERFACE
     FUNCTION socket_connect(sockfd, addr, addrlen) BIND(C, name="connect")
      USE ISO_C_BINDING
         INTEGER(KIND=C_INT), VALUE  :: sockfd
         TYPE(C_PTR), INTENT(IN), VALUE :: addr
         INTEGER(KIND=C_SIZE_T), INTENT(IN), VALUE  :: addrlen
         INTEGER(KIND=C_INT)  :: socket_connect
     END FUNCTION
   END INTERFACE

   INTERFACE
      FUNCTION socket_write(sid, data, dcount) BIND(C, name="write")
         USE ISO_C_BINDING
            INTEGER(KIND=C_INT), INTENT(IN), VALUE     :: sid
            TYPE(C_PTR), INTENT(IN), VALUE     :: data
            INTEGER(KIND=C_SIZE_T), INTENT(IN), VALUE  :: dcount
            INTEGER(KIND=C_SIZE_T)     :: socket_write
      END FUNCTION
   END INTERFACE

   INTERFACE 
      FUNCTION socket_read(sid, data, dcount) BIND(C, name="read")
         USE ISO_C_BINDING
            INTEGER(KIND=C_INT), INTENT(IN), VALUE     :: sid
            TYPE(C_PTR), INTENT(IN), VALUE     :: data
            INTEGER(KIND=C_SIZE_T), INTENT(IN), VALUE  :: dcount
            INTEGER(KIND=C_SIZE_T)     :: socket_read
      END FUNCTION
   END INTERFACE

   INTERFACE
      SUBROUTINE memcpy(dout, din, dcount) BIND(C, name="memcpy")
         USE ISO_C_BINDING
            TYPE(C_PTR), INTENT(IN), VALUE     :: din, dout
            INTEGER(KIND=C_SIZE_T), INTENT(IN), VALUE  :: dcount
      END SUBROUTINE
   END INTERFACE

   INTERFACE
      TYPE(C_PTR) FUNCTION memset(s, c, n) BIND(C, name="memset")
         USE ISO_C_BINDING
            TYPE(C_PTR), INTENT(IN), VALUE     :: s
            INTEGER(KIND=C_INT), INTENT(IN), VALUE :: c
            INTEGER(KIND=C_SIZE_T), INTENT(IN), VALUE  :: n
      END FUNCTION
   END INTERFACE


   ! writebuffer interfaces
   INTERFACE writebuffer
      MODULE PROCEDURE writebuffer_s, &
                       writebuffer_d, writebuffer_dv, &
                       writebuffer_i

   END INTERFACE 

   ! readbuffer interfaces
   INTERFACE readbuffer
      MODULE PROCEDURE readbuffer_s, &
                       readbuffer_dv, readbuffer_d, &
                       readbuffer_i

   END INTERFACE 
   CONTAINS

   SUBROUTINE fstr2cstr(fstr, cstr, plen)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fstr
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: cstr(:)
      INTEGER, INTENT(IN), OPTIONAL :: plen

      INTEGER i,n
      IF (PRESENT(plen)) THEN
         n = plen
         DO i=1,n
            cstr(i) = fstr(i:i)
         ENDDO
      ELSE
         n = LEN_TRIM(fstr)
         DO i=1,n
            cstr(i) = fstr(i:i)
         ENDDO
         cstr(n+1) = C_NULL_CHAR
      END IF
   END SUBROUTINE

   SUBROUTINE open_socket(psockfd, inet, port, host)      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: inet, port
      INTEGER, INTENT(OUT) :: psockfd
      CHARACTER(LEN=1024), INTENT(IN) :: host

      INTEGER :: ai_err
      CHARACTER(LEN=256) :: service
      CHARACTER(LEN=1,KIND=C_CHAR), TARGET :: cservice(256)
      CHARACTER(LEN=1,KIND=C_CHAR), TARGET :: chost(1024)

      TYPE(addrinfo), TARGET :: hints
      TYPE(addrinfo), TARGET :: res
      TYPE(sockaddr_un), TARGET :: addrun
      TYPE(C_PTR) :: ptr

      CALL fstr2cstr(host, chost)
      IF (INET>0) THEN
      ! creates an internet socket

      ! fetches information on the host

         ptr = memset(c_loc(hints), 0, SIZEOF(hints))
         hints%ai_socktype = SOCK_STREAM
         hints%ai_family = AF_INET
         hints%ai_flags = AI_PASSIVE

         WRITE(service,'(I10)') port
         service=ADJUSTL(TRIM(service))
         CALL fstr2cstr(service,cservice)

         ai_err = getaddrinfo(c_loc(chost(1)), c_loc(cservice(1)), &
                    c_loc(hints), ptr)
         IF (ai_err < 0) THEN
            WRITE(6,*) "Error fetching host data. Wrong host name?"
            STOP " ENDED "
         ENDIF

         CALL memcpy(c_loc(res), ptr, sizeof(res))
         WRITE(6,*) "pointer", res
         psockfd = socket_make(res%ai_family, res%ai_socktype, res%ai_protocol)
         IF (psockfd < 0)  THEN 
            WRITE(6,*) "Error opening socket"
            STOP " ENDED "
         ENDIF

         ai_err = socket_connect(psockfd, res%ai_addr, res%ai_addrlen)
         IF (ai_err < 0) THEN
            WRITE(6,*) "Error opening INET socket: wrong port or server unreachable"
            STOP " ENDED "
         ENDIF

         CALL freeaddrinfo(ptr)
      ELSE
         ! creates an unix socket
         ptr = memset(c_loc(addrun), 0, SIZEOF(addrun))

         addrun%sun_family = AF_UNIX
         CALL fstr2cstr("/tmp/ipi_"//host, addrun%sun_path) 

         psockfd = socket_make(AF_UNIX, SOCK_STREAM, 0)

         ai_err = socket_connect(psockfd, c_loc(addrun), sizeof(addrun))
         IF (ai_err < 0) THEN
            WRITE(6,*) "Could not open UNIX socket. Non-existing path?"
            STOP " ENDED "
         ENDIF
      END IF
   END SUBROUTINE


   ! Set of wrappers to socket_write. Write data to a socket.
   ! Args:
   !  psockfd: The id of the socket that will be written to.
   !   data: The data to be written to the socket (different kinds)
   !   plen: The length of the data (not for single element writes)

   SUBROUTINE writebuffer_s(psockfd, fstring, plen)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: psockfd
      INTEGER, INTENT(IN) ::  plen
      CHARACTER(LEN=*), TARGET, INTENT(IN)  :: fstring
      INTEGER(8) :: nwrite, nlen
      CHARACTER(LEN=1,KIND=C_CHAR), TARGET :: cstring(plen)
      INTEGER i,n

      n = plen
      DO i = 1,n
         cstring(i) = fstring(i:i)
      ENDDO
      nlen = plen
      nwrite = socket_write(psockfd, c_loc(cstring(1)), nlen)
      IF (nwrite/=nlen) THEN
         WRITE(6,*) "Error in writing to socket buffer"
         STOP " ENDED "
      ENDIF
   END SUBROUTINE

   SUBROUTINE writebuffer_d(psockfd, fdata)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: psockfd
      DOUBLE PRECISION, TARGET, INTENT(IN)  :: fdata
      INTEGER(8) :: nwrite, nlen

      nlen = 8
      nwrite = socket_write(psockfd, c_loc(fdata), nlen)
      IF (nwrite/=nlen) THEN
         WRITE(6,*) "Error in writing to socket buffer"
         STOP " ENDED "
      ENDIF
   END SUBROUTINE

   SUBROUTINE writebuffer_i(psockfd, fdata)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: psockfd
      INTEGER, TARGET, INTENT(IN)  :: fdata
      INTEGER(8) :: nwrite, nlen

      nlen = 4
      nwrite = socket_write(psockfd, c_loc(fdata), nlen)
      IF (nwrite/=nlen) THEN
         WRITE(6,*) "Error in writing to socket buffer"
         STOP " ENDED "
      ENDIF
   END SUBROUTINE

   SUBROUTINE writebuffer_dv(psockfd, fdata, plen)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: psockfd, plen
      DOUBLE PRECISION, TARGET, INTENT(IN)  :: fdata(plen)
      INTEGER(8) :: nwrite, nlen

      nlen = 8*plen
      nwrite = socket_write(psockfd, c_loc(fdata(1)), nlen)
      IF (nwrite/=nlen) THEN
         WRITE(6,*) "Error in writing to socket buffer"
         STOP " ENDED "
      ENDIF
   END SUBROUTINE


   ! Set of wrappers to socket_read. Read data from a socket.
   ! Args:
   !  psockfd: The id of the socket we will read from.
   !   data: The data to be read from the socket (different kinds)
   !   plen: The length of the data (not for single element reads)
   ! NB: we always need to read to a c_str, since socket_read can return
   ! before the read has completed. Then, the c_str data must be copied to
   ! the output variable

   SUBROUTINE readbuffer_cstr(psockfd, cstring, plen)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: psockfd
      INTEGER, INTENT(IN) ::  plen
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT), TARGET :: cstring(plen)
      INTEGER(8) :: nread, nlen, n

      nlen = plen

      nread = socket_read(psockfd, c_loc(cstring(1)), nlen)
      n = nread
      DO WHILE(nread>0 .AND. n<nlen)
         nread = socket_read(psockfd, c_loc(cstring(n+1)), nlen-n)
         n = n + nread
      ENDDO

      IF (n<nlen) THEN
         WRITE(6,*) "Error in reading from socket"
         STOP " ENDED "
      ENDIF               
   END SUBROUTINE

   SUBROUTINE readbuffer_s(psockfd, fstring, plen)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: psockfd
      INTEGER, INTENT(IN) ::  plen
      CHARACTER(LEN=*), TARGET, INTENT(OUT)  :: fstring
      INTEGER(8) :: n, i
      CHARACTER(LEN=1,KIND=C_CHAR), TARGET :: cstring(plen)

      CALL readbuffer_cstr(psockfd, cstring, plen)

      n = plen
      DO i = 1,n
         fstring(i:i) = cstring(i)
      ENDDO
   END SUBROUTINE

   SUBROUTINE readbuffer_dv(psockfd, fdata, plen)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: psockfd, plen
      DOUBLE PRECISION, TARGET, INTENT(OUT)  :: fdata(plen)
      INTEGER(8) :: n
      CHARACTER(LEN=1,KIND=C_CHAR), TARGET :: cstring(plen*8)

      CALL readbuffer_cstr(psockfd, cstring, plen*8)
      n = plen*8
      CALL memcpy(c_loc(fdata(1)), c_loc(cstring(1)), n)
   END SUBROUTINE

   SUBROUTINE readbuffer_d(psockfd, fdata)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: psockfd
      DOUBLE PRECISION, TARGET, INTENT(OUT)  :: fdata
      INTEGER(8) :: n
      CHARACTER(LEN=1,KIND=C_CHAR), TARGET :: cstring(8)

      CALL readbuffer_cstr(psockfd, cstring, 8)
      n = 8
      CALL memcpy(c_loc(fdata), c_loc(cstring(1)), n)
   END SUBROUTINE

   SUBROUTINE readbuffer_i(psockfd, fdata)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: psockfd
      INTEGER, TARGET, INTENT(OUT)  :: fdata
      INTEGER(8) :: n
      CHARACTER(LEN=1,KIND=C_CHAR), TARGET :: cstring(4)

      CALL readbuffer_cstr(psockfd, cstring, 4)
      n = 4
      CALL memcpy(c_loc(fdata), c_loc(cstring(1)), n)
   END SUBROUTINE
   END MODULE
