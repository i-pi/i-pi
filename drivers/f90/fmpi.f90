! MPI transport for the i-PI ffmpi forcefield (compiled-driver pilot).
!
! Copyright (C) 2024, i-PI developers
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

      MODULE FMPI
         USE MPI
         IMPLICIT NONE

         ! message tags (must match forcefields.py and driver.py)
         INTEGER, PARAMETER :: TAG_HELLO = 1, TAG_INIT = 2, TAG_WORK = 3
         INTEGER, PARAMETER :: TAG_RESULT = 4, TAG_EXTRA = 5, TAG_EXIT = 6
         ! intra-group broadcast control flags (driver-internal)
         INTEGER, PARAMETER :: GROUP_WORK = 0, GROUP_EXIT = 1

      CONTAINS

         ! Announces this group to i-PI: the root sends its world rank, the
         ! group size and the world ranks of all its members.
         SUBROUTINE mpi_send_hello(world, server, world_rank, group, ff_name)
            INTEGER, INTENT(IN) :: world, server, world_rank, group
            CHARACTER(LEN=*), INTENT(IN) :: ff_name
            INTEGER :: gsize, grank, ierr
            INTEGER, ALLOCATABLE :: members(:), hello(:)

            CALL MPI_Comm_size(group, gsize, ierr)
            CALL MPI_Comm_rank(group, grank, ierr)
            ALLOCATE(members(gsize))
            CALL MPI_Allgather(world_rank, 1, MPI_INTEGER, members, 1, &
                               MPI_INTEGER, group, ierr)
            IF (grank == 0) THEN
               ALLOCATE(hello(2 + gsize))
               hello(1) = world_rank
               hello(2) = gsize
               hello(3:) = members
               CALL MPI_Send(hello, 2 + gsize, MPI_INTEGER, server, &
                             TAG_HELLO, world, ierr)
               DEALLOCATE(hello)
               ! the forcefield name follows as a byte message (empty = untagged)
               CALL MPI_Send(ff_name, LEN_TRIM(ff_name), MPI_BYTE, server, &
                             TAG_HELLO, world, ierr)
            END IF
            DEALLOCATE(members)
         END SUBROUTINE

         ! Sends the packed numeric results, then the per-structure 'extra'
         ! strings as one byte message: [int32 length][bytes] per structure.
         SUBROUTINE mpi_send_results(world, server, nat, batch, numeric, &
                                     extra_msg, nextra)
            INTEGER, INTENT(IN) :: world, server, nat, batch, nextra
            DOUBLE PRECISION, INTENT(IN) :: numeric(:)
            CHARACTER(LEN=*), INTENT(IN) :: extra_msg
            INTEGER :: ierr

            CALL MPI_Send(numeric, batch * (1 + 3 * nat + 9), &
                          MPI_DOUBLE_PRECISION, server, TAG_RESULT, world, ierr)
            CALL MPI_Send(extra_msg, nextra, MPI_BYTE, server, TAG_EXTRA, &
                          world, ierr)
         END SUBROUTINE

      END MODULE FMPI
