! Copyright 2016
! United States Government as represented by the National Oceanic and
! Atmospheric Administration.
! No copyright is claimed in the United States under Title 17, U.S.Code.
! All Other Rights Reserved.
!
! Copyright 2016 Regents of the University of Colorado. All rights
! reserved.
! Copyright 2016 Regents of the Colorado State University. All rights
! reserved.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! globals: Contains global, static data. 
!          Perhaps things like ims, ime, ips, ipe, ihe could be added here to avoid
!          giant argument lists?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module globals
#ifndef SERIAL
  use mpi
#endif
  implicit none

  integer :: myRank = 0    ! MPI rank (always zero in serial mode)
  integer :: nPEs = 1      ! number of MPI tasks (always one in serial mode)
  logical :: GPUrun
  logical :: parallelBuild = .false.
#ifdef SERIAL
  character(len=6), parameter :: mynode = 'serial'
#else
  character(len=MPI_MAX_PROCESSOR_NAME) :: mynode = 'not_set_yet'
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! init_globals: Set global data for convenient use anywhere in NIM.
!               In serial mode this routine currently does nothing, since mynode is
!               a parameter already set to 'serial'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_globals ()
    integer :: resultlen   ! length returned from mpi_get_processor_name
    integer :: ierr        ! error return from MPI call

#ifndef SERIAL
    parallelBuild=.true.
!sms$comm_size(npes)
!sms$comm_rank(myrank)
    call mpi_get_processor_name (mynode, resultlen, ierr)
    if (ierr /= 0) then
      write(6,*)'init_globals: failure to get processor name'
      return
    end if
#endif

    return
  end subroutine init_globals
end module globals
