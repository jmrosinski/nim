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


subroutine taskinfo_comm(compute_tasks, cpu_cores_per_node,             &
                         max_compute_tasks_per_node,                    &
                         omp_threads_per_compute_task, num_write_tasks, &
                         max_write_tasks_per_node, root_own_node,       &
                         icosio_debugmsg_on, max_compute_tasks_per_mic, &
                         max_accelerators_per_node,                     &
                         omp_threads_per_mic_mpi_task, comm)

! If this code is executed in an MPI program, package task info in a
! user-defined type, and broadcast from the root to the other tasks using 
! MPI communicator "comm".  Otherwise do nothing.  
! See taskinfo.F90 for an explanation of the other arguments.   
!
! This is an MPI program iff cpp token "SERIAL" is undefined.  
!
! Note that direct MPI calls must be used here because 
! this routine may be called before SMS setup is complete.  

#ifndef SERIAL
  use mpi
#endif

  implicit none

  ! these INOUT arguments are intent(in) on the root and intent(out) elsewhere
  integer, intent(inout) :: compute_tasks
  integer, intent(inout) :: cpu_cores_per_node
  integer, intent(inout) :: max_compute_tasks_per_node
  integer, intent(inout) :: omp_threads_per_compute_task
  integer, intent(inout) :: num_write_tasks
  integer, intent(inout) :: max_write_tasks_per_node
  integer, intent(inout) :: max_compute_tasks_per_mic
  integer, intent(inout) :: max_accelerators_per_node
  integer, intent(inout) :: omp_threads_per_mic_mpi_task

  logical, intent(inout) :: root_own_node
  logical, intent(inout) :: icosio_debugmsg_on

  integer, intent(in   ) :: comm                 ! an MPI intracommunicator

  integer :: ierr           ! mpi return code
  integer :: me             ! my rank
  integer :: iarr(9)        ! array holds integers to be broadcast (the integer inout args above)
  logical :: larr(2)        ! array holds logicals to be broadcast (the logical inout args above)

!SMS$IGNORE BEGIN

#ifndef SERIAL

  call mpi_comm_rank(comm, me, ierr)
  if (ierr.ne.0) then
    write (*,'(a,i0)') 'task_info: MPI_Comm_Rank returned ',ierr
    call flush(6)
    call mpi_abort (comm, 999, ierr)
    stop 999 ! In case mpi_abort fails
  endif

  if (me.eq.0) then

! Populate integer and logical arrays with values to be broadcast

    iarr(1) = compute_tasks
    iarr(2) = cpu_cores_per_node
    iarr(3) = max_compute_tasks_per_node
    iarr(4) = omp_threads_per_compute_task
    iarr(5) = num_write_tasks
    iarr(6) = max_write_tasks_per_node
    iarr(7) = max_compute_tasks_per_mic
    iarr(8) = omp_threads_per_mic_mpi_task
    iarr(9) = max_accelerators_per_node

    larr(1) = root_own_node
    larr(2) = icosio_debugmsg_on
  endif

! Broadcast the integer task configuration info

  call mpi_bcast (iarr, size(iarr), MPI_INTEGER, 0, comm, ierr)
  if (ierr.ne.0) then
    write (*,'(a,i0)') 'taskinfo_comm bcast of integers: MPI_Bcast returned ',ierr
    call flush(6)
    call mpi_abort (comm, 999, ierr)
    stop 999 ! In case mpi_abort fails
  endif

! Broadcast the logical task configuration info

  call mpi_bcast (larr, size(larr), MPI_LOGICAL, 0, comm, ierr)
  if (ierr.ne.0) then
    write (*,'(a,i0)') 'taskinfo_comm bcast of logicals: MPI_Bcast returned ',ierr
    call flush(6)
    call mpi_abort (comm, 999, ierr)
    stop 999 ! In case mpi_abort fails
  endif

! Unpack conf to individual values (only rquired on non-root tasks)

  compute_tasks                = iarr(1)
  cpu_cores_per_node           = iarr(2)
  max_compute_tasks_per_node   = iarr(3)
  omp_threads_per_compute_task = iarr(4)
  num_write_tasks              = iarr(5)
  max_write_tasks_per_node     = iarr(6)
  max_compute_tasks_per_mic    = iarr(7)
  omp_threads_per_mic_mpi_task = iarr(8)
  max_accelerators_per_node    = iarr(9)

  root_own_node                = larr(1)
  icosio_debugmsg_on           = larr(2)
#endif

  return

!SMS$IGNORE END

end subroutine taskinfo_comm
