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

program get_task_info

  use taskinfo, only:task_info

  implicit none

  integer, parameter :: unitno = 11        ! unit number for namelist read
  logical :: i_read
  logical :: write_yaml
  integer :: compute_tasks                 ! number of compute tasks
  integer :: cpu_cores_per_node            ! number of cores per node
  integer :: max_compute_tasks_per_node    ! number of MPI tasks per node
  integer :: omp_threads_per_compute_task  ! number of OpenMP threads
  integer :: num_write_tasks               ! number of write tasks
  integer :: max_write_tasks_per_node      ! max write tasks per node
  integer :: max_compute_tasks_per_mic
  integer :: max_accelerators_per_node
  integer :: omp_threads_per_mic_mpi_task
  logical :: root_own_node                 ! .false. if compute root shares node
  logical :: icosio_debugmsg_on            ! ignored in this program
  integer :: total_cpu_cores               ! total number of CPU cores
  integer :: total_mpi_tasks               ! total number of MPI tasks
  integer :: total_nodes                   ! total number of nodes
  integer :: mpi_tasks_per_node            ! number of MPI tasks per node
  integer :: write_nodes                   ! number of write nodes
  integer :: do_nothing_tasks              ! number of do-nothing tasks
  integer :: ret

  i_read = .true.
  write_yaml = .true.
  ! TODO:  Wrap this call in a simpler API since INTENT(OUT) values other 
  ! TODO:  than ret are never used here.  
  call task_info(unitno, i_read, write_yaml, compute_tasks,                   &
                 cpu_cores_per_node, max_compute_tasks_per_node,              &
                 omp_threads_per_compute_task, num_write_tasks,               &
                 max_write_tasks_per_node, root_own_node, icosio_debugmsg_on, &
                 max_compute_tasks_per_mic, omp_threads_per_mic_mpi_task,     &
                 total_cpu_cores, total_mpi_tasks, total_nodes,               &
                 mpi_tasks_per_node, write_nodes, do_nothing_tasks,           &
                 max_accelerators_per_node, ret)
  if (ret/=0) then
    write(6,*) 'get_task_info: task_info failed, ret=', ret
    stop 999
  endif
end program get_task_info
