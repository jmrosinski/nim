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

module core_setup
  use taskinfo, only: task_info
  use units, only: getunit, returnunit
  use globals,      only: init_globals
  use print_taskinfo,only: print_nodeinfo
#ifndef SERIAL
  use mpi
#endif

  implicit none
  save
  private

#include <gptl.inc>

  integer, parameter :: badret = 1 ! Return code for MPI_Abort to pass to environment

! Public functions

  public :: core_setup_nim, iam_compute_root, iam_write_root
#ifdef SERIAL
  integer, parameter :: MPI_COMM_NULL=0
  integer, parameter :: mpi_max_processor_name=80
  integer, parameter :: mpi_max_error_string=80
  integer, parameter :: mpi_status_size=100
#endif

! Public data

! TODO:  Why are there so many dup-ed communicators?  Untangle logic and 
! TODO:  minimize duplicates...  
  ! MPI intra-communicator for NIM compute or write tasks.
  integer, public :: my_comm=MPI_COMM_NULL
   ! MPI inter-communicator between NIM compute and write tasks.
  integer, public :: write_intercomm=MPI_COMM_NULL
  ! MPI intra-communicator for all available MPI tasks.
  integer :: total_comm=MPI_COMM_NULL
  ! rank of this task in mpi_comm_in or mpi_comm_world
  integer :: rank_comm_in

! TODO:  Ensure use_write_tasks, iam_nim_task, iam_write task dont get used 
! TODO:  before core_setup_nim is called -- make them functions...  
  ! .true. iff write tasks are enabled (init to false)
  logical, public :: use_write_tasks = .false.
  ! .true. iff I am a NIM compute task (init to true)
  logical, public :: iam_nim_task = .true.
  ! .true. iff i am a NIM write task (init to false)
  logical, public :: iam_write_task = .false.
  ! flag to pass to icosio to turn on/off debug messages
  logical, public :: icosio_debugmsg_on = .false.
  ! .true. iff core_setup_nim has finished without error
  logical :: core_setup_done = .false.

  integer,public  :: max_compute_tasks_per_node
  integer,public  :: max_accelerators_per_node

  ! arguments to task_info(), see taskinfo.F90 for descriptions
  integer :: compute_tasks
  integer, public :: cpu_cores_per_node
  integer :: omp_threads_per_compute_task
  integer :: num_write_tasks
  integer :: max_write_tasks_per_node
  integer :: max_compute_tasks_per_mic
  integer :: omp_threads_per_mic_mpi_task
  logical :: root_own_node
  integer :: total_cpu_cores
  integer :: total_mpi_tasks
  integer :: total_nodes
  integer :: mpi_tasks_per_node
  integer :: write_nodes
  integer :: do_nothing_tasks
  integer, pointer :: tasklist_compute(:)
  integer, pointer :: tasklist_write(:)
  integer, pointer :: tasklist_donothing(:)
  integer, pointer :: nodelist(:)

CONTAINS


!*********************************************************************
  subroutine core_setup_nim(mpi_comm_in)
!
! SUMMARY:  
!       When MPI is used, set up communicators for NIM compute tasks and 
!       optional write tasks.  
!       When MPI is not used, return immediately with serial defaults.  
! NOTE:  Executable SMS directives must not be placed before 
! NOTE:  this routine is called.  
! NOTE:  This includes writes or prints without !SMS$IGNORE because they 
! NOTE:  cause SMS to generate code.
!
! ARGUMENTS:  
!       Optional argument mpi_comm_in is an MPI communicator that may be 
!       passed in to restrict NIM to a subset of MPI_COMM_WORLD. 
!
! DETAILS:  
!       Split the MPI communicator and create intercommunicators.  
!     MPI tasks may be split into two or three groups depending upon 
!     settings in NIMnamelist.  The first group will contain the 
!     "compute" tasks responsible for model computations.  The 
!     optional second group will contain "write" tasks responsible 
!     for writing output to disk allowing overlap of computation with 
!     output.  If needed, a third group will contain "do-nothing" tasks 
!     that are neither compute tasks nor write tasks.  The use of 
!     "do-nothing" tasks allows the compute root to optionally live 
!     on its own node (if more memory is needed).  The "do-nothing" 
!     tasks also allow write tasks to be mapped to nodes optionally 
!     leaving some cores unused (if more memory is needed or to allow 
!     concurrent writes of output files from multiple nodes which is 
!     very fast on some machines).  Use of "do-nothing" tasks avoids 
!     dependence upon site-specific features of batch queuing systems 
!     to map MPI tasks onto nodes (e.g. mpich $MACHINE_FILE).  
!       Details of how the MPI communicator is split are described in 
!     taskinfo.F90.  
!*********************************************************************

!TODO:  MPI communicators created here are never freed.  Fix.  

    IMPLICIT NONE

! Arguments
    INTEGER, OPTIONAL, INTENT(IN) :: mpi_comm_in

! Local declarations
    integer, parameter :: tag = 998          ! tag for MPI_Recv
    integer, parameter :: intercomm_tag = 15 ! tag for intercommnicator creation

    LOGICAL :: initialized     ! whether mpi_init has been called
    logical :: nodename_ok
    INTEGER :: size_comm_in    ! size of mpi_comm_in or mpi_comm_world
    INTEGER :: color           ! input arg to mpi_comm_split()--must be >0

    integer :: status(mpi_status_size) ! returned from mpi_recv
    integer :: ignore          ! return code from mpi routines (ignored)
    integer :: ierr            ! return code from mpi routines (not ignored)
    integer :: i,n,nn          ! loop indices
    integer :: unitno          ! unit number for namelist read

    character(len=mpi_max_processor_name) :: mynode ! node name
    character(len=mpi_max_error_string)   :: estring ! error string
    integer :: resultlen       ! length return from MPI routines
    integer :: tmp_comm_in     ! input communicator (MPI_COMM_WORLD or derivative)
    integer :: dup_comm        ! dup of tmp_comm_in
    integer :: comm_active     ! dup_comm modified to contain only "active" 
                               ! (compute or write) tasks except on the 
                               ! "do-nothing" tasks where comm_active contains 
                               ! only "do-nothing" tasks.  
    integer :: rank_active     ! Rank in comm_active
    integer :: write_root_active_rank    ! rank of write root in comm_active
    integer :: wrar            ! task-local temporary copy of above value
    integer :: compute_root_active_rank  ! rank of compute root in comm_active
    integer :: crar            ! task-local temporary copy of above value
    real*8 :: tcore_setup_nim
    integer :: ret
!$$$DEBUG BEGIN
    !integer :: medbg,sizedbg
!$$$DEBUG END
    logical :: i_read
    logical :: write_yaml

! Dynamic arrays

    ! MPI node names
    character(len=mpi_max_processor_name), allocatable :: nodename(:)

! Return early iff this is a serial run (i.e. MPI_INIT has not been called)

    initialized = .false.
#ifndef SERIAL
    call mpi_initialized (initialized, ierr)
#endif
    if (.not.initialized) then
! this is a serial run
      iam_nim_task = .true.
      iam_write_task = .false.
      core_setup_done = .true.
      unitno = getunit()
      if (unitno < 0) then
!SMS$IGNORE BEGIN
        print *,'ERROR in core_setup_nim: returned unit number < 0'
        call flush(6)
        stop
!SMS$IGNORE END
      endif
      i_read = .true.
      write_yaml = .false.
      call task_info(unitno, i_read, write_yaml, compute_tasks,         &
                     cpu_cores_per_node, max_compute_tasks_per_node,    &
                     omp_threads_per_compute_task, num_write_tasks,     &
                     max_write_tasks_per_node, root_own_node,           &
                     icosio_debugmsg_on, max_compute_tasks_per_mic,     &
                     omp_threads_per_mic_mpi_task,                      &
                     total_cpu_cores, total_mpi_tasks, total_nodes,     &
                     mpi_tasks_per_node, write_nodes, do_nothing_tasks, &
                     max_accelerators_per_node, ret)
      if (ret/=0) then
!SMS$IGNORE BEGIN
        print *,'ERROR in core_setup_nim: task_info() returned ',ret
        call flush(6)
        stop
!SMS$IGNORE END
      endif
      call returnunit (unitno)
      return
    endif

! NOTE:  If we reach this point, this is an MPI run.  

    ret = gptlstart ('core_setup_nim')

! use passed-in communicator if present, otherwise use MPI_COMM_WORLD
#ifndef SERIAL
    tmp_comm_in = MPI_COMM_WORLD
#endif
    IF ( PRESENT( mpi_comm_in ) ) THEN
      tmp_comm_in = mpi_comm_in
    ENDIF
! dup is needed to make sends and receives below safe
#ifndef SERIAL
    CALL MPI_COMM_DUP (tmp_comm_in, dup_comm, ierr)
    IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
      PRINT *,'ERROR in core_setup_nim:  MPI_COMM_DUP(tmp_comm_in) returned ',ierr
      call flush(6)
      CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
      STOP
!SMS$IGNORE END
    ENDIF
#endif

#ifndef SERIAL
    call mpi_comm_size (dup_comm, size_comm_in, ierr)
    if (ierr /= 0) then
!SMS$IGNORE BEGIN
      call mpi_error_string (ierr, estring, resultlen, ignore)
      write(6,*)'Error in mpi_comm_size:', estring(1:resultlen)
      call flush (6)
      call mpi_abort (mpi_comm_world, badret, ignore)
!SMS$IGNORE END
    endif
    call mpi_comm_rank (dup_comm, rank_comm_in, ierr)
    if (ierr /= 0) then
!SMS$IGNORE BEGIN
      call mpi_error_string (ierr, estring, resultlen, ignore)
      write(6,*)'Error in mpi_comm_rank:', estring(1:resultlen)
      call flush (6)
      call mpi_abort (mpi_comm_world, badret, ignore)
!SMS$IGNORE END
    endif
#endif

    unitno = getunit()
    if (unitno < 0) then
!SMS$IGNORE BEGIN
      print *,'ERROR in core_setup_nim: returned unit number < 0'
      call flush(6)
#ifndef SERIAL
      call mpi_abort (mpi_comm_world, badret, ignore)
#endif
      stop
!SMS$IGNORE END
    endif
#ifndef SERIAL
    i_read = (rank_comm_in == 0)
    NULLIFY(tasklist_compute)
    NULLIFY(tasklist_write)
    NULLIFY(tasklist_donothing)
    NULLIFY(nodelist)
    write_yaml = .false.
    call task_info(unitno, i_read, write_yaml, compute_tasks,         &
                   cpu_cores_per_node, max_compute_tasks_per_node,    &
                   omp_threads_per_compute_task, num_write_tasks,     &
                   max_write_tasks_per_node, root_own_node,           &
                   icosio_debugmsg_on, max_compute_tasks_per_mic,     &
                   omp_threads_per_mic_mpi_task,                      &
                   total_cpu_cores, total_mpi_tasks, total_nodes,     &
                   mpi_tasks_per_node, write_nodes, do_nothing_tasks, &
                   max_accelerators_per_node, ret, dup_comm,          &
                   tasklist_compute,tasklist_write,tasklist_donothing,&
                   nodelist)
    if (ret/=0) then
!SMS$IGNORE BEGIN
      print *,'ERROR in core_setup_nim: task_info() returned ',ret
      call flush(6)
      call mpi_abort (mpi_comm_world, badret, ignore)
      stop
!SMS$IGNORE END
    endif
#endif
    call returnunit (unitno)

    ! verify size of MPI communicator
    if (size_comm_in /= total_mpi_tasks) then
!SMS$IGNORE BEGIN
      print *,'ERROR in core_setup_nim: size of MPI communicator is incorrect (',size_comm_in,' but expected ',total_mpi_tasks,')'
      call flush(6)
#ifndef SERIAL
      call mpi_abort (mpi_comm_world, badret, ignore)
#endif
      stop
!SMS$IGNORE END
    endif

    iam_nim_task = .false.
    do n = 1,size(tasklist_compute,1)
      if (tasklist_compute(n)==rank_comm_in) then
        iam_nim_task = .true.
      endif
    enddo
    iam_write_task = .false.
    do n = 1,size(tasklist_write,1)
      if (tasklist_write(n)==rank_comm_in) then
        iam_write_task = .true.
      endif
    enddo
!TODO:  remove this after debugging is finished
    if (iam_nim_task.and.iam_write_task) then
!SMS$IGNORE BEGIN
      print *,'ERROR in core_setup_nim: internal error in task distribution'
      call flush(6)
#ifndef SERIAL
      call mpi_abort (mpi_comm_world, badret, ignore)
#endif
      stop
!SMS$IGNORE END
    endif
    use_write_tasks = (num_write_tasks>0)

#ifdef _OPENMP
    !JR Added requirement for root_own_node to be true when threading enabled 
    !JR to avoid the problem on zeus where "omplace" is either not used, or 
    !JR used incorrectly, in a way that causes all threads or all MPI tasks to 
    !JR be placed on a single core of a node.  This can result in a disastrous 
    !JR 10X slowdown!  Not making zeus-specific because the same problem could 
    !JR occur on other machines.  
    ! Note that this problem can occur for an OpenMP build even when 
    ! OMP_NUM_THREADS==1, so the error check must be made here.  
!TBH:  Commented out this check to allow TACC symmetric-mode runs with OpenMP
!TBH:  turned on and root_own_node==.false.  If any vendor has trouble with 
!TBH:  this issue we want to know.  
!TBH:  Do not commit this to the NIM trunk.  
!    if (.not. root_own_node) then
!!SMS$IGNORE BEGIN
!      print *,'ERROR in core_setup_nim: OPENMP requires root_own_node = .true.'
!      call flush(6)
!#ifndef SERIAL
!      call mpi_abort (mpi_comm_world, badret, ignore)
!#endif
!      stop
!!SMS$IGNORE END
!    endif
#else
    ! This check only makes sense for host-only runs!  
    if (max_accelerators_per_node==0) then
      ! omp_threads_per_compute_task>1 only allowed for OpenMP builds
      if (omp_threads_per_compute_task>1) then
!SMS$IGNORE BEGIN
        print *,'ERROR in core_setup_nim: omp_threads_per_compute_task>1 only allowed for OPENMP builds'
        call flush(6)
#ifndef SERIAL
        call mpi_abort (mpi_comm_world, badret, ignore)
#endif
        stop
!SMS$IGNORE END
      endif
    endif
#endif

    ! Use nodelist to verify that node names differ between tasks.  
    ! that are supposed to be mapped to different nodes.  
    ! get my node name
    mynode = ' '
#ifndef SERIAL
    call mpi_get_processor_name (mynode, resultlen, ierr)
    if (ierr /= 0) then
!SMS$IGNORE BEGIN
      call mpi_error_string (ierr, estring, resultlen, ignore)
      write(6,*)'Error from mpi_get_processor_name:', estring(1:resultlen)
      call flush(6)
      call mpi_abort (mpi_comm_world, badret, ignore)
      stop
!SMS$IGNORE END
    endif
#endif

    allocate (nodename(0:size_comm_in-1))
    nodename(0) = mynode
#ifndef SERIAL
    call MPI_Gather(mynode,mpi_max_processor_name,mpi_character,nodename,mpi_max_processor_name,mpi_character,0,dup_comm,ierr)
    if (ierr /= 0) then
!SMS$IGNORE BEGIN
      call mpi_error_string (ierr, estring, resultlen, ignore)
      write(6,*)'Error in mpi_recv:', estring(1:resultlen)
      call flush (6)
      call mpi_abort (mpi_comm_world, badret, ignore)
!SMS$IGNORE END
    endif
#endif

    ! Root task checks that node names change when node number changes
!TODO:  No check for repeated node names -- need to add this?  
    if (rank_comm_in == 0) then
      do n=1,UBOUND(nodelist,1)
        if (nodelist(n-1)==nodelist(n)) then
          if (max_compute_tasks_per_mic == 0) then
            nodename_ok = (trim(nodename(n-1))==trim(nodename(n)))
          else
            ! When running in symmetric mode at TACC, nodename may change 
            ! within a node between MPI tasks on host and MIC.  
            nodename_ok = .true.
          endif
        else
          nodename_ok = (trim(nodename(n-1))/=trim(nodename(n)))
        endif
        if (.not.nodename_ok) then
!SMS$IGNORE BEGIN
          print *,'WARNING in core_setup_nim: internal error, node names and node numbers do not change together'
          print *,n,nodelist(n-1),nodelist(n),trim(nodename(n-1)),' ',trim(nodename(n))
          call flush(6)
#ifndef SERIAL
!JR Cannot get 1 MPI on root, then others on MIC, to work without commenting this out
!          call mpi_abort (mpi_comm_world, badret, ignore)
#endif
!          stop
!SMS$IGNORE END
        endif
      enddo
    endif
    deallocate(nodename)

! Split "do-nothing" tasks from other tasks.  
! Yes, two splits are needed so intercommunicators can be correctly 
! constructed.  
! Define the new communicator comm_active.  
! tasks with color==1 will be in the "do-nothing" group
    color = 1
    if (iam_nim_task.or.iam_write_task) color = 0
#ifndef SERIAL
    CALL MPI_COMM_SPLIT (dup_comm, color, rank_comm_in, comm_active, ierr)
    if (ierr /= 0) then
!SMS$IGNORE BEGIN
      call mpi_error_string (ierr, estring, resultlen, ignore)
      write(6,*)'Error in mpi_comm_split:', estring(1:resultlen)
      call flush (6)
      call mpi_abort (mpi_comm_world, badret, ignore)
!SMS$IGNORE END
    endif
#endif

    if (iam_nim_task.or.iam_write_task) then

#ifndef SERIAL
      call mpi_comm_rank (comm_active, rank_active, ierr)
      if (ierr /= 0) then
!SMS$IGNORE BEGIN
        call mpi_error_string (ierr, estring, resultlen, ignore)
        write(6,*)'Error in mpi_comm_rank:', estring(1:resultlen)
        call flush (6)
        call mpi_abort (mpi_comm_world, badret, ignore)
!SMS$IGNORE END
      endif
#endif
!SMS$IGNORE BEGIN
!JR Namelist hasnt yet been read so control this print with ifdef instead
#ifndef SUPPRESS_WHOIAM_PRINT
      if (iam_nim_task) then
        write(6,'(a,i0,a,i0,a,a)') 'Compute task: rank_comm_in ',rank_comm_in, &
             ' = comm_active rank ', rank_active, ' is running on ', trim(mynode)
      else 
        write(6,'(a,i0,a,i0,a,a)') 'Write task: rank_comm_in ',rank_comm_in, &
             ' = comm_active rank ', rank_active, ' is running on ', trim(mynode)
      endif
#endif
!SMS$IGNORE END
      ! find rank of first compute task ("compute root") in comm_active and 
      ! broadcast to all active tasks for use later to create MPI 
      ! intercommunicators
      if (tasklist_compute(1)==rank_comm_in) then
        crar = rank_active
      else
        crar = -1
      endif
#ifndef SERIAL
      call mpi_allreduce(crar,compute_root_active_rank,1,mpi_integer, &
                         mpi_max,comm_active,ierr)
      if (ierr /= 0) then
!SMS$IGNORE BEGIN
        call mpi_error_string (ierr, estring, resultlen, ignore)
        write(6,*)'Error in mpi_allreduce:', estring(1:resultlen)
        call flush (6)
        call mpi_abort (mpi_comm_world, badret, ignore)
!SMS$IGNORE END
      endif
#endif
      if (use_write_tasks) then
        ! find rank of first write task ("write root") in comm_active and 
        ! broadcast to all active tasks for use later to create MPI 
        ! intercommunicators
        if (tasklist_write(1)==rank_comm_in) then
          wrar = rank_active
        else
          wrar = -1
        endif
#ifndef SERIAL
        call mpi_allreduce(wrar,write_root_active_rank,1,mpi_integer, &
                           mpi_max,comm_active,ierr)
        if (ierr /= 0) then
!SMS$IGNORE BEGIN
          call mpi_error_string (ierr, estring, resultlen, ignore)
          write(6,*)'Error in mpi_allreduce:', estring(1:resultlen)
          call flush (6)
          call mpi_abort (mpi_comm_world, badret, ignore)
!SMS$IGNORE END
        endif
#endif
      endif

    else   ! we are a "do-nothing" task

      my_comm = comm_active   ! used later in !SMS$SET_COMMUNICATOR
!!SMS$IGNORE BEGIN
!      write(6,'(a,i0,a,a)') 'Do-nothing task: rank_comm_in ',rank_comm_in, &
!        ' is running on ', trim(mynode)
!!SMS$IGNORE END

    endif

! Split the remaining MPI communicator into "compute" and "write" 
! sections and create intercommunicators.  
! Not surprisingly, do-nothing tasks skip this.  
    finish_core_setup: if (iam_nim_task.or.iam_write_task) then

!TODO:  Why dup comm_active?  
#ifndef SERIAL
      CALL MPI_COMM_DUP (comm_active, total_comm, ierr)
      IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
        PRINT *,'ERROR in core_setup_nim:  MPI_COMM_DUP returned ',ierr
        call flush(6)
        CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
        STOP
!SMS$IGNORE END
      ENDIF
#endif

! set up MPI communicators for NIM compute and optional write tasks
      IF (.not.use_write_tasks) THEN
! my_comm = total_comm
#ifndef SERIAL
        CALL MPI_COMM_DUP (total_comm, my_comm, ierr)
        IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in core_setup_nim:  MPI_COMM_DUP returned ',ierr
          call flush(6)
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
          STOP
!SMS$IGNORE END
        ENDIF
#endif
      ELSE   ! num_write_tasks > 0
! Split the MPI communicator and create intercommunicators for 
! NIM compute and write tasks.  
! See comments above for details.  
        IF (iam_nim_task) THEN
          color = 0
        ELSE
          color = 1
        ENDIF

#ifndef SERIAL
        CALL MPI_COMM_SPLIT (total_comm, color, rank_active, my_comm, ierr)
        IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in core_setup_nim:  MPI_COMM_SPLIT returned ',ierr
          call flush(6)
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
          STOP
!SMS$IGNORE END
        ENDIF
#endif
        IF (iam_nim_task) THEN
          ! call from NIM compute tasks: "remote leader" is 
          ! write_root_active_rank
#ifndef SERIAL
          CALL MPI_INTERCOMM_CREATE (my_comm, 0, total_comm, &
                                     write_root_active_rank, &
                                     intercomm_tag,          &
                                     write_intercomm, ierr)
          IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
            PRINT *,'ERROR in core_setup_nim:  MPI_INTERCOMM_CREATE returned ',ierr
            call flush(6)
            CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
            STOP
!SMS$IGNORE END
          ENDIF
#endif
        ELSE IF (iam_write_task) THEN
          ! call from NIM write tasks: "remote leader" is 
          ! compute_root_active_rank
!TODO:  extend to multiple groups of write tasks
#ifndef SERIAL
          CALL MPI_INTERCOMM_CREATE (my_comm, 0, total_comm,   &
                                     compute_root_active_rank, &
                                     intercomm_tag,            &
                                     write_intercomm, ierr)
          IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
            PRINT *,'ERROR in core_setup_nim:  MPI_INTERCOMM_CREATE returned ',ierr
            call flush(6)
            CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
            STOP
!SMS$IGNORE END
          ENDIF
#endif
        ELSE
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in core_setup_nim: entered no-mans land where I am neither ', &
            'nim task nor write task nor do-nothing task'
          call flush(6)
#ifndef SERIAL
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
#endif
          STOP
!SMS$IGNORE END
        ENDIF

      ENDIF

!$$$DEBUG BEGIN
#ifndef SERIAL
!SMS$IGNORE BEGIN
!IF (iam_nim_task) THEN
!  call mpi_comm_rank(my_comm,medbg,ignore)
!  call mpi_comm_size(my_comm,sizedbg,ignore)
!  print *,'DEBUG core CT ',medbg,': SIZE(my_comm) = ',sizedbg
!  CALL MPI_BARRIER (my_comm, ignore)
!endif
!IF (iam_write_task) THEN
!  call mpi_comm_rank(my_comm,medbg,ignore)
!  call mpi_comm_size(my_comm,sizedbg,ignore)
!  print *,'DEBUG core WT ',medbg,': SIZE(my_comm) = ',sizedbg
!  CALL MPI_BARRIER (my_comm, ignore)
!endif
!SMS$IGNORE END
#endif
!$$$DEBUG END

    endif finish_core_setup


! TODO Remove references, even in comments, to SMS library internals!

! Finally, tell SMS what communicator to use.  
! All tasks do this passing in their respective intra-communicators.  
! Execution of the SMS SET_COMMUNICATOR directive on all tasks ensures that 
! SMS_INIT() is called properly on the write tasks and on the do-nothing 
! tasks.  This allows these tasks to call sms__stop() without assertion 
! errors.  As a side-effect, duplicate messages are printed by the roots of 
! the write and do-nothing groups during sms__stop().
! NOTE:  Executable SMS directives must not be placed before 
! NOTE:  this directive.  
! NOTE:  This includes writes or prints without !SMS$ignore because they 
! NOTE:  cause SMS to generate code.
!SMS$SET_COMMUNICATOR( my_comm )

    call init_globals()
    call print_nodeinfo(iam_nim_task,iam_write_task)

    ! Do-nothing tasks wait here for other tasks to call MPI_FINALIZE
    ! Note that do-nothing tasks only exist if SMS and MPI are in-use.  
    if ((.not.iam_nim_task).and.(.not.iam_write_task)) then
! NOTE:  In FIM, the do-nothing tasks must return from this routine for NEMS 
! NOTE:  runs due to the ESMF requirement that all tasks in a VM participate 
! NOTE:  in component creation, even tasks that are not used.  If NIM is 
! NOTE:  required to adopt NEMS this behavior will have to be changed to match 
! NOTE:  FIM.  Until then, do-nothing tasks will *not* return from this 
! NOTE:  routine.    
!SMS$INSERT call sms__stop
    endif

! print time only from NIM compute root
    ret = gptlstop ('core_setup_nim')
    ret = gptlget_wallclock ('core_setup_nim', 0, tcore_setup_nim)  ! The "0" is thread number
    IF (iam_nim_task) THEN
      PRINT *,'core_setup_nim time =', tcore_setup_nim
      PRINT"(' Number of Write Tasks:'               ,I24,' processors')",num_write_tasks
    ENDIF

    core_setup_done = .true.

    return
  end subroutine core_setup_nim

  logical function iam_compute_root()
    logical, save :: first_call = .true.
    logical, save :: save_result
    logical       :: mpiinitdone
    integer       :: myrank, ierr, ignore
    if (first_call) then
      mpiinitdone=.false.
      save_result = .false.
#ifndef SERIAL
      call mpi_initialized(mpiinitdone,ierr)
#endif
      if (.not.mpiinitdone) then ! serial run
        save_result=.true.
      else
        if (iam_nim_task) then
! figure out who is the "root" of the compute tasks
#ifndef SERIAL
          CALL MPI_COMM_RANK (my_comm, myrank, ierr)
          IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
            PRINT *,'ERROR in iam_compute_root:  MPI_COMM_RANK returned ',ierr
            call flush(6)
            CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
            STOP
!SMS$IGNORE END
          ENDIF
#endif
          save_result = (myrank == 0)
        endif
      endif
      first_call = .false.
    endif
    iam_compute_root = save_result
  end function iam_compute_root

!TODO:  remove duplication with iam_compute_root
  logical function iam_write_root()
    logical, save :: first_call = .true.
    logical, save :: save_result
    integer       :: myrank, ierr, ignore
    if (first_call) then
      save_result = .false.
      if (.NOT.iam_nim_task) then
! figure out who is the "root" of the write tasks
#ifndef SERIAL
        CALL MPI_COMM_RANK (my_comm, myrank, ierr)
        IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in iam_write_root:  MPI_COMM_RANK returned ',ierr
          call flush(6)
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
          STOP
!SMS$IGNORE END
        ENDIF
#endif
        save_result = (myrank == 0)
      endif
      first_call = .false.
    endif
    iam_write_root = save_result
  end function iam_write_root

end module core_setup

!**************************************************************************
!
!       Nonhydrostatic Icosahedral Model (NIM)
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)!
!
!**************************************************************************
