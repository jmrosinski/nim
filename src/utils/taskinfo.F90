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

module taskinfo

  use namelistdata, only: GetNprocs

  implicit none

  logical, save :: write_yaml_mod = .false.

contains

  subroutine task_info(unitno, i_read, write_yaml, compute_tasks_out,     &
                       cpu_cores_per_node_out,                            &
                       max_compute_tasks_per_node_out,                    &
                       omp_threads_per_compute_task_out,                  &
                       num_write_tasks_out, max_write_tasks_per_node_out, &
                       root_own_node_out, icosio_debugmsg_on_out,         &
                       max_compute_tasks_per_mic_out,                     &
                       omp_threads_per_mic_mpi_task_out,                  &
                       total_cpu_cores_out, total_mpi_tasks_out,          &
                       total_nodes_out, mpi_tasks_per_node_out,           &
                       write_nodes_out, do_nothing_tasks_out,             &
                       max_accelerators_per_node_out, status_out, comm_in,&
                       tasklist_compute_out, tasklist_write_out,          &
                       tasklist_donothing_out, nodelist_out)

! Read task and thread layout information from namelist, compute distribution 
! of MPI tasks and OpenMP threads across nodes, and optionally broadcast 
! results to other MPI tasks.  
!
! Note that for parallel (MPI) builds, every CPU core is assigned to either 
! an MPI task or an OpenMP thread.  MPI tasks are divided into three 
! categories:  
!   "compute tasks"     perform model computations
!   "write tasks"       perform concurrent writes of output files enabling 
!                       overlap of model computation with disk writes
!   "do-nothing tasks"  fill in unused CPU cores to avoid site-specific 
!                       methods of mapping MPI tasks to nodes
! The compute "root" (first compute task) is always assigned to the first node 
! due to parallel file system issues observed at some sites.  When write tasks 
! are used they always follow the last compute node.  This simplifies task 
! distribution.  Finally, the same number of MPI tasks are assigned to every 
! node to conform to requirements of the more limited batch systems.  See 
! examples below for more details.  
!
! The following namelist settings (minus the "_out" postfix) are returned 
! from this call:  
!   compute_tasks_out                Total number of MPI compute tasks, or 1 
!                                    if this is a serial run.  
!   cpu_cores_per_node_out           Number of CPU cores per node on this 
!                                    machine.  If deliberate undersubscription 
!                                    or oversubscription of cores is desired, 
!                                    set cpu_cores_per_node smaller or larger 
!                                    than the number of CPU cores per node on
!                                    this machine.  
!   max_compute_tasks_per_node_out   Maximum number of MPI compute tasks to 
!                                    place on a single node.  This is also a 
!                                    cap on max_write_tasks_per_node, see below.
!   omp_threads_per_compute_task_out Number of OMP threads per compute task.  
!                                    This is ignored when OMP is disabled.  
!   num_write_tasks_out              Total number of MPI tasks to dedicate to 
!                                    disk writes.  
!   max_write_tasks_per_node_out     Maximum number of MPI write tasks to place 
!                                    on a single node.  Write tasks do not use 
!                                    OpenMP threads.  Note that this will be 
!                                    reduced to <=(max_compute_tasks_per_node 
!                                    - max_compute_tasks_per_mic)
!                                    when max_write_tasks_per_node is larger 
!                                    than this value.  It is 
!                                    OK to have max_write_tasks_per_node < 
!                                    this value because extra 
!                                    "do-nothing" tasks sharing a node with a 
!                                    write task are not expected to 
!                                    significantly slow down the write task.  
!                                    Similarly, "do-nothing" tasks on the root 
!                                    node have been shown to not slow down 
!                                    OpenMP threads when root_own_node==.true..
!                                    It is possible that "do-nothing" tasks 
!                                    could consume cycles on some machine, so 
!                                    minimizing their impact on the compute 
!                                    nodes seems wise.  
!   root_own_node_out                When .true., the compute root process is 
!                                    placed on its own node.  
!   icosio_debugmsg_on_out           If .true., icosio library will print 
!                                    verbose debugging messages.  
!   max_compute_tasks_per_mic_out    Number of compute tasks per MIC (Intel 
!                                    Xeon Phi) for symmetric mode ONLY.  When 
!                                    max_compute_tasks_per_mic>0 then MICs are 
!                                    assumed to be present.  Currently, MIC 
!                                    support is limited to TACC Stampede which 
!                                    has 1 MIC per node and assigns MPI tasks 
!                                    to the host first then to the MIC on each 
!                                    node before moving to the next node.  
!                                    Each node is required to use 
!                                    max_compute_tasks_per_node MPI tasks with 
!                                    max_compute_tasks_per_mic running on the 
!                                    MIC and (max_compute_tasks_per_node - 
!                                    max_compute_tasks_per_mic) MPI tasks 
!                                    running on the host CPUs.  For example, 
!                                    if max_compute_tasks_per_node=5 and 
!                                    max_compute_tasks_per_mic=3 and 
!                                    compute_tasks=10 then MPI tasks will be 
!                                    assigned to host ("H") and MIC ("M") as 
!                                    follows:  
!                                      Node      0 0 0 0 0 1 1 1 1 1
!                                      Host/MIC  H H M M M H H M M M
!                                      MPI rank  0 1 2 3 4 5 6 7 8 9
!                                    When MIC is used, total_nodes = 
!                                    compute_tasks/max_compute_tasks_per_node
!                                    and a remainder is not allowed.  
!                                    Write tasks will not be mapped to MIC 
!                                    MPI tasks.  Root task will not be mapped 
!                                    to a MIC MPI task.  
!   omp_threads_per_mic_mpi_task_out Number of OpenMP threads to assign to 
!                                    each MIC MPI task.  This must be postive 
!                                    if max_compute_tasks_per_mic>0.  There is  
!                                    no cap on its maximum value (to allow 
!                                    intentional oversubscription of MIC 
!                                    cores) so use with care.  
!
! The following computed values are also returned from this call:  
!   total_cpu_cores_out              Total number of compute cores needed by 
!                                    this run.  This includes compute tasks, 
!                                    write tasks, and any "do-nothing" tasks.  
!   total_mpi_tasks_out              Total number of MPI tasks needed by this 
!                                    run.  This should be passed to the 
!                                    mpirun/mpiexec command when starting the 
!                                    job.  
!   total_nodes_out                  Total number of nodes needed by this run.  
!                                    Some batch systems need this to start a 
!                                    job.  
!   mpi_tasks_per_node_out           Number of MPI tasks to map to each node.  
!                                    This is a convenience 
!                                    (total_mpi_tasks/total_nodes).  
!                                    Note that "do-nothing" tasks are used to 
!                                    ensure that the same number of MPI tasks 
!                                    are mapped to each node.  
!   write_nodes_out                  Number of nodes devoted to write tasks.  
!   do_nothing_tasks_out             Number of "do-nothing" MPI tasks.  
!
! The following optional lists of MPI tasks may also be returned from this 
! call.  
!   tasklist_compute_out
!     If present, optional argument tasklist_compute_out must be an 
!     unassociated pointer to a rank-1 integer array.  It will be pointed to 
!     an array filled with MPI ranks of "compute" tasks.  The caller must not 
!     deallocate tasklist_compute_out.  If this is a serial run, 
!     tasklist_compute_out will point to a size-1 array and the value of its 
!     only element will be set to 0.  
!   tasklist_write_out
!     If present, optional argument tasklist_write_out must be an 
!     unassociated pointer to a rank-1 integer array.  It will be pointed to 
!     an array filled with MPI ranks of "write" tasks.  If the number of write 
!     tasks requested is 0, tasklist_write_out will point to an empty array.  
!     The caller must not deallocate tasklist_write_out.  If this is a serial 
!     run, tasklist_write_out will be allocated with size 0.  
!   tasklist_donothing_out
!     If present, optional argument tasklist_donothing_out must be an 
!     unassociated pointer to a rank-1 integer array.  It will be pointed to 
!     an array filled with MPI ranks of "do nothing" tasks.  If the number of 
!     "do nothing" tasks is 0, tasklist_donothing_out will point to an empty 
!     array.  The caller must not deallocate tasklist_donothing_out.  If this 
!     is a serial run, tasklist_donothing_out will be allocated with size 0.  
! The following optional list of nodes may also be returned from this call.  
!   nodelist_out
!     If present, optional argument nodelist_out must be an unassociated 
!     pointer to a rank-1 integer array.  It will be pointed to an array 
!     filled with the node numbers associated with each MPI rank.  Nodes 
!     numbers start at 0.  Bounds of nodelist_out will be 
!     (0:total_mpi_tasks-1).  This argument is useful for diagnosing bad MPI 
!     task distributions.  The caller must not deallocate nodelist_out.  If 
!     this is a serial run, nodelist_out will point to an empty array.  
!
! Self-consistency checks of namelist settings are made here.  
! 
! The caller must set "i_read" to .true. iff the caller will read the namelist 
! file using logical unit number "unitno" and perform other I/O.  "unitno" is 
! ignored if i_read==.false. .    
!
! If write_yaml and i_read are both .true., namelist settings and returned 
! values will be written to stdout in YAML format.  This feature is used 
! for automated testing and by the NIM run automation.  
! Note that error messages printed by this routine must also be valid YAML 
! so they can also be parsed by the test automation.  
!
! If optional MPI communicator argument comm_in is present, then this routine 
! is assumed to have been called from within an MPI program and will call 
! routine taskinfo_comm() to broadcast values from the root to other MPI 
! tasks.  Note that the number of MPI tasks used and their ranks are deduced 
! from the namelist, *not* from optional argument comm_in.  
!
! If this is a serial build, the namelist setting for ComputeTasks must be set 
! to 1 (or this routine will return with and error message).  The namelist 
! settings for max_compute_tasks_per_node, num_write_tasks, 
! max_write_tasks_per_node, root_own_node, and max_compute_tasks_per_mic are 
! *ignored* and serial defaults are used for intent(out) arguments derived 
! from these values.  The namelist settings for omp_threads_per_compute_task 
! and omp_threads_per_mic_mpi_task are ignored unless OpenMP was enabled 
! during a serial build.  

! Upon return, status_out/=0 indicates an error has occurred.  The caller must 
! abort using appropriate means (i.e. "stop", "mpi_abort()", etc.).  
!
! EXAMPLES:  
!   See taskinfocases.yaml for examples.  
!

    implicit none

    integer, intent(in)  :: unitno
    logical, intent(in)  :: i_read
    logical, intent(in)  :: write_yaml
    integer, intent(out) :: compute_tasks_out
    integer, intent(out) :: cpu_cores_per_node_out
    integer, intent(out) :: max_compute_tasks_per_node_out
    integer, intent(out) :: omp_threads_per_compute_task_out
    integer, intent(out) :: num_write_tasks_out
    integer, intent(out) :: max_write_tasks_per_node_out
    logical, intent(out) :: root_own_node_out
    integer, intent(out) :: max_compute_tasks_per_mic_out
    integer, intent(out) :: max_accelerators_per_node_out
    integer, intent(out) :: omp_threads_per_mic_mpi_task_out
    logical, intent(out) :: icosio_debugmsg_on_out
    integer, intent(out) :: status_out
    integer, intent(out) :: total_cpu_cores_out
    integer, intent(out) :: total_mpi_tasks_out
    integer, intent(out) :: total_nodes_out
    integer, intent(out) :: mpi_tasks_per_node_out
    integer, intent(out) :: write_nodes_out
    integer, intent(out) :: do_nothing_tasks_out
    integer, intent(in ),optional :: comm_in
    integer, optional, pointer :: tasklist_compute_out(:)
    integer, optional, pointer :: tasklist_write_out(:)
    integer, optional, pointer :: tasklist_donothing_out(:)
    integer, optional, pointer :: nodelist_out(:)

    integer       :: istatus
    logical, save :: alreadyReadTaskInfo = .false. ! initial value
    integer, parameter :: UNSET = -9999
    integer :: ComputeTasks_nl
    integer :: cpu_cores_per_node_nl
    integer :: max_compute_tasks_per_node_nl
    integer :: omp_threads_per_compute_task_nl
    integer :: num_write_tasks_nl
    integer :: max_write_tasks_per_node_nl
    logical :: root_own_node_nl
    integer :: max_compute_tasks_per_mic_nl
    integer :: max_accelerators_per_node_nl
    integer :: omp_threads_per_mic_mpi_task_nl
    integer :: compute_nodes
    integer :: n,inode,icompute,iwrite,idonothing,nodecount
    integer :: mpirank,tasksthisnode,tasksleft
    integer :: max_host_tasks_per_node,testcomputenodes

! Namelist variables and default values
    integer, save :: compute_tasks=-1                 ! bad value
    integer, save :: cpu_cores_per_node=UNSET         ! required value
    integer, save :: max_compute_tasks_per_node=UNSET ! required value
    integer, save :: num_write_tasks=0                ! no write tasks
    integer, save :: max_write_tasks_per_node=0       ! no write tasks
    integer, save :: omp_threads_per_compute_task=1   ! 1 OpenMP thread
    logical, save :: icosio_debugmsg_on=.false.       ! no debug prints
    logical, save :: root_own_node=.false.            ! root shares a node
    integer, save :: max_compute_tasks_per_mic=0
    integer, save :: max_accelerators_per_node=0
    integer, save :: omp_threads_per_mic_mpi_task=240
    integer, save :: total_cpu_cores
    integer, save :: total_mpi_tasks
    integer, save :: total_nodes
    integer, save :: mpi_tasks_per_node
    integer, save :: write_nodes
    integer, save :: do_nothing_tasks
    ! These arrays are always allocated (i.e. regardless of presence of 
    ! optional arguments) the first time this routine is called and re-used 
    ! during subsequent calls.  
    integer, pointer, save :: tasklist_compute(:) => NULL()
    integer, pointer, save :: tasklist_write(:) => NULL()
    integer, pointer, save :: tasklist_donothing(:) => NULL()
    integer, pointer, save :: nodelist(:) => NULL()

    namelist /TASKnamelist/ cpu_cores_per_node, max_compute_tasks_per_node,   &
                            num_write_tasks, max_write_tasks_per_node,        &
                            omp_threads_per_compute_task, icosio_debugmsg_on, &
                            root_own_node, max_compute_tasks_per_mic,         &
                            max_accelerators_per_node,                        &
                            omp_threads_per_mic_mpi_task

    status_out = 0
    write_yaml_mod = write_yaml

    if (.not.alreadyReadTaskInfo) then

      if (i_read) then
        ! read compute_tasks from namelist
        call GetNprocs (compute_tasks, istatus)

        if (istatus/=0) then
          write(6,*) 'task_info: GetNprocs failed, istatus=', istatus
          call flush(6)
          status_out = istatus
          return
        endif
        ! save original namelist values for diagnostics
        ComputeTasks_nl = compute_tasks

        open (unitno,file='NIMnamelist',status='old',action='read',iostat=istatus)
        if (istatus/=0) then
          write (*,'(a,i0,a,i0)') 'task_info: failed to open namelist file on unit ',unitno
          call flush(6)
          status_out = istatus
          return
        endif

        read (unitno,TASKnamelist,iostat=istatus)
        if (istatus/=0) then
          write (*,'(a,i0,a,i0)') 'task_info: failed to read TASKnamelist on unit ',unitno
          call flush(6)
          status_out = istatus
          return
        endif

        close(unitno)

        ! save original namelist values for diagnostics
        cpu_cores_per_node_nl = cpu_cores_per_node
        max_compute_tasks_per_node_nl = max_compute_tasks_per_node
        omp_threads_per_compute_task_nl = omp_threads_per_compute_task
        num_write_tasks_nl = num_write_tasks
        max_write_tasks_per_node_nl = max_write_tasks_per_node
        root_own_node_nl = root_own_node
        max_compute_tasks_per_mic_nl = max_compute_tasks_per_mic
        max_accelerators_per_node_nl = max_accelerators_per_node
        omp_threads_per_mic_mpi_task_nl = omp_threads_per_mic_mpi_task

        ! perform self-consistency checks on namelist settings
        if (cpu_cores_per_node == UNSET) then
          call yaml_error('cpu_cores_per_node not found in namelist')
          status_out = 2
          return
        endif
        if (cpu_cores_per_node <= 0) then
          call yaml_error('cpu_cores_per_node must be positive')
          write(6,*)'got ', cpu_cores_per_node
          status_out = 2
          return
        endif
#ifdef _OPENMP
        if (omp_threads_per_compute_task <= 0) then
          call yaml_error('omp_threads_per_compute_task must be positive')
          status_out = 2
          write(6,*)'got ', omp_threads_per_compute_task
          return
        endif
#else
        ! override if not an OpenMP build
        omp_threads_per_compute_task = 1
#endif

      endif

#ifdef SERIAL

      ! if this was a serial *build*

      ! perform self-consistency checks on namelist settings
      if (compute_tasks /= 1) then
        call yaml_error('ComputeTasks must be 1 for a serial build')
        write(6,*)'got ', compute_tasks
        status_out = 2
        return
      endif

#ifdef _OPENMP
      if (omp_threads_per_compute_task > cpu_cores_per_node) then
        call yaml_error('omp_threads_per_compute_task exceeds cpu_cores_per_node')
        write(6,*)'got ', omp_threads_per_compute_task, ' and ', cpu_cores_per_node
        status_out = 2
        return
      endif
#endif
      total_cpu_cores = cpu_cores_per_node

      ! override remaining namelist settings and use serial defaults
      max_compute_tasks_per_node = 1
      num_write_tasks = 0
      max_write_tasks_per_node = 0
      root_own_node = .true.
      max_compute_tasks_per_mic = 0
      icosio_debugmsg_on = .false.
      total_mpi_tasks = 1
      total_nodes = 1
      mpi_tasks_per_node = 1
      write_nodes = 0
      do_nothing_tasks = 0
      ! always allocate and compute tasklist_compute, tasklist_write, and 
      ! tasklist_donothing because some future call may ask for them
      allocate(tasklist_compute(1))
      tasklist_compute(1) = 0
      allocate(tasklist_write(0))
      allocate(tasklist_donothing(0))
      ! always allocate and compute nodelist because some future call may 
      ! ask for it
      allocate(nodelist(0))

#else

      ! if this was not a serial *build*

      if (i_read) then

        ! check for required namelist settings
        if (max_compute_tasks_per_node == UNSET) then
          call yaml_error('max_compute_tasks_per_node not found in namelist')
          status_out = 2
          return
        endif

        ! perform self-consistency checks on namelist settings
        if (compute_tasks <= 0) then
          call yaml_error('ComputeTasks must be positive')
          write(6,*)'got ', compute_tasks
          status_out = 2
          return
        endif
        if (max_compute_tasks_per_node <= 0) then
          call yaml_error('max_compute_tasks_per_node must be positive')
          write(6,*)'got ', max_compute_tasks_per_node
          status_out = 2
          return
        endif
        if (max_compute_tasks_per_mic == 0) then
          if (max_compute_tasks_per_node > cpu_cores_per_node) then
            call yaml_error('max_compute_tasks_per_node > cpu_cores_per_node')
            write(6,*)'got ', max_compute_tasks_per_node, ' and ', cpu_cores_per_node
            status_out = 2
            return
          endif
        endif
        if (num_write_tasks < 0) then
          call yaml_error('num_write_tasks must be non-negative')
          write(6,*)'got ', num_write_tasks
          status_out = 2
          return
        else if (num_write_tasks > 0) then
          ! only check max_write_tasks_per_node if num_write_tasks>0
          if (max_write_tasks_per_node <= 0) then
            call yaml_error('max_write_tasks_per_node must be positive when num_write_tasks > 0')
            write(6,*)'got ', max_write_tasks_per_node, ' and ', num_write_tasks
            status_out = 2
            return
          endif
          if (max_write_tasks_per_node > cpu_cores_per_node) then
            call yaml_error('max_write_tasks_per_node > cpu_cores_per_node')
            write(6,*)'got ', max_write_tasks_per_node, ' and ', cpu_cores_per_node
            status_out = 2
            return
          endif
        else  ! num_write_tasks==0
          ! set to 0 when write tasks are not used
          max_write_tasks_per_node = 0
        endif
        if (max_compute_tasks_per_mic < 0) then
          call yaml_error('max_compute_tasks_per_mic must be non-negative')
          write(6,*)'got ', max_compute_tasks_per_mic
          status_out = 2
          return
        endif
        if (max_compute_tasks_per_mic > 0) then
          if (omp_threads_per_mic_mpi_task < 0) then
            call yaml_error('omp_threads_per_mic_mpi_task must be non-negative')
            write(6,*)'got ', omp_threads_per_mic_mpi_task
            status_out = 2
            return
          endif
        endif

        if (max_compute_tasks_per_mic == 0) then
          ! make automatic adjustments to namelist settings if needed
          ! reduce max_compute_tasks_per_node if needed for small node counts
          if (max_compute_tasks_per_node>compute_tasks) then
            max_compute_tasks_per_node=compute_tasks
          endif
          if (root_own_node.and.(compute_tasks>1)) then
            if (max_compute_tasks_per_node>=compute_tasks) then
              max_compute_tasks_per_node=compute_tasks-1
            endif
          endif
          max_host_tasks_per_node = max_compute_tasks_per_node
        else
          ! We are running on MICs at TACC
          ! In symmetric mode at TACC, at least one MPI rank per node must 
          ! be a host task.  max_compute_tasks_per_node includes 
          ! max_compute_tasks_per_mic.  
          ! In native mode at TACC, all MPI ranks on a node are MIC ranks.  
          if (max_compute_tasks_per_mic>max_compute_tasks_per_node) then
            call yaml_error('max_compute_tasks_per_mic must be less than max_compute_tasks_per_node')
            write(6,*)'got ', max_compute_tasks_per_mic, ' and ', max_compute_tasks_per_node
            status_out = 2
            return
          endif
          max_host_tasks_per_node = max_compute_tasks_per_node - &
                                    max_accelerators_per_node*max_compute_tasks_per_mic
          ! In native mode at TACC, write tasks are not allowed.  
          if ((max_host_tasks_per_node==0).and.(num_write_tasks>0)) then
            call yaml_error('write tasks not allowed with MIC native mode')
            status_out = 2
            return
          endif
        endif
        ! reduce max_write_tasks_per_node if needed
        ! Fit write tasks on the fewest possible nodes such that 
        ! max_write_tasks_per_node<=max_host_tasks_per_node and write tasks 
        ! are distributed across the nodes as evenly as possible.
        if (max_write_tasks_per_node > max_host_tasks_per_node) then
          ! Reduce max_write_tasks_per_node to <=max_host_tasks_per_node to 
          ! avoid oversubscription of cores (i.e. we dont want "do-nothing" 
          ! tasks to compete for CPU resources with OpenMP threads).  
          max_write_tasks_per_node = max_host_tasks_per_node
        endif
        ! further reduce max_write_tasks_per_node if doing so distributes 
        ! write tasks more evenly among the write nodes
        call get_num_nodes(num_write_tasks, max_write_tasks_per_node, &
                           write_nodes)
        if (write_nodes>0) then
          max_write_tasks_per_node = num_write_tasks/write_nodes
          if ((max_write_tasks_per_node*write_nodes)<num_write_tasks) then
            max_write_tasks_per_node = max_write_tasks_per_node + 1
          endif
        else
          max_write_tasks_per_node = 0
        endif

#ifdef _OPENMP
        if (max_compute_tasks_per_mic == 0) then
          ! additional self-consistency checks after automatic adjustments
          if ((max_compute_tasks_per_node*omp_threads_per_compute_task) &
              > cpu_cores_per_node) then
            call yaml_error('max_compute_tasks_per_node*omp_threads_per_compute_task exceeds cpu_cores_per_node')
            write(6,*)'got ', max_compute_tasks_per_node*omp_threads_per_compute_task, ' and ', &
                      cpu_cores_per_node
            status_out = 2
            return
          endif
        endif
#endif

      endif ! i_read

      if (present(comm_in)) then
        ! broadcast namelist settings to all compute tasks
        call taskinfo_comm (compute_tasks, cpu_cores_per_node,              &
                            max_compute_tasks_per_node,                     &
                            omp_threads_per_compute_task, num_write_tasks,  &
                            max_write_tasks_per_node, root_own_node,        &
                            icosio_debugmsg_on, max_compute_tasks_per_mic,  &
                            max_accelerators_per_node,                      &
                            omp_threads_per_mic_mpi_task, comm_in)
      endif

      ! all compute tasks have all needed namelist data at this point
      ! all compute derived values below
      if (max_compute_tasks_per_mic == 0) then
        ! mpi_tasks_per_node
        mpi_tasks_per_node = cpu_cores_per_node / omp_threads_per_compute_task
        ! this assumes max_write_tasks_per_node<=max_compute_tasks_per_node
        mpi_tasks_per_node = min(mpi_tasks_per_node, max_compute_tasks_per_node)
      else
        mpi_tasks_per_node = max_compute_tasks_per_node
      endif
      ! write_nodes (needed on all tasks where i_read == .FALSE.)
      call get_num_nodes(num_write_tasks, max_write_tasks_per_node, &
                         write_nodes)
      ! total_nodes and do_nothing_tasks
      if (root_own_node) then
        call get_num_nodes(compute_tasks-1, max_compute_tasks_per_node, &
                           compute_nodes)
        compute_nodes = compute_nodes + 1
      else
        call get_num_nodes(compute_tasks, max_compute_tasks_per_node, &
                           compute_nodes)
      endif
      total_nodes = write_nodes + compute_nodes
      do_nothing_tasks = (total_nodes*mpi_tasks_per_node) - &
                         (compute_tasks + num_write_tasks)
      ! total_cpu_cores
      total_cpu_cores = total_nodes * cpu_cores_per_node
      ! total_mpi_tasks

      total_mpi_tasks = num_write_tasks + compute_tasks + do_nothing_tasks

      ! always allocate and compute tasklist_compute, tasklist_write, and 
      ! tasklist_donothing because some future call may ask for them
      allocate(tasklist_compute(compute_tasks))
      allocate(tasklist_write(num_write_tasks))
      allocate(tasklist_donothing(do_nothing_tasks))
      ! always allocate and compute nodelist because some future call may 
      ! ask for it
      allocate(nodelist(0:total_mpi_tasks-1))
      ! start with MPI rank 0...
      mpirank = 0
      ! indices into the tasklist_* arrays
      icompute = 1
      iwrite = 1
      idonothing = 1
      nodecount = 0
      ! start node numbering from 0
      inode = 0
      ! compute nodes go first...
      if (root_own_node) then
        ! optionally put compute root on its own node
        call add_node('compute', 1, mpi_tasks_per_node, &
                      icompute, tasklist_compute,       &
                      idonothing, tasklist_donothing,   &
                      inode, nodelist, mpirank,         &
                      status_out)
        if (status_out/=0) return
        nodecount = 1
      endif
      do n=(nodecount+1),compute_nodes
        tasksleft = compute_tasks-icompute+1
        tasksthisnode = min(max_compute_tasks_per_node, tasksleft)
        call add_node('compute', tasksthisnode, mpi_tasks_per_node, &
                      icompute, tasklist_compute,                   &
                      idonothing, tasklist_donothing,               &
                      inode, nodelist, mpirank, status_out)
        if (status_out/=0) return
      enddo
      ! write nodes go last...
      do n=1,write_nodes
        tasksleft = num_write_tasks-iwrite+1
        tasksthisnode = min(max_write_tasks_per_node, tasksleft)
        call add_node('write', tasksthisnode, mpi_tasks_per_node, &
                      iwrite, tasklist_write,                     &
                      idonothing, tasklist_donothing,             &
                      inode, nodelist, mpirank, status_out)
        if (status_out/=0) return
      enddo

#endif

      alreadyReadTaskInfo=.true.

    endif ! .not.alreadyReadTaskInfo

    compute_tasks_out                = compute_tasks
    cpu_cores_per_node_out           = cpu_cores_per_node
    icosio_debugmsg_on_out           = icosio_debugmsg_on
    max_write_tasks_per_node_out     = max_write_tasks_per_node
    num_write_tasks_out              = num_write_tasks
    root_own_node_out                = root_own_node
    max_compute_tasks_per_mic_out    = max_compute_tasks_per_mic
    max_accelerators_per_node_out    = max_accelerators_per_node
    omp_threads_per_mic_mpi_task_out = omp_threads_per_mic_mpi_task
    max_compute_tasks_per_node_out   = max_compute_tasks_per_node
    omp_threads_per_compute_task_out = omp_threads_per_compute_task
    total_cpu_cores_out              = total_cpu_cores
    total_mpi_tasks_out              = total_mpi_tasks
    total_nodes_out                  = total_nodes
    mpi_tasks_per_node_out           = mpi_tasks_per_node
    write_nodes_out                  = write_nodes
    do_nothing_tasks_out             = do_nothing_tasks
    if (present(tasklist_compute_out)) then
      tasklist_compute_out => tasklist_compute
    endif
    if (present(tasklist_write_out)) then
      tasklist_write_out => tasklist_write
    endif
    if (present(tasklist_donothing_out)) then
      tasklist_donothing_out => tasklist_donothing
    endif
    if (present(nodelist_out)) then
      nodelist_out => nodelist
    endif

    ! YAML prints mimic documentation block for automated unit testing
    if (write_yaml_mod) then
      if (i_read) then
        write(6,'(a)') '---'
#ifdef SERIAL
        write(6,'(a)') 'build: serial'
#else
        write(6,'(a)') 'build: parallel'
#endif
        write(6,'(a)') 'namelist_settings:'
        write(6,'(a,i0,a,i0)') '  ComputeTasks: ',ComputeTasks_nl
        write(6,'(a,i0,a,i0)') '  cpu_cores_per_node: ',cpu_cores_per_node_nl
        write(6,'(a,i0,a,i0)') '  max_compute_tasks_per_node: ',max_compute_tasks_per_node_nl
        write(6,'(a,i0,a,i0)') '  omp_threads_per_compute_task: ',omp_threads_per_compute_task_nl
        write(6,'(a,i0,a,i0)') '  num_write_tasks: ',num_write_tasks_nl
        write(6,'(a,i0,a,i0)') '  max_write_tasks_per_node: ',max_write_tasks_per_node_nl
        write(6,'(a,i0,a,i0)') '  max_compute_tasks_per_mic: ',max_compute_tasks_per_mic_nl
        write(6,'(a,i0,a,i0)') '  max_accelerators_per_node: ',max_accelerators_per_node_nl
        write(6,'(a,i0,a,i0)') '  omp_threads_per_mic_mpi_task: ',omp_threads_per_mic_mpi_task_nl
        write(6,'(a,l1,a,l1)') '  root_own_node: ',root_own_node_nl
        write(6,'(a)') 'returned_values:'
        write(6,'(a,i0)') '  actual_max_compute_tasks_per_node: ',max_compute_tasks_per_node
        write(6,'(a,i0)') '  actual_max_write_tasks_per_node: ',max_write_tasks_per_node
        write(6,'(a,i0)') '  total_cpu_cores: ',total_cpu_cores
        write(6,'(a,i0)') '  total_mpi_tasks: ',total_mpi_tasks
        write(6,'(a,i0)') '  total_nodes: ',total_nodes
        write(6,'(a,i0)') '  mpi_tasks_per_node: ',mpi_tasks_per_node
        write(6,'(a,i0)') '  write_nodes: ',write_nodes
        write(6,'(a,i0)') '  do_nothing_tasks: ',do_nothing_tasks
        call printlist('  tasklist_compute',tasklist_compute)
        call printlist('  tasklist_write',tasklist_write)
        call printlist('  tasklist_donothing',tasklist_donothing)
        call printlist('  nodelist',nodelist)
        write(6,'(a)') '---'
        call flush(6)
      endif
    endif

    return
  end subroutine task_info


  ! compute the number of nodes needed
  subroutine get_num_nodes(num_tasks, max_tasks_per_node, nodes)
    integer, intent(in   ) :: num_tasks
    integer, intent(in   ) :: max_tasks_per_node
    integer, intent(inout) :: nodes
    if (num_tasks > 0) then
      nodes = num_tasks/max_tasks_per_node
      if ((nodes*max_tasks_per_node)<num_tasks) then
        nodes = nodes + 1
      endif
    else
      nodes = 0
    endif
    return
  end subroutine get_num_nodes

  ! add current MPI rank to tasklist.  Also add to the node list.  
  subroutine add_task(listname, itask, mpirank, tasklist, & 
                      inode, nodelist, status_out)
    character(len=*), intent(in   ) :: listname
    integer,          intent(inout) :: itask
    integer,          intent(inout) :: mpirank
    integer,          intent(inout) :: tasklist(:)
    integer,          intent(in   ) :: inode
    integer,          intent(inout) :: nodelist(0:)
    integer,          intent(  out) :: status_out
    status_out = 0
    ! bounds-check
    if (itask>ubound(tasklist,1)) then
      call yaml_error(trim(listname),' INTERNAL ERROR: task list too small')
      status_out = 3
      return
    endif
    ! add MPI rank to task list
    tasklist(itask) = mpirank
    ! add node for this MPI rank to node list
    ! bounds-check
    if (mpirank>ubound(nodelist,1)) then
      call yaml_error('INTERNAL ERROR: nodelist too small')
      status_out = 3
      return
    endif
    nodelist(mpirank) = inode
    ! increment index and rank
    itask = itask + 1
    mpirank = mpirank + 1
    return
  end subroutine add_task

  ! add a node by adding tasks to either the compute or write task list and 
  ! to the do-nothing task list.  Also add to the node list.  
  subroutine add_node(nodename, tasksthisnode, mpi_tasks_per_node,     &
                      itask, tasklist, idonothing, tasklist_donothing, &
                      inode, nodelist, mpirank, status_out)
    character(len=*), intent(in   ) :: nodename
    integer,          intent(in   ) :: tasksthisnode
    integer,          intent(in   ) :: mpi_tasks_per_node
    integer,          intent(inout) :: itask
    integer,          intent(inout) :: tasklist(:)
    integer,          intent(inout) :: idonothing
    integer,          intent(inout) :: tasklist_donothing(:)
    integer,          intent(inout) :: inode
    integer,          intent(inout) :: nodelist(0:)
    integer,          intent(inout) :: mpirank
    integer,          intent(  out) :: status_out
    integer :: i
    status_out = 0
    do i=1,tasksthisnode
      call add_task(trim(nodename), itask, mpirank, &
                    tasklist, inode, nodelist,      &
                    status_out)
      if (status_out/=0) return
    enddo
    do i=tasksthisnode+1,mpi_tasks_per_node
      call add_task('donothing', idonothing, mpirank,    &
                    tasklist_donothing, inode, nodelist, &
                    status_out)
      if (status_out/=0) return
    enddo
    ! increment node index
    inode = inode + 1
    return
  end subroutine add_node

  ! Print array values
  ! More than two consecutive values are abbreviated using a dash.  For 
  ! example,  
  !   (3,4,5,6,7) becomes (3-7)
  !   (10,11,13)  becomes (10-11,13)
  !   (3,4)       remains (3,4).  
  subroutine printlist(listname,list)
    character(len=*), intent(in) :: listname
    integer,          intent(in) :: list(:)
    integer :: i,last_consecutive_value
    integer :: consecutive
    consecutive = 0
    write(6,'(2a)',advance='no') trim(listname),': ('
    do i=1,size(list)
      if (i>1) then
        if (list(i)==list(i-1)+1) then
          ! consecutive values
          consecutive = consecutive + 1
          last_consecutive_value = list(i)
        else
          if (consecutive>0) then
            if (consecutive>1) then
              ! dash separator
              write(6,'(a)',advance='no') '-'
            else
              ! comma separator
              write(6,'(a)',advance='no') ','
            endif
            write(6,'(i0)',advance='no') last_consecutive_value
          endif
          consecutive = 0
          ! comma separator
          write(6,'(a)',advance='no') ','
          write(6,'(i0)',advance='no') list(i)
        endif
      else
        write(6,'(i0)',advance='no') list(i)
      endif
    enddo
    if (consecutive>0) then
      if (consecutive>1) then
        ! dash separator
        write(6,'(a)',advance='no') '-'
      else
        ! comma separator
        write(6,'(a)',advance='no') ','
      endif
      write(6,'(i0)',advance='no') last_consecutive_value
    endif
    write(6,'(a)') ')'
    return
  end subroutine printlist

  ! Write error message, possibly in YAML format used by test suite
  !TODO:  A better approach *might* be to split YAML into a separate file  
  !TODO:  but this is also not ideal since Fortran95 does not handle ARGV 
  !TODO:  portably...  
  subroutine yaml_error(msg, msg2)
    character(len=*), intent(in)           :: msg
    character(len=*), intent(in), optional :: msg2
    if (write_yaml_mod) then
      write(6,'(a)') '---'
      write(6,'(a)') 'returned_values:'
    endif
    if (present(msg2)) then
      write(6,'(a,a,a)') '  ERROR: task_info ', msg, msg2
    else
      write(6,'(a,a)') '  ERROR: task_info ', msg
    endif
    ! this stops YAML processing so any text printed after this need not be 
    ! in YAML format
    if (write_yaml_mod) then
      write(6,'(a)') '---'
    endif
    call flush(6)
  end subroutine yaml_error

end module taskinfo

