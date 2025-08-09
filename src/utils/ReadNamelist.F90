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

module ReadNamelist
!******************************************************************************
! All NIM-specific namelist data lives here, along with the routine that reads 
! them.  
! Model-independent namelist data lives in src/utils/.  
!******************************************************************************
  use kinds, only: rt
  implicit none

  public
  save

! QUEUEnamelist
  character(8)   :: MaxQueueTime = '00:10:00' ! Run time for the complete job (HH:MM:SS)
  character(128) :: DataDir    = '/lfs0/projects/nim/DATGFS/'
  character(128) :: RestartDir = 'G5_K32_PHYSICS_P1_1111'

! CNTLnamelist
  integer :: numday              = 5     ! number of integration days
  integer :: numhour             = 0     ! number of integration hours
  integer :: glvl                = 5     ! The grid level
  integer :: gtype               = 2     ! Grid type: Standard recursive(0),Modified recursive (2), Modified great circle (3)
  integer :: SubdivNum(12)       = 2     ! The subdivsion number for each recursive refinement
  integer, parameter :: nz = 96    ! Number of vertical levels
  character(2)  :: ArchvTimeUnit = 'hr'  !  ts:timestep; hr:hour dy:day
  integer :: itsbeg              = 1     ! Timestep starting index
  integer :: RestartBegin        = 0     ! Begin restart if .ne.0
  integer :: ForecastLength      = 100   
  integer :: ArchvIntvl          = 10    ! Archive interval (in ArchvTimeUnit) to do output
  integer :: minmaxPrintInterval = 10    ! Interval to print out MAXs and MINs
  integer :: PrintIpnDiag        = -1    ! ipn at which to print diagnostics (-1 means no print)
  character(12) :: physics       ='NONE' ! Physical package: 'GRIMS' or 'GFS' or 'NONE' for no physics
  logical :: GravityWaveDrag     =.true. ! True means calculate gravity wave drag
  logical :: OutputByHour        =.false.! TRUE = output by hour, FALSE = output by step number
  character(len=12) :: yyyymmddhhmm = "200707170000"      ! Forecast initial time
  real(4) :: pertlim = 0._rt            ! Perturbation to apply to initial temperature (default zero)
  !TBH:  Flag to run "**" operations on CPU for bitwise-exact comparison with 
  !TBH:  GPU.  Running "**" on CPU is slower but exact.  No effect if GPU is 
  !TBH:  not used.  
  logical :: powDouble            =.false. ! True means execute "**" in double precision
  integer :: vdmint_combine      =3  ! This options controls reuse of amtx1.  
                                     ! 3:  Combine u,v into one call (vdmintv) 
                                     !     and combine rp,tp,trp into a 2nd 
                                     !     call (vdmints3).  Default behavior.  
                                     ! 5:  Combine u,v,rp,tp,trp into 1 call 
                                     !     (vdmintsv5).  Maximizes re-use but 
                                     !     does not work on Fermi GPU.  
  integer :: tiles_per_thread = 1    ! Multiplies OMP_NUM_THREADS in GRIMS to obtain num_tiles
  integer :: num_tiles               ! number of WRF-style tiles to use in GRIMS

  logical :: pin_to_single_core = .false.  ! true means pin each MPI rank to a single core, false means pin to a socket
  logical :: root_on_socket1    = .true.   ! true means place the root process on socket 1 (on FGE socket 1 has the Mellanox card)
  logical :: writeOutput = .true.          ! false means do not write any output (do not call output)
  logical :: preExchangeBarrier = .false.  ! true means call a barrier before each exchange
  logical :: postExchangeBarrier = .false. ! true means call a barrier after each exchange
  logical :: OutputBarrier = .false.       ! true means call a barrier before output
  logical :: dyn_phy_barrier = .false.
  logical :: read_amtx = .true.            ! false means compute amtx1 and amtx2 internally, instead of read in from a file
  logical :: TimeInitExchanges = .false.   ! true means output the SMS exchange timings before the main loop
  logical :: zeroSMStimers = .false.       ! true means zero the SMS timers before the main loop starts.
! /POSTnamelist/
  integer ::  numvars       = 12
  character(128)   :: var_list      = "urZZ vrZZ trpZ wrZZ pZZZ tZZZ rZZZ qvZZ stZZ uZZZ vZZZ wZZZ"
  character(1)   :: vert_cord     = 'S'            ! S,Z,P
  character(1)   :: projection    = 'G'
  real           :: center_lat    = 90.
  real           :: center_lon    =  0.
  integer        :: gptx          = 64
  integer        :: gpty          = 64
  real           :: xlen          = 6.e6
  real           :: ylen          = 6.e6

! /PLTVARnamelist/
  real ::  xmin                     = 0.
  real ::  xmax                     = 0.
  real ::  xinc                     = 0.
  character(20)   :: pltvar_list    = "uZZZ 3d XYP0025 00  "

! /PLTICOSnamelist/
  character(128) :: ginfofile ="/scratch1/portfolios/BMC/nim/DATADIR/"
  character(256) :: datafile ="/scratch1/portfolios/BMC/nim/DATADIR/"
  integer :: grid_level = 5
  character(128) :: var_name = "rn2d"    ! name of the variable
  integer :: nvlvls = 1                  ! number of vertical levels of the dataset
  integer :: level = 1                   ! the level of the data to plot
  character(2) :: proj = 'OR'            ! plot projection 
  real    :: latc = 27.344               ! US
  real    :: lonc = -80.0                ! US
  real    :: extent(2) = 180.0           ! the extent of domain for CE proj.
  integer :: map_vis = 1                 ! plot map
  integer :: cell_vis = 0                ! plot voronoi cell
  integer :: ll_vis = 1                  ! plot lat lon lines
  integer :: ipn_label = 0               ! plot ipn index label
  integer :: print_version = 0           ! create graphic file for printing

contains

  subroutine readnl (retval)
!******************************************************************************
! Read NIM namelists. Check validity of input. Print an error message and
! return a non-zero error code if anything bad happens. Since this routine
! may be called from parallel code, no "stop" or "abort" statements are allowed
!
! Jim Rosinski, June 2011
!******************************************************************************

    integer, intent(out) :: retval

    integer :: ioerr
    logical, save :: already_read_namelists = .false.

    namelist /QUEUEnamelist/ MaxQueueTime, DataDir, RestartDir
    namelist /CNTLnamelist/ glvl, gtype, SubdivNum, ArchvTimeUnit, itsbeg, &
      RestartBegin, ForecastLength, ArchvIntvl, minmaxPrintInterval,       &
      PrintIpnDiag, physics, GravityWaveDrag, yyyymmddhhmm,                &
      writeOutput, pertlim, powDouble,outputBarrier,pin_to_single_core,    &
      root_on_socket1, preExchangeBarrier,postExchangeBarrier,             &
      tiles_per_thread, dyn_phy_barrier, read_amtx, vdmint_combine,        &
      TimeInitExchanges,zeroSMStimers
    namelist /POSTnamelist/ numvars,var_list,vert_cord, projection, &
                            center_lat, center_lon, gptx, gpty, xlen, ylen
    namelist /PLTVARnamelist/ xmin, xmax, xinc, pltvar_list
    namelist /PLTICOSnamelist/ ginfofile, datafile, grid_level, var_name, &
    nvlvls, level, proj, latc, lonc, extent, map_vis, cell_vis, ll_vis, &
    ipn_label, print_version

#ifdef _OPENMP
    integer, external :: omp_get_max_threads
#endif

    if (.not. already_read_namelists) then

      retval=opennl()
      if (retval.ne.0) then
        return
      else
        read (unit=10, NML=QUEUEnamelist, iostat=ioerr)
        if (ioerr /= 0) then
          write(6,*)'readnl: error reading QUEUEnamelist from ./NIMnamelist'
          retval = -1
          return
        endif
      endif
      close (10)

      retval=opennl()
      if (retval.ne.0) then
        return
      else
        read (unit=10, NML=CNTLnamelist, iostat=ioerr)
        if (ioerr /= 0) then
          write(6,*)'readnl: error reading CNTLnamelist from ./NIMnamelist'
          retval = -1
          return
        endif
      endif
      close(10)
      call ToUpper(physics)

      retval=opennl()
      if (retval.ne.0) then
        return
      else
        read (unit=10, NML=POSTnamelist, iostat=ioerr)
        if (ioerr /= 0) then
          write(6,*)'readnl: error reading POSTnamelist from ./NIMnamelist'
          retval = -1
          return
        endif
      endif
      close (10)

      retval=opennl()
      if (retval.ne.0) then
        return
      else
        read (unit=10, NML=PLTVARnamelist, iostat=ioerr)
        if (ioerr /= 0) then
          write(6,*)'readnl: error reading PLTVARnamelist from ./NIMnamelist'
          retval = -1
          return
        endif
      endif
      close (10)

      retval=opennl()
      if (retval.ne.0) then
        return
      else
        read (unit=10, NML=PLTICOSnamelist, iostat=ioerr)
        if (ioerr /= 0) then
          write(6,*)'readnl: error reading PLTICOSnamelist from ./NIMnamelist'
          retval = -1
          return
        endif
      endif
      close (10)

      already_read_namelists = .true.

#ifdef _OPENMP
      num_tiles = tiles_per_thread * omp_get_max_threads()
#else
      num_tiles = tiles_per_thread
#endif
    endif

    retval = 0

  end subroutine readnl

  integer function opennl()
    implicit none
    character(len=13) :: nlfile='./NIMnamelist'
    opennl=0 !Work around because ppp does not proadcast opennl
    open (unit=10,file=nlfile,status='old',action='read',iostat=opennl)
    if (opennl.ne.0) then
      write (6,'(a,a)') 'readnl: error opening ',nlfile
    endif
  end function opennl

end module ReadNamelist
