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

!********************************************************************
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
!	This module specifies control variables for NIM
!  	A. E. MacDonald		October 11, 2004
!  	J. LEE         		September,  2005
!  	J. LEE         		December ,  2008
!********************************************************************

module module_control
use kinds, only: rt
use ReadNamelist
implicit none
save
!
!  MODEL GRID LEVEL: glvl

!  Grid levels for the globe are:

!       glvl    Number of grid points   Linear Scale (km)
!       0               12                7,071
!       1               42                3,779
!       2               162               1,174
!       3               642               891
!       4               2562              446
!       5               10,242            223
!       6               40,962            111
!       7               163,842           56
!       8               655,362           28
!       9               2,621,442         14
!       10              10,485,762        7
!       11              41,943,042        3.5
!       12              167,772,162       1.75

integer , parameter :: nrkl                = 4     ! number of Runge-Kutta time levels
real(rt), parameter :: zt                  = 3.e4_rt
integer , parameter :: nob                 = 9
integer , parameter :: nbf                 = 9
integer , parameter :: nvars               = 6     ! r,ur,vr,wr,tr,trp
integer , parameter :: ntr1                = 3     ! # of tracer in group 1 where 1:qv, 2:qw, 3:O3
integer , parameter :: VarNameLen          = 4
integer , parameter :: filenamelen         = 14
logical , parameter :: PartialStep         =.false.! Do a partial time step
logical , parameter :: FullStep            =.true. ! Do a full time-step. 
logical             :: noOverlap                   ! SMS variable fetched in init: TRUE means do not overlap calc with comm
logical             :: doBarrier                   ! used in preExchange and postExchange
integer             :: numphr                      ! # of time steps/hr
integer             :: nip                         ! # of icosahedral points
integer             :: ips                         ! Start of patch on a given processor
integer             :: ipe                         ! End of patch on a given processor
integer,allocatable :: i_start(:)                  ! tile start indices
integer,allocatable :: i_end(:)                    ! tile end indices
integer             :: ihs                         ! Start of halo on a given processor
integer             :: ihe                         ! End of halo on a given processor
integer             :: ims                         ! Start of memory on a given processor
integer             :: ime                         ! End of memory on a given processor
integer             :: ibx
integer             :: haloSize                    ! Halo size
integer             :: itsend
real(rt)            :: dt
real(rt)            :: cn0
real(rt)            :: cn1
real(rt)            :: dz
real(rt)            :: rhomin
real(rt)            :: qvmin,qwmin,qrmin
real(rt)            :: scl
real(rt)            :: scz
integer             :: ArchvStep                   ! archive interval in time steps
real(rt), parameter :: hrs_in_month   = 730._rt    ! length of month in hrs(= 24*365/12)
integer             :: step_per_unit 
integer             :: nts

contains

subroutine control
integer :: rc

if(nob /= 9) then
  print*,'Error in module_control: nob must be equal to 9',nob
  print*,'9 is hard coded in vdmints.F90 and vdmintv.F90'
  stop
endif
if(nbf /= 9) then
  print*,'Error in module_control: nbf must be equal to 9',nbf
  print*,'9 is hard coded in vdmints.F90 and vdmintv.F90'
  stop
endif

call readnl (rc)
if (rc /= 0) then
  print*,'Error in module_control: bad return from readnl'
  stop
end if

nip    = 10*(2**glvl)**2+2  ! # of icosahedral points
ims    = 1
ime    = nip
ips    = 1
ipe    = nip
!dt     = 900.0/2**(glvl-5)
dt     = 3600.0/2**(glvl-3)
cn0    = 0.4*dt
cn1    = (dt-cn0)
ibx    = nz*nip      ! # of icosahedral boxes
dz     = zt/nz
rhomin = 1.e-5
qvmin  = 1.e-6
qwmin  = 0.  ! 2nd exp, 1.e-9*1.e-3  ! 1st exp: 1.e-9
qrmin  = 0.  ! 2nd exp, 1.e-9*1.e-3  ! 1st exp: 1.e-9
scl    =  2.**glvl/7.68e6
scz    =  1./zt
numphr = 3600./dt+0.5       ! # of time steps/hr

    if (ArchvTimeUnit == 'ts') then
      step_per_unit = 1
!KY    else if (ArchvTimeUnit == 'mi') then
!KY      step_per_unit = nint(60./dt)
    else if (ArchvTimeUnit == 'hr') then
      step_per_unit = nint(3600./dt)
    else if (ArchvTimeUnit == 'dy') then
      step_per_unit = nint(86400./dt)
!KY    else if (ArchvTimeUnit == 'mo') then
!KY      step_per_unit = nint(hrs_in_month*3600./dt)
    else
      write (*,'(a,a)') 'ERROR in module_control unrecognized output time unit:',ArchvTimeUnit
      stop
    end if

    ArchvStep    = ArchvIntvl * step_per_unit

    nts = nint( ForecastLength * step_per_unit + 0.)
!    nts = nint( ForecastLength * step_per_unit + 0.4)

itsend = nts

end subroutine control

end module module_control
