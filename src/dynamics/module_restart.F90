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

module module_restart

   integer            :: itatm,itsfc
   real*8,parameter   :: rad_inthr=3.0   ! from lw_inthr in interface
   logical,parameter  :: dbg=.true.
   character(len=2)   :: RestartTimeUnit

contains

!-------------------------------------------------------------------------------
   subroutine get_restart(index,fn)
!-------------------------------------------------------------------------------
   use module_control
   use module_constants
   use module_variables
!
   implicit none
!
   integer            :: index,it,its2time,request_timestep,itmp
   integer            :: ios           ! return from I/O call
   character(len=64)  :: fn            ! file name
   character(len=12)  :: fnam
   logical            :: rstfile_e
   logical,parameter  :: dbg_minor=.false.
!-------------------------------------------------------------------------------
!
! check time unit of previous run
!
!sms$ignore begin
   if (index.eq.1) then
     fnam='old.rst_atm_'
     write(6,*) '  Check availability of old.rst_atm file'
   else if (index.eq.2) then
     fnam='old.rst_sfc_'
     write(6,*) '  Check availability of old.rst_sfc file' 
   endif
   do it = 0,10000
     call get_fname(fnam,12,it,fn)
     inquire( file=trim(fn), exist=rstfile_e )
     if ( rstfile_e ) then
       open (30, file=trim(fn), form="unformatted", action='read', iostat=ios)
       read(30) RestartTimeUnit
       read(30) itmp
       rewind(30)
       close(30)
       if( dbg_minor ) then
         write(6,*) 'RestartTimeUnit,itsout in restart.F90 = ', &
                     RestartTimeUnit,itmp
       endif
       goto 777
     end if
   enddo
777 continue
!
! Request time check
!
  if      (ArchvTimeUnit == 'ts') then
    request_timestep = int(RestartBegin)
  else if (ArchvTimeUnit == 'mi') then
    request_timestep = int(RestartBegin*60./dt)
  else if (ArchvTimeUnit == 'hr') then
    request_timestep = int(RestartBegin*3600./dt)
  else if (ArchvTimeUnit == 'dy') then
    request_timestep = int(RestartBegin*86400./dt)
  else if (ArchvTimeUnit == 'mo') then
    request_timestep = int(RestartBegin*hrs_in_month*3600./dt)
  end if
  if (index.eq.1 .and. dbg) then
    write(6,*) '  - You reqeusted restart from',RestartBegin,ArchvTimeUnit
    write(6,*) '   (which is',request_timestep,' timestep)'
  endif
!
! Get restart filename according to requested timestep 
!
  do its=request_timestep,0,-1
    if (mod(dble(its*dt),3600.*rad_inthr) == 0.) then
      if      (RestartTimeUnit == 'ts') then
        its2time=its
      else if (RestartTimeUnit == 'mi') then
        its2time=int(its*dt/60.)
      else if (RestartTimeUnit == 'hr') then
        its2time=int(its*dt/3600.)
      else if (RestartTimeUnit == 'dy') then
        its2time=int(its*dt/86400.)
      else if (RestartTimeUnit == 'mo') then
        its2time=int(its*dt/(hrs_in_month*3600.))
      end if
      call get_fname(fnam,12,its2time,fn)
      inquire( file=trim(fn), exist=rstfile_e )
      if ( rstfile_e ) then
        if (dbg) then
          write(6,*) '  - Here it is! Restart from',its2time,RestartTimeUnit
        endif
        goto 555
      else if ( .not.rstfile_e .and. its.eq.request_timestep .and. dbg) then
        write(6,*) '  - Searching proper restart file ...'
      else if ( .not.rstfile_e .and. its == 0 ) then
        write(6,*) '  - Ooops!! No files for restart exist!!!!!'
        stop
      endif
    else
      if (its.eq.request_timestep .and. dbg) then
        write(6,*) '  - ... Searching proper restart file'
      endif
    endif
  enddo
555 continue
!-------------------------------------------------------------------------------
   return
!sms$ignore end
   end subroutine get_restart
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine get_fname(fnamin,nchar,itime,fnamout)
!-------------------------------------------------------------------------------
   implicit none
!
   integer,intent(in)               :: nchar,itime
   character(len=nchar),intent(in)  :: fnamin
   character(len=64)   ,intent(out) :: fnamout
   logical,parameter                :: dbg_minor=.false.
!-------------------------------------------------------------------------------
!sms$ignore begin
   if (itime.ge.0.and.itime.lt.10) then
     write(fnamout,'(a,a,i1.1)') fnamin,'00000',itime
   else if (itime.ge.10.and.itime.lt.100) then
     write(fnamout,'(a,a,i2.2)') fnamin,'0000',itime
   else if (itime.ge.100.and.itime.lt.1000) then
     write(fnamout,'(a,a,i3.3)') fnamin,'000',itime
   else if (itime.ge.1000.and.itime.lt.10000) then
     write(fnamout,'(a,a,i4.4)') fnamin,'00',itime
   else if (itime.ge.10000.and.itime.lt.100000) then
     write(fnamout,'(a,a,i5.5)') fnamin,'0',itime
   else if (itime.ge.100000.and.itime.lt.1000000) then
     write(fnamout,'(a,a,i6.6)') fnamin,'',itime
   endif
   if( dbg_minor ) then
     write(6,*) 'get_fname : ', fnamout
   endif
!-------------------------------------------------------------------------------
   return
!sms$ignore end
   end subroutine get_fname

end module module_restart
