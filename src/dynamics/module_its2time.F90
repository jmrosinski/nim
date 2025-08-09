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

module module_its2time
  use module_control,only:ArchvTimeUnit,dt,hrs_in_month
  implicit none
!  real         (out):: its2time

contains
  function its2time(its)
    implicit none
    integer,intent(in):: its
    real              :: its2time

!sms$ignore begin
! Changes in this block may require changes in dyn_init.F90 as well
  if (ArchvTimeUnit.eq.'ts') then
    its2time=float(its)
!KY  else if (ArchvTimeUnit.eq.'mi') then
!KY    its2time=its*dt/60.
  else if (ArchvTimeUnit.eq.'hr') then
    its2time=its*dt/3600.
  else if (ArchvTimeUnit.eq.'dy') then
    its2time=its*dt/86400.
!KY  else if (ArchvTimeUnit.eq.'mo') then
!KY    its2time=its*dt/(hrs_in_month*3600.)
  else
    write (*,'(a,a)') 'ERROR in its2time: unrecognized output time unit: ',ArchvTimeUnit
    stop
  endif
!sms$ignore end
end function its2time

end module module_its2time
