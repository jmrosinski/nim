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

program GetQueueTime
  use namelistdata, only: GetMaxQueueTime
  implicit none
  character(8) :: QueueTime
  integer :: ret

  call GetMaxQueueTime (QueueTime, ret)
  if (ret == 0) then
    write(6,'(a)') QueueTime

  else
    write(6,'(a)') 'GetQueueTime: failure from GetMaxQueueTime'
  end if
end program GetQueueTime
