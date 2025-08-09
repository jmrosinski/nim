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

program GetRestart
  use namelistdata, only: GetRestartDir,ReturnRESTARTBEGIN
  implicit none
  integer :: RestartBegin,ret
  character(128) :: RestartDir

  call GetRestartDir (RestartDir, ret)
  if (ret == 0) then
    write (6,'(a,a)') 'RestartDir:', RestartDir
  else
    write(6,'(a)') 'GetRestart: failure from GetRestartDir'
  end if
!
  call ReturnRESTARTBEGIN (RestartBegin, ret)
  if (ret == 0) then
    write(6,'(a,i0)') 'RestartBegin:', RestartBegin
  else
    write(6,'(a)') 'GetRestart: failure from ReturnRESTARTBEGIN'
  end if
end program GetRestart
