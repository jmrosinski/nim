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

program GetPath
!KY  use namelistdata, only: GetDataGrd
  use namelistdata, only: GetDataDir
  implicit none
  integer :: ret
!KY  character(128) :: DataGrd
  character(128) :: DataDir

!KY  call GetDataGrd (DataGrd, ret)
!KY  if (ret == 0) then
!KY    write(6,'(a)') DataGrd
  call GetDataDir (DataDir, ret)
  if (ret == 0) then
    write(6,'(a)') DataDir
  else
    write(6,'(a)') 'GetPath: failure from GetDataGrd'
  end if
end program GetPath
