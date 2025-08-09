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

program GetComputeTasks
  use namelistdata, only: GetNprocs
  implicit none
  integer :: nprocs, ret

  call GetNprocs (nprocs, ret)
  if (ret == 0) then
    write(6,"(i0)") nprocs
  else
    write(6,*)'GetComputeTasks: failure from GetNprocs'
  end if
end program GetComputeTasks
