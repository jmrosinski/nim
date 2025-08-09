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

program GetPhysics
  use namelistdata, only: ReturnPhysics
  implicit none
  integer       :: ret
  character(12) :: physics ! Physical package: GRIMS or GFS or none for no physics


  call ReturnPhysics(physics, ret)
  if (ret == 0) then
    write(6,"(A)") physics
  else
    write(6,*)'GetPhysics: failure from ReturnPhysics'
  end if
end program GetPhysics

