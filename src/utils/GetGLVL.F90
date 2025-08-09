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

program GetGLVL
  use namelistdata, only: ReturnGLVL
  implicit none
  integer :: glvl, ret

  call ReturnGLVL(glvl, ret)
  if (ret == 0) then
    write(6,"(i0)") glvl
  else
    write(6,*)'GetGLVL: failure from ReturnGLVL'
  end if
end program GetGLVL

