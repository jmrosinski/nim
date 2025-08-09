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

subroutine read1Di4(varName,unit,dim1,var)
implicit none

character(*),intent(in)  :: varName
integer     ,intent(in)  :: unit
integer     ,intent(in)  ::     dim1
integer     ,intent(out) :: var(dim1)
integer                  :: ierr

read(unit, iostat=ierr) var
if (ierr.ne.0) then
  write (*,'(a,a)') 'read1Di4: read error reading ',varName
  stop
endif

end subroutine read1Di4
