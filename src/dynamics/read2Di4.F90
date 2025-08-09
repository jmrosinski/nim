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

subroutine read2Di4(varName,unit,dim1,dim2,var)
implicit none

character(*),intent(in)  :: varName
integer     ,intent(in)  :: unit
integer     ,intent(in)  ::     dim1,dim2
integer     ,intent(out) :: var(dim1,dim2)
integer                  :: i,ierr

do i=1,dim1
  read(unit, iostat=ierr) var(i,:)
  if (ierr.ne.0) then
    write (*,'(a,a)') 'read2Di4: read error reading ',varName
    stop
  endif
enddo

end subroutine read2Di4
