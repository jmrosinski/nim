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

subroutine ToUpper(str)
character(*), intent(in out) :: str
integer :: i
 
do i = 1, len(str)
  select case(str(i:i))
    case("a":"z")
      str(i:i) = achar(iachar(str(i:i))-32)
  end select
end do 
end subroutine ToUpper
