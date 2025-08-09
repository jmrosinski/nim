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

! Substitutes for non-standard system intrinsic subroutines.  

module module_sys_share
  implicit none
end module module_sys_share

! For IBM #define IBMFLUSH and build this proxy.  
! IBM uses flush() instead of call flush()
#ifdef IBMFLUSH
subroutine flush(lun)
  implicit none
  integer,intent(in)::lun
  flush(lun) ! for IBM, generalize later if needed
end subroutine flush
#endif /* IBMFLUSH */

subroutine stub_to_satisfy_linkers_that_cant_handle_empty_doto_files
end subroutine stub_to_satisfy_linkers_that_cant_handle_empty_doto_files
