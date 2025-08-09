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

! Real kinds for NIM dynamics computations.  

module kinds
implicit none

save

integer, parameter, public :: sp = selected_real_kind(6)    ! REAL*4
integer, parameter, public :: dp = selected_real_kind(15)   ! REAL*8

#ifdef DOUBLEPREC
integer, parameter, public :: rt = dp
#else
integer, parameter, public :: rt = sp
#endif
 
end module kinds

