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

module infnan

  implicit none

  private
  public :: inf, nan, negint

! TODO: add big endian ifdef for e.g. IBM
! TODO: Should NaN be signaling or non-signaling?

  real,    parameter :: inf = z'7f800000'   ! Infinity
  real,    parameter :: nan = z'ffc00000'   ! NaN
  integer, parameter :: negint = -999       ! Bad integer value

end module infnan
