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


#ifdef SOLVEITHLS_RECIPROCAL
#define _OP_ *
#else
#define _OP_ /
#endif

SUBROUTINE solveiThLS3(nob,nbf,b1,b2,b3,amtx)

use kinds, only: rt
!use ReadNamelist, only: nz
IMPLICIT NONE
!MWG: workaround due to F2C bug
  integer, parameter     :: nz=96

integer, intent(IN)    :: nbf,nob
real(rt), INTENT(INOUT):: b1(NZ,nob),b2(NZ,nob),b3(NZ,nob)
real(rt), intent(IN)   :: amtx(NZ,nbf,nob+1)
real(rt)               :: wk
INTEGER                :: k
#ifdef FINEGRAINED_TIMING
  integer :: ret
  integer :: handle = 0
#include <gptl.inc>

  ret = gptlstart_handle ('solveiThLS3', handle)
#endif
!$acc routine(solveiThLS3) vector
!$acc loop vector
!dir$ vector aligned
  do k=1,NZ
      wk =                b1(k,1) &
          + amtx(k,2,1) * b1(k,2) &
          + amtx(k,3,1) * b1(k,3) &
          + amtx(k,4,1) * b1(k,4) &
          + amtx(k,5,1) * b1(k,5) &
          + amtx(k,6,1) * b1(k,6) &
          + amtx(k,7,1) * b1(k,7) &
          + amtx(k,8,1) * b1(k,8) &
          + amtx(k,9,1) * b1(k,9)
      b1(k,1) = b1(k,1) - amtx(k,1,10) * wk
      b1(k,2) = b1(k,2) - amtx(k,1,10) * wk * amtx(k,2,1)
      b1(k,3) = b1(k,3) - amtx(k,1,10) * wk * amtx(k,3,1)
      b1(k,4) = b1(k,4) - amtx(k,1,10) * wk * amtx(k,4,1)
      b1(k,5) = b1(k,5) - amtx(k,1,10) * wk * amtx(k,5,1)
      b1(k,6) = b1(k,6) - amtx(k,1,10) * wk * amtx(k,6,1)
      b1(k,7) = b1(k,7) - amtx(k,1,10) * wk * amtx(k,7,1)
      b1(k,8) = b1(k,8) - amtx(k,1,10) * wk * amtx(k,8,1)
      b1(k,9) = b1(k,9) - amtx(k,1,10) * wk * amtx(k,9,1)
      wk =                b2(k,1) &
          + amtx(k,2,1) * b2(k,2) &
          + amtx(k,3,1) * b2(k,3) &
          + amtx(k,4,1) * b2(k,4) &
          + amtx(k,5,1) * b2(k,5) &
          + amtx(k,6,1) * b2(k,6) &
          + amtx(k,7,1) * b2(k,7) &
          + amtx(k,8,1) * b2(k,8) &
          + amtx(k,9,1) * b2(k,9)
      b2(k,1) = b2(k,1) - amtx(k,1,10) * wk
      b2(k,2) = b2(k,2) - amtx(k,1,10) * wk * amtx(k,2,1)
      b2(k,3) = b2(k,3) - amtx(k,1,10) * wk * amtx(k,3,1)
      b2(k,4) = b2(k,4) - amtx(k,1,10) * wk * amtx(k,4,1)
      b2(k,5) = b2(k,5) - amtx(k,1,10) * wk * amtx(k,5,1)
      b2(k,6) = b2(k,6) - amtx(k,1,10) * wk * amtx(k,6,1)
      b2(k,7) = b2(k,7) - amtx(k,1,10) * wk * amtx(k,7,1)
      b2(k,8) = b2(k,8) - amtx(k,1,10) * wk * amtx(k,8,1)
      b2(k,9) = b2(k,9) - amtx(k,1,10) * wk * amtx(k,9,1)
      wk =                b3(k,1) &
          + amtx(k,2,1) * b3(k,2) &
          + amtx(k,3,1) * b3(k,3) &
          + amtx(k,4,1) * b3(k,4) &
          + amtx(k,5,1) * b3(k,5) &
          + amtx(k,6,1) * b3(k,6) &
          + amtx(k,7,1) * b3(k,7) &
          + amtx(k,8,1) * b3(k,8) &
          + amtx(k,9,1) * b3(k,9)
      b3(k,1) = b3(k,1) - amtx(k,1,10) * wk
      b3(k,2) = b3(k,2) - amtx(k,1,10) * wk * amtx(k,2,1)
      b3(k,3) = b3(k,3) - amtx(k,1,10) * wk * amtx(k,3,1)
      b3(k,4) = b3(k,4) - amtx(k,1,10) * wk * amtx(k,4,1)
      b3(k,5) = b3(k,5) - amtx(k,1,10) * wk * amtx(k,5,1)
      b3(k,6) = b3(k,6) - amtx(k,1,10) * wk * amtx(k,6,1)
      b3(k,7) = b3(k,7) - amtx(k,1,10) * wk * amtx(k,7,1)
      b3(k,8) = b3(k,8) - amtx(k,1,10) * wk * amtx(k,8,1)
      b3(k,9) = b3(k,9) - amtx(k,1,10) * wk * amtx(k,9,1)
      wk =                b1(k,2) &
          + amtx(k,3,2) * b1(k,3) &
          + amtx(k,4,2) * b1(k,4) &
          + amtx(k,5,2) * b1(k,5) &
          + amtx(k,6,2) * b1(k,6) &
          + amtx(k,7,2) * b1(k,7) &
          + amtx(k,8,2) * b1(k,8) &
          + amtx(k,9,2) * b1(k,9)
      b1(k,2) = b1(k,2) - amtx(k,2,10) * wk
      b1(k,3) = b1(k,3) - amtx(k,2,10) * wk * amtx(k,3,2)
      b1(k,4) = b1(k,4) - amtx(k,2,10) * wk * amtx(k,4,2)
      b1(k,5) = b1(k,5) - amtx(k,2,10) * wk * amtx(k,5,2)
      b1(k,6) = b1(k,6) - amtx(k,2,10) * wk * amtx(k,6,2)
      b1(k,7) = b1(k,7) - amtx(k,2,10) * wk * amtx(k,7,2)
      b1(k,8) = b1(k,8) - amtx(k,2,10) * wk * amtx(k,8,2)
      b1(k,9) = b1(k,9) - amtx(k,2,10) * wk * amtx(k,9,2)
      wk =                b2(k,2) &
          + amtx(k,3,2) * b2(k,3) &
          + amtx(k,4,2) * b2(k,4) &
          + amtx(k,5,2) * b2(k,5) &
          + amtx(k,6,2) * b2(k,6) &
          + amtx(k,7,2) * b2(k,7) &
          + amtx(k,8,2) * b2(k,8) &
          + amtx(k,9,2) * b2(k,9)
      b2(k,2) = b2(k,2) - amtx(k,2,10) * wk
      b2(k,3) = b2(k,3) - amtx(k,2,10) * wk * amtx(k,3,2)
      b2(k,4) = b2(k,4) - amtx(k,2,10) * wk * amtx(k,4,2)
      b2(k,5) = b2(k,5) - amtx(k,2,10) * wk * amtx(k,5,2)
      b2(k,6) = b2(k,6) - amtx(k,2,10) * wk * amtx(k,6,2)
      b2(k,7) = b2(k,7) - amtx(k,2,10) * wk * amtx(k,7,2)
      b2(k,8) = b2(k,8) - amtx(k,2,10) * wk * amtx(k,8,2)
      b2(k,9) = b2(k,9) - amtx(k,2,10) * wk * amtx(k,9,2)
      wk =                b3(k,2) &
          + amtx(k,3,2) * b3(k,3) &
          + amtx(k,4,2) * b3(k,4) &
          + amtx(k,5,2) * b3(k,5) &
          + amtx(k,6,2) * b3(k,6) &
          + amtx(k,7,2) * b3(k,7) &
          + amtx(k,8,2) * b3(k,8) &
          + amtx(k,9,2) * b3(k,9)
      b3(k,2) = b3(k,2) - amtx(k,2,10) * wk
      b3(k,3) = b3(k,3) - amtx(k,2,10) * wk * amtx(k,3,2)
      b3(k,4) = b3(k,4) - amtx(k,2,10) * wk * amtx(k,4,2)
      b3(k,5) = b3(k,5) - amtx(k,2,10) * wk * amtx(k,5,2)
      b3(k,6) = b3(k,6) - amtx(k,2,10) * wk * amtx(k,6,2)
      b3(k,7) = b3(k,7) - amtx(k,2,10) * wk * amtx(k,7,2)
      b3(k,8) = b3(k,8) - amtx(k,2,10) * wk * amtx(k,8,2)
      b3(k,9) = b3(k,9) - amtx(k,2,10) * wk * amtx(k,9,2)
      wk =                b1(k,3) &
          + amtx(k,4,3) * b1(k,4) &
          + amtx(k,5,3) * b1(k,5) &
          + amtx(k,6,3) * b1(k,6) &
          + amtx(k,7,3) * b1(k,7) &
          + amtx(k,8,3) * b1(k,8) &
          + amtx(k,9,3) * b1(k,9)
      b1(k,3) = b1(k,3) - amtx(k,3,10) * wk
      b1(k,4) = b1(k,4) - amtx(k,3,10) * wk * amtx(k,4,3)
      b1(k,5) = b1(k,5) - amtx(k,3,10) * wk * amtx(k,5,3)
      b1(k,6) = b1(k,6) - amtx(k,3,10) * wk * amtx(k,6,3)
      b1(k,7) = b1(k,7) - amtx(k,3,10) * wk * amtx(k,7,3)
      b1(k,8) = b1(k,8) - amtx(k,3,10) * wk * amtx(k,8,3)
      b1(k,9) = b1(k,9) - amtx(k,3,10) * wk * amtx(k,9,3)
      wk =                b2(k,3) &
          + amtx(k,4,3) * b2(k,4) &
          + amtx(k,5,3) * b2(k,5) &
          + amtx(k,6,3) * b2(k,6) &
          + amtx(k,7,3) * b2(k,7) &
          + amtx(k,8,3) * b2(k,8) &
          + amtx(k,9,3) * b2(k,9)
      b2(k,3) = b2(k,3) - amtx(k,3,10) * wk
      b2(k,4) = b2(k,4) - amtx(k,3,10) * wk * amtx(k,4,3)
      b2(k,5) = b2(k,5) - amtx(k,3,10) * wk * amtx(k,5,3)
      b2(k,6) = b2(k,6) - amtx(k,3,10) * wk * amtx(k,6,3)
      b2(k,7) = b2(k,7) - amtx(k,3,10) * wk * amtx(k,7,3)
      b2(k,8) = b2(k,8) - amtx(k,3,10) * wk * amtx(k,8,3)
      b2(k,9) = b2(k,9) - amtx(k,3,10) * wk * amtx(k,9,3)
      wk =                b3(k,3) &
          + amtx(k,4,3) * b3(k,4) &
          + amtx(k,5,3) * b3(k,5) &
          + amtx(k,6,3) * b3(k,6) &
          + amtx(k,7,3) * b3(k,7) &
          + amtx(k,8,3) * b3(k,8) &
          + amtx(k,9,3) * b3(k,9)
      b3(k,3) = b3(k,3) - amtx(k,3,10) * wk
      b3(k,4) = b3(k,4) - amtx(k,3,10) * wk * amtx(k,4,3)
      b3(k,5) = b3(k,5) - amtx(k,3,10) * wk * amtx(k,5,3)
      b3(k,6) = b3(k,6) - amtx(k,3,10) * wk * amtx(k,6,3)
      b3(k,7) = b3(k,7) - amtx(k,3,10) * wk * amtx(k,7,3)
      b3(k,8) = b3(k,8) - amtx(k,3,10) * wk * amtx(k,8,3)
      b3(k,9) = b3(k,9) - amtx(k,3,10) * wk * amtx(k,9,3)
      wk =                b1(k,4) &
          + amtx(k,5,4) * b1(k,5) &
          + amtx(k,6,4) * b1(k,6) &
          + amtx(k,7,4) * b1(k,7) &
          + amtx(k,8,4) * b1(k,8) &
          + amtx(k,9,4) * b1(k,9)
      b1(k,4) = b1(k,4) - amtx(k,4,10) * wk
      b1(k,5) = b1(k,5) - amtx(k,4,10) * wk * amtx(k,5,4)
      b1(k,6) = b1(k,6) - amtx(k,4,10) * wk * amtx(k,6,4)
      b1(k,7) = b1(k,7) - amtx(k,4,10) * wk * amtx(k,7,4)
      b1(k,8) = b1(k,8) - amtx(k,4,10) * wk * amtx(k,8,4)
      b1(k,9) = b1(k,9) - amtx(k,4,10) * wk * amtx(k,9,4)
      wk =                b2(k,4) &
          + amtx(k,5,4) * b2(k,5) &
          + amtx(k,6,4) * b2(k,6) &
          + amtx(k,7,4) * b2(k,7) &
          + amtx(k,8,4) * b2(k,8) &
          + amtx(k,9,4) * b2(k,9)
      b2(k,4) = b2(k,4) - amtx(k,4,10) * wk
      b2(k,5) = b2(k,5) - amtx(k,4,10) * wk * amtx(k,5,4)
      b2(k,6) = b2(k,6) - amtx(k,4,10) * wk * amtx(k,6,4)
      b2(k,7) = b2(k,7) - amtx(k,4,10) * wk * amtx(k,7,4)
      b2(k,8) = b2(k,8) - amtx(k,4,10) * wk * amtx(k,8,4)
      b2(k,9) = b2(k,9) - amtx(k,4,10) * wk * amtx(k,9,4)
      wk =                b3(k,4) &
          + amtx(k,5,4) * b3(k,5) &
          + amtx(k,6,4) * b3(k,6) &
          + amtx(k,7,4) * b3(k,7) &
          + amtx(k,8,4) * b3(k,8) &
          + amtx(k,9,4) * b3(k,9)
      b3(k,4) = b3(k,4) - amtx(k,4,10) * wk
      b3(k,5) = b3(k,5) - amtx(k,4,10) * wk * amtx(k,5,4)
      b3(k,6) = b3(k,6) - amtx(k,4,10) * wk * amtx(k,6,4)
      b3(k,7) = b3(k,7) - amtx(k,4,10) * wk * amtx(k,7,4)
      b3(k,8) = b3(k,8) - amtx(k,4,10) * wk * amtx(k,8,4)
      b3(k,9) = b3(k,9) - amtx(k,4,10) * wk * amtx(k,9,4)
      wk =                b1(k,5) &
          + amtx(k,6,5) * b1(k,6) &
          + amtx(k,7,5) * b1(k,7) &
          + amtx(k,8,5) * b1(k,8) &
          + amtx(k,9,5) * b1(k,9)
      b1(k,5) = b1(k,5) - amtx(k,5,10) * wk
      b1(k,6) = b1(k,6) - amtx(k,5,10) * wk * amtx(k,6,5)
      b1(k,7) = b1(k,7) - amtx(k,5,10) * wk * amtx(k,7,5)
      b1(k,8) = b1(k,8) - amtx(k,5,10) * wk * amtx(k,8,5)
      b1(k,9) = b1(k,9) - amtx(k,5,10) * wk * amtx(k,9,5)
      wk =                b2(k,5) &
          + amtx(k,6,5) * b2(k,6) &
          + amtx(k,7,5) * b2(k,7) &
          + amtx(k,8,5) * b2(k,8) &
          + amtx(k,9,5) * b2(k,9)
      b2(k,5) = b2(k,5) - amtx(k,5,10) * wk
      b2(k,6) = b2(k,6) - amtx(k,5,10) * wk * amtx(k,6,5)
      b2(k,7) = b2(k,7) - amtx(k,5,10) * wk * amtx(k,7,5)
      b2(k,8) = b2(k,8) - amtx(k,5,10) * wk * amtx(k,8,5)
      b2(k,9) = b2(k,9) - amtx(k,5,10) * wk * amtx(k,9,5)
      wk =                b3(k,5) &
          + amtx(k,6,5) * b3(k,6) &
          + amtx(k,7,5) * b3(k,7) &
          + amtx(k,8,5) * b3(k,8) &
          + amtx(k,9,5) * b3(k,9)
      b3(k,5) = b3(k,5) - amtx(k,5,10) * wk
      b3(k,6) = b3(k,6) - amtx(k,5,10) * wk * amtx(k,6,5)
      b3(k,7) = b3(k,7) - amtx(k,5,10) * wk * amtx(k,7,5)
      b3(k,8) = b3(k,8) - amtx(k,5,10) * wk * amtx(k,8,5)
      b3(k,9) = b3(k,9) - amtx(k,5,10) * wk * amtx(k,9,5)
      wk =                b1(k,6) &
          + amtx(k,7,6) * b1(k,7) &
          + amtx(k,8,6) * b1(k,8) &
          + amtx(k,9,6) * b1(k,9)
      b1(k,6) = b1(k,6) - amtx(k,6,10) * wk
      b1(k,7) = b1(k,7) - amtx(k,6,10) * wk * amtx(k,7,6)
      b1(k,8) = b1(k,8) - amtx(k,6,10) * wk * amtx(k,8,6)
      b1(k,9) = b1(k,9) - amtx(k,6,10) * wk * amtx(k,9,6)
      wk =                b2(k,6) &
          + amtx(k,7,6) * b2(k,7) &
          + amtx(k,8,6) * b2(k,8) &
          + amtx(k,9,6) * b2(k,9)
      b2(k,6) = b2(k,6) - amtx(k,6,10) * wk
      b2(k,7) = b2(k,7) - amtx(k,6,10) * wk * amtx(k,7,6)
      b2(k,8) = b2(k,8) - amtx(k,6,10) * wk * amtx(k,8,6)
      b2(k,9) = b2(k,9) - amtx(k,6,10) * wk * amtx(k,9,6)
      wk =                b3(k,6) &
          + amtx(k,7,6) * b3(k,7) &
          + amtx(k,8,6) * b3(k,8) &
          + amtx(k,9,6) * b3(k,9)
      b3(k,6) = b3(k,6) - amtx(k,6,10) * wk
      b3(k,7) = b3(k,7) - amtx(k,6,10) * wk * amtx(k,7,6)
      b3(k,8) = b3(k,8) - amtx(k,6,10) * wk * amtx(k,8,6)
      b3(k,9) = b3(k,9) - amtx(k,6,10) * wk * amtx(k,9,6)
      wk =                b1(k,7) &
          + amtx(k,8,7) * b1(k,8) &
          + amtx(k,9,7) * b1(k,9)
      b1(k,7) = b1(k,7) - amtx(k,7,10) * wk
      b1(k,8) = b1(k,8) - amtx(k,7,10) * wk * amtx(k,8,7)
      b1(k,9) = b1(k,9) - amtx(k,7,10) * wk * amtx(k,9,7)
      wk =                b2(k,7) &
          + amtx(k,8,7) * b2(k,8) &
          + amtx(k,9,7) * b2(k,9)
      b2(k,7) = b2(k,7) - amtx(k,7,10) * wk
      b2(k,8) = b2(k,8) - amtx(k,7,10) * wk * amtx(k,8,7)
      b2(k,9) = b2(k,9) - amtx(k,7,10) * wk * amtx(k,9,7)
      wk =                b3(k,7) &
          + amtx(k,8,7) * b3(k,8) &
          + amtx(k,9,7) * b3(k,9)
      b3(k,7) = b3(k,7) - amtx(k,7,10) * wk
      b3(k,8) = b3(k,8) - amtx(k,7,10) * wk * amtx(k,8,7)
      b3(k,9) = b3(k,9) - amtx(k,7,10) * wk * amtx(k,9,7)
      wk =                b1(k,8) &
          + amtx(k,9,8) * b1(k,9)
      b1(k,8) = b1(k,8) - amtx(k,8,10) * wk
      b1(k,9) = b1(k,9) - amtx(k,8,10) * wk * amtx(k,9,8)
      wk =                b2(k,8) &
          + amtx(k,9,8) * b2(k,9)
      b2(k,8) = b2(k,8) - amtx(k,8,10) * wk
      b2(k,9) = b2(k,9) - amtx(k,8,10) * wk * amtx(k,9,8)
      wk =                b3(k,8) &
          + amtx(k,9,8) * b3(k,9)
      b3(k,8) = b3(k,8) - amtx(k,8,10) * wk
      b3(k,9) = b3(k,9) - amtx(k,8,10) * wk * amtx(k,9,8)
      wk =                b1(k,9)
      b1(k,9) = b1(k,9) - amtx(k,9,10) * wk
      wk =                b2(k,9)
      b2(k,9) = b2(k,9) - amtx(k,9,10) * wk
      wk =                b3(k,9)
      b3(k,9) = b3(k,9) - amtx(k,9,10) * wk

      b1(k,9) = b1(k,9) _OP_ amtx(k,9,9)
        b1(k,1) = b1(k,1) - b1(k,9)*amtx(k,1,9)
        b1(k,2) = b1(k,2) - b1(k,9)*amtx(k,2,9)
        b1(k,3) = b1(k,3) - b1(k,9)*amtx(k,3,9)
        b1(k,4) = b1(k,4) - b1(k,9)*amtx(k,4,9)
        b1(k,5) = b1(k,5) - b1(k,9)*amtx(k,5,9)
        b1(k,6) = b1(k,6) - b1(k,9)*amtx(k,6,9)
        b1(k,7) = b1(k,7) - b1(k,9)*amtx(k,7,9)
        b1(k,8) = b1(k,8) - b1(k,9)*amtx(k,8,9)
      b2(k,9) = b2(k,9) _OP_ amtx(k,9,9)
        b2(k,1) = b2(k,1) - b2(k,9)*amtx(k,1,9)
        b2(k,2) = b2(k,2) - b2(k,9)*amtx(k,2,9)
        b2(k,3) = b2(k,3) - b2(k,9)*amtx(k,3,9)
        b2(k,4) = b2(k,4) - b2(k,9)*amtx(k,4,9)
        b2(k,5) = b2(k,5) - b2(k,9)*amtx(k,5,9)
        b2(k,6) = b2(k,6) - b2(k,9)*amtx(k,6,9)
        b2(k,7) = b2(k,7) - b2(k,9)*amtx(k,7,9)
        b2(k,8) = b2(k,8) - b2(k,9)*amtx(k,8,9)
      b3(k,9) = b3(k,9) _OP_ amtx(k,9,9)
        b3(k,1) = b3(k,1) - b3(k,9)*amtx(k,1,9)
        b3(k,2) = b3(k,2) - b3(k,9)*amtx(k,2,9)
        b3(k,3) = b3(k,3) - b3(k,9)*amtx(k,3,9)
        b3(k,4) = b3(k,4) - b3(k,9)*amtx(k,4,9)
        b3(k,5) = b3(k,5) - b3(k,9)*amtx(k,5,9)
        b3(k,6) = b3(k,6) - b3(k,9)*amtx(k,6,9)
        b3(k,7) = b3(k,7) - b3(k,9)*amtx(k,7,9)
        b3(k,8) = b3(k,8) - b3(k,9)*amtx(k,8,9)
      b1(k,8) = b1(k,8) _OP_ amtx(k,8,8)
        b1(k,1) = b1(k,1) - b1(k,8)*amtx(k,1,8)
        b1(k,2) = b1(k,2) - b1(k,8)*amtx(k,2,8)
        b1(k,3) = b1(k,3) - b1(k,8)*amtx(k,3,8)
        b1(k,4) = b1(k,4) - b1(k,8)*amtx(k,4,8)
        b1(k,5) = b1(k,5) - b1(k,8)*amtx(k,5,8)
        b1(k,6) = b1(k,6) - b1(k,8)*amtx(k,6,8)
        b1(k,7) = b1(k,7) - b1(k,8)*amtx(k,7,8)
      b2(k,8) = b2(k,8) _OP_ amtx(k,8,8)
        b2(k,1) = b2(k,1) - b2(k,8)*amtx(k,1,8)
        b2(k,2) = b2(k,2) - b2(k,8)*amtx(k,2,8)
        b2(k,3) = b2(k,3) - b2(k,8)*amtx(k,3,8)
        b2(k,4) = b2(k,4) - b2(k,8)*amtx(k,4,8)
        b2(k,5) = b2(k,5) - b2(k,8)*amtx(k,5,8)
        b2(k,6) = b2(k,6) - b2(k,8)*amtx(k,6,8)
        b2(k,7) = b2(k,7) - b2(k,8)*amtx(k,7,8)
      b3(k,8) = b3(k,8) _OP_ amtx(k,8,8)
        b3(k,1) = b3(k,1) - b3(k,8)*amtx(k,1,8)
        b3(k,2) = b3(k,2) - b3(k,8)*amtx(k,2,8)
        b3(k,3) = b3(k,3) - b3(k,8)*amtx(k,3,8)
        b3(k,4) = b3(k,4) - b3(k,8)*amtx(k,4,8)
        b3(k,5) = b3(k,5) - b3(k,8)*amtx(k,5,8)
        b3(k,6) = b3(k,6) - b3(k,8)*amtx(k,6,8)
        b3(k,7) = b3(k,7) - b3(k,8)*amtx(k,7,8)
      b1(k,7) = b1(k,7) _OP_ amtx(k,7,7)
        b1(k,1) = b1(k,1) - b1(k,7)*amtx(k,1,7)
        b1(k,2) = b1(k,2) - b1(k,7)*amtx(k,2,7)
        b1(k,3) = b1(k,3) - b1(k,7)*amtx(k,3,7)
        b1(k,4) = b1(k,4) - b1(k,7)*amtx(k,4,7)
        b1(k,5) = b1(k,5) - b1(k,7)*amtx(k,5,7)
        b1(k,6) = b1(k,6) - b1(k,7)*amtx(k,6,7)
      b2(k,7) = b2(k,7) _OP_ amtx(k,7,7)
        b2(k,1) = b2(k,1) - b2(k,7)*amtx(k,1,7)
        b2(k,2) = b2(k,2) - b2(k,7)*amtx(k,2,7)
        b2(k,3) = b2(k,3) - b2(k,7)*amtx(k,3,7)
        b2(k,4) = b2(k,4) - b2(k,7)*amtx(k,4,7)
        b2(k,5) = b2(k,5) - b2(k,7)*amtx(k,5,7)
        b2(k,6) = b2(k,6) - b2(k,7)*amtx(k,6,7)
      b3(k,7) = b3(k,7) _OP_ amtx(k,7,7)
        b3(k,1) = b3(k,1) - b3(k,7)*amtx(k,1,7)
        b3(k,2) = b3(k,2) - b3(k,7)*amtx(k,2,7)
        b3(k,3) = b3(k,3) - b3(k,7)*amtx(k,3,7)
        b3(k,4) = b3(k,4) - b3(k,7)*amtx(k,4,7)
        b3(k,5) = b3(k,5) - b3(k,7)*amtx(k,5,7)
        b3(k,6) = b3(k,6) - b3(k,7)*amtx(k,6,7)
      b1(k,6) = b1(k,6) _OP_ amtx(k,6,6)
        b1(k,1) = b1(k,1) - b1(k,6)*amtx(k,1,6)
        b1(k,2) = b1(k,2) - b1(k,6)*amtx(k,2,6)
        b1(k,3) = b1(k,3) - b1(k,6)*amtx(k,3,6)
        b1(k,4) = b1(k,4) - b1(k,6)*amtx(k,4,6)
        b1(k,5) = b1(k,5) - b1(k,6)*amtx(k,5,6)
      b2(k,6) = b2(k,6) _OP_ amtx(k,6,6)
        b2(k,1) = b2(k,1) - b2(k,6)*amtx(k,1,6)
        b2(k,2) = b2(k,2) - b2(k,6)*amtx(k,2,6)
        b2(k,3) = b2(k,3) - b2(k,6)*amtx(k,3,6)
        b2(k,4) = b2(k,4) - b2(k,6)*amtx(k,4,6)
        b2(k,5) = b2(k,5) - b2(k,6)*amtx(k,5,6)
      b3(k,6) = b3(k,6) _OP_ amtx(k,6,6)
        b3(k,1) = b3(k,1) - b3(k,6)*amtx(k,1,6)
        b3(k,2) = b3(k,2) - b3(k,6)*amtx(k,2,6)
        b3(k,3) = b3(k,3) - b3(k,6)*amtx(k,3,6)
        b3(k,4) = b3(k,4) - b3(k,6)*amtx(k,4,6)
        b3(k,5) = b3(k,5) - b3(k,6)*amtx(k,5,6)
      b1(k,5) = b1(k,5) _OP_ amtx(k,5,5)
        b1(k,1) = b1(k,1) - b1(k,5)*amtx(k,1,5)
        b1(k,2) = b1(k,2) - b1(k,5)*amtx(k,2,5)
        b1(k,3) = b1(k,3) - b1(k,5)*amtx(k,3,5)
        b1(k,4) = b1(k,4) - b1(k,5)*amtx(k,4,5)
      b2(k,5) = b2(k,5) _OP_ amtx(k,5,5)
        b2(k,1) = b2(k,1) - b2(k,5)*amtx(k,1,5)
        b2(k,2) = b2(k,2) - b2(k,5)*amtx(k,2,5)
        b2(k,3) = b2(k,3) - b2(k,5)*amtx(k,3,5)
        b2(k,4) = b2(k,4) - b2(k,5)*amtx(k,4,5)
      b3(k,5) = b3(k,5) _OP_ amtx(k,5,5)
        b3(k,1) = b3(k,1) - b3(k,5)*amtx(k,1,5)
        b3(k,2) = b3(k,2) - b3(k,5)*amtx(k,2,5)
        b3(k,3) = b3(k,3) - b3(k,5)*amtx(k,3,5)
        b3(k,4) = b3(k,4) - b3(k,5)*amtx(k,4,5)
      b1(k,4) = b1(k,4) _OP_ amtx(k,4,4)
        b1(k,1) = b1(k,1) - b1(k,4)*amtx(k,1,4)
        b1(k,2) = b1(k,2) - b1(k,4)*amtx(k,2,4)
        b1(k,3) = b1(k,3) - b1(k,4)*amtx(k,3,4)
      b2(k,4) = b2(k,4) _OP_ amtx(k,4,4)
        b2(k,1) = b2(k,1) - b2(k,4)*amtx(k,1,4)
        b2(k,2) = b2(k,2) - b2(k,4)*amtx(k,2,4)
        b2(k,3) = b2(k,3) - b2(k,4)*amtx(k,3,4)
      b3(k,4) = b3(k,4) _OP_ amtx(k,4,4)
        b3(k,1) = b3(k,1) - b3(k,4)*amtx(k,1,4)
        b3(k,2) = b3(k,2) - b3(k,4)*amtx(k,2,4)
        b3(k,3) = b3(k,3) - b3(k,4)*amtx(k,3,4)
      b1(k,3) = b1(k,3) _OP_ amtx(k,3,3)
        b1(k,1) = b1(k,1) - b1(k,3)*amtx(k,1,3)
        b1(k,2) = b1(k,2) - b1(k,3)*amtx(k,2,3)
      b2(k,3) = b2(k,3) _OP_ amtx(k,3,3)
        b2(k,1) = b2(k,1) - b2(k,3)*amtx(k,1,3)
        b2(k,2) = b2(k,2) - b2(k,3)*amtx(k,2,3)
      b3(k,3) = b3(k,3) _OP_ amtx(k,3,3)
        b3(k,1) = b3(k,1) - b3(k,3)*amtx(k,1,3)
        b3(k,2) = b3(k,2) - b3(k,3)*amtx(k,2,3)
      b1(k,2) = b1(k,2) _OP_ amtx(k,2,2)
        b1(k,1) = b1(k,1) - b1(k,2)*amtx(k,1,2)
      b2(k,2) = b2(k,2) _OP_ amtx(k,2,2)
        b2(k,1) = b2(k,1) - b2(k,2)*amtx(k,1,2)
      b3(k,2) = b3(k,2) _OP_ amtx(k,2,2)
        b3(k,1) = b3(k,1) - b3(k,2)*amtx(k,1,2)
      b1(k,1) = b1(k,1) _OP_ amtx(k,1,1)
      b2(k,1) = b2(k,1) _OP_ amtx(k,1,1)
      b3(k,1) = b3(k,1) _OP_ amtx(k,1,1)
  enddo !k-loop
#ifdef FINEGRAINED_TIMING
  ret = gptlstop_handle ('solveiThLS3', handle)
#endif
END SUBROUTINE solveiThLS3

!**************************************************************************
!
!       Nonhydrostatic Icosahedral Model (NIM)
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)!
!
!**************************************************************************
