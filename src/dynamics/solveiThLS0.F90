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

SUBROUTINE solveiThLS0(nob,nbf,b,amtx,ipn,ips,startk)
!**************************************************************************
!
!       This file contains the classes that implement the least square fitting. 
!       using spline wavelet as basis functions
!       Ning Wang, Feb.. 2009
!       Jacques Middlecoff April 2009 optimized
!       Jacques Middlecoff July 2009 Put K on the inside for the GPU
!       Tom Hendeson       August 2010 bitwise-exact optimizations for GPU
!
!**************************************************************************
use kinds, only: rt
!MWG workaround due to an F2C generation error
!use ReadNamelist, only: nz
IMPLICIT NONE
! MWG: F2C workaround, declare nz explicitly
INTEGER,PARAMETER      ::  nz = 96

INTEGER, INTENT(IN)    :: nbf,nob,ipn,ips
real(rt), INTENT(INOUT):: b   (    NZ      ,nob      )
real(rt), INTENT(IN)   :: amtx(0:NZ-1,nob,nob+1)
integer, intent(in)    :: startk
real(rt)               :: wk
INTEGER                :: k

!DIR$ ASSUME_ALIGNED b:64

#ifdef FINEGRAINED_TIMING
#include <gptl.inc>
  integer :: ret
  integer :: handle = 0
  ret = gptlstart_handle ('solveiThLS0', handle)
#endif
!$acc routine(solveiThLS0) vector

!JR Need to verify that this vectorizes
!$acc loop vector private(wk)
  do k=startk,min(NZ-1,startk+32-1)
      wk =                b(k,1) &
          + amtx(k,2,1) * b(k,2) &
          + amtx(k,3,1) * b(k,3) &
          + amtx(k,4,1) * b(k,4) &
          + amtx(k,5,1) * b(k,5) &
          + amtx(k,6,1) * b(k,6) &
          + amtx(k,7,1) * b(k,7) &
          + amtx(k,8,1) * b(k,8) &
          + amtx(k,9,1) * b(k,9)
      b(k,1) = b(k,1) - amtx(k,1,10) * wk
      b(k,2) = b(k,2) - amtx(k,1,10) * wk * amtx(k,2,1)
      b(k,3) = b(k,3) - amtx(k,1,10) * wk * amtx(k,3,1)
      b(k,4) = b(k,4) - amtx(k,1,10) * wk * amtx(k,4,1)
      b(k,5) = b(k,5) - amtx(k,1,10) * wk * amtx(k,5,1)
      b(k,6) = b(k,6) - amtx(k,1,10) * wk * amtx(k,6,1)
      b(k,7) = b(k,7) - amtx(k,1,10) * wk * amtx(k,7,1)
      b(k,8) = b(k,8) - amtx(k,1,10) * wk * amtx(k,8,1)
      b(k,9) = b(k,9) - amtx(k,1,10) * wk * amtx(k,9,1)
      wk =                b(k,2) &
          + amtx(k,3,2) * b(k,3) &
          + amtx(k,4,2) * b(k,4) &
          + amtx(k,5,2) * b(k,5) &
          + amtx(k,6,2) * b(k,6) &
          + amtx(k,7,2) * b(k,7) &
          + amtx(k,8,2) * b(k,8) &
          + amtx(k,9,2) * b(k,9)
      b(k,2) = b(k,2) - amtx(k,2,10) * wk
      b(k,3) = b(k,3) - amtx(k,2,10) * wk * amtx(k,3,2)
      b(k,4) = b(k,4) - amtx(k,2,10) * wk * amtx(k,4,2)
      b(k,5) = b(k,5) - amtx(k,2,10) * wk * amtx(k,5,2)
      b(k,6) = b(k,6) - amtx(k,2,10) * wk * amtx(k,6,2)
      b(k,7) = b(k,7) - amtx(k,2,10) * wk * amtx(k,7,2)
      b(k,8) = b(k,8) - amtx(k,2,10) * wk * amtx(k,8,2)
      b(k,9) = b(k,9) - amtx(k,2,10) * wk * amtx(k,9,2)
      wk =                b(k,3) &
          + amtx(k,4,3) * b(k,4) &
          + amtx(k,5,3) * b(k,5) &
          + amtx(k,6,3) * b(k,6) &
          + amtx(k,7,3) * b(k,7) &
          + amtx(k,8,3) * b(k,8) &
          + amtx(k,9,3) * b(k,9)
      b(k,3) = b(k,3) - amtx(k,3,10) * wk
      b(k,4) = b(k,4) - amtx(k,3,10) * wk * amtx(k,4,3)
      b(k,5) = b(k,5) - amtx(k,3,10) * wk * amtx(k,5,3)
      b(k,6) = b(k,6) - amtx(k,3,10) * wk * amtx(k,6,3)
      b(k,7) = b(k,7) - amtx(k,3,10) * wk * amtx(k,7,3)
      b(k,8) = b(k,8) - amtx(k,3,10) * wk * amtx(k,8,3)
      b(k,9) = b(k,9) - amtx(k,3,10) * wk * amtx(k,9,3)
      wk =                b(k,4) &
          + amtx(k,5,4) * b(k,5) &
          + amtx(k,6,4) * b(k,6) &
          + amtx(k,7,4) * b(k,7) &
          + amtx(k,8,4) * b(k,8) &
          + amtx(k,9,4) * b(k,9)
      b(k,4) = b(k,4) - amtx(k,4,10) * wk
      b(k,5) = b(k,5) - amtx(k,4,10) * wk * amtx(k,5,4)
      b(k,6) = b(k,6) - amtx(k,4,10) * wk * amtx(k,6,4)
      b(k,7) = b(k,7) - amtx(k,4,10) * wk * amtx(k,7,4)
      b(k,8) = b(k,8) - amtx(k,4,10) * wk * amtx(k,8,4)
      b(k,9) = b(k,9) - amtx(k,4,10) * wk * amtx(k,9,4)
      wk =                b(k,5) &
          + amtx(k,6,5) * b(k,6) &
          + amtx(k,7,5) * b(k,7) &
          + amtx(k,8,5) * b(k,8) &
          + amtx(k,9,5) * b(k,9)
      b(k,5) = b(k,5) - amtx(k,5,10) * wk
      b(k,6) = b(k,6) - amtx(k,5,10) * wk * amtx(k,6,5)
      b(k,7) = b(k,7) - amtx(k,5,10) * wk * amtx(k,7,5)
      b(k,8) = b(k,8) - amtx(k,5,10) * wk * amtx(k,8,5)
      b(k,9) = b(k,9) - amtx(k,5,10) * wk * amtx(k,9,5)
      wk =                b(k,6) &
          + amtx(k,7,6) * b(k,7) &
          + amtx(k,8,6) * b(k,8) &
          + amtx(k,9,6) * b(k,9)
      b(k,6) = b(k,6) - amtx(k,6,10) * wk
      b(k,7) = b(k,7) - amtx(k,6,10) * wk * amtx(k,7,6)
      b(k,8) = b(k,8) - amtx(k,6,10) * wk * amtx(k,8,6)
      b(k,9) = b(k,9) - amtx(k,6,10) * wk * amtx(k,9,6)
      wk =                b(k,7) &
          + amtx(k,8,7) * b(k,8) &
          + amtx(k,9,7) * b(k,9)
      b(k,7) = b(k,7) - amtx(k,7,10) * wk
      b(k,8) = b(k,8) - amtx(k,7,10) * wk * amtx(k,8,7)
      b(k,9) = b(k,9) - amtx(k,7,10) * wk * amtx(k,9,7)
      wk =                b(k,8) &
          + amtx(k,9,8) * b(k,9)
      b(k,8) = b(k,8) - amtx(k,8,10) * wk
      b(k,9) = b(k,9) - amtx(k,8,10) * wk * amtx(k,9,8)
      wk =                b(k,9)
      b(k,9) = b(k,9) - amtx(k,9,10) * wk

      b(k,9) = b(k,9) _OP_ amtx(k,9,9)
        b(k,1) = b(k,1) - b(k,9)*amtx(k,1,9)
        b(k,2) = b(k,2) - b(k,9)*amtx(k,2,9)
        b(k,3) = b(k,3) - b(k,9)*amtx(k,3,9)
        b(k,4) = b(k,4) - b(k,9)*amtx(k,4,9)
        b(k,5) = b(k,5) - b(k,9)*amtx(k,5,9)
        b(k,6) = b(k,6) - b(k,9)*amtx(k,6,9)
        b(k,7) = b(k,7) - b(k,9)*amtx(k,7,9)
        b(k,8) = b(k,8) - b(k,9)*amtx(k,8,9)
      b(k,8) = b(k,8) _OP_ amtx(k,8,8)
        b(k,1) = b(k,1) - b(k,8)*amtx(k,1,8)
        b(k,2) = b(k,2) - b(k,8)*amtx(k,2,8)
        b(k,3) = b(k,3) - b(k,8)*amtx(k,3,8)
        b(k,4) = b(k,4) - b(k,8)*amtx(k,4,8)
        b(k,5) = b(k,5) - b(k,8)*amtx(k,5,8)
        b(k,6) = b(k,6) - b(k,8)*amtx(k,6,8)
        b(k,7) = b(k,7) - b(k,8)*amtx(k,7,8)
      b(k,7) = b(k,7) _OP_ amtx(k,7,7)
        b(k,1) = b(k,1) - b(k,7)*amtx(k,1,7)
        b(k,2) = b(k,2) - b(k,7)*amtx(k,2,7)
        b(k,3) = b(k,3) - b(k,7)*amtx(k,3,7)
        b(k,4) = b(k,4) - b(k,7)*amtx(k,4,7)
        b(k,5) = b(k,5) - b(k,7)*amtx(k,5,7)
        b(k,6) = b(k,6) - b(k,7)*amtx(k,6,7)
      b(k,6) = b(k,6) _OP_ amtx(k,6,6)
        b(k,1) = b(k,1) - b(k,6)*amtx(k,1,6)
        b(k,2) = b(k,2) - b(k,6)*amtx(k,2,6)
        b(k,3) = b(k,3) - b(k,6)*amtx(k,3,6)
        b(k,4) = b(k,4) - b(k,6)*amtx(k,4,6)
        b(k,5) = b(k,5) - b(k,6)*amtx(k,5,6)
      b(k,5) = b(k,5) _OP_ amtx(k,5,5)
        b(k,1) = b(k,1) - b(k,5)*amtx(k,1,5)
        b(k,2) = b(k,2) - b(k,5)*amtx(k,2,5)
        b(k,3) = b(k,3) - b(k,5)*amtx(k,3,5)
        b(k,4) = b(k,4) - b(k,5)*amtx(k,4,5)
      b(k,4) = b(k,4) _OP_ amtx(k,4,4)
        b(k,1) = b(k,1) - b(k,4)*amtx(k,1,4)
        b(k,2) = b(k,2) - b(k,4)*amtx(k,2,4)
        b(k,3) = b(k,3) - b(k,4)*amtx(k,3,4)
      b(k,3) = b(k,3) _OP_ amtx(k,3,3)
        b(k,1) = b(k,1) - b(k,3)*amtx(k,1,3)
        b(k,2) = b(k,2) - b(k,3)*amtx(k,2,3)
      b(k,2) = b(k,2) _OP_ amtx(k,2,2)
        b(k,1) = b(k,1) - b(k,2)*amtx(k,1,2)
      b(k,1) = b(k,1) _OP_ amtx(k,1,1)
  enddo !k-loop
#ifdef FINEGRAINED_TIMING
  ret = gptlstop ('solveiThLS0', handle)
#endif
END SUBROUTINE solveiThLS0

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
