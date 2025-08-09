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


SUBROUTINE invert_A_diag(nz,nob,ims,ime,ips,ipe,amtx1,amtx2)
!**************************************************************************
!
!       Invert diagonal elements of "A matrix" variants to avoid 
!       division in solveiThLS*() routines as an optimization.  
!
!**************************************************************************
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nz,nob,ims,ime,ips,ipe
  REAL   , INTENT(INOUT) :: amtx1(nz  ,nob,nob+1,ims:ime)
  REAL   , INTENT(INOUT) :: amtx2(nz-1,nob,nob+1,ims:ime)
  INTEGER                :: k,ipn

print *,'Inverting diagonal elements of amtx1 and amtx2...'  
call flush(6)

!$OMP PARALLEL DO SCHEDULE(runtime)
do ipn=ips,ipe
  do k=1,nz
      amtx1(k,9,9,ipn) = 1.0 / amtx1(k,9,9,ipn)
      amtx1(k,8,8,ipn) = 1.0 / amtx1(k,8,8,ipn)
      amtx1(k,7,7,ipn) = 1.0 / amtx1(k,7,7,ipn)
      amtx1(k,6,6,ipn) = 1.0 / amtx1(k,6,6,ipn)
      amtx1(k,5,5,ipn) = 1.0 / amtx1(k,5,5,ipn)
      amtx1(k,4,4,ipn) = 1.0 / amtx1(k,4,4,ipn)
      amtx1(k,3,3,ipn) = 1.0 / amtx1(k,3,3,ipn)
      amtx1(k,2,2,ipn) = 1.0 / amtx1(k,2,2,ipn)
      amtx1(k,1,1,ipn) = 1.0 / amtx1(k,1,1,ipn)
  enddo !k-loop
  do k=1,nz-1
      amtx2(k,9,9,ipn) = 1.0 / amtx2(k,9,9,ipn)
      amtx2(k,8,8,ipn) = 1.0 / amtx2(k,8,8,ipn)
      amtx2(k,7,7,ipn) = 1.0 / amtx2(k,7,7,ipn)
      amtx2(k,6,6,ipn) = 1.0 / amtx2(k,6,6,ipn)
      amtx2(k,5,5,ipn) = 1.0 / amtx2(k,5,5,ipn)
      amtx2(k,4,4,ipn) = 1.0 / amtx2(k,4,4,ipn)
      amtx2(k,3,3,ipn) = 1.0 / amtx2(k,3,3,ipn)
      amtx2(k,2,2,ipn) = 1.0 / amtx2(k,2,2,ipn)
      amtx2(k,1,1,ipn) = 1.0 / amtx2(k,1,1,ipn)
  enddo !k-loop
enddo ! ipn loop
!$OMP END PARALLEL DO
END SUBROUTINE invert_A_diag

