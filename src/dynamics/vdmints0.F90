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

subroutine vdmints0(w,sedgvar,bedgvar,nvars,nvi,amtx2,nbf,npp,ims,ime,nob,&
                   nprox,prox,xyc,zc,z,ca4k,ca4p,ips,ipe)
!**********************************************************************
!
!       Three-dimensional interpolation to obtain edge values 
!       from those at centers
!       Jin Lee                  February, 2009
!       Jacques Middlecoff July 2009 Put K on the inside for the GPU
!       Jacques Middlecoff July 2009 Changed sedgvar, bedgvar to 4D to stop copy
!       Jacques Middlecoff July 2009 Put variables in calling sequence.
!
!**********************************************************************

use kinds, only: rt
use ReadNamelist, only: nz
use gptl
use gptl_acc
implicit none
integer,intent (IN)    :: nvi,npp,ims,ime,nvars,nob,nbf,ips,ipe
integer,intent (IN)    :: nprox(           ims:ime) ! holds number of proximity points
integer,intent (IN)    :: prox (       npp,ims:ime) ! holds index of proximity points
real(rt),intent (IN)   :: xyc  (     2,npp,ims:ime)
real(rt),intent (IN)   :: zc   (  0:NZ,npp,ims:ime)
real(rt),intent (IN)   :: z    (  0:NZ,    ims:ime)
real(rt),intent (IN)   :: ca4k (    nz            )
real(rt),intent (IN)   :: ca4p (    nz            )
real(rt),intent (IN)   :: w    (  0:NZ,    ims:ime)
real(rt),intent (INOUT):: sedgvar (  1:NZ,npp,nvars,ims:ime)
real(rt),intent (INOUT):: bedgvar (  0:NZ,    ims:ime,nvars)
real(rt),intent (IN)   :: amtx2(0:nz-1,(nob+1)*nbf,ims:ime)
real(rt)               :: tgtc
real(rt)               :: rhs  (NZ,9)
real(rt)               :: Tgt  (NZ,npp)
INTEGER                :: ipn,k,isn,isp,ism,i
real(rt)               :: xyc1, xyc2  ! xyc temporaries

integer                :: mythread  ! thread number when FINEGRAINED_TIMING defined
integer                :: ret       ! GPTL return code when FINEGRAINED_TIMING defined
logical, save          :: first = .true.
integer, save          :: vdmints0_handle, ipn_handle, solvei_handle
integer                :: startk

#ifdef FINEGRAINED_TIMING
#include <gptl.inc>
  integer, external :: omp_get_thread_num
  ret = gptlstart ('vdmints0_outside')
#endif

!$acc routine(solveiThLS0) vector

  if (first) then
    first = .false.
!$acc parallel private(ret) copyout(vdmints0_handle, ipn_handle, solvei_handle)
    ret = gptlinit_handle_gpu ('vdmints0', vdmints0_handle)
    ret = gptlinit_handle_gpu ('vdmints0_ipn', ipn_handle)
    ret = gptlinit_handle_gpu ('vdmints0_solvei', solvei_handle)
!$acc end parallel
  end if

!$acc parallel private(ret) copyin(vdmints0_handle)
  ret = gptlstart_gpu (vdmints0_handle)
!$acc end parallel

!DIR$ ASSUME_ALIGNED ca4k:64, ca4p:64, sedgvar:64
!$acc parallel private(ret) num_workers(PAR_WRK) vector_length(VEC_LEN) copyin(ipn_handle,solvei_handle)
!$acc loop gang private(rhs,Tgt)
!$OMP PARALLEL DO PRIVATE(mythread,ret,k,rhs,i,isn,tgtc,tgt,isp,ism,xyc1,xyc2) SCHEDULE(runtime)
do ipn=ips,ipe
  ret = gptlstart_gpu (ipn_handle)
!$acc cache(rhs,tgt)
#ifdef FINEGRAINED_TIMING
  mythread = omp_get_thread_num ()
  if (mythread == 0) then
    ret = gptlstart_handle ('vdmints0_inside', handle)
  end if
#endif
!$acc loop worker vector
  do k=2,NZ
    rhs(k,1) = w(k-1,prox(1         ,ipn)) - w(k-1,ipn)
    rhs(k,2) = w(k-1,prox(2         ,ipn)) - w(k-1,ipn)
    rhs(k,3) = w(k-1,prox(3         ,ipn)) - w(k-1,ipn)
    rhs(k,4) = w(k-1,prox(5         ,ipn)) - w(k-1,ipn)
    rhs(k,5) = w(k  ,prox(2         ,ipn)) - w(k-1,ipn)
    rhs(k,6) = w(k  ,prox(3         ,ipn)) - w(k-1,ipn)
    rhs(k,7) = w(k  ,prox(4         ,ipn)) - w(k-1,ipn)
    rhs(k,8) = w(k  ,prox(nprox(ipn),ipn)) - w(k-1,ipn)
    rhs(k,9) = .5_rt*w(k-1,ipn) + .5_rt*w(k,ipn) - w(k-1,ipn)
  enddo !k-loop
!$acc loop vector
  do i=1,nob
    rhs(1,i)=0._rt
  enddo

  ret = gptlstart_gpu (solvei_handle)
!$acc loop worker
  do startk=1,NZ-1,32
    CALL solveiThLS0(nob,nbf,rhs,amtx2(0,1,ipn),ipn,ips,startk)
  end do
  ret = gptlstop_gpu (solvei_handle)

!JR Defining temporary variables xyc1 and xyc2 prevents ifort from swapping the k and isn loops below
!$acc loop worker private(xyc1,xyc2)
  do isn = 1,nprox(ipn)
    xyc1 = xyc(1,isn,ipn)
    xyc2 = xyc(2,isn,ipn)
!$acc loop vector private(tgtc)
    do k=2,NZ ! Tgt(nz) is calculated for ks0=0
      tgtc = ( (.5_rt*(zc(k-1,isn,ipn) + zc(k,isn,ipn)))-z(k-1,ipn) )
      Tgt(k-1,isn) = xyc1*( rhs(k,1)*xyc1  &
                         +rhs(k,4)*xyc2  &
                         +rhs(k,5)*tgtc+rhs(k,7)) &
                 + xyc2*( rhs(k,2)*xyc2  &
                         +rhs(k,6)*tgtc+rhs(k,8)) &
                 + tgtc*( rhs(k,3)*tgtc+rhs(k,9)) &
                 + w(k-1,ipn)
    enddo !k-loop
  end do  ! isn-loop

!$acc loop worker private(isp,ism)
  do isn = 1,nprox(ipn)
     isp=mod(isn,nprox(ipn))+1
     ism=isn-1
     if(ism.eq.0) ism=nprox(ipn)
!$acc loop vector
     do k = 2,NZ-1
       sedgvar(k,isn,nvi,ipn) = .5_rt*( ca4k(k-1)*(Tgt(k-1,isn)+Tgt(k-1 ,isp)) &
                               +     ca4p(k-1)*(Tgt(k  ,isn)+Tgt(k   ,isp))  )
     end do !  k -loop
     sedgvar( 1,isn,nvi,ipn) = (z(1,ipn)-z(3,ipn))/(z(2,ipn)-z(3,ipn))*sedgvar(2,isn,nvi,ipn) &
                              +(z(1,ipn)-z(2,ipn))/(z(3,ipn)-z(2,ipn))*sedgvar(3,isn,nvi,ipn) 
     sedgvar(NZ,isn,nvi,ipn) = 2._rt*sedgvar(nz-1,isn,nvi,ipn)-sedgvar(nz-2,isn,nvi,ipn)
  end do  ! isn-loop

!$acc loop vector
  do k=0,NZ-1
    bedgvar(k,ipn,nvi) = .5_rt*w(k,ipn) + .5_rt*w(k+1,ipn)
  end do
  bedgvar(NZ,ipn,nvi) = .5_rt*w(NZ,ipn) + .5_rt*w(NZ,ipn)
#ifdef FINEGRAINED_TIMING
  if (mythread == 0) then
    ret = gptlstop_handle ('vdmints0_inside', handle)
  end if
#endif
  ret = gptlstop_gpu (ipn_handle)
enddo !ipn-loop
!$acc end parallel

!$OMP END PARALLEL DO
#ifdef FINEGRAINED_TIMING
ret = gptlstop ('vdmints0_outside')
#endif

!$acc parallel private(ret) copyin(vdmints0_handle)
ret = gptlstop_gpu (vdmints0_handle)
!$acc end parallel

return
end subroutine vdmints0

!**********************************************************************
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
!**********************************************************************
#include "solveiThLS0.F90"
