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

subroutine vdmintv(u,v,sedgvar,bedgvar,nvars,amtx1,nbf,npp,ims,ime,nob,&
                   nprox,prox,cs,sn,xyc,zc,zm,ca4k,ca4p,z,ips,ipe)

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

use gptl
#ifdef _OPENACC
use gptl_acc
#endif
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent (IN)    :: npp,ims,ime,nvars,nob,nbf,ips,ipe
integer,intent (IN)    :: nprox  (         ims:ime) ! holds number of proximity points
integer,intent (IN)    :: prox   (     npp,ims:ime) ! holds index of proximity points
real(rt),intent (IN)   :: cs     (     npp,ims:ime) ! cosine transform function
real(rt),intent (IN)   :: sn     (     npp,ims:ime) ! sine transform function
real(rt),intent (IN)   :: xyc    (   2,npp,ims:ime)
real(rt),intent (IN)   :: zc     (0:NZ,npp,ims:ime)
real(rt),intent (IN)   :: zm     (  NZ,    ims:ime)
real(rt),intent (IN)   :: ca4k   (  NZ            )
real(rt),intent (IN)   :: ca4p   (  NZ            )
real(rt),intent (IN)   :: u      (1:NZ,    ims:ime)
real(rt),intent (IN)   :: v      (1:NZ,    ims:ime)
real(rt),intent (IN)   :: z      (0:NZ,    ims:ime)
real(rt),intent (INOUT):: sedgvar(1:NZ,npp,nvars,ims:ime)
real(rt),intent (INOUT):: bedgvar(0:NZ,    ims:ime,nvars)
real(rt),intent (IN)   :: amtx1  (1:NZ,nbf,nob+1,ims:ime)
real(rt)               :: rhsu(NZ,nob)
real(rt)               :: rhsv(NZ,nob)
real(rt)               :: Tgtu(NZ,npp)
real(rt)               :: Tgtv(NZ,npp)
INTEGER                :: ipn,k,isn,isp,n,ipp1,ipp2,ipp3,ipp4,ipp5,ippn
real(rt)               :: xyc1, xyc2  ! xyc temporaries
integer                :: startk
integer, save          :: vdmintv_handle, ipn_handle, solvei_handle

!DIR$ ASSUME_ALIGNED zm:64, ca4k:64, ca4p:64, u:64, v:64, sedgvar:64, amtx1:64
!$acc routine(solveiThLS2) vector

integer :: mythread  ! thread number when FINEGRAINED_TIMING defined
integer :: ret       ! GPTL return code when FINEGRAINED_TIMING defined

logical, save :: first = .true.

#ifdef _OPENACC
if (first) then
  first = .false.
!$acc parallel private(ret) copyout(vdmintv_handle, ipn_handle, solvei_handle)
  ret = gptlinit_handle_gpu ('vdmintv',        vdmintv_handle)
  ret = gptlinit_handle_gpu ('vdmintv_ipn',    ipn_handle)
  ret = gptlinit_handle_gpu ('vdmintv_solvei', solvei_handle)
!$acc end parallel
end if
!$acc parallel private(ret) copyin(vdmintv_handle)
  ret = gptlstart_gpu (vdmintv_handle)
!$acc end parallel
#endif
  
#ifdef FINEGRAINED_TIMING
#include <gptl.inc>
  integer :: handle = 0
  integer, external :: omp_get_thread_num
  ret = gptlstart ('vdmintv_outside')
#endif

!$acc parallel num_workers(PAR_WRK) vector_length(VEC_LEN) copyin(ipn_handle, solvei_handle)
!$acc loop gang worker private(Tgtu,Tgtv,rhsu,rhsv)
!$OMP PARALLEL DO PRIVATE(mythread,ret,n,ipp1,ipp2,ipp3,ipp4,ipp5,ippn,k, &
!$OMP                     rhsu,rhsv,isn,tgtu,tgtv,isp,xyc1,xyc2) &
!$OMP             SCHEDULE(runtime)
do ipn=ips,ipe
!$acc cache(tgtu,tgtv)
#ifdef _OPENACC
  ret = gptlstart_gpu (ipn_handle)
#endif
#ifdef FINEGRAINED_TIMING
  mythread = omp_get_thread_num ()
  if (mythread == 0) then
    ret = gptlstart_handle ('vdmintv_inside', handle)
  end if
#endif
  n    = nprox(ipn)
  ipp1 = prox(1,ipn)
  ipp2 = prox(2,ipn)
  ipp3 = prox(3,ipn)
  ipp4 = prox(4,ipn)
  ipp5 = prox(5,ipn)
  ippn = prox(n,ipn)
!$acc loop vector 
  do k=1,NZ-1
    rhsu(k,1) = cs(1,ipn)*u(k  ,ipp1)+sn(1,ipn)*v(k  ,ipp1) - u(k,ipn)
    rhsu(k,2) = cs(2,ipn)*u(k  ,ipp2)+sn(2,ipn)*v(k  ,ipp2) - u(k,ipn)
    rhsu(k,3) = cs(3,ipn)*u(k  ,ipp3)+sn(3,ipn)*v(k  ,ipp3) - u(k,ipn)
    rhsu(k,4) = cs(5,ipn)*u(k  ,ipp5)+sn(5,ipn)*v(k  ,ipp5) - u(k,ipn)
    rhsu(k,5) = cs(2,ipn)*u(k+1,ipp2)+sn(2,ipn)*v(k+1,ipp2) - u(k,ipn)
    rhsu(k,6) = cs(3,ipn)*u(k+1,ipp3)+sn(3,ipn)*v(k+1,ipp3) - u(k,ipn)
    rhsu(k,7) = cs(4,ipn)*u(k+1,ipp4)+sn(4,ipn)*v(k+1,ipp4) - u(k,ipn)
    rhsu(k,8) = cs(n,ipn)*u(k+1,ippn)+sn(n,ipn)*v(k+1,ippn) - u(k,ipn)
    rhsu(k,9) = ca4k(k)*u(k,ipn)+ca4p(k)*u(k+1,ipn) - u(k,ipn)

    rhsv(k,1) =-sn(1,ipn)*u(k  ,ipp1)+cs(1,ipn)*v(k  ,ipp1) - v(k,ipn)
    rhsv(k,2) =-sn(2,ipn)*u(k  ,ipp2)+cs(2,ipn)*v(k  ,ipp2) - v(k,ipn)
    rhsv(k,3) =-sn(3,ipn)*u(k  ,ipp3)+cs(3,ipn)*v(k  ,ipp3) - v(k,ipn)
    rhsv(k,4) =-sn(5,ipn)*u(k  ,ipp5)+cs(5,ipn)*v(k  ,ipp5) - v(k,ipn)
    rhsv(k,5) =-sn(2,ipn)*u(k+1,ipp2)+cs(2,ipn)*v(k+1,ipp2) - v(k,ipn)
    rhsv(k,6) =-sn(3,ipn)*u(k+1,ipp3)+cs(3,ipn)*v(k+1,ipp3) - v(k,ipn)
    rhsv(k,7) =-sn(4,ipn)*u(k+1,ipp4)+cs(4,ipn)*v(k+1,ipp4) - v(k,ipn)
    rhsv(k,8) =-sn(n,ipn)*u(k+1,ippn)+cs(n,ipn)*v(k+1,ippn) - v(k,ipn)
    rhsv(k,9) = ca4k(k)*v(k,ipn)+ca4p(k)*v(k+1,ipn) - v(k,ipn)
  enddo
  k=nz-1
  rhsu(k+1,1) = cs(1,ipn)*u(k  ,ipp1)+sn(1,ipn)*v(k  ,ipp1) - u(k,ipn)
  rhsu(k+1,2) = cs(2,ipn)*u(k  ,ipp2)+sn(2,ipn)*v(k  ,ipp2) - u(k,ipn)
  rhsu(k+1,3) = cs(3,ipn)*u(k  ,ipp3)+sn(3,ipn)*v(k  ,ipp3) - u(k,ipn)
  rhsu(k+1,4) = cs(5,ipn)*u(k  ,ipp5)+sn(5,ipn)*v(k  ,ipp5) - u(k,ipn)
  rhsu(k+1,5) = cs(2,ipn)*u(k+1,ipp2)+sn(2,ipn)*v(k+1,ipp2) - u(k,ipn)
  rhsu(k+1,6) = cs(3,ipn)*u(k+1,ipp3)+sn(3,ipn)*v(k+1,ipp3) - u(k,ipn)
  rhsu(k+1,7) = cs(4,ipn)*u(k+1,ipp4)+sn(4,ipn)*v(k+1,ipp4) - u(k,ipn)
  rhsu(k+1,8) = cs(n,ipn)*u(k+1,ippn)+sn(n,ipn)*v(k+1,ippn) - u(k,ipn)
  rhsu(k+1,9) = ca4k(k)*u(k,ipn)+ca4p(k)*u(k+1,ipn) - u(k,ipn)

  rhsv(k+1,1) =-sn(1,ipn)*u(k  ,ipp1)+cs(1,ipn)*v(k  ,ipp1) - v(k,ipn)
  rhsv(k+1,2) =-sn(2,ipn)*u(k  ,ipp2)+cs(2,ipn)*v(k  ,ipp2) - v(k,ipn)
  rhsv(k+1,3) =-sn(3,ipn)*u(k  ,ipp3)+cs(3,ipn)*v(k  ,ipp3) - v(k,ipn)
  rhsv(k+1,4) =-sn(5,ipn)*u(k  ,ipp5)+cs(5,ipn)*v(k  ,ipp5) - v(k,ipn)
  rhsv(k+1,5) =-sn(2,ipn)*u(k+1,ipp2)+cs(2,ipn)*v(k+1,ipp2) - v(k,ipn)
  rhsv(k+1,6) =-sn(3,ipn)*u(k+1,ipp3)+cs(3,ipn)*v(k+1,ipp3) - v(k,ipn)
  rhsv(k+1,7) =-sn(4,ipn)*u(k+1,ipp4)+cs(4,ipn)*v(k+1,ipp4) - v(k,ipn)
  rhsv(k+1,8) =-sn(n,ipn)*u(k+1,ippn)+cs(n,ipn)*v(k+1,ippn) - v(k,ipn)
  rhsv(k+1,9) = ca4k(k)*v(k,ipn)+ca4p(k)*v(k+1,ipn) - v(k,ipn)

#ifdef _OPENACC
  ret = gptlstart_gpu(solvei_handle)
#endif
  CALL solveiThLS2(nob,nbf,rhsu,rhsv,amtx1(1,1,1,ipn),ipn,ips)
#ifdef _OPENACC
  ret = gptlstop_gpu(solvei_handle)
#endif
  
!JR Defining temporary variables xyc1 and xyc2 prevents ifort from swapping the loops below
!$acc loop seq
  do isn = 1,nprox(ipn)
    xyc1 = xyc(1,isn,ipn)
    xyc2 = xyc(2,isn,ipn)
!$acc loop vector 
    do k=1,NZ-1 
      Tgtu(k,isn) = xyc1*( rhsu(k,1)*xyc1  &
                    +rhsu(k,4)*xyc2  &
                    +rhsu(k,5)*( zc(k,isn,ipn)-zm(k,ipn) )+rhsu(k,7)) &
                    +xyc2*( rhsu(k,2)*xyc2  &
                    +rhsu(k,6)*( zc(k,isn,ipn)-zm(k,ipn) )+rhsu(k,8)) &
                    +( zc(k,isn,ipn)-zm(k,ipn) ) *( rhsu(k,3)*( zc(k,isn,ipn)-zm(k,ipn) )+rhsu(k,9)) &
                    +u(k,ipn)
      Tgtv(k,isn) = xyc1*( rhsv(k,1)*xyc1  &
                    +rhsv(k,4)*xyc2  &
                    +rhsv(k,5)*( zc(k,isn,ipn)-zm(k,ipn) )+rhsv(k,7)) &
                    +xyc2*( rhsv(k,2)*xyc2  &
                    +rhsv(k,6)*( zc(k,isn,ipn)-zm(k,ipn) )+rhsv(k,8)) &
                    +( zc(k,isn,ipn)-zm(k,ipn) ) *( rhsv(k,3)*( zc(k,isn,ipn)-zm(k,ipn) )+rhsv(k,9)) &
                    +v(k,ipn)
    enddo
!TBH:  eliminate this?  Never used...  
    Tgtu(nz,isn) =  2._rt*Tgtu(nz-1,isn) - Tgtu(nz-2,isn)
    Tgtv(nz,isn) =  2._rt*Tgtv(nz-1,isn) - Tgtv(nz-2,isn)
  end do  ! isn-loop
!$acc loop seq
  do isn = 1,nprox(ipn)
    isp=mod(isn,nprox(ipn))+1
!$acc loop vector 
    do k = 2,NZ-1
      sedgvar(k,isn,1,ipn)=.25_rt*( (Tgtu(k-1,isn)+Tgtu(k-1,isp)) &
                          +       (  Tgtu(k  ,isn)+Tgtu(k  ,isp)) )
      sedgvar(k,isn,2,ipn)=.25_rt*( (Tgtv(k-1,isn)+Tgtv(k-1,isp)) &
                          +       (  Tgtv(k  ,isn)+Tgtv(k  ,isp)) )
    end do !  k -loop
    sedgvar( 1,isn,1,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,1,ipn) &
                         +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,1,ipn)
    sedgvar( 1,isn,2,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,2,ipn) &
                         +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,2,ipn)
    sedgvar(nz,isn,1,ipn)=2._rt*sedgvar(nz-1,isn,1,ipn)-sedgvar(nz-2,isn,1,ipn)
    sedgvar(nz,isn,2,ipn)=2._rt*sedgvar(nz-1,isn,2,ipn)-sedgvar(nz-2,isn,2,ipn)
  end do   ! isn-loop

!$acc loop vector 
  do k=1,NZ-1
    bedgvar(k,ipn,1)=ca4k(k)* u(k,ipn)+ca4p(k)* u(k+1,ipn)
    bedgvar(k,ipn,2)=ca4k(k)* v(k,ipn)+ca4p(k)* v(k+1,ipn)
  end do
  bedgvar(nz,ipn,1)=ca4k(NZ)* u(nz,ipn)+ca4p(NZ)* u(nz,ipn)
  bedgvar(nz,ipn,2)=ca4k(NZ)* v(nz,ipn)+ca4p(NZ)* v(nz,ipn)
  bedgvar( 0, ipn,1)= (z(0,ipn)-z(2,ipn))/(z(1,ipn)-z(2,ipn))*bedgvar(1,ipn,1) &
                    + (z(0,ipn)-z(1,ipn))/(z(2,ipn)-z(1,ipn))*bedgvar(2,ipn,1)
  bedgvar( 0, ipn,2)= (z(0,ipn)-z(2,ipn))/(z(1,ipn)-z(2,ipn))*bedgvar(1,ipn,2) &
                    + (z(0,ipn)-z(1,ipn))/(z(2,ipn)-z(1,ipn))*bedgvar(2,ipn,2)
#ifdef FINEGRAINED_TIMING
  if (mythread == 0) then
    ret = gptlstop_handle ('vdmintv_inside', handle)
  end if
#endif

#ifdef _OPENACC
  ret = gptlstop_gpu (ipn_handle)
#endif
enddo !ipn-loop
!$acc end parallel
!$OMP END PARALLEL DO

#ifdef _OPENACC
!$acc parallel private(ret) copyin(vdmintv_handle)
ret = gptlstop_gpu (vdmintv_handle)
!$acc end parallel
#endif

#ifdef FINEGRAINED_TIMING
ret = gptlstop ('vdmintv_outside')
#endif

return
end subroutine vdmintv

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

# include "solveiThLS2.F90"
