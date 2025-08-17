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

subroutine vdmints3(rp,tp,trp,teb,sedgvar,bedgvar,nvars,nvi,amtx1,nbf,npp,ims,ime,nob,&
                    nprox,prox,xyc,zc,zm,ca4k,ca4p,ips,ipe)
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
use gptl_acc
use kinds, only: rt
use ReadNamelist, only: nz
use openacc
implicit none
integer ,intent (IN)   :: nvi,npp,ims,ime,nvars,nob,nbf,ips,ipe
integer ,intent (IN)   :: nprox(           ims:ime) ! holds number of proximity points
integer ,intent (IN)   :: prox (       npp,ims:ime) ! holds index of proximity points
real(rt),intent (IN)   :: xyc  (     2,npp,ims:ime)
real(rt),intent (IN)   :: zc   (  0:NZ,npp,ims:ime)
real(rt),intent (IN)   :: zm   (    NZ,    ims:ime)
real(rt),intent (IN)   :: ca4k (    NZ            )
real(rt),intent (IN)   :: ca4p (    NZ            )
real(rt),intent (IN)   :: rp   (    NZ,    ims:ime) ! nvi
real(rt),intent (IN)   :: tp   (    NZ,    ims:ime) ! nvi+1
real(rt),intent (IN)   :: trp  (    NZ,    ims:ime) ! nvi+2
real(rt),intent (IN)   :: teb  (    NZ,npp,ims:ime)
real(rt),intent (INOUT):: sedgvar (  1:NZ,npp,nvars,ims:ime)
real(rt),intent (INOUT):: bedgvar (  0:NZ,    ims:ime,nvars)
real(rt),intent (IN)   :: amtx1   (  1:NZ,(nob+1)*nbf,ims:ime)
real(rt)               :: tgtc
real(rt)               :: rhs1 (NZ,9)
real(rt)               :: rhs2 (NZ,9)
real(rt)               :: rhs3 (NZ,9)
real(rt)               :: Tgt1 (NZ,npp)
real(rt)               :: Tgt2 (NZ,npp)
real(rt)               :: Tgt3 (NZ,npp)
INTEGER                :: ipn,k,isn,isp,nv2,nv3
real(rt)               :: xyc1, xyc2  ! xyc temporaries
integer                :: startk

integer :: mythread  ! thread number when FINEGRAINED_TIMING defined
integer :: ret
integer, save :: vdmints3_handle, ipn_handle, kloop1_handle, kloop2_handle, k3_handle, k4_handle
integer, save :: isn1_handle, scalar_handle, solvei_handle, isn2_handle
logical, save :: first = .true.

!$acc routine(solveiThLS3) vector

#define INNER_TIMERS

if (first) then
  first = .false.
  ret = gptlstart('vdmints3_sethandles'//char(0))
!$acc parallel private(ret) copyout(vdmints3_handle, ipn_handle, kloop1_handle, kloop2_handle, &
!$acc&                              k3_handle, k4_handle, isn1_handle, isn2_handle, &  
!$acc&                              scalar_handle, solvei_handle)
  ret = gptlinit_handle_gpu ('vdmints3',        vdmints3_handle)
  ret = gptlinit_handle_gpu ('vdmints3_ipn',    ipn_handle)
  ret = gptlinit_handle_gpu ('vdmints3_kloop1', kloop1_handle)
  ret = gptlinit_handle_gpu ('vdmints3_kloop2', kloop2_handle)
  ret = gptlinit_handle_gpu ('vd3_k3',          k3_handle)
  ret = gptlinit_handle_gpu ('vd3_k4',          k4_handle)
  ret = gptlinit_handle_gpu ('vdmints3_isn1',   isn1_handle)
  ret = gptlinit_handle_gpu ('vdmints3_isn2',   isn2_handle)
  ret = gptlinit_handle_gpu ('vdmints3_scalar', scalar_handle)
  ret = gptlinit_handle_gpu ('vdmints3_solvei', solvei_handle)
  ret = gptlstart_gpu (vdmints3_handle)
!$acc end parallel
  ret = gptlstop('vdmints3_sethandles'//char(0))
else
!$acc parallel private(ret) copyin(vdmints3_handle)
  ret = gptlstart_gpu (vdmints3_handle)
!$acc end parallel
end if

nv2=nvi+1
nv3=nvi+2
#ifdef INNER_TIMERS
!$acc parallel private(ret) num_gangs(10242) num_workers(3) vector_length(32), &
!$acc&         copyin(ipn_handle, kloop1_handle, kloop2_handle, k3_handle, isn1_handle, &
!$acc&                isn2_handle, scalar_handle, solvei_handle)
#else
!$acc parallel private(ret) num_workers(PAR_WRK) vector_length(VEC_LEN)
#endif

!$acc loop gang private(rhs1,rhs2,rhs3,Tgt1,Tgt2,Tgt3)

!$OMP PARALLEL DO PRIVATE(mythread,ret,k,rhs1,rhs2,rhs3,isn,tgtc,tgt1,tgt2, &
!$OMP                     tgt3,isp,xyc1,xyc2) SCHEDULE(runtime)
do ipn=ips,ipe
#ifdef INNER_TIMERS
  ret = gptlstart_gpu (ipn_handle)
  ret = gptlstart_gpu (kloop1_handle)
#endif
!!$acc cache(rhs1,rhs2,rhs3)
!$acc cache(rhs1,rhs2)
!$acc loop worker
  do startk=1,NZ-1,32
!$acc loop vector
  do k=startk,min(NZ-1,startk+32-1)

    rhs1(k,1) = rp(k  ,prox(1         ,ipn)) - rp(k,ipn)
    rhs1(k,2) = rp(k  ,prox(2         ,ipn)) - rp(k,ipn)
    rhs1(k,3) = rp(k  ,prox(3         ,ipn)) - rp(k,ipn)
    rhs1(k,4) = rp(k  ,prox(5         ,ipn)) - rp(k,ipn)
    rhs1(k,5) = rp(k+1,prox(2         ,ipn)) - rp(k,ipn)
    rhs1(k,6) = rp(k+1,prox(3         ,ipn)) - rp(k,ipn)
    rhs1(k,7) = rp(k+1,prox(4         ,ipn)) - rp(k,ipn)
    rhs1(k,8) = rp(k+1,prox(nprox(ipn),ipn)) - rp(k,ipn)
    rhs1(k,9) = ca4k(k)*rp(k,ipn)+ca4p(k)*rp(k+1,ipn) - rp(k,ipn)
    rhs2(k,1) = tp(k  ,prox(1         ,ipn)) - tp(k,ipn)
    rhs2(k,2) = tp(k  ,prox(2         ,ipn)) - tp(k,ipn)
    rhs2(k,3) = tp(k  ,prox(3         ,ipn)) - tp(k,ipn)
    rhs2(k,4) = tp(k  ,prox(5         ,ipn)) - tp(k,ipn)
    rhs2(k,5) = tp(k+1,prox(2         ,ipn)) - tp(k,ipn)
    rhs2(k,6) = tp(k+1,prox(3         ,ipn)) - tp(k,ipn)
    rhs2(k,7) = tp(k+1,prox(4         ,ipn)) - tp(k,ipn)
    rhs2(k,8) = tp(k+1,prox(nprox(ipn),ipn)) - tp(k,ipn)
    rhs2(k,9) = ca4k(k)*tp(k,ipn)+ca4p(k)*tp(k+1,ipn) - tp(k,ipn)
    rhs3(k,1) = trp(k  ,prox(1         ,ipn)) - trp(k,ipn)
    rhs3(k,2) = trp(k  ,prox(2         ,ipn)) - trp(k,ipn)
    rhs3(k,3) = trp(k  ,prox(3         ,ipn)) - trp(k,ipn)
    rhs3(k,4) = trp(k  ,prox(5         ,ipn)) - trp(k,ipn)
    rhs3(k,5) = trp(k+1,prox(2         ,ipn)) - trp(k,ipn)
    rhs3(k,6) = trp(k+1,prox(3         ,ipn)) - trp(k,ipn)
    rhs3(k,7) = trp(k+1,prox(4         ,ipn)) - trp(k,ipn)
    rhs3(k,8) = trp(k+1,prox(nprox(ipn),ipn)) - trp(k,ipn)
    rhs3(k,9) = ca4k(k)*trp(k,ipn)+ca4p(k)*trp(k+1,ipn) - trp(k,ipn)
  enddo !k-loop
  enddo !k-loop
#ifdef INNER_TIMERS
  ret = gptlstop_gpu (kloop1_handle)
  ret = gptlstart_gpu(scalar_handle)
#endif
  k=NZ-1
  rhs1(k+1,1) = rp(k  ,prox(1         ,ipn)) - rp(k,ipn)
  rhs1(k+1,2) = rp(k  ,prox(2         ,ipn)) - rp(k,ipn)
  rhs1(k+1,3) = rp(k  ,prox(3         ,ipn)) - rp(k,ipn)
  rhs1(k+1,4) = rp(k  ,prox(5         ,ipn)) - rp(k,ipn)
  rhs1(k+1,5) = rp(k+1,prox(2         ,ipn)) - rp(k,ipn)
  rhs1(k+1,6) = rp(k+1,prox(3         ,ipn)) - rp(k,ipn)
  rhs1(k+1,7) = rp(k+1,prox(4         ,ipn)) - rp(k,ipn)
  rhs1(k+1,8) = rp(k+1,prox(nprox(ipn),ipn)) - rp(k,ipn)
  rhs1(k+1,9) = ca4k(k)*rp(k,ipn)+ca4p(k)*rp(k+1,ipn) - rp(k,ipn)
  rhs2(k+1,1) = tp(k  ,prox(1         ,ipn)) - tp(k,ipn)
  rhs2(k+1,2) = tp(k  ,prox(2         ,ipn)) - tp(k,ipn)
  rhs2(k+1,3) = tp(k  ,prox(3         ,ipn)) - tp(k,ipn)
  rhs2(k+1,4) = tp(k  ,prox(5         ,ipn)) - tp(k,ipn)
  rhs2(k+1,5) = tp(k+1,prox(2         ,ipn)) - tp(k,ipn)
  rhs2(k+1,6) = tp(k+1,prox(3         ,ipn)) - tp(k,ipn)
  rhs2(k+1,7) = tp(k+1,prox(4         ,ipn)) - tp(k,ipn)
  rhs2(k+1,8) = tp(k+1,prox(nprox(ipn),ipn)) - tp(k,ipn)
  rhs2(k+1,9) = ca4k(k)*tp(k,ipn)+ca4p(k)*tp(k+1,ipn) - tp(k,ipn)
  rhs3(k+1,1) = trp(k  ,prox(1         ,ipn)) - trp(k,ipn)
  rhs3(k+1,2) = trp(k  ,prox(2         ,ipn)) - trp(k,ipn)
  rhs3(k+1,3) = trp(k  ,prox(3         ,ipn)) - trp(k,ipn)
  rhs3(k+1,4) = trp(k  ,prox(5         ,ipn)) - trp(k,ipn)
  rhs3(k+1,5) = trp(k+1,prox(2         ,ipn)) - trp(k,ipn)
  rhs3(k+1,6) = trp(k+1,prox(3         ,ipn)) - trp(k,ipn)
  rhs3(k+1,7) = trp(k+1,prox(4         ,ipn)) - trp(k,ipn)
  rhs3(k+1,8) = trp(k+1,prox(nprox(ipn),ipn)) - trp(k,ipn)
  rhs3(k+1,9) = ca4k(k)*trp(k,ipn)+ca4p(k)*trp(k+1,ipn) - trp(k,ipn)    
#ifdef INNER_TIMERS
  ret = gptlstop_gpu (scalar_handle)
  ret = gptlstart_gpu(solvei_handle)
#endif
!$acc loop worker
  do startk=1,NZ,32
    CALL solveiThLS3(nob,nbf,rhs1,rhs2,rhs3,amtx1(1,1,ipn),startk)
  end do
#ifdef INNER_TIMERS
  ret = gptlstop_gpu(solvei_handle)
  ret = gptlstart_gpu(isn1_handle)
#endif

!JR Defining temporary variables xyc1 and xyc2 prevents ifort from swapping the loops below
!$acc loop worker private(xyc1,xyc2)
  do isn = 1,nprox(ipn)
    xyc1 = xyc(1,isn,ipn)
    xyc2 = xyc(2,isn,ipn)
!$acc loop vector private(tgtc)
    do k=1,NZ-1
!      ret = gptlstart_gpu(kloop2_handle)
      tgtc = ( zc(k,isn,ipn)-zm(k,ipn) )
      Tgt1(k,isn) =xyc1*( rhs1(k,1)*xyc1  &
                         +rhs1(k,4)*xyc2  &
                         +rhs1(k,5)*tgtc+rhs1(k,7)) &
                 + xyc2*( rhs1(k,2)*xyc2  &
                         +rhs1(k,6)*tgtc+rhs1(k,8)) &
                 + tgtc*( rhs1(k,3)*tgtc+rhs1(k,9)) &
                 + rp(k,ipn)
      Tgt2(k,isn) =xyc1*( rhs2(k,1)*xyc1  &
                         +rhs2(k,4)*xyc2  &
                         +rhs2(k,5)*tgtc+rhs2(k,7)) &
                 + xyc2*( rhs2(k,2)*xyc2  &
                         +rhs2(k,6)*tgtc+rhs2(k,8)) &
                 + tgtc*( rhs2(k,3)*tgtc+rhs2(k,9)) &
                 + tp(k,ipn)
      Tgt3(k,isn) =xyc1*( rhs3(k,1)*xyc1  &
                         +rhs3(k,4)*xyc2  &
                         +rhs3(k,5)*tgtc+rhs3(k,7)) &
                 + xyc2*( rhs3(k,2)*xyc2  &
                         +rhs3(k,6)*tgtc+rhs3(k,8)) &
                 + tgtc*( rhs3(k,3)*tgtc+rhs3(k,9)) &
                 + trp(k,ipn)
!      ret = gptlstop_gpu(kloop2_handle)
    enddo !k-loop
  end do  ! isn-loop

#ifdef INNER_TIMERS
  ret = gptlstop_gpu(isn1_handle)
  ret = gptlstart_gpu(isn2_handle)
#endif

!$acc loop worker private(isp)
  do isn = 1,nprox(ipn)
    isp=mod(isn,nprox(ipn))+1
!$acc loop vector
    do k = 2,NZ-1
!      ret = gptlstart_gpu(k3_handle)
      sedgvar(k,isn,nvi,ipn)=.5_rt*(.5_rt*(Tgt1(k-1,isn)+Tgt1(k-1 ,isp)) &
                         +          .5_rt*(Tgt1(k  ,isn)+Tgt1(k   ,isp)))
      sedgvar(k,isn,nv2,ipn)=.5_rt*(.5_rt*(Tgt2(k-1,isn)+Tgt2(k-1 ,isp)) &
                         +          .5_rt*(Tgt2(k  ,isn)+Tgt2(k   ,isp)))
      sedgvar(k,isn,nv3,ipn)=.5_rt*(.5_rt*(Tgt3(k-1,isn)+Tgt3(k-1 ,isp)) &
                         +          .5_rt*(Tgt3(k  ,isn)+Tgt3(k   ,isp)))
!      ret = gptlstop_gpu(k3_handle)
    end do

#ifdef INNER_TIMERS
    ret = gptlstart_gpu (scalar_handle)
#endif
    sedgvar( 1,isn,nvi,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,nvi,ipn) &
                           +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,nvi,ipn)
    sedgvar( 1,isn,nv2,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,nv2,ipn) &
                           +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,nv2,ipn)
    sedgvar( 1,isn,nv3,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,nv3,ipn) &
                           +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,nv3,ipn)
    sedgvar(NZ,isn,nvi,ipn)=2._rt*sedgvar(NZ-1,isn,nvi,ipn)-sedgvar(NZ-2,isn,nvi,ipn)
    sedgvar(NZ,isn,nv2,ipn)=2._rt*sedgvar(NZ-1,isn,nv2,ipn)-sedgvar(NZ-2,isn,nv2,ipn)
    sedgvar(NZ,isn,nv3,ipn)=2._rt*sedgvar(NZ-1,isn,nv3,ipn)-sedgvar(NZ-2,isn,nv3,ipn)
#ifdef INNER_TIMERS
    ret = gptlstop_gpu (scalar_handle)
#endif
!$acc loop vector
    do k = 1,NZ
      sedgvar(k,isn,5,ipn) = teb(k,isn,ipn)+sedgvar(k,isn,5,ipn)
    enddo
  end do  ! isn-loop

#ifdef INNER_TIMERS
  ret = gptlstop_gpu(isn2_handle)
  ret = gptlstart_gpu(k4_handle)
#endif

!$acc loop vector
  do k=1,NZ-1
    bedgvar(k,ipn,nvi)=ca4k(k)* rp(k,ipn)+ca4p(k)* rp(k+1,ipn)
    bedgvar(k,ipn,nv2)=ca4k(k)* tp(k,ipn)+ca4p(k)* tp(k+1,ipn)
    bedgvar(k,ipn,nv3)=ca4k(k)*trp(k,ipn)+ca4p(k)*trp(k+1,ipn)
  end do
#ifdef INNER_TIMERS
  ret = gptlstop_gpu(k4_handle)
  ret = gptlstart_gpu(scalar_handle)
#endif
  bedgvar( 0, ipn,nvi)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*bedgvar(1,ipn,nvi) &
                      +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*bedgvar(2,ipn,nvi) 
  bedgvar( 0, ipn,nv2)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*bedgvar(1,ipn,nv2) &
                      +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*bedgvar(2,ipn,nv2) 
  bedgvar( 0, ipn,nv3)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*bedgvar(1,ipn,nv3) &
                      +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*bedgvar(2,ipn,nv3) 
  bedgvar(NZ,ipn,nvi)=ca4k(NZ)* rp(NZ,ipn)+ca4p(NZ)* rp(NZ,ipn)
  bedgvar(NZ,ipn,nv2)=ca4k(NZ)* tp(NZ,ipn)+ca4p(NZ)* tp(NZ,ipn)
  bedgvar(NZ,ipn,nv3)=ca4k(NZ)*trp(NZ,ipn)+ca4p(NZ)*trp(NZ,ipn)

#ifdef INNER_TIMERS
  ret = gptlstop_gpu(scalar_handle)
  ret = gptlstop_gpu(ipn_handle)
#endif
enddo !ipn-loop
!$acc end parallel
!$OMP END PARALLEL DO

!$acc parallel private(ret) copyin(vdmints3_handle)
ret = gptlstop_gpu (vdmints3_handle)
!$acc end parallel

return
end subroutine vdmints3

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
#include "solveiThLS3.F90"
