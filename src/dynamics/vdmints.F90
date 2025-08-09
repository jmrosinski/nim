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

subroutine vdmints(trc1,sedgvar,bedgvar,nvars,amtx1,nbf,npp,ims,ime,nob,ntr1,&
                   nprox,prox,xyc,zc,zm,z,ca4k,ca4p,ips,ipe)
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
implicit none
integer ,intent (IN)   :: npp,ims,ime,nvars,nob,ntr1,nbf,ips,ipe
integer ,intent (IN)   :: nprox(          ims:ime) ! holds number of proximity points
integer ,intent (IN)   :: prox (     npp, ims:ime) ! holds index of proximity points
real(rt),intent (IN)   :: xyc  (   2,npp, ims:ime)
real(rt),intent (IN)   :: zc   (0:NZ,npp, ims:ime)
real(rt),intent (IN)   :: zm   (  NZ,     ims:ime)
real(rt),intent (IN)   :: z    (0:NZ,     ims:ime)
real(rt),intent (IN)   :: ca4k (  nz             )
real(rt),intent (IN)   :: ca4p (  nz             )
real(rt),intent (IN)   :: trc1 (  NZ,ntr1,ims:ime)
real(rt),intent (INOUT):: sedgvar (  1:NZ,npp,nvars,ims:ime)
real(rt),intent (INOUT):: bedgvar (  0:NZ,    ims:ime,nvars)
real(rt),intent (IN)   :: amtx1   (  1:NZ,(nob+1)*nbf,ims:ime)
real(rt)               :: tgtc
real(rt)               :: rhs  (NZ,9)
real(rt)               :: Tgt  (NZ,npp)
INTEGER                :: ipn,k,isn,isp,ism,nvi
real(rt)               :: xyc1, xyc2  ! xyc temporaries

!$acc routine(solveiThLS) vector
!$acc parallel num_workers(PAR_WRK) vector_length(VEC_LEN)

!$OMP PARALLEL DO PRIVATE(k,rhs,isn,tgtc,tgt,isp,ism,xyc1,xyc2) SCHEDULE(runtime)
!$acc loop gang private(rhs,tgt)
do ipn=ips,ipe
!$acc cache(Tgt,rhs)
!$acc loop seq
do nvi=1,ntr1
!$acc loop vector
  do k=1,NZ-1
    rhs(k,1) = trc1(k  ,nvi,prox(1         ,ipn)) - trc1(k,nvi,ipn)
    rhs(k,2) = trc1(k  ,nvi,prox(2         ,ipn)) - trc1(k,nvi,ipn)
    rhs(k,3) = trc1(k  ,nvi,prox(3         ,ipn)) - trc1(k,nvi,ipn)
    rhs(k,4) = trc1(k  ,nvi,prox(5         ,ipn)) - trc1(k,nvi,ipn)
    rhs(k,5) = trc1(k+1,nvi,prox(2         ,ipn)) - trc1(k,nvi,ipn)
    rhs(k,6) = trc1(k+1,nvi,prox(3         ,ipn)) - trc1(k,nvi,ipn)
    rhs(k,7) = trc1(k+1,nvi,prox(4         ,ipn)) - trc1(k,nvi,ipn)
    rhs(k,8) = trc1(k+1,nvi,prox(nprox(ipn),ipn)) - trc1(k,nvi,ipn)
    rhs(k,9) = ca4k(k)*trc1(k,nvi,ipn)+ca4p(k)*trc1(k+1,nvi,ipn)-trc1(k,nvi,ipn)
  enddo !k-loop
  k = nz-1
  rhs(k+1,1) = trc1(k  ,nvi,prox(1         ,ipn)) - trc1(k,nvi,ipn)
  rhs(k+1,2) = trc1(k  ,nvi,prox(2         ,ipn)) - trc1(k,nvi,ipn)
  rhs(k+1,3) = trc1(k  ,nvi,prox(3         ,ipn)) - trc1(k,nvi,ipn)
  rhs(k+1,4) = trc1(k  ,nvi,prox(5         ,ipn)) - trc1(k,nvi,ipn)
  rhs(k+1,5) = trc1(k+1,nvi,prox(2         ,ipn)) - trc1(k,nvi,ipn)
  rhs(k+1,6) = trc1(k+1,nvi,prox(3         ,ipn)) - trc1(k,nvi,ipn)
  rhs(k+1,7) = trc1(k+1,nvi,prox(4         ,ipn)) - trc1(k,nvi,ipn)
  rhs(k+1,8) = trc1(k+1,nvi,prox(nprox(ipn),ipn)) - trc1(k,nvi,ipn)
  rhs(k+1,9) = ca4k(k)*trc1(k,nvi,ipn)+ca4p(k)*trc1(k+1,nvi,ipn)-trc1(k,nvi,ipn)
  CALL solveiThLS(nob,nbf,rhs,amtx1(1,1,ipn),ipn,ips)
!JR Defining temporary variables xyc1 and xyc2 prevents ifort from swapping the loops below
  do isn = 1,nprox(ipn)
    xyc1 = xyc(1,isn,ipn)
    xyc2 = xyc(2,isn,ipn)
!$acc loop vector
    do k=1,NZ-1
      tgtc = ( zc(k,isn,ipn)-zm(k,ipn) )
      Tgt(k,isn) = xyc1*( rhs(k,1)*xyc1 &
                         +rhs(k,4)*xyc2 &
                         +rhs(k,5)*tgtc+rhs(k,7)) &
                 + xyc2*( rhs(k,2)*xyc2 &
                         +rhs(k,6)*tgtc+rhs(k,8)) &
                 + tgtc*( rhs(k,3)*tgtc+rhs(k,9)) &
                 + trc1(k,nvi,ipn)
    enddo !k-loop
  end do  ! isn-loop
!$acc loop seq
  do isn = 1,nprox(ipn)
    isp=mod(isn,nprox(ipn))+1
    ism=isn-1
    if(ism.eq.0) ism=nprox(ipn)
!$acc loop vector
    do k = 2,NZ-1
      sedgvar(k,isn,nvi,ipn) = .5_rt*( .5_rt*(Tgt(k-1 ,isn)+Tgt(k-1 ,isp)) &
                            +  .5_rt*(     Tgt( k  ,isn)+Tgt( k  ,isp))  )
    end do !  k -loop
    sedgvar( 1,isn,nvi,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,nvi,ipn) &
                           +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,nvi,ipn) 
    sedgvar(NZ,isn,nvi,ipn)=2._rt*sedgvar(nz-1,isn,nvi,ipn)-sedgvar(nz-2,isn,nvi,ipn)
  end do  ! isn-loop

!$acc loop vector
  do k=1,NZ-1
    bedgvar(k,ipn,nvi)=ca4k(k)* trc1(k,nvi,ipn)+ca4p(k)* trc1(k+1,nvi,ipn)
  end do
  bedgvar( 0, ipn,nvi) = (zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*bedgvar(1,ipn,nvi) &
                      +  (zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*bedgvar(2,ipn,nvi) 
  bedgvar(NZ,ipn,nvi)=ca4k(nz)* trc1(NZ,nvi,ipn)+ca4p(nz)* trc1(NZ,nvi,ipn)
enddo ! nvi loop
enddo !ipn-loop
!$OMP END PARALLEL DO
!$acc end parallel

return
end subroutine vdmints

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
#include "solveiThLS.F90"
