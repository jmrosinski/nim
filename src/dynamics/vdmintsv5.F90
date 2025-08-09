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

subroutine vdmintsv5(u,v,rp,tp,trp,teb,sedgvar,bedgvar,nvars,nvi,amtx1,  &
                     nbf,npp,ims,ime,nob,nprox,prox,cs,sn,xyc,zc,zm,z,&
                     ca4k,ca4p,ips,ipe)
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

use ReadNamelist, only: nz
implicit none
integer,intent (IN)    :: nvi,npp,ims,ime,nvars,nob,nbf,ips,ipe
integer,intent (IN)    :: nprox(           ims:ime) ! holds number of proximity points
integer,intent (IN)    :: prox (       npp,ims:ime) ! holds index of proximity points
real   ,intent (IN)    :: cs     (     npp,ims:ime) ! cosine transform function
real   ,intent (IN)    :: sn     (     npp,ims:ime) ! sine transform function
real   ,intent (IN)    :: xyc  (     2,npp,ims:ime)
real   ,intent (IN)    :: zc   (  0:nz,npp,ims:ime)
real   ,intent (IN)    :: zm   (    nz,    ims:ime)
real   ,intent (IN)    :: z    (  0:nz,    ims:ime)
real   ,intent (IN)    :: ca4k (    nz            )
real   ,intent (IN)    :: ca4p (    nz            )
real   ,intent (IN)    :: u      (1:nz,    ims:ime) ! 1
real   ,intent (IN)    :: v      (1:nz,    ims:ime) ! 2
real   ,intent (IN)    :: rp   (    nz,    ims:ime) ! nvi
real   ,intent (IN)    :: tp   (    nz,    ims:ime) ! nvi+1
real   ,intent (IN)    :: trp  (    nz,    ims:ime) ! nvi+2
real   ,intent (IN)    :: teb  (    nz,npp,ims:ime)
real   ,intent (INOUT) :: sedgvar (  1:nz,npp,nvars,ims:ime)
real   ,intent (INOUT) :: bedgvar (  0:nz,    ims:ime,nvars)
real   ,intent (IN)    :: amtx1   (  1:nz,(nob+1)*nbf,ims:ime)
REAL                   :: tgtc
REAL                   :: rhs1 (nz,9)
REAL                   :: rhs2 (nz,9)
REAL                   :: rhs3 (nz,9)
REAL                   :: Tgt1 (nz,npp)
REAL                   :: Tgt2 (nz,npp)
REAL                   :: Tgt3 (nz,npp)
REAL                   :: rhsu(nz,nob)
REAL                   :: rhsv(nz,nob)
REAL                   :: Tgtu(nz,npp)
REAL                   :: Tgtv(nz,npp)
INTEGER                :: ipn,k,isn,isp,ism,nv2,nv3,n,ipp1,ipp2,ipp3,ipp4,ipp5,ippn
real                   :: xyc1, xyc2  ! xyc temporaries

nv2=nvi+1
nv3=nvi+2

!$OMP PARALLEL DO PRIVATE(k,rhsu,rhsv,rhs1,rhs2,rhs3,isn,tgtu,tgtv,tgtc, &
!$OMP                     tgt1,tgt2,tgt3,isp,ism,xyc1,xyc2,n,ipp1,ipp2,ipp3, &
!$OMP                     ipp4,ipp5,ippn) SCHEDULE(runtime)
do ipn=ips,ipe
    n    = nprox(ipn)
    ipp1 = prox(1,ipn)
    ipp2 = prox(2,ipn)
    ipp3 = prox(3,ipn)
    ipp4 = prox(4,ipn)
    ipp5 = prox(5,ipn)
    ippn = prox(n,ipn)
  do k=1,nz-1

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

    rhs1(k,1) = rp(k  ,ipp1) - rp(k,ipn)
    rhs1(k,2) = rp(k  ,ipp2) - rp(k,ipn)
    rhs1(k,3) = rp(k  ,ipp3) - rp(k,ipn)
    rhs1(k,4) = rp(k  ,ipp5) - rp(k,ipn)
    rhs1(k,5) = rp(k+1,ipp2) - rp(k,ipn)
    rhs1(k,6) = rp(k+1,ipp3) - rp(k,ipn)
    rhs1(k,7) = rp(k+1,ipp4) - rp(k,ipn)
    rhs1(k,8) = rp(k+1,ippn) - rp(k,ipn)
    rhs1(k,9) = ca4k(k)*rp(k,ipn)+ca4p(k)*rp(k+1,ipn) - rp(k,ipn)
    rhs2(k,1) = tp(k  ,ipp1) - tp(k,ipn)
    rhs2(k,2) = tp(k  ,ipp2) - tp(k,ipn)
    rhs2(k,3) = tp(k  ,ipp3) - tp(k,ipn)
    rhs2(k,4) = tp(k  ,ipp5) - tp(k,ipn)
    rhs2(k,5) = tp(k+1,ipp2) - tp(k,ipn)
    rhs2(k,6) = tp(k+1,ipp3) - tp(k,ipn)
    rhs2(k,7) = tp(k+1,ipp4) - tp(k,ipn)
    rhs2(k,8) = tp(k+1,ippn) - tp(k,ipn)
    rhs2(k,9) = ca4k(k)*tp(k,ipn)+ca4p(k)*tp(k+1,ipn) - tp(k,ipn)
    rhs3(k,1) = trp(k  ,ipp1) - trp(k,ipn)
    rhs3(k,2) = trp(k  ,ipp2) - trp(k,ipn)
    rhs3(k,3) = trp(k  ,ipp3) - trp(k,ipn)
    rhs3(k,4) = trp(k  ,ipp5) - trp(k,ipn)
    rhs3(k,5) = trp(k+1,ipp2) - trp(k,ipn)
    rhs3(k,6) = trp(k+1,ipp3) - trp(k,ipn)
    rhs3(k,7) = trp(k+1,ipp4) - trp(k,ipn)
    rhs3(k,8) = trp(k+1,ippn) - trp(k,ipn)
    rhs3(k,9) = ca4k(k)*trp(k,ipn)+ca4p(k)*trp(k+1,ipn) - trp(k,ipn)
  enddo !k-loop
!TBH:  eliminate this?  Never used...  
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

    rhs1(k+1,1) = rp(k  ,ipp1) - rp(k,ipn)
    rhs1(k+1,2) = rp(k  ,ipp2) - rp(k,ipn)
    rhs1(k+1,3) = rp(k  ,ipp3) - rp(k,ipn)
    rhs1(k+1,4) = rp(k  ,ipp5) - rp(k,ipn)
    rhs1(k+1,5) = rp(k+1,ipp2) - rp(k,ipn)
    rhs1(k+1,6) = rp(k+1,ipp3) - rp(k,ipn)
    rhs1(k+1,7) = rp(k+1,ipp4) - rp(k,ipn)
    rhs1(k+1,8) = rp(k+1,ippn) - rp(k,ipn)
  rhs1(k+1,9) = ca4k(k)*rp(k,ipn)+ca4p(k)*rp(k+1,ipn) - rp(k,ipn)
    rhs2(k+1,1) = tp(k  ,ipp1) - tp(k,ipn)
    rhs2(k+1,2) = tp(k  ,ipp2) - tp(k,ipn)
    rhs2(k+1,3) = tp(k  ,ipp3) - tp(k,ipn)
    rhs2(k+1,4) = tp(k  ,ipp5) - tp(k,ipn)
    rhs2(k+1,5) = tp(k+1,ipp2) - tp(k,ipn)
    rhs2(k+1,6) = tp(k+1,ipp3) - tp(k,ipn)
    rhs2(k+1,7) = tp(k+1,ipp4) - tp(k,ipn)
    rhs2(k+1,8) = tp(k+1,ippn) - tp(k,ipn)
  rhs2(k+1,9) = ca4k(k)*tp(k,ipn)+ca4p(k)*tp(k+1,ipn) - tp(k,ipn)
    rhs3(k+1,1) = trp(k  ,ipp1) - trp(k,ipn)
    rhs3(k+1,2) = trp(k  ,ipp2) - trp(k,ipn)
    rhs3(k+1,3) = trp(k  ,ipp3) - trp(k,ipn)
    rhs3(k+1,4) = trp(k  ,ipp5) - trp(k,ipn)
    rhs3(k+1,5) = trp(k+1,ipp2) - trp(k,ipn)
    rhs3(k+1,6) = trp(k+1,ipp3) - trp(k,ipn)
    rhs3(k+1,7) = trp(k+1,ipp4) - trp(k,ipn)
    rhs3(k+1,8) = trp(k+1,ippn) - trp(k,ipn)
  rhs3(k+1,9) = ca4k(k)*trp(k,ipn)+ca4p(k)*trp(k+1,ipn) - trp(k,ipn)

    CALL solveiThLS5(nob,nbf,rhsu,rhsv,rhs1,rhs2,rhs3,amtx1(1,1,ipn))

!JR Defining temporary variables xyc1 and xyc2 prevents ifort from swapping the loops below
  do isn = 1,nprox(ipn)
    xyc1 = xyc(1,isn,ipn)
    xyc2 = xyc(2,isn,ipn)
    do k=1,nz-1
      tgtc = ( zc(k,isn,ipn)-zm(k,ipn) )
      Tgtu(k,isn) =xyc1*( rhsu(k,1)*xyc1  &
                         +rhsu(k,4)*xyc2  &
                         +rhsu(k,5)*tgtc+rhsu(k,7)) &
                 + xyc2*( rhsu(k,2)*xyc2  &
                         +rhsu(k,6)*tgtc+rhsu(k,8)) &
                 + tgtc*( rhsu(k,3)*tgtc+rhsu(k,9)) &
                 + u(k,ipn)
      Tgtv(k,isn) =xyc1*( rhsv(k,1)*xyc1  &
                         +rhsv(k,4)*xyc2  &
                         +rhsv(k,5)*tgtc+rhsv(k,7)) &
                 + xyc2*( rhsv(k,2)*xyc2  &
                         +rhsv(k,6)*tgtc+rhsv(k,8)) &
                 + tgtc*( rhsv(k,3)*tgtc+rhsv(k,9)) &
                 + v(k,ipn)
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
    enddo !k-loop
!TBH:  eliminate this?  Never used...  
    Tgtu(nz,isn) =  2.*Tgtu(nz-1,isn) - Tgtu(nz-2,isn)
    Tgtv(nz,isn) =  2.*Tgtv(nz-1,isn) - Tgtv(nz-2,isn)
  end do  ! isn-loop

  do isn = 1,nprox(ipn)
    isp=mod(isn,nprox(ipn))+1
    ism=isn-1
    if(ism.eq.0) ism=nprox(ipn)
    do k = 2,nz-1
      sedgvar(k,isn,  1,ipn)=.25*(  (Tgtu(k-1,isn)+Tgtu(k-1,isp)) &
                          +          (Tgtu(k  ,isn)+Tgtu(k  ,isp)) )
      sedgvar(k,isn,  2,ipn)=.25*(  (Tgtv(k-1,isn)+Tgtv(k-1,isp)) &
                          +          (Tgtv(k  ,isn)+Tgtv(k  ,isp)) )
      sedgvar(k,isn,nvi,ipn)=.5*(.5*(Tgt1(k-1,isn)+Tgt1(k-1 ,isp)) &
                         +       .5*(Tgt1(k  ,isn)+Tgt1(k   ,isp)))
      sedgvar(k,isn,nv2,ipn)=.5*(.5*(Tgt2(k-1,isn)+Tgt2(k-1 ,isp)) &
                         +       .5*(Tgt2(k  ,isn)+Tgt2(k   ,isp)))
      sedgvar(k,isn,nv3,ipn)=.5*(.5*(Tgt3(k-1,isn)+Tgt3(k-1 ,isp)) &
                         +       .5*(Tgt3(k  ,isn)+Tgt3(k   ,isp)))
    end do !  k -loop
    sedgvar( 1,isn,  1,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,  1,ipn) &
                           +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,  1,ipn)
    sedgvar( 1,isn,  2,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,  2,ipn) &
                           +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,  2,ipn)
    sedgvar( 1,isn,nvi,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,nvi,ipn) &
                           +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,nvi,ipn)
    sedgvar( 1,isn,nv2,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,nv2,ipn) &
                           +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,nv2,ipn)
    sedgvar( 1,isn,nv3,ipn)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*sedgvar(2,isn,nv3,ipn) &
                           +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*sedgvar(3,isn,nv3,ipn)
    sedgvar(nz,isn,  1,ipn)=2.*sedgvar(nz-1,isn,  1,ipn)-sedgvar(nz-2,isn,  1,ipn)
    sedgvar(nz,isn,  2,ipn)=2.*sedgvar(nz-1,isn,  2,ipn)-sedgvar(nz-2,isn,  2,ipn)
    sedgvar(nz,isn,nvi,ipn)=2.*sedgvar(nz-1,isn,nvi,ipn)-sedgvar(nz-2,isn,nvi,ipn)
    sedgvar(nz,isn,nv2,ipn)=2.*sedgvar(nz-1,isn,nv2,ipn)-sedgvar(nz-2,isn,nv2,ipn)
    sedgvar(nz,isn,nv3,ipn)=2.*sedgvar(nz-1,isn,nv3,ipn)-sedgvar(nz-2,isn,nv3,ipn)
!$acc loop vector
    do k = 1,nz
      sedgvar(k,isn,5,ipn) = teb(k,isn,ipn)+sedgvar(k,isn,5,ipn)
    enddo
  end do  ! isn-loop

  do k=1,nz-1
    bedgvar(k,ipn,  1)=ca4k(k)*  u(k,ipn)+ca4p(k)*  u(k+1,ipn)
    bedgvar(k,ipn,  2)=ca4k(k)*  v(k,ipn)+ca4p(k)*  v(k+1,ipn)
    bedgvar(k,ipn,nvi)=ca4k(k)* rp(k,ipn)+ca4p(k)* rp(k+1,ipn)
    bedgvar(k,ipn,nv2)=ca4k(k)* tp(k,ipn)+ca4p(k)* tp(k+1,ipn)
    bedgvar(k,ipn,nv3)=ca4k(k)*trp(k,ipn)+ca4p(k)*trp(k+1,ipn)
  end do
    bedgvar( 0, ipn,  1)=( z(0,ipn)- z(2,ipn))/( z(1,ipn)- z(2,ipn))*bedgvar(1,ipn,  1) &
                        +( z(0,ipn)- z(1,ipn))/( z(2,ipn)- z(1,ipn))*bedgvar(2,ipn,  1)
    bedgvar( 0, ipn,  2)=( z(0,ipn)- z(2,ipn))/( z(1,ipn)- z(2,ipn))*bedgvar(1,ipn,  2) &
                        +( z(0,ipn)- z(1,ipn))/( z(2,ipn)- z(1,ipn))*bedgvar(2,ipn,  2)
  bedgvar( 0, ipn,nvi)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*bedgvar(1,ipn,nvi) &
                      +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*bedgvar(2,ipn,nvi) 
  bedgvar( 0, ipn,nv2)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*bedgvar(1,ipn,nv2) &
                      +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*bedgvar(2,ipn,nv2) 
  bedgvar( 0, ipn,nv3)=(zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*bedgvar(1,ipn,nv3) &
                      +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*bedgvar(2,ipn,nv3) 
  bedgvar(nz,ipn,  1)=ca4k(nz)*  u(nz,ipn)+ca4p(nz)*  u(nz,ipn)
  bedgvar(nz,ipn,  2)=ca4k(nz)*  v(nz,ipn)+ca4p(nz)*  v(nz,ipn)
  bedgvar(nz,ipn,nvi)=ca4k(nz)* rp(nz,ipn)+ca4p(nz)* rp(nz,ipn)
  bedgvar(nz,ipn,nv2)=ca4k(nz)* tp(nz,ipn)+ca4p(nz)* tp(nz,ipn)
  bedgvar(nz,ipn,nv3)=ca4k(nz)*trp(nz,ipn)+ca4p(nz)*trp(nz,ipn)
enddo !ipn-loop
!$OMP END PARALLEL DO

return
end subroutine vdmintsv5

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
#include "solveiThLS5.F90"
