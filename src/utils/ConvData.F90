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

program ConvData
!
! This program is created temporarily to convert g3d.dat & amx.dat 
! from old format to new formate which allows high resolution runs. 
! 10/31/12   Jin Lee, Ka Yee Wong
!
implicit none
integer, parameter :: glvl=4
integer, parameter :: nz=32
!
integer :: nip,npp,k,ipn,isn,i,nob,nbf
real   ,allocatable :: lat(:),lon(:)
real   ,allocatable :: z(:,:),zm(:,:)
real   ,allocatable :: prox(:,:),nprox(:),proxs(:,:)
real   ,allocatable :: area(:)
real   ,allocatable :: cs(:,:),sn(:,:)
real   ,allocatable :: xyc(:,:,:),xys(:,:,:),xyb(:,:,:)
real   ,allocatable :: zc(:,:,:),zcs(:,:,:)
real   ,allocatable :: nvecs(:,:,:,:),nvecb(:,:,:)
real   ,allocatable :: sa(:,:,:),vol(:,:),saoprx(:,:,:)
real   ,allocatable :: ca4k(:),ca4p(:),caw4tk(:),caw4tp(:),ct4w(:)
real   ,allocatable :: amtx1(:,:,:)
real   ,allocatable :: amtx2(:,:,:)
!
nip = 10*(2**glvl)**2+2
npp = 6
nob = 9
nbf = 9
allocate(lat(nip),lon(nip)) 
allocate(z(0:nz,nip),zm(nz,nip))
allocate(prox(npp,nip),nprox(nip),proxs(npp,nip))
allocate(area(nip))
allocate(cs(npp,nip),sn(npp,nip))
allocate(xyc(2,npp,nip),xys(2,npp,nip),xyb(2,npp,nip))
allocate(zc(0:nz,npp,nip),zcs(nz,npp,nip))
allocate(nvecs(nz,npp,nip,2),nvecb(0:nz,nip,3))
allocate(sa(nz,npp,nip),vol(nz,nip),saoprx(nz,npp,nip))
allocate(ca4k(1:nz),ca4p(1:nz),caw4tk(1:nz),caw4tp(1:nz),ct4w(1:nz))
allocate(amtx1(  nz  ,(nob+1)*nbf,nip))
allocate(amtx2(0:nz-1,(nob+1)*nbf,nip))
!
open(unit=28,file="/scratch1/portfolios/BMC/nim/DATADIR/DATAtmp/CTL/G3D_G04K032.dat.old", form="unformatted")
print *,'  Reading g3d.dat'
read(28) lat,lon,z,zm,prox,nprox,proxs,area,cs,sn,xyc,xys,xyb, &
          zc,zcs,nvecs,nvecb,sa,vol,saoprx,ca4k,ca4p,caw4tk,caw4tp,ct4w
close(28)
print*,'reading g3d.dat completed!'
print *,'Start writing g3d.dat'
open(unit=28,file="./g3d.dat",form="unformatted")
write(28) lat
write(28) lon
do k=0,nz
  write(28) z(k,:)
enddo
do k=1,nz
  write(28) zm(k,:)
enddo
do isn=1,npp
  write(28) prox(isn,:) 
enddo
write(28) nprox
do isn=1,npp
  write(28)  proxs(isn,:)
enddo
write(28) area
do isn=1,npp
  write(28) cs(isn,:)
enddo
do isn=1,npp
  write(28) sn(isn,:)
enddo
do isn=1,npp
  do i=1,2
    write(28) xyc(i,isn,:)
  enddo
enddo
do isn=1,npp
  do i=1,2
    write(28) xys(i,isn,:)
  enddo
enddo
do isn=1,npp
  do i=1,2
    write(28) xyb(i,isn,:)
  enddo
enddo
do isn=1,npp
  do k=0,nz
    write(28) zc(k,isn,:)
  enddo
enddo
do isn=1,npp
  do k=1,nz
    write(28) zcs(k,isn,:)
  enddo
 enddo
do i=1,2
  do isn=1,npp
    do k=1,nz
      write(28) nvecs(k,isn,:,i)
    enddo
  enddo
enddo
do i=1,3
  do k=0,nz
    write(28) nvecb(k,:,i)
  enddo
enddo
do isn=1,npp
  do k=1,nz
    write(28) sa(k,isn,:)
  enddo
enddo
do k=1,nz
  write(28) vol(k,:)
enddo
do isn=1,npp
  do k=1,nz
    write(28) saoprx(k,isn,:)
  enddo
enddo
write(28) ca4k,ca4p,caw4tk,caw4tp,ct4w
close(28)
print*,'Writing g3d.dat completed!'

open(unit=28,file="/scratch1/portfolios/BMC/nim/DATADIR/DATAtmp/CTL/AM_G04K032.dat.old", form="unformatted")
print *,'  Reading amx.dat'
read(28) amtx1,amtx2
close(28)
print *,'  Reading amx.dat completed'
print *,'Start writing amx.dat'
open(unit=28,file="./amx.dat",form="unformatted")
do i=1,(nob+1)*nbf
  do k=1,nz
     write(28) amtx1(k,i,:)
  enddo
enddo

do i=1,(nob+1)*nbf
  do k=0,nz-1
     write(28) amtx2(k,i,:)
  enddo
enddo
close(28)
print *,'Writing amx.dat completed'

end program ConvData
