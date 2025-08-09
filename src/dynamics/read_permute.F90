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


module permute_io

  implicit none

  interface read_permute
    module procedure read_permute_r41d
    module procedure read_permute_r81d
    module procedure read_permute_i1d
    module procedure read_permute_r42d
    module procedure read_permute_r82d
    module procedure read_permute_i2d
    module procedure read_permute_r43d
    module procedure read_permute_r83d
    module procedure read_permute_i3d
  end interface read_permute

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read_permute: reads in, permutes, and writes out to distributed memory the
!               contents of a single distributed array.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!TODO:  use cpp to DRY this code!  

subroutine read_permute_r41d (unitno, outarr, FileName, ArrayName)
  use module_control  , only: nip
  use module_variables, only: globalperm

  implicit none
#include <gptl.inc>

  integer, intent(in) :: unitno
  character(len=*), intent(in) :: FileName
  character(len=*), intent(in) :: ArrayName
!SMS$DISTRIBUTE(dh,1) BEGIN
  real*4, intent(out) :: outarr(:)
!SMS$DISTRIBUTE END

  real*4,allocatable :: fromdisk(:)
  integer :: ipn
  integer :: ret,ierr

  ret = gptlstart ('read_permute')
!SMS$SERIAL (<outarr,OUT>:default=ignore) BEGIN
  allocate(fromdisk(nip))
  read(unitno, iostat=ierr) fromdisk
  if (ierr.ne.0) then
    write (*,'(4a)') 'read_permute: error reading ',trim(ArrayName),' from ',trim(Filename)
    stop
  endif
! Currently require VECTOR ALWAYS to vectorize this loop, but it doesnt help much.
! So leave as is for now since always runs on CPU.
  do ipn=1,nip
    outarr(ipn) = fromdisk(globalperm(ipn))
  enddo
  deallocate(fromdisk)
!SMS$SERIAL END
  ret = gptlstop ('read_permute')
  return
end subroutine read_permute_r41d


subroutine read_permute_r81d (unitno, outarr, FileName, ArrayName)
  use module_control  , only: nip
  use module_variables, only: globalperm

  implicit none
#include <gptl.inc>

  integer, intent(in) :: unitno
  character(len=*), intent(in) :: FileName
  character(len=*), intent(in) :: ArrayName
!SMS$DISTRIBUTE(dh,1) BEGIN
  real*8, intent(out) :: outarr(:)
!SMS$DISTRIBUTE END

  real*4,allocatable :: fromdisk(:)
  integer :: ipn
  integer :: ret,ierr

  ret = gptlstart ('read_permute')
!SMS$SERIAL (<outarr,OUT>:default=ignore) BEGIN
  allocate(fromdisk(nip))
  read(unitno, iostat=ierr) fromdisk
  if (ierr.ne.0) then
    write (*,'(4a)') 'read_permute: error reading ',trim(ArrayName),' from ',trim(Filename)
    stop
  endif
! Currently require VECTOR ALWAYS to vectorize this loop, but it doesnt help much.
! So leave as is for now since always runs on CPU.
  do ipn=1,nip
    outarr(ipn) = fromdisk(globalperm(ipn))
  enddo
  deallocate(fromdisk)
!SMS$SERIAL END
  ret = gptlstop ('read_permute')
  return
end subroutine read_permute_r81d


subroutine read_permute_i1d (unitno, outarr, FileName, ArrayName)
  use module_control  , only: nip
  use module_variables, only: globalperm

  implicit none
#include <gptl.inc>

  integer, intent(in) :: unitno
  character(len=*), intent(in) :: FileName
  character(len=*), intent(in) :: ArrayName
!SMS$DISTRIBUTE(dh,1) BEGIN
  integer, intent(out) :: outarr(:)
!SMS$DISTRIBUTE END

  integer,allocatable :: fromdisk(:)
  integer :: ipn
  integer :: ret,ierr

  ret = gptlstart ('read_permute')
!SMS$SERIAL (<outarr,OUT>:default=ignore) BEGIN
  allocate(fromdisk(nip))
  read(unitno, iostat=ierr) fromdisk
  if (ierr.ne.0) then
    write (*,'(4a)') 'read_permute: error reading ',trim(ArrayName),' from ',trim(Filename)
    stop
  endif
! Currently require VECTOR ALWAYS to vectorize this loop, but it doesnt help much.
! So leave as is for now since always runs on CPU.
  do ipn=1,nip
    outarr(ipn) = fromdisk(globalperm(ipn))
  enddo
  deallocate(fromdisk)
!SMS$SERIAL END
  ret = gptlstop ('read_permute')
  return
end subroutine read_permute_i1d


subroutine read_permute_r42d (unitno, outarr, FileName, ArrayName)
  use module_control  , only: nip
  use module_variables, only: globalperm

  implicit none
#include <gptl.inc>

  integer, intent(in) :: unitno
  character(len=*), intent(in) :: FileName
  character(len=*), intent(in) :: ArrayName
!SMS$DISTRIBUTE(dh,2) BEGIN
  real*4, intent(out) :: outarr(:,:)
!SMS$DISTRIBUTE END

  real*4,allocatable :: fromdisk(:,:)
  integer :: ipn, k, nz
  integer :: ret,ierr
  ret = gptlstart ('read_permute')
  nz = size(outarr,1)
!SMS$SERIAL (<outarr,OUT>:default=ignore) BEGIN
  allocate(fromdisk(nz,nip))
  read(unitno, iostat=ierr) fromdisk
  if (ierr.ne.0) then
    write (*,'(4a)') 'read_permute: error reading ',trim(ArrayName),' from ',trim(Filename)
    stop
  endif
! Currently require VECTOR ALWAYS to vectorize this loop, but it doesnt help much.
! So leave as is for now since always runs on CPU.
  do ipn=1,nip
    do k=1,nz
      outarr(k,ipn) = fromdisk(k,globalperm(ipn))
    enddo
  enddo
  deallocate(fromdisk)
!SMS$SERIAL END
  ret = gptlstop ('read_permute')
  return
end subroutine read_permute_r42d


subroutine read_permute_r82d (unitno, outarr, FileName, ArrayName)
  use module_control  , only: nip
  use module_variables, only: globalperm

  implicit none
#include <gptl.inc>

  integer, intent(in) :: unitno
  character(len=*), intent(in) :: FileName
  character(len=*), intent(in) :: ArrayName
!SMS$DISTRIBUTE(dh,2) BEGIN
  real*8, intent(out) :: outarr(:,:)
!SMS$DISTRIBUTE END

  real*4,allocatable :: fromdisk(:,:)
  integer :: ipn, k, nz
  integer :: ret,ierr

  ret = gptlstart ('read_permute')
  nz = size(outarr,1)
!SMS$SERIAL (<outarr,OUT>:default=ignore) BEGIN
  allocate(fromdisk(nz,nip))
  read(unitno, iostat=ierr) fromdisk
  if (ierr.ne.0) then
    write (*,'(4a)') 'read_permute: error reading ',trim(ArrayName),' from ',trim(Filename)
    stop
  endif
! Currently require VECTOR ALWAYS to vectorize this loop, but it doesnt help much.
! So leave as is for now since always runs on CPU.
  do ipn=1,nip
    do k=1,nz
      outarr(k,ipn) = fromdisk(k,globalperm(ipn))
    enddo
  enddo
  deallocate(fromdisk)
!SMS$SERIAL END
  ret = gptlstop ('read_permute')
  return
end subroutine read_permute_r82d


subroutine read_permute_i2d (unitno, outarr, FileName, ArrayName)
  use module_control  , only: nip
  use module_variables, only: globalperm

  implicit none
#include <gptl.inc>

  integer, intent(in) :: unitno
  character(len=*), intent(in) :: FileName
  character(len=*), intent(in) :: ArrayName
!SMS$DISTRIBUTE(dh,2) BEGIN
  integer, intent(out) :: outarr(:,:)
!SMS$DISTRIBUTE END

  integer,allocatable :: fromdisk(:,:)
  integer :: ipn, k, nz
  integer :: ret,ierr

  ret = gptlstart ('read_permute')
  nz = size(outarr,1)
!SMS$SERIAL (<outarr,OUT>:default=ignore) BEGIN
  allocate(fromdisk(nz,nip))
  read(unitno, iostat=ierr) fromdisk
  if (ierr.ne.0) then
    write (*,'(4a)') 'read_permute: error reading ',trim(ArrayName),' from ',trim(Filename)
    stop
  endif
! Currently require VECTOR ALWAYS to vectorize this loop, but it doesnt help much.
! So leave as is for now since always runs on CPU.
  do ipn=1,nip
    do k=1,nz
      outarr(k,ipn) = fromdisk(k,globalperm(ipn))
    enddo
  enddo
  deallocate(fromdisk)
!SMS$SERIAL END
  ret = gptlstop ('read_permute')
  return
end subroutine read_permute_i2d


subroutine read_permute_r43d (unitno, outarr, FileName, ArrayName)
  use module_control  , only: nip
  use module_variables, only: globalperm

  implicit none
#include <gptl.inc>

  integer, intent(in) :: unitno
  character(len=*), intent(in) :: FileName
  character(len=*), intent(in) :: ArrayName
!SMS$DISTRIBUTE(dh,3) BEGIN
  real*4, intent(out) :: outarr(:,:,:)
!SMS$DISTRIBUTE END

  real*4,allocatable :: fromdisk(:,:,:)
  integer :: ipn, k, nz, isn, npp
  integer :: ret,ierr

  ret = gptlstart ('read_permute')
  nz = size(outarr,1)
  npp = size(outarr,2)
!SMS$SERIAL (<outarr,OUT>:default=ignore) BEGIN
  allocate(fromdisk(nz,npp,nip))
  read(unitno, iostat=ierr) fromdisk
  if (ierr.ne.0) then
    write (*,'(4a)') 'read_permute: error reading ',trim(ArrayName),' from ',trim(Filename)
    stop
  endif
! Currently require VECTOR ALWAYS to vectorize this loop, but it doesnt help much.
! So leave as is for now since always runs on CPU.
  do ipn=1,nip
    do isn=1,npp
      do k=1,nz
        outarr(k,isn,ipn) = fromdisk(k,isn,globalperm(ipn))
      enddo
    enddo
  enddo
  deallocate(fromdisk)
!SMS$SERIAL END
  ret = gptlstop ('read_permute')
  return
end subroutine read_permute_r43d


subroutine read_permute_r83d (unitno, outarr, FileName, ArrayName)
  use module_control  , only: nip
  use module_variables, only: globalperm

  implicit none
#include <gptl.inc>

  integer, intent(in) :: unitno
  character(len=*), intent(in) :: FileName
  character(len=*), intent(in) :: ArrayName
!SMS$DISTRIBUTE(dh,3) BEGIN
  real*8, intent(out) :: outarr(:,:,:)
!SMS$DISTRIBUTE END

  real*4,allocatable :: fromdisk(:,:,:)
  integer :: ipn, k, nz, isn, npp
  integer :: ret,ierr

  ret = gptlstart ('read_permute')
  nz = size(outarr,1)
  npp = size(outarr,2)
!SMS$SERIAL (<outarr,OUT>:default=ignore) BEGIN
  allocate(fromdisk(nz,npp,nip))
  read(unitno, iostat=ierr) fromdisk
  if (ierr.ne.0) then
    write (*,'(4a)') 'read_permute: error reading ',trim(ArrayName),' from ',trim(Filename)
    stop
  endif
! Currently require VECTOR ALWAYS to vectorize this loop, but it doesnt help much.
! So leave as is for now since always runs on CPU.
  do ipn=1,nip
    do isn=1,npp
      do k=1,nz
        outarr(k,isn,ipn) = fromdisk(k,isn,globalperm(ipn))
      enddo
    enddo
  enddo
  deallocate(fromdisk)
!SMS$SERIAL END
  ret = gptlstop ('read_permute')
  return
end subroutine read_permute_r83d


subroutine read_permute_i3d (unitno, outarr, FileName, ArrayName)
  use module_control  , only: nip
  use module_variables, only: globalperm

  implicit none
#include <gptl.inc>

  integer, intent(in) :: unitno
  character(len=*), intent(in) :: FileName
  character(len=*), intent(in) :: ArrayName
!SMS$DISTRIBUTE(dh,3) BEGIN
  integer, intent(out) :: outarr(:,:,:)
!SMS$DISTRIBUTE END

  integer,allocatable :: fromdisk(:,:,:)
  integer :: ipn, k, nz, isn, npp
  integer :: ret,ierr

  ret = gptlstart ('read_permute')
  nz = size(outarr,1)
  npp = size(outarr,2)
!SMS$SERIAL (<outarr,OUT>:default=ignore) BEGIN
  allocate(fromdisk(nz,npp,nip))
  read(unitno, iostat=ierr) fromdisk
  if (ierr.ne.0) then
    write (*,'(4a)') 'read_permute: error reading ',trim(ArrayName),' from ',trim(Filename)
    stop
  endif
! Currently require VECTOR ALWAYS to vectorize this loop, but it doesnt help much.
! So leave as is for now since always runs on CPU.
  do ipn=1,nip
    do isn=1,npp
      do k=1,nz
        outarr(k,isn,ipn) = fromdisk(k,isn,globalperm(ipn))
      enddo
    enddo
  enddo
  deallocate(fromdisk)
!SMS$SERIAL END
  ret = gptlstop ('read_permute')
  return
end subroutine read_permute_i3d


end module permute_io

