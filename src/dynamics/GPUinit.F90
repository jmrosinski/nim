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

! Routine to initialize the GPU
! Author:  Jacques Middlecoff
! Date:  September 2010 
! For Fortran this routine does nothing except return GPUrun=0 and error=0
! See dynamics/cuda/GPUinit.cu for the real routine.

subroutine GPUinit(myRank,max_accelerators_per_node,GPUrun,error)
#ifdef _OPENACC
!SMS$insert use mpi
use openacc
#endif
integer, intent(IN ) :: myRank
integer, intent(IN ) :: max_accelerators_per_node
logical, intent(OUT) :: GPUrun
integer, intent(OUT) :: error
integer              :: devicenum,ngpus

error = 0
#ifdef _OPENACC
  GPUrun = .true.
  call acc_init(ACC_DEVICE_NVIDIA)
  ngpus = acc_get_num_devices(ACC_DEVICE_NVIDIA)
  if (ngpus .eq. 0) then
    print *,'No GPUs found on this system.  Exiting'
    error = -5
    return
  endif
!  if (max_accelerators_per_node > ngpus ) then
!    print *,'max_accelerators_per_node > ngpus.  Exiting',max_accelerators_per_node,ngpus
!    error = -7
!    return
!  endif
! devicenum = mod(myRank,max_accelerators_per_node) + 4
  devicenum = mod(myRank,max_accelerators_per_node)
  call acc_set_device_num(devicenum,ACC_DEVICE_NVIDIA)
  print *,'GPUinit: rank,device   ',myRank,acc_get_device_num(ACC_DEVICE_NVIDIA)
#else
  GPUrun = .false.
#endif

return
end subroutine GPUinit
