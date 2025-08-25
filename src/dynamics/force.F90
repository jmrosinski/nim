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

subroutine force(ims,ime,ips,ipe,vol,f,ur,vr,fu0,fv0,fw0,tfur,tfvr,tfwr)
!*********************************************************************
!Calculates NIM dynamical forcing terms
!
!       Nonhydrostatic Icosahedral Model (NIM)
!  Original subroutine force in QNH:  A.E. MacDonald, Jin Lee -1991
!       Jin Lee                  December, 2008
!  Jacques Middlecoff July 2009 Put K on the inside for the GPU
!  Jacques Middlecoff July 2009 Put variables in calling sequence.
!       Mark Govett August 2010 GPU parallelization
!  Jacques Middlecoff Oct. 2010 Removed SMS
!  Jacques Middlecoff July 2016 Split force into preForce and force.
!*********************************************************************

use gptl
#ifdef _OPENACC
use gptl_acc
#endif
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer ,intent(IN   ) ::           ims,ime,ips,ipe
real(rt),intent(IN   ) :: vol (  NZ,ims:ime) ! volume for each control volume
real(rt),intent(IN   ) :: f   (     ims:ime) ! Coriolis parameters
real(rt),intent(IN   ) :: ur  (  NZ,ims:ime) ! u*r, flux variable u times density
real(rt),intent(IN   ) :: vr  (  NZ,ims:ime) ! v*r, flux variable v times density
real(rt),intent(IN   ) :: fu0 (  NZ,ims:ime)
real(rt),intent(IN   ) :: fv0 (  NZ,ims:ime)
real(rt),intent(IN   ) :: fw0 (  NZ,ims:ime)
real(rt),intent(INOUT) :: tfur(  NZ,ims:ime) ! tendency function of ur
real(rt),intent(INOUT) :: tfvr(  NZ,ims:ime) ! tendency function of vr
real(rt),intent(INOUT) :: tfwr(0:NZ,ims:ime) ! tendency function of wr
integer                :: ipn,k
integer                :: ret, ret2
integer                :: force_handle
integer                :: force_ipn_handle

#ifdef _OPENACC
!$acc parallel private(ret) copyout(force_handle,force_ipn_handle)
  ret = gptlinit_handle_gpu ('force', force_handle)
  ret = gptlinit_handle_gpu ('force_ipn', force_ipn_handle)
  ret = gptlstart_gpu (force_handle)
!$acc end parallel
#endif

!$OMP PARALLEL DO PRIVATE(k)
!$acc parallel private(ret2) num_workers(PAR_WRK) vector_length(VEC_LEN) &
!$acc&         copyin (force_ipn_handle)  
!$acc loop gang
do ipn=ips,ipe
#ifdef _OPENACC
  ret2 = gptlstart_gpu(force_ipn_handle)
#endif
 !$acc loop vector
  do k=1,NZ
    tfur(k,ipn) = tfur(k,ipn) + fu0(k,ipn)/vol(k,ipn) + f(ipn)*vr(k,ipn)
    tfvr(k,ipn) = tfvr(k,ipn) + fv0(k,ipn)/vol(k,ipn) - f(ipn)*ur(k,ipn)
    tfwr(k,ipn) = tfwr(k,ipn) + fw0(k,ipn)/vol(k,ipn)
  end do
#ifdef _OPENACC
  ret2 = gptlstop_gpu(force_ipn_handle)
#endif
end do
!$acc end parallel
!$OMP END PARALLEL DO

#ifdef _OPENACC
!$acc parallel private(ret) copyin (force_handle)
ret = gptlstop_gpu (force_handle)
!$acc end parallel
#endif

return
end subroutine force

!*********************************************************************
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)!
!
!*********************************************************************
