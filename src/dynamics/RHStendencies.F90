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

subroutine RHStendencies(ims,ime,ips,ipe,spncef,su,sv,tfr,tfur,tfvr,tfwr,tftr)
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent (IN ) :: ims,ime,ips,ipe
real(rt),intent (IN) :: spncef(0:nz)   
real(rt),intent (IN) :: su (  nz,ims:ime)   
real(rt),intent (IN) :: sv (  nz,ims:ime)   
real(rt),intent (OUT):: tfr (  nz,ims:ime)  ! forcing of density (r)
real(rt),intent (OUT):: tfur(  nz,ims:ime)  ! forcing of ur
real(rt),intent (OUT):: tfvr(  nz,ims:ime)  ! forcing of vr
real(rt),intent (OUT):: tfwr(0:nz,ims:ime)  ! forcing of wr
real(rt),intent (OUT):: tftr(  nz,ims:ime)  ! forcing of tr
integer              :: k,ipn

!$acc parallel vector_length(96)
!$OMP PARALLEL DO PRIVATE(k)
!$acc loop gang
do ipn=ips,ipe
  tfwr(0,ipn) = 0._rt
!$acc loop vector
  do k=1,nz
    tfur(k,ipn) = (1._rt-.5_rt*(spncef(k-1)+spncef(k)))*su(k,ipn)
    tfvr(k,ipn) = (1._rt-.5_rt*(spncef(k-1)+spncef(k)))*sv(k,ipn)
    tftr(k,ipn) = 0._rt
    tfr(k,ipn)  = 0._rt
    tfwr(k,ipn) = 0._rt
  end do
end do
!$acc end parallel
!$OMP END PARALLEL DO
return

end subroutine RHStendencies

!**************************************************************************
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
!**************************************************************************
