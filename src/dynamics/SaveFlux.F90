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

subroutine SaveFlux(ims,ime,ips,ipe,ur,vr,wr,trp, r,urs,vrs,wrs,trs, rs)
! Jacques Middlecoff October 2010 Removed SMS
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN ) :: ims,ime,ips,ipe
real(rt),intent(IN ):: ur (  NZ,ims:ime) ! u*r, flux variable u times density
real(rt),intent(IN ):: vr (  NZ,ims:ime) ! v*r, flux variable v times density
real(rt),intent(IN ):: wr (0:NZ,ims:ime) ! w*r, flux variable w times density
real(rt),intent(IN ):: trp(  NZ,ims:ime) ! flux pot temp prime ( = tr - trb)	
real(rt),intent(IN )::  r (  NZ,ims:ime) ! rho (density) prime (kg/m**3) (= r - rb)
real(rt),intent(OUT):: urs(  NZ,ims:ime)	
real(rt),intent(OUT):: vrs(  NZ,ims:ime)	
real(rt),intent(OUT):: wrs(0:NZ,ims:ime)	
real(rt),intent(OUT):: trs(  NZ,ims:ime)	
real(rt),intent(OUT)::  rs(  NZ,ims:ime)	
integer :: k,ipn

!$OMP PARALLEL DO PRIVATE(k)
!$acc parallel vector_length(96)
!$acc loop gang
do ipn=ips,ipe
  wrs(0,ipn) = wr(0,ipn)
!$acc loop vector
  do k=1,NZ
    urs(k,ipn) =  ur(k,ipn)
    vrs(k,ipn) =  vr(k,ipn)
    wrs(k,ipn) =  wr(k,ipn)
    trs(k,ipn) = trp(k,ipn)
     rs(k,ipn) =   r(k,ipn)
  end do
end do
!$acc end parallel
!$OMP END PARALLEL DO
return
end subroutine SaveFlux

!***************************************************************************
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
!***************************************************************************
